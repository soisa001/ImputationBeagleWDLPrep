version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    String? disk_type
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

workflow GLIMPSE2Summarize {
    input {
        File panel_vcf              # split to biallelic
        File panel_vcf_idx
        File imputed_vcf            # split to biallelic, variants in same order as in panel
        File imputed_vcf_idx
        File population_tsv         # should contain columns: Sample name, Population code (e.g., igsr_samples.tsv tables generated from https://www.internationalgenome.org/data/)
        String output_prefix
        File? summarize_wheelhouse  # pre-staged pip wheelhouse (cp311 manylinux); perimeter blocks the conda/PyPI install
    }

    call SummarizeAndPlot { input:
        panel_vcf = panel_vcf,
        panel_vcf_idx = panel_vcf_idx,
        imputed_vcf = imputed_vcf,
        imputed_vcf_idx = imputed_vcf_idx,
        population_tsv = population_tsv,
        summarize_wheelhouse = summarize_wheelhouse,
        output_prefix = output_prefix
    }

    output {
        File summarize_pearson_tsv = SummarizeAndPlot.summarize_pearson_tsv
        Array[File] summarize_plots = SummarizeAndPlot.summarize_plots
    }
}

task SummarizeAndPlot {
    input {
        File panel_vcf
        File panel_vcf_idx
        File imputed_vcf
        File imputed_vcf_idx
        File population_tsv
        String output_prefix
        File? summarize_wheelhouse

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 20 + ceil(size(panel_vcf, "GiB") + size(imputed_vcf, "GiB"))

    command <<<
        set -euxo pipefail
        
        export MPLBACKEND=Agg     # headless matplotlib

        # Install deps for variant streaming, stats, and plotting. The perimeter blocks anaconda/PyPI
        # from the task, so install OFFLINE from a pre-staged pip wheelhouse (cp311 manylinux) if one
        # is provided; otherwise fall back to an online pip install.
        # matplotlib pinned <3.10: seaborn 0.13.2's boxplot legend (_configure_legend) raises
        # "UnboundLocalError: boxprops" with matplotlib >=3.11. Pinning here (not just in the wheelhouse
        # build) means a stale/mixed wheelhouse fails LOUDLY at install instead of silently grabbing 3.11.
        if [ -n "~{summarize_wheelhouse}" ]; then
            mkdir -p wheelhouse && tar -xzf ~{summarize_wheelhouse} -C wheelhouse
            pip install --no-index --find-links=wheelhouse cyvcf2 pandas numpy 'matplotlib<3.10' seaborn scipy
        else
            pip install cyvcf2 pandas numpy 'matplotlib<3.10' seaborn scipy
        fi
        python -c 'import matplotlib; assert tuple(map(int, matplotlib.__version__.split(".")[:2])) < (3,10), "matplotlib "+matplotlib.__version__+" breaks seaborn 0.13.2 boxplot; rebuild the wheelhouse (pull the matplotlib<3.10 fix)"'

        python - ~{panel_vcf} \
                 ~{imputed_vcf} \
                 ~{population_tsv} \
                 ~{output_prefix} <<-'EOF'
        import sys
        import time
        import numpy as np
        import pandas as pd
        import cyvcf2
        import matplotlib.pyplot as plt
        import seaborn as sns
        import matplotlib.colors
        from scipy.stats import pearsonr
        from datetime import timedelta
        import warnings

        warnings.filterwarnings('ignore', category=RuntimeWarning)

        panel_vcf_path = sys.argv[1]
        imputed_vcf_path = sys.argv[2]
        population_tsv_path = sys.argv[3]
        output_prefix = sys.argv[4]

        # 1. Initialize VCF readers
        panel_vcf = cyvcf2.VCF(panel_vcf_path)
        imputed_vcf = cyvcf2.VCF(imputed_vcf_path)

        panel_samples = np.array(panel_vcf.samples)
        target_samples = np.array(imputed_vcf.samples)
        num_p_samples = len(panel_samples)
        num_c_samples = len(target_samples)

        # 2. Global Storage for Metrics
        altlen_all = []
        panel_af_all, target_af_all = [], []

        panel_hwe = {'hom_ref': [], 'het': [], 'hom_alt': []}
        target_hwe = {'hom_ref': [], 'het': [], 'hom_alt': []}

        panel_mean_alt_alleles_all, target_mean_alt_alleles_all = [], []

        # Pre-allocate per-sample accumulators (fast 1D numpy integer arrays)
        p_het_all = np.zeros(num_p_samples, dtype=int)
        p_hom_ref_all = np.zeros(num_p_samples, dtype=int)
        p_hom_alt_all = np.zeros(num_p_samples, dtype=int)
        p_het_ins = np.zeros(num_p_samples, dtype=int)
        p_het_del = np.zeros(num_p_samples, dtype=int)

        c_het_all = np.zeros(num_c_samples, dtype=int)
        c_hom_ref_all = np.zeros(num_c_samples, dtype=int)
        c_hom_alt_all = np.zeros(num_c_samples, dtype=int)
        c_het_ins = np.zeros(num_c_samples, dtype=int)
        c_het_del = np.zeros(num_c_samples, dtype=int)

        # Pre-allocated boolean buffers for the loop (avoids millions of memory allocations)
        p_is_het = np.empty(num_p_samples, dtype=bool)
        p_is_hom_ref = np.empty(num_p_samples, dtype=bool)
        p_is_hom_alt = np.empty(num_p_samples, dtype=bool)

        c_is_het = np.empty(num_c_samples, dtype=bool)
        c_is_hom_ref = np.empty(num_c_samples, dtype=bool)
        c_is_hom_alt = np.empty(num_c_samples, dtype=bool)

        # 3. Stream variants with a (CHROM,POS,REF,ALT) merge-join.
        #    The site-matched panel and the imputed output are both coordinate-sorted, but at a
        #    multiallelic position they need not carry the same allele SET or ORDER (e.g. panel
        #    chr1:10626 A>AG vs imputed chr1:10626 A>AGGCGCAG), so a positional zip() desyncs and
        #    raises "VCFs are out of sync". Instead we walk both as position-blocks and pair records
        #    by exact (REF,ALT) within each position -- order-independent, tolerant of allele-set
        #    differences, memory-bounded (one position block per file at a time). Only matched sites
        #    contribute to the AF correlation (the only places where it is meaningful anyway);
        #    panel-only / target-only sites are counted and reported.
        print("Streaming variants (merge-join on CHROM,POS,REF,ALT)...", flush=True)
        variants_processed = 0
        panel_only = 0
        target_only = 0
        start_time = time.time()

        # contig rank for cross-file position ordering (assumes a shared reference)
        try:
            chrom_rank = {c: i for i, c in enumerate(imputed_vcf.seqnames)}
        except Exception:
            chrom_rank = {}
        def pos_order(chrom, pos):
            return (chrom_rank.get(chrom, 1 << 30), pos)

        def pos_blocks(vcf):
            # yield ((CHROM, POS), [variants at that position]); input is coordinate-sorted
            block = []
            key = None
            for v in vcf:
                k = (v.CHROM, v.POS)
                if key is not None and k != key:
                    yield key, block
                    block = []
                block.append(v)
                key = k
            if block:
                yield key, block

        def allele_key(v):
            return (v.REF, v.ALT[0] if v.ALT else ".")

        pit = pos_blocks(panel_vcf)
        cit = pos_blocks(imputed_vcf)
        pb = next(pit, None)
        cb = next(cit, None)

        while pb is not None and cb is not None:
            p_key, p_block = pb
            c_key, c_block = cb
            if p_key == c_key:
                p_amap = {}
                for v in p_block:
                    p_amap.setdefault(allele_key(v), v)
                c_amap = {}
                for v in c_block:
                    c_amap.setdefault(allele_key(v), v)

                for ak, c_var in c_amap.items():
                    p_var = p_amap.get(ak)
                    if p_var is None:
                        target_only += 1
                        continue

                    alts = p_var.ALT
                    if not alts:
                        continue

                    altlen = len(alts[0]) - len(p_var.REF)
                    altlen_all.append(altlen)

                    is_ins = altlen >= 50
                    is_del = altlen <= -50

                    # --- Panel Extraction ---
                    p_gt = p_var.gt_types
                    np.equal(p_gt, 1, out=p_is_het)
                    np.equal(p_gt, 0, out=p_is_hom_ref)
                    np.equal(p_gt, 3, out=p_is_hom_alt)

                    p_het_all += p_is_het
                    p_hom_ref_all += p_is_hom_ref
                    p_hom_alt_all += p_is_hom_alt
                    if is_ins: p_het_ins += p_is_het
                    elif is_del: p_het_del += p_is_het

                    p_n_het = p_var.num_het
                    p_n_hom_ref = p_var.num_hom_ref
                    p_n_hom_alt = p_var.num_hom_alt

                    panel_hwe['hom_ref'].append(p_n_hom_ref)
                    panel_hwe['het'].append(p_n_het)
                    panel_hwe['hom_alt'].append(p_n_hom_alt)

                    p_valid = (p_n_hom_ref + p_n_het + p_n_hom_alt) * 2
                    p_alt_count = p_n_het + 2 * p_n_hom_alt
                    panel_af_all.append(p_alt_count / p_valid if p_valid > 0 else 0.0)
                    panel_mean_alt_alleles_all.append(p_alt_count / num_p_samples)

                    # --- Target Extraction ---
                    c_gt = c_var.gt_types
                    np.equal(c_gt, 1, out=c_is_het)
                    np.equal(c_gt, 0, out=c_is_hom_ref)
                    np.equal(c_gt, 3, out=c_is_hom_alt)

                    c_het_all += c_is_het
                    c_hom_ref_all += c_is_hom_ref
                    c_hom_alt_all += c_is_hom_alt
                    if is_ins: c_het_ins += c_is_het
                    elif is_del: c_het_del += c_is_het

                    c_n_het = c_var.num_het
                    c_n_hom_ref = c_var.num_hom_ref
                    c_n_hom_alt = c_var.num_hom_alt

                    target_hwe['hom_ref'].append(c_n_hom_ref)
                    target_hwe['het'].append(c_n_het)
                    target_hwe['hom_alt'].append(c_n_hom_alt)

                    c_valid = (c_n_hom_ref + c_n_het + c_n_hom_alt) * 2
                    c_alt_count = c_n_het + 2 * c_n_hom_alt
                    target_af_all.append(c_alt_count / c_valid if c_valid > 0 else 0.0)
                    target_mean_alt_alleles_all.append(c_alt_count / num_c_samples)

                    variants_processed += 1
                    if variants_processed % 10000 == 0:
                        elapsed_secs = time.time() - start_time
                        elapsed_str = str(timedelta(seconds=int(elapsed_secs)))
                        print(f"Processed {variants_processed:,} matched records... [Elapsed: {elapsed_str}] [Location: {p_var.CHROM}:{p_var.POS}]", flush=True)

                panel_only += sum(1 for ak in p_amap if ak not in c_amap)
                pb = next(pit, None)
                cb = next(cit, None)
            elif pos_order(*p_key) < pos_order(*c_key):
                panel_only += len(p_block)
                pb = next(pit, None)
            else:
                target_only += len(c_block)
                cb = next(cit, None)

        # drain remainders (unmatched tails) for the report
        while pb is not None:
            panel_only += len(pb[1]); pb = next(pit, None)
        while cb is not None:
            target_only += len(cb[1]); cb = next(cit, None)

        elapsed_secs = time.time() - start_time
        elapsed_str = str(timedelta(seconds=int(elapsed_secs)))
        print(f"\nFinished: {variants_processed:,} matched sites in {elapsed_str}. "
              f"(panel-only: {panel_only:,}, target-only: {target_only:,})", flush=True)
        if variants_processed == 0:
            raise ValueError("No (CHROM,POS,REF,ALT)-matched sites between the panel and imputed VCFs "
                             "-- check the site-matched panel build (chrom naming / allele representation).")

        # Assign sample stats back into original dict structure for plotting
        sample_stats = {
            'panel': {
                'het_all': p_het_all, 'hom_ref_all': p_hom_ref_all, 'hom_alt_all': p_hom_alt_all,
                'het_ins': p_het_ins, 'het_del': p_het_del
            },
            'target': {
                'het_all': c_het_all, 'hom_ref_all': c_hom_ref_all, 'hom_alt_all': c_hom_alt_all,
                'het_ins': c_het_ins, 'het_del': c_het_del
            }
        }

        # Convert globals to numpy arrays for fast indexing
        panel_af = np.array(panel_af_all)
        target_af = np.array(target_af_all)
        altlen = np.array(altlen_all)
        is_sv_ins = altlen >= 50
        is_sv_del = altlen <= -50

        # 4. Pearson Correlations (TSV Output)
        print("Calculating Pearson correlations...", flush=True)
        pearson_records = []

        # Helper to safely calculate pearsonr against edge-case 0-variance bins
        def safe_pearson(x, y):
            if len(x) > 1 and np.std(x) > 0 and np.std(y) > 0:
                return pearsonr(x, y)[0]
            return np.nan

        r_all = safe_pearson(panel_af, target_af)
        pearson_records.append({"VARIANT_TYPE": "ALL", "PEARSON_R": r_all})

        if np.sum(is_sv_ins) > 1:
            r_ins = safe_pearson(panel_af[is_sv_ins], target_af[is_sv_ins])
            pearson_records.append({"VARIANT_TYPE": "SV_INS", "PEARSON_R": r_ins})

        if np.sum(is_sv_del) > 1:
            r_del = safe_pearson(panel_af[is_sv_del], target_af[is_sv_del])
            pearson_records.append({"VARIANT_TYPE": "SV_DEL", "PEARSON_R": r_del})

        is_sv = is_sv_ins | is_sv_del
        if np.sum(is_sv) > 1:
            r_sv = safe_pearson(panel_af[is_sv], target_af[is_sv])
            pearson_records.append({"VARIANT_TYPE": "SV", "PEARSON_R": r_sv})

        pd.DataFrame(pearson_records).to_csv(f'{output_prefix}.pearson.tsv', sep='\t', index=False)

        # 5. AF Hist2D Plots
        print("Generating AF Hist2D Plots...", flush=True)
        def plot_hist2d(p_af, c_af, title, outfile):
            if len(p_af) == 0: return
            plt.figure()
            plt.hist2d(p_af, c_af, bins=np.linspace(0, 1, 50), norm=matplotlib.colors.LogNorm())
            plt.title(title)
            plt.xlabel('AoU+HPRC2+HGSVC3 allele frequency')
            plt.ylabel('Target allele frequency')
            plt.gca().set_aspect('equal')
            cbar = plt.colorbar()
            cbar.set_label('Number of variants', rotation=270, labelpad=10)
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')
            plt.close()

        plot_hist2d(panel_af, target_af, 'All variants', f'{output_prefix}-AF-all')
        plot_hist2d(panel_af[is_sv_ins], target_af[is_sv_ins], 'SV-length insertions', f'{output_prefix}-AF-SV-ins')
        plot_hist2d(panel_af[is_sv_del], target_af[is_sv_del], 'SV-length deletions', f'{output_prefix}-AF-SV-del')

        # 6. Sample Metrics & Boxplots
        print("Generating population boxplots...", flush=True)
        population_df = pd.read_csv(population_tsv_path, sep='\t')

        pop_color_dict = {
            'ESN': '#ffcd00', 'GWD': '#ffb900', 'LWK': '#cc9933', 'MSL': '#e1b919', 'YRI': '#ffb933',
            'ACB': '#ff9900', 'ASW': '#cc6600', 'CLM': '#cc3333', 'MXL': '#e10033', 'PEL': '#ff0000',
            'PUR': '#cc3300', 'CDX': '#339900', 'CHB': '#adcd00', 'CHS': '#01ff00', 'JPT': '#008b00',
            'KHV': '#00cc33', 'CEU': '#0000ff', 'GBR': '#00c5cd', 'FIN': '#00ebff', 'IBS': '#6495ed',
            'TSI': '#00008b', 'BEB': '#8b008b', 'GIH': '#9400d3', 'ITU': '#b03060', 'PJL': '#e11289',
            'STU': '#ff00ff'
        }

        def make_results_df(samples, stats):
            df = pd.DataFrame()
            df['Sample'] = samples
            df['Heterozygous variants per sample'] = stats['het_all']
            df['Homozygous reference variants per sample'] = stats['hom_ref_all']
            df['Homozygous alternate variants per sample'] = stats['hom_alt_all']
            df['Heterozygous SV-length insertions per sample'] = stats['het_ins']
            df['Heterozygous SV-length deletions per sample'] = stats['het_del']

            df = pd.merge(df, population_df[['Sample name', 'Population code']], left_on='Sample', right_on='Sample name')
            return df.rename(columns={'Population code': 'Population'})

        panel_results_df = make_results_df(panel_samples, sample_stats['panel'])
        target_results_df = make_results_df(target_samples, sample_stats['target'])

        def plot_boxplot(df, x_col, title, outfile, xlim):
            # Only plot populations that actually have data. order=pop_color_dict.keys() would include
            # empty categories (e.g. when the auto one-population TSV uses code "ALL", which is not in
            # pop_color_dict) -> seaborn 0.13.2 draws zero boxes and leaves 'boxprops' unbound in
            # _configure_legend (UnboundLocalError). Restrict order/palette to present populations.
            present = list(pd.unique(df['Population'].dropna()))
            order = [p for p in pop_color_dict if p in present] + [p for p in present if p not in pop_color_dict]
            if not order:
                print(f"  [boxplot] no populations with data for '{x_col}'; skipping {outfile}", flush=True)
                return
            fallback = ['#7f7f7f', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
            palette = [pop_color_dict.get(p, fallback[i % len(fallback)]) for i, p in enumerate(order)]
            plt.figure(figsize=(4, 6))
            sns.boxplot(data=df[df['Population'].isin(order)], x=x_col, y='Population', hue='Population',
                        order=order, hue_order=order, palette=palette)
            plt.title(title)
            plt.xlim(xlim)
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')
            plt.close()

        # Guard each boxplot: a plotting failure must NOT discard the (expensive) streaming results --
        # pearson.tsv and the other PDFs are already written, so log and continue on any plot error.
        _boxplots = [
            (panel_results_df,  'Heterozygous variants per sample',            'HPRC2+HGSVC3 in AoU+HPRC2+HGSVC3', f'{output_prefix}-panel-het-all',    [0, 6E4]),
            (target_results_df, 'Heterozygous variants per sample',            'Target',                          f'{output_prefix}-target-het-all',   [0, 6E4]),
            (panel_results_df,  'Heterozygous SV-length insertions per sample','HPRC2+HGSVC3 in AoU+HPRC2+HGSVC3', f'{output_prefix}-panel-het-SV-ins', [0, 500]),
            (target_results_df, 'Heterozygous SV-length insertions per sample','Target',                          f'{output_prefix}-target-het-SV-ins',[0, 500]),
            (panel_results_df,  'Heterozygous SV-length deletions per sample', 'HPRC2+HGSVC3 in AoU+HPRC2+HGSVC3', f'{output_prefix}-panel-het-SV-del', [0, 500]),
            (target_results_df, 'Heterozygous SV-length deletions per sample', 'Target',                          f'{output_prefix}-target-het-SV-del',[0, 500]),
        ]
        for _df, _x, _t, _o, _xl in _boxplots:
            try:
                plot_boxplot(_df, _x, _t, _o, _xl)
            except Exception as _e:
                print(f"  [boxplot] WARNING: failed for '{_x}' ({_o}): {_e}", flush=True)
                plt.close('all')

        # 7. ALT Length Weighted Histogram
        print("Generating ALT Length Histogram...", flush=True)
        bins = list(np.linspace(-10000, -100, 397)) + [-75, -50, -25, -1, -0.1, 0.1, 1, 25, 50, 75] + list(np.linspace(100, 10000, 397))
        plt.figure()
        plt.hist(altlen, bins=bins, label='AoU+HPRC2+HGSVC3', log=True, histtype='step', alpha=0.5, weights=panel_mean_alt_alleles_all)
        plt.hist(altlen, bins=bins, label='Target', log=True, histtype='step', alpha=0.5, weights=target_mean_alt_alleles_all)
        plt.xscale('symlog', linthresh=100)
        plt.xticks([-1E4, -1E3] + list(np.linspace(-100, 100, 9)) + [1E3, 1E4],
                   labels=['$-10^4$', '$-10^3$', '$-10^2$'] + ['', '$-50$', '', '$0$', '', '$50$', ''] + ['$10^2$', '$10^3$', '$10^4$'])
        plt.ylabel('Number of ALT alleles per sample')
        plt.xlabel('ALT length - REF length (bp)')
        plt.legend()
        plt.savefig(f'{output_prefix}-alt-alleles-per-sample-hist.pdf', bbox_inches='tight')
        plt.close()

        # 8. De Finetti Plots
        print("Generating De Finetti Plots...", flush=True)
        ternary_to_cartesian = lambda a, b, c: (0.5 * (2 * b + c) / (a + b + c + 1E-10), 0.5 * np.sqrt(3) * c / (a + b + c + 1E-10))

        def calc_hwe_ternary(x, m=1, f=1):
            return np.array([1 - f, 0, 0]) + f * np.array([m * (1 - x)**2 + (1 - m) * (1 - x), m * x**2 + (1 - m) * x, 2 * m * x * (1 - x)])

        def make_de_finetti_ax():
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.set_xlim([-0.1, 1.1])
            ax.set_ylim([-0.1, np.sqrt(3) / 2 + 0.1])
            ax.set_aspect(np.sqrt(3) / 2)
            ax.axis('off')
            ep = 0.02
            ax.plot([-2 * ep / np.sqrt(3), 1 + 2 * ep / np.sqrt(3)], [-ep, -ep], lw=3, c='k')
            ax.plot([-2 * ep / np.sqrt(3), 0.5], [-ep, 0.5 * np.sqrt(3) + ep], lw=3, c='k')
            ax.plot([1 + 2 * ep / np.sqrt(3), 0.5], [-ep, 0.5 * np.sqrt(3) + ep], lw=3, c='k')
            ax.text(-4 * ep / np.sqrt(3), -4 * ep, 'HOM\nREF ', fontsize=16, ha='right')
            ax.text(0.5, 0.5 * np.sqrt(3) + 3 * ep, 'HET', fontsize=16, ha='center')
            ax.text(1 + 4 * ep / np.sqrt(3), -4 * ep, 'HOM\n ALT', fontsize=16, ha='left')
            return fig, ax

        def plot_de_finetti(hom_ref_arr, hom_alt_arr, het_arr, title, outfile, gridsize=70):
            if len(hom_ref_arr) == 0: return
            fig, ax = make_de_finetti_ax()

            # 1. Combine arrays into a 2D matrix
            counts_arr = np.column_stack((hom_ref_arr, hom_alt_arr, het_arr))

            # 2. Pre-aggregate identical variant counts
            unique_counts, point_weights = np.unique(counts_arr, axis=0, return_counts=True)

            # 3. Calculate X/Y coordinates ONLY for unique combinations
            x_ternary_v, y_ternary_v = ternary_to_cartesian(unique_counts[:,0], unique_counts[:,1], unique_counts[:,2])

            # 4. Pass the point_weights directly to hexbin using 'C'
            hb = ax.hexbin(
                x_ternary_v, 
                y_ternary_v, 
                C=point_weights, 
                reduce_C_function=np.sum, 
                gridsize=gridsize, 
                extent=[0, 1, 0, np.sqrt(3) / 2], 
                norm=matplotlib.colors.LogNorm()
            )

            cbar = plt.colorbar(hb, ax=ax, shrink=0.5)
            cbar.set_label('Number of variants', rotation=270, labelpad=10)

            x_values = np.linspace(0, 1, 50)
            cart_values = np.array([ternary_to_cartesian(*calc_hwe_ternary(x)) for x in x_values])
            ax.plot(cart_values[:, 0], cart_values[:, 1], c='C1', ls='solid', lw=3)

            ax.text(0.5, -0.3, title, fontsize=18, ha='center')
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')
            plt.close()

        # Slice HWE dictionaries for subsets
        def slice_hwe(hwe_dict, mask):
            mask = np.array(mask)
            return {k: np.array(v)[mask] for k, v in hwe_dict.items()}

        p_hwe_ins = slice_hwe(panel_hwe, is_sv_ins)
        p_hwe_del = slice_hwe(panel_hwe, is_sv_del)
        c_hwe_ins = slice_hwe(target_hwe, is_sv_ins)
        c_hwe_del = slice_hwe(target_hwe, is_sv_del)

        plot_de_finetti(panel_hwe['hom_ref'], panel_hwe['hom_alt'], panel_hwe['het'], 
                        'AoU+HPRC2+HGSVC3\nAll variants', f'{output_prefix}-panel-hwe-all')
        plot_de_finetti(p_hwe_ins['hom_ref'], p_hwe_ins['hom_alt'], p_hwe_ins['het'], 
                        'AoU+HPRC2+HGSVC3\nSV-length insertions', f'{output_prefix}-panel-hwe-SV-ins')
        plot_de_finetti(p_hwe_del['hom_ref'], p_hwe_del['hom_alt'], p_hwe_del['het'], 
                        'AoU+HPRC2+HGSVC3\nSV-length deletions', f'{output_prefix}-panel-hwe-SV-del')

        plot_de_finetti(target_hwe['hom_ref'], target_hwe['hom_alt'], target_hwe['het'], 
                        'Target\nAll variants', f'{output_prefix}-target-hwe-all', gridsize=80)
        plot_de_finetti(c_hwe_ins['hom_ref'], c_hwe_ins['hom_alt'], c_hwe_ins['het'], 
                        'Target\nSV-length insertions', f'{output_prefix}-target-hwe-SV-ins', gridsize=80)
        plot_de_finetti(c_hwe_del['hom_ref'], c_hwe_del['hom_alt'], c_hwe_del['het'], 
                        'Target\nSV-length deletions', f'{output_prefix}-target-hwe-SV-del', gridsize=80)

        print("All tasks completed successfully.", flush=True)
        EOF
    >>>

    output {
        File summarize_pearson_tsv = "~{output_prefix}.pearson.tsv"
        Array[File] summarize_plots = glob("*.pdf")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "python:3.11-slim"     # pinned cp311 -> matches the staged manylinux wheelhouse
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

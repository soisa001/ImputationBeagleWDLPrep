version 1.0

workflow GLIMPSE2Concordance {
    input {
        File panel_vcf          # split to biallelic
        File panel_vcf_idx
        File imputed_vcf        # split to biallelic
        File imputed_vcf_idx
        File trh_bed
        File trh_bed_idx
        String region
        String output_prefix
        File? concordance_binary   # pre-staged GLIMPSE2_concordance_static (perimeter blocks the github wget)
        File? metrics_wheelhouse   # cyvcf2+numpy wheelhouse for ExactGenotypeMetrics (perimeter blocks PyPI)
    }

    Array[String] trh_bins = ["outTRH", "inTRH"]
    Array[String] length_bins = ["SV_DEL", "DEL", "SNP", "INS", "SV_INS"]

    call AnnotateImputed { input:
        imputed_vcf = imputed_vcf,
        imputed_vcf_idx = imputed_vcf_idx,
        trh_bed = trh_bed,
        trh_bed_idx = trh_bed_idx,
        region = region,
        output_prefix = output_prefix
    }

    scatter (trh_bin in trh_bins) {
        scatter (length_bin in length_bins) {
            call FilterAndConcordance { input:
                annotated_bcf = AnnotateImputed.annotated_vcf,
                annotated_bcf_index = AnnotateImputed.annotated_vcf_idx,
                panel_vcf = panel_vcf,
                panel_vcf_idx = panel_vcf_idx,
                trh_bin = trh_bin,
                length_bin = length_bin,
                region = region,
                concordance_binary = concordance_binary,
                output_prefix = output_prefix + "." + trh_bin + "." + length_bin
            }
        }
    }

    Array[File] all_rsquare_grp_files = flatten(flatten(FilterAndConcordance.rsquare_grp_files))
    Array[File] all_rsquare_spl_files = flatten(flatten(FilterAndConcordance.rsquare_spl_files))
    Array[File] all_error_grp_files = flatten(flatten(FilterAndConcordance.error_grp_files))
    Array[File] all_error_spl_files = flatten(flatten(FilterAndConcordance.error_spl_files))
    Array[File] all_error_cal_files = flatten(flatten(FilterAndConcordance.error_cal_files))

    call PlotResults { input:
        rsquare_grp_files = all_rsquare_grp_files,
        error_spl_files = all_error_spl_files,
        output_prefix = output_prefix,
        panel_name = basename(panel_vcf),
        imputed_name = basename(imputed_vcf)
    }

    # EXACT per-sample precision/recall/F1 from the full 3x3 genotype confusion matrix, by streaming
    # truth (panel) vs imputed directly (the error.spl-derived per_sample_metrics can only give an
    # UPPER bound on precision because GLIMPSE groups mismatches by truth class, hiding het<->hom-alt swaps).
    call ExactGenotypeMetrics { input:
        panel_vcf = panel_vcf,
        panel_vcf_idx = panel_vcf_idx,
        annotated_imputed_vcf = AnnotateImputed.annotated_vcf,       # has GT:DS:GP + INFO/TRH
        annotated_imputed_vcf_idx = AnnotateImputed.annotated_vcf_idx,
        output_prefix = output_prefix,
        metrics_wheelhouse = metrics_wheelhouse
    }

    output {
        Array[File] concordance_plots = [PlotResults.r2_plot_inTRH, PlotResults.r2_plot_outTRH, PlotResults.nrd_plot_inTRH, PlotResults.nrd_plot_outTRH, PlotResults.recall_fpr_plot_inTRH, PlotResults.recall_fpr_plot_outTRH]
        File per_sample_metrics_tsv = PlotResults.per_sample_metrics_tsv
        File exact_per_sample_metrics_tsv = ExactGenotypeMetrics.exact_per_sample_metrics_tsv
        Array[Array[File]] concordance_results = [all_rsquare_grp_files, all_rsquare_spl_files, all_error_grp_files, all_error_spl_files, all_error_cal_files]
    }
}

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

task AnnotateImputed {
    input {
        File imputed_vcf
        File imputed_vcf_idx
        File trh_bed
        File trh_bed_idx
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 3 * ceil(size(imputed_vcf, "GiB") + size(trh_bed, "GiB"))

    command <<<
        set -euox pipefail
        
        bcftools annotate ~{imputed_vcf} \
            --threads $(nproc) \
            -r ~{region} \
            -a ~{trh_bed} -c CHROM,FROM,TO -m +TRH \
            --write-index=csi -Ob -o ~{output_prefix}.imputed.annotated.bcf
    >>>

    output {
        File annotated_vcf = "~{output_prefix}.imputed.annotated.bcf"
        File annotated_vcf_idx = "~{output_prefix}.imputed.annotated.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
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

task FilterAndConcordance {
    input {
        File annotated_bcf
        File annotated_bcf_index
        File panel_vcf
        File panel_vcf_idx
        String trh_bin
        String length_bin
        String region
        String output_prefix
        File? concordance_binary

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(annotated_bcf, "GiB") + size(panel_vcf, "GiB"))

    command <<<
        set -euox pipefail
        
        # Determine TRH Expression
        if [ "~{trh_bin}" == "outTRH" ]; then 
            TRH_EXP="TRH!=1"
        else 
            TRH_EXP="TRH==1"
        fi

        # Determine Length Expression
        if [ "~{length_bin}" == "SV_DEL" ]; then LEN_EXP="(STRLEN(ALT)-STRLEN(REF) <= -50)"; fi
        if [ "~{length_bin}" == "DEL" ]; then LEN_EXP="((-50 < STRLEN(ALT)-STRLEN(REF)) && (STRLEN(ALT)-STRLEN(REF) <= -1))"; fi
        if [ "~{length_bin}" == "SNP" ]; then LEN_EXP="((STRLEN(REF) == 1) && (STRLEN(ALT) == 1))"; fi
        if [ "~{length_bin}" == "INS" ]; then LEN_EXP="((STRLEN(REF) != 1) && (0 <= STRLEN(ALT)-STRLEN(REF)) && (STRLEN(ALT)-STRLEN(REF) < 50))"; fi
        if [ "~{length_bin}" == "SV_INS" ]; then LEN_EXP="(50 <= STRLEN(ALT)-STRLEN(REF))"; fi

        echo "Filtering with: $TRH_EXP & $LEN_EXP"

        bcftools view ~{annotated_bcf} \
            -i "$TRH_EXP & $LEN_EXP" \
            --threads $(nproc) \
            --write-index=csi -Ob -o ~{output_prefix}.bcf

        echo "~{region} ~{panel_vcf} ~{panel_vcf} ~{output_prefix}.bcf" > ~{output_prefix}.concordance-input.txt

        # Use a pre-staged binary if provided (the VPC-SC perimeter blocks the github wget); else fetch.
        if [ -n "~{concordance_binary}" ]; then
            cp ~{concordance_binary} GLIMPSE2_concordance_static
        else
            wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_concordance_static
        fi
        chmod +x GLIMPSE2_concordance_static

        ./GLIMPSE2_concordance_static \
            --min-tar-gp 0.0 0.9 \
            --gt-val \
            --use-alt-af \
            --out-r2-per-site \
            --bins 0.00001 0.00002 0.00005 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99 0.995 0.999 1.0 \
            --input ~{output_prefix}.concordance-input.txt \
            --threads $(nproc) \
            --output ~{output_prefix}.concordance-result
    >>>

    output {
        Array[File] rsquare_grp_files = glob("*.rsquare.grp.txt.gz")
        Array[File] rsquare_spl_files = glob("*.rsquare.spl.txt.gz")
        Array[File] error_grp_files = glob("*.error.grp.txt.gz")
        Array[File] error_spl_files = glob("*.error.spl.txt.gz")
        Array[File] error_cal_files = glob("*.error.cal.txt.gz")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
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

task PlotResults {
    input {
        Array[File] rsquare_grp_files
        Array[File] error_spl_files
        String output_prefix
        String panel_name
        String imputed_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + ceil(size(rsquare_grp_files, "GiB") + size(error_spl_files, "GiB"))

    command <<<
        set -euox pipefail

        cat << 'EOF' > plot_script.py
        import sys
        import os
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns

        rsquare_files = sys.argv[1].split(',')
        error_files = sys.argv[2].split(',')
        panel_name = sys.argv[3]
        imputed_name = sys.argv[4]

        # 1. Load Dosage R2 Data
        r2_vs_af_df_values = []
        for filepath in rsquare_files:
            if not filepath.strip(): continue
            filename = os.path.basename(filepath)
            
            parts = filename.split('_GPfilt_')
            prefix_parts = parts[0].split('.')
            length_bin = prefix_parts[-2]
            trh_bin = prefix_parts[-3]
            min_tar_gp = float(parts[1].split('.rsquare')[0])

            df = pd.read_csv(filepath, sep=' ', comment='#', names=['AF_BIN_INDEX', 'AF_BIN_COUNT', 'AF_BIN_MEAN', 'R2_GT', 'R2_DS'])
            for _, row in df.iterrows():
                if row['AF_BIN_COUNT'] > 0:
                    r2_vs_af_df_values.append([trh_bin, length_bin, min_tar_gp, row['AF_BIN_COUNT'], row['AF_BIN_MEAN'], row['R2_DS']])

        r2_vs_af_df = pd.DataFrame(r2_vs_af_df_values, columns=['TRH_BIN', 'LENGTH_BIN', 'MIN_TAR_GP', 'AF_BIN_COUNT', 'AF_BIN_MEAN', 'R2_DS'])

        # 2. Load Error/Concordance Rate Data
        sample_df_values = []
        for filepath in error_files:
            if not filepath.strip(): continue
            filename = os.path.basename(filepath)
            
            parts = filename.split('_GPfilt_')
            prefix_parts = parts[0].split('.')
            length_bin = prefix_parts[-2]
            trh_bin = prefix_parts[-3]
            min_tar_gp = float(parts[1].split('.error')[0])

            cols = 'GCsV id sample_name #val_gt_RR #val_gt_RA #val_gt_AA #filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared'.split(' ')
            df = pd.read_csv(filepath, sep=' ', comment='#', names=cols)
            
            for _, row in df.iterrows():
                if row['GCsV'] != 'GCsV': continue
                sample_df_values.append([trh_bin, length_bin, min_tar_gp, row['sample_name'], float(row['non_reference_discordanc_rate_percent']), float(row['imputed_ds_rsquared'])])

        sample_df = pd.DataFrame(sample_df_values, columns=['TRH_BIN', 'LENGTH_BIN', 'MIN_TAR_GP', 'sample_name', 'non_reference_discordanc_rate_percent', 'imputed_ds_rsquared'])

        # 3. Calculate metrics for title
        num_samples = sample_df['sample_name'].nunique()
        title_metadata = f"Panel: {panel_name}\nTarget: {imputed_name}\nEvaluated samples: {num_samples}"

        # 4. Generate Dosage R2 Plots
        for trh_bin in ['outTRH', 'inTRH']:
            fig, ax = plt.subplots(1, 5, figsize=(10, 2))
            for i, length_bin in enumerate(['SV_DEL', 'DEL', 'SNP', 'INS', 'SV_INS']):
                ax2 = ax[i].twinx()
                for min_tar_gp in [0.0, 0.9]:
                    label = 'unfiltered' if min_tar_gp == 0.0 else f'GP > {min_tar_gp}'
                    x = (r2_vs_af_df['TRH_BIN'] == trh_bin) & (r2_vs_af_df['LENGTH_BIN'] == length_bin) & (r2_vs_af_df['MIN_TAR_GP'] == min_tar_gp)
                    
                    if not r2_vs_af_df[x].empty:
                        ax[i].plot(r2_vs_af_df[x]['AF_BIN_MEAN'], r2_vs_af_df[x]['R2_DS'], label=label, 
                                   ls={'unfiltered': 'solid', 'GP > 0.9': 'dotted'}[label], color='C0')
                        ax2.plot(r2_vs_af_df[x]['AF_BIN_MEAN'], r2_vs_af_df[x]['AF_BIN_COUNT'], label=label, 
                                 ls={'unfiltered': 'solid', 'GP > 0.9': 'dotted'}[label], color='C1')

                ax[i].set_xscale('log')
                ax[i].set_xlim([1E-4, 1])
                ax[i].set_ylim([0, 1])
                ax2.set_yscale('log')
                ax2.set_ylim([10**2, 10**9])
                ax2.set_yticks([10**j for j in range(2, 10)])

                if i == 0:
                    ax[i].set_ylabel('$r^2_{dosage}$', fontsize=14, color='C0')
                    ax2.set_yticklabels([])
                elif i == 4:
                    ax2.set_ylabel('number of variants', fontsize=14, color='C1', rotation=270, va='bottom')
                    ax[i].set_yticklabels([])
                else:
                    ax[i].set_yticklabels([])
                    ax2.set_yticklabels([])

                if i == 2:
                    trh_tag = {'outTRH': 'non-TR/homopolymer', 'inTRH': 'TR/homopolymer'}[trh_bin]
                    ax[i].set_title(f"{title_metadata}\n{trh_tag}\n", fontsize=10)
                    ax[i].set_xlabel(f'panel allele frequency\n\n{length_bin}\nALT length - REF length (bp)', fontsize=12)
                    ax[i].legend(loc='lower right', fontsize=8)
                else:
                    length_bin_label = {'SV_DEL': '(-inf, -50]', 'DEL': '(-50, -1]', 'SNP': 'SNP', 'INS': '[0, 50)', 'SV_INS': '[50, inf)'}[length_bin]
                    ax[i].set_xlabel(f'\n\n{length_bin_label}', fontsize=12)

            plt.savefig(f'~{output_prefix}.{trh_bin}.r2.png', bbox_inches='tight')
            plt.close()

        # 5. Generate Error Plots
        for trh_bin in ['outTRH', 'inTRH']:
            fig, ax = plt.subplots(1, 1, figsize=(6, 3))
            trh_tag = {'outTRH': 'non-TR/homopolymer', 'inTRH': 'TR/homopolymer'}[trh_bin]
            ax.set_title(f"{title_metadata}\n{trh_tag}", fontsize=10)
            
            plt_df_values = []
            for length_bin in ['SV_DEL', 'DEL', 'SNP', 'INS', 'SV_INS']:
                length_bin_label = {'SV_DEL': '(-inf, -50]', 'DEL': '(-50, -1]', 'SNP': 'SNP', 'INS': '[0, 50)', 'SV_INS': '[50, inf)'}[length_bin]
                for min_tar_gp in [0.0, 0.9]:
                    min_tar_gp_label = 'unfiltered' if min_tar_gp == 0.0 else f'GP > {min_tar_gp}'
                    x = (sample_df['TRH_BIN'] == trh_bin) & (sample_df['LENGTH_BIN'] == length_bin) & (sample_df['MIN_TAR_GP'] == min_tar_gp)
                    bin_df = sample_df[x]
                    
                    for s in range(bin_df.shape[0]):
                        plt_df_values.append([length_bin_label, min_tar_gp_label, 1 - 0.01 * bin_df['non_reference_discordanc_rate_percent'].values[s]])

            if plt_df_values:
                plt_df = pd.DataFrame(plt_df_values, columns=['LENGTH_BIN_TEXT', 'MIN_TAR_GP_TEXT', 'non_reference_discordanc_rate_percent'])
                sns.boxplot(data=plt_df, x='LENGTH_BIN_TEXT', y='non_reference_discordanc_rate_percent', hue='MIN_TAR_GP_TEXT', ax=ax)
                ax.set_xlabel('ALT length - REF length (bp)', fontsize=8)
                ax.set_ylabel('non-reference concordance rate', fontsize=8)
                ax.set_ylim([0, 1.01])
                plt.legend(loc='lower center', fontsize=8)
                plt.tight_layout()
                plt.savefig(f'~{output_prefix}.{trh_bin}.nrd.png', bbox_inches='tight')
            plt.close()

        # 6. Per-sample precision/recall + FP/FN from the genotype confusion counts.
        #    GLIMPSE reports matches/mismatches grouped by the TRUTH genotype class, so:
        #      - false positives (truth hom-ref called as carrying an alt) are EXACT (mismatch_RR);
        #      - non-ref recall / overall concordance / NRC are EXACT;
        #      - a clean false-NEGATIVE rate is NOT recoverable: het/hom-alt mismatches lump
        #        "called hom-ref" (a true miss) with "called the other non-ref genotype" (a
        #        genotype swap, still a detected variant). We report carrier_error_rate as the
        #        UPPER bound on the miss/FN rate, and PPV only vs hom-ref FPs (swaps excluded).
        #        The full 3x3 (and thus exact precision/FN) would need the genotypes, not these
        #        per-truth-class aggregates -- see README.
        print("Computing per-sample concordance metrics...", flush=True)
        # field order per GLIMPSE2 concordance .error.spl (call_set_writing.cpp), tag GCsV = SNP+indel
        count_cols = ('tag idx sample_name val_RR val_RA val_AA filtered_gp '
                      'mRR mRA mAA xRR xRA xAA '
                      'fp_rate_pct het_mm_rate_pct homalt_mm_rate_pct '
                      'nrd_pct best_gt_r2 ds_r2').split()
        metric_rows = []
        for filepath in error_files:
            if not filepath.strip(): continue
            filename = os.path.basename(filepath)
            parts = filename.split('_GPfilt_')
            prefix_parts = parts[0].split('.')
            length_bin = prefix_parts[-2]
            trh_bin = prefix_parts[-3]
            min_tar_gp = float(parts[1].split('.error')[0])

            df = pd.read_csv(filepath, sep=' ', comment='#', names=count_cols)
            df = df[df['tag'] == 'GCsV']
            for _, r in df.iterrows():
                val_RR, val_RA, val_AA = float(r['val_RR']), float(r['val_RA']), float(r['val_AA'])
                mRR, mRA, mAA = float(r['mRR']), float(r['mRA']), float(r['mAA'])
                xRR, xRA, xAA = float(r['xRR']), float(r['xRA']), float(r['xAA'])

                nonref_true  = val_RA + val_AA           # condition-positive: truth carries an alt
                nonref_match = mRA + mAA                 # exact-GT recovery of those
                nonref_mm    = xRA + xAA                 # truth-carrier miscalls (miss OR swap)
                total_gt     = val_RR + val_RA + val_AA

                nan = float('nan')
                recall        = nonref_match / nonref_true if nonref_true > 0 else nan   # exact-GT sensitivity
                het_recall    = mRA / val_RA if val_RA > 0 else nan
                homalt_recall = mAA / val_AA if val_AA > 0 else nan
                fp_rate       = xRR / val_RR if val_RR > 0 else nan                      # truth hom-ref -> alt (exact)
                nrc           = nonref_match / (nonref_match + xRR + nonref_mm) if (nonref_match + xRR + nonref_mm) > 0 else nan
                gt_conc       = (mRR + mRA + mAA) / total_gt if total_gt > 0 else nan    # overall exact-match (incl hom-ref)
                carrier_err   = nonref_mm / nonref_true if nonref_true > 0 else nan      # UPPER bound on miss/FN rate
                ppv_vs_homref = nonref_match / (nonref_match + xRR) if (nonref_match + xRR) > 0 else nan

                metric_rows.append([
                    trh_bin, length_bin, min_tar_gp, r['sample_name'],
                    int(val_RR), int(val_RA), int(val_AA),
                    int(mRR), int(mRA), int(mAA), int(xRR), int(xRA), int(xAA),
                    recall, het_recall, homalt_recall, fp_rate, nrc, gt_conc, carrier_err, ppv_vs_homref,
                ])

        metrics_df = pd.DataFrame(metric_rows, columns=[
            'TRH_BIN', 'LENGTH_BIN', 'MIN_TAR_GP', 'sample_name',
            'n_truth_RR', 'n_truth_RA', 'n_truth_AA',
            'match_RR', 'match_RA', 'match_AA', 'mismatch_RR_FP', 'mismatch_RA', 'mismatch_AA',
            'nonref_recall', 'het_recall', 'homalt_recall', 'false_pos_rate', 'nonref_concordance',
            'overall_gt_concordance', 'carrier_error_rate_upperbound', 'ppv_vs_homref',
        ])
        metrics_df.to_csv(f'~{output_prefix}.per_sample_metrics.tsv', sep='\t', index=False)
        print(f"Wrote {metrics_df.shape[0]} per-sample metric rows -> ~{output_prefix}.per_sample_metrics.tsv", flush=True)

        # 7. Recall + FP-rate boxplots (mirror the NRC plot layout: length bins x GP filter)
        length_label = {'SV_DEL': '(-inf, -50]', 'DEL': '(-50, -1]', 'SNP': 'SNP', 'INS': '[0, 50)', 'SV_INS': '[50, inf)'}
        panels = [('nonref_recall', 'non-ref recall (sensitivity)', [0, 1.01]),
                  ('false_pos_rate', 'false-positive rate (hom-ref -> alt)', None)]
        for trh_bin in ['outTRH', 'inTRH']:
            sub = metrics_df[metrics_df['TRH_BIN'] == trh_bin]
            trh_tag = {'outTRH': 'non-TR/homopolymer', 'inTRH': 'TR/homopolymer'}[trh_bin]
            fig, axes = plt.subplots(1, 2, figsize=(11, 3))
            for ax, (col, ylab, ylim) in zip(axes, panels):
                rows = []
                for length_bin in ['SV_DEL', 'DEL', 'SNP', 'INS', 'SV_INS']:
                    for min_tar_gp in [0.0, 0.9]:
                        gp_lab = 'unfiltered' if min_tar_gp == 0.0 else f'GP > {min_tar_gp}'
                        m = (sub['LENGTH_BIN'] == length_bin) & (sub['MIN_TAR_GP'] == min_tar_gp)
                        for v in sub[m][col].values:
                            rows.append([length_label[length_bin], gp_lab, v])
                if rows:
                    bx = pd.DataFrame(rows, columns=['LENGTH_BIN_TEXT', 'GP', 'val'])
                    sns.boxplot(data=bx, x='LENGTH_BIN_TEXT', y='val', hue='GP', ax=ax)
                    ax.legend(loc='best', fontsize=7)
                ax.set_xlabel('ALT length - REF length (bp)', fontsize=8)
                ax.set_ylabel(ylab, fontsize=8)
                if ylim is not None:
                    ax.set_ylim(ylim)
            fig.suptitle(f"{title_metadata}\n{trh_tag}", fontsize=9)
            plt.tight_layout()
            plt.savefig(f'~{output_prefix}.{trh_bin}.recall_fpr.png', bbox_inches='tight')
            plt.close()
        EOF

        python3 plot_script.py "~{sep=',' rsquare_grp_files}" "~{sep=',' error_spl_files}" "~{panel_name}" "~{imputed_name}"
    >>>

    output {
        File r2_plot_inTRH = "~{output_prefix}.inTRH.r2.png"
        File r2_plot_outTRH = "~{output_prefix}.outTRH.r2.png"
        File nrd_plot_inTRH = "~{output_prefix}.inTRH.nrd.png"
        File nrd_plot_outTRH = "~{output_prefix}.outTRH.nrd.png"
        File recall_fpr_plot_inTRH = "~{output_prefix}.inTRH.recall_fpr.png"
        File recall_fpr_plot_outTRH = "~{output_prefix}.outTRH.recall_fpr.png"
        File per_sample_metrics_tsv = "~{output_prefix}.per_sample_metrics.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "jupyter/scipy-notebook:latest"
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

task ExactGenotypeMetrics {
    input {
        File panel_vcf                    # truth (full panel, multi-sample), biallelic
        File panel_vcf_idx
        File annotated_imputed_vcf        # imputed + INFO/TRH (GT:DS:GP), biallelic
        File annotated_imputed_vcf_idx
        String output_prefix
        File? metrics_wheelhouse          # cyvcf2 + numpy wheelhouse (offline; perimeter blocks PyPI)

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 20 + 3 * ceil(size(panel_vcf, "GiB") + size(annotated_imputed_vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # cyvcf2 (bundles htslib) + numpy, offline from the staged wheelhouse if provided.
        if [ -n "~{metrics_wheelhouse}" ]; then
            mkdir -p wheelhouse && tar -xzf ~{metrics_wheelhouse} -C wheelhouse
            pip install --no-index --find-links=wheelhouse cyvcf2 numpy
        else
            pip install cyvcf2 numpy
        fi

        cat << 'EOF' > exact_metrics.py
        import sys, time
        import numpy as np
        import cyvcf2

        panel_path, imputed_path, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]

        GP_THRESH = [0.0, 0.9]
        LEN_NAMES = ['SV_DEL', 'DEL', 'SNP', 'INS', 'SV_INS']
        TRH_NAMES = ['outTRH', 'inTRH']
        CLS_LUT = np.array([0, 1, -1, 2], dtype=np.int64)   # cyvcf2 gt_types 0,1,2,3 -> RR,RA,missing,AA

        imp = cyvcf2.VCF(imputed_path)
        imp_samples = list(imp.samples)
        n = len(imp_samples)
        imp_index = {s: i for i, s in enumerate(imp_samples)}

        # subset the (huge) panel to the imputed samples at the htslib level -> only those GTs decoded
        panel = cyvcf2.VCF(panel_path, samples=imp_samples)
        panel_samples = list(panel.samples)                 # subset, in panel's file order
        psamp = np.array([imp_index[s] for s in panel_samples], dtype=np.int64)  # panel col -> canonical (imputed) idx
        sys.stderr.write(f"samples: imputed={n}, panel-matched={len(panel_samples)}\n"); sys.stderr.flush()
        if len(panel_samples) == 0:
            sys.exit("ERROR: no overlapping samples between imputed and panel")

        # M[sample, trh(2), length(5), gp(2), truth(3), call(3)]
        M = np.zeros((n, 2, 5, 2, 3, 3), dtype=np.int64)

        def length_idx(ref, alt):
            d = len(alt) - len(ref)                          # replicate FilterAndConcordance's LEN_EXP
            if d <= -50: return 0                            # SV_DEL
            if -50 < d <= -1: return 1                       # DEL
            if len(ref) == 1 and len(alt) == 1: return 2     # SNP
            if len(ref) != 1 and 0 <= d < 50: return 3       # INS
            if d >= 50: return 4                             # SV_INS
            return -1

        def pos_blocks(vcf):
            block, key = [], None
            for v in vcf:
                k = (v.CHROM, v.POS)
                if key is not None and k != key:
                    yield key, block; block = []
                block.append(v); key = k
            if block: yield key, block

        def akey(v):
            return (v.REF, v.ALT[0] if v.ALT else '.')

        try:
            chrom_rank = {c: i for i, c in enumerate(imp.seqnames)}
        except Exception:
            chrom_rank = {}
        def order(k):
            return (chrom_rank.get(k[0], 1 << 30), k[1])

        pit, cit = pos_blocks(panel), pos_blocks(imp)
        pb, cb = next(pit, None), next(cit, None)
        matched, t0 = 0, time.time()
        while pb is not None and cb is not None:
            pk, pblk = pb; ck, cblk = cb
            if pk == ck:
                pmap = {}
                for v in pblk: pmap.setdefault(akey(v), v)
                for cv in cblk:
                    pv = pmap.get(akey(cv))
                    if pv is None: continue
                    li = length_idx(cv.REF, cv.ALT[0] if cv.ALT else '')
                    if li < 0: continue
                    trh = cv.INFO.get('TRH')
                    ti = 1 if (trh is not None and int(trh) == 1) else 0

                    truth_canon = np.full(n, 2, dtype=np.int64)   # default UNKNOWN -> excluded
                    truth_canon[psamp] = pv.gt_types              # panel GTs mapped to canonical order
                    t_cls = CLS_LUT[truth_canon]
                    q_cls = CLS_LUT[cv.gt_types]                  # imputed already canonical
                    gp = cv.format('GP')
                    maxgp = gp.max(axis=1) if gp is not None else np.ones(n, dtype=np.float32)

                    vi = np.where((t_cls >= 0) & (q_cls >= 0))[0]
                    if vi.size:
                        M[:, ti, li, 0] += np.bincount(vi * 9 + t_cls[vi] * 3 + q_cls[vi],
                                                       minlength=n * 9).reshape(n, 3, 3)
                        vf = vi[maxgp[vi] >= 0.9]
                        if vf.size:
                            M[:, ti, li, 1] += np.bincount(vf * 9 + t_cls[vf] * 3 + q_cls[vf],
                                                           minlength=n * 9).reshape(n, 3, 3)
                    matched += 1
                    if matched % 200000 == 0:
                        sys.stderr.write(f"  matched {matched:,} sites ... {int(time.time()-t0)}s\n"); sys.stderr.flush()
                pb, cb = next(pit, None), next(cit, None)
            elif order(pk) < order(ck):
                pb = next(pit, None)
            else:
                cb = next(cit, None)
        sys.stderr.write(f"done: {matched:,} matched sites in {int(time.time()-t0)}s\n")

        def safe(a, b):
            return a / b if b > 0 else float('nan')

        hdr = ['TRH_BIN', 'LENGTH_BIN', 'MIN_TAR_GP', 'sample_name',
               'n_RR_RR', 'n_RR_RA', 'n_RR_AA', 'n_RA_RR', 'n_RA_RA', 'n_RA_AA', 'n_AA_RR', 'n_AA_RA', 'n_AA_AA',
               'TP_carrier', 'FP', 'FN_carrier', 'TN', 'precision', 'recall', 'F1',
               'TP_exact', 'precision_exact', 'recall_exact', 'F1_exact', 'overall_gt_conc', 'nonref_conc']

        def row(tname, lname, g, sname, m):
            # m is a 3x3 truth(RR,RA,AA) x call(RR,RA,AA) confusion matrix
            n00, n01, n02 = int(m[0, 0]), int(m[0, 1]), int(m[0, 2])
            n10, n11, n12 = int(m[1, 0]), int(m[1, 1]), int(m[1, 2])
            n20, n21, n22 = int(m[2, 0]), int(m[2, 1]), int(m[2, 2])
            tot = n00 + n01 + n02 + n10 + n11 + n12 + n20 + n21 + n22
            truth_nonref = n10 + n11 + n12 + n20 + n21 + n22
            called_nonref = n01 + n11 + n21 + n02 + n12 + n22
            # carrier level (a het<->hom-alt swap still detects the carrier -> TP)
            TPc = n11 + n12 + n21 + n22
            FPc = n01 + n02
            FNc = n10 + n20
            TN = n00
            prec_c, rec_c = safe(TPc, TPc + FPc), safe(TPc, TPc + FNc)
            f1_c = safe(2 * TPc, 2 * TPc + FPc + FNc)
            # genotype-exact (a swap is both a wrong call and a missed correct call)
            TPg = n11 + n22
            prec_g, rec_g = safe(TPg, called_nonref), safe(TPg, truth_nonref)
            f1_g = safe(2 * TPg, called_nonref + truth_nonref)
            return [tname, lname, g, sname,
                    n00, n01, n02, n10, n11, n12, n20, n21, n22,
                    TPc, FPc, FNc, TN, prec_c, rec_c, f1_c,
                    TPg, prec_g, rec_g, f1_g,
                    safe(n00 + n11 + n22, tot), safe(n11 + n22, tot - n00)]

        import csv
        with open(f'{out_prefix}.exact_per_sample_metrics.tsv', 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t'); w.writerow(hdr)
            # top-line: micro-average over all samples (pool confusion counts, then compute rates).
            # sample_name='ALL' rows come first so the aggregate is easy to grep out.
            Msum = M.sum(axis=0)   # [trh, length, gp, 3, 3]
            for ti, tname in enumerate(TRH_NAMES):
                for li, lname in enumerate(LEN_NAMES):
                    for gi, g in enumerate(GP_THRESH):
                        w.writerow(row(tname, lname, g, 'ALL', Msum[ti, li, gi]))
            for s in range(n):
                for ti, tname in enumerate(TRH_NAMES):
                    for li, lname in enumerate(LEN_NAMES):
                        for gi, g in enumerate(GP_THRESH):
                            w.writerow(row(tname, lname, g, imp_samples[s], M[s, ti, li, gi]))
        sys.stderr.write("wrote exact_per_sample_metrics.tsv (with sample_name='ALL' micro-average rows)\n")
        EOF

        python exact_metrics.py ~{panel_vcf} ~{annotated_imputed_vcf} ~{output_prefix}
    >>>

    output {
        File exact_per_sample_metrics_tsv = "~{output_prefix}.exact_per_sample_metrics.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "python:3.11-slim"
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

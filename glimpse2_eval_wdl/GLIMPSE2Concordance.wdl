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

    output {
        Array[File] concordance_plots = [PlotResults.r2_plot_inTRH, PlotResults.r2_plot_outTRH, PlotResults.nrd_plot_inTRH, PlotResults.nrd_plot_outTRH]
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

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_concordance_static
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
        EOF

        python3 plot_script.py "~{sep=',' rsquare_grp_files}" "~{sep=',' error_spl_files}" "~{panel_name}" "~{imputed_name}"
    >>>

    output {
        File r2_plot_inTRH = "~{output_prefix}.inTRH.r2.png"
        File r2_plot_outTRH = "~{output_prefix}.outTRH.r2.png"
        File nrd_plot_inTRH = "~{output_prefix}.inTRH.nrd.png"
        File nrd_plot_outTRH = "~{output_prefix}.outTRH.nrd.png"
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

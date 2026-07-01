version 1.0

# vcfdist phasing + precision/recall evaluation of the holdout198 imputation.
#
# Adapted from the Broad starter for this pipeline's reality:
#   - ONE multi-sample eval VCF (the imputed output) and ONE multi-sample truth VCF
#     (the 198's full-panel genotypes) -- both carry the same 198 sample ids -- instead of
#     a per-sample array of truth VCFs (so the truth is localized once per subset, not 198x).
#   - every workflow input is a SCALAR (File/String/Boolean) or a newline-delimited FILE
#     (sample names, bed list, labels) read with read_lines(), so the launcher can pass them
#     through the proven `wb` batch-input-CSV / column-mapping mechanism (which is scalar-only).
#   - SummarizeEvaluations installs pandas OFFLINE from a pre-staged pip wheelhouse (the VPC-SC
#     perimeter blocks PyPI from tasks) and also aggregates the PHASING metrics (switch/flip),
#     not just precision/recall.
#
# vcfdist is representation-aware (it aligns query vs truth), so the imputed output and the
# truth do NOT need a matched allele representation -- popped or un-popped both work as eval.

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

struct VcfdistOutputs {
    File summary_vcf
    File precision_recall_summary_tsv
    File precision_recall_tsv
    File query_tsv
    File truth_tsv
    File phasing_summary_tsv
    File switchflips_tsv
    File superclusters_tsv
    File phase_blocks_tsv
}

workflow VcfdistEvaluation {
    input {
        File eval_vcf                       # multi-sample imputed output (popped or un-popped); phased
        File eval_vcf_idx
        File truth_vcf                      # multi-sample truth (the 198 full-panel genotypes); phased
        File? truth_vcf_idx

        File sample_names_file              # eval sample ids, one per line (the 198)
        File? truth_sample_names_file       # truth ids if they differ (e.g. "syndip"); default = sample_names_file

        File? confident_regions_bed         # optional shared confident-regions BED (applied to both sides)

        String region                       # e.g. "chr1"
        File reference_fasta
        File reference_fasta_fai

        Boolean do_naively_phase = false    # convert / to | in the EVAL (only) if it is unphased

        File vcfdist_bed_list               # analysis/stratification BEDs (gs:// or local), one per line
        File labels_list                    # label per stratification, parallel to vcfdist_bed_list
        # NOT String? : Cromwell makes each scatter body a sub-workflow, for which a referenced
        # optional-without-default String becomes a *required* sub-workflow input and fails to look up
        # when omitted. A defaulted String is always defined, so it propagates into the scatter cleanly.
        String vcfdist_extra_args = ""
        Int vcfdist_mem_gb = 64              # vcfdist RAM per shard; 16 OOMs on the dense popped rep

        File? pip_wheelhouse                 # pandas wheelhouse for SummarizeEvaluations (offline pip)
    }

    Array[String] eval_sample_names  = read_lines(sample_names_file)
    Array[String] truth_sample_names = read_lines(select_first([truth_sample_names_file, sample_names_file]))
    Array[File]   vcfdist_bed_files  = read_lines(vcfdist_bed_list)
    Array[String] labels_per_stratification = read_lines(labels_list)

    scatter (i in range(length(eval_sample_names))) {
        call SubsetSampleFromVcf as SubsetSampleFromVcfEval { input:
            vcf = eval_vcf,
            vcf_idx = eval_vcf_idx,
            original_sample_name = eval_sample_names[i],
            sample_name = eval_sample_names[i],
            region = region,
            bed_file = confident_regions_bed,
            reference_fasta_fai = reference_fasta_fai,
            do_naively_phase = do_naively_phase
        }

        call SubsetSampleFromVcf as SubsetSampleFromVcfTruth { input:
            vcf = truth_vcf,
            vcf_idx = truth_vcf_idx,
            original_sample_name = truth_sample_names[i],
            sample_name = eval_sample_names[i],     # rename truth sample to match eval
            region = region,
            bed_file = confident_regions_bed,
            reference_fasta_fai = reference_fasta_fai,
            do_naively_phase = false
        }
    }

    scatter (j in range(length(vcfdist_bed_files))) {
        scatter (i in range(length(eval_sample_names))) {
            call Vcfdist { input:
                sample_name = eval_sample_names[i],
                eval_vcf = SubsetSampleFromVcfEval.single_sample_vcf[i],
                truth_vcf = SubsetSampleFromVcfTruth.single_sample_vcf[i],
                bed_file = vcfdist_bed_files[j],
                reference_fasta = reference_fasta,
                extra_args = vcfdist_extra_args,
                mem_gb = vcfdist_mem_gb
            }
        }
    }

    call SummarizeEvaluations { input:
        labels_per_vcf = labels_per_stratification,
        vcfdist_outputs_per_vcf_and_sample = Vcfdist.outputs,
        pip_wheelhouse = pip_wheelhouse
    }

    output {
        # stratification x sample
        Array[Array[VcfdistOutputs]] vcfdist_summary = Vcfdist.outputs
        File evaluation_summary_tsv = SummarizeEvaluations.evaluation_summary_tsv
    }
}

task SubsetSampleFromVcf {
    input {
        File vcf
        File? vcf_idx
        String original_sample_name
        String sample_name
        String region
        File? bed_file
        File reference_fasta_fai
        Boolean do_naively_phase = false

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size([vcf, bed_file], "GiB")) + 10

    command <<<
        set -euxo pipefail

        # localize the (multi-sample) VCF + its index into the cwd under a matching basename so bcftools
        # finds the index; create one if none was provided (avoids the fragile Terra "##idx##" hint).
        B="$(basename ~{vcf})"
        ln -s ~{vcf} "./${B}"
        IDX="~{if defined(vcf_idx) then vcf_idx else ''}"
        if [ -n "${IDX}" ]; then ln -s "${IDX}" "./$(basename "${IDX}")"; else bcftools index "./${B}"; fi

        # subset to one sample over the region (must combine -r region with -T bed to intersect properly)
        bcftools view "./${B}" \
            -s ~{original_sample_name} \
            -r ~{region} \
            ~{"-T " + bed_file} \
            ~{if do_naively_phase
                then "-Ou | bcftools +setGT -Oz -o " + sample_name + ".subset.vcf.gz -- -t a -n p"
                else "-Oz -o " + sample_name + ".subset.vcf.gz"}

        # rename the (single) sample to the eval id and refresh contig lines from the .fai
        echo ~{sample_name} > sample_name.txt
        bcftools reheader ~{sample_name}.subset.vcf.gz \
            -s sample_name.txt \
            --fai ~{reference_fasta_fai} \
            -o ~{sample_name}.subset.reheadered.vcf.gz
        bcftools index -t ~{sample_name}.subset.reheadered.vcf.gz
    >>>

    output {
        File single_sample_vcf = "~{sample_name}.subset.reheadered.vcf.gz"
        File single_sample_vcf_tbi = "~{sample_name}.subset.reheadered.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "LOCAL",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
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

task Vcfdist {
    input {
        String sample_name
        File eval_vcf
        File truth_vcf
        File bed_file
        File reference_fasta
        String extra_args = ""       # defaulted (not String?) so the scatter sub-workflow always resolves it
        Int verbosity = 1
        Int cpu = 1
        Int mem_gb = 64              # 16 OOM-killed during wavefront clustering on the dense popped rep;
                                     # pair with a supercluster cap (-s) in extra_args to bound it

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size([eval_vcf, truth_vcf, reference_fasta], "GiB")) + 10

    command <<<
        set -euxo pipefail

        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{reference_fasta} \
            -b ~{bed_file} \
            -v ~{verbosity} \
            ~{extra_args}

        for tsv in $(ls *.tsv); do mv "$tsv" ~{sample_name}."$tsv"; done
        mv summary.vcf ~{sample_name}.summary.vcf
    >>>

    output {
        VcfdistOutputs outputs = object {
            summary_vcf: "~{sample_name}.summary.vcf",
            precision_recall_summary_tsv: "~{sample_name}.precision-recall-summary.tsv",
            precision_recall_tsv: "~{sample_name}.precision-recall.tsv",
            query_tsv: "~{sample_name}.query.tsv",
            truth_tsv: "~{sample_name}.truth.tsv",
            phasing_summary_tsv: "~{sample_name}.phasing-summary.tsv",
            switchflips_tsv: "~{sample_name}.switchflips.tsv",
            superclusters_tsv: "~{sample_name}.superclusters.tsv",
            phase_blocks_tsv: "~{sample_name}.phase-blocks.tsv"
        }
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpu,
        mem_gb:             mem_gb,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "timd1/vcfdist:v2.6.4"
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

task SummarizeEvaluations {
    input {
        Array[String] labels_per_vcf
        Array[Array[VcfdistOutputs]] vcfdist_outputs_per_vcf_and_sample
        File? pip_wheelhouse

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        # pandas: the perimeter blocks PyPI from the task, so install OFFLINE from a pre-staged wheelhouse
        # (cp311 manylinux, matching python:3.11-slim) when provided; else fall back to an online install.
        if [ -n "~{pip_wheelhouse}" ]; then
            mkdir -p wheelhouse && tar -xzf ~{pip_wheelhouse} -C wheelhouse
            pip install --no-index --find-links=wheelhouse pandas
        else
            pip install pandas
        fi

        cat << 'EOF' > summarize.py
        import argparse
        import json
        import pandas as pd

        # vcfdist precision-recall-summary.tsv: 2-level index (VAR_TYPE, FILTER); we read VAR_TYPE x 'NONE'.
        PR_COLS = ['TRUTH_TP', 'QUERY_TP', 'TRUTH_FN', 'QUERY_FP', 'PREC', 'RECALL', 'F1_SCORE']
        VAR_TYPES = ['SNP', 'INDEL', 'SV']

        def pr_metrics_for_sample(path):
            out = {}
            try:
                df = pd.read_csv(path, sep='\t', index_col=[0, 1])
            except Exception as e:
                print(f"  [pr] could not read {path}: {e}", flush=True)
                return out
            for vt in VAR_TYPES:
                try:
                    row = df.loc[vt, 'NONE']
                    cols = [c for c in PR_COLS if c in row.index]
                    out.update(row[cols].add_prefix(f'{vt}_').to_dict())
                except Exception:
                    pass
            return out

        def phasing_metrics_for_sample(outs):
            # Average all NUMERIC columns of phasing-summary.tsv (robust to vcfdist column naming),
            # plus count switch/flip events and phase blocks. Prefix PHASE_.
            out = {}
            try:
                ph = pd.read_csv(outs['phasing_summary_tsv'], sep='\t')
                num = ph.select_dtypes(include='number')
                if not num.empty:
                    out.update(num.mean(axis=0).add_prefix('PHASE_').to_dict())
            except Exception as e:
                print(f"  [phase] could not read phasing-summary: {e}", flush=True)
            for key, label in [('switchflips_tsv', 'PHASE_N_SWITCHFLIP_EVENTS'),
                               ('phase_blocks_tsv', 'PHASE_N_PHASE_BLOCKS')]:
                try:
                    out[label] = float(max(0, sum(1 for _ in open(outs[key])) - 1))  # rows minus header
                except Exception:
                    pass
            return out

        def summarize_over_samples(vcfdist_outputs_per_sample):
            per_sample = {}
            for s, outs in enumerate(vcfdist_outputs_per_sample):
                rec = {}
                rec.update(pr_metrics_for_sample(outs['precision_recall_summary_tsv']))
                rec.update(phasing_metrics_for_sample(outs))
                per_sample[s] = rec
            return pd.DataFrame.from_dict(per_sample, orient='index').mean(axis=0)

        def summarize(labels_per_vcf_txt, vcfdist_outputs_per_vcf_and_sample_json):
            with open(labels_per_vcf_txt) as f:
                labels = f.read().splitlines()
            with open(vcfdist_outputs_per_vcf_and_sample_json) as f:
                vcfdist_outputs_per_vcf_and_sample = json.load(f)

            summary = {}
            for i, label in enumerate(labels):
                summary[label] = {'NUM_VCFDIST_SAMPLES': len(vcfdist_outputs_per_vcf_and_sample[i])}
                summary[label].update(summarize_over_samples(vcfdist_outputs_per_vcf_and_sample[i]).to_dict())
            pd.DataFrame.from_dict(summary, orient='index').to_csv(
                'evaluation_summary.tsv', sep='\t', float_format='%.4f')

        def main():
            p = argparse.ArgumentParser()
            p.add_argument('--labels_per_vcf_txt', type=str)
            p.add_argument('--vcfdist_outputs_per_vcf_and_sample_json', type=str)
            a = p.parse_args()
            summarize(a.labels_per_vcf_txt, a.vcfdist_outputs_per_vcf_and_sample_json)

        if __name__ == '__main__':
            main()
        EOF

        python3 summarize.py \
            --labels_per_vcf_txt ~{write_lines(labels_per_vcf)} \
            --vcfdist_outputs_per_vcf_and_sample_json ~{write_json(vcfdist_outputs_per_vcf_and_sample)}
    >>>

    output {
        File evaluation_summary_tsv = "evaluation_summary.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        disk_type:          "HDD",
        preemptible_tries:  2,
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

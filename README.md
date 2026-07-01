# aou-lrwgs-imputebeagle-holdout

Leave-N-out (chr1) benchmark of **ACAF short-read → AoU lrWGS Phase-2 panel**
imputation with Beagle, on Verily Workbench / All of Us controlled tier.

198 held-out AoU samples' **real ACAF short-read SNVs** are projected onto the
panel's bubble-allele representation, imputed against the **leaveout** panel
(bref3), the imputed bubbles are popped back to constituent variants, and the
result is scored against those samples' **full-panel long-read genotypes** — i.e.
empirical ACAF→panel imputation accuracy, not an optimistic panel-derived target.

---

## Pipeline at a glance

There is **one main script** (imputation) and **two independent eval scripts**
that run afterwards. Each script stages its inputs, registers its WDL with `wb`,
and launches it on the workspace Cromwell.

| stage | script | what it does | key outputs |
|-------|--------|--------------|-------------|
| **① Impute** (main) | `prep_imputebeagle_holdout198.sh` | build leaveout bref3 + ACAF target scaffold → project to the panel's bubble rep → run the Beagle WDL (`ImputeBeagleWithPop`) → pop bubbles to atomic variants | `…imputed.vcf.gz` (un-popped / bubble.split) and `…imputed.popped.vcf.gz` (popped / atomic) |
| **② Eval — accuracy** | `prep_eval_holdout198.sh` | `GLIMPSE2Concordance` (dosage r² & NRD + per-sample recall/FP) and `GLIMPSE2Summarize` (AF / HWE / ALT-length) vs the panel truth | r²/NRD/recall PNGs, `pearson.tsv`, `per_sample_metrics.tsv` |
| **② Eval — phasing** | `prep_vcfdist_holdout198.sh` | `vcfdist` per-sample vs the 198 truth: **phasing** (switch/flip) + precision/recall (representation-aware) | per-sample TSVs + `evaluation_summary.tsv` |
| **(opt) Validate** | `debug_check_outputs.sh` | coherence + representation checks on the popped / un-popped outputs | pass/fail report |

The two evals are independent — run either, both, or neither, in any order, once
imputation has produced its outputs.

---

## Layout

```
prep_imputebeagle_holdout198.sh     # ① MAIN: build ref+target, project, register + run the imputation WDL
prep_eval_holdout198.sh             # ② EVAL (accuracy): GLIMPSE2Concordance + GLIMPSE2Summarize
prep_vcfdist_holdout198.sh          # ② EVAL (phasing):  vcfdist per-sample
debug_check_outputs.sh              # (opt) validate the imputed outputs (representation coherence)
run_summarize_local_once.sh         # (opt) one-off LOCAL re-run of GLIMPSE2Summarize (skips Cromwell)

project_to_panel_rep.rs             # FAST minimal-rep projection (rustc, pure-std; the default projector)
project_to_panel_rep.py             # reference projection (Python); the .rs is byte-identical to this
test_projection_equivalence.sh      # asserts the Rust projector == the Python (adversarial + fuzz)
test_bubble_representation_matching.py   # optional: quantify exact vs min-rep scaffold reach

imputationbeagle_wdl_flat/          # ① imputation workflow bundle (import closure = these 4)
    ImputeBeagleWithPop.wdl         #   the workflow
    ImputationBeagleStructs.wdl
    ImputationTasks.wdl
    ImputationBeagleTasks.wdl       #   contains the Rust pop task (PopAndMarginalizeCollisions)

glimpse2_eval_wdl/                  # ② accuracy eval workflows (single-file each)
    GLIMPSE2Concordance.wdl
    GLIMPSE2Summarize.wdl
    GLIMPSE2_concordance_static     #   vendored binary (perimeter blocks the WDL's github wget)

vcfdist_eval_wdl/                   # ② phasing eval workflow (single-file)
    VcfdistEvaluation.wdl

pop_glimpse2_rust/                  # pop engine source (built once -> static-musl binary)
    pop-glimpse2.rs
    Cargo.toml
minrep_equiv/                       # optional: prove the projection's matching == the Rust minrep
scripts/build_pop_binary.sh         # standalone build of the static-musl pop-glimpse2 binary
```

Run every script **from the repo root** — they auto-detect `./imputationbeagle_wdl_flat`,
`./glimpse2_eval_wdl`, `./vcfdist_eval_wdl`, and `./project_to_panel_rep.*`.

---

## Prerequisites (on the workbench VM)

- CLIs: `bcftools`, `tabix`, `bgzip`, `gsutil`, `gcloud`, `plink2`, `java`, `python3`, `wb`, `rustc` (pop/projector builds).
- The workspace bucket resolves via `wb resource resolve --name rw-migration-aou-rw-f178dfde`.
- The VM needs PyPI egress (it builds the eval pip wheelhouses); the Cromwell **tasks** do not (perimeter-blocked — see *Perimeter / offline deps* below).

## Controlled-access data (supplied on the VM, never committed)

- **Truth** — `truth.aou.chr1.bcf` (+`.csi`): the 198's full-panel genotypes. Used to derive the 198 id list and as the eval truth. Default `gs://cloned-rw-migration-aou-rw-f178dfde-wb-sharp-papaya-7463/vcf/truth.aou.chr1.bcf`; override `TRUTH_AOU=`.
- **ACAF target** — pulled from the chromosome-sharded ACAF **PGEN** (`<contig>.pgen`+`.pvar[.zst]`+`.psam`) with `plink2 --keep <198> --export vcf`; override `AOU_PGEN_GS=`. One `≤~100 GB` file/chromosome is downloaded to local SSD, used, then removed (a single chromosome fits in ~500 GB local disk).
- **Panels** — leaveout / sites-only / id-split / popped panels are pulled from `gs://rw-long-reads-transfer-2026-06-17/…` automatically.

`.gitignore` blocks all VCF/BCF/BED/binary artifacts.

---

## How to run

### ① Impute (main)

```bash
# Idempotent: reruns skip already-staged artifacts. The prep BUILDS the static-musl
# pop-glimpse2 binary itself (the in-task cargo build can't reach crates.io in the perimeter);
# reuses ~/pop-build. Pass POP_BINARY_LOCAL=<path> to skip the build.
nohup bash prep_imputebeagle_holdout198.sh > holdout198_prep.log 2>&1 &
tail -f holdout198_prep.log
#   watch for: ">> built pop-glimpse2 locally: …", "target samples: 198 …",
#              "[project] in(alleles)=… exact=… recovered=… dropped=…"

wb workflow job list | head        # monitor to completion
```
On success, two outputs land under `…/imputebeagle-${RUN_ID}-run/` (default `RUN_ID=holdout198_v2`):
- `aou_holdout198_chr1.imputed.vcf.gz` — un-popped (bubble.split)
- `aou_holdout198_chr1.imputed.popped.vcf.gz` — popped (constituent/atomic variants)

### ② Evaluate (after imputation completes)

Two independent evals. Both auto-discover the imputed output under the run namespace,
and both take `POPPED=true` to score the popped output instead of the un-popped one.

**Accuracy — GLIMPSE2 (dosage r² / NRD / AF / HWE):**
```bash
nohup bash prep_eval_holdout198.sh > holdout198_eval.log 2>&1 &          # un-popped vs bubble.split panel
POPPED=true bash prep_eval_holdout198.sh                                  # popped vs popped panel
```

**Phasing + precision/recall — vcfdist:**
```bash
MAX_SAMPLES=5 bash prep_vcfdist_holdout198.sh                             # cheap test on 5 samples first
MAX_SAMPLES=0 bash prep_vcfdist_holdout198.sh                             # all 198
POPPED=true MAX_SAMPLES=0 bash prep_vcfdist_holdout198.sh                 # score the popped output
```

### (optional) Validate the imputed outputs
```bash
USER_PROJECT=<billing-project> bash debug_check_outputs.sh               # requester-pays panel reads
```

### After editing any eval WDL
Re-run that eval once with `FORCE_WF_REGISTER=true` so the edited WDL is re-uploaded and
re-registered (e.g. `POPPED=true FORCE_WF_REGISTER=true bash prep_eval_holdout198.sh`).

---

## Single chromosome only (chr1)

The whole pipeline is **per-contig** and fully supports having only chr1 done — nothing
assumes a genome-wide (chr1–22) VCF:

- **Imputation** scatters over `CONTIGS` (default `chr1`); add contigs to do more.
- **vcfdist** loads the genome-wide reference FASTA but only evaluates the `-b` region
  (a whole-`chr1` BED built from the `.fai`); it processes whatever contigs the query/truth
  actually contain (your run logged *"All contig checks passed"* on chr1-only input).
- **Both eval summaries aggregate across *samples*, not chromosomes** — vcfdist's
  `SummarizeEvaluations` and GLIMPSE2's plots bin by AF / variant-length / TR-context, never
  by a fixed chromosome set, and every per-sample read is guarded so a class with no data on
  chr1 is skipped, not fatal. vcfdist produces per-sample **TSVs** (no plots); GLIMPSE2
  produces the plots.

If you later feed genome-wide VCFs to an *external* plotting/reporting tool that assumes all
chromosomes, restrict it to chr1 — but none of the scripts here need that.

---

## Key env overrides

### Imputation (`prep_imputebeagle_holdout198.sh`)
| var | default | meaning |
|-----|---------|---------|
| `RUN_ID` | `holdout198_v2` | run namespace. Bump for a **fresh run**: changes the gs:// staging dir, the registered workflow name, and the output path together. Use the same `RUN_ID` for the eval scripts. |
| `CONTIGS` | `chr1` | contig(s) to impute |
| `TRUTH_AOU` | `…/vcf/truth.aou.chr1.bcf` | 198 full-panel genotypes (id source + eval truth) |
| `AOU_PGEN_GS` | acaf_threshold `pgen/` dir | ACAF PGEN source (`<contig>.pgen`+`.pvar[.zst]`+`.psam`) |
| `TARGET_FILTER` | `PASS,.` | site FILTER kept on the plink2 export (`PASS,.` keeps PASS + unfiltered) |
| `ENABLE_POP` | `true` | also emit the popped output |
| `POP_BINARY_LOCAL` | (unset → prep builds it) | prebuilt pop-glimpse2 binary; else built from `pop_glimpse2_rust/` |
| `PROJECT_BIN` | (built from `.rs`) | prebuilt Rust projector; else `project_to_panel_rep.rs` is built with `rustc -O` |
| `THREADS` / `TARGET_CONCURRENT` | `nproc` / `true` | plink2+norm threads / overlap the ACAF pull with the panel prep |

### Evals (`prep_eval_holdout198.sh`, `prep_vcfdist_holdout198.sh`)
| var | default | meaning |
|-----|---------|---------|
| `RUN_ID` / `REGION` | `holdout198_v2` / `chr1` | must match the imputation run |
| `POPPED` | `false` | score the popped output (`.popped`/`-popped`-namespaced) vs the popped panel |
| `FORCE_WF_REGISTER` | `false` | re-upload + re-register the WDL (use after editing it) |
| `IMPUTED_VCF` | (auto-discover) | override the imputed output path |
| `MAX_SAMPLES` (vcfdist) | `5` | samples to score; `0` = all 198 |
| `VCFDIST_MEM_GB` (vcfdist) | `16` | RAM per shard |
| `VCFDIST_EXTRA_ARGS` (vcfdist) | `--cluster size 100 --max-supercluster-size 10000` | bound clustering memory on the dense popped rep |
| `REF_FASTA` (vcfdist) | public GRCh38 | reference; a chr1-subset cuts the per-shard localization ~12× |
| `SUMMARIZE_WHEELHOUSE` / `PIP_WHEELHOUSE` | (built on the VM) | pre-staged pip wheelhouse (if the VM lacks PyPI egress) |

---

## Implementation notes

### Imputation & the ACAF scaffold
- The ACAF target is **hard-call GTs only** (ACAF carries no PL), so it's a GT scaffold, not PL-based.
  Multiallelics are split (`bcftools norm -m -any`, and the projection splits internally too) with GT
  recoding; only SNV-equivalent alleles are scaffolded — indels/MNVs are imputed, not scaffolded.
- The projection runs the Rust `project_to_panel_rep.rs` by default (`rustc -O`, cached in
  `PROJECT_BUILD_DIR`) — a pure-std port **byte-identical** to `project_to_panel_rep.py`
  (`bash test_projection_equivalence.sh`), ~10–50× faster. Output is `bcftools sort`ed (temp on data disk via `-T`).
- The pop task (`PopAndMarginalizeCollisions`) no longer runs `bcftools annotate` (memory heavy):
  `pop-glimpse2` reads the bubble ID straight from the sites-only panel via a streaming, memory-bounded
  `(CHROM,POS,REF,ALT)` join (ported from upstream `pop-glimpse2-joint-opt.rs`, minus mimalloc). Verified
  identical to the old annotate path.

### GLIMPSE2 accuracy eval
- Concordance emits, alongside the NRC boxplot, **per-sample precision/recall + FP metrics** from GLIMPSE's
  `.error.spl` confusion counts: `<prefix>.per_sample_metrics.tsv` (raw `n_truth_*`/`match_*`/`mismatch_*`
  counts + `nonref_recall`, `false_pos_rate`, `nonref_concordance` (=1−NRD), `overall_gt_concordance`, …)
  and `<prefix>.{in,out}TRH.recall_fpr.png`. **Caveat:** GLIMPSE groups mismatches by the *truth* class, so
  hom-ref→alt false positives are exact but a clean false-negative rate is not recoverable (a truth-het
  mismatch could be a miss or a het↔hom-alt swap). `carrier_error_rate_upperbound` bounds the miss rate;
  `ppv_vs_homref` is precision vs hom-ref FPs only. An exact 3×3 would need re-streaming the genotypes.
- `GLIMPSE2Summarize` streams panel vs imputed with a **`(CHROM,POS,REF,ALT)` merge-join** (not a positional
  `zip()`), so a different allele set/order at a multiallelic position (e.g. `chr1:10626 A>AG` vs `A>AGGCGCAG`)
  no longer raises "VCFs out of sync" — it pairs by exact allele per position and skips/reports the rest.

### vcfdist phasing eval
- [vcfdist](https://github.com/TimD1/vcfdist) is **representation-aware** (it realigns query vs truth), so the
  imputed output and the truth need no matched allele rep — popped or un-popped both work. It scores per sample:
  the eval and truth multi-sample VCFs are split per sample (`SubsetSampleFromVcf`), then `vcfdist` runs once per
  sample; `SummarizeEvaluations` aggregates the per-sample TSVs into `evaluation_summary.tsv`.
- Adapted from the Broad starter: a single multi-sample truth VCF (the 198 share one truth) instead of a
  per-sample array, and every workflow input is a scalar or a `read_lines` file so the `wb` batch-input
  mechanism applies.
- The dense **popped** representation can make vcfdist's superclusters blow up memory/runtime, so the prep
  defaults to `--cluster size 100 --max-supercluster-size 10000` at `VCFDIST_MEM_GB=16`. If a shard still OOMs,
  raise `VCFDIST_MEM_GB` and/or tighten `--max-supercluster-size` (floor = largest-variant+2 = 5002).

### Perimeter / offline deps (VPC-SC blocks anaconda/PyPI/github from Cromwell tasks)
- **GLIMPSE2_concordance** binary: vendored (`glimpse2_eval_wdl/GLIMPSE2_concordance_static`), staged by the prep,
  passed as `concordance_binary` (the WDL's github `wget` is blocked).
- **GLIMPSE2Summarize / vcfdist SummarizeEvaluations** python deps: the prep builds a pip wheelhouse (cp311
  manylinux) on the VM, stages it, and the task installs offline with `pip install --no-index --find-links`
  (docker `python:3.11-slim`). The staged wheelhouse path is content-addressed by the package set, so changing
  the pins rebuilds it rather than reusing a stale tarball. matplotlib is pinned `<3.10`, and the Summarize
  boxplot only orders populations that have data (seaborn 0.13.2 raises `boxprops` on empty groups). Override
  with `SUMMARIZE_WHEELHOUSE=gs://…` / `PIP_WHEELHOUSE=gs://…` if the VM lacks PyPI egress.

### Misc
- Sample-id namespace must match between panel/truth and ACAF; the prep errors if 0 of the 198 are found.
- `run_summarize_local_once.sh` re-runs GLIMPSE2Summarize **locally** (venv with the correct pins) when you just
  need the plots without another Cromwell round-trip; it re-streams (~1h20m) but hands you the outputs directly.

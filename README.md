# aou-lrwgs-imputebeagle-holdout

Leave-N-out (chr1) benchmark of **ACAF short-read → AoU lrWGS Phase-2 panel**
imputation with Beagle, on Verily Workbench / All of Us controlled tier.

The 198 held-out AoU samples' **real ACAF short-read SNVs** are projected onto the
panel's bubble-allele representation, imputed against the **leaveout** panel
(bref3), the imputed bubbles are popped to constituent variants, and the result is
scored against those samples' **full-panel** long-read genotypes — i.e. empirical
ACAF→panel imputation accuracy, not an optimistic panel-derived target.

## Layout

```
prep_imputebeagle_holdout198.sh     # STEP 1: reference (bref3) + ACAF target (split+project) + register + run
prep_eval_holdout198.sh             # STEP 3: GLIMPSE2Concordance + GLIMPSE2Summarize on the imputed output
project_to_panel_rep.rs             # FAST minimal-rep projection (rustc, pure-std; the default projector)
project_to_panel_rep.py             # reference projection (Python); the .rs is byte-identical to this
test_projection_equivalence.sh      # asserts the Rust projector == the Python (adversarial + fuzz)
test_bubble_representation_matching.py   # optional: quantify exact vs min-rep scaffold reach

imputationbeagle_wdl_flat/          # main workflow bundle (import closure = these 4)
    ImputeBeagleWithPop.wdl
    ImputationBeagleStructs.wdl
    ImputationTasks.wdl
    ImputationBeagleTasks.wdl       # contains the Rust pop task (PopAndMarginalizeCollisions)

glimpse2_eval_wdl/                  # eval workflows (single-file each)
    GLIMPSE2Concordance.wdl
    GLIMPSE2Summarize.wdl

vcfdist_eval_wdl/                   # vcfdist phasing + precision/recall eval (single-file)
    VcfdistEvaluation.wdl
prep_vcfdist_holdout198.sh          # STEP 4: register + run vcfdist per-sample vs the 198 truth

pop_glimpse2_rust/                  # pop engine source (build once -> static musl binary)
    pop-glimpse2.rs
    Cargo.toml

minrep_equiv/                       # optional: prove the Python projection's matching == the Rust minrep
    minrep_cli.rs
    check_minrep_equivalence.py

scripts/
    build_pop_binary.sh             # builds the static-musl pop-glimpse2 binary (no sudo/docker)
```

**Run the prep scripts from the repo root** — they auto-detect
`./imputationbeagle_wdl_flat`, `./project_to_panel_rep.py`, and
`./glimpse2_eval_wdl`.

## Prerequisites on the workbench app

- `bcftools`, `tabix`, `bgzip`, `gsutil`, `gcloud`, `plink2`, `java`, `python3`, `wb` CLI, `rustc` (for the pop binary).
- The clone workspace bucket resolves via `wb resource resolve --name rw-migration-aou-rw-f178dfde`.

## Data NOT in this repo (controlled access — supply on the VM)

- `truth.aou.chr1.bcf` (+ `.csi`): the 198's full-panel genotypes (your `check_holdout_panel.sh` Step A).
  Used to derive the 198 id list and as the eval truth. Defaults to `gs://cloned-rw-migration-aou-rw-f178dfde-wb-sharp-papaya-7463/vcf/truth.aou.chr1.bcf`; override with `TRUTH_AOU=`.
- ACAF genotypes. Pulled from the **chromosome-sharded ACAF PGEN** with `plink2` (`--keep` the
  198 → export VCF), from `gs://vwb-aou-datasets-controlled/.../acaf_threshold/pgen` (override
  `AOU_PGEN_GS`); expects `<contig>.pgen` + `.pvar[.zst]` + `.psam`. plink2 subsets the 198 from
  the packed binary genotypes (no 245k-column text parse) and PGEN preserves REF/ALT/indels for
  the min-rep projection. One `<=~100GB` file per chromosome is downloaded to local SSD, used,
  then removed (so a single chromosome fits in ~500GB local disk).
- Leaveout / sites-only / id-split panels are pulled from `gs://rw-long-reads-transfer-2026-06-17/…` automatically.

`.gitignore` blocks all VCF/BCF/BED/binary artifacts so controlled data is never committed.

## Run order

```bash
# (optional) prove the projection's matching logic == the Rust minrep
#   (the checker auto-builds minrep_cli from minrep_cli.rs if rustc is present)
( cd minrep_equiv && python3 check_minrep_equivalence.py )

# 1) prep + launch imputation (idempotent; reruns skip staged artifacts).
#    The prep BUILDS the static-musl pop-glimpse2 binary itself from pop_glimpse2_rust/
#    (the in-task cargo build can't reach crates.io inside the perimeter). Reuses an
#    existing build at ~/pop-build. Pass POP_BINARY_LOCAL=<path> to skip the build.
nohup bash prep_imputebeagle_holdout198.sh > holdout198_prep.log 2>&1 &
tail -f holdout198_prep.log
#   watch for: ">> building pop-glimpse2 (static ...)" then ">> built pop-glimpse2 locally: ...",
#              "target samples: 198 ...",
#              "[project] in(alleles)=... exact=... recovered=... dropped=..."

# 2) monitor the workflow to completion
wb workflow job list | head
#   on success two outputs land under .../imputebeagle-${RUN_ID}-run/ (default RUN_ID=holdout198_v2):
#     aou_holdout198_chr1.imputed.vcf.gz         (un-popped, bubble.split)
#     aou_holdout198_chr1.imputed.popped.vcf.gz  (popped, constituent variants)

# 3) eval the un-popped output vs full-panel truth (no required inputs)
nohup bash prep_eval_holdout198.sh > holdout198_eval.log 2>&1 &
tail -f holdout198_eval.log

# 3b) popped eval: scores the popped (constituent) output vs the popped FULL panel
#     (panel_popped_vcf). Staging + outputs are namespaced (.popped / -popped) so it
#     does not clobber the un-popped eval.
# POPPED=true bash prep_eval_holdout198.sh
```

The first imputation run incurs a one-time rust toolchain install (user-local rustup) + pop
build; later runs reuse `~/pop-build`. `scripts/build_pop_binary.sh` is an optional standalone
pre-build (e.g. to build once before launching); the prep does the same automatically.

## Key env overrides

| var | default | meaning |
|-----|---------|---------|
| `RUN_ID` | `holdout198_v2` | run namespace. Bump it for a **fresh run**: changes the gs:// staging dir (`imp_<RUN_ID>/`), the registered workflow name (`ImputeBeagleWithPop_<RUN_ID>`), and the run output path (`imputebeagle-<RUN_ID>-run`) together, so nothing reuses previously staged/cached artifacts. Set the same `RUN_ID` for `prep_eval_holdout198.sh`. |
| `POP_BINARY_LOCAL` | (unset → prep builds it) | prebuilt pop-glimpse2 binary; if unset the prep builds a static-musl one from `pop_glimpse2_rust/` |
| `POP_BUILD_LOCAL` | `true` | build the pop binary locally in the prep (set `false` to fall back to the in-task source build) |
| `POP_BUILD_DIR` | `~/pop-build` | where the pop binary is built/cached (idempotent reuse) |
| `TRUTH_AOU` | `gs://cloned-rw-migration-aou-rw-f178dfde-wb-sharp-papaya-7463/vcf/truth.aou.chr1.bcf` | 198 full-panel genotypes (id source + eval truth) |
| `AOU_SAMPLES` | (derive from TRUTH_AOU) | explicit 198 id list (local/gs) |
| `AOU_PGEN_GS` | acaf_threshold `pgen/` gs:// dir | ACAF PGEN source dir; expects `<contig>.pgen` + `.pvar[.zst]` + `.psam` |
| `TARGET_FILTER` | `PASS,.` | site FILTER kept (bcftools `-f`) on the plink2 export; `PASS,.` keeps PASS + unfiltered (`.`), drops LowQual/ExcessHet/… Set `PASS` for strict, empty to disable |
| `THREADS` | `nproc` | plink2 export + bcftools norm threads |
| `PROJECT_BIN` | (built from `.rs`) | prebuilt Rust projector; if unset, `project_to_panel_rep.rs` is built with `rustc -O` and used (byte-identical to the Python, ~10-50x faster). No rustc/.rs → falls back to the Python |
| `PROJECT_BUILD_DIR` | `~/project-build` | where the Rust projector is built/cached (idempotent reuse) |
| `PROJECT_LOCAL` | (auto-detect) | path to project_to_panel_rep.py (the fallback projector) |
| `CONTIGS` | `chr1` | contig(s) |
| `ENABLE_POP` | `true` | emit popped output |
| `TARGET_CONCURRENT` | `true` | overlap the ACAF PGEN pull+export with the panel prep + bref3 build |

## Notes

- The ACAF target comes from the chromosome-sharded ACAF **PGEN** via `plink2 --keep <198> --export vcf`
  (hard-call GTs only — ACAF carries no PL, so this is a GT scaffold, not PL-based). Multiallelic sites
  are split (`bcftools norm -m -any`, and the projection splits internally too) with GT recoding; only
  SNV-equivalent alleles (true + padded) are scaffolded — indels/MNVs are imputed, not scaffolded.
- The ACAF PGEN pull+export runs concurrently with the panel prep + bref3 build by default
  (`TARGET_CONCURRENT=true`). Set `TARGET_CONCURRENT=false` for the old serial order.
- The projection runs the Rust `project_to_panel_rep.rs` by default (built once with `rustc -O`,
  cached in `PROJECT_BUILD_DIR`); it is a pure-std port that is **byte-identical** to
  `project_to_panel_rep.py` (verify with `bash test_projection_equivalence.sh`) but ~10-50x faster.
  The output is still piped through `bcftools sort` (temp on the data disk via `-T`).
- The bubble-pop task (`PopAndMarginalizeCollisions`) no longer runs `bcftools annotate` (memory
  heavy). `pop-glimpse2` now reads the bubble ID straight from the sites-only panel (its 2nd
  positional arg) via a streaming, memory-bounded (CHROM,POS,REF,ALT) join, ported from upstream
  `pop-glimpse2-joint-opt.rs` (minus mimalloc, to keep the static-musl build pure-Rust). Verified to
  produce identical popped output to the old annotate-based path.
- The perimeter blocks the Concordance WDL's github `wget` of `GLIMPSE2_concordance_static`, so the
  binary is vendored (`glimpse2_eval_wdl/GLIMPSE2_concordance_static`); the eval prep stages it to the
  bucket and passes it as the `concordance_binary` input. After updating the eval WDL, re-run the eval
  with `FORCE_WF_REGISTER=true` once so the new WDL is re-uploaded + re-registered.
- The Concordance WDL also emits **per-sample precision/recall + false-positive metrics** alongside the
  NRC boxplot, derived from GLIMPSE's per-sample `.error.spl` confusion counts (no extra passes over the
  data, so the binning matches the NRC plot exactly):
    - `<prefix>.per_sample_metrics.tsv` — one row per (sample × TRH bin × length bin × GP filter) with the
      raw confusion counts (`n_truth_RR/RA/AA`, `match_*`, `mismatch_*`) and derived rates:
      `nonref_recall` (exact-GT sensitivity for alt carriers), `het_recall`, `homalt_recall`,
      `false_pos_rate` (truth hom-ref called as carrying an alt — "called but not in truth"),
      `nonref_concordance` (= 1 − NRD, the boxplot metric), `overall_gt_concordance` (exact match incl.
      hom-ref — the inflated "raw GT match rate"), `carrier_error_rate_upperbound`, and `ppv_vs_homref`.
    - `<prefix>.{in,out}TRH.recall_fpr.png` — recall + FP-rate boxplots in the same layout as the NRC plot.
  **Caveat (by design of GLIMPSE's output):** GLIMPSE reports matches/mismatches grouped by the *truth*
  genotype class, so hom-ref→alt **false positives are exact**, but a clean **false-negative** rate is
  *not* recoverable — a truth-het scored as a mismatch could be "called hom-ref" (a true miss) or "called
  hom-alt" (a genotype swap, still a detected variant), and the aggregates don't separate them.
  `carrier_error_rate_upperbound` (= 1 − `nonref_recall`) bounds the miss/FN rate from above, and
  `ppv_vs_homref` is precision counting only hom-ref false positives (swaps excluded), not a full PPV. An
  exact 3×3 (and thus exact precision/FN) would require re-streaming the genotypes rather than reading
  GLIMPSE's per-truth-class counts.
- The perimeter also blocks the Summarize WDL's `conda install cyvcf2 pandas numpy matplotlib seaborn
  scipy`, so the eval prep builds a pip wheelhouse (cp311 manylinux) on the notebook VM, stages it to
  the bucket, and passes it as the `summarize_wheelhouse` input; the task installs offline with
  `pip install --no-index --find-links` (docker pinned to `python:3.11-slim` to match the cp311 wheels,
  `MPLBACKEND=Agg` for headless plotting). Override with `SUMMARIZE_WHEELHOUSE=gs://...` to supply a
  prebuilt wheelhouse (e.g. if this VM lacks PyPI egress). As with Concordance, re-run the eval with
  `FORCE_WF_REGISTER=true` once after updating the WDL.
- GLIMPSE2Summarize streams the panel and imputed VCFs with a **(CHROM,POS,REF,ALT) merge-join** rather
  than a positional `zip()`. The popped output and the popped panel can carry a different allele set or
  order at a multiallelic position (e.g. `chr1:10626 A>AG` vs `A>AGGCGCAG`), which used to raise
  "VCFs are out of sync"; the merge-join pairs records by exact allele within each position, skips and
  reports panel-only / imputed-only sites, and is order-independent and memory-bounded. The eval's
  site-matched panel no longer needs an exact record-count match with the imputed output.
- **vcfdist phasing eval** (`vcfdist_eval_wdl/VcfdistEvaluation.wdl`, launched by
  `prep_vcfdist_holdout198.sh`) scores the imputation per-sample against the 198's full-panel truth with
  [vcfdist](https://github.com/TimD1/vcfdist), which is **representation-aware** (it aligns query vs
  truth), so it reports **phasing** accuracy (switch / flip errors, phase blocks) alongside
  precision/recall — no allele-rep matching needed, and either the popped or un-popped output works as
  the eval VCF (default: un-popped, the direct phased Beagle haplotypes; `POPPED=true` for the popped
  output). Adapted from the Broad starter for this pipeline: a single multi-sample truth VCF instead of
  a per-sample array (the 198 share one truth), and every workflow input is a scalar or a
  newline-delimited file (`read_lines`) so the proven `wb` batch-input mechanism applies. The
  `SummarizeEvaluations` task installs pandas offline from a staged wheelhouse (perimeter blocks PyPI)
  and aggregates `phasing-summary.tsv` (all numeric columns) + switch/flip event and phase-block counts
  into `evaluation_summary.tsv`. Run: `bash prep_vcfdist_holdout198.sh` (set `MAX_SAMPLES=5` for a cheap
  test; `REF_FASTA=` to override the default public GRCh38; `CONF_BED=` for confident-region restriction;
  `FORCE_WF_REGISTER=true` after editing the WDL). vcfdist localizes the reference FASTA per sample, so a
  chr1-subset reference cuts cost on the 198-way scatter.
- `debug_check_outputs.sh` validates the imputation outputs end-to-end: that the popped output is truly
  popped (biallelic atomic, single-token `INFO/ID`, `FORMAT=GT:DS:GP`, `DS == gp1+2*gp2`); that each
  popped allele is **minimal/parsimonious** (no trimmable flanking bases -> bubbles really collapsed to
  simple atomic constituents, with a bubble-vs-popped allele-length comparison); that its sites are a
  subset of `panel_popped_vcf` / the id-split panel; that the un-popped output's sites are a subset of
  the bubble.split leaveout panel; and it checks the Summarize merge-join alignment. Read-only; streams
  panel sites (set `USER_PROJECT` for requester-pays reads).
- Sample-id namespace must match between panel/truth and ACAF; the prep errors if 0 of the 198 are
  found and reports the count otherwise.

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

# 3b) later: popped eval (truth side should become the id-split panel; current wiring is for un-popped)
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
- The eval WDLs fetch tools at runtime (concordance binary via wget, cyvcf2 via conda); if the
  VPC-SC perimeter blocks that, switch them to a prebuilt/offline approach.
- Sample-id namespace must match between panel/truth and ACAF; the prep errors if 0 of the 198 are
  found and reports the count otherwise.

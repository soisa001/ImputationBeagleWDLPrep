#!/usr/bin/env bash
# =============================================================================
# ImputeBeagleWithPop (pop OFF) — AoU lrWGS Phase-2 HELD-OUT eval imputation.
#
# Leave-N-out design (chr1) per the handoff:
#   REFERENCE = leaveout panel (held-out 490 removed)         -> bref3
#   TARGET    = the SAME 198 held-out AoU samples' REAL short-read (ACAF) biallelic
#               SNVs, PROJECTED onto the leaveout panel's bubble-allele representation
#               (minimal-rep matching, the GT analog of extract-bubble-PLs), UNPHASED
#                                          <- ACAF shards (gsutil cp from gs://, or mount)
#   We impute the projected ACAF SNV scaffold against the leaveout reference,
#   recovering indels/SVs, then score vs the FULL panel (truth) with
#   GLIMPSE2Concordance / Summarize (see prep_eval_holdout198.sh). This measures
#   EMPIRICAL ACAF->panel imputation accuracy, not the optimistic panel-derived target.
#   TRUTH_AOU (full-panel genotypes of the 198) is used here only to derive the id list.
#
# Pop is ON: the run emits BOTH the un-popped imputed_multi_sample_vcf
# (<OUT_BASE>.imputed.vcf.gz, bubble.split) AND the popped imputed_popped_vcf
# (<OUT_BASE>.imputed.popped.vcf.gz). The eval scores the un-popped output first
# (it matches the bubble.split panel truth); the popped output is for a later run.
#
# Reuses the validated bref3/refsrc-cleaning, genetic-map staging, inputs/CSV,
# and register/run machinery from prep_imputebeagle_aou.sh.
# =============================================================================
set -euo pipefail

# --- quiet, log-friendly wrappers --------------------------------------------
gsutil() { command gsutil -q "$@"; }
wb()     { TERM=dumb NO_COLOR=1 command wb "$@"; }
sanitize() { LC_ALL=C sed -E $'s/\x1b\\[[0-9;?]*[ -/]*[@-~]//g; s/\x1b\\][^\x07]*(\x07|\x1b\\\\)//g' \
             | LC_ALL=C tr -d '\000' | LC_ALL=C tr '\r' '\n'; }

# ============================ CONFIG =========================================
BUCKET_ID="${BUCKET_ID:-rw-migration-aou-rw-f178dfde}"
if [ -z "${BUCKET:-}" ]; then
  _B="$(wb resource resolve --name "${BUCKET_ID}" 2>/dev/null | tr -d '[:space:]')"
  case "${_B}" in
    gs://*) BUCKET="${_B}" ;;
    "")     BUCKET="gs://${BUCKET_ID}" ;;
    *)      BUCKET="gs://${_B}" ;;
  esac
fi
PROJECT="${GOOGLE_CLOUD_PROJECT:-$(gcloud config get-value project 2>/dev/null || true)}"
AOU_USER_PROJECT="${AOU_USER_PROJECT:-${PROJECT}}"   # billing project for requester-pays reads

# Run namespace. Bump RUN_ID (or pass RUN_ID=...) for a FRESH run that does NOT reuse
# previously staged artifacts or call-cached workflow outputs: it changes the gs://
# staging dir, the registered workflow name, and the run output path together.
RUN_ID="${RUN_ID:-holdout198_v2}"
WORK="${BUCKET}/imp_${RUN_ID}"                   # separate workspace dir from the production run
LOCAL="${HOME}/imp_${RUN_ID}"

CONTIGS=(${CONTIGS:-chr1})                        # held-out eval is chr1

# REFERENCE = leaveout panel (held-out removed). gs:// -> used in place.
PANEL_SRC_TMPL="${PANEL_SRC_TMPL-}"
[ -n "${PANEL_SRC_TMPL}" ] || PANEL_SRC_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_leaveout_vcf/aou_lr_phase2_v1.{contig}.bubble.split.leaveout.bcf'

# TRUTH source for the 198 (their full-panel genotypes; from check_holdout_panel.sh Step A).
# Used ONLY to build the sparse target. local path or gs://. Index (.csi/.tbi) expected alongside.
TRUTH_AOU="${TRUTH_AOU:-gs://cloned-rw-migration-aou-rw-f178dfde-wb-sharp-papaya-7463/vcf/truth.aou.chr1.bcf}"

# ---- ACAF (short-read) target source ----------------------------------------
# The 198 target is built from the SAME 198 held-out AoU samples' REAL short-read
# (ACAF) calls, then PROJECTED onto the leaveout panel's bubble-allele
# representation (minimal-rep matching, the GT analog of extract-bubble-PLs).
# Reference = leaveout panel (bref3); truth (eval) = full panel. This measures
# empirical ACAF->panel imputation, not the optimistic panel-derived target.
AOU_VCF_MOUNT="${AOU_VCF_MOUNT:-${HOME}/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/vcf}"
AOU_VCF_GS="${AOU_VCF_GS:-gs://vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/vcf}"   # for the mount-missing hint only
# 198 held-out AoU ids: local/gs file (one id per line). If unset, derived from TRUTH_AOU's sample list.
AOU_SAMPLES="${AOU_SAMPLES:-}"
# representation projection script (runs locally over the mount); auto-detected at ~ or ./ if unset.
PROJECT_LOCAL="${PROJECT_LOCAL:-}"

REF_PREFIX="${WORK}/ref/aou_lr_phase2_v1_leaveout"   # bref3 = <prefix>.<contig>.bref3 (+ .unique_variants)
MAPS_DIR="${WORK}/ref/genetic_maps/"             # trailing slash REQUIRED by the WDL
TARGET_DIR="${WORK}/target/"                     # gs:// dir (trailing slash)
MAPS_ZIP_URL="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"
BREF3_JAR_URL="${BREF3_JAR_URL:-https://faculty.washington.edu/browning/beagle/bref3.17Dec24.224.jar}"

N_TARGET="${N_TARGET:-198}"
TARGET_BASE="${TARGET_BASE:-aou_holdout198}"
OUT_BASE="${OUT_BASE:-aou_holdout198_${CONTIGS[0]}}"   # workflow output_basename -> <OUT_BASE>.imputed.vcf.gz

# ---- workflow registration + run (shares the ImputeBeagleWithPop WDL; pop OFF) ----
RUN_WORKFLOW="${RUN_WORKFLOW:-true}"
WF_NAME="${WF_NAME:-ImputeBeagleWithPop_${RUN_ID}}"   # distinct name -> fresh registration, no call-cache collision
MAIN_WDL="${MAIN_WDL:-ImputeBeagleWithPop.wdl}"
FORCE_WF_REGISTER="${FORCE_WF_REGISTER:-false}"
RUN_TAG="${RUN_TAG:-${RUN_ID}}"
OUTPUT_PATH="${OUTPUT_PATH:-imputebeagle-${RUN_TAG}-run}"   # runs land at gs://<bucket>/<OUTPUT_PATH>/ImputeBeagleWithPop/<uuid>/
USE_BATCH_CSV="${USE_BATCH_CSV:-true}"
SAMPLE_CHUNK_SIZE="${SAMPLE_CHUNK_SIZE:-1000}"   # 198 <= 1000 -> single sample chunk
PHASE_MEM_GB="${PHASE_MEM_GB:-60}"
IMPUTE_MEM_GB="${IMPUTE_MEM_GB:-65}"
ERROR_COUNT_OVERRIDE="${ERROR_COUNT_OVERRIDE:-0}"

# ---- bubble-pop / collision postprocess (ON: run produces BOTH outputs) ------
# With pop ON the workflow emits the un-popped imputed_multi_sample_vcf
# (<OUT_BASE>.imputed.vcf.gz, bubble.split = what the eval scores first) AND the
# popped imputed_popped_vcf (<OUT_BASE>.imputed.popped.vcf.gz, for the later run).
ENABLE_POP="${ENABLE_POP:-true}"
SITES_ONLY_TMPL="${SITES_ONLY_TMPL:-}"
[ -n "${SITES_ONLY_TMPL}" ] || SITES_ONLY_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_sites_only_vcf/aou_lr_phase2_v1.{contig}.bubble.split.sites.bcf'
ID_SPLIT_TMPL="${ID_SPLIT_TMPL:-}"
[ -n "${ID_SPLIT_TMPL}" ] || ID_SPLIT_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_id_split_vcf_gz/aou_lr_phase2_v1.{contig}.id.split.vcf.gz'
# pop-glimpse2 Rust engine: stage the source (.rs) + Cargo.toml (built in-task with cargo),
# OR a prebuilt binary (skips the build). Defaults look for the files next to this script.
POP_RS_LOCAL="${POP_RS_LOCAL:-}"      # Rust source; auto-resolve repo layout if unset
[ -n "${POP_RS_LOCAL}" ] || for d in ./pop_glimpse2_rust/pop-glimpse2.rs ./pop-glimpse2.rs "${HOME}/pop_glimpse2_rust/pop-glimpse2.rs"; do [ -s "$d" ] && { POP_RS_LOCAL="$d"; break; }; done
POP_CARGO_LOCAL="${POP_CARGO_LOCAL:-}"   # Cargo.toml (package pop-glimpse2); auto-resolve repo layout if unset
[ -n "${POP_CARGO_LOCAL}" ] || for d in ./pop_glimpse2_rust/Cargo.toml ./Cargo.toml "${HOME}/pop_glimpse2_rust/Cargo.toml"; do [ -s "$d" ] && { POP_CARGO_LOCAL="$d"; break; }; done
POP_BINARY_LOCAL="${POP_BINARY_LOCAL:-}"             # optional prebuilt binary (skips the build entirely)
POP_RS_URL="${POP_RS_URL:-}"                         # optional: fetch .rs from a URL instead of local
POP_CARGO_URL="${POP_CARGO_URL:-}"                   # optional: fetch Cargo.toml from a URL
POP_DOCKER="${POP_DOCKER:-us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-rust:v1}"  # needs bcftools + cargo (Rust toolchain)
# Build the pop-glimpse2 binary LOCALLY in this prep (static musl, no sudo/docker). The in-task
# cargo build can't reach crates.io inside the VPC-SC perimeter, so we build here on the notebook
# VM (which can) and ship the binary. Idempotent: reuses an existing build. Set POP_BUILD_LOCAL=false
# to skip and fall back to the in-task source build (only works if the task VM has crates.io egress).
POP_BUILD_LOCAL="${POP_BUILD_LOCAL:-true}"
POP_BUILD_DIR="${POP_BUILD_DIR:-${HOME}/pop-build}"
POP_BUILD_TARGET="${POP_BUILD_TARGET:-x86_64-unknown-linux-musl}"

# WDL bundle (same flat bundle as the production launcher).
WDL_SRC_GCS="${WDL_SRC_GCS:-gs://longreadsphase2imputation/wdl/imputation_beagle/}"
WDL_REL="wdl/imputebeagle_withpop"
WDL_GCS="${BUCKET}/${WDL_REL}/"

PANEL_THREADS="${PANEL_THREADS:-8}"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 8)}"
# target build runs LOCALLY: PARALLEL shards x SHARD_THREADS each (~THREADS). With
# TARGET_PULL_MODE=copy each bracketed shard is gsutil-cp'd to local SSD first (sliced
# parallel download, no FUSE), so peak local disk ~= PARALLEL * shard_size (each shard is
# removed right after its part is written). Lower PARALLEL if local disk is tight.
SHARD_THREADS="${SHARD_THREADS:-1}"   # threads per shard; PARALLEL=THREADS/SHARD_THREADS shards at once -> one shard per thread
PARALLEL="${PARALLEL:-$(( THREADS/SHARD_THREADS > 0 ? THREADS/SHARD_THREADS : 1 ))}"
MERGE_THREADS="${MERGE_THREADS:-${THREADS}}"
TARGET_CONCURRENT="${TARGET_CONCURRENT:-true}"   # overlap ACAF shard pulls with panel prep + bref3 build
# Shard pull transport (the slow part): 'copy' = cp each bracketed shard to local SSD then
# bcftools locally (fast; needs only AOU_VCF_GS + local disk, NOT the mount); 'mount' = read
# in place over the gcsfuse mount (AOU_VCF_MOUNT). Sample-subsetting reads the whole region
# regardless, and the bracketed shards are ~whole-file for the contig, so copy (parallel
# object download) beats gcsfuse streaming + .tbi random reads.
TARGET_PULL_MODE="${TARGET_PULL_MODE:-copy}"
# copy-mode prefetch buffer: a downloader pulls shards a little ahead of the reader so bcftools
# never waits on the network. DL_PARALLEL concurrent downloads; at most PREFETCH shards held
# downloaded-but-unread -> peak local disk ~= (PREFETCH + DL_PARALLEL) * shard_size.
DL_PARALLEL="${DL_PARALLEL:-${PARALLEL}}"
PREFETCH="${PREFETCH:-$(( PARALLEL + 2 ))}"
# keep only ACAF sites whose FILTER passes this list (bcftools -f). Default "PASS,." keeps
# PASS plus unfiltered (FILTER=".") sites (here "." implies PASS); other filters (LowQual,
# ExcessHet, ...) are dropped. Set TARGET_FILTER="PASS" to require an explicit PASS, or
# TARGET_FILTER= (empty) to disable site filtering.
TARGET_FILTER="${TARGET_FILTER-PASS,.}"
IDEMPOTENT="${IDEMPOTENT:-true}"

# ============================ preflight ======================================
for t in gsutil gcloud curl unzip awk wb python3 bcftools tabix bgzip java; do
  command -v "$t" >/dev/null 2>&1 || { echo "MISSING TOOL: $t"; exit 1; }
done
[ -n "${PROJECT}" ] || { echo "PROJECT empty (set GOOGLE_CLOUD_PROJECT)"; exit 1; }
mkdir -p "${LOCAL}"
gs_exists() { gsutil -q stat "$1" 2>/dev/null; }

GS_RETRIES="${GS_RETRIES:-5}"
retry() {
  local a=1
  until "$@"; do
    [ "$a" -ge "${GS_RETRIES}" ] && { echo ">> gave up after ${a} tries: $*" >&2; return 1; }
    echo ">> retry ${a}/${GS_RETRIES} (sleep $((a*10))s): $*" >&2
    sleep "$((a*10))"; a=$((a+1))
  done
}

# --- ACAF (requester-pays) transport helpers: copy-to-local instead of the gcsfuse mount ---
acaf_cp()  { retry gcloud storage cp --billing-project="${AOU_USER_PROJECT}" "$1" "$2"; }  # gs:// -> local (parallel striped dl)
acaf_cat() { gsutil -u "${AOU_USER_PROJECT}" cat "$1"; }            # gs:// -> stdout (streamed; caller may early-close)
acaf_ls()  { gsutil -u "${AOU_USER_PROJECT}" ls "$1"; }            # list a gs:// prefix (full gs:// urls)

PANEL_CANON=""; PANEL_LOCAL=""
ensure_panel_local() {                                 # download leaveout PANEL_CANON -> PANEL_LOCAL once
  [ -s "${PANEL_LOCAL}" ] && return 0
  echo ">> downloading reference (leaveout) locally: ${PANEL_CANON} -> ${PANEL_LOCAL}" >&2
  retry gsutil -u "${AOU_USER_PROJECT}" cp "${PANEL_CANON}" "${PANEL_LOCAL}.part"
  mv "${PANEL_LOCAL}.part" "${PANEL_LOCAL}"
}

# ============================ per-contig functions ===========================
resolve_panel() {                          # $1=contig -> PANEL_CANON (gs://) + PANEL_LOCAL (lazy)
  local c="$1" src; src="${PANEL_SRC_TMPL//\{contig\}/$c}"
  [[ "$src" == gs://* ]] || { echo "ERROR: panel must be gs://; got: $src"; exit 1; }
  gs_exists "$src" || { echo "ERROR: leaveout panel not found: $src"; exit 1; }
  PANEL_LOCAL="${LOCAL}/leaveout.${c}.bcf"
  local H="" a=1
  while :; do
    H="$(gsutil -u "${AOU_USER_PROJECT}" cat "$src" 2>/dev/null | bcftools view -h - 2>/dev/null | sed '/^#CHROM/q')"
    grep -q '^#CHROM' <<<"$H" && break
    [ "$a" -ge "${GS_RETRIES}" ] && { echo "ERROR: could not read panel header after ${a} tries: $src"; exit 1; }
    echo ">> [$c] panel header read retry ${a}/${GS_RETRIES}" >&2; sleep "$((a*5))"; a=$((a+1))
  done
  if grep -q "^##contig=<ID=${c}[,>]" <<<"$H"; then
    PANEL_CANON="$src"
  else
    echo ">> [$c] contig not chr-named; downloading + renaming locally" >&2
    PANEL_CANON="$src"; ensure_panel_local
    local act cm="${LOCAL}/rename_${c}.txt" ren="${LOCAL}/leaveout.${c}.vcf.gz"
    act="$(bcftools view -H "${PANEL_LOCAL}" 2>/dev/null | head -1 | cut -f1)"
    printf '%s\t%s\n' "$act" "$c" > "$cm"
    bcftools annotate --rename-chrs "$cm" -Oz -o "$ren" "${PANEL_LOCAL}"
    rm -f "${PANEL_LOCAL}"; PANEL_LOCAL="$ren"
  fi
}

prep_panel() {                             # $1=contig -> refsrc + unique_variants (cleaned leaveout)
  local c="$1"
  local REFSRC="${WORK}/ref/refsrc_noinfo.${c}.vcf.gz"
  local UNIQV="${REF_PREFIX}.${c}.unique_variants"
  if gs_exists "$REFSRC" && gs_exists "$UNIQV"; then echo ">> [$c] refsrc+uniqv exist; skipping panel prep"; return; fi
  echo ">> [$c] panel prep (local): strip INFO + norm -d exact + drop REF==ALT -> refsrc + unique_variants"
  ensure_panel_local
  local RL="${LOCAL}/refsrc_noinfo.${c}.vcf.gz" UL="${LOCAL}/${c}.unique_variants"
  bcftools annotate -x INFO -Ou "${PANEL_LOCAL}" \
    | bcftools norm -d exact --threads "${PANEL_THREADS}" -Ov - \
    | awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {n=split($5,a,","); bad=0; for(i=1;i<=n;i++) if(a[i]==$4){bad=1;break} if(bad){d++} else print} END{if(d) print "  [refsrc] dropped "d" REF==ALT record(s)" > "/dev/stderr"}' \
    | bcftools view -Oz -o "${RL}" --threads "${PANEL_THREADS}" -
  bcftools index -t --threads "${THREADS}" "${RL}"
  bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "${RL}" | LC_ALL=C sort -u --parallel="${THREADS}" -S 4G | sed '/^$/d' > "${UL}"
  gsutil cp "${RL}" "${REFSRC}"; gsutil cp "${RL}.tbi" "${REFSRC}.tbi"; gsutil cp "${UL}" "${UNIQV}"
}

build_bref3() {                            # $1=contig -> bref3 from cleaned leaveout refsrc
  local c="$1" BREF3="${REF_PREFIX}.${c}.bref3" REFSRC="${WORK}/ref/refsrc_noinfo.${c}.vcf.gz"
  if gs_exists "$BREF3"; then echo ">> [$c] bref3 exists"; return; fi
  local JAR="${LOCAL}/bref3.jar"
  [ -s "$JAR" ] || { echo ">> downloading bref3 jar"; curl -fsSL -o "$JAR" "${BREF3_JAR_URL}"; }
  local BL="${LOCAL}/${c}.bref3" RL="${LOCAL}/refsrc_noinfo.${c}.vcf.gz"
  [ -s "$RL" ] || gsutil cp "${REFSRC}" "${RL}"
  echo ">> [$c] bref3 (local): java -jar bref3 refsrc(cleaned) -> bref3"
  java -jar "$JAR" "${RL}" > "${BL}"
  [ -s "${BL}" ] || { echo ">> [$c] ERROR: bref3 build produced empty output"; exit 1; }
  gsutil cp "${BL}" "${BREF3}"
}

stage_map() {                              # $1=contig
  local c="$1" GMAP="${MAPS_DIR}plink.${c}.GRCh38.withchr.map"
  gs_exists "$GMAP" && { echo ">> [$c] map exists"; return; }
  echo ">> [$c] staging genetic map"
  [ -f "${LOCAL}/maps.zip" ] || curl -fL --retry 6 --retry-all-errors --retry-delay 10 "${MAPS_ZIP_URL}" -o "${LOCAL}/maps.zip"
  ( cd "${LOCAL}" && unzip -o maps.zip >/dev/null )
  local SRC; SRC="$(find "${LOCAL}" \( -name "plink.${c}.GRCh38.map" -o -name "plink.${c#chr}.GRCh38.map" \) -print | sort | head -1)"
  [ -n "$SRC" ] || { echo "no map for $c in zip"; exit 1; }
  awk 'BEGIN{OFS="\t"}{ if($1 !~ /^chr/) $1="chr"$1; print }' "$SRC" > "${LOCAL}/${c}.withchr.map"
  local MC; MC="$(head -1 "${LOCAL}/${c}.withchr.map" | cut -f1)"
  [ "$MC" = "$c" ] || { echo "ERROR: map for $c resolved to chrom '$MC' (wrong source file: $SRC)"; exit 1; }
  gsutil cp "${LOCAL}/${c}.withchr.map" "$GMAP"
}

# ---- ACAF shard ordering (first-contig order via tabix -l) for binary search ----
corder() {                                 # chr name -> numeric order key
  local x="${1#chr}"
  case "$x" in
    X) echo 23;; Y) echo 24;; M|MT) echo 25;;
    ''|*[!0-9]*) echo 99;;
    *) echo "$x";;
  esac
}
# first contig of a shard. mount: from the .tbi (tabix -l). copy: stream just the first data
# record's CHROM over the network (no index / full-file read needed). memoized per ref.
declare -A _FC_CACHE
shard_first_contig() {                     # $1 = shard ref (mount path or gs:// url)
  local s="$1"
  if [ -n "${_FC_CACHE[$s]+x}" ]; then printf '%s' "${_FC_CACHE[$s]}"; return; fi
  local fc
  if [ "${TARGET_PULL_MODE}" = "copy" ]; then
    fc="$( set +o pipefail; acaf_cat "$s" 2>/dev/null | bcftools view -H 2>/dev/null | head -1 | cut -f1 )"
  else
    fc="$( tabix -l "$s" 2>/dev/null | head -1 )"
  fi
  _FC_CACHE[$s]="$fc"; printf '%s' "$fc"
}
_firstord() { corder "$(shard_first_contig "${SHARDS[$1]}")"; }   # uses global SHARDS
_lb() {                                    # lower_bound: first shard with first-contig order >= $1
  local tgt="$1" lo=0 hi=${#SHARDS[@]} m fo
  while [ "$lo" -lt "$hi" ]; do
    m=$(( (lo + hi) / 2 )); fo="$(_firstord "$m")"
    if [ "$fo" -lt "$tgt" ]; then lo=$((m + 1)); else hi="$m"; fi
  done
  echo "$lo"
}

resolve_project_script() {                 # -> path to project_to_panel_rep.py (local), or empty
  local d
  for d in "${PROJECT_LOCAL}" "${HOME}/project_to_panel_rep.py" "./project_to_panel_rep.py" "${PWD}/project_to_panel_rep.py"; do
    [ -n "$d" ] && [ -s "$d" ] && { echo "$d"; return; }
  done
  echo ""
}

resolve_198_samples() {                    # $1=truth-local -> ${LOCAL}/holdout198.samples (the held-out AoU ids)
  local TA="$1" S="${LOCAL}/holdout198.samples"
  if [ -s "$S" ]; then echo "$S"; return; fi
  if [ -n "${AOU_SAMPLES}" ]; then
    if [[ "${AOU_SAMPLES}" == gs://* ]]; then retry gsutil -u "${AOU_USER_PROJECT}" cp "${AOU_SAMPLES}" "$S"
    else [ -s "${AOU_SAMPLES}" ] || { echo "ERROR: AOU_SAMPLES not found: ${AOU_SAMPLES}" >&2; exit 1; }; cp "${AOU_SAMPLES}" "$S"; fi
  else
    bcftools query -l "${TA}" > "$S"        # the 198 in the eval truth == the held-out AoU set
  fi
  sed -e 's/\r$//' -e '/^$/d' "$S" | awk '!seen[$0]++' > "${S}.tmp" && mv "${S}.tmp" "$S"
  [ -s "$S" ] || { echo "ERROR: empty 198 sample list" >&2; exit 1; }
  echo "$S"
}

pull_acaf_target() {                       # $1=contig -> raw ACAF biallelic target (NETWORK phase; safe to run in background)
  local c="$1"
  local TGT="${TARGET_DIR}${TARGET_BASE}_${c}.snps.vcf.gz"
  if [ "${IDEMPOTENT}" = "true" ] && gs_exists "${TGT}" && gs_exists "${TGT}.tbi"; then
    echo ">> [$c] target exists: ${TGT} (skip ACAF pull)"; return
  fi
  local ACAF="${LOCAL}/${TARGET_BASE}_${c}.acaf.snps.vcf.gz"
  if [ -s "$ACAF" ] && [ "$(bcftools query -l "$ACAF" 2>/dev/null | wc -l)" -gt 0 ]; then
    echo ">> [$c] ACAF concat exists, reusing: ${ACAF}"; return
  fi
  # The 198 held-out AoU samples' REAL short-read (ACAF) biallelic calls are extracted from the
  # ACAF shards (multiallelics split, GTs recoded). This phase depends only on TRUTH_AOU (for the id
  # list) and the ACAF shards (copied from gs:// by default) -- NOT on the panel prep / bref3
  # build -- so it runs concurrently
  # with them. The projection onto panel representation happens later in project_target.

  # --- resolve the 198 ids (default: TRUTH_AOU's sample list == the held-out set) ---
  local TA
  if [[ "${TRUTH_AOU}" == gs://* ]]; then
    TA="${LOCAL}/truth.aou.${c}.bcf"
    [ -s "${TA}" ] || retry gsutil -u "${AOU_USER_PROJECT}" cp "${TRUTH_AOU}" "${TA}"
    gs_exists "${TRUTH_AOU}.csi" && { [ -s "${TA}.csi" ] || gsutil cp "${TRUTH_AOU}.csi" "${TA}.csi"; }
  else
    [ -s "${TRUTH_AOU}" ] || { echo "ERROR: TRUTH_AOU not found: ${TRUTH_AOU}"; exit 1; }
    TA="${TRUTH_AOU}"
  fi
  [ -s "${TA}.csi" ] || [ -s "${TA}.tbi" ] || bcftools index "${TA}"
  local SAMP_LOC; SAMP_LOC="$(resolve_198_samples "$TA")"
  echo ">> [$c] target samples: $(wc -l < "$SAMP_LOC") held-out AoU ids"

  # --- discover ACAF shards (copy: list gs://; mount: ls the gcsfuse mount) ---
  SHARDS=()
  if [ "${TARGET_PULL_MODE}" = "copy" ]; then
    mapfile -t SHARDS < <(acaf_ls "${AOU_VCF_GS}/" 2>/dev/null | grep -E '\.vcf\.bgz$' | sort)
    if [ "${#SHARDS[@]}" -eq 0 ]; then
      echo "ERROR: no *.vcf.bgz under ${AOU_VCF_GS} (requester-pays; billing project=${AOU_USER_PROJECT}). Check AOU_VCF_GS / AOU_USER_PROJECT, or use TARGET_PULL_MODE=mount."
      exit 1
    fi
    echo ">> [$c] ${#SHARDS[@]} ACAF shards listed at ${AOU_VCF_GS} (pull mode: copy)"
  else
    [ -d "${AOU_VCF_MOUNT}" ] || { echo "ERROR: ACAF mount not found: ${AOU_VCF_MOUNT} (set AOU_VCF_MOUNT to the acaf_threshold vcf/ dir, or use TARGET_PULL_MODE=copy)"; exit 1; }
    mapfile -t SHARDS < <(ls -U -1 "${AOU_VCF_MOUNT}" 2>/dev/null | grep -E '\.vcf\.bgz$' | sed "s#^#${AOU_VCF_MOUNT}/#" | sort)
    if [ "${#SHARDS[@]}" -eq 0 ]; then
      echo "ERROR: no *.vcf.bgz under ${AOU_VCF_MOUNT}"
      gsutil -u "${AOU_USER_PROJECT}" ls "${AOU_VCF_GS}/" 2>/dev/null | grep -qE '\.vcf\.bgz$' \
        && echo "  -> shards exist in ${AOU_VCF_GS} (requester-pays); set TARGET_PULL_MODE=copy, or remount with --billing-project + --implicit-dirs and set AOU_VCF_MOUNT."
      exit 1
    fi
    echo ">> [$c] ${#SHARDS[@]} ACAF shards visible under the mount (pull mode: mount)"
  fi

  # --- bracket shards for contig c (first-contig order via tabix -l) ---
  local want; want="$(corder "$c")"
  local L R; L="$(_lb "$want")"; R="$(_lb $((want + 1)))"
  local lo_i=$(( L - 1 < 0 ? 0 : L - 1 )) hi_i=$(( R < ${#SHARDS[@]} ? R : ${#SHARDS[@]} - 1 ))
  echo ">> [$c] shard bracket [${lo_i}..${hi_i}] of ${#SHARDS[@]}"
  local sel=() i
  for (( i=lo_i; i<=hi_i; i++ )); do sel+=("$i"); done

  # --- per-shard parallel extract: 198 samples, SPLIT multiallelics -> biallelic (idempotent .tmp -> mv) ---
  # bcftools norm -m -any decomposes multiallelic ACAF sites into biallelic records and recodes GTs
  # (a het-of-two-alts 1/2 -> 1/0 in one record, 0/1 in the other; biallelic-valid). We deliberately do
  # NOT pre-filter -v snps here: padded SNVs from a split (e.g. AT>AG) stay len>1 and would be dropped by
  # -v snps, but the projection's minimal-rep reduces them to a SNV and matches. (pipefail: a mid-stream
  # read error fails the whole pipe and the until-loop retries.)
  local PARTS="${LOCAL}/target_parts_${c}"; mkdir -p "$PARTS"

  _extract_part() {                          # $1=shard-index $2=reader-path -> writes PARTS part (idempotent)
    local i="$1" reader="$2" part a=1 filt=()
    part="${PARTS}/part_$(printf '%010d' "$i").vcf.gz"
    [ -n "${TARGET_FILTER}" ] && filt=(-f "${TARGET_FILTER}")   # keep only PASS sites (site FILTER)
    until bcftools view "${filt[@]}" -r "$c" -S "$SAMP_LOC" --force-samples --threads "${SHARD_THREADS}" "$reader" -Ou \
            | bcftools norm -m -any --threads "${SHARD_THREADS}" -Oz -o "${part}.tmp"; do
      rm -f "${part}.tmp"
      [ "$a" -ge "${GS_RETRIES}" ] && { echo ">> [$c] shard idx $i extract FAILED after ${a} tries" >&2; return 1; }
      echo ">> [$c] shard idx $i extract retry ${a}/${GS_RETRIES}; sleep $((a*5))s" >&2
      sleep "$((a*5))"; a=$((a+1))
    done
    mv "${part}.tmp" "$part"
  }

  if [ "${TARGET_PULL_MODE}" = "copy" ]; then
    # --- prefetch buffer: a downloader pulls shards (gcloud storage cp) a little AHEAD of the
    #     reader into a bounded buffer, so bcftools never waits on the network. The bucket .tbi/.csi
    #     index is downloaded alongside each shard and used for the -r region query. Backpressure
    #     keeps at most PREFETCH downloaded-but-unread shards on disk; DL_PARALLEL download at once
    #     -> peak local disk ~= (PREFETCH + DL_PARALLEL) * shard_size (freed right after each read). ---
    local BUF="${LOCAL}/acaf_buf_${c}"; rm -rf "$BUF"; mkdir -p "$BUF"

    # producer: download shard + its index ahead of the reader (in order, DL_PARALLEL concurrent)
    (
      for i in "${sel[@]}"; do
        [ -s "${PARTS}/part_$(printf '%010d' "$i").vcf.gz" ] && continue   # already extracted (resume)
        # backpressure: don't run more than PREFETCH ahead of the reader
        while [ "$(( $(ls "${BUF}"/*.dl.ok 2>/dev/null | wc -l) - $(ls "${BUF}"/*.rd.ok 2>/dev/null | wc -l) ))" -ge "${PREFETCH}" ]; do sleep 0.5; done
        while [ "$(jobs -rp | wc -l)" -ge "${DL_PARALLEL}" ]; do wait -n 2>/dev/null || break; done
        dl="${BUF}/$(printf '%010d' "$i")"; shard="${SHARDS[$i]}"
        (
          set +e
          acaf_cp "$shard" "${dl}.vcf.bgz.part" && mv "${dl}.vcf.bgz.part" "${dl}.vcf.bgz"; rc=$?
          if [ "$rc" -eq 0 ]; then                                  # fetch the bucket index (.tbi or .csi)
            if   acaf_cp "${shard}.tbi" "${dl}.idx.part" 2>/dev/null; then mv "${dl}.idx.part" "${dl}.vcf.bgz.tbi"
            elif acaf_cp "${shard}.csi" "${dl}.idx.part" 2>/dev/null; then mv "${dl}.idx.part" "${dl}.vcf.bgz.csi"
            else tabix -f -p vcf "${dl}.vcf.bgz"; rc=$?; fi
          fi
          [ "$rc" -eq 0 ] && touch "${dl}.dl.ok" || touch "${dl}.dl.fail"
        ) &
      done
      wait
    ) &
    local DL_PID=$!; disown "$DL_PID" 2>/dev/null || true   # keep it out of the consumer's `wait -n`

    # consumer: extract each shard as soon as its download lands, then free its buffer slot
    local running=0 i dl
    for i in "${sel[@]}"; do
      [ -s "${PARTS}/part_$(printf '%010d' "$i").vcf.gz" ] && continue
      dl="${BUF}/$(printf '%010d' "$i")"
      while [ ! -e "${dl}.dl.ok" ]; do
        [ -e "${dl}.dl.fail" ] && { echo ">> [$c] shard idx $i download FAILED" >&2; kill "$DL_PID" 2>/dev/null; exit 1; }
        kill -0 "$DL_PID" 2>/dev/null || { [ -e "${dl}.dl.ok" ] || { echo ">> [$c] downloader exited before shard idx $i" >&2; exit 1; }; }
        sleep 0.5
      done
      (
        _extract_part "$i" "${dl}.vcf.bgz" || { rm -f "${dl}".* ; touch "${dl}.rd.ok"; exit 1; }
        rm -f "${dl}.vcf.bgz" "${dl}.vcf.bgz.tbi" "${dl}.vcf.bgz.csi"   # free disk
        touch "${dl}.rd.ok"                                            # release a buffer slot
      ) &
      running=$(( running + 1 ))
      if [ "$running" -ge "${PARALLEL}" ]; then wait -n 2>/dev/null || wait; running=$(( running - 1 )); fi
    done
    wait
    while kill -0 "$DL_PID" 2>/dev/null; do sleep 0.2; done    # producer is done once all reads completed
    rm -rf "$BUF"
  else
    # mount mode: read shards in place over the gcsfuse mount (original behavior)
    local running=0 i
    for i in "${sel[@]}"; do
      [ -s "${PARTS}/part_$(printf '%010d' "$i").vcf.gz" ] && continue
      ( _extract_part "$i" "${SHARDS[$i]}" || exit 1 ) &
      running=$(( running + 1 ))
      if [ "$running" -ge "${PARALLEL}" ]; then wait -n 2>/dev/null || wait; running=$(( running - 1 )); fi
    done
    wait
  fi

  # --- collect parts in order, concat to the raw ACAF target ---
  local PARTLIST=()
  for i in "${sel[@]}"; do
    part="${PARTS}/part_$(printf '%010d' "$i").vcf.gz"
    [ -s "$part" ] || { echo "ERROR: missing part for shard index $i ($part); rerun to resume"; exit 1; }
    PARTLIST+=("$part")
  done
  bcftools concat --threads "${MERGE_THREADS}" -Oz -o "${ACAF}.tmp" "${PARTLIST[@]}" && mv "${ACAF}.tmp" "$ACAF"
  local NS; NS="$(bcftools query -l "$ACAF" | wc -l)"
  [ "$NS" -gt 0 ] || { echo "ERROR: 0 of the requested ids found in ACAF (id namespace mismatch between panel/truth and ACAF? provide AOU_SAMPLES with ACAF-matching ids)"; exit 1; }
  echo ">> [$c] ACAF SNVs (pre-projection): $(bcftools index -n "$ACAF" 2>/dev/null || echo '?') records x ${NS} samples"
}

project_target() {                         # $1=contig -> PROJECT ACAF onto panel rep -> final target (needs refsrc + ACAF concat)
  local c="$1"
  local TGT="${TARGET_DIR}${TARGET_BASE}_${c}.snps.vcf.gz"
  if [ "${IDEMPOTENT}" = "true" ] && gs_exists "${TGT}" && gs_exists "${TGT}.tbi"; then
    echo ">> [$c] target exists: ${TGT}"; return
  fi
  local ACAF="${LOCAL}/${TARGET_BASE}_${c}.acaf.snps.vcf.gz"
  [ -s "$ACAF" ] || { echo "ERROR: [$c] ACAF concat missing: ${ACAF} (pull phase did not complete)"; exit 1; }

  # --- projection script (runs locally) ---
  local PROJ; PROJ="$(resolve_project_script)"
  [ -n "$PROJ" ] || { echo "ERROR: project_to_panel_rep.py not found (set PROJECT_LOCAL=/path/to/project_to_panel_rep.py)"; exit 1; }
  echo ">> [$c] projection script: ${PROJ}"

  # --- panel sites TSV = the cleaned leaveout refsrc actually in the bref3 (produced by prep_panel) ---
  local RL="${LOCAL}/refsrc_noinfo.${c}.vcf.gz" REFSRC="${WORK}/ref/refsrc_noinfo.${c}.vcf.gz"
  [ -s "$RL" ] || retry gsutil cp "${REFSRC}" "${RL}"
  local PSITES="${LOCAL}/panel_sites.${c}.tsv.gz"
  if [ ! -s "$PSITES" ]; then
    echo ">> [$c] building panel sites TSV from leaveout refsrc"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${RL}" | bgzip -c > "${PSITES}.tmp" && mv "${PSITES}.tmp" "${PSITES}"
  fi

  # --- PROJECT onto panel bubble-allele representation, unphase, sort ---
  local TL="${LOCAL}/${TARGET_BASE}_${c}.snps.vcf.gz"
  echo ">> [$c] projecting ACAF SNVs onto panel representation (minimal-rep matching)"
  bcftools view "$ACAF" -Ov \
    | python3 "$PROJ" "${PSITES}" \
    | bcftools sort -Oz -o "${TL}"
  bcftools index -t "${TL}"
  echo ">> [$c] target (projected): $(bcftools query -l "${TL}" | wc -l) samples x $(bcftools index -n "${TL}") SNVs"
  gsutil cp "${TL}" "${TGT}"; gsutil cp "${TL}.tbi" "${TGT}.tbi"
}

# ============================ build reference + target =======================
for c in "${CONTIGS[@]}"; do
  echo "===================== ${c} ====================="
  resolve_panel "$c"
  if [ "${TARGET_CONCURRENT}" = "true" ]; then
    # ACAF shard pulls (network-bound) overlap the CPU-bound panel prep + bref3 build (the slow
    # part that was blocking the copies). The pull needs only TRUTH_AOU + the mount, not the panel.
    echo ">> [$c] launching ACAF target pull in background; building reference concurrently"
    pull_acaf_target "$c" & _acaf_pid=$!
    prep_panel "$c"
    build_bref3 "$c"
    stage_map "$c"
    if ! wait "${_acaf_pid}"; then echo "ERROR: [$c] ACAF target pull failed (see log above)"; exit 1; fi
    project_target "$c"
  else
    prep_panel "$c"; build_bref3 "$c"; stage_map "$c"
    pull_acaf_target "$c"; project_target "$c"
  fi
done

# ---- stage the pop-glimpse2 Rust engine + resolve per-contig pop inputs -----
POP_JSON=""
if [ "${ENABLE_POP}" = "true" ]; then
  if [ "${#CONTIGS[@]}" -ne 1 ]; then
    echo ">> WARN: pop wired for a single contig; using ${CONTIGS[0]} pop panels only"
  fi
  stage_pop() {   # $1=local $2=url $3=dest-gcs ; rc=1 if neither local nor url provided
    gs_exists "$3" && return 0
    if [ -n "$1" ] && [ -s "$1" ]; then echo ">> staging $(basename "$3") (local) -> $3"; gsutil cp "$1" "$3"
    elif [ -n "$2" ]; then echo ">> downloading $(basename "$3") -> $3"; curl -fsSL -o "${LOCAL}/$(basename "$3")" "$2" || { echo "ERROR: cannot fetch $2"; exit 1; }; gsutil cp "${LOCAL}/$(basename "$3")" "$3"
    else return 1; fi
  }
  build_pop_local() {   # echoes built static-musl binary path on stdout (progress on stderr); empty on failure
    [ "${POP_BUILD_LOCAL}" = "true" ] || { echo ""; return; }
    local bin="${POP_BUILD_DIR}/target/${POP_BUILD_TARGET}/release/pop-glimpse2"
    if [ -x "$bin" ]; then echo "$bin"; return; fi                       # idempotent: reuse existing build
    local rs="" cargo=""
    for d in "${POP_RS_LOCAL}" "${HOME}/pop_glimpse2_rust/pop-glimpse2.rs" "./pop_glimpse2_rust/pop-glimpse2.rs" "./pop-glimpse2.rs"; do
      [ -n "$d" ] && [ -s "$d" ] && { rs="$d"; break; }; done
    for d in "${POP_CARGO_LOCAL}" "${HOME}/pop_glimpse2_rust/Cargo.toml" "./pop_glimpse2_rust/Cargo.toml" "./Cargo.toml"; do
      [ -n "$d" ] && [ -s "$d" ] && { cargo="$d"; break; }; done
    [ -n "$rs" ] && [ -n "$cargo" ] || { echo "" ; return; }            # no source -> can't build locally
    if ! command -v cargo >/dev/null 2>&1; then                          # user-local rustup (no sudo)
      if [ -s "$HOME/.cargo/env" ]; then . "$HOME/.cargo/env"
      else
        command -v curl >/dev/null 2>&1 || { echo "" >&2; echo ""; return; }
        echo ">> installing rustup (user-local) for the pop build..." >&2
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y >&2 || { echo ""; return; }
        . "$HOME/.cargo/env"
      fi
    fi
    rustup target add "${POP_BUILD_TARGET}" >&2 2>/dev/null || true
    rm -rf "${POP_BUILD_DIR}"; mkdir -p "${POP_BUILD_DIR}/src/bin" "${POP_BUILD_DIR}/.cargo"
    cp "$rs" "${POP_BUILD_DIR}/src/bin/pop-glimpse2.rs"
    cp "$cargo" "${POP_BUILD_DIR}/Cargo.toml"
    printf '\n[workspace]\n' >> "${POP_BUILD_DIR}/Cargo.toml"            # don't climb to a stray ~/Cargo.toml
    printf '[target.%s]\nlinker = "rust-lld"\n' "${POP_BUILD_TARGET}" > "${POP_BUILD_DIR}/.cargo/config.toml"
    echo ">> building pop-glimpse2 (static ${POP_BUILD_TARGET})..." >&2
    ( cd "${POP_BUILD_DIR}" && rm -f Cargo.lock && cargo build --release --target "${POP_BUILD_TARGET}" ) >&2 || { echo ""; return; }
    [ -x "$bin" ] && echo "$bin" || echo ""
  }
  # If no prebuilt binary was provided, BUILD one locally (static musl) from the Rust source.
  if [ -z "${POP_BINARY_LOCAL}" ]; then
    _pop_built="$(build_pop_local || true)"
    if [ -n "${_pop_built}" ] && [ -s "${_pop_built}" ]; then
      POP_BINARY_LOCAL="${_pop_built}"
      echo ">> built pop-glimpse2 locally: ${POP_BINARY_LOCAL}"
    elif [ "${POP_BUILD_LOCAL}" = "true" ]; then
      echo ">> WARN: local pop build unavailable (no rustc/source?); falling back to in-task source build (needs crates.io egress)"
    fi
  fi
  if [ -n "${POP_BINARY_LOCAL}" ] && [ -s "${POP_BINARY_LOCAL}" ]; then
    POP_BIN_GCS="${WORK}/prep/$(basename "${POP_BINARY_LOCAL}")"
    gs_exists "${POP_BIN_GCS}" || { echo ">> staging prebuilt pop-glimpse2 binary -> ${POP_BIN_GCS}"; gsutil cp "${POP_BINARY_LOCAL}" "${POP_BIN_GCS}"; }
    POP_ENGINE_JSON=",
  \"ImputeBeagleWithPop.pop_glimpse2_binary\": \"${POP_BIN_GCS}\""
    POP_ENGINE_INLINE=",pop_glimpse2_binary=${POP_BIN_GCS}"
    echo ">> pop engine: prebuilt binary ${POP_BIN_GCS}"
  else
    POP_RS_GCS="${WORK}/prep/pop-glimpse2.rs"
    POP_CARGO_GCS="${WORK}/prep/Cargo.toml"
    stage_pop "${POP_RS_LOCAL}" "${POP_RS_URL}" "${POP_RS_GCS}" \
      || { echo "ERROR: pop Rust source not found. Set POP_RS_LOCAL=<pop-glimpse2.rs> (or POP_RS_URL=, or POP_BINARY_LOCAL=)."; exit 1; }
    stage_pop "${POP_CARGO_LOCAL}" "${POP_CARGO_URL}" "${POP_CARGO_GCS}" \
      || { echo "ERROR: Cargo.toml not found. Set POP_CARGO_LOCAL=<Cargo.toml> (or POP_CARGO_URL=)."; exit 1; }
    POP_ENGINE_JSON=",
  \"ImputeBeagleWithPop.pop_glimpse2_script\": \"${POP_RS_GCS}\",
  \"ImputeBeagleWithPop.pop_glimpse2_cargo_toml\": \"${POP_CARGO_GCS}\""
    POP_ENGINE_INLINE=",pop_glimpse2_script=${POP_RS_GCS},pop_glimpse2_cargo_toml=${POP_CARGO_GCS}"
    echo ">> pop engine: build from source ${POP_RS_GCS} + ${POP_CARGO_GCS}"
  fi
  SITES_ONLY="${SITES_ONLY_TMPL//\{contig\}/${CONTIGS[0]}}"
  ID_SPLIT="${ID_SPLIT_TMPL//\{contig\}/${CONTIGS[0]}}"
  for f in "${SITES_ONLY}" "${SITES_ONLY}.csi" "${ID_SPLIT}" "${ID_SPLIT}.tbi"; do
    gs_exists "$f" || echo ">> WARN: pop input not found: $f"
  done
  POP_JSON=",
  \"ImputeBeagleWithPop.panel_bubble_split_sites_only_vcf\": \"${SITES_ONLY}\",
  \"ImputeBeagleWithPop.panel_bubble_split_sites_only_vcf_idx\": \"${SITES_ONLY}.csi\",
  \"ImputeBeagleWithPop.panel_id_split_vcf_gz\": \"${ID_SPLIT}\",
  \"ImputeBeagleWithPop.panel_id_split_vcf_gz_tbi\": \"${ID_SPLIT}.tbi\"${POP_ENGINE_JSON},
  \"ImputeBeagleWithPop.pop_docker\": \"${POP_DOCKER}\",
  \"ImputeBeagleWithPop.pop_region\": \"${CONTIGS[0]}\""
  echo ">> pop enabled: sites=${SITES_ONLY} id_split=${ID_SPLIT} region=${CONTIGS[0]}"
fi

# ============================ inputs: JSON + batch CSV + column map ===========
SCS="${SAMPLE_CHUNK_SIZE}"
CJ=""; for c in "${CONTIGS[@]}"; do CJ="${CJ}\"${c}\","; done; CJ="[${CJ%,}]"
mkdir -p "${LOCAL}/batch"
IN_JSON="${LOCAL}/batch/${TARGET_BASE}.inputs.json"
IN_CSV="${LOCAL}/batch/${TARGET_BASE}.inputs.csv"
IN_COLS="${LOCAL}/batch/${TARGET_BASE}.columns.json"
ECO_JSON=""
[ -n "${ERROR_COUNT_OVERRIDE}" ] && ECO_JSON=",
  \"ImputeBeagleWithPop.error_count_override\": ${ERROR_COUNT_OVERRIDE}"
cat > "${IN_JSON}" <<JSON
{
  "ImputeBeagleWithPop.multi_sample_vcf": "${TARGET_DIR}${TARGET_BASE}_${CONTIGS[0]}.snps.vcf.gz",
  "ImputeBeagleWithPop.ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
  "ImputeBeagleWithPop.contigs": ${CJ},
  "ImputeBeagleWithPop.reference_panel_path_prefix": "${REF_PREFIX}",
  "ImputeBeagleWithPop.genetic_maps_path": "${MAPS_DIR}",
  "ImputeBeagleWithPop.output_basename": "${OUT_BASE}",
  "ImputeBeagleWithPop.sample_chunk_size": ${SCS},
  "ImputeBeagleWithPop.beagle_phase_memory_in_gb": ${PHASE_MEM_GB},
  "ImputeBeagleWithPop.beagle_impute_memory_in_gb": ${IMPUTE_MEM_GB}${ECO_JSON}${POP_JSON}
}
JSON
python3 - "${IN_JSON}" "${IN_CSV}" "${IN_COLS}" <<'PY'
import json, csv, sys
d = json.load(open(sys.argv[1]))
header, row, cols = [], [], {}
for k, v in d.items():
    col = k.split(".", 1)[1] if "." in k else k
    cols[k] = col; header.append(col)
    row.append(json.dumps(v) if isinstance(v, (list, dict, bool)) else str(v))
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f); w.writerow(header); w.writerow(row)
json.dump(cols, open(sys.argv[3], "w"), indent=2)
PY
BATCH_GCS="${WORK}/batch/"
BATCH_CSV_REL="${BATCH_GCS#${BUCKET}/}${TARGET_BASE}.inputs.csv"
COLS_URI="${BATCH_GCS}${TARGET_BASE}.columns.json"
gsutil cp "${IN_JSON}" "${IN_CSV}" "${IN_COLS}" "${BATCH_GCS}"

# ============================ WDL bundle =====================================
if [ -z "${WDL_SRC_LOCAL:-}" ]; then
  for d in "${HOME}/imputationbeagle_wdl_flat" "./imputationbeagle_wdl_flat" "${PWD}/imputationbeagle_wdl_flat"; do
    if [ -f "${d}/${MAIN_WDL}" ]; then WDL_SRC_LOCAL="$d"; echo ">> auto-detected local WDL bundle: ${d}"; break; fi
  done
fi
if [ -n "${WDL_SRC_LOCAL:-}" ]; then
  echo ">> uploading WDLs ${WDL_SRC_LOCAL%/}/*.wdl -> ${WDL_GCS}"
  gsutil cp "${WDL_SRC_LOCAL%/}/"*.wdl "${WDL_GCS}"
elif gs_exists "${WDL_GCS}${MAIN_WDL}"; then
  echo ">> WDLs already staged in ${WDL_GCS}; no copy needed"
elif [ "${WDL_SRC_GCS%/}" != "${WDL_GCS%/}" ]; then
  echo ">> copying WDLs ${WDL_SRC_GCS}*.wdl -> ${WDL_GCS}"
  gsutil cp "${WDL_SRC_GCS}"*.wdl "${WDL_GCS}" \
    || { echo "ERROR: cannot copy WDLs from ${WDL_SRC_GCS} (cross-perimeter? VPC-SC 403)."; echo "       put the bundle on this VM at ~/imputationbeagle_wdl_flat/ (or set WDL_SRC_LOCAL=<dir>) and rerun."; exit 1; }
else
  echo ">> WDLs already in this workspace bucket (${WDL_GCS}); no copy needed"
fi

# ============================ register + run =================================
if [ "${RUN_WORKFLOW}" = "true" ]; then
  miss=0
  for w in "${MAIN_WDL}" ImputationBeagleStructs.wdl ImputationBeagleTasks.wdl ImputationTasks.wdl; do
    gs_exists "${WDL_GCS}${w}" || { echo "ERROR: missing ${WDL_GCS}${w}"; miss=1; }
  done
  [ "$miss" -eq 0 ] || { echo "       set WDL_SRC_LOCAL=<dir with the .wdl bundle> to upload them from this VM."; exit 1; }
  if [ "${FORCE_WF_REGISTER}" = "true" ] && wb workflow list 2>/dev/null | grep -qw "${WF_NAME}"; then
    echo ">> FORCE_WF_REGISTER: deleting stale ${WF_NAME}"
    wb workflow delete --workflow="${WF_NAME}" --quiet 2>/dev/null || true
  fi
  if wb workflow list 2>/dev/null | grep -qw "${WF_NAME}"; then
    echo ">> workflow ${WF_NAME} already registered"
  else
    echo ">> registering workflow ${WF_NAME} -> ${WDL_GCS}${MAIN_WDL}"
    wb workflow create --bucket-id="${BUCKET_ID}" \
      --path="${WDL_REL}/${MAIN_WDL}" --workflow="${WF_NAME}"
  fi
  echo ">> launching ${WF_NAME} (pop=${ENABLE_POP}): output-path=${OUTPUT_PATH}"
  RUN_LOG="${LOCAL}/wb_run_${OUT_BASE}.log"
  if [ "${USE_BATCH_CSV}" = "true" ]; then
    wb workflow job run \
      --workflow="${WF_NAME}" \
      --output-bucket-id="${BUCKET_ID}" \
      --output-path="${OUTPUT_PATH}" \
      --batch-input-bucket-id="${BUCKET_ID}" \
      --batch-input-csv-path="${BATCH_CSV_REL}" \
      --column-mapping-uri="${COLS_URI}" \
      --read-from-cache --write-to-cache 2>&1 | sanitize | tee "${RUN_LOG}" \
    || { echo ">> batch-CSV run failed (inputs at ${BATCH_GCS}). Try USE_BATCH_CSV=false."; exit 1; }
  else
    INPUTS="multi_sample_vcf=${TARGET_DIR}${TARGET_BASE}_${CONTIGS[0]}.snps.vcf.gz"
    INPUTS="${INPUTS},ref_dict=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    INPUTS="${INPUTS},contigs=${CJ}"
    INPUTS="${INPUTS},reference_panel_path_prefix=${REF_PREFIX}"
    INPUTS="${INPUTS},genetic_maps_path=${MAPS_DIR}"
    INPUTS="${INPUTS},output_basename=${OUT_BASE}"
    INPUTS="${INPUTS},sample_chunk_size=${SCS}"
    INPUTS="${INPUTS},beagle_phase_memory_in_gb=${PHASE_MEM_GB}"
    INPUTS="${INPUTS},beagle_impute_memory_in_gb=${IMPUTE_MEM_GB}"
    [ -n "${ERROR_COUNT_OVERRIDE}" ] && INPUTS="${INPUTS},error_count_override=${ERROR_COUNT_OVERRIDE}"
    if [ "${ENABLE_POP}" = "true" ]; then
      INPUTS="${INPUTS},panel_bubble_split_sites_only_vcf=${SITES_ONLY},panel_bubble_split_sites_only_vcf_idx=${SITES_ONLY}.csi"
      INPUTS="${INPUTS},panel_id_split_vcf_gz=${ID_SPLIT},panel_id_split_vcf_gz_tbi=${ID_SPLIT}.tbi"
      INPUTS="${INPUTS}${POP_ENGINE_INLINE},pop_docker=${POP_DOCKER},pop_region=${CONTIGS[0]}"
    fi
    wb workflow job run \
      --workflow="${WF_NAME}" \
      --output-bucket-id="${BUCKET_ID}" \
      --output-path="${OUTPUT_PATH}" \
      --inputs="${INPUTS}" \
      --read-from-cache --write-to-cache 2>&1 | sanitize | tee "${RUN_LOG}" \
    || { echo ">> inline --inputs run failed; see 'wb workflow job run --help'."; exit 1; }
  fi
  JOB_ID="$(grep -m1 'Job ID:' "${RUN_LOG}" 2>/dev/null | sed 's/.*Job ID:[[:space:]]*//' | tr -d '[:space:]' || true)"
  [ -n "${JOB_ID}" ] && echo ">> job id: ${JOB_ID}"
  JOB_OUT="$(grep -m1 'Output bucket path:' "${RUN_LOG}" 2>/dev/null | sed 's/.*Output bucket path:[[:space:]]*//' | tr -d '[:space:]' || true)"
  [ -n "${JOB_OUT}" ] && echo ">> job output dir: ${JOB_OUT}"
fi

echo "================= DONE ================="
echo "reference bref3:  ${REF_PREFIX}.${CONTIGS[0]}.bref3   (from leaveout panel)"
echo "scaffold target:  ${TARGET_DIR}${TARGET_BASE}_${CONTIGS[0]}.snps.vcf.gz"
echo "workflow:         ${WF_NAME} (pop=${ENABLE_POP})  ->  ${BUCKET}/${OUTPUT_PATH}/ImputeBeagleWithPop/<uuid>/"
echo "un-popped (eval now):  <...>/${OUT_BASE}.imputed.vcf.gz         = imputed_multi_sample_vcf"
[ "${ENABLE_POP}" = "true" ] && echo "popped    (eval later): <...>/${OUT_BASE}.imputed.popped.vcf.gz  = imputed_popped_vcf"
echo
echo "Once the run finishes, score it with:"
echo "  IMPUTED_VCF=<the .imputed.vcf.gz above> bash prep_eval_holdout198.sh"

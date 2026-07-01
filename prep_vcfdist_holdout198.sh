#!/usr/bin/env bash
# =============================================================================
# vcfdist phasing + precision/recall eval prep -- holdout198 imputation (chr1).
#
# Registers + runs vcfdist_eval_wdl/VcfdistEvaluation.wdl on the holdout Beagle
# output, per-sample, vs the 198's full-panel truth. vcfdist is representation-
# aware (it aligns query vs truth), so it scores PHASING (switch/flip errors,
# phase blocks) as well as precision/recall -- no allele-rep matching needed.
#
# Defaults to the UN-POPPED output (the direct, phased Beagle haplotypes); set
# POPPED=true to score the popped output instead. All workflow inputs are passed
# as scalars (sample list / bed list / labels are newline-delimited FILES read
# with read_lines) so the proven `wb` batch-input-CSV mechanism applies.
#
# Perimeter: SummarizeEvaluations installs pandas OFFLINE from a pip wheelhouse
# staged here; the vcfdist + bcftools images are pulled normally.
# =============================================================================
set -euo pipefail

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
AOU_USER_PROJECT="${AOU_USER_PROJECT:-${PROJECT}}"

RUN_ID="${RUN_ID:-holdout198_v2}"
WORK="${BUCKET}/imp_${RUN_ID}"
LOCAL="${HOME}/imp_${RUN_ID}"
VCFDIST_WORK="${WORK}/vcfdist"                    # gs:// staging for vcfdist inputs
REGION="${REGION:-chr1}"
OUT_BASE="${OUT_BASE:-aou_holdout198_${REGION}}"
HOLDOUT_OUTPUT_PATH="${HOLDOUT_OUTPUT_PATH:-imputebeagle-${RUN_ID}-run}"

# ---- eval VCF = the holdout imputed output (un-popped by default) ----
IMPUTED_VCF="${IMPUTED_VCF:-}"
POPPED="${POPPED:-false}"
EVAL_TAG=""; EVAL_PATHTAG=""
if [ "${POPPED}" = "true" ]; then EVAL_TAG=".popped"; EVAL_PATHTAG="-popped"; fi

# ---- truth VCF = the 198 full-panel genotypes (phased); id source too ----
TRUTH_AOU="${TRUTH_AOU:-${BUCKET}/vcf/truth.aou.${REGION}.bcf}"
TRUTH_AOU_IDX="${TRUTH_AOU_IDX:-${TRUTH_AOU}.csi}"
# if truth ids differ from eval ids (e.g. dipcall "syndip"), provide a file of truth ids (one/line)
TRUTH_SAMPLES_FILE="${TRUTH_SAMPLES_FILE:-}"

# ---- reference (vcfdist needs the FASTA + .fai; GRCh38, chr-named) ----
REF_FASTA="${REF_FASTA:-gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta}"
REF_FAI="${REF_FAI:-${REF_FASTA}.fai}"

# ---- analysis regions / stratifications ----
# CONF_BED: optional confident-regions BED (gs:// bed[.gz]); applied to BOTH sides at subset time.
CONF_BED="${CONF_BED:-}"
# VCFDIST_BEDS / VCFDIST_LABELS: parallel comma-lists of analysis BEDs + labels for vcfdist's -b.
# If unset, a single whole-${REGION} BED (built from the .fai) is staged with label "${REGION}".
VCFDIST_BEDS="${VCFDIST_BEDS:-}"
VCFDIST_LABELS="${VCFDIST_LABELS:-}"

DO_NAIVELY_PHASE="${DO_NAIVELY_PHASE:-false}"     # convert / to | in the EVAL only (Beagle output is already phased)
MAX_SAMPLES="${MAX_SAMPLES:-5}"                   # test cap: 5 samples. Set MAX_SAMPLES=0 for all 198.
# vcfdist memory has TWO independent knobs:
#   1) clustering: cap supercluster size (--cluster size <min_gap> + --max-supercluster-size) so a dense
#      popped region can't build a giant supercluster.
#   2) precision/recall: vcfdist's --max-ram defaults to 64 GB, so the P/R phase parallelizes to fill
#      ~64 GB (threads-per-tier = max-ram/tier). MUST set --max-ram <= the task RAM or it OOMs there
#      regardless of clustering. We set it to VCFDIST_MEM_GB/2 to leave headroom for the base (reference
#      FASTA + variants). Override VCFDIST_EXTRA_ARGS / VCFDIST_MEM_GB to tune the speed/RAM tradeoff.
VCFDIST_MEM_GB="${VCFDIST_MEM_GB:-32}"
VCFDIST_EXTRA_ARGS="${VCFDIST_EXTRA_ARGS:---cluster size 100 --max-supercluster-size 10000 --max-ram $((VCFDIST_MEM_GB / 2))}"

# ---- pandas wheelhouse for SummarizeEvaluations (perimeter blocks PyPI in-task) ----
PIP_WHEELHOUSE="${PIP_WHEELHOUSE:-}"              # gs:// to a prebuilt wheelhouse.tar.gz (skips the build)
PIP_PKGS="${PIP_PKGS:-pandas}"
PIP_PY_VER="${PIP_PY_VER:-311}"

# ---- registration ----
WF_VCFDIST="${WF_VCFDIST:-VcfdistEvaluation}"
WDL_REL="wdl/vcfdist_eval"
WDL_GCS="${BUCKET}/${WDL_REL}/"
FORCE_WF_REGISTER="${FORCE_WF_REGISTER:-false}"

if [ -z "${WDL_SRC_LOCAL:-}" ]; then
  for d in "${HOME}/vcfdist_eval_wdl" "./vcfdist_eval_wdl" "${PWD}/vcfdist_eval_wdl" "${HOME}"; do
    [ -f "${d}/${WF_VCFDIST}.wdl" ] && { WDL_SRC_LOCAL="$d"; break; }
  done
fi

# ============================ preflight ======================================
for t in gsutil gcloud wb bcftools tabix bgzip python3; do
  command -v "$t" >/dev/null 2>&1 || { echo "MISSING TOOL: $t"; exit 1; }
done
[ -n "${PROJECT}" ] || { echo "PROJECT empty (set GOOGLE_CLOUD_PROJECT)"; exit 1; }
mkdir -p "${LOCAL}"
gs_exists() { gsutil -q stat "$1" 2>/dev/null; }
GS_RETRIES="${GS_RETRIES:-5}"
retry() { local a=1; until "$@"; do [ "$a" -ge "${GS_RETRIES}" ] && { echo ">> gave up: $*" >&2; return 1; }; echo ">> retry ${a} (sleep $((a*10))s): $*" >&2; sleep "$((a*10))"; a=$((a+1)); done; }

# list sample names from a VCF/BCF header (local or gs://). The awk `exit` closes the pipe at
# #CHROM (SIGPIPE upstream); the subshell disables pipefail so the early close is not a failure.
vcf_samples() {
  local src="$1"
  ( set +o pipefail
    if [[ "${src}" == gs://* ]]; then gsutil cat "${src}"; else cat "${src}"; fi \
      | bcftools view -h /dev/stdin 2>/dev/null \
      | awk -F'\t' '/^#CHROM/{ for(i=10;i<=NF;i++) print $i; exit }'
  )
}

# ============================ resolve eval (imputed) output ==================
if [ -z "${IMPUTED_VCF}" ]; then
  ROOT="${BUCKET}/${HOLDOUT_OUTPUT_PATH}"
  if [ "${POPPED}" = "true" ]; then PAT="/${OUT_BASE}\.imputed\.popped\.vcf\.gz$"; else PAT="/${OUT_BASE}\.imputed\.vcf\.gz$"; fi
  echo ">> discovering eval output (${PAT}) under ${ROOT}/"
  IMPUTED_VCF="$(gsutil ls -r "${ROOT}/" 2>/dev/null | grep -E "${PAT}" | sort | tail -1 || true)"
  [ -n "${IMPUTED_VCF}" ] || { echo "ERROR: could not find ${PAT} under ${ROOT}/. Pass IMPUTED_VCF=gs://.../<file>."; exit 1; }
fi
IMPUTED_IDX="${IMPUTED_IDX:-${IMPUTED_VCF}.tbi}"
gs_exists "${IMPUTED_VCF}" || { echo "ERROR: eval VCF not found: ${IMPUTED_VCF}"; exit 1; }
gs_exists "${IMPUTED_IDX}" || { echo "ERROR: eval index not found: ${IMPUTED_IDX}"; exit 1; }
gs_exists "${TRUTH_AOU}"   || { echo "ERROR: truth VCF not found: ${TRUTH_AOU} (set TRUTH_AOU=)"; exit 1; }
gs_exists "${TRUTH_AOU_IDX}" || { echo ">> WARN: truth index ${TRUTH_AOU_IDX} not found; vcfdist subset will create one"; TRUTH_AOU_IDX=""; }
gs_exists "${REF_FASTA}"   || { echo "ERROR: reference FASTA not found: ${REF_FASTA} (set REF_FASTA=)"; exit 1; }
gs_exists "${REF_FAI}"     || { echo "ERROR: reference .fai not found: ${REF_FAI} (set REF_FAI=)"; exit 1; }
echo ">> eval (imputed): ${IMPUTED_VCF}"
echo ">> truth         : ${TRUTH_AOU}"
echo ">> reference     : ${REF_FASTA}"

# ============================ build + stage inputs ==========================
# sample names: the eval (imputed) carries exactly the 198 held-out ids.
SN_LOCAL="${LOCAL}/vcfdist.sample_names.txt"
vcf_samples "${IMPUTED_VCF}" > "${SN_LOCAL}"
NS_ALL="$(wc -l < "${SN_LOCAL}")"
[ "${NS_ALL}" -gt 0 ] || { echo "ERROR: read 0 sample names from ${IMPUTED_VCF}"; exit 1; }
if [ "${MAX_SAMPLES}" -gt 0 ] && [ "${MAX_SAMPLES}" -lt "${NS_ALL}" ]; then
  head -n "${MAX_SAMPLES}" "${SN_LOCAL}" > "${SN_LOCAL}.head" && mv "${SN_LOCAL}.head" "${SN_LOCAL}"
  echo ">> MAX_SAMPLES=${MAX_SAMPLES}: scoring ${MAX_SAMPLES} of ${NS_ALL} samples"
fi
NS="$(wc -l < "${SN_LOCAL}")"
# sanity: how many of these are present in the truth?
TRUTH_SN_LOCAL="${LOCAL}/vcfdist.truth_sample_names.txt"
vcf_samples "${TRUTH_AOU}" > "${TRUTH_SN_LOCAL}" || true
NMATCH="$(LC_ALL=C comm -12 <(LC_ALL=C sort -u "${SN_LOCAL}") <(LC_ALL=C sort -u "${TRUTH_SN_LOCAL}") | wc -l)"
echo ">> samples: ${NS} eval ids; ${NMATCH} also present in truth"
[ "${NMATCH}" -gt 0 ] || { echo "ERROR: none of the eval sample ids are in the truth (${TRUTH_AOU}); id namespace mismatch."; exit 1; }
[ "${NMATCH}" -eq "${NS}" ] || echo ">> WARN: ${NS} eval ids but only ${NMATCH} in truth; missing ids will fail their subset task."

SN_GCS="${VCFDIST_WORK}/sample_names.txt"
gsutil cp "${SN_LOCAL}" "${SN_GCS}"

# analysis BED list + labels
BEDLIST_LOCAL="${LOCAL}/vcfdist.bed_list.txt"; : > "${BEDLIST_LOCAL}"
LABELS_LOCAL="${LOCAL}/vcfdist.labels.txt";    : > "${LABELS_LOCAL}"
if [ -n "${VCFDIST_BEDS}" ]; then
  IFS=',' read -r -a _beds   <<< "${VCFDIST_BEDS}"
  IFS=',' read -r -a _labels <<< "${VCFDIST_LABELS:-}"
  for k in "${!_beds[@]}"; do
    gs_exists "${_beds[$k]}" || { echo "ERROR: VCFDIST_BEDS entry not found: ${_beds[$k]}"; exit 1; }
    echo "${_beds[$k]}" >> "${BEDLIST_LOCAL}"
    echo "${_labels[$k]:-strat${k}}" >> "${LABELS_LOCAL}"
  done
else
  # single whole-${REGION} BED built from the reference .fai
  FAI_L="${LOCAL}/$(basename "${REF_FAI}")"
  [ -s "${FAI_L}" ] || retry gsutil cp "${REF_FAI}" "${FAI_L}"
  REGION_BED_L="${LOCAL}/vcfdist.${REGION}.bed"
  awk -v c="${REGION}" 'BEGIN{OFS="\t"} $1==c{ print $1, 0, $2 }' "${FAI_L}" > "${REGION_BED_L}"
  [ -s "${REGION_BED_L}" ] || { echo "ERROR: ${REGION} not found in ${REF_FAI}; set VCFDIST_BEDS/REF_FAI."; exit 1; }
  REGION_BED_GCS="${VCFDIST_WORK}/${REGION}.bed"
  gsutil cp "${REGION_BED_L}" "${REGION_BED_GCS}"
  echo "${REGION_BED_GCS}" >> "${BEDLIST_LOCAL}"
  echo "${REGION}"         >> "${LABELS_LOCAL}"
fi
BEDLIST_GCS="${VCFDIST_WORK}/bed_list.txt"; gsutil cp "${BEDLIST_LOCAL}" "${BEDLIST_GCS}"
LABELS_GCS="${VCFDIST_WORK}/labels.txt";    gsutil cp "${LABELS_LOCAL}" "${LABELS_GCS}"
echo ">> stratifications: $(wc -l < "${BEDLIST_LOCAL}")  ($(paste -d= "${LABELS_LOCAL}" "${BEDLIST_LOCAL}" | tr '\n' ' '))"

# truth sample-names file (only if the user supplied differing truth ids)
TRUTH_SN_GCS=""
if [ -n "${TRUTH_SAMPLES_FILE}" ]; then
  if [[ "${TRUTH_SAMPLES_FILE}" == gs://* ]]; then TRUTH_SN_GCS="${TRUTH_SAMPLES_FILE}";
  else TRUTH_SN_GCS="${VCFDIST_WORK}/truth_sample_names.txt"; gsutil cp "${TRUTH_SAMPLES_FILE}" "${TRUTH_SN_GCS}"; fi
fi

# pandas wheelhouse (cp311 manylinux) for SummarizeEvaluations
ensure_pip_wheelhouse() {
  if [[ "${PIP_WHEELHOUSE}" == gs://* ]]; then
    gs_exists "${PIP_WHEELHOUSE}" || { echo "ERROR: PIP_WHEELHOUSE not found: ${PIP_WHEELHOUSE}"; exit 1; }
    echo ">> pip wheelhouse (user gs://): ${PIP_WHEELHOUSE}"; return
  fi
  local tag; tag="$(printf '%s' "${PIP_PKGS} cp${PIP_PY_VER}" | cksum | cut -d' ' -f1)"
  local gcs="${VCFDIST_WORK}/bin/pandas_wheelhouse.${tag}.tar.gz"
  if gs_exists "${gcs}"; then echo ">> pip wheelhouse already staged: ${gcs}"; PIP_WHEELHOUSE="${gcs}"; return; fi
  local whd="${LOCAL}/vcfdist_wheelhouse"; rm -rf "${whd}"; mkdir -p "${whd}"
  echo ">> building pip wheelhouse (cp${PIP_PY_VER} manylinux): ${PIP_PKGS}"
  python3 -m pip download --only-binary=:all: \
      --python-version "${PIP_PY_VER}" --implementation cp --abi "cp${PIP_PY_VER}" \
      --platform manylinux_2_17_x86_64 --platform manylinux2014_x86_64 \
      -d "${whd}" ${PIP_PKGS} \
    || { echo "ERROR: pip download failed (PyPI egress?). Supply PIP_WHEELHOUSE=gs://.../pandas_wheelhouse.tar.gz."; exit 1; }
  local tb="${LOCAL}/pandas_wheelhouse.tar.gz"; tar -czf "${tb}" -C "${whd}" .
  echo ">> staging pip wheelhouse ($(du -h "${tb}" | cut -f1), $(ls "${whd}" | wc -l) wheels) -> ${gcs}"
  gsutil cp "${tb}" "${gcs}"; PIP_WHEELHOUSE="${gcs}"
}
ensure_pip_wheelhouse

# ============================ register WDL ==================================
register_wf() {
  local name="$1" file="$2"
  if [ "${FORCE_WF_REGISTER}" = "true" ] || ! gs_exists "${WDL_GCS}${file}"; then
    [ -n "${WDL_SRC_LOCAL:-}" ] && [ -f "${WDL_SRC_LOCAL%/}/${file}" ] \
      || { echo "ERROR: ${file} not in the bucket and no local copy found; set WDL_SRC_LOCAL=<dir with ${file}>."; exit 1; }
    echo ">> uploading ${WDL_SRC_LOCAL%/}/${file} -> ${WDL_GCS}"
    gsutil cp "${WDL_SRC_LOCAL%/}/${file}" "${WDL_GCS}"
  fi
  if [ "${FORCE_WF_REGISTER}" = "true" ] && wb workflow list 2>/dev/null | grep -qw "${name}"; then
    echo ">> FORCE_WF_REGISTER: deleting stale ${name}"; wb workflow delete --workflow="${name}" --quiet 2>/dev/null || true
  fi
  if wb workflow list 2>/dev/null | grep -qw "${name}"; then
    echo ">> workflow ${name} already registered"
  else
    echo ">> registering ${name} -> ${WDL_GCS}${file}"
    wb workflow create --bucket-id="${BUCKET_ID}" --path="${WDL_REL}/${file}" --workflow="${name}"
  fi
}
register_wf "${WF_VCFDIST}" "${WF_VCFDIST}.wdl"

# ============================ run ===========================================
run_wf() {                                 # $1=workflow $2=output_path $3=inputs(csv key=val)
  local name="$1" opath="$2" inputs="$3"
  local log="${LOCAL}/wb_vcfdist_${name}.log"
  echo ">> launching ${name}: output-path=${opath}"
  local bdir="${LOCAL}/batch_vcfdist"; mkdir -p "${bdir}"
  local injson="${bdir}/${name}.inputs.json" incsv="${bdir}/${name}.inputs.csv" incols="${bdir}/${name}.columns.json"
  WF="${name}" python3 - "${inputs}" "${injson}" "${incsv}" "${incols}" <<'PY'
import json, csv, os, sys
wf = os.environ["WF"]
pairs = [kv.split("=", 1) for kv in sys.argv[1].split(",") if kv]
d = {f"{wf}.{k}": v for k, v in pairs}
json.dump(d, open(sys.argv[2], "w"), indent=2)
header = [k.split(".", 1)[1] for k in d]
with open(sys.argv[3], "w", newline="") as f:
    w = csv.writer(f); w.writerow(header); w.writerow(list(d.values()))
json.dump({k: k.split(".", 1)[1] for k in d}, open(sys.argv[4], "w"), indent=2)
PY
  local bgcs="${VCFDIST_WORK}/batch/"
  local csv_rel="${bgcs#${BUCKET}/}${name}.inputs.csv"
  local cols_uri="${bgcs}${name}.columns.json"
  gsutil cp "${injson}" "${incsv}" "${incols}" "${bgcs}"
  wb workflow job run \
    --workflow="${name}" \
    --output-bucket-id="${BUCKET_ID}" \
    --output-path="${opath}" \
    --batch-input-bucket-id="${BUCKET_ID}" \
    --batch-input-csv-path="${csv_rel}" \
    --column-mapping-uri="${cols_uri}" \
    --read-from-cache --write-to-cache 2>&1 | sanitize | tee "${log}" \
  || { echo ">> ${name} run failed; see ${log}"; return 1; }
  local jid; jid="$(grep -m1 'Job ID:' "${log}" 2>/dev/null | sed 's/.*Job ID:[[:space:]]*//' | tr -d '[:space:]' || true)"
  [ -n "$jid" ] && echo ">> ${name} job id: ${jid}"
}

# assemble inputs (all scalar: gs:// paths / strings / bools)
IN="eval_vcf=${IMPUTED_VCF},eval_vcf_idx=${IMPUTED_IDX}"
IN="${IN},truth_vcf=${TRUTH_AOU}"
[ -n "${TRUTH_AOU_IDX}" ] && IN="${IN},truth_vcf_idx=${TRUTH_AOU_IDX}"
IN="${IN},sample_names_file=${SN_GCS}"
[ -n "${TRUTH_SN_GCS}" ] && IN="${IN},truth_sample_names_file=${TRUTH_SN_GCS}"
[ -n "${CONF_BED}" ] && { gs_exists "${CONF_BED}" || { echo "ERROR: CONF_BED not found: ${CONF_BED}"; exit 1; }; IN="${IN},confident_regions_bed=${CONF_BED}"; }
IN="${IN},region=${REGION}"
IN="${IN},reference_fasta=${REF_FASTA},reference_fasta_fai=${REF_FAI}"
IN="${IN},do_naively_phase=${DO_NAIVELY_PHASE}"
IN="${IN},vcfdist_bed_list=${BEDLIST_GCS},labels_list=${LABELS_GCS}"
[ -n "${VCFDIST_EXTRA_ARGS}" ] && IN="${IN},vcfdist_extra_args=${VCFDIST_EXTRA_ARGS}"
IN="${IN},vcfdist_mem_gb=${VCFDIST_MEM_GB}"
IN="${IN},pip_wheelhouse=${PIP_WHEELHOUSE}"

run_wf "${WF_VCFDIST}" "vcfdist-holdout198${EVAL_PATHTAG}" "${IN}"

echo
echo "vcfdist out: ${BUCKET}/vcfdist-holdout198${EVAL_PATHTAG}/${WF_VCFDIST}/<uuid>/  (per-sample TSVs + evaluation_summary.tsv)"
echo "  phasing metrics: <sample>.phasing-summary.tsv / .switchflips.tsv / .phase-blocks.tsv"
echo "  aggregate      : evaluation_summary.tsv (PHASE_* + SNP/INDEL/SV precision/recall, averaged over samples)"

#!/usr/bin/env bash
# =============================================================================
# GLIMPSE2 eval prep — score the 198-AoU held-out imputation (chr1).
#
# Registers + runs the two eval WDLs against the holdout Beagle output
# (ImputeBeagleWithPop.imputed_multi_sample_vcf, un-popped = bubble.split):
#
#   Step C (PRIMARY)  GLIMPSE2Concordance : dosage r2 by AF bin, NRD, by length
#                     class x TR/homopolymer context, at 2 GP thresholds.
#                     panel_vcf = FULL panel (truth + freq, matched by name;
#                     intersects internally -> robust to site differences).
#
#   Step D (SANITY)   GLIMPSE2Summarize   : AF Pearson, HWE/de-Finetti, burden,
#                     ALT-length. Zips panel+imputed and HARD-FAILS on any
#                     CHROM/POS/REF/ALT mismatch -> panel MUST be subset to the
#                     imputed (variant-only) sites first (built here).
#
# Both eval WDLs are single-file and fetch tools at RUNTIME (Concordance wgets
# GLIMPSE2_concordance_static v2.0.1; Summarize conda-installs cyvcf2 etc.), so
# their tasks need network egress.
#
# Both inputs are OPTIONAL: TRH_BED (Step C) defaults to the GIAB TR+homopolymer
# stratification; POP_TSV (Step D) defaults to a single population ("ALL").
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

# Must match the imputation prep's RUN_ID so the eval reads the fresh run's outputs.
RUN_ID="${RUN_ID:-holdout198_v2}"
WORK="${BUCKET}/imp_${RUN_ID}"
LOCAL="${HOME}/imp_${RUN_ID}"
EVAL_WORK="${WORK}/eval"                          # gs:// staging for site-matched panel
REGION="${REGION:-chr1}"
OUT_PREFIX="${OUT_PREFIX:-aou_holdout198_chr1}"

# ---- the holdout imputed output (Step B) ----
# Pass IMPUTED_VCF explicitly, or let us discover it under the holdout run namespace.
# Default scores the UN-POPPED output (imputed_multi_sample_vcf, matches the bubble.split
# panel truth). Set POPPED=true later to score the popped output (imputed_popped_vcf).
IMPUTED_VCF="${IMPUTED_VCF:-}"
POPPED="${POPPED:-false}"
HOLDOUT_OUTPUT_PATH="${HOLDOUT_OUTPUT_PATH:-imputebeagle-${RUN_ID}-run}"
OUT_BASE="${OUT_BASE:-aou_holdout198_${REGION}}"  # the imputation output_basename -> <OUT_BASE>.imputed[.popped].vcf.gz

# POPPED=true scores the popped (constituent/atomic) output against the POPPED full panel; the default
# scores the un-popped (bubble.split) output against the bubble.split panel. EVAL_TAG namespaces every
# eval staging artifact + output path so the two evals never overwrite each other.
EVAL_TAG=""; EVAL_PATHTAG=""
if [ "${POPPED}" = "true" ]; then EVAL_TAG=".popped"; EVAL_PATHTAG="-popped"; fi
OUT_PREFIX="${OUT_PREFIX}${EVAL_TAG}"

# ---- truth/freq panel (FULL panel; bubble.split for un-popped, popped for POPPED) ----
FULL_SRC_TMPL="${FULL_SRC_TMPL-}"
[ -n "${FULL_SRC_TMPL}" ] || FULL_SRC_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.{contig}.bubble.split.bcf'
FULL_POPPED_SRC_TMPL="${FULL_POPPED_SRC_TMPL-}"
[ -n "${FULL_POPPED_SRC_TMPL}" ] || FULL_POPPED_SRC_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_popped_vcf/aou_lr_phase2_v1.{contig}.popped.bcf'
if [ "${POPPED}" = "true" ]; then
  PANEL_FULL="${PANEL_FULL:-${FULL_POPPED_SRC_TMPL//\{contig\}/${REGION}}}"
else
  PANEL_FULL="${PANEL_FULL:-${FULL_SRC_TMPL//\{contig\}/${REGION}}}"
fi
PANEL_FULL_IDX="${PANEL_FULL_IDX:-${PANEL_FULL}.csi}"

# ---- Step C: tandem-repeat / homopolymer stratification BED (OPTIONAL) ----
# TRH_BED may be a gs:// bed.gz (bgzipped+tabixed) OR an http(s) URL. If unset, the GIAB
# v3.6 all-TR+homopolymer (slop5) stratification is downloaded, indexed, and staged.
TRH_BED="${TRH_BED:-}"
TRH_BED_IDX="${TRH_BED_IDX:-${TRH_BED}.tbi}"
TRH_BED_URL="${TRH_BED_URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz}"

# ---- Step D open deps ----
ENABLE_SUMMARIZE="${ENABLE_SUMMARIZE:-true}"
POP_TSV="${POP_TSV:-}"                             # optional igsr-style TSV ("Sample name","Population code"). If unset, a one-population TSV (code "ALL") is auto-generated.
PANEL_SUMM_VCF="${PANEL_SUMM_VCF:-}"               # optional: pre-built site-matched panel (skips the heavy build)
PANEL_SUMM_IDX="${PANEL_SUMM_IDX:-${PANEL_SUMM_VCF}.tbi}"

# ---- registration ----
WF_CONC="${WF_CONC:-GLIMPSE2Concordance}"
WF_SUMM="${WF_SUMM:-GLIMPSE2Summarize}"
WDL_EVAL_REL="wdl/glimpse2_eval"                  # <<< registration dir for the eval WDLs
WDL_EVAL_GCS="${BUCKET}/${WDL_EVAL_REL}/"
FORCE_WF_REGISTER="${FORCE_WF_REGISTER:-false}"
ENABLE_CONCORDANCE="${ENABLE_CONCORDANCE:-true}"

# auto-detect the eval WDL bundle on this VM (the two single-file WDLs)
if [ -z "${WDL_EVAL_SRC_LOCAL:-}" ]; then
  for d in "${HOME}/glimpse2_eval_wdl" "./glimpse2_eval_wdl" "${PWD}/glimpse2_eval_wdl" "${HOME}"; do
    if [ -f "${d}/${WF_CONC}.wdl" ] && [ -f "${d}/${WF_SUMM}.wdl" ]; then WDL_EVAL_SRC_LOCAL="$d"; break; fi
  done
fi

# ============================ preflight ======================================
for t in gsutil gcloud wb bcftools tabix bgzip curl; do
  command -v "$t" >/dev/null 2>&1 || { echo "MISSING TOOL: $t"; exit 1; }
done
[ -n "${PROJECT}" ] || { echo "PROJECT empty (set GOOGLE_CLOUD_PROJECT)"; exit 1; }
mkdir -p "${LOCAL}"
gs_exists() { gsutil -q stat "$1" 2>/dev/null; }
GS_RETRIES="${GS_RETRIES:-5}"
retry() { local a=1; until "$@"; do [ "$a" -ge "${GS_RETRIES}" ] && { echo ">> gave up: $*" >&2; return 1; }; echo ">> retry ${a} (sleep $((a*10))s): $*" >&2; sleep "$((a*10))"; a=$((a+1)); done; }

# list sample names from a VCF/BCF (local or gs://), reading only the header. The awk `exit`
# closes the pipe at #CHROM (SIGPIPE upstream) so only the header is transferred; the subshell
# disables pipefail so that early close is not treated as a pipeline failure.
vcf_samples() {
  local f="$1"
  ( set +o pipefail
    if [[ "$f" == gs://* ]]; then
      case "$f" in
        *.bcf) gsutil cat "$f" 2>/dev/null | bcftools view -h - 2>/dev/null ;;
        *)     gsutil cat "$f" 2>/dev/null | zcat 2>/dev/null ;;
      esac | awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) print $i; exit}'
    else
      bcftools query -l "$f"
    fi
  )
}

# Resolve TRH_BED (+idx) to gs:// paths the WDL can localize. Accepts a gs:// bed.gz (with idx),
# an http(s) URL, or empty (-> TRH_BED_URL default). URLs are downloaded, ensured bgzipped +
# coordinate-sorted + tabixed, and staged to the bucket (idempotent).
ensure_trh_bed() {
  if [[ "${TRH_BED}" == gs://* ]]; then
    gs_exists "${TRH_BED}"     || { echo "ERROR: TRH_BED not found: ${TRH_BED}"; exit 1; }
    gs_exists "${TRH_BED_IDX}" || { echo "ERROR: TRH_BED_IDX not found: ${TRH_BED_IDX} (need a tabix .tbi for the bed)."; exit 1; }
    echo ">> TRH bed (user gs://): ${TRH_BED}"; return
  fi
  local url="${TRH_BED}"; [[ "${url}" == http*://* ]] || url="${TRH_BED_URL}"
  local base; base="$(basename "${url}")"
  local gcs="${EVAL_WORK}/trh/${base}"
  if gs_exists "${gcs}" && gs_exists "${gcs}.tbi"; then
    echo ">> TRH bed already staged: ${gcs}"
  else
    echo ">> staging TRH bed from ${url}"
    mkdir -p "${LOCAL}/trh"
    local L="${LOCAL}/trh/${base}"
    [ -s "${L}" ] || curl -fL --retry 6 --retry-all-errors --retry-delay 10 -o "${L}" "${url}" \
      || { echo "ERROR: cannot download TRH bed from ${url} (or supply your own: TRH_BED=gs://...bed.gz TRH_BED_IDX=gs://...bed.gz.tbi)."; exit 1; }
    if ! tabix -f -p bed "${L}" 2>/dev/null; then
      echo ">> re-bgzipping + coordinate-sorting TRH bed"
      local S="${LOCAL}/trh/sorted.${base}"
      zcat -f "${L}" | grep -v '^#' | LC_ALL=C sort -k1,1 -k2,2n | bgzip > "${S}"
      mv "${S}" "${L}"; tabix -f -p bed "${L}"
    fi
    gsutil cp "${L}" "${gcs}"; gsutil cp "${L}.tbi" "${gcs}.tbi"
  fi
  TRH_BED="${gcs}"; TRH_BED_IDX="${gcs}.tbi"
}

# ============================ resolve imputed output =========================
if [ -z "${IMPUTED_VCF}" ]; then
  ROOT="${BUCKET}/${HOLDOUT_OUTPUT_PATH}"
  if [ "${POPPED}" = "true" ]; then
    PAT="/${OUT_BASE}\.imputed\.popped\.vcf\.gz$"; echo ">> discovering POPPED output under ${ROOT}/"
  else
    PAT="/${OUT_BASE}\.imputed\.vcf\.gz$";        echo ">> discovering un-popped output under ${ROOT}/"
  fi
  IMPUTED_VCF="$(gsutil ls -r "${ROOT}/" 2>/dev/null | grep -E "${PAT}" | sort | tail -1 || true)"
  [ -n "${IMPUTED_VCF}" ] || { echo "ERROR: could not find a match for ${PAT} under ${ROOT}/."; echo "       pass it explicitly: IMPUTED_VCF=gs://.../<file> bash $0"; exit 1; }
fi
IMPUTED_IDX="${IMPUTED_IDX:-${IMPUTED_VCF}.tbi}"
gs_exists "${IMPUTED_VCF}" || { echo "ERROR: imputed VCF not found: ${IMPUTED_VCF}"; exit 1; }
gs_exists "${IMPUTED_IDX}" || { echo "ERROR: imputed index not found: ${IMPUTED_IDX}"; exit 1; }
echo ">> imputed : ${IMPUTED_VCF}"
echo ">> panel   : ${PANEL_FULL}"

# ============================ register eval WDLs =============================
register_eval_wf() {                       # $1=workflow name, $2=wdl filename
  local name="$1" file="$2"
  gs_exists "${WDL_EVAL_GCS}${file}" || {
    [ -n "${WDL_EVAL_SRC_LOCAL:-}" ] && [ -f "${WDL_EVAL_SRC_LOCAL%/}/${file}" ] \
      || { echo "ERROR: ${file} not in this workspace bucket and no local copy found."; echo "       set WDL_EVAL_SRC_LOCAL=<dir with ${WF_CONC}.wdl + ${WF_SUMM}.wdl>."; exit 1; }
    echo ">> uploading ${WDL_EVAL_SRC_LOCAL%/}/${file} -> ${WDL_EVAL_GCS}"
    gsutil cp "${WDL_EVAL_SRC_LOCAL%/}/${file}" "${WDL_EVAL_GCS}"
  }
  if [ "${FORCE_WF_REGISTER}" = "true" ] && wb workflow list 2>/dev/null | grep -qw "${name}"; then
    echo ">> FORCE_WF_REGISTER: deleting stale ${name}"; wb workflow delete --workflow="${name}" --quiet 2>/dev/null || true
  fi
  if wb workflow list 2>/dev/null | grep -qw "${name}"; then
    echo ">> workflow ${name} already registered"
  else
    echo ">> registering ${name} -> ${WDL_EVAL_GCS}${file}"
    wb workflow create --bucket-id="${BUCKET_ID}" --path="${WDL_EVAL_REL}/${file}" --workflow="${name}"
  fi
}
[ "${ENABLE_CONCORDANCE}" = "true" ] && register_eval_wf "${WF_CONC}" "${WF_CONC}.wdl"
[ "${ENABLE_SUMMARIZE}"   = "true" ] && register_eval_wf "${WF_SUMM}" "${WF_SUMM}.wdl"

run_wf() {                                 # $1=workflow $2=output_path $3=inputs(csv key=val)
  local name="$1" opath="$2" inputs="$3"
  local log="${LOCAL}/wb_eval_${name}.log"     # separate decl: a single `local` expands all RHS (incl ${name}) before assigning, tripping set -u
  echo ">> launching ${name}: output-path=${opath}"

  # `wb workflow job run --inputs` expects JSON, not key=val CSV. Use the same proven mechanism as
  # the imputation prep: build a fully-qualified inputs JSON -> 1-row batch CSV (bare columns) +
  # column-mapping JSON, stage them, and run with --batch-input-csv-path / --column-mapping-uri.
  local bdir="${LOCAL}/batch_eval"; mkdir -p "${bdir}"
  local injson="${bdir}/${name}.inputs.json"
  local incsv="${bdir}/${name}.inputs.csv"
  local incols="${bdir}/${name}.columns.json"
  WF="${name}" python3 - "${inputs}" "${injson}" "${incsv}" "${incols}" <<'PY'
import json, csv, os, sys
wf = os.environ["WF"]
pairs = [kv.split("=", 1) for kv in sys.argv[1].split(",") if kv]
d = {f"{wf}.{k}": v for k, v in pairs}                       # fully-qualified inputs (all string-valued)
json.dump(d, open(sys.argv[2], "w"), indent=2)
header = [k.split(".", 1)[1] for k in d]                     # bare column names
with open(sys.argv[3], "w", newline="") as f:
    w = csv.writer(f); w.writerow(header); w.writerow(list(d.values()))
json.dump({k: k.split(".", 1)[1] for k in d}, open(sys.argv[4], "w"), indent=2)  # qualified -> bare
PY

  local bgcs="${EVAL_WORK}/batch/"
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

# ============================ Step C: Concordance ============================
if [ "${ENABLE_CONCORDANCE}" = "true" ]; then
  ensure_trh_bed
  C_IN="panel_vcf=${PANEL_FULL},panel_vcf_idx=${PANEL_FULL_IDX}"
  C_IN="${C_IN},imputed_vcf=${IMPUTED_VCF},imputed_vcf_idx=${IMPUTED_IDX}"
  C_IN="${C_IN},trh_bed=${TRH_BED},trh_bed_idx=${TRH_BED_IDX}"
  C_IN="${C_IN},region=${REGION},output_prefix=${OUT_PREFIX}"
  run_wf "${WF_CONC}" "glimpse2-concordance-holdout198${EVAL_PATHTAG}" "${C_IN}"
fi

# ============================ Step D: Summarize ==============================
if [ "${ENABLE_SUMMARIZE}" = "true" ]; then
  if [ -n "${POP_TSV}" ]; then
    gs_exists "${POP_TSV}" || { echo "ERROR: POP_TSV not found: ${POP_TSV}"; exit 1; }
  fi

  # Summarize zips panel+imputed by record and hard-fails on any site mismatch. The imputed
  # output is variant-records-only (all-hom-ref-in-198 sites dropped), so it is a SUBSET of
  # the full-panel sites. Build a site-matched panel = full panel restricted to the imputed
  # sites (allele-aware), same position order, all panel samples retained.
  if [ -z "${PANEL_SUMM_VCF}" ]; then
    PANEL_SUMM_VCF="${EVAL_WORK}/panel.${REGION}${EVAL_TAG}.site_matched.vcf.gz"
    PANEL_SUMM_IDX="${PANEL_SUMM_VCF}.tbi"
    if gs_exists "${PANEL_SUMM_VCF}" && gs_exists "${PANEL_SUMM_IDX}"; then
      echo ">> site-matched panel exists: ${PANEL_SUMM_VCF}"
    else
      echo ">> building site-matched panel (full panel restricted to imputed sites)"
      IMP_L="${LOCAL}/imputed.${REGION}${EVAL_TAG}.vcf.gz"
      [ -s "${IMP_L}" ]      || retry gsutil cp "${IMPUTED_VCF}" "${IMP_L}"
      [ -s "${IMP_L}.tbi" ]  || retry gsutil cp "${IMPUTED_IDX}" "${IMP_L}.tbi"
      SITES="${LOCAL}/imp.sites.${REGION}${EVAL_TAG}.tsv.gz"
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${IMP_L}" | bgzip > "${SITES}"
      tabix -s1 -b2 -e2 "${SITES}"
      echo ">> imputed sites: $(zcat "${SITES}" | wc -l)"
      FULL_L="${LOCAL}/full_panel.${REGION}${EVAL_TAG}.bcf"
      if [ ! -s "${FULL_L}" ]; then
        retry gsutil -u "${AOU_USER_PROJECT}" cp "${PANEL_FULL}" "${FULL_L}.part"
        mv "${FULL_L}.part" "${FULL_L}"
      fi
      PS_L="${LOCAL}/panel.${REGION}${EVAL_TAG}.site_matched.vcf.gz"
      bcftools view -T "${SITES}" -m2 -M2 "${FULL_L}" -Oz -o "${PS_L}"
      tabix -p vcf "${PS_L}"
      # guard: record counts should match (imputed subset of panel sites)
      NI="$(bcftools index -n "${IMP_L}")"; NP="$(bcftools index -n "${PS_L}")"
      echo ">> site-matched panel records: ${NP}  (imputed: ${NI})"
      [ "${NP}" = "${NI}" ] || echo ">> WARNING: panel(${NP}) != imputed(${NI}) record count; Summarize zip may fail. Inspect before relying on Step D."
      gsutil cp "${PS_L}" "${PANEL_SUMM_VCF}"; gsutil cp "${PS_L}.tbi" "${PANEL_SUMM_IDX}"
    fi
  fi
  gs_exists "${PANEL_SUMM_VCF}" || { echo "ERROR: site-matched panel missing: ${PANEL_SUMM_VCF}"; exit 1; }

  # population_tsv: if not supplied, treat ALL samples as one population ("ALL"). The WDL input
  # is mandatory (it merges sample -> Population code), so we synthesize a one-pop map covering
  # every panel + imputed sample. NB: the per-population boxplots key color/order to the 1000G
  # code set, so a single "ALL" code leaves those boxplots empty; the AF-Pearson, AF-hist2d,
  # HWE/de-Finetti and ALT-length outputs do not use the map and are unaffected.
  if [ -z "${POP_TSV}" ]; then
    POP_TSV="${EVAL_WORK}/pop.all.${REGION}.tsv"
    if gs_exists "${POP_TSV}"; then
      echo ">> one-pop TSV exists: ${POP_TSV}"
    else
      echo ">> generating one-population TSV (all samples -> 'ALL')"
      PT_L="${LOCAL}/pop.all.${REGION}.tsv"
      { vcf_samples "${PANEL_SUMM_VCF}"; vcf_samples "${IMPUTED_VCF}"; } \
        | LC_ALL=C sort -u | sed '/^$/d' \
        | awk 'BEGIN{FS=OFS="\t"; print "Sample name","Population code"} {print $1,"ALL"}' > "${PT_L}"
      echo ">> one-pop TSV samples: $(( $(wc -l < "${PT_L}") - 1 ))"
      [ "$(( $(wc -l < "${PT_L}") - 1 ))" -gt 0 ] || { echo "ERROR: could not read any sample names for the one-pop TSV (check ${PANEL_SUMM_VCF} / ${IMPUTED_VCF})."; exit 1; }
      gsutil cp "${PT_L}" "${POP_TSV}"
    fi
  fi

  S_IN="panel_vcf=${PANEL_SUMM_VCF},panel_vcf_idx=${PANEL_SUMM_IDX}"
  S_IN="${S_IN},imputed_vcf=${IMPUTED_VCF},imputed_vcf_idx=${IMPUTED_IDX}"
  S_IN="${S_IN},population_tsv=${POP_TSV},output_prefix=${OUT_PREFIX}"
  run_wf "${WF_SUMM}" "glimpse2-summarize-holdout198${EVAL_PATHTAG}" "${S_IN}"
fi

echo "================= DONE ================="
echo "eval WDLs registered from: ${WDL_EVAL_GCS}"
echo "  ${WF_CONC} -> ${WDL_EVAL_REL}/${WF_CONC}.wdl"
echo "  ${WF_SUMM} -> ${WDL_EVAL_REL}/${WF_SUMM}.wdl"
[ "${ENABLE_CONCORDANCE}" = "true" ] && echo "Concordance out: ${BUCKET}/glimpse2-concordance-holdout198${EVAL_PATHTAG}/${WF_CONC}/<uuid>/  (r2/nrd PNGs + rsquare/error tables)"
[ "${ENABLE_SUMMARIZE}"   = "true" ] && echo "Summarize out:   ${BUCKET}/glimpse2-summarize-holdout198${EVAL_PATHTAG}/${WF_SUMM}/<uuid>/  (pearson.tsv + AF/HWE/altlen PDFs)"

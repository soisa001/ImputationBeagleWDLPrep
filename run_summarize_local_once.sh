#!/usr/bin/env bash
# =============================================================================
# One-off LOCAL re-run of GLIMPSE2 Summarize for the failed popped run.
#
# The Cromwell SummarizeAndPlot task completed the (CHROM,POS,REF,ALT) merge-join
# (2,899,122 matched, 0 target-only) but then crashed in the population boxplot
# (seaborn 0.13.2 + matplotlib 3.11 'boxprops' bug) -- BEFORE any output was
# written/delocalized, so nothing persisted. This runs the SAME Summarize script
# locally (on the notebook, which has PyPI egress) in a venv with matplotlib<3.10,
# then uploads the results. It re-streams once (~1h20m) -- unavoidable, the task's
# in-memory arrays are gone -- but it is a single local step, no Cromwell/wb.
#
# Inputs default to the exact paths from the failed run; override via env if needed.
#   bash run_summarize_local_once.sh
# =============================================================================
set -euo pipefail

BUCKET_GS="${BUCKET_GS:-gs://cloned-rw-migration-aou-rw-f178dfde-wb-sharp-papaya-7463}"
PANEL_SUMM_VCF="${PANEL_SUMM_VCF:-${BUCKET_GS}/imp_holdout198_v2/eval/panel.chr1.popped.site_matched.vcf.gz}"
IMPUTED_VCF="${IMPUTED_VCF:-${BUCKET_GS}/imputebeagle-holdout198_v2-run/ImputeBeagleWithPop_holdout198_v2_job_2026-06-26T085204336474180Z/ImputeBeagleWithPop/d66bb401-a0e1-4489-bf58-92932d2e8cd2/call-PopAndMarginalizeCollisions/attempt-2/aou_holdout198_chr1.imputed.popped.vcf.gz}"
POP_TSV="${POP_TSV:-${BUCKET_GS}/imp_holdout198_v2/eval/pop.all.chr1.tsv}"
OUT_PREFIX="${OUT_PREFIX:-aou_holdout198_chr1.popped}"
OUT_GCS="${OUT_GCS:-$(dirname "${POP_TSV}")/summarize_manual}"   # where to upload pearson.tsv + *.pdf

WORK="${WORK:-${HOME}/summarize_local_once}"
VENV="${VENV:-${WORK}/venv}"
PY_PKGS="${PY_PKGS:-cyvcf2 pandas numpy matplotlib<3.10 seaborn scipy}"

# locate the Summarize WDL (source of the exact script to run)
WDL="${WDL:-}"
if [ -z "${WDL}" ]; then
  for d in "./glimpse2_eval_wdl" "${PWD}/glimpse2_eval_wdl" "${HOME}/glimpse2_eval_wdl" "$(dirname "$0")/glimpse2_eval_wdl"; do
    [ -f "${d}/GLIMPSE2Summarize.wdl" ] && { WDL="${d}/GLIMPSE2Summarize.wdl"; break; }
  done
fi
[ -n "${WDL}" ] && [ -f "${WDL}" ] || { echo "ERROR: GLIMPSE2Summarize.wdl not found (run from the repo root or set WDL=)"; exit 1; }

for t in gsutil python3; do command -v "$t" >/dev/null || { echo "MISSING TOOL: $t"; exit 1; }; done
mkdir -p "${WORK}"; cd "${WORK}"
echo ">> work dir: ${WORK}"
echo ">> panel  : ${PANEL_SUMM_VCF}"
echo ">> imputed: ${IMPUTED_VCF}"
echo ">> pop tsv: ${POP_TSV}"
echo ">> prefix : ${OUT_PREFIX}"

# ---- venv with the FIXED pins (matplotlib<3.10) ----
if [ ! -x "${VENV}/bin/python" ]; then
  echo ">> creating venv ${VENV} with: ${PY_PKGS}"
  python3 -m venv "${VENV}"
  "${VENV}/bin/pip" install -q --upgrade pip
  # shellcheck disable=SC2086
  "${VENV}/bin/pip" install -q ${PY_PKGS}
else
  echo ">> reusing venv ${VENV}"
fi
"${VENV}/bin/python" - <<'PY'
import matplotlib, seaborn, cyvcf2, pandas, numpy, scipy
assert tuple(int(x) for x in matplotlib.__version__.split('.')[:2]) < (3,10), matplotlib.__version__
print(f">> deps OK: matplotlib {matplotlib.__version__} | seaborn {seaborn.__version__} | cyvcf2 {cyvcf2.__version__}")
PY

# ---- download inputs (cyvcf2 needs local files; iterates sequentially, no index needed) ----
dl() { local src="$1" dst="$2"; [ -s "${dst}" ] && { echo "   (have ${dst})"; return; }; echo ">> downloading ${src}"; gsutil -q cp "${src}" "${dst}"; }
PANEL_L="${WORK}/$(basename "${PANEL_SUMM_VCF}")"
IMP_L="${WORK}/$(basename "${IMPUTED_VCF}")"
TSV_L="${WORK}/$(basename "${POP_TSV}")"
dl "${PANEL_SUMM_VCF}" "${PANEL_L}"
dl "${IMPUTED_VCF}"    "${IMP_L}"
dl "${POP_TSV}"        "${TSV_L}"

# ---- extract the EXACT Summarize python from the WDL (the merge-join version) ----
python3 - "${WDL}" "${WORK}/summarize.py" <<'PY'
import re, sys
src = open(sys.argv[1]).read()
m = re.search(r"\n(        import sys\n.*?)\n        EOF", src, re.S)
assert m, "could not extract the embedded Summarize python from the WDL"
body = "\n".join(ln[8:] if ln.startswith("        ") else ln for ln in m.group(1).splitlines())
open(sys.argv[2], "w").write(body)
print(">> extracted summarize.py")
PY

# ---- run it (same arg order as the WDL: panel, imputed, population_tsv, output_prefix) ----
echo ">> running Summarize locally (this re-streams ${PANEL_SUMM_VCF##*/} vs the imputed; ~1h20m)..."
MPLBACKEND=Agg "${VENV}/bin/python" "${WORK}/summarize.py" \
    "${PANEL_L}" "${IMP_L}" "${TSV_L}" "${OUT_PREFIX}"

echo ">> outputs written locally in ${WORK}:"
ls -la "${OUT_PREFIX}.pearson.tsv" ./*.pdf 2>/dev/null || true

# ---- upload ----
if [ -n "${OUT_GCS}" ]; then
  echo ">> uploading to ${OUT_GCS}/"
  gsutil -q cp "${OUT_PREFIX}.pearson.tsv" "${OUT_GCS}/" || true
  gsutil -q cp ./*.pdf "${OUT_GCS}/" || true
  echo ">> done. Results: ${OUT_GCS}/"
else
  echo ">> OUT_GCS empty -> results kept local only (${WORK})"
fi

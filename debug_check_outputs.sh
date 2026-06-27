#!/usr/bin/env bash
# debug_check_outputs.sh -- coherence + representation checks for the holdout198 imputation outputs.
#
# Verifies that:
#   POPPED output  (aou_holdout198_<contig>.imputed.popped.vcf.gz)
#     - is actually POPPED (atomic constituent representation: biallelic, single-token INFO/ID,
#       FORMAT=GT:DS:GP, DS == GP-implied dosage);
#     - every (CHROM,POS,REF,ALT) exists in the POPPED full panel (panel_popped_vcf) and in the
#       id-split (atomic) panel -- i.e. the popped representation matches the popped panel;
#     - reproduces the GLIMPSE2Summarize zip() alignment against the site-matched panel and reports
#       the first desync (the "VCFs are out of sync" failure mode).
#   UN-POPPED output (aou_holdout198_<contig>.imputed.vcf.gz)
#     - is the bubble.split representation (biallelic, FORMAT=GT:DS:GP);
#     - every (CHROM,POS,REF,ALT) exists in the bubble.split sites-only LEAVEOUT panel (the bref3 the
#       posteriors come from) -- i.e. the un-popped representation matches the un-popped panel; and
#       reports overlap with the FULL bubble.split panel (the un-popped eval truth).
#
# Read-only: pulls nothing into the project, downloads only the (small) imputed outputs locally and
# STREAMS panel sites (no full panel copy). Panels are controlled-access; set USER_PROJECT for billing.
#
# Usage:
#   bash debug_check_outputs.sh                 # discover both outputs under the run namespace, check
#   POPPED_OUT=gs://.../x.popped.vcf.gz UNPOPPED_OUT=gs://.../x.vcf.gz bash debug_check_outputs.sh
#
set -euo pipefail

# ----------------------------- config (mirrors the prep scripts) -----------------------------
RUN_ID="${RUN_ID:-holdout198_v2}"
REGION="${REGION:-chr1}"
BUCKET_ID="${BUCKET_ID:-cloned-rw-migration-aou-rw-f178dfde-wb-sharp-papaya-7463}"
_B="${BUCKET:-}"
case "${_B}" in
  gs://*) BUCKET="${_B}" ;;
  "")     BUCKET="gs://${BUCKET_ID}" ;;
  *)      BUCKET="gs://${_B}" ;;
esac
HOLDOUT_OUTPUT_PATH="${HOLDOUT_OUTPUT_PATH:-imputebeagle-${RUN_ID}-run}"
OUT_BASE="${OUT_BASE:-aou_holdout198_${REGION}}"
USER_PROJECT="${USER_PROJECT:-${AOU_USER_PROJECT:-}}"   # billing project for requester-pays panel reads

# Panels (same templates/defaults as prep_imputebeagle_holdout198.sh / prep_eval_holdout198.sh).
# NOTE: assign brace-containing defaults via a separate single-quoted statement -- a bare '}' inside a
# ${VAR:-default} default value prematurely closes the expansion and mangles the {contig} placeholder.
POPPED_PANEL_TMPL="${FULL_POPPED_SRC_TMPL:-}"
[ -n "${POPPED_PANEL_TMPL}" ] || POPPED_PANEL_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_popped_vcf/aou_lr_phase2_v1.{contig}.popped.bcf'
IDSPLIT_PANEL_TMPL="${ID_SPLIT_TMPL:-}"
[ -n "${IDSPLIT_PANEL_TMPL}" ] || IDSPLIT_PANEL_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_id_split_vcf_gz/aou_lr_phase2_v1.{contig}.id.split.vcf.gz'
BUBBLE_LEAVEOUT_TMPL="${SITES_ONLY_TMPL:-}"
[ -n "${BUBBLE_LEAVEOUT_TMPL}" ] || BUBBLE_LEAVEOUT_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_sites_only_leaveout_vcf/aou_lr_phase2_v1.{contig}.bubble.split.sites.leaveout.bcf'
BUBBLE_FULL_TMPL="${FULL_SRC_TMPL:-}"
[ -n "${BUBBLE_FULL_TMPL}" ] || BUBBLE_FULL_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.{contig}.bubble.split.bcf'

POPPED_PANEL="${POPPED_PANEL_TMPL//\{contig\}/${REGION}}"
IDSPLIT_PANEL="${IDSPLIT_PANEL_TMPL//\{contig\}/${REGION}}"
BUBBLE_LEAVEOUT="${BUBBLE_LEAVEOUT_TMPL//\{contig\}/${REGION}}"
BUBBLE_FULL="${BUBBLE_FULL_TMPL//\{contig\}/${REGION}}"

WORK="${WORK:-${TMPDIR:-/tmp}/holdout198_debug}"
mkdir -p "${WORK}"
SAMPLE_N="${SAMPLE_N:-200000}"   # how many records to spot-check for DS/GP coherence (0 = all)

fail=0; warn=0
red(){ printf '\033[31m%s\033[0m\n' "$*"; }
grn(){ printf '\033[32m%s\033[0m\n' "$*"; }
ylw(){ printf '\033[33m%s\033[0m\n' "$*"; }
PASS(){ grn   "  PASS  $*"; }
WARN(){ ylw   "  WARN  $*"; warn=$((warn+1)); }
FAIL(){ red   "  FAIL  $*"; fail=$((fail+1)); }
hdr(){ printf '\n=== %s ===\n' "$*"; }

retry(){ local n=0; until "$@"; do n=$((n+1)); [ "$n" -ge 4 ] && return 1; sleep $((2**n)); done; }

# stream a (possibly requester-pays, possibly local) VCF/BCF body through bcftools query of CHROM/POS/REF/ALT
gsargs(){ if [ -n "${USER_PROJECT}" ]; then echo "-u ${USER_PROJECT}"; fi; }
panel_keys(){  # $1 = gs:// or local src ; writes sorted-unique CHROM\tPOS\tREF\tALT to stdout-cache $2
  local src="$1" out="$2"
  if [ -s "${out}" ]; then echo "  (cached: ${out} -> $(wc -l < "${out}") keys)" >&2; return; fi
  echo "  extracting sites from ${src} ..." >&2
  if [[ "${src}" == gs://* ]]; then
    # shellcheck disable=SC2046
    gsutil $(gsargs) cat "${src}" | bcftools view -G -Ou - | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -
  else
    bcftools view -G -Ou "${src}" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -
  fi | LC_ALL=C sort -u > "${out}"
  echo "  -> $(wc -l < "${out}") keys" >&2
}

# count keys in A not present in B (A,B are sorted-unique key files); echoes the count, writes misses to $3
missing_keys(){ LC_ALL=C comm -23 "$1" "$2" | tee "${3:-/dev/null}" | wc -l; }

# ----------------------------- resolve the two imputed outputs -----------------------------
discover(){ # $1 = regex pattern under the run dir
  gsutil ls -r "${BUCKET}/${HOLDOUT_OUTPUT_PATH}/" 2>/dev/null | grep -E "$1" | sort | tail -1 || true
}
POPPED_OUT="${POPPED_OUT:-}"
UNPOPPED_OUT="${UNPOPPED_OUT:-}"
[ -n "${POPPED_OUT}"   ] || POPPED_OUT="$(discover "/${OUT_BASE}\.imputed\.popped\.vcf\.gz$")"
[ -n "${UNPOPPED_OUT}" ] || UNPOPPED_OUT="$(discover "/${OUT_BASE}\.imputed\.vcf\.gz$")"

echo "run namespace : ${BUCKET}/${HOLDOUT_OUTPUT_PATH}"
echo "region        : ${REGION}"
echo "popped out    : ${POPPED_OUT:-<not found>}"
echo "un-popped out : ${UNPOPPED_OUT:-<not found>}"
echo "popped panel  : ${POPPED_PANEL}"
echo "id-split panel: ${IDSPLIT_PANEL}"
echo "bubble leaveout: ${BUBBLE_LEAVEOUT}"
echo "bubble full    : ${BUBBLE_FULL}"
echo "workdir       : ${WORK}"

dl(){ # $1 src $2 dst  (small imputed outputs only)
  [ -s "$2" ] && return 0
  retry gsutil $(gsargs) cp "$1" "$2"
}

# =========================================================================================
#                                   POPPED OUTPUT
# =========================================================================================
if [ -n "${POPPED_OUT}" ]; then
  hdr "POPPED output coherence: ${POPPED_OUT##*/}"
  PL="${WORK}/popped.${REGION}.vcf.gz"
  dl "${POPPED_OUT}" "${PL}"; dl "${POPPED_OUT}.tbi" "${PL}.tbi" || retry bcftools index -t "${PL}"
  NREC="$(bcftools index -n "${PL}")"
  echo "  records: ${NREC}"

  # C1: FORMAT must be GT:DS:GP
  FMTS="$( ( set +o pipefail; bcftools view -H "${PL}" | cut -f9 | LC_ALL=C sort -u ) )"
  if [ "${FMTS}" = "GT:DS:GP" ]; then PASS "FORMAT == GT:DS:GP"; else FAIL "FORMAT not uniformly GT:DS:GP -> {${FMTS//$'\n'/,}}"; fi

  # C2: biallelic (atomic) -- no multiallelic ALT
  NMULTI="$( ( set +o pipefail; bcftools query -f '%ALT\n' "${PL}" | grep -c ',' ) || true )"
  if [ "${NMULTI:-0}" -eq 0 ]; then PASS "all records biallelic (popped/atomic)"; else FAIL "${NMULTI} records have multiallelic ALT (not popped)"; fi

  # C3: INFO/ID is a single atomic id -- a bubble id (comma/colon-joined) means NOT popped
  NBUBBLE="$( ( set +o pipefail; bcftools query -f '%ID\n' "${PL}" | grep -c '[,:]' ) || true )"
  if [ "${NBUBBLE:-0}" -eq 0 ]; then PASS "all INFO/ID single-token (no residual bubble ids)"; else FAIL "${NBUBBLE} records still carry a multi-id bubble ID (not popped)"; fi

  # C4: DS == GP-implied dosage (gp1 + 2*gp2), within rounding (|.|<=0.01). SAMPLE_N = #genotype lines.
  echo "  checking DS vs GP dosage (first ${SAMPLE_N} genotype rows)..."
  BADDS="$( ( set +o pipefail
    bcftools query -f '[%DS\t%GP\n]' "${PL}" | { [ "${SAMPLE_N}" -gt 0 ] && head -n "${SAMPLE_N}" || cat; }
  ) | awk -F'\t' '
      { n=split($2,g,","); if(n<3) next;
        ds=$1+0; dose=g[2]+2*g[3];
        d=ds-dose; if(d<0)d=-d;
        tot=g[1]+g[2]+g[3];
        if(d>0.01) bad++;
        if(tot<0.97||tot>1.03) badsum++;
        c++ }
      END{ printf "%d %d %d", bad+0, badsum+0, c+0 }' )"
  read -r b_ds b_sum n_ds <<<"${BADDS}"
  if [ "${n_ds:-0}" -gt 0 ] && [ "${b_ds:-0}" -eq 0 ]; then PASS "DS == gp1+2*gp2 for all ${n_ds} sampled genotypes"; \
    elif [ "${n_ds:-0}" -eq 0 ]; then WARN "DS/GP spot-check produced no rows (query path?)"; \
    else FAIL "${b_ds}/${n_ds} sampled genotypes have DS != gp1+2*gp2 (>0.01)"; fi
  [ "${b_sum:-0}" -eq 0 ] || WARN "${b_sum}/${n_ds} sampled GP vectors do not sum to ~1"

  # keys of the popped output
  POUT_KEYS="${WORK}/popped_out.keys"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${PL}" | LC_ALL=C sort -u > "${POUT_KEYS}"
  echo "  popped output unique keys: $(wc -l < "${POUT_KEYS}")"

  # C5: every popped key is in the POPPED panel
  PP_KEYS="${WORK}/popped_panel.keys"; panel_keys "${POPPED_PANEL}" "${PP_KEYS}"
  M="$(missing_keys "${POUT_KEYS}" "${PP_KEYS}" "${WORK}/popped_not_in_poppedpanel.keys")"
  if [ "${M}" -eq 0 ]; then PASS "all popped keys present in panel_popped_vcf (representation matches popped panel)"; \
    else FAIL "${M} popped keys NOT in panel_popped_vcf (see ${WORK}/popped_not_in_poppedpanel.keys)"; fi

  # C6: every popped key is in the ID-SPLIT (atomic) panel
  ID_KEYS="${WORK}/idsplit_panel.keys"; panel_keys "${IDSPLIT_PANEL}" "${ID_KEYS}"
  M2="$(missing_keys "${POUT_KEYS}" "${ID_KEYS}" "${WORK}/popped_not_in_idsplit.keys")"
  if [ "${M2}" -eq 0 ]; then PASS "all popped keys present in id-split panel (atomic source)"; \
    else WARN "${M2} popped keys NOT in id-split panel (see ${WORK}/popped_not_in_idsplit.keys)"; fi

  # C8: check the GLIMPSE2Summarize (CHROM,POS,REF,ALT) merge-join against the site-matched panel.
  #     The eval builds the panel via `bcftools view -T <imputed sites> -m2 -M2` (position- not
  #     allele-aware), so at a multiallelic position the panel keeps MORE alleles than the imputed
  #     carries -> the panel record count is legitimately larger. The Summarize merge-join pairs by
  #     exact (REF,ALT), so what matters is: (a) every imputed key is present in the site-matched
  #     panel (-> all imputed sites get scored), and (b) how many panel-only alleles get skipped.
  hdr "POPPED Summarize-alignment check (merge-join on CHROM,POS,REF,ALT)"
  SITES_T="${WORK}/popped_out.sites.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${PL}" | bgzip > "${SITES_T}.gz"
  tabix -s1 -b2 -e2 "${SITES_T}.gz"
  # panel restricted to the imputed sites (what the eval's -T produces), as sorted-unique keys
  SM_KEYS="${WORK}/site_matched_panel.keys"
  if [[ "${POPPED_PANEL}" == gs://* ]]; then
    # shellcheck disable=SC2046
    gsutil $(gsargs) cat "${POPPED_PANEL}" \
      | bcftools view -T "${SITES_T}.gz" -m2 -M2 - 2>/dev/null \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' - | LC_ALL=C sort -u > "${SM_KEYS}" || true
  else
    bcftools view -T "${SITES_T}.gz" -m2 -M2 "${POPPED_PANEL}" \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' - | LC_ALL=C sort -u > "${SM_KEYS}" || true
  fi
  NSMP="$(wc -l < "${SM_KEYS}")"; NIMP="$(wc -l < "${POUT_KEYS}")"
  UNMATCHED="$(missing_keys "${POUT_KEYS}" "${SM_KEYS}" "${WORK}/imputed_not_in_sitematched.keys")"
  MATCHED=$(( NIMP - UNMATCHED )); PANELONLY=$(( NSMP - MATCHED ))
  echo "  site-matched panel unique keys: ${NSMP}   imputed unique keys: ${NIMP}"
  echo "  merge-join: matched=${MATCHED}  imputed-only(skipped)=${UNMATCHED}  panel-only(skipped)=${PANELONLY}"
  if [ "${UNMATCHED}" -eq 0 ]; then
    PASS "merge-join aligns ALL ${MATCHED} imputed sites; ${PANELONLY} panel-only alleles skipped (expected: -T keeps extra multiallelic alleles). Summarize will NOT desync."
  else
    WARN "${UNMATCHED} imputed keys absent from the site-matched panel (merge-join skips them as target-only; see ${WORK}/imputed_not_in_sitematched.keys)"
  fi
else
  WARN "no popped output found/given -- skipping popped checks"
fi

# =========================================================================================
#                                   UN-POPPED OUTPUT
# =========================================================================================
if [ -n "${UNPOPPED_OUT}" ]; then
  hdr "UN-POPPED output coherence: ${UNPOPPED_OUT##*/}"
  UL="${WORK}/unpopped.${REGION}.vcf.gz"
  dl "${UNPOPPED_OUT}" "${UL}"; dl "${UNPOPPED_OUT}.tbi" "${UL}.tbi" || retry bcftools index -t "${UL}"
  NRECU="$(bcftools index -n "${UL}")"
  echo "  records: ${NRECU}"

  # U1: FORMAT GT:DS:GP
  FMTSU="$( ( set +o pipefail; bcftools view -H "${UL}" | cut -f9 | LC_ALL=C sort -u ) )"
  if [ "${FMTSU}" = "GT:DS:GP" ]; then PASS "FORMAT == GT:DS:GP"; else FAIL "FORMAT not uniformly GT:DS:GP -> {${FMTSU//$'\n'/,}}"; fi

  # U2: biallelic (bubble.split is split to biallelic)
  NMU="$( ( set +o pipefail; bcftools query -f '%ALT\n' "${UL}" | grep -c ',' ) || true )"
  if [ "${NMU:-0}" -eq 0 ]; then PASS "all records biallelic (bubble.split)"; else FAIL "${NMU} multiallelic ALT records (expected split)"; fi

  UOUT_KEYS="${WORK}/unpopped_out.keys"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${UL}" | LC_ALL=C sort -u > "${UOUT_KEYS}"
  echo "  un-popped output unique keys: $(wc -l < "${UOUT_KEYS}")"

  # U3: every un-popped key is in the bubble.split sites-only LEAVEOUT panel (the bref3 source)
  BL_KEYS="${WORK}/bubble_leaveout.keys"; panel_keys "${BUBBLE_LEAVEOUT}" "${BL_KEYS}"
  MU="$(missing_keys "${UOUT_KEYS}" "${BL_KEYS}" "${WORK}/unpopped_not_in_leaveout.keys")"
  if [ "${MU}" -eq 0 ]; then PASS "all un-popped keys present in bubble.split LEAVEOUT panel (representation matches the un-popped panel)"; \
    else FAIL "${MU} un-popped keys NOT in the bubble.split leaveout panel (see ${WORK}/unpopped_not_in_leaveout.keys)"; fi

  # U4: overlap with the FULL bubble.split panel (the un-popped eval truth) -- informational
  BF_KEYS="${WORK}/bubble_full.keys"; panel_keys "${BUBBLE_FULL}" "${BF_KEYS}"
  MF="$(missing_keys "${UOUT_KEYS}" "${BF_KEYS}" "${WORK}/unpopped_not_in_full.keys")"
  NOV=$(( $(wc -l < "${UOUT_KEYS}") - MF ))
  if [ "${MF}" -eq 0 ]; then PASS "all un-popped keys present in the FULL bubble.split panel (eval truth)"; \
    else WARN "${MF} un-popped keys NOT in the FULL bubble.split panel (${NOV} overlap); the un-popped eval scores only the overlap"; fi
else
  WARN "no un-popped output found/given -- skipping un-popped checks"
fi

# =========================================================================================
hdr "SUMMARY"
echo "  failures: ${fail}   warnings: ${warn}"
if [ "${fail}" -eq 0 ]; then grn "ALL COHERENCE/REPRESENTATION CHECKS PASSED"; else red "${fail} CHECK(S) FAILED -- see above"; fi
exit "$(( fail > 0 ? 1 : 0 ))"

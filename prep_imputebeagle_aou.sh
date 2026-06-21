#!/usr/bin/env bash
# =============================================================================
# ImputationBeagle — AoU OUT-OF-PANEL imputation (chr21 demo), start to finish.
#
# Reference  : AoU lrWGS Phase 2 v1 panel (~10k, SNV+INDEL+SV), built into bref3.
# Target     : 1000 AoU srWGS (ACAF) samples that are NOT in the panel.
# Evaluation : no per-sample truth for out-of-panel samples, so instead of
#              concordance we run SANITY CHECKS on the imputed SVs vs the panel:
#                - ALT allele frequency (imputed 1000 vs panel)
#                - HWE p-value distribution (imputed vs panel)
#                - Beagle DR2 quality distribution
#                - AF bias vs panel MAF
#
# Pipeline:
#   1. panel -> refsrc (INFO strip + norm -d exact), unique_variants,
#      unique_variants, bref3, genetic map        [host: bcftools + java, streamed]
#   2. fixed sample list: first N AoU samples (from an ACAF shard header) not in
#      the panel sample set                        [host: header parse]
#   3. discover chr21 shards (tabix -l binary search over the mount) + fixed sample
#      list, then build the target LOCALLY: subset each chr21 shard in parallel over
#      the mount (bcftools view -r, biallelic SNP, sample subset; parts cached/atomic),
#      concat + index, upload to the bucket           [host: multithreaded bcftools]
#   4. inputs CSV + column-map JSON, upload WDLs, register, run (cache on, chunk=1000)
#   5. wait for the imputed VCF, run the sanity checks
#
# Everything runs on the host: panel prep streams the gs:// panel through bcftools
# (+ java for bref3); the shard/target work reads the gcsfuse-mounted ACAF bucket
# (byte-range reads, samples-subset), all multithreaded. No dsub.
# =============================================================================
set -euo pipefail

# --- quiet, log-friendly wrappers (avoid spinner/ANSI/NUL noise in nohup logs) ---
gsutil() { command gsutil -q "$@"; }                 # -q: drop progress bars/info (keeps errors + stdout)
wb()     { TERM=dumb NO_COLOR=1 command wb "$@"; }    # force plain output (no JLine spinner/ANSI)
# strip ANSI/OSC escapes, NUL, and CR-redraws from a captured stream
sanitize() { LC_ALL=C sed -E $'s/\x1b\\[[0-9;?]*[ -/]*[@-~]//g; s/\x1b\\][^\x07]*(\x07|\x1b\\\\)//g' \
             | LC_ALL=C tr -d '\000' | LC_ALL=C tr '\r' '\n'; }

# ============================ CONFIG =========================================
# BUCKET_ID is the workspace bucket *resource* name (see `wb resource list`). In a cloned
# workspace it physically resolves to gs://cloned-<id>-<project>, NOT gs://<id>. gsutil writes
# must target that same physical bucket or wb/the workflow engine won't find what we upload.
BUCKET_ID="${BUCKET_ID:-rw-migration-aou-rw-f178dfde}"   # resource name (NOT necessarily the gs:// name)
if [ -z "${BUCKET:-}" ]; then
  _B="$(wb resource resolve --name "${BUCKET_ID}" 2>/dev/null | tr -d '[:space:]')"
  case "${_B}" in
    gs://*) BUCKET="${_B}" ;;                            # resolved physical URL (clone-aware)
    "")     BUCKET="gs://${BUCKET_ID}" ;;                # resolve unavailable -> assume literal
    *)      BUCKET="gs://${_B}" ;;                       # resolved bare name
  esac
fi
WORK="${BUCKET}/imp_phase1"
LOCAL="${HOME}/imp_phase1"

# AoU lrWGS Phase 2 v1 10k panel, bubble-decomposed + split BCF. gs:// -> used in place.
# Default on its own line: a {contig} brace inside ${VAR:-default} would truncate it.
PANEL_SRC_TMPL="${PANEL_SRC_TMPL-}"
[ -n "${PANEL_SRC_TMPL}" ] || PANEL_SRC_TMPL='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.{contig}.bubble.split.bcf'

CONTIGS=(${CONTIGS:-chr21})                       # chr21 only for this demo

REF_PREFIX="${WORK}/ref/aou_lr_phase2_v1_bubble"  # bref3 = <prefix>.<contig>.bref3 (distinct from the popped run)
MAPS_DIR="${WORK}/ref/genetic_maps/"              # trailing slash REQUIRED by the WDL
TARGET_DIR="${WORK}/target/"                      # gs:// dir (trailing slash)
SAMP_DIR="${WORK}/samples/"                       # gs:// dir (trailing slash)
MAPS_ZIP_URL="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"

# ---- AoU out-of-panel target ----
AOU_VCF_MOUNT="${AOU_VCF_MOUNT:-${HOME}/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/vcf}"
AOU_VCF_GS="${AOU_VCF_GS:-gs://vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/vcf}"
N_TARGET="${N_TARGET:-1000}"                       # # samples to impute
# Candidate samples: one research_id per line, header 'research_id'. The target is
# the first N_TARGET of these that are AVAILABLE in the ACAF callset (and, if
# EXCLUDE_PANEL, not in the panel). gs:// is downloaded; local path used as-is.
TARGET_TSV="${TARGET_TSV:-selected_samples_100k_95_sexmatched_SNOMED.tsv}"
EXCLUDE_PANEL="${EXCLUDE_PANEL:-true}"             # keep the out-of-panel premise (drop any TSV id in the panel)
TARGET_BASE="${TARGET_BASE:-aou_oop}"             # out-of-panel target basename
OUT_BASE="${OUT_BASE:-aou_oop_${N_TARGET}_${CONTIGS[0]}}"   # workflow output_basename


# ---- workflow registration + run (Workbench `wb`) ----
RUN_WORKFLOW="${RUN_WORKFLOW:-true}"              # false = stage inputs only
WF_NAME="${WF_NAME:-imputation-beagle}"
FORCE_WF_REGISTER="${FORCE_WF_REGISTER:-false}"   # true = delete+recreate the workflow so its WDL path matches BUCKET
RUN_TAG="${RUN_TAG:-aou1000bub}"                  # tag baked into the output namespace (distinct from the popped run)
OUTPUT_PATH="${OUTPUT_PATH:-imputebeagle-${RUN_TAG}-run}"   # runs land at gs://<bucket>/<OUTPUT_PATH>/ImputationBeagle/<uuid>/
USE_BATCH_CSV="${USE_BATCH_CSV:-true}"            # file-based inputs (CSV + column-mapping JSON)
SAMPLE_CHUNK_SIZE="${SAMPLE_CHUNK_SIZE:-1000}"    # WDL's designed sample-chunk size
PHASE_MEM_GB="${PHASE_MEM_GB:-60}"                # beagle_phase_memory_in_gb  (WDL default 40)
IMPUTE_MEM_GB="${IMPUTE_MEM_GB:-65}"              # beagle_impute_memory_in_gb (WDL default 45)
ERROR_COUNT_OVERRIDE="${ERROR_COUNT_OVERRIDE:-0}" # 0 = bypass chunk-QC fail gate (impute all chunks); set empty to enforce QC

# ---- sanity eval ----
RUN_EVAL="${RUN_EVAL:-true}"
EVAL_POLL_SECS="${EVAL_POLL_SECS:-120}"
EVAL_MAX_POLLS="${EVAL_MAX_POLLS:-180}"           # 180 x 120s = 6h

# WDL bundle: the canonical WDLs live in another (controlled-tier) workspace bucket. We copy them
# into THIS workspace's bucket so `wb workflow create` (which reads --path within --bucket-id =
# this workspace) can register them. All other inputs are read cross-workspace in place; only the
# WDLs must be local to register.
WDL_SRC_GCS="${WDL_SRC_GCS:-gs://longreadsphase2imputation/wdl/imputation_beagle/}"   # source (cross-workspace, read)
WDL_REL="wdl/imputation_beagle"
WDL_GCS="${BUCKET}/${WDL_REL}/"                   # local workspace dest (registration reads here)

# clean-slate controls
IDEMPOTENT="${IDEMPOTENT:-true}"  # keep all prior artifacts (cleaned panel/refsrc, bref3, maps) AND the run namespace; reuse via gs_exists
CLEAN="${CLEAN:-false}"           # delete stale panel-derived inputs for the CONTIGS (ignored when IDEMPOTENT=true)
CLEAN_RUNS="${CLEAN_RUNS:-false}" # delete this run's prior output namespace (ignored when IDEMPOTENT=true)
if [ "${IDEMPOTENT}" = "true" ]; then CLEAN=false; CLEAN_RUNS=false; fi

# ---- local host execution (no dsub) ----
PROJECT="${PROJECT:-${GOOGLE_CLOUD_PROJECT:-$(gcloud config get-value project 2>/dev/null || true)}}"
AOU_USER_PROJECT="${AOU_USER_PROJECT:-${PROJECT}}"   # billing project for requester-pays AoU bucket (empty-mount probe)
BREF3_JAR_URL="${BREF3_JAR_URL:-https://faculty.washington.edu/browning/beagle/bref3.17Dec24.224.jar}"
THREADS="${THREADS:-16}"                             # local VM parallelism budget (this VM = 16 vCPU)
PANEL_THREADS="${PANEL_THREADS:-${THREADS}}"         # bcftools threads for panel norm + sort --parallel
# target build runs LOCALLY over the gcsfuse mount: PARALLEL shards x SHARD_THREADS each (~= THREADS)
SHARD_THREADS="${SHARD_THREADS:-2}"                  # bcftools --threads per shard subset (compression)
PARALLEL="${PARALLEL:-$(( THREADS/SHARD_THREADS > 0 ? THREADS/SHARD_THREADS : 1 ))}"   # concurrent shard subsets
MERGE_THREADS="${MERGE_THREADS:-${THREADS}}"         # bcftools concat/index threads

# ============================ preflight ======================================
for t in gsutil gcloud curl unzip awk wb python3 bcftools tabix java; do
  command -v "$t" >/dev/null 2>&1 || { echo "MISSING TOOL: $t"; exit 1; }
done
python3 -c "import cyvcf2,numpy,pandas,scipy,matplotlib" 2>/dev/null \
  || { echo ">> installing python deps for sanity eval"; pip install --quiet cyvcf2 numpy pandas scipy matplotlib; }
[ -n "${PROJECT}" ] || { echo "PROJECT empty (set GOOGLE_CLOUD_PROJECT)"; exit 1; }
[ -d "${AOU_VCF_MOUNT}" ] || { echo "ERROR: ACAF mount not found: ${AOU_VCF_MOUNT}"; echo "       set AOU_VCF_MOUNT to the gcsfuse vcf/ dir."; exit 1; }
mkdir -p "${LOCAL}"
gs_exists() { gsutil -q stat "$1" 2>/dev/null; }

# Robust GCS access: large reads (panel BCFs, ACAF shards over gcsfuse) can be reset mid-stream
# ("Connection reset by peer"). retry() re-runs a command with linear backoff; ensure_panel_local()
# downloads the panel once with a resumable, retried cp so prep/bref3/target read a local file.
GS_RETRIES="${GS_RETRIES:-5}"
retry() {                                              # retry <cmd...> : up to GS_RETRIES, linear backoff
  local a=1
  until "$@"; do
    [ "$a" -ge "${GS_RETRIES}" ] && { echo ">> gave up after ${a} tries: $*" >&2; return 1; }
    echo ">> retry ${a}/${GS_RETRIES} (sleep $((a*10))s): $*" >&2
    sleep "$((a*10))"; a=$((a+1))
  done
}
ensure_panel_local() {                                 # download PANEL_CANON -> PANEL_LOCAL once (lazy, resumable)
  [ -s "${PANEL_LOCAL}" ] && return 0
  echo ">> downloading panel locally (resumable cp): ${PANEL_CANON} -> ${PANEL_LOCAL}" >&2
  retry gsutil -u "${AOU_USER_PROJECT}" cp "${PANEL_CANON}" "${PANEL_LOCAL}.part"
  mv "${PANEL_LOCAL}.part" "${PANEL_LOCAL}"             # atomic: final name appears only when complete
}
# ============================ helper scripts (written once) ===================
mkdir -p "${LOCAL}"
SANITY_PY="${LOCAL}/sanity_imputation.py"


cat > "${SANITY_PY}" <<'PYEOF'
# Sanity checks for out-of-panel imputed SVs (no per-sample truth): compare the
# imputed cohort's SV ALT-AF / HWE / DR2 against the reference panel.
import os, subprocess, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from cyvcf2 import VCF
try: from scipy.stats import chi2
except Exception: chi2 = None

IMP   = os.environ["IMPUTED_VCF"]
PANEL = os.environ["PANEL_VCF"]                    # refsrc (uniquified panel) for AF/HWE
CONTIG= os.environ.get("CONTIG", "chr21")
OUT_BASE = os.environ.get("OUT_BASE", "aou_oop")
RESDIR = os.path.expanduser(os.environ.get("RESDIR", "~/imp_sanity"))
FS = 26
os.makedirs(RESDIR, exist_ok=True)

def loc(p):
    if p.startswith("gs://"):
        lp = os.path.join(RESDIR, os.path.basename(p))
        if not os.path.exists(lp): subprocess.run(["gsutil", "cp", p, lp], check=True)
        return lp
    return p
def is_sv(ref, alt): return alt.startswith("<") or max(len(ref), len(alt)) >= 50
def svtype(alt): return alt.strip("<>").split(":")[0] if alt.startswith("<") else "OTHER"
def gtm(v): return np.array([r[:2] for r in v.genotypes], dtype=int)
def sv_alleles(v): return [(j, a) for j, a in enumerate(v.ALT, start=1) if is_sv(v.REF, a)]
def af_hwe(G, j):
    a, b = G[:, 0], G[:, 1]; miss = (a < 0) | (b < 0)
    aa, bb = a[~miss], b[~miss]; N = len(aa)
    if N < 1: return (np.nan,) * 5
    alt = (aa == j).astype(int) + (bb == j).astype(int)
    n0 = int(np.sum(alt == 0)); n1 = int(np.sum(alt == 1)); n2 = int(np.sum(alt == 2))
    if (n0 + n1 + n2) < 1: return (np.nan,) * 5
    p = (2 * n2 + n1) / (2 * (n0 + n1 + n2))
    hwe = np.nan
    if chi2 is not None and 0 < p < 1:
        e0 = (1 - p) ** 2 * N; e1 = 2 * p * (1 - p) * N; e2 = p ** 2 * N
        if min(e0, e1, e2) > 0:
            x2 = (n0 - e0) ** 2 / e0 + (n1 - e1) ** 2 / e1 + (n2 - e2) ** 2 / e2
            hwe = float(chi2.sf(x2, 1))
    return p, n0, n1, n2, hwe

P = VCF(loc(PANEL)); pan = {}
for v in P:
    if v.CHROM != CONTIG: continue
    G = gtm(v)
    for j, a in sv_alleles(v):
        af, _, _, _, hwe = af_hwe(G, j)
        if not np.isnan(af): pan[f"{v.CHROM}:{v.POS}:{v.REF}:{a}"] = (af, hwe, svtype(a), a)
print(f"panel SVs: {len(pan)}", flush=True)

I = VCF(loc(IMP)); rows = []; matched = 0
for v in I:
    if v.CHROM != CONTIG: continue
    svs = sv_alleles(v)
    if not svs: continue
    G = gtm(v); dr2 = v.INFO.get("DR2")
    for j, a in svs:
        k = f"{v.CHROM}:{v.POS}:{v.REF}:{a}"
        af, n0, n1, n2, hwe = af_hwe(G, j)
        if np.isnan(af): continue
        col = j - 1
        if dr2 is None: d = np.nan
        elif isinstance(dr2, (tuple, list, np.ndarray)): d = float(dr2[col]) if col < len(dr2) else np.nan
        else: d = float(dr2)
        paf, phwe, svt, alt = pan.get(k, (np.nan, np.nan, svtype(a), a))
        if k in pan: matched += 1
        rows.append((v.CHROM, int(v.POS), k, alt, svt, af, paf, d, hwe, phwe, n0, n1, n2))
df = pd.DataFrame(rows, columns=["chrom","pos","key","alt","svtype","imp_af","panel_af",
                                 "dr2","imp_hwe_p","panel_hwe_p","n_homref","n_het","n_homalt"])
df.to_csv(f"{RESDIR}/{OUT_BASE}.sv_sanity.tsv", sep="\t", index=False)
print(f"imputed SVs: {len(df)} | matched to panel: {matched} ({matched}/{len(pan)} panel SVs)", flush=True)
if len(df) == 0:
    print("No imputed SVs found.", flush=True); raise SystemExit(0)

m = df.dropna(subset=["imp_af", "panel_af"])
r = float(np.corrcoef(m.imp_af, m.panel_af)[0, 1]) if len(m) > 1 else np.nan
print(f"\n=== {OUT_BASE} sanity ===")
print(f"AF corr imputed-vs-panel: r={r:.4f} r2={r*r:.4f} | mean|dAF|={np.mean(np.abs(m.imp_af-m.panel_af)):.4f} (n={len(m)})")
if df.dr2.notna().any():
    print(f"DR2: mean={df.dr2.mean():.4f} median={df.dr2.median():.4f} frac>=0.8={np.mean(df.dr2.dropna()>=0.8):.3f}")
hh = df.imp_hwe_p.dropna()
if len(hh): print(f"imputed HWE: frac p<1e-6 = {np.mean(hh<1e-6):.4f} (n={len(hh)})")

# ---- DR2 (and AF concordance) stratified by panel-MAF bin ----
# bin by panel AF where the SV matched the panel, else by the imputed AF; report per bin.
df["_afsrc"] = df.panel_af.where(df.panel_af.notna(), df.imp_af)
df["_maf"]   = np.minimum(df["_afsrc"], 1 - df["_afsrc"])
df["_absdaf"] = (df.imp_af - df.panel_af).abs()
maf_bins = [0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5]
maf_lab  = ["<0.5%", "0.5-1%", "1-5%", "5-10%", "10-20%", "20-50%"]
df["_mafbin"] = pd.cut(df["_maf"].clip(upper=0.5), bins=maf_bins, labels=maf_lab, include_lowest=True)
gb = df.groupby("_mafbin", observed=False)
binned = pd.DataFrame({
    "n":              gb.size(),
    "mean_dr2":       gb["dr2"].mean(),
    "median_dr2":     gb["dr2"].median(),
    "frac_dr2_ge0.8": gb["dr2"].apply(lambda s: float((s.dropna() >= 0.8).mean()) if s.notna().any() else np.nan),
    "mean_abs_dAF":   gb["_absdaf"].mean(),
})
binned.to_csv(f"{RESDIR}/{OUT_BASE}.dr2_by_maf.tsv", sep="\t")
print("\nDR2 / AF-concordance by panel-MAF bin (panel AF where matched, else imputed):")
print(binned.to_string(float_format=lambda x: f"{x:.4f}"))

def nlp(s): s = s.dropna().clip(lower=1e-300); return -np.log10(s)
fig, ax = plt.subplots(2, 3, figsize=(33, 18))
ax[0,0].scatter(m.panel_af, m.imp_af, s=14, alpha=0.35, color="#534AB7")
ax[0,0].plot([0,1],[0,1], "--", lw=2, color="#A32D2D")
ax[0,0].set_title(f"SV ALT AF: imputed vs panel (r={r:.3f})", fontsize=FS+4)
ax[0,0].set_xlabel("Panel ALT AF", fontsize=FS); ax[0,0].set_ylabel("Imputed ALT AF", fontsize=FS)
ax[0,0].set_xlim(0,1); ax[0,0].set_ylim(0,1); ax[0,0].tick_params(labelsize=FS-2)
if df.dr2.notna().any(): ax[0,1].hist(df.dr2.dropna(), bins=30, color="#0F6E56")
ax[0,1].set_title("Imputed SV Beagle DR2", fontsize=FS+4)
ax[0,1].set_xlabel("DR2", fontsize=FS); ax[0,1].set_ylabel("SV count", fontsize=FS); ax[0,1].tick_params(labelsize=FS-2)
ax[1,0].hist(nlp(df.imp_hwe_p), bins=40, alpha=0.6, label="imputed", color="#C25E00")
if df.panel_hwe_p.notna().any(): ax[1,0].hist(nlp(df.panel_hwe_p), bins=40, alpha=0.6, label="panel", color="#444444")
ax[1,0].axvline(6, ls="--", color="red"); ax[1,0].legend(fontsize=FS-4)
ax[1,0].set_title("SV HWE -log10(p)", fontsize=FS+4)
ax[1,0].set_xlabel("-log10 HWE p", fontsize=FS); ax[1,0].set_ylabel("SV count", fontsize=FS); ax[1,0].tick_params(labelsize=FS-2)
m2 = m.copy(); m2["pmaf"] = np.minimum(m2.panel_af, 1 - m2.panel_af); m2["dAF"] = m2.imp_af - m2.panel_af
ax[1,1].scatter(m2.pmaf, m2.dAF, s=14, alpha=0.35, color="#1A6FB0")
ax[1,1].axhline(0, ls="--", color="#A32D2D")
ax[1,1].set_title("AF bias vs panel MAF", fontsize=FS+4)
ax[1,1].set_xlabel("Panel MAF", fontsize=FS); ax[1,1].set_ylabel("Imputed - panel AF", fontsize=FS); ax[1,1].tick_params(labelsize=FS-2)
# new: DR2 stratified by panel-MAF bin (mean DR2, and frac DR2>=0.8)
xb = np.arange(len(binned.index)); xl = list(binned.index.astype(str))
ax[0,2].bar(xb, binned["mean_dr2"].values, color="#0F6E56")
ax[0,2].set_xticks(xb); ax[0,2].set_xticklabels(xl, rotation=30, ha="right")
ax[0,2].set_ylim(0,1); ax[0,2].set_title("Mean DR2 by panel-MAF bin", fontsize=FS+4)
ax[0,2].set_xlabel("Panel MAF bin", fontsize=FS); ax[0,2].set_ylabel("Mean DR2", fontsize=FS); ax[0,2].tick_params(labelsize=FS-2)
for xx, nn in zip(xb, binned["n"].values):
    if nn > 0: ax[0,2].text(xx, 0.02, f"n={int(nn)}", ha="center", va="bottom", fontsize=FS-10, rotation=90)
ax[1,2].bar(xb, binned["frac_dr2_ge0.8"].values, color="#534AB7")
ax[1,2].set_xticks(xb); ax[1,2].set_xticklabels(xl, rotation=30, ha="right")
ax[1,2].set_ylim(0,1); ax[1,2].set_title("Frac DR2>=0.8 by panel-MAF bin", fontsize=FS+4)
ax[1,2].set_xlabel("Panel MAF bin", fontsize=FS); ax[1,2].set_ylabel("Frac DR2>=0.8", fontsize=FS); ax[1,2].tick_params(labelsize=FS-2)
plt.tight_layout(); fig.savefig(f"{RESDIR}/{OUT_BASE}.sanity.png", dpi=150, bbox_inches="tight")
print(f"\nwrote {RESDIR}/{OUT_BASE}.sv_sanity.tsv\nwrote {RESDIR}/{OUT_BASE}.dr2_by_maf.tsv\nwrote {RESDIR}/{OUT_BASE}.sanity.png", flush=True)
PYEOF

# ============================ clean slate ====================================
if [ "${IDEMPOTENT}" = "true" ]; then
  echo ">> IDEMPOTENT=true: keeping prior panel/refsrc/bref3/maps and the run namespace; reusing whatever already exists (set IDEMPOTENT=false to allow CLEAN/CLEAN_RUNS)"
fi
if [ "${CLEAN}" = "true" ]; then
  echo ">> CLEAN: deleting stale panel-derived inputs for: ${CONTIGS[*]}"
  for c in "${CONTIGS[@]}"; do
    gsutil -m rm -f \
      "${REF_PREFIX}.${c}.bref3" \
      "${REF_PREFIX}.${c}.unique_variants" \
      "${WORK}/ref/refsrc_noinfo.${c}.vcf.gz" \
      "${WORK}/ref/refsrc_noinfo.${c}.vcf.gz.tbi" 2>/dev/null || true
  done
fi
if [ "${CLEAN_RUNS}" = "true" ]; then
  echo ">> CLEAN: deleting this run's namespace ${BUCKET}/${OUTPUT_PATH}/ (make sure none are running)"
  gsutil -m rm -rf "${BUCKET}/${OUTPUT_PATH}/" 2>/dev/null || true
fi

# ============================ per-contig functions ===========================
PANEL_CANON=""; PANEL_LOCAL=""    # set by resolve_panel

resolve_panel() {                          # $1=contig -> sets PANEL_CANON (gs://) + PANEL_LOCAL (path; downloaded lazily)
  local c="$1" src; src="${PANEL_SRC_TMPL//\{contig\}/$c}"
  [[ "$src" == gs://* ]] || { echo "ERROR: AoU panel must be gs://; got: $src"; exit 1; }
  gs_exists "$src" || { echo "ERROR: panel not found: $src"; exit 1; }
  PANEL_LOCAL="${LOCAL}/panel.${c}.bcf"

  # header-only read (small) to test chr-naming; retry on transient reset (content-checked, since
  # `view -h` closes the pipe early and SIGPIPEs gsutil, so we can't trust the exit code).
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
    # not chr-named: download once and rename locally (no full-file stream)
    echo ">> [$c] contig not chr-named; downloading + renaming locally" >&2
    PANEL_CANON="$src"; ensure_panel_local
    local act cm="${LOCAL}/rename_${c}.txt" ren="${LOCAL}/panel.${c}.vcf.gz"
    act="$(bcftools view -H "${PANEL_LOCAL}" 2>/dev/null | head -1 | cut -f1)"
    printf '%s\t%s\n' "$act" "$c" > "$cm"
    bcftools annotate --rename-chrs "$cm" -Oz -o "$ren" "${PANEL_LOCAL}"
    rm -f "${PANEL_LOCAL}"; PANEL_LOCAL="$ren"          # subsequent reads use the renamed local file
  fi
}

prep_panel() {                             # $1=contig -> refsrc + unique_variants (uses PANEL_CANON)
  local c="$1"
  local REFSRC="${WORK}/ref/refsrc_noinfo.${c}.vcf.gz"
  local UNIQV="${REF_PREFIX}.${c}.unique_variants"
  if gs_exists "$REFSRC" && gs_exists "$UNIQV"; then echo ">> [$c] refsrc+uniqv exist; skipping panel prep"; return; fi
  echo ">> [$c] panel prep (local): strip INFO + norm -d exact -> refsrc + unique_variants"
  ensure_panel_local
  # Convert the BCF to bref3-ready VCF.gz. Strip INFO (bref3 ignores it; also drops the
  # malformed 'AC Number=A' header warning and trims per-record INFO). norm -d exact drops
  # exact-duplicate markers. (No symbolic-<SV> ALTs in this panel, so no ID-folding needed.)
  local RL="${LOCAL}/refsrc_noinfo.${c}.vcf.gz" UL="${LOCAL}/${c}.unique_variants"
  # Strip INFO; norm -d exact drops exact-duplicate markers. Also drop records where an ALT
  # equals REF (degenerate bubble traversals): Beagle carries them into the imputed VCF and
  # GATK/htsjdk (SelectVariants, GatherVcfs) reject the resulting duplicate allele. The panel
  # is bubble-split (biallelic), so dropping these non-variant records is lossless.
  bcftools annotate -x INFO -Ou "${PANEL_LOCAL}" \
    | bcftools norm -d exact --threads "${PANEL_THREADS}" -Ov - \
    | awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {n=split($5,a,","); bad=0; for(i=1;i<=n;i++) if(a[i]==$4){bad=1;break} if(bad){d++} else print} END{if(d) print "  [refsrc] dropped "d" REF==ALT record(s)" > "/dev/stderr"}' \
    | bcftools view -Oz -o "${RL}" --threads "${PANEL_THREADS}" -
  bcftools index -t --threads "${THREADS}" "${RL}"
  bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "${RL}" | LC_ALL=C sort -u --parallel="${THREADS}" -S 4G | sed '/^$/d' > "${UL}"
  gsutil cp "${RL}" "${REFSRC}"; gsutil cp "${RL}.tbi" "${REFSRC}.tbi"; gsutil cp "${UL}" "${UNIQV}"
}

build_bref3() {                            # $1=contig
  local c="$1" BREF3="${REF_PREFIX}.${c}.bref3" REFSRC="${WORK}/ref/refsrc_noinfo.${c}.vcf.gz"
  if gs_exists "$BREF3"; then echo ">> [$c] bref3 exists"; return; fi
  local JAR="${LOCAL}/bref3.jar"
  [ -s "$JAR" ] || { echo ">> downloading bref3 jar"; curl -fsSL -o "$JAR" "${BREF3_JAR_URL}"; }
  local BL="${LOCAL}/${c}.bref3"

  # Build bref3 from the prepped, CLEANED refsrc (INFO-stripped, deduped, REF==ALT dropped),
  # NOT the raw panel: the raw panel still carries the degenerate ALT==REF records that make
  # the imputed VCF fail GATK downstream. bref3 reads vcf.gz natively (no bcftools needed).
  local RL="${LOCAL}/refsrc_noinfo.${c}.vcf.gz"
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
  # guard: the staged map MUST be this contig (catches a wrong-file pick immediately)
  local MC; MC="$(head -1 "${LOCAL}/${c}.withchr.map" | cut -f1)"
  [ "$MC" = "$c" ] || { echo "ERROR: map for $c resolved to chrom '$MC' (wrong source file: $SRC)"; exit 1; }
  gsutil cp "${LOCAL}/${c}.withchr.map" "$GMAP"
}

corder() {                                 # chr name -> numeric order key
  local x="${1#chr}"
  case "$x" in
    X) echo 23;; Y) echo 24;; M|MT) echo 25;;
    ''|*[!0-9]*) echo 99;;
    *) echo "$x";;
  esac
}
_firstord() { corder "$(tabix -l "${SHARDS[$1]}" 2>/dev/null | head -1)"; }   # uses global SHARDS
_lb() {                                    # lower_bound: first shard with first-contig order >= $1
  local tgt="$1" lo=0 hi=${#SHARDS[@]} m fo
  while [ "$lo" -lt "$hi" ]; do
    m=$(( (lo + hi) / 2 )); fo="$(_firstord "$m")"
    if [ "$fo" -lt "$tgt" ]; then lo=$((m + 1)); else hi="$m"; fi
  done
  echo "$lo"
}

build_target_aou() {                       # $1=contig (host bcftools/tabix over the mount)
  local c="$1"
  local TGT="${TARGET_DIR}${TARGET_BASE}_${c}.snps.vcf.gz"
  if gs_exists "$TGT"; then echo ">> [$c] AoU target exists: $TGT"; return; fi
  local SAMP_GS="${SAMP_DIR}target_${N_TARGET}.txt" SAMP_LOC="${LOCAL}/target_${N_TARGET}.txt"
  local PANEL_SAMP="${LOCAL}/panel_samples.txt" AOU_SAMP="${LOCAL}/aou_samples.txt"

  SHARDS=()                                  # global (so _firstord/_lb see it)
  # readdir (ls -U: no stat, no shell glob -> no ARG_MAX) then filter by name; prefix + sort.
  mapfile -t SHARDS < <(ls -U -1 "${AOU_VCF_MOUNT}" 2>/dev/null | grep -E '\.vcf\.bgz$' | sed "s#^#${AOU_VCF_MOUNT}/#" | sort)
  if [ "${#SHARDS[@]}" -eq 0 ]; then
    local _N_ALL; _N_ALL="$(ls -U -1 "${AOU_VCF_MOUNT}" 2>/dev/null | wc -l)"
    echo "ERROR: no *.vcf.bgz under ${AOU_VCF_MOUNT} (readdir saw ${_N_ALL} entries)"
    if [ "${_N_ALL}" -gt 0 ]; then
      echo "  -> dir is non-empty but nothing matches '*.vcf.bgz'. First entries seen:"
      ls -U -1 "${AOU_VCF_MOUNT}" 2>/dev/null | head -5 | sed 's/^/       /'
      echo "     Check the real extension/subdir and set AOU_VCF_MOUNT accordingly."
    elif gsutil -u "${AOU_USER_PROJECT}" ls "${AOU_VCF_GS}/" 2>/dev/null | grep -qE '\.vcf\.bgz$'; then
      echo "  -> shards DO exist in ${AOU_VCF_GS} (requester-pays); the gcsfuse mount isn't exposing them."
      echo "     Remount with a billing project + implicit dirs, then set AOU_VCF_MOUNT:"
      echo "       gcsfuse --billing-project=${AOU_USER_PROJECT} --implicit-dirs -o ro vwb-aou-datasets-controlled \"\$HOME/aou_acaf_rp_mount\""
      echo "       AOU_VCF_MOUNT=\"\$HOME/aou_acaf_rp_mount/v8/wgs/short_read/snpindel/acaf_threshold/vcf\" ./prep_imputebeagle_aou.sh"
    else
      echo "  -> also not found via gsutil -u ${AOU_USER_PROJECT} ls ${AOU_VCF_GS}/ ; check AOU_VCF_GS / AOU_USER_PROJECT / permissions."
    fi
    exit 1
  fi
  echo ">> [$c] ${#SHARDS[@]} ACAF shards visible under the mount"

  # ---- candidate TSV: mirror local ./ -> bucket, and bucket -> local (future runs) ----
  local TSV_BN TSV_GS TSV_LOC=""
  TSV_BN="$(basename "${TARGET_TSV}")"; TSV_GS="${SAMP_DIR}${TSV_BN}"
  if [[ "${TARGET_TSV}" == gs://* ]]; then
    TSV_LOC="${LOCAL}/${TSV_BN}"; gsutil cp "${TARGET_TSV}" "${TSV_LOC}"
    gs_exists "${TSV_GS}" || gsutil cp "${TSV_LOC}" "${TSV_GS}"
  elif [ -s "${TARGET_TSV}" ]; then
    TSV_LOC="${TARGET_TSV}"
    gsutil cp "${TARGET_TSV}" "${TSV_GS}"; echo ">> [$c] TSV written to bucket: ${TSV_GS}"
  elif gs_exists "${TSV_GS}"; then
    TSV_LOC="${LOCAL}/${TSV_BN}"; gsutil cp "${TSV_GS}" "${TSV_LOC}"; echo ">> [$c] TSV pulled from bucket: ${TSV_GS}"
  fi

  # ---- fixed sample list: first N TSV ids available in ACAF (reuse from bucket if staged) ----
  # "available" = present in the ACAF callset; if EXCLUDE_PANEL, also not in the panel.
  gs_exists "$SAMP_GS" && gsutil cp "$SAMP_GS" "$SAMP_LOC" 2>/dev/null || true
  if [ ! -s "$SAMP_LOC" ]; then
    [ -s "${TSV_LOC}" ] || { echo "ERROR: candidate TSV not found (local ./${TSV_BN}, bucket ${TSV_GS}, or TARGET_TSV=${TARGET_TSV})"; exit 1; }
    bcftools query -l "${SHARDS[0]}" > "$AOU_SAMP"          # ACAF sample set (no header)
    local ALLOWED="${LOCAL}/allowed_samples.txt"
    if [ "${EXCLUDE_PANEL}" = "true" ]; then
      if [ ! -s "$PANEL_SAMP" ]; then
        ensure_panel_local
        bcftools view -h "${PANEL_LOCAL}" 2>/dev/null | sed -n '/^#CHROM/{p;q}' | cut -f10- | tr '\t' '\n' > "$PANEL_SAMP"
      fi
      awk 'NR==FNR{p[$0];next} !($0 in p)' "$PANEL_SAMP" "$AOU_SAMP" > "$ALLOWED"   # ACAF minus panel
    else
      cp "$AOU_SAMP" "$ALLOWED"
    fi
    # walk the TSV IN ORDER, keep ids that are allowed, take the first N (skip 'research_id' header)
    awk 'NR==FNR{a[$1]=1;next} (FNR==1 && $1=="research_id"){next} ($1 in a){print $1}' \
      "$ALLOWED" "$TSV_LOC" | head -n "${N_TARGET}" > "$SAMP_LOC"
    [ -s "$SAMP_LOC" ] || { echo "ERROR: 0 TSV ids available in ACAF (id namespace mismatch?)"; exit 1; }
    gsutil cp "$SAMP_LOC" "$SAMP_GS"
    echo ">> [$c] sample list written to bucket: ${SAMP_GS}"
  else
    echo ">> [$c] sample list reused from bucket: ${SAMP_GS}"
  fi
  local NTSV NACAF NSEL
  NTSV="$( [ -s "$TSV_LOC" ] && awk 'END{print NR-1}' "$TSV_LOC" || echo NA )"   # minus header
  NACAF="$(bcftools query -l "${SHARDS[0]}" | wc -l)"
  NSEL="$(wc -l < "$SAMP_LOC")"
  echo ">> [$c] candidates(TSV)=${NTSV} acaf=${NACAF} exclude_panel=${EXCLUDE_PANEL} -> target=${NSEL}"

  # ---- binary search shards for contig c (first-contig order via tabix -l) ----
  local want; want="$(corder "$c")"
  local L R; L="$(_lb "$want")"; R="$(_lb $((want + 1)))"
  local lo_i=$(( L - 1 < 0 ? 0 : L - 1 )) hi_i=$(( R < ${#SHARDS[@]} ? R : ${#SHARDS[@]} - 1 ))
  echo ">> [$c] shard bracket [${lo_i}..${hi_i}] of ${#SHARDS[@]}"

  # ---- membership-verified chr21 shard indices ----
  local sel=() i
  for i in $(seq "$lo_i" "$hi_i"); do
    tabix -l "${SHARDS[$i]}" 2>/dev/null | grep -qx "$c" && sel+=("$i")   # drop boundary/other-contig shards
  done
  local N="${#sel[@]}"
  [ "$N" -gt 0 ] || { echo "ERROR: no $c shards selected"; exit 1; }
  echo ">> [$c] ${N} ${c} shards; local subset PARALLEL=${PARALLEL} x ${SHARD_THREADS} threads each"

  # ---- per-shard subset over the mount, parallel with a concurrency cap; parts cached (resumable) ----
  # atomic via .tmp -> mv so an interrupted part is never mistaken for complete.
  local PARTS="${LOCAL}/target_parts_${c}"; mkdir -p "$PARTS"
  local running=0 part
  for i in "${sel[@]}"; do
    part="${PARTS}/part_$(printf '%010d' "$i").vcf.gz"
    [ -s "$part" ] && continue                                  # idempotent: keep finished parts
    (
      a=1
      until bcftools view -r "$c" -v snps -m2 -M2 -S "$SAMP_LOC" --force-samples \
              --threads "${SHARD_THREADS}" "${SHARDS[$i]}" -Oz -o "${part}.tmp"; do
        rm -f "${part}.tmp"
        [ "$a" -ge "${GS_RETRIES}" ] && { echo ">> [$c] shard idx $i FAILED after ${a} tries" >&2; exit 1; }
        echo ">> [$c] shard idx $i retry ${a}/${GS_RETRIES} (mount/GCS read); sleep $((a*10))s" >&2
        sleep "$((a*10))"; a=$((a+1))
      done
      mv "${part}.tmp" "$part"
    ) &
    running=$(( running + 1 ))
    if [ "$running" -ge "${PARALLEL}" ]; then wait -n 2>/dev/null || wait; running=$(( running - 1 )); fi
  done
  wait

  # ---- collect parts in genomic order (verify all present), merge + index, upload ----
  local PARTLIST=()
  for i in "${sel[@]}"; do
    part="${PARTS}/part_$(printf '%010d' "$i").vcf.gz"
    [ -s "$part" ] || { echo "ERROR: missing part for shard index $i ($part); rerun to resume"; exit 1; }
    PARTLIST+=("$part")
  done
  local COMB="${LOCAL}/${TARGET_BASE}_${c}.snps.vcf.gz"
  bcftools concat --threads "${MERGE_THREADS}" -Oz -o "$COMB" "${PARTLIST[@]}"
  bcftools index -t --threads "${MERGE_THREADS}" "$COMB"
  gsutil cp "$COMB" "$TGT"
  gsutil cp "${COMB}.tbi" "${TGT}.tbi"
  echo ">> [$c] AoU target staged: ${TGT}  ($(bcftools index -n "$COMB" 2>/dev/null || echo '?') records)"
}

# ============================ run ============================================
echo ">> project=${PROJECT}  bucket=${BUCKET}  mount=${AOU_VCF_MOUNT}  N_TARGET=${N_TARGET}"
for c in "${CONTIGS[@]}"; do
  echo "================= ${c} ================="
  resolve_panel "$c"
  echo ">> [$c] panel: ${PANEL_CANON}"
  prep_panel "$c"
  build_bref3 "$c"
  stage_map "$c"
  build_target_aou "$c"
  # free this contig's local intermediates (everything needed is in the bucket now)
  rm -rf "${PANEL_LOCAL}" \
         "${LOCAL}/refsrc_noinfo.${c}.vcf.gz" "${LOCAL}/refsrc_noinfo.${c}.vcf.gz.tbi" \
         "${LOCAL}/${c}.unique_variants" "${LOCAL}/${c}.bref3" \
         "${LOCAL}/target_parts_${c}" "${LOCAL}/${TARGET_BASE}_${c}.snps.vcf.gz" "${LOCAL}/${TARGET_BASE}_${c}.snps.vcf.gz.tbi"
done

# ============================ inputs: JSON + batch CSV + column map ===========
SCS="${SAMPLE_CHUNK_SIZE}"
CJ=""; for c in "${CONTIGS[@]}"; do CJ="${CJ}\"${c}\","; done; CJ="[${CJ%,}]"
mkdir -p "${LOCAL}/batch"
IN_JSON="${LOCAL}/batch/${TARGET_BASE}.inputs.json"
IN_CSV="${LOCAL}/batch/${TARGET_BASE}.inputs.csv"
IN_COLS="${LOCAL}/batch/${TARGET_BASE}.columns.json"
ECO_JSON=""
[ -n "${ERROR_COUNT_OVERRIDE}" ] && ECO_JSON=",
  \"ImputationBeagle.error_count_override\": ${ERROR_COUNT_OVERRIDE}"
cat > "${IN_JSON}" <<JSON
{
  "ImputationBeagle.multi_sample_vcf": "${TARGET_DIR}${TARGET_BASE}_${CONTIGS[0]}.snps.vcf.gz",
  "ImputationBeagle.ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
  "ImputationBeagle.contigs": ${CJ},
  "ImputationBeagle.reference_panel_path_prefix": "${REF_PREFIX}",
  "ImputationBeagle.genetic_maps_path": "${MAPS_DIR}",
  "ImputationBeagle.output_basename": "${OUT_BASE}",
  "ImputationBeagle.sample_chunk_size": ${SCS},
  "ImputationBeagle.beagle_phase_memory_in_gb": ${PHASE_MEM_GB},
  "ImputationBeagle.beagle_impute_memory_in_gb": ${IMPUTE_MEM_GB}${ECO_JSON}
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
BATCH_GCS="${WORK}/batch/"                         # gs:// dir (trailing slash)
BATCH_CSV_REL="${BATCH_GCS#${BUCKET}/}${TARGET_BASE}.inputs.csv"
COLS_URI="${BATCH_GCS}${TARGET_BASE}.columns.json"
gsutil cp "${IN_JSON}" "${IN_CSV}" "${IN_COLS}" "${BATCH_GCS}"

# ============================ WDL bundle =====================================
# Get the canonical WDLs into THIS workspace's bucket so we can register them.
# NOTE: cross-perimeter reads (e.g. from another workspace's bucket in a different
# VPC-SC perimeter) are blocked by org policy, so prefer uploading from a LOCAL
# dir on this VM (a within-perimeter write to ${BUCKET}, which is always allowed).
#   WDL_SRC_LOCAL=/path/to/imputationbeagle_wdl_flat   (dir holding the 4 .wdl)
# Auto-detect a local bundle on this VM so the within-perimeter upload is the DEFAULT; the
# cross-perimeter gs:// copy (WDL_SRC_GCS) is only used if no local bundle is found.
if [ -z "${WDL_SRC_LOCAL:-}" ]; then
  for d in "${HOME}/imputationbeagle_wdl_flat" "./imputationbeagle_wdl_flat" "${PWD}/imputationbeagle_wdl_flat"; do
    if [ -f "${d}/ImputationBeagle.wdl" ]; then WDL_SRC_LOCAL="$d"; echo ">> auto-detected local WDL bundle: ${d}"; break; fi
  done
fi
if [ -n "${WDL_SRC_LOCAL:-}" ]; then
  echo ">> uploading WDLs ${WDL_SRC_LOCAL%/}/*.wdl -> ${WDL_GCS}"
  gsutil cp "${WDL_SRC_LOCAL%/}/"*.wdl "${WDL_GCS}"
elif [ "${WDL_SRC_GCS%/}" != "${WDL_GCS%/}" ]; then
  echo ">> copying WDLs ${WDL_SRC_GCS}*.wdl -> ${WDL_GCS}"
  gsutil cp "${WDL_SRC_GCS}"*.wdl "${WDL_GCS}"
else
  echo ">> WDLs already in this workspace bucket (${WDL_GCS}); no copy needed"
fi

# ============================ register + run =================================
if [ "${RUN_WORKFLOW}" = "true" ]; then
  miss=0
  for w in ImputationBeagle.wdl ImputationBeagleStructs.wdl ImputationBeagleTasks.wdl ImputationTasks.wdl; do
    gs_exists "${WDL_GCS}${w}" || { echo "ERROR: missing ${WDL_GCS}${w}"; miss=1; }
  done
  [ "$miss" -eq 0 ] || { echo "       all 4 WDLs must be present in ${WDL_GCS} (set WDL_SRC_LOCAL=<dir with the 4 .wdl> to upload them from this VM, or WDL_SRC_GCS to a same-perimeter bucket)."; exit 1; }
  if [ "${FORCE_WF_REGISTER}" = "true" ] && wb workflow list 2>/dev/null | grep -qw "${WF_NAME}"; then
    echo ">> FORCE_WF_REGISTER: deleting stale ${WF_NAME} so its WDL path is rebuilt"
    wb workflow delete --workflow="${WF_NAME}" --quiet 2>/dev/null || true
  fi
  if wb workflow list 2>/dev/null | grep -qw "${WF_NAME}"; then
    echo ">> workflow ${WF_NAME} already registered (rerun with FORCE_WF_REGISTER=true if its WDL path is stale)"
  else
    echo ">> registering workflow ${WF_NAME} -> ${WDL_GCS}ImputationBeagle.wdl"
    wb workflow create --bucket-id="${BUCKET_ID}" \
      --path="${WDL_REL}/ImputationBeagle.wdl" --workflow="${WF_NAME}"
  fi
  echo ">> launching ${WF_NAME}: output-path=${OUTPUT_PATH}, sample_chunk_size=${SCS}"
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
    || { echo ">> batch-CSV run failed (inputs at ${BATCH_GCS}). Try USE_BATCH_CSV=false ./prep_imputebeagle_aou.sh"; exit 1; }
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
  # the THIS-job output dir (timestamped, single contig) — scope the imputed-VCF search to it
  # so a stale *.imputed.vcf.gz from an earlier run in the same namespace can't be picked up.
  JOB_OUT="$(grep -m1 'Output bucket path:' "${RUN_LOG}" 2>/dev/null | sed 's/.*Output bucket path:[[:space:]]*//' | tr -d '[:space:]' || true)"
  [ -n "${JOB_OUT}" ] && echo ">> job output dir: ${JOB_OUT}"
fi

# ============================ wait for output + sanity =======================
if [ "${RUN_WORKFLOW}" = "true" ] && [ "${RUN_EVAL}" = "true" ]; then
  echo ">> waiting for ${OUT_BASE}.imputed.vcf.gz under ${JOB_OUT:-${BUCKET}/${OUTPUT_PATH}}/ (poll ${EVAL_POLL_SECS}s x ${EVAL_MAX_POLLS})"

  # locate the final merged imputed VCF for THIS job only. Scope to JOB_OUT (the
  # timestamped, single-contig run dir) so a stale *.imputed.vcf.gz from an earlier
  # run in the same namespace (e.g. a previous chr21) can never be matched.
  find_imp() {
    local root="${JOB_OUT:-${BUCKET}/${OUTPUT_PATH}}" hit
    hit="$(gsutil ls -r "${root%/}/" 2>/dev/null | grep -E "/${OUT_BASE}\.imputed\.vcf\.gz$" | head -1)"
    [ -n "$hit" ] || hit="$(gsutil ls -r "${root%/}/" 2>/dev/null | grep -E '\.imputed\.vcf\.gz$' | grep -vE 'no_overlaps|hom_ref' | head -1)"
    printf '%s' "$hit"
  }

  # aggregate job state via `wb workflow job describe` (works for single + batch job ids).
  # defaults to RUNNING on any parse uncertainty so we never abort a healthy run.
  job_state() {
    local id="${1:-}" out f o s p
    [ -n "$id" ] || { echo UNKNOWN; return; }
    out="$(wb workflow job describe --job-id="$id" 2>/dev/null)" || { echo UNKNOWN; return; }
    if grep -q 'Batch job summary' <<<"$out"; then               # batch: decide from run counts
      f=$(awk -F':[[:space:]]*' '/Failed runs:/{print $2;exit}'   <<<"$out"); f="${f//[!0-9]/}"; : "${f:=0}"
      o=$(awk -F':[[:space:]]*' '/Ongoing runs:/{print $2;exit}'  <<<"$out"); o="${o//[!0-9]/}"; : "${o:=0}"
      s=$(awk -F':[[:space:]]*' '/Starting runs:/{print $2;exit}' <<<"$out"); s="${s//[!0-9]/}"; : "${s:=0}"
      p=$(awk -F':[[:space:]]*' '/Pending runs:/{print $2;exit}'  <<<"$out"); p="${p//[!0-9]/}"; : "${p:=0}"
      if [ "$((o+s+p))" -eq 0 ]; then { [ "$f" -gt 0 ] && echo FAILED || echo DONE; }; else echo RUNNING; fi
      return
    fi
    case "$(awk -F':[[:space:]]*' '/^Status:/{print $2;exit}' <<<"$out" | tr '[:lower:]' '[:upper:]')" in
      *FAIL*|*CANCEL*|*ABORT*|*ERROR*) echo FAILED ;;
      *SUCC*|*COMPLETE*|*DONE*)        echo DONE ;;
      *)                               echo RUNNING ;;
    esac
  }

  IMP="" ST="RUNNING"
  for i in $(seq 1 "${EVAL_MAX_POLLS}"); do
    IMP="$(find_imp)"
    [ -n "$IMP" ] && { ST=DONE; break; }
    ST="$(job_state "${JOB_ID:-}")"
    if [ "$ST" = FAILED ]; then
      echo ">> job FAILED — stopping poll. Status detail:"
      wb workflow job describe --job-id="${JOB_ID:-}" 2>/dev/null \
        | grep -E 'Status|Failed runs|Output bucket' | sed 's/^/     /'
      [ "${USE_BATCH_CSV}" = "true" ] \
        && echo "     per-run: wb workflow job list --batch-job-id=${JOB_ID:-} --status=FAILED  (then 'job describe' a failed id)"
      break
    fi
    if [ "$ST" = DONE ]; then            # job finished -> stop polling and locate the output
      echo ">> job state DONE; locating imputed VCF under ${BUCKET}/${OUTPUT_PATH}/"
      IMP="$(find_imp)"
      break
    fi
    echo ">> [poll ${i}/${EVAL_MAX_POLLS}] not ready (state=${ST}); sleeping ${EVAL_POLL_SECS}s"
    sleep "${EVAL_POLL_SECS}"
  done

  if [ -n "$IMP" ]; then
    echo ">> found imputed VCF: ${IMP}"
    IMPUTED_VCF="${IMP}" PANEL_VCF="${WORK}/ref/refsrc_noinfo.${CONTIGS[0]}.vcf.gz" \
      CONTIG="${CONTIGS[0]}" OUT_BASE="${OUT_BASE}" RESDIR="${LOCAL}/sanity" \
      python3 "${SANITY_PY}"
    echo ">> sanity outputs in ${LOCAL}/sanity/"
  elif [ "$ST" = FAILED ]; then
    echo ">> workflow failed before producing ${OUT_BASE}.imputed.vcf.gz; skipping sanity (debug with the commands above)."
  elif [ "$ST" = DONE ]; then
    echo ">> job DONE but no *.imputed.vcf.gz located under ${BUCKET}/${OUTPUT_PATH}/ ; list it with:"
    echo "   gsutil ls -r ${BUCKET}/${OUTPUT_PATH}/ | grep -E '\\.imputed\\.vcf\\.gz\$'"
  else
    echo ">> imputed VCF not found within timeout (last state=${ST}). Check: wb workflow job describe --job-id=${JOB_ID:-<id>}"
    echo "   then run sanity later with:"
    echo "   IMPUTED_VCF=<gs://...> PANEL_VCF=${WORK}/ref/refsrc_noinfo.${CONTIGS[0]}.vcf.gz CONTIG=${CONTIGS[0]} OUT_BASE=${OUT_BASE} RESDIR=${LOCAL}/sanity python3 ${SANITY_PY}"
  fi
fi

echo
echo "================= DONE (sample_chunk_size=${SCS}) ================="
echo "target:        ${TARGET_DIR}${TARGET_BASE}_${CONTIGS[0]}.snps.vcf.gz  (N=${N_TARGET} AoU, out of panel)"
echo "candidate TSV: ${SAMP_DIR}$(basename "${TARGET_TSV}")"
echo "sample list:   ${SAMP_DIR}target_${N_TARGET}.txt"
echo "ref bref3:     ${REF_PREFIX}.${CONTIGS[0]}.bref3"
echo "inputs:        CSV ${BATCH_GCS}${TARGET_BASE}.inputs.csv | cols ${COLS_URI}"
echo "WDL bundle:    ${WDL_GCS}"
if [ "${RUN_WORKFLOW}" = "true" ]; then
  echo "workflow:      ${WF_NAME}  ->  output ${BUCKET}/${OUTPUT_PATH}/"
  echo "monitor:       wb workflow job list ; wb workflow job describe --job-id=<id>"
  echo "sanity:        ${LOCAL}/sanity/${OUT_BASE}.sanity.png  +  .sv_sanity.tsv"
fi

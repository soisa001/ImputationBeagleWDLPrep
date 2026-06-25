# =============================================================================
# Representation matching test: ACAF SNV scaffold  vs  panel bubble alleles.
#
# Quantifies how many panel bubble-split alleles the short-read (ACAF) scaffold
# can reach under (a) CURRENT exact CHROM:POS:REF:ALT matching, the way Beagle
# joins target<->reference markers, vs (b) MINIMAL-REPRESENTATION matching, the
# normalization extract-bubble-PLs applies (trim shared suffix, then prefix).
#
# The delta = bubble alleles the panel encodes with flanking context (e.g.
# GAT/GGT) that a minimal short-read SNV (A/G at +1) only matches after
# normalization. True multi-change MNV paths (GA/CG) are NOT recovered by
# min-rep alone (they need atomization) and are reported separately so the
# ceiling of this approach is explicit.
# =============================================================================
import subprocess, collections
import matplotlib.pyplot as plt

# ----------------------------- config ----------------------------------------
REGION = "chr1:1000000-3000000"            # test window (indexed reads; keep modest)
# PANEL = the bubble.split reference Beagle imputes against (biallelic bubble alleles).
PANEL  = f"{__import__('os').path.expanduser('~')}/imp_holdout198/leaveout.chr1.bcf"
# QUERY = the ACAF / short-read calls that form the scaffold (one shard, or a pre-cut region VCF).
QUERY  = "/path/to/acaf_or_target.chr1.vcf.gz"   # <-- set to your ACAF shard / target VCF
# bcftools must be able to read + index-seek these (local files, or gs:// if your build supports it).

# --------------------------- minimal representation --------------------------
def minrep(pos, ref, alt):
    """Parsimonious left-aligned key: trim shared suffix, then shared prefix (advancing pos).
    Mirrors get_minimal_representation in extract-bubble-PLs. Symbolic <...> alleles pass through."""
    pos = int(pos)
    if alt[:1] == "<" or ref in (".", "") or alt == ".":
        return (pos, ref, alt)
    r, a = ref, alt
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    while r and a and r[0] == a[0]:
        r, a = r[1:], a[1:]; pos += 1
    return (pos, r, a)

def kind(ref, alt):
    """Classify the (already min-rep) allele."""
    if alt[:1] == "<":            return "symbolic"
    if len(ref) == 1 and len(alt) == 1:   return "SNV"
    if len(ref) <= 1 or len(alt) <= 1:    return "indel"       # one side empty/1bp after trim
    return "MNV/complex"          # both sides >1bp after trim -> needs atomization, not just min-rep

# --------------------------- red-green self-test ------------------------------
_cases = [  # (panel pos,ref,alt), (acaf pos,ref,alt), expect_naive, expect_minrep
    ((1000,"A","G"),    (1000,"A","G"),   True,  True),    # exact SNV
    ((1000,"GAT","GGT"),(1001,"A","G"),   False, True),    # padded SNV  -> recovered
    ((1000,"GAT","GT"), (1000,"GA","G"),  False, True),    # padded DEL  -> recovered
    ((1000,"GA","GTA"), (1000,"G","GT"),  False, True),    # padded INS  -> recovered
    ((1000,"GA","CG"),  (1000,"G","C"),   False, False),   # true MNV    -> NOT recovered
]
for (pp,pr,pa),(qp,qr,qa),en,em in _cases:
    assert ((pp,pr,pa)==(qp,qr,qa)) == en
    assert (minrep(pp,pr,pa)==minrep(qp,qr,qa)) == em
print("self-test OK (padded SNV/indel recovered; true MNV not)\n")

# ----------------------------- load records ----------------------------------
def fetch(path, region, snv_only=False):
    """Return list of (chrom,pos,ref,alt) for each ALT allele in the region."""
    cmd = ["bcftools", "view", "-r", region]
    if snv_only:
        cmd += ["-v", "snps", "-m2", "-M2"]
    cmd += [path]
    q = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n"],
                         stdin=q.stdout, capture_output=True, text=True)
    q.stdout.close(); q.wait()
    recs = []
    for line in out.stdout.splitlines():
        c, p, r, a = line.split("\t")
        for alt in a.split(","):                 # bubble.split is biallelic; split defensively
            recs.append((c, int(p), r, alt))
    return recs

print(f"reading panel : {PANEL}")
panel = fetch(PANEL, REGION, snv_only=False)         # all bubble alleles (SNV/indel/MNV)
print(f"reading scaffold (ACAF SNVs): {QUERY}")
acaf  = fetch(QUERY, REGION, snv_only=True)          # the short-read SNV scaffold

# scaffold key sets (the ACAF side)
acaf_naive  = {(c, p, r, a) for (c, p, r, a) in acaf}
acaf_minrep = {(c,) + minrep(p, r, a) for (c, p, r, a) in acaf}

# --------------------------- compare matching --------------------------------
# panel allele is "reachable" if the scaffold carries a matching marker.
matched_naive, matched_minrep, recovered = [], [], []
for (c, p, r, a) in panel:
    n_hit = (c, p, r, a) in acaf_naive
    m_hit = ((c,) + minrep(p, r, a)) in acaf_minrep
    if n_hit: matched_naive.append((c, p, r, a))
    if m_hit: matched_minrep.append((c, p, r, a))
    if m_hit and not n_hit: recovered.append((c, p, r, a))

n_panel = len(panel)
n_naive = len(matched_naive)
n_mr    = len(matched_minrep)
n_rec   = len(recovered)

# residual ceiling: panel alleles that stay unmatched AND are MNV/complex after min-rep
unmatched = [(c,p,r,a) for (c,p,r,a) in panel if ((c,)+minrep(p,r,a)) not in acaf_minrep]
mnv_residual = [t for t in unmatched if kind(*minrep(t[1],t[2],t[3])[1:]) == "MNV/complex"]

print(f"\nregion {REGION}")
print(f"panel bubble alleles            : {n_panel}")
print(f"reachable, exact (current)      : {n_naive}  ({100*n_naive/max(n_panel,1):.2f}%)")
print(f"reachable, minimal-rep          : {n_mr}  ({100*n_mr/max(n_panel,1):.2f}%)")
print(f"RECOVERED by min-rep            : {n_rec}  (+{100*n_rec/max(n_naive,1):.2f}% over current)")
print(f"still-unmatched MNV/complex     : {len(mnv_residual)}  (would need atomization, not min-rep)")

# breakdown of what got recovered, by allele type (on the min-rep allele, so padded SNVs read as SNV)
by_kind = collections.Counter(kind(*minrep(p, r, a)[1:]) for (_, p, r, a) in recovered)
print("\nrecovered, by panel allele type :", dict(by_kind))

print("\nexamples (panel raw -> min-rep, matched ACAF SNV):")
seen = 0
for (c, p, r, a) in recovered:
    mp = minrep(p, r, a)
    print(f"  {c}:{p} {r}>{a:<12} -> {c}:{mp[0]} {mp[1] or '-'}>{mp[2] or '-'}")
    seen += 1
    if seen >= 8: break

# ------------------------------- plot ----------------------------------------
plt.rcParams.update({"font.size": 26})
fig, ax = plt.subplots(figsize=(10, 8))
bars = ax.bar(["exact\n(current)", "minimal-rep\n(extract-bubble-PLs)"],
              [n_naive, n_mr], color=["#9aa0a6", "#1a73e8"], width=0.6)
ax.bar_label(bars, fmt="%d", fontsize=26, padding=6)
ax.set_ylabel("panel bubble alleles reached", fontsize=26)
ax.set_title(f"Scaffold reach over {REGION}\n+{n_rec} recovered ({100*n_rec/max(n_naive,1):.1f}%)",
             fontsize=26)
ax.tick_params(labelsize=26)
plt.tight_layout()
plt.show()

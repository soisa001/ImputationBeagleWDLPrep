#!/usr/bin/env python3
# =============================================================================
# Project a target VCF onto the PANEL's allele representation via minimal-
# representation matching, so panel bubble alleles that encode an isolated SNV
# with flanking context (e.g. GAT/GGT) are reachable by a minimal short-read
# call (A/G at +1). GT-scaffold analog of what extract-bubble-PLs does for GLIMPSE2.
#
# MULTIALLELIC-AWARE: a multiallelic input record is decomposed into one
# biallelic record per ALT, with GTs recoded exactly like `bcftools norm -m -any`
# (focal ALT -> 1, REF -> 0, any other ALT -> 0, missing -> .). So feeding raw
# multiallelic ACAF directly is fine; splitting upstream with `bcftools norm
# -m -any` first is also fine (idempotent). Each per-ALT allele is then reduced
# to minimal rep and looked up in the panel's minimal-rep -> raw map: if found,
# emit at the PANEL's raw (POS,REF,ALT) carrying the recoded GTs (GT-only,
# UNPHASED); else drop (a marker absent from the reference is unusable by Beagle,
# and non-SNV alts like deletions fall away since the map holds SNV-minrep keys).
# Exact matches are preferred so already-aligned SNVs are emitted unchanged.
#
# Usage (either works; do NOT use -m2 -M2/-v snps or multiallelics are lost):
#   bcftools view <acaf.vcf[.gz]> | project_to_panel_rep.py <panel_sites.tsv[.gz]> > projected.vcf
#   bcftools view <acaf.vcf[.gz]> | bcftools norm -m -any | project_to_panel_rep.py <panel_sites.tsv[.gz]> > projected.vcf
#   (panel_sites = CHROM<TAB>POS<TAB>REF<TAB>ALT, e.g. bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n')
# Output is NOT position-sorted (padded panel alleles sit at lower POS) -> pipe through `bcftools sort`.
# =============================================================================
import sys, gzip

def minrep(pos, ref, alt):
    pos = int(pos)
    if alt[:1] == "<":
        return (pos, ref, alt)
    r, a = ref, alt
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    while r and a and r[0] == a[0]:
        r, a = r[1:], a[1:]; pos += 1
    return (pos, r, a)

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def recode_gt(g, k):
    """Decompose a (possibly multiallelic) genotype onto the focal ALT index k:
    0->0, k->1, any other ALT->0, missing->. ; output UNPHASED. Mirrors the GT
    recoding of `bcftools norm -m -any` (a het-of-two-alts 1/2 -> 1/0 for k=1,
    0/1 for k=2). For a biallelic record (k=1) this is a pass-through+unphase."""
    out = []
    for a in g.replace("|", "/").split("/"):
        if a in (".", ""):
            out.append(".")
        else:
            try:
                out.append("1" if int(a) == k else "0")
            except ValueError:
                out.append(".")
    return "/".join(out)

def main():
    if len(sys.argv) < 2:
        sys.stderr.write("usage: ... | project_to_panel_rep.py <panel_sites.tsv[.gz]>\n"); sys.exit(1)

    # --- build panel minimal-rep -> raw map (SNV-minrep keys only; a SNV target can
    #     only match a panel allele whose minimal rep is that same SNV). Prefer the
    #     exact representation when several panel alleles share a minrep. ---
    pmap = {}
    with opener(sys.argv[1]) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            c, p, r, a = parts[0], parts[1], parts[2], parts[3]
            for alt in a.split(","):
                mp, mr, ma = minrep(p, r, alt)
                if len(mr) != 1 or len(ma) != 1:        # keep only SNV-minrep keys
                    continue
                key = (c, mp, mr, ma)
                raw = (c, int(p), r, alt)
                cur = pmap.get(key)
                if cur is None or (raw[1] == mp and raw[2] == mr and raw[3] == ma):
                    pmap[key] = raw                     # prefer exact (raw == minrep)
    sys.stderr.write(f"[project] panel SNV-minrep keys: {len(pmap)}\n")

    n_in = n_exact = n_recov = n_drop = 0
    emitted = set()
    w = sys.stdout.write
    for line in sys.stdin:
        if line[:1] == "#":
            w(line); continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        c, p, ref, alt_field = f[0], f[1], f[3], f[4]
        fmt = f[8].split(":")
        gi = fmt.index("GT") if "GT" in fmt else 0
        raw_gt = [s.split(":")[gi] if s not in (".", "") else "./." for s in f[9:]]
        alts = alt_field.split(",")                       # MULTIALLELIC-aware: one biallelic record per ALT
        for k, alt in enumerate(alts, start=1):
            n_in += 1
            praw = pmap.get((c,) + minrep(p, ref, alt))
            if praw is None:
                n_drop += 1; continue
            key = (praw[0], praw[1], praw[2], praw[3])
            if key in emitted:                            # dedupe collisions onto the same panel allele
                continue
            emitted.add(key)
            # recode GTs for this focal ALT k (0->0, k->1, any other ALT->0, missing->.) UNPHASED.
            # This matches `bcftools norm -m -any`; for a biallelic record (k=1) it is a pass-through.
            gts = [recode_gt(g, k) for g in raw_gt]
            if praw[1] == int(p) and praw[2] == ref and praw[3] == alt:
                n_exact += 1
            else:
                n_recov += 1
            w("\t".join([praw[0], str(praw[1]), ".", praw[2], praw[3], ".", ".", ".", "GT", *gts]) + "\n")

    sys.stderr.write(f"[project] in(alleles)={n_in} exact={n_exact} recovered={n_recov} dropped(no panel match)={n_drop}\n")

if __name__ == "__main__":
    main()

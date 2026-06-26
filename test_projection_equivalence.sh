#!/usr/bin/env bash
# =============================================================================
# Differential equivalence test: the Rust projector (project_to_panel_rep.rs)
# vs the reference Python (project_to_panel_rep.py). Builds the Rust binary,
# generates an adversarial + randomized (derived-from-panel) corpus, runs both,
# and asserts byte-identical stdout AND the [project] stderr summary.
#
#   bash test_projection_equivalence.sh
#
# Exit 0 = identical on every case; non-zero = a mismatch (diff printed).
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${HERE}/project_to_panel_rep.py"
RS="${HERE}/project_to_panel_rep.rs"
[ -s "$PY" ] && [ -s "$RS" ] || { echo "ERROR: need project_to_panel_rep.py and .rs next to this script"; exit 1; }
command -v rustc   >/dev/null 2>&1 || { echo "ERROR: rustc not found"; exit 2; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found"; exit 2; }

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
echo ">> building rust projector"
rustc -O "$RS" -o "${TMP}/project_bin"

gen_and_check() {   # $1=label  (reads PSITES from $TMP/p.tsv, VCF from $TMP/in.vcf)
  python3 "$PY" "${TMP}/p.tsv" < "${TMP}/in.vcf" > "${TMP}/py.out" 2> "${TMP}/py.err"
  "${TMP}/project_bin" "${TMP}/p.tsv" < "${TMP}/in.vcf" > "${TMP}/rs.out" 2> "${TMP}/rs.err"
  local emit; emit="$(grep -vc '^#' "${TMP}/py.out" || true)"
  if diff -q "${TMP}/py.out" "${TMP}/rs.out" >/dev/null \
     && diff <(grep '^\[project\]' "${TMP}/py.err") <(grep '^\[project\]' "${TMP}/rs.err") >/dev/null; then
    echo "   [$1] IDENTICAL (${emit} emitted records; $(grep 'in(alleles)' "${TMP}/py.err"))"
  else
    echo "   [$1] MISMATCH:"; diff "${TMP}/py.out" "${TMP}/rs.out" | head -40; exit 1
  fi
}

# --- Case 1: adversarial (exact, padded-recovery, dedup, indel-drop, multiallelic, phased, missing) ---
printf 'chr1\t1000\tA\tG\nchr1\t1000\tGAT\tGGT\nchr1\t1001\tA\tG\nchr1\t2000\tC\tCAG\nchr1\t3000\tA\tT,G\nchr1\t3500\tA\tT,G\nchr1\t4000\tACGT\tACGT\nchr1\t5000\tN\tN\n' > "${TMP}/p.tsv"
cat > "${TMP}/in.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	1000	.	GAT	GGT	.	.	.	GT	0/1	0/0	.
chr1	1001	.	A	G	.	.	.	GT	0/1	1/1	0/0
chr1	2000	.	C	CAG	.	.	.	GT	0/1	0/0	1/1
chr1	3000	.	A	T	.	.	.	GT	1/1	0/1	./.
chr1	3000	.	A	G	.	.	.	GT	0/1	1|0	0/0
chr1	3500	.	A	T,G	.	.	.	GT	0/1	1/2	2/2
chr1	9999	.	A	C	.	.	.	GT	0/1	0/0	1/1
EOF
echo ">> running cases"
gen_and_check "adversarial"

# --- Cases 2-4: randomized, input derived from the panel min-reps so matching is exercised ---
for SEED in 1 17 101; do
python3 - "$SEED" "${TMP}/p.tsv" "${TMP}/in.vcf" <<'PY'
import random, sys
seed=int(sys.argv[1]); random.seed(seed); B="ACGT"
def minrep(pos,r,a):
    r=list(r); a=list(a)
    while len(r)>1 and len(a)>1 and r[-1]==a[-1]: r.pop(); a.pop()
    while len(r)>1 and len(a)>1 and r[0]==a[0]: r.pop(0); a.pop(0); pos+=1
    return pos,"".join(r),"".join(a)
panel=[]
for _ in range(4000):
    pos=random.randint(100,99000); r0=random.choice(B); a0=random.choice([b for b in B if b!=r0]); k=random.random()
    if k<0.5: panel.append((pos,r0,a0))
    elif k<0.8:
        pre="".join(random.choice(B) for _ in range(random.randint(0,2))); suf="".join(random.choice(B) for _ in range(random.randint(0,2)))
        panel.append((pos-len(pre),pre+r0+suf,pre+a0+suf))
    else: panel.append((pos,r0,r0+"".join(random.choice(B) for _ in range(random.randint(1,3)))))
with open(sys.argv[2],"w") as f:
    i=0
    while i<len(panel):
        pos,r,a=panel[i]
        if i+1<len(panel) and panel[i+1][0]==pos and panel[i+1][1]==r and random.random()<0.3: a=a+","+panel[i+1][2]; i+=1
        f.write(f"chr1\t{pos}\t{r}\t{a}\n"); i+=1
def gt():
    if random.random()<0.05: return "."
    al=lambda: random.choice(["0","1","2","."]); return al()+random.choice(["/","|"])+al()
rows=[]
for pos,r,a in panel:
    a1=a.split(",")[0]; mp,mr,ma=minrep(pos,r,a1); c=random.random()
    if c<0.5: rows.append((mp,mr,ma))
    elif c<0.8: rows.append((pos,r,a1))
    else: rows.append((random.randint(100,99000),random.choice(B),random.choice(B)))
random.shuffle(rows)
NS=4
with open(sys.argv[3],"w") as f:
    f.write('##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(f"S{i}" for i in range(NS))+"\n")
    for pos,r,a in rows: f.write(f"chr1\t{pos}\t.\t{r}\t{a}\t.\t.\t.\tGT\t"+"\t".join(gt() for _ in range(NS))+"\n")
PY
gen_and_check "random.seed=${SEED}"
done

echo "RESULT: rust projector is byte-identical to the python on every case."

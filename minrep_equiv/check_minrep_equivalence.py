#!/usr/bin/env python3
# =============================================================================
# Differential equivalence check for the minimal-representation matching shared
# by project_to_panel_rep.py (Python) and extract-bubble-PLs (Rust).
#
# It does NOT claim the two PROGRAMS are equivalent (they aren't: PL vs GT,
# windowed buffer, gVCF fallback). It verifies the ONE thing that must agree for
# matching to be correct: get_minimal_representation, i.e. whether a given
# (POS,REF,ALT) reduces to the same key on both sides. If the minreps agree on
# every input, then "alleles match iff minreps equal" makes identical decisions
# in both tools.
#
# Usage:
#   rustc -O minrep_cli.rs -o minrep_cli           # build the Rust side (VM has rustc)
#   python3 check_minrep_equivalence.py ./minrep_cli           # real differential test
#   python3 check_minrep_equivalence.py --self-test            # harness sanity (no Rust)
#
# Exit 0 = all cases agree; non-zero = mismatches (first 20 printed).
# =============================================================================
import sys, os, random, subprocess, shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
try:
    from project_to_panel_rep import minrep            # the REAL function under test
except Exception:
    # fallback: same source, inlined (keep identical to project_to_panel_rep.minrep)
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

# byte-faithful re-implementation of the Rust fn, used ONLY for --self-test
# (so the harness can be exercised without a Rust toolchain). It is a separate
# transcription on purpose: agreement here checks the harness, not the tools.
def _rust_ref(pos, ref, alt):
    pos = int(pos)
    r = ref.encode(); a = alt.encode()
    if len(a) > 0 and a[0:1] == b"<":
        return (pos, ref, alt)
    while len(r) and len(a) and r[-1:] == a[-1:]:
        r = r[:-1]; a = a[:-1]
    while len(r) and len(a) and r[0:1] == a[0:1]:
        r = r[1:]; a = a[1:]; pos += 1
    return (pos, r.decode(), a.decode())

def gen_cases(n, seed=0):
    """Adversarial + random corpus of (pos, ref, alt) with non-empty ref/alt."""
    rng = random.Random(seed)
    bases = "ACGT"
    # fixed adversarial cases first (the ones that trip naive implementations)
    fixed = [
        (1000, "A", "G"),            # plain SNV
        (1000, "GAT", "GGT"),        # SNV padded both? (shared suffix T, shared prefix G) -> 1001 A/G
        (1000, "GAT", "GT"),         # deletion of A -> 1001 A/-
        (1000, "GA", "GTA"),         # insertion -> 1001 -/T
        (1000, "GA", "CG"),          # true MNV (no shared flanks)
        (1000, "ACGT", "ACGT"),      # ref==alt -> empties both
        (1000, "AAAA", "AA"),        # repeat-collapse deletion
        (1000, "TA", "TATA"),        # tandem insertion
        (1000, "N", "N"),            # N handling
        (1000, "A", "<DEL>"),        # symbolic -> passthrough
        (1000, "ACGT", "<INS:ME:L1>"),
        (1000, "CAG", "C"),          # left-anchored deletion -> 1001 AG/-
        (1000, "C", "CAG"),          # left-anchored insertion -> 1001 -/AG
        (1000, "GAAG", "GAG"),       # internal-looking del, shared prefix+suffix
        (1000, "ATAT", "AT"),        # suffix-collapse
        (1, "A", "T"), (999999999, "G", "C"),   # pos bounds
    ]
    cases = list(fixed)
    for _ in range(n):
        pos = rng.randint(1, 250_000_000)
        lr = rng.randint(1, 8); la = rng.randint(1, 8)
        # bias toward shared flanks to stress the trimming
        pre = "".join(rng.choice(bases) for _ in range(rng.randint(0, 3)))
        suf = "".join(rng.choice(bases) for _ in range(rng.randint(0, 3)))
        core_r = "".join(rng.choice(bases) for _ in range(lr))
        core_a = "".join(rng.choice(bases) for _ in range(la))
        ref = pre + core_r + suf
        alt = pre + core_a + suf
        if rng.random() < 0.05:                 # occasional symbolic alt
            alt = "<DEL>"
        cases.append((pos, ref, alt))
    return cases

def py_minrep_str(pos, ref, alt):
    p, r, a = minrep(pos, ref, alt)
    return f"{p}\t{r or '-'}\t{a or '-'}"

def main():
    args = sys.argv[1:]
    self_test = "--self-test" in args
    rust_bin = next((a for a in args if not a.startswith("-")), None)
    n = 200_000
    for a in args:
        if a.startswith("--n="):
            n = int(a.split("=", 1)[1])

    cases = gen_cases(n)
    inp = "".join(f"{p}\t{r}\t{a}\n" for (p, r, a) in cases)

    # reference side
    if self_test:
        ref_lines = [f"{p}\t{r or '-'}\t{a or '-'}" for (p, r, a) in
                     (_rust_ref(p_, r_, a_) for (p_, r_, a_) in cases)]
        ref_name = "byte-reference (Python; harness self-test only)"
    else:
        # auto-build ./minrep_cli from minrep_cli.rs if no binary was given and rustc is present
        if not rust_bin:
            here = os.path.dirname(os.path.abspath(__file__))
            cand = os.path.join(here, "minrep_cli")
            src = os.path.join(here, "minrep_cli.rs")
            if not os.path.exists(cand) and os.path.exists(src) and shutil.which("rustc"):
                print(f"building {cand} from minrep_cli.rs ...")
                b = subprocess.run(["rustc", "-O", src, "-o", cand], capture_output=True, text=True)
                if b.returncode != 0:
                    print("rustc build failed:\n", b.stderr); sys.exit(2)
            if os.path.exists(cand):
                rust_bin = cand
        if not rust_bin or not os.path.exists(rust_bin):
            print("ERROR: no Rust binary. Either install rustc (the checker will build minrep_cli.rs),\n"
                  "  pass a path:   python3 check_minrep_equivalence.py ./minrep_cli\n"
                  "  or self-test:  python3 check_minrep_equivalence.py --self-test")
            sys.exit(2)
        proc = subprocess.run([rust_bin], input=inp, capture_output=True, text=True)
        if proc.returncode != 0:
            print("Rust binary failed:\n", proc.stderr); sys.exit(2)
        ref_lines = proc.stdout.splitlines()
        ref_name = f"Rust extract-bubble-PLs minrep ({rust_bin})"

    py_lines = [py_minrep_str(p, r, a) for (p, r, a) in cases]

    assert len(ref_lines) == len(py_lines) == len(cases), \
        f"line count mismatch: rust={len(ref_lines)} py={len(py_lines)} cases={len(cases)}"

    mism = []
    for (p, r, a), rl, pl in zip(cases, ref_lines, py_lines):
        if rl != pl:
            mism.append((f"{p}\t{r}\t{a}", rl, pl))

    print(f"compared {len(cases)} cases: Python minrep  vs  {ref_name}")
    if not mism:
        print("RESULT: IDENTICAL on every case  ->  the matching decision is equivalent.")
        sys.exit(0)
    print(f"RESULT: {len(mism)} MISMATCH(es). First {min(20, len(mism))}:")
    for inp_s, rl, pl in mism[:20]:
        print(f"  in[{inp_s}]  rust[{rl}]  py[{pl}]")
    sys.exit(1)

if __name__ == "__main__":
    main()

// minrep_cli.rs — differential-test harness.
// Wraps the EXACT get_minimal_representation function copied verbatim from the
// extract-bubble-PLs source (do not edit the function body). Reads lines of
// "POS<TAB>REF<TAB>ALT" on stdin; writes "POS<TAB>REF<TAB>ALT" of the minimal
// representation on stdout. Empty REF/ALT after trimming is printed as "-" so
// the diff against the Python side is unambiguous.
//
// Build (on the VM, rustc present):  rustc -O minrep_cli.rs -o minrep_cli
// (no crates needed; pure std)

use std::io::{self, BufRead, Write, BufWriter};

// ===== verbatim from extract-bubble-PLs (src/main.rs) — keep byte-identical =====
fn get_minimal_representation<'a>(mut pos: i64, mut r: &'a [u8], mut a: &'a [u8]) -> (i64, &'a [u8], &'a [u8]) {
    if !a.is_empty() && a[0] == b'<' {
        return (pos, r, a);
    }
    while !r.is_empty() && !a.is_empty() && r.last() == a.last() {
        r = &r[..r.len() - 1];
        a = &a[..a.len() - 1];
    }
    while !r.is_empty() && !a.is_empty() && r.first() == a.first() {
        r = &r[1..];
        a = &a[1..];
        pos += 1;
    }
    (pos, r, a)
}
// ================================================================================

fn main() {
    let stdin = io::stdin();
    let mut out = BufWriter::new(io::stdout().lock());
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        if line.is_empty() { continue; }
        let mut it = line.splitn(3, '\t');
        let pos: i64 = it.next().unwrap().parse().unwrap();
        let r = it.next().unwrap().as_bytes();
        let a = it.next().unwrap().as_bytes();
        let (p, mr, ma) = get_minimal_representation(pos, r, a);
        let mr = if mr.is_empty() { b"-" as &[u8] } else { mr };
        let ma = if ma.is_empty() { b"-" as &[u8] } else { ma };
        out.write_all(p.to_string().as_bytes()).unwrap();
        out.write_all(b"\t").unwrap();
        out.write_all(mr).unwrap();
        out.write_all(b"\t").unwrap();
        out.write_all(ma).unwrap();
        out.write_all(b"\n").unwrap();
    }
}

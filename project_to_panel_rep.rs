// =============================================================================
// project_to_panel_rep.rs -- Rust port of project_to_panel_rep.py (GT-only).
//
// Projects a hard-GT target VCF onto the PANEL's allele representation via
// minimal-representation matching, byte-faithful to the Python tool:
//   bcftools view <acaf.vcf[.gz]> | project_to_panel_rep <panel_sites.tsv> > projected.vcf
//   (panel_sites = CHROM<TAB>POS<TAB>REF<TAB>ALT, plain text; gz not supported here)
// Output is NOT position-sorted -> pipe through `bcftools sort`. Pure std (build:
// rustc -O project_to_panel_rep.rs -o project_to_panel_rep).
// =============================================================================
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

// ===== verbatim from extract-bubble-PLs / minrep_cli (keep byte-identical) =====
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

fn parse_i64(b: &[u8]) -> Option<i64> {
    std::str::from_utf8(b).ok()?.trim().parse().ok()
}

// Recode a GT string onto focal ALT index k (mirrors `bcftools norm -m -any`):
// 0->0, k->1, any other ALT->0, missing->. ; always UNPHASED ('/').
fn recode_gt(g: &[u8], k: i64, out: &mut Vec<u8>) {
    out.clear();
    let n = g.len();
    if n == 0 {
        out.push(b'.');
        return;
    }
    let mut first = true;
    let mut i = 0usize;
    loop {
        let mut j = i;
        while j < n && g[j] != b'/' && g[j] != b'|' {
            j += 1;
        }
        let tok = &g[i..j];
        if !first {
            out.push(b'/');
        }
        first = false;
        if tok.is_empty() || tok == b"." {
            out.push(b'.');
        } else {
            match parse_i64(tok) {
                Some(v) => out.push(if v == k { b'1' } else { b'0' }),
                None => out.push(b'.'),
            }
        }
        if j >= n {
            break;
        }
        i = j + 1;
    }
}

struct Interner {
    map: HashMap<Vec<u8>, u32>,
}
impl Interner {
    fn new() -> Self {
        Interner { map: HashMap::new() }
    }
    fn intern(&mut self, c: &[u8]) -> u32 {
        if let Some(&id) = self.map.get(c) {
            return id;
        }
        let id = self.map.len() as u32;
        self.map.insert(c.to_vec(), id);
        id
    }
    fn get(&self, c: &[u8]) -> Option<u32> {
        self.map.get(c).copied()
    }
}

#[derive(PartialEq, Eq, Hash)]
struct Key {
    chr: u32,
    pos: i64,
    r: u8,
    a: u8,
}
struct Val {
    pos: i64,
    r: Box<[u8]>,
    a: Box<[u8]>,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: ... | {} <panel_sites.tsv>", args[0]);
        std::process::exit(1);
    }

    // --- build panel minimal-rep -> raw map (SNV-minrep keys only; prefer exact) ---
    let mut interner = Interner::new();
    let mut pmap: HashMap<Key, Val> = HashMap::new();
    {
        let f = File::open(&args[1]).unwrap_or_else(|e| {
            eprintln!("ERROR opening {}: {}", args[1], e);
            std::process::exit(1);
        });
        let mut fh = BufReader::new(f);
        let mut buf: Vec<u8> = Vec::new();
        loop {
            buf.clear();
            if fh.read_until(b'\n', &mut buf).unwrap() == 0 {
                break;
            }
            while matches!(buf.last(), Some(b'\n') | Some(b'\r')) {
                buf.pop();
            }
            if buf.is_empty() || buf[0] == b'#' {
                continue;
            }
            // split first 4 tab fields: CHROM POS REF ALT
            let mut t = [0usize; 3];
            let mut nt = 0;
            for (i, &b) in buf.iter().enumerate() {
                if b == b'\t' {
                    if nt < 3 {
                        t[nt] = i;
                    }
                    nt += 1;
                    if nt == 3 {
                        break;
                    }
                }
            }
            if nt < 3 {
                continue;
            }
            let c = &buf[..t[0]];
            let p = match parse_i64(&buf[t[0] + 1..t[1]]) {
                Some(v) => v,
                None => continue,
            };
            let r = &buf[t[1] + 1..t[2]];
            // ALT field runs to the next tab (col 4) or end of line
            let alt_end = buf[t[2] + 1..].iter().position(|&b| b == b'\t').map(|x| t[2] + 1 + x).unwrap_or(buf.len());
            let alt_field = &buf[t[2] + 1..alt_end];
            let chr_id = interner.intern(c);
            for alt in alt_field.split(|&b| b == b',') {
                let (mp, mr, ma) = get_minimal_representation(p, r, alt);
                if mr.len() != 1 || ma.len() != 1 {
                    continue; // keep only SNV-minrep keys
                }
                let key = Key { chr: chr_id, pos: mp, r: mr[0], a: ma[0] };
                let is_exact = p == mp && r == mr && alt == ma;
                match pmap.entry(key) {
                    Entry::Vacant(slot) => {
                        slot.insert(Val { pos: p, r: r.into(), a: alt.into() });
                    }
                    Entry::Occupied(mut slot) => {
                        if is_exact {
                            slot.insert(Val { pos: p, r: r.into(), a: alt.into() });
                        }
                    }
                }
            }
        }
    }
    eprintln!("[project] panel SNV-minrep keys: {}", pmap.len());

    // --- stream stdin VCF, project matched alleles onto panel coords ---
    let stdin = io::stdin();
    let mut rdr = stdin.lock();
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let mut emitted: HashSet<(u32, i64, Box<[u8]>, Box<[u8]>)> = HashSet::new();
    let (mut n_in, mut n_exact, mut n_recov, mut n_drop): (u64, u64, u64, u64) = (0, 0, 0, 0);
    let mut line: Vec<u8> = Vec::new();
    let mut gtbuf: Vec<u8> = Vec::new();

    loop {
        line.clear();
        if rdr.read_until(b'\n', &mut line).unwrap() == 0 {
            break;
        }
        if !line.is_empty() && line[0] == b'#' {
            out.write_all(&line).unwrap(); // passthrough (keeps newline)
            continue;
        }
        let mut end = line.len();
        while end > 0 && (line[end - 1] == b'\n' || line[end - 1] == b'\r') {
            end -= 1;
        }
        let rec = &line[..end];
        if rec.is_empty() {
            continue;
        }
        // index the first 9 tabs (CHROM..FORMAT, then samples)
        let mut tabs = [0usize; 9];
        let mut nt = 0;
        for (i, &b) in rec.iter().enumerate() {
            if b == b'\t' {
                if nt < 9 {
                    tabs[nt] = i;
                }
                nt += 1;
                if nt == 9 {
                    break;
                }
            }
        }
        if nt < 9 {
            continue; // need 9 fixed cols + >=1 sample
        }
        let c = &rec[..tabs[0]];
        let p = match parse_i64(&rec[tabs[0] + 1..tabs[1]]) {
            Some(v) => v,
            None => continue,
        };
        let ref_b = &rec[tabs[2] + 1..tabs[3]];
        let alt_field = &rec[tabs[3] + 1..tabs[4]];
        let fmt = &rec[tabs[7] + 1..tabs[8]];
        let samples = &rec[tabs[8] + 1..];
        let gi = fmt.split(|&b| b == b':').position(|tg| tg == b"GT").unwrap_or(0);
        let chr_id = interner.get(c);

        for (k0, alt) in alt_field.split(|&b| b == b',').enumerate() {
            let k = (k0 + 1) as i64;
            n_in += 1;
            let (mp, mr, ma) = get_minimal_representation(p, ref_b, alt);
            if mr.len() != 1 || ma.len() != 1 {
                n_drop += 1;
                continue;
            }
            let cid = match chr_id {
                Some(id) => id,
                None => {
                    n_drop += 1;
                    continue;
                }
            };
            let praw = match pmap.get(&Key { chr: cid, pos: mp, r: mr[0], a: ma[0] }) {
                Some(v) => v,
                None => {
                    n_drop += 1;
                    continue;
                }
            };
            let ekey = (cid, praw.pos, praw.r.clone(), praw.a.clone());
            if emitted.contains(&ekey) {
                continue; // dedupe collisions onto the same panel allele
            }
            emitted.insert(ekey);
            if praw.pos == p && praw.r.as_ref() == ref_b && praw.a.as_ref() == alt {
                n_exact += 1;
            } else {
                n_recov += 1;
            }
            // emit at panel raw coords, GT-only: CHROM POS . REF ALT . . . GT gts...
            out.write_all(c).unwrap();
            write!(out, "\t{}\t.\t", praw.pos).unwrap();
            out.write_all(&praw.r).unwrap();
            out.write_all(b"\t").unwrap();
            out.write_all(&praw.a).unwrap();
            out.write_all(b"\t.\t.\t.\tGT").unwrap();
            for sample in samples.split(|&b| b == b'\t') {
                let s_gt: &[u8] = if sample.is_empty() || sample == b"." {
                    &b"./."[..]
                } else if gi == 0 {
                    match sample.iter().position(|&b| b == b':') {
                        Some(p2) => &sample[..p2],
                        None => sample,
                    }
                } else {
                    sample.split(|&b| b == b':').nth(gi).unwrap_or(&b"."[..])
                };
                recode_gt(s_gt, k, &mut gtbuf);
                out.write_all(b"\t").unwrap();
                out.write_all(&gtbuf).unwrap();
            }
            out.write_all(b"\n").unwrap();
        }
    }
    out.flush().unwrap();
    eprintln!(
        "[project] in(alleles)={} exact={} recovered={} dropped(no panel match)={}",
        n_in, n_exact, n_recov, n_drop
    );
}

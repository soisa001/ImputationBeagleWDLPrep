use flate2::read::MultiGzDecoder; 
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, BufWriter};

struct Record {
    info_map: HashMap<String, String>,
    atomic_ids: HashSet<String>,
    gts: Vec<String>,
    gps: Vec<Vec<f32>>,
}

#[derive(Clone)]
struct PhasedDist {
    p0: f32, 
    p1: f32, 
}

fn smart_open(filename: &str) -> Box<dyn BufRead> {
    let file = File::open(filename).expect("Cannot open file");
    if filename.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file))) 
    } else {
        Box::new(BufReader::new(file))
    }
}

fn parse_info(info_str: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for item in info_str.split(';') {
        if let Some((k, v)) = item.split_once('=') {
            map.insert(k.to_string(), v.to_string());
        }
    }
    map
}

fn format_float(val: f32) -> String {
    let s = format!("{:.3}", val);
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

fn process_group(
    group_lines: &[String],
    id_buffer: &HashMap<String, (u32, String, String, String)>,
    max_alleles: usize,
    out_handle: &mut impl Write,
) {
    if group_lines.is_empty() { return; }

    let parsed_lines: Vec<Vec<&str>> = group_lines.iter()
        .map(|line| line.trim_end().split('\t').collect())
        .collect();

    if parsed_lines[0].len() <= 9 {
        panic!("Error: VCF does not contain sample columns. This script requires sample Genotype (GT) and Probability (GP) columns to project joint distributions.");
    }

    let num_samples = parsed_lines[0].len() - 9;
    let chrom = parsed_lines[0][0].to_string();
    let num_alleles = parsed_lines.len();

    let mut records: Vec<Record> = Vec::with_capacity(num_alleles);
    let mut all_atomic_ids = HashSet::new();

    for fields in &parsed_lines {
        let info_map = parse_info(fields[7]);
        
        let mut atomic_ids = HashSet::new();
        if let Some(id_str) = info_map.get("ID") {
            let replaced_id = id_str.replace(',', ":");
            for j in replaced_id.split(':').map(|s| s.trim()) {
                if id_buffer.contains_key(j) {
                    atomic_ids.insert(j.to_string());
                    all_atomic_ids.insert(j.to_string());
                }
            }
        }

        let fmt: Vec<&str> = fields[8].split(':').collect();
        let gt_idx = fmt.iter().position(|&x| x == "GT");
        let gp_idx = fmt.iter().position(|&x| x == "GP");

        let mut gts = Vec::with_capacity(num_samples);
        let mut gps = Vec::with_capacity(num_samples);

        for s in 0..num_samples {
            let s_vals: Vec<&str> = fields[9 + s].split(':').collect();
            
            let gt = if let Some(idx) = gt_idx {
                if idx < s_vals.len() && s_vals[idx] != "." { s_vals[idx].to_string() } else { "0|0".to_string() }
            } else { "0|0".to_string() };

            let mut sample_gp = vec![1.0, 0.0, 0.0];
            if let Some(idx) = gp_idx {
                if idx < s_vals.len() && s_vals[idx] != "." {
                    let parsed: Vec<f32> = s_vals[idx].split(',')
                        .filter_map(|x| x.parse::<f32>().ok())
                        .collect();
                    if parsed.len() == 3 { sample_gp = parsed; }
                }
            }

            gts.push(gt);
            gps.push(sample_gp);
        }

        records.push(Record {
            info_map,
            atomic_ids,
            gts,
            gps,
        });
    }

    // Step 1: Natively derive Hap0 and Hap1 marginals directly from GP arrays
    let mut hap_probs: Vec<Vec<(f32, f32)>> = vec![vec![(0.0, 0.0); num_alleles]; num_samples];
    
    for s in 0..num_samples {
        for a in 0..num_alleles {
            let gt = &records[a].gts[s];
            let gp = &records[a].gps[s];
            
            let gp1 = gp[1];
            let gp2 = gp[2];
            
            let (p0, p1) = if gt.starts_with("1|0") {
                (gp2 + gp1, gp2)
            } else if gt.starts_with("0|1") {
                (gp2, gp2 + gp1)
            } else {
                (gp2 + (gp1 / 2.0), gp2 + (gp1 / 2.0))
            };
            
            // Apply un-rounding clamp to strictly prevent division by zero in the odds ratio
            // 1e-5 represents a maximum internal confidence cap of 99.999%
            hap_probs[s][a] = (p0.clamp(1e-5, 1.0 - 1e-5), p1.clamp(1e-5, 1.0 - 1e-5));
        }
    }

    // Step 2 & 3: Truncate, convert to Odds, and Normalize
    let mut sample_atomic_dists: Vec<HashMap<String, PhasedDist>> = vec![HashMap::new(); num_samples];

    for s in 0..num_samples {
        let mut allele_scores: Vec<(usize, f32)> = Vec::with_capacity(num_alleles);
        for a in 0..num_alleles {
            let score = hap_probs[s][a].0 + hap_probs[s][a].1;
            allele_scores.push((a, score));
        }

        allele_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let m = allele_scores.len().min(max_alleles);
        let top_alleles: Vec<usize> = allele_scores.iter().take(m).map(|x| x.0).collect();

        let mut z0 = 1.0_f32; // Weight for Reference allele is always 1.0
        let mut z1 = 1.0_f32;
        
        let mut w0_map = HashMap::new();
        let mut w1_map = HashMap::new();

        for &a in &top_alleles {
            let (p0, p1) = hap_probs[s][a];
            let w0 = p0 / (1.0 - p0);
            let w1 = p1 / (1.0 - p1);
            
            w0_map.insert(a, w0);
            w1_map.insert(a, w1);
            
            z0 += w0;
            z1 += w1;
        }

        // Map the normalized haplotypic probabilities down to the atomic variants
        for &a in &top_alleles {
            let norm_p0 = w0_map[&a] / z0;
            let norm_p1 = w1_map[&a] / z1;

            for atomic_id in &records[a].atomic_ids {
                let entry = sample_atomic_dists[s].entry(atomic_id.clone()).or_insert(PhasedDist { p0: 0.0, p1: 0.0 });
                entry.p0 += norm_p0;
                entry.p1 += norm_p1;
            }
        }
    }

    let mut sorted_atomic_vars: Vec<_> = all_atomic_ids.into_iter().collect();
    
    // STRICT SORTING FIX: Sort by POS, then REF, then ALT to ensure deterministic output
    sorted_atomic_vars.sort_by(|a, b| {
        let data_a = id_buffer.get(a);
        let data_b = id_buffer.get(b);
        
        match (data_a, data_b) {
            (Some(da), Some(db)) => {
                da.0.cmp(&db.0)                      // 1. Compare POS (coord)
                    .then_with(|| da.2.cmp(&db.2))   // 2. Compare REF string
                    .then_with(|| da.3.cmp(&db.3))   // 3. Compare ALT string
            }
            _ => std::cmp::Ordering::Equal,
        }
    });

    // Step 4: Derive Final VCF Fields
    for assigned_id in sorted_atomic_vars {
        let var_data = &id_buffer[&assigned_id];
        let coord = var_data.0;
        
        let t_rec = records.iter().find(|r| r.atomic_ids.contains(&assigned_id)).unwrap();
        let mut new_info = vec![format!("ID={}", assigned_id)];
        for k in ["MA", "UK", "RAF", "AF", "INFO"] {
            if let Some(v) = t_rec.info_map.get(k) {
                new_info.push(format!("{}={}", k, v));
            }
        }

        let mut vcf_line = vec![
            chrom.clone(),
            coord.to_string(),
            var_data.1.clone(), 
            var_data.2.clone(), 
            var_data.3.clone(), 
            ".".to_string(),
            ".".to_string(),
            new_info.join(";"),
            "GT:DS:GP".to_string(),
        ];

        for s in 0..num_samples {
            let empty_dist = PhasedDist { p0: 0.0, p1: 0.0 };
            let dist = sample_atomic_dists[s].get(&assigned_id).unwrap_or(&empty_dist);
            
            let p0 = dist.p0.clamp(0.0, 1.0);
            let p1 = dist.p1.clamp(0.0, 1.0);

            let hap0_gt = if p0 > 0.5 { "1" } else { "0" };
            let hap1_gt = if p1 > 0.5 { "1" } else { "0" };
            let gt = format!("{}|{}", hap0_gt, hap1_gt);

            let ds_str = format_float(p0 + p1);

            let gp0_raw = (1.0 - p0) * (1.0 - p1);
            let gp1_raw = p0 * (1.0 - p1) + (1.0 - p0) * p1;
            let gp2_raw = p0 * p1;
            
            let mut v0 = (gp0_raw * 1000.0).round() as i32;
            let mut v1 = (gp1_raw * 1000.0).round() as i32;
            let mut v2 = (gp2_raw * 1000.0).round() as i32;
            
            let diff = 1000 - (v0 + v1 + v2);
            if diff != 0 {
                if v0 >= v1 && v0 >= v2 { v0 += diff; }
                else if v1 >= v0 && v1 >= v2 { v1 += diff; }
                else { v2 += diff; }
            }
            let gp_str = format!("{},{},{}", 
                format_float(v0 as f32 / 1000.0), 
                format_float(v1 as f32 / 1000.0), 
                format_float(v2 as f32 / 1000.0)
            );

            vcf_line.push(format!("{}:{}:{}", gt, ds_str, gp_str));
        }
        writeln!(out_handle, "{}", vcf_line.join("\t")).unwrap();
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: cat <multiallelic VCF> | {} <biallelic ID VCF> [max_alleles]", args[0]);
        std::process::exit(1);
    }
    
    let vcf_path = &args[1];
    let max_alleles: usize = if args.len() > 2 {
        args[2].parse().unwrap_or(10)
    } else {
        10
    };

    let mut id_iter = smart_open(vcf_path).lines().peekable();
    
    let stdout = io::stdout();
    let mut out_handle = BufWriter::new(stdout.lock());

    let stdin = io::stdin();
    let mut current_pos: Option<String> = None;
    let mut current_chrom: String = String::new();
    let mut group: Vec<String> = Vec::new(); 
    let mut records_processed: usize = 0;

    let mut id_buffer: HashMap<String, (u32, String, String, String)> = HashMap::new();
    let mut active_id_chrom = String::new();

    eprintln!("Starting Phased Joint-Distribution projection (Max Alleles: {})...", max_alleles);
    for line_result in stdin.lock().lines() {
        let line = line_result.unwrap();
        if line.starts_with('#') {
            if !line.contains("INFO=<ID=AK") && !line.contains("FORMAT=<ID=GL") && !line.contains("FORMAT=<ID=KC") {
                writeln!(out_handle, "{}", line).unwrap();
            }
            continue;
        }

        let fields: Vec<&str> = line.splitn(3, '\t').collect();
        if fields.len() < 2 { continue; }
        
        let chrom = fields[0].to_string();
        let pos = fields[1].to_string();

        records_processed += 1;
        if records_processed % 10_000 == 0 {
            eprintln!("Processed {} input records... (Currently at {}:{})", records_processed, chrom, pos);
        }

        if current_pos.is_none() {
            current_pos = Some(pos.clone());
            current_chrom = chrom.clone();
        }

        if pos != *current_pos.as_ref().unwrap() || chrom != current_chrom {
            if current_chrom != active_id_chrom {
                id_buffer.clear();
                active_id_chrom = current_chrom.clone();
                
                while let Some(Ok(peek_line)) = id_iter.peek() {
                    if peek_line.starts_with('#') {
                        id_iter.next();
                        continue;
                    }
                    let peek_fields: Vec<&str> = peek_line.splitn(3, '\t').collect();
                    let peek_chrom = peek_fields[0];
                    
                    if peek_chrom != active_id_chrom {
                        if id_buffer.is_empty() {
                            id_iter.next(); 
                            continue;
                        } else {
                            break; 
                        }
                    }
                    
                    let pop_line = id_iter.next().unwrap().unwrap();
                    let pop_fields: Vec<&str> = pop_line.split('\t').collect();
                    let pop_pos: u32 = pop_fields[1].parse().unwrap_or(0);
                    let orig_id = pop_fields[2].to_string();
                    let ref_seq = pop_fields[3].to_string();
                    let alt_seq = pop_fields[4].to_string();
                    
                    for item in pop_fields[7].split(';') {
                        if let Some(id_val) = item.strip_prefix("ID=") {
                            id_buffer.insert(id_val.to_string(), (pop_pos, orig_id, ref_seq, alt_seq));
                            break;
                        }
                    }
                }
            }

            process_group(&group, &id_buffer, max_alleles, &mut out_handle);
            group.clear();
            
            current_pos = Some(pos);
            current_chrom = chrom;
        }
        group.push(line);
    }

    // note that this might not be safe for edge cases (single-variant in the VCF, or last record being a single variant on a new chromosome)
    if !group.is_empty() {
        process_group(&group, &id_buffer, max_alleles, &mut out_handle);
    }
    
    out_handle.flush().unwrap();
    eprintln!("Finished! Processed a total of {} input records.", records_processed);
}

//! Parse AlphaFold Server (AlphaFold 3) results.
//! af3_report <directory of unpacked results> <reference.fa>
//!
//! For each job prints:
//!   - overall metrics (ipTM, pTM, ranking score, clashes, fraction disordered);
//!   - per-chain-pair ipTM matrix with human-readable chain names;
//!   - inter-chain contacts (residues closer than 4.5 Å) for pairs involving TMCO6/ZNF511.
//! Chain names are assigned by best-match of chain sequence against reference.fa;
//! when the match is imperfect the number of mismatches is appended to the name.
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

fn three_to_one(s: &str) -> char {
    match s {
        "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D', "CYS" => 'C', "GLN" => 'Q',
        "GLU" => 'E', "GLY" => 'G', "HIS" => 'H', "ILE" => 'I', "LEU" => 'L', "LYS" => 'K',
        "MET" => 'M', "PHE" => 'F', "PRO" => 'P', "SER" => 'S', "THR" => 'T', "TRP" => 'W',
        "TYR" => 'Y', "VAL" => 'V', "A" => 'a', "C" => 'c', "G" => 'g', "U" => 'u', _ => 'X',
    }
}

struct Atom {
    chain: String,
    seq: i64,
    comp: String,
    xyz: [f64; 3],
    plddt: f64,
}

/// mmCIF atom_site loop parser: columns are looked up by name, not by position.
fn parse_cif(path: &Path) -> Vec<Atom> {
    let txt = fs::read_to_string(path).unwrap();
    let mut cols: Vec<String> = Vec::new();
    let mut out = Vec::new();
    let mut in_loop = false;
    for line in txt.lines() {
        let t = line.trim();
        if t.starts_with("_atom_site.") {
            cols.push(t.trim_start_matches("_atom_site.").to_string());
            in_loop = true;
            continue;
        }
        if in_loop {
            if t.starts_with('#') || t.starts_with("loop_") || t.is_empty() {
                if !out.is_empty() {
                    break;
                }
                continue;
            }
            let f: Vec<&str> = t.split_whitespace().collect();
            if f.len() < cols.len() {
                continue;
            }
            let g = |n: &str| cols.iter().position(|c| c == n).map(|i| f[i]).unwrap_or("");
            if g("group_PDB") != "ATOM" && g("group_PDB") != "HETATM" {
                continue;
            }
            if g("pdbx_PDB_model_num") != "1" && !g("pdbx_PDB_model_num").is_empty() {
                continue;
            }
            let num = |n: &str| g(n).parse::<f64>().unwrap_or(0.0);
            out.push(Atom {
                chain: g("auth_asym_id").to_string(),
                seq: g("label_seq_id").parse().unwrap_or(0),
                comp: g("label_comp_id").to_string(),
                xyz: [num("Cartn_x"), num("Cartn_y"), num("Cartn_z")],
                plddt: num("B_iso_or_equiv"),
            });
        }
    }
    out
}

fn read_fasta(p: &str) -> Vec<(String, String)> {
    fs::read_to_string(p).unwrap().split('>').skip(1).map(|r| {
        let mut l = r.lines();
        let h = l.next().unwrap();
        let name = h.split('|').nth(2).and_then(|x| x.split('_').next()).unwrap_or(h).to_string();
        (name, l.collect::<String>())
    }).collect()
}

/// Minimum number of mismatches when lengths are equal; otherwise None.
fn mismatches(a: &str, b: &str) -> Option<usize> {
    (a.len() == b.len()).then(|| a.bytes().zip(b.bytes()).filter(|(x, y)| x != y).count())
}

fn jnum(v: &serde_json::Value, k: &str) -> String {
    match &v[k] {
        serde_json::Value::Number(n) => format!("{:.3}", n.as_f64().unwrap()),
        serde_json::Value::Bool(b) => b.to_string(),
        serde_json::Value::Null => "—".into(),
        x => x.to_string(),
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let refs = read_fasta(&a[2]);
    // collect jobs: key is the filename prefix up to _summary_confidences
    let mut jobs: BTreeMap<String, (Option<PathBuf>, Option<PathBuf>)> = BTreeMap::new();
    let mut stack = vec![PathBuf::from(&a[1])];
    while let Some(d) = stack.pop() {
        for e in fs::read_dir(&d).unwrap().flatten() {
            let p = e.path();
            if p.is_dir() {
                stack.push(p);
                continue;
            }
            // take only the best model (rank 0)
            let n = p.file_name().unwrap().to_string_lossy().to_string();
            if let Some(key) = n.strip_suffix("_model_0.cif") {
                jobs.entry(key.to_string()).or_default().0 = Some(p);
            } else if let Some(key) = n.strip_suffix("_summary_confidences_0.json") {
                jobs.entry(key.to_string()).or_default().1 = Some(p);
            }
        }
    }
    for (job, (cif, conf)) in &jobs {
        println!("\n===== {job}");
        let Some(cif) = cif else {
            println!("  no .cif found");
            continue;
        };
        let atoms = parse_cif(cif);
        // sequence of each chain
        let mut chains: BTreeMap<String, Vec<(i64, String)>> = BTreeMap::new();
        for at in &atoms {
            let e = chains.entry(at.chain.clone()).or_default();
            if e.last().map(|x| x.0) != Some(at.seq) {
                e.push((at.seq, at.comp.clone()));
            }
        }
        let mut name: BTreeMap<String, String> = BTreeMap::new();
        for (c, res) in &chains {
            let s: String = res.iter().map(|r| three_to_one(&r.1)).collect();
            let best = refs.iter().filter_map(|(n, rs)| mismatches(&s, rs).map(|m| (m, n))).min();
            name.insert(c.clone(), match best {
                Some((0, n)) => n.clone(),
                Some((m, n)) => format!("{n}[{m} sub]"),
                None => format!("chain {c} ({} res.)", s.len()),
            });
        }
        println!("  chains: {}", name.iter().map(|(c, n)| format!("{c}={n}")).collect::<Vec<_>>().join(", "));

        if let Some(conf) = conf {
            let v: serde_json::Value = serde_json::from_str(&fs::read_to_string(conf).unwrap()).unwrap();
            println!("  ipTM {} | pTM {} | ranking {} | clash {} | disordered {}",
                jnum(&v, "iptm"), jnum(&v, "ptm"), jnum(&v, "ranking_score"), jnum(&v, "has_clash"), jnum(&v, "fraction_disordered"));
            // chain ids are taken from the JSON file, not from the structural order
            let ids: Vec<String> = v["chain_ids"].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap_or("?").to_string()).collect())
                .unwrap_or_else(|| chains.keys().cloned().collect());
            if let Some(m) = v["chain_pair_iptm"].as_array() {
                let pae = v["chain_pair_pae_min"].as_array();
                let mut rows: Vec<(f64, String)> = Vec::new();
                for (i, row) in m.iter().enumerate() {
                    for (j, x) in row.as_array().unwrap().iter().enumerate() {
                        if j > i {
                            let unk = format!("chain {}", ids[i]);
                            let ni = name.get(&ids[i]).unwrap_or(&unk).clone();
                            let nj = name.get(&ids[j]).cloned().unwrap_or_else(|| format!("chain {}", ids[j]));
                            let iptm = x.as_f64().unwrap_or(f64::NAN);
                            let pm = pae.and_then(|p| p[i][j].as_f64()).unwrap_or(f64::NAN);
                            rows.push((iptm, format!("    {ni} — {nj}\tipTM {iptm:.2}\tPAE_min {pm:.1}")));
                        }
                    }
                }
                rows.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
                println!("  chain pairs (descending ipTM):");
                for (_, s) in &rows {
                    println!("{s}");
                }
            }
        }
        // inter-chain contacts
        let keys: Vec<String> = chains.keys().cloned().collect();
        for i in 0..keys.len() {
            for j in (i + 1)..keys.len() {
                let (ca, cb) = (&keys[i], &keys[j]);
                // by default only pairs involving TMCO6/ZNF511; ALL_PAIRS=1 to show all
                let all = std::env::var("ALL_PAIRS").is_ok();
                let want = |c: &String| all || name[c].starts_with("TMCO6") || name[c].starts_with("ZN511");
                if !want(ca) && !want(cb) {
                    continue;
                }
                let aa: Vec<&Atom> = atoms.iter().filter(|x| &x.chain == ca).collect();
                let bb: Vec<&Atom> = atoms.iter().filter(|x| &x.chain == cb).collect();
                let mut ra: BTreeMap<i64, (String, f64)> = BTreeMap::new();
                let mut rb: BTreeMap<i64, (String, f64)> = BTreeMap::new();
                for x in &aa {
                    for y in &bb {
                        let d = (0..3).map(|k| (x.xyz[k] - y.xyz[k]).powi(2)).sum::<f64>();
                        if d < 20.25 {
                            ra.entry(x.seq).or_insert((x.comp.clone(), x.plddt));
                            rb.entry(y.seq).or_insert((y.comp.clone(), y.plddt));
                        }
                    }
                }
                if ra.is_empty() {
                    continue;
                }
                let fmt = |m: &BTreeMap<i64, (String, f64)>| m.iter().map(|(k, v)| format!("{}{k}", three_to_one(&v.0))).collect::<Vec<_>>().join(" ");
                let pl = |m: &BTreeMap<i64, (String, f64)>| m.values().map(|v| v.1).sum::<f64>() / m.len() as f64;
                println!("  contact {} ({} res., pLDDT {:.0}) — {} ({} res., pLDDT {:.0})", name[ca], ra.len(), pl(&ra), name[cb], rb.len(), pl(&rb));
                println!("    {}: {}", name[ca], fmt(&ra));
                println!("    {}: {}", name[cb], fmt(&rb));
            }
        }
    }
}

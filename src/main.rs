//! Профиль данных для малоизученных (Tdark) белков: публикации (NCBI gene2pubmed), DepMap (Chronos),
//! BioPlex 3.0 (293T + HCT116), локализация HPA. Пишет TSV и печатает лучших кандидатов.
use flate2::read::GzDecoder;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (depmap, g2p) = (&a[1], &a[2]);
    // 1. Tdark
    let tdark: Vec<serde_json::Value> = serde_json::from_str(&std::fs::read_to_string("data/pharos_tdark.json").unwrap()).unwrap();
    let tsyms: HashSet<String> = tdark.iter().map(|t| t["sym"].as_str().unwrap_or("").to_string()).collect();
    let tname: HashMap<String, String> = tdark.iter().map(|t| (t["sym"].as_str().unwrap_or("").to_string(), t["name"].as_str().unwrap_or("").to_string())).collect();
    // 2. символ -> GeneID
    let mut sym2id: HashMap<String, String> = HashMap::new();
    let mut band: HashMap<String, String> = HashMap::new();
    for l in BufReader::new(GzDecoder::new(std::fs::File::open("data/Homo_sapiens.gene_info.gz").unwrap())).lines().skip(1) {
        let l = l.unwrap(); let f: Vec<&str> = l.split('\t').collect();
        band.insert(f[2].to_string(), f[7].split('.').next().unwrap_or("").to_string());
        if tsyms.contains(f[2]) { sym2id.insert(f[2].to_string(), f[1].to_string()); }
    }
    let ids: HashSet<String> = sym2id.values().cloned().collect();
    // 3. публикации
    let mut pubs: HashMap<String, u32> = HashMap::new();
    for l in BufReader::new(GzDecoder::new(std::fs::File::open(g2p).unwrap())).lines().skip(1) {
        let l = l.unwrap(); let mut f = l.split('\t');
        if f.next() != Some("9606") { continue; }
        let g = f.next().unwrap();
        if ids.contains(g) { *pubs.entry(g.to_string()).or_default() += 1; }
    }
    // 4. DepMap: колонки «SYM (GeneID)»
    let mut lines = BufReader::new(std::fs::File::open(depmap).unwrap()).lines();
    let header: Vec<String> = lines.next().unwrap().unwrap().split(',').map(String::from).collect();
    let all_syms: Vec<String> = header.iter().skip(1).map(|h| h.split(' ').next().unwrap().to_string()).collect();
    let mut mat: Vec<Vec<f64>> = vec![vec![]; all_syms.len()]; // ген -> значения по линиям (NaN сохраняем)
    for l in lines {
        let l = l.unwrap(); let f: Vec<&str> = l.split(',').collect();
        for j in 0..all_syms.len() { mat[j].push(f[j + 1].parse::<f64>().unwrap_or(f64::NAN)); }
    }
    let want: Vec<(usize, String)> = all_syms.iter().enumerate().filter(|(_, s)| tsyms.contains(*s)).map(|(i, s)| (i, s.clone())).collect();
    let vals: Vec<Vec<f64>> = want.iter().map(|(i, _)| mat[*i].iter().cloned().filter(|x| !x.is_nan()).collect()).collect();
    // Стандартизованная матрица для ко-зависимости (NaN -> 0 после центрирования), только гены с ≥ 500 линиями
    let nl = mat[0].len();
    let z: Vec<Option<Vec<f32>>> = mat.iter().map(|v| {
        let ok: Vec<f64> = v.iter().cloned().filter(|x| !x.is_nan()).collect();
        if ok.len() < 500 { return None; }
        let m = ok.iter().sum::<f64>() / ok.len() as f64;
        let sd = (ok.iter().map(|x| (x - m).powi(2)).sum::<f64>() / ok.len() as f64).sqrt();
        if sd < 0.05 { return None; }
        Some(v.iter().map(|x| if x.is_nan() { 0.0 } else { ((x - m) / sd) as f32 }).collect())
    }).collect();
    use rayon::prelude::*;
    let multicopy = |s: &str| ["OR", "KRTAP", "GOLGA6", "GOLGA8", "AGAP", "POTE", "SPDYE", "PRAMEF", "NUTM2", "LRRC37A", "CTAGE", "TBC1D3", "USP17L",
        "FAM90A", "NBPF", "SPATA31", "ANKRD20A", "ANKRD30", "ANKRD36", "DEFB", "LCE", "GAGE", "XAGE", "PAGE", "MAGE", "SSX", "CT45", "CT47", "RGPD",
        "NPIPA", "NPIPB", "TP53TG3", "FAM230", "CCDC144", "TRIM64", "TRIM43", "TRIM49", "TRIM51", "TRIM61", "ZNF658", "SMN", "SERF1", "GTF2IRD2",
        "NCF1", "HYDIN", "FRG2", "DUX4", "MUC", "HSFX", "HSFY", "RBMY", "TSPY", "VCX", "VCY", "PRY", "CDY", "DAZ", "SPANX", "OFD1P", "FAM25", "OCLM"]
        .iter().any(|p| s.starts_with(p) && s[p.len()..].chars().next().map_or(false, |c| c.is_ascii_digit() || c.is_ascii_uppercase()));
    let coess: HashMap<String, Vec<(String, f32)>> = want.par_iter().filter_map(|(i, s)| {
        let zi = z[*i].as_ref()?;
        let mut best: Vec<(f32, usize)> = (0..all_syms.len()).filter(|j| *j != *i && band.get(&all_syms[*j]) != band.get(s) && z[*j].is_some())
            .map(|j| (zi.iter().zip(z[j].as_ref().unwrap()).map(|(a, b)| a * b).sum::<f32>() / nl as f32, j)).collect();
        best.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        Some((s.clone(), best.iter().take(5).map(|(r, j)| (all_syms[*j].clone(), *r)).collect()))
    }).collect();
    let mut dm: HashMap<String, (usize, f64, f64, f64)> = HashMap::new(); // n, среднее, доля < -0.5, минимум
    for (k, (_, s)) in want.iter().enumerate() {
        let v = &vals[k]; if v.is_empty() { continue; }
        let n = v.len();
        dm.insert(s.clone(), (n, v.iter().sum::<f64>() / n as f64, v.iter().filter(|x| **x < -0.5).count() as f64 / n as f64, v.iter().cloned().fold(f64::MAX, f64::min)));
    }
    // 5. BioPlex
    let mut partners: HashMap<String, HashSet<String>> = HashMap::new();
    for f in ["data/BioPlex_293T_Network_10K_Dec_2019.tsv", "data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv"] {
        for l in std::fs::read_to_string(f).unwrap().lines().skip(1) {
            let c: Vec<&str> = l.split('\t').collect();
            let (x, y) = (c[4].trim_matches('"'), c[5].trim_matches('"'));
            partners.entry(x.to_string()).or_default().insert(y.to_string());
            partners.entry(y.to_string()).or_default().insert(x.to_string());
        }
    }
    // 6. HPA
    let mut hpa: HashMap<String, (String, String)> = HashMap::new();
    for l in std::fs::read_to_string("data/subcellular_location.tsv").unwrap().lines().skip(1) {
        let c: Vec<&str> = l.split('\t').collect();
        hpa.insert(c[1].to_string(), (c[2].to_string(), c[3].to_string()));
    }
    // Сводка
    let mut out = String::from("sym\tname\tpubs\tdepmap_n\tdepmap_mean\tdepmap_frac_lt_-0.5\tdepmap_min\tbioplex_partners\thpa_reliability\thpa_location\tcoessential_top5\tscore\n");
    let mut rows = vec![];
    for s in &tsyms {
        let p = sym2id.get(s).and_then(|g| pubs.get(g)).cloned().unwrap_or(0);
        let d = dm.get(s).cloned();
        let bp = partners.get(s).map(|x| x.len()).unwrap_or(0);
        let (rel, loc) = hpa.get(s).cloned().unwrap_or_default();
        // Очки «насыщенности данными»
        let mut score = 0.0;
        if let Some((n, mean, frac, _)) = d { if n >= 500 { if frac >= 0.05 { score += 2.0; } if mean < -0.3 { score += 1.0; } } }
        if bp >= 3 { score += 1.0; } if bp >= 8 { score += 1.0; }
        let ce = coess.get(s).cloned().unwrap_or_default();
        let top_r = ce.first().map(|x| x.1).unwrap_or(0.0);
        if top_r >= 0.25 { score += 1.0; } if top_r >= 0.4 { score += 2.0; }
        if multicopy(s) { score = -1.0; }
        let ce_s = ce.iter().map(|(g, r)| format!("{g}:{r:.2}")).collect::<Vec<_>>().join(",");
        if rel == "Enhanced" || rel == "Supported" { score += 1.0; }
        let (n, mean, frac, mn) = d.unwrap_or((0, f64::NAN, f64::NAN, f64::NAN));
        out += &format!("{s}\t{}\t{p}\t{n}\t{mean:.3}\t{frac:.3}\t{mn:.2}\t{bp}\t{rel}\t{loc}\t{ce_s}\t{score}\n", tname[s]);
        rows.push((score, p, s.clone(), n, mean, frac, bp, rel, ce_s));
    }
    std::fs::write("data/tdark_profile.tsv", out).unwrap();
    rows.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
    println!("Tdark: {}, с GeneID: {}, в DepMap: {}, в BioPlex: {}, в HPA: {}", tsyms.len(), sym2id.len(), dm.len(),
        tsyms.iter().filter(|s| partners.contains_key(*s)).count(), tsyms.iter().filter(|s| hpa.contains_key(*s)).count());
    println!("{:<11} {:>4} {:>4} {:>6} {:>5} {:>3} {:<9} {}", "ген", "очки", "стат", "DM ср.", "<-.5", "BP", "HPA", "ко-зависимость (топ-3)");
    for (sc, p, s, _, mean, frac, bp, rel, ce) in rows.iter().filter(|r| r.1 <= 15).take(40) {
        let ce3 = ce.split(',').take(3).collect::<Vec<_>>().join(" ");
        println!("{s:<11} {sc:>4} {p:>4} {mean:>6.2} {:>4.0}% {bp:>3} {rel:<9} {ce3}", frac * 100.0);
    }
}

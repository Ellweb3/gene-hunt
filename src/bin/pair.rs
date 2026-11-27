//! Разбор пары TMCO6–ZNF511 по DepMap Chronos.
//! pair <gene_effect.csv> <ach_disease.tsv> <out_prefix>
//! 1) корреляции пары с панелью генов LSm/сплайсинга/распада мРНК;
//! 2) типы опухолей, где пара нужнее (Welch t против остальных);
//! 3) геномные корреляты среднего профиля пары (модуль).
use std::collections::HashMap;
use std::fs;
use std::io::Write;

fn pearson(a: &[f64], b: &[f64]) -> (f64, usize) {
    let pairs: Vec<(f64, f64)> = a.iter().zip(b).filter(|(x, y)| x.is_finite() && y.is_finite()).map(|(x, y)| (*x, *y)).collect();
    let n = pairs.len();
    if n < 50 {
        return (f64::NAN, n);
    }
    let mx = pairs.iter().map(|p| p.0).sum::<f64>() / n as f64;
    let my = pairs.iter().map(|p| p.1).sum::<f64>() / n as f64;
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (x, y) in &pairs {
        sxy += (x - mx) * (y - my);
        sxx += (x - mx).powi(2);
        syy += (y - my).powi(2);
    }
    (sxy / (sxx * syy).sqrt(), n)
}

fn mean_var(v: &[f64]) -> (f64, f64, usize) {
    let n = v.len();
    let m = v.iter().sum::<f64>() / n as f64;
    let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n as f64 - 1.0).max(1.0);
    (m, var, n)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (ge, dis, out) = (&a[1], &a[2], &a[3]);
    let txt = fs::read_to_string(ge).expect("gene_effect");
    let mut lines = txt.lines();
    let header: Vec<String> = lines.next().unwrap().split(',').skip(1).map(|s| s.split(" (").next().unwrap().trim_matches('"').to_string()).collect();
    let g = header.len();
    let mut ids = Vec::new();
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); g];
    for l in lines {
        let mut it = l.split(',');
        ids.push(it.next().unwrap().to_string());
        for (j, v) in it.enumerate() {
            cols[j].push(v.parse::<f64>().unwrap_or(f64::NAN));
        }
    }
    let idx: HashMap<&str, usize> = header.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let t = &cols[idx["TMCO6"]];
    let z = &cols[idx["ZNF511"]];
    let pairavg: Vec<f64> = t.iter().zip(z).map(|(x, y)| (x + y) / 2.0).collect();

    // 1. панель
    let panel = [
        "LSM1", "LSM2", "LSM3", "LSM4", "LSM5", "LSM6", "LSM7", "LSM8", "LSM10", "LSM11", "LSM12", "LSM14A", "PATL1", "DDX6",
        "SART3", "PRPF3", "PRPF4", "PRPF31", "PRPF8", "USB1", "SNRNP200", "EFTUD2", "XRN1", "DCP2", "TUT1", "MEPCE", "LARP7",
        "KPNA1", "KPNA2", "KPNA3", "KPNA4", "KPNA6", "KPNB1", "CSE1L", "TMCO6", "ZNF511",
    ];
    let mut f = fs::File::create(format!("{out}_panel.tsv")).unwrap();
    writeln!(f, "gene\tmean\tr_TMCO6\tr_ZNF511\tn").unwrap();
    for p in panel {
        if let Some(&j) = idx.get(p) {
            let (rt, n) = pearson(t, &cols[j]);
            let (rz, _) = pearson(z, &cols[j]);
            let fin: Vec<f64> = cols[j].iter().copied().filter(|x| x.is_finite()).collect();
            let m = fin.iter().sum::<f64>() / fin.len() as f64;
            writeln!(f, "{p}\t{m:.3}\t{rt:.3}\t{rz:.3}\t{n}").unwrap();
        }
    }

    // 2. типы опухолей
    let dmap: HashMap<String, String> = fs::read_to_string(dis).unwrap().lines().filter_map(|l| {
        let v: Vec<&str> = l.split('\t').collect();
        (v.len() >= 3).then(|| (v[0].to_string(), v[2].to_string()))
    }).collect();
    let mut groups: HashMap<String, Vec<f64>> = HashMap::new();
    let mut all = Vec::new();
    for (i, id) in ids.iter().enumerate() {
        if !pairavg[i].is_finite() {
            continue;
        }
        all.push(pairavg[i]);
        if let Some(d) = dmap.get(id) {
            groups.entry(d.clone()).or_default().push(pairavg[i]);
        }
    }
    let (ma, va, na) = mean_var(&all);
    let mut rows: Vec<(f64, String, f64, usize, f64)> = groups.iter().filter(|(_, v)| v.len() >= 8).map(|(d, v)| {
        let (m, var, n) = mean_var(v);
        let tstat = (m - ma) / (var / n as f64 + va / na as f64).sqrt();
        let frac = v.iter().filter(|x| **x < -0.5).count() as f64 / n as f64;
        (tstat, d.clone(), m, n, frac)
    }).collect();
    rows.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    let mut f = fs::File::create(format!("{out}_lineage.tsv")).unwrap();
    writeln!(f, "disease\tn\tmean_pair\tfrac_lt_-0.5\tt_vs_all\t(all mean {ma:.3}, n {na})").unwrap();
    for r in &rows {
        writeln!(f, "{}\t{}\t{:.3}\t{:.2}\t{:.2}", r.1, r.3, r.2, r.4, r.0).unwrap();
    }
    // самые зависимые линии
    let mut order: Vec<usize> = (0..ids.len()).filter(|&i| pairavg[i].is_finite()).collect();
    order.sort_by(|&x, &y| pairavg[x].partial_cmp(&pairavg[y]).unwrap());
    let mut f = fs::File::create(format!("{out}_toplines.tsv")).unwrap();
    writeln!(f, "ach\tdisease\tTMCO6\tZNF511").unwrap();
    for &i in order.iter().take(40) {
        writeln!(f, "{}\t{}\t{:.2}\t{:.2}", ids[i], dmap.get(&ids[i]).map(|s| s.as_str()).unwrap_or("NA"), t[i], z[i]).unwrap();
    }

    // 3. модуль: корреляты среднего профиля
    let mut cor: Vec<(f64, &str)> = (0..g).filter_map(|j| {
        let (r, _) = pearson(&pairavg, &cols[j]);
        r.is_finite().then(|| (r, header[j].as_str()))
    }).collect();
    cor.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
    let mut f = fs::File::create(format!("{out}_module.tsv")).unwrap();
    writeln!(f, "gene\tr_pairavg").unwrap();
    for (r, s) in cor.iter().take(60) {
        writeln!(f, "{s}\t{r:.3}").unwrap();
    }
    for (r, s) in cor.iter().rev().take(15) {
        writeln!(f, "{s}\t{r:.3}").unwrap();
    }
    eprintln!("lines {}, genes {g}", ids.len());
}

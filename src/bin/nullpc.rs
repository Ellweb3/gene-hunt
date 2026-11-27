//! Control for technical artifacts in co-dependency analysis.
//! nullpc <gene_effect.csv> <K> <geneA> <geneB> [extra genes...]
//! The null set consists of single-copy olfactory receptors (OR\d…) that are not expressed
//! in cancer cell lines: their shared structure equals screening noise. We extract K principal
//! components of the null set, regress them out of gene profiles, and compute partial correlations.
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;

/// Z-score over measured cell lines; missing values → 0 (mean) so all genes share the same grid.
fn standardize(v: &[f64]) -> Vec<f64> {
    let fin: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    let n = fin.len() as f64;
    let m = fin.iter().sum::<f64>() / n;
    let sd = (fin.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
    let z: Vec<f64> = v.iter().map(|x| if x.is_finite() { (x - m) / sd } else { 0.0 }).collect();
    // renormalize so that dot/(N-1) remains a correlation
    let s = (z.iter().map(|x| x * x).sum::<f64>() / (z.len() as f64 - 1.0)).sqrt();
    z.iter().map(|x| x / s).collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn r(a: &[f64], b: &[f64]) -> f64 {
    let (a, b) = (standardize(a), standardize(b));
    dot(&a, &b) / (a.len() as f64 - 1.0)
}

/// Residual after regression on orthonormal components.
fn residual(v: &[f64], pcs: &[Vec<f64>]) -> Vec<f64> {
    let mut res = standardize(v);
    for p in pcs {
        let c = dot(&res, p);
        for (x, y) in res.iter_mut().zip(p) {
            *x -= c * y;
        }
    }
    res
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let k: usize = a[2].parse().unwrap();
    let txt = fs::read_to_string(&a[1]).unwrap();
    let mut lines = txt.lines();
    let header: Vec<String> = lines.next().unwrap().split(',').skip(1).map(|s| s.split(" (").next().unwrap().to_string()).collect();
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); header.len()];
    for l in lines {
        for (j, v) in l.split(',').skip(1).enumerate() {
            cols[j].push(v.parse::<f64>().unwrap_or(f64::NAN));
        }
    }
    // only genes with enough non-missing values, so everything is computed on the same cell lines
    let full: Vec<usize> = (0..header.len()).filter(|&j| cols[j].iter().filter(|x| x.is_finite()).count() >= 1000).collect();
    let idx: HashMap<&str, usize> = full.iter().map(|&j| (header[j].as_str(), j)).collect();
    // null set: OR<digit>…, excluding multi-copy subfamilies (keep all OR genes but
    // require |mean| < 0.15 — a strong mean shift indicates a multi-copy locus)
    let null: Vec<usize> = full.iter().copied().filter(|&j| {
        let h = &header[j];
        h.starts_with("OR") && h[2..].chars().next().is_some_and(|c| c.is_ascii_digit()) && {
            let fin: Vec<f64> = cols[j].iter().copied().filter(|x| x.is_finite()).collect();
            let m = fin.iter().sum::<f64>() / fin.len() as f64;
            m.abs() < 0.15
        }
    }).collect();
    let n = cols[0].len();
    eprintln!("lines {n}, full genes {}, null OR genes {}", full.len(), null.len());
    let z: Vec<Vec<f64>> = null.iter().map(|&j| standardize(&cols[j])).collect();
    // PCs in cell-line space: power iteration on X·Xᵀ without forming the matrix explicitly
    let mut pcs: Vec<Vec<f64>> = Vec::new();
    let mut total = 0.0;
    for zz in &z {
        total += dot(zz, zz);
    }
    for c in 0..k {
        let mut v: Vec<f64> = (0..n).map(|i| ((i * 7919 + c * 104729) % 1000) as f64 / 1000.0 - 0.5).collect();
        let mut lam = 0.0;
        for _ in 0..300 {
            // w = X Xᵀ v
            let mut w = vec![0.0; n];
            for zz in &z {
                let s = dot(zz, &v);
                for (wi, zi) in w.iter_mut().zip(zz) {
                    *wi += s * zi;
                }
            }
            for p in &pcs {
                let s = dot(&w, p);
                for (wi, pi) in w.iter_mut().zip(p) {
                    *wi -= s * pi;
                }
            }
            let nn = dot(&w, &w).sqrt();
            lam = nn;
            v = w.iter().map(|x| x / nn).collect();
        }
        eprintln!("PC{} explains {:.1}% of null variance", c + 1, 100.0 * lam / total);
        pcs.push(v);
    }
    let ga = idx[a[3].as_str()];
    let gb = idx[a[4].as_str()];
    println!("raw r({},{}) = {:.3}", a[3], a[4], r(&cols[ga], &cols[gb]));
    for kk in [1, 3, 5, k] {
        let ra = residual(&cols[ga], &pcs[..kk]);
        let rb = residual(&cols[gb], &pcs[..kk]);
        println!("after {kk} null PCs: r = {:.3}", r(&ra, &rb));
    }
    let ra = residual(&cols[ga], &pcs);
    let rb = residual(&cols[gb], &pcs);
    for g in &a[5..] {
        if let Some(&j) = idx.get(g.as_str()) {
            let rg = residual(&cols[j], &pcs);
            println!("{g}: raw r with {} {:.3} / {} {:.3}; after PCs {:.3} / {:.3}", a[3], r(&cols[ga], &cols[j]), a[4], r(&cols[gb], &cols[j]), r(&ra, &rg), r(&rb, &rg));
        }
    }
    // background: distribution of r between random null-gene pairs before/after correction
    let zr: Vec<Vec<f64>> = null.par_iter().map(|&j| residual(&cols[j], &pcs)).collect();
    let mut raw = Vec::new();
    let mut adj = Vec::new();
    for i in 0..null.len() {
        for j in (i + 1)..null.len() {
            raw.push(dot(&z[i], &z[j]) / (n as f64 - 1.0));
            adj.push(r(&zr[i], &zr[j]));
        }
    }
    for (name, v) in [("raw", &mut raw), ("adj", &mut adj)] {
        v.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let q = |p: f64| v[((v.len() - 1) as f64 * p) as usize];
        println!("null OR pairs {name}: median {:.3}, 99% {:.3}, 99.9% {:.3}, max {:.3}", q(0.5), q(0.99), q(0.999), v[v.len() - 1]);
    }
    // genome-wide background after correction: rank of partner gene among all genes
    let all: Vec<(f64, &str)> = full.par_iter().filter(|&&j| j != ga).map(|&j| (r(&ra, &residual(&cols[j], &pcs)), header[j].as_str())).collect();
    let mut all = all;
    all.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
    println!("top partners of {} after correction:", a[3]);
    for (rv, g) in all.iter().take(10) {
        println!("  {g}\t{rv:.3}");
    }
}

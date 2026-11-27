//! AlphaFold model profile: pLDDT by segment, helices (from Cα geometry), hydrophobicity
//! (Kyte–Doolittle, window 19), and residue accessibility (Cα neighbor count within 10 Å).
//! plddt <model.pdb> [start-end ...]
use std::fs;

fn kd(c: u8) -> f64 {
    match c {
        b'I' => 4.5, b'V' => 4.2, b'L' => 3.8, b'F' => 2.8, b'C' => 2.5, b'M' => 1.9, b'A' => 1.8,
        b'G' => -0.4, b'T' => -0.7, b'S' => -0.8, b'W' => -0.9, b'Y' => -1.3, b'P' => -1.6, b'H' => -3.2,
        b'E' | b'Q' | b'D' | b'N' => -3.5, b'K' => -3.9, b'R' => -4.5, _ => 0.0,
    }
}

fn aa1(s: &str) -> u8 {
    let t = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"];
    let o = b"ARNDCQEGHILKMFPSTWYV";
    t.iter().position(|x| *x == s).map(|i| o[i]).unwrap_or(b'X')
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let txt = fs::read_to_string(&a[1]).unwrap();
    let mut seq = Vec::new();
    let mut pl = Vec::new();
    let mut ca: Vec<[f64; 3]> = Vec::new();
    for l in txt.lines().filter(|l| l.starts_with("ATOM") && &l[12..16] == " CA ") {
        seq.push(aa1(&l[17..20]));
        pl.push(l[60..66].trim().parse::<f64>().unwrap());
        let p = |i: usize| l[i..i + 8].trim().parse::<f64>().unwrap();
        ca.push([p(30), p(38), p(46)]);
    }
    let n = seq.len();
    let d = |i: usize, j: usize| ((0..3).map(|k| (ca[i][k] - ca[j][k]).powi(2)).sum::<f64>()).sqrt();
    // α-helix: d(i,i+3) ≈ 5.0–5.6 Å and d(i,i+4) ≈ 5.8–6.6 Å
    let helix: Vec<bool> = (0..n).map(|i| {
        let ok = |s: usize| s + 4 < n && (4.6..5.8).contains(&d(s, s + 3)) && (5.6..6.9).contains(&d(s, s + 4));
        (i.saturating_sub(3)..=i).any(ok)
    }).collect();
    let nb: Vec<usize> = (0..n).map(|i| (0..n).filter(|&j| j != i && d(i, j) < 10.0).count()).collect();
    println!("length {n}");
    println!("sequence {}", String::from_utf8_lossy(&seq));
    // 25-residue blocks
    println!("block\tmean_pLDDT\thelix%\tKDmax19\tmean_neighbors");
    for s in (0..n).step_by(25) {
        let e = (s + 25).min(n);
        let m = pl[s..e].iter().sum::<f64>() / (e - s) as f64;
        let h = helix[s..e].iter().filter(|x| **x).count() as f64 / (e - s) as f64;
        let kdmax = (s..e).filter(|&i| i + 19 <= n).map(|i| seq[i..i + 19].iter().map(|c| kd(*c)).sum::<f64>() / 19.0).fold(f64::MIN, f64::max);
        let nbm = nb[s..e].iter().sum::<usize>() as f64 / (e - s) as f64;
        println!("{}-{}\t{:.0}\t{:.0}\t{:.2}\t{:.1}", s + 1, e, m, 100.0 * h, kdmax, nbm);
    }
    for r in &a[2..] {
        let (x, y) = r.split_once('-').unwrap();
        let (x, y): (usize, usize) = (x.parse().unwrap(), y.parse().unwrap());
        let seg = &seq[x - 1..y];
        let m = pl[x - 1..y].iter().sum::<f64>() / (y - x + 1) as f64;
        let kdm = seg.iter().map(|c| kd(*c)).sum::<f64>() / seg.len() as f64;
        let nbm = nb[x - 1..y].iter().sum::<usize>() as f64 / seg.len() as f64;
        let charged = seg.iter().filter(|c| b"KRDE".contains(c)).count();
        println!("{r}: {} pLDDT {m:.0}, KD {kdm:.2}, charged {charged}, CA neighbors {nbm:.1}", String::from_utf8_lossy(seg));
    }
}

//! Does the outward-facing surface of a helix look like lipid or like solvent?
//!
//! lipid_face <model.pdb> <start-end> ...
//!
//! Rationale: burial alone does not distinguish a transmembrane helix from a buried helix of
//! a soluble solenoid — in a multi-pass membrane protein the TM helices pack against each
//! other and are just as buried. What should differ is the chemistry of the EXPOSED face:
//! a TM helix presents hydrophobic residues to the lipid, whereas a core helix of a soluble
//! protein presents polar residues to the solvent.
//!
//! Exposure is approximated by the Calpha neighbour count within 10 A; residues in the most
//! exposed tertile of the whole chain are treated as surface. For each requested segment the
//! tool reports the mean Kyte-Doolittle hydrophobicity of its surface residues, and the same
//! quantity for the rest of the chain as an internal reference.
use std::fs;

fn kd(c: char) -> f64 {
    match c {
        'I' => 4.5, 'V' => 4.2, 'L' => 3.8, 'F' => 2.8, 'C' => 2.5, 'M' => 1.9, 'A' => 1.8,
        'G' => -0.4, 'T' => -0.7, 'S' => -0.8, 'W' => -0.9, 'Y' => -1.3, 'P' => -1.6, 'H' => -3.2,
        'E' | 'Q' | 'D' | 'N' => -3.5, 'K' => -3.9, 'R' => -4.5, _ => 0.0,
    }
}

fn aa1(s: &str) -> char {
    let t = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"];
    let o = "ARNDCQEGHILKMFPSTWYV";
    t.iter().position(|x| *x == s).map(|i| o.as_bytes()[i] as char).unwrap_or('X')
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let txt = fs::read_to_string(&a[1]).unwrap();
    let mut ca: Vec<[f64; 3]> = Vec::new();
    let mut seq: Vec<char> = Vec::new();
    for l in txt.lines().filter(|l| l.starts_with("ATOM") && &l[12..16] == " CA ") {
        let p = |i: usize| l[i..i + 8].trim().parse::<f64>().unwrap();
        ca.push([p(30), p(38), p(46)]);
        seq.push(aa1(&l[17..20]));
    }
    let n = ca.len();
    let nb: Vec<usize> = (0..n).map(|i| (0..n).filter(|&j| j != i && (0..3).map(|k| (ca[i][k] - ca[j][k]).powi(2)).sum::<f64>() < 100.0).count()).collect();
    // exposure threshold: most exposed third of the chain
    let mut sorted = nb.clone();
    sorted.sort_unstable();
    let thr = sorted[n / 3];
    let exposed: Vec<bool> = nb.iter().map(|&x| x <= thr).collect();

    let segs: Vec<(usize, usize)> = a[2..].iter().map(|r| {
        let (x, y) = r.split_once('-').unwrap();
        (x.parse::<usize>().unwrap(), y.parse::<usize>().unwrap())
    }).collect();
    let in_seg = |i: usize| segs.iter().any(|&(x, y)| i + 1 >= x && i + 1 <= y);

    let stat = |idx: Vec<usize>| -> (f64, usize) {
        let e: Vec<usize> = idx.into_iter().filter(|&i| exposed[i]).collect();
        let m = if e.is_empty() { f64::NAN } else { e.iter().map(|&i| kd(seq[i])).sum::<f64>() / e.len() as f64 };
        (m, e.len())
    };

    println!("chain length {n}, exposure threshold <= {thr} neighbours");
    println!("segment\tsurface_res\tmean_KD_of_surface");
    for &(x, y) in &segs {
        let (m, c) = stat((x - 1..y).collect());
        println!("{x}-{y}\t{c}\t{m:+.2}");
    }
    let (m, c) = stat((0..n).filter(|&i| in_seg(i)).collect());
    println!("ALL_SEGMENTS\t{c}\t{m:+.2}");
    let (m2, c2) = stat((0..n).filter(|&i| !in_seg(i)).collect());
    println!("REST_OF_CHAIN\t{c2}\t{m2:+.2}");
    println!("difference (segments - rest): {:+.2}", m - m2);
}

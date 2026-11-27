//! Sequence identity to two reference sequences (SW, fraction of matching aligned columns,
//! and reference coverage). ident <seqs.fa> <refA.fa> <refB.fa> <acc_list.txt>
use gene_hunt::align::{encode, sw_align};
use std::collections::HashSet;
use std::fs;

fn fasta(p: &str) -> Vec<(String, String, Vec<u8>)> {
    fs::read_to_string(p).unwrap().split('>').skip(1).map(|r| {
        let mut l = r.lines();
        let h = l.next().unwrap().to_string();
        let acc = h.split('|').nth(1).unwrap_or("").to_string();
        (acc, h, l.collect::<String>().into_bytes())
    }).collect()
}

fn ident(r: &[u8], s: &[u8]) -> (f64, f64) {
    let al = sw_align(&encode(r), &encode(s));
    let cols = al.qa.iter().zip(&al.sa).filter(|(x, y)| **x != b'-' && **y != b'-').count();
    let same = al.qa.iter().zip(&al.sa).filter(|(x, y)| x == y && **x != b'-').count();
    (same as f64 / cols.max(1) as f64, (al.q.1 - al.q.0) as f64 / r.len() as f64)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let seqs = fasta(&a[1]);
    let ra = &fasta(&a[2])[0];
    let rb = &fasta(&a[3])[0];
    let want: HashSet<String> = fs::read_to_string(&a[4]).unwrap().lines().map(String::from).collect();
    println!("acc\tid_{}\tcov\tid_{}\tcov\tname", ra.0, rb.0);
    for (acc, h, s) in seqs.iter().filter(|x| want.contains(&x.0)) {
        let (ia, ca) = ident(&ra.2, s);
        let (ib, cb) = ident(&rb.2, s);
        let os = h.split("OS=").nth(1).and_then(|x| x.split(" OX=").next()).unwrap_or("");
        println!("{acc}\t{ia:.2}\t{ca:.2}\t{ib:.2}\t{cb:.2}\t{os}");
    }
}

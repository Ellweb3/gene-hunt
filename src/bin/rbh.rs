//! Best hit against a reference panel (e.g., all human importins-α + TMCO6).
//! rbh <seqs.fa> <refs.fa> <acc_list.txt>
//! Prints for each sequence the best reference hit, the gap to the second hit (in bits), and
//! a per-group summary from *.species (if a 4th argument is supplied).
use gene_hunt::align::{bits, encode, sw_best};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;

fn fasta(p: &str) -> Vec<(String, String, Vec<u8>)> {
    fs::read_to_string(p).unwrap().split('>').skip(1).map(|r| {
        let mut l = r.lines();
        let h = l.next().unwrap().to_string();
        let acc = h.split('|').nth(1).unwrap_or("").to_string();
        let name = h.split('|').nth(2).and_then(|x| x.split(' ').next()).unwrap_or("").to_string();
        (acc, name, l.collect::<String>().into_bytes())
    }).collect()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let seqs = fasta(&a[1]);
    let refs: Vec<(String, Vec<u8>)> = fasta(&a[2]).into_iter().map(|(_, n, s)| (n, encode(&s))).collect();
    let want: HashSet<String> = fs::read_to_string(&a[3]).unwrap().lines().map(String::from).collect();
    let group: HashMap<String, String> = a.get(4).map(|p| fs::read_to_string(p).unwrap().lines().skip(1).map(|l| {
        let v: Vec<&str> = l.split('\t').collect();
        (v[0].to_string(), v[2].to_string())
    }).collect()).unwrap_or_default();
    let res: Vec<(String, String, f64)> = seqs.par_iter().filter(|x| want.contains(&x.0)).map(|(acc, _, s)| {
        let e = encode(s);
        let mut sc: Vec<(f64, &str)> = refs.iter().map(|(n, r)| (bits(sw_best(r, &e).0), n.as_str())).collect();
        sc.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
        (acc.clone(), sc[0].1.to_string(), sc[0].0 - sc[1].0)
    }).collect();
    let mut summ: BTreeMap<(String, String), (usize, f64)> = BTreeMap::new();
    for (acc, best, margin) in &res {
        let g = group.get(acc).cloned().unwrap_or_else(|| "?".into());
        let e = summ.entry((g, best.clone())).or_insert((0, 0.0));
        e.0 += 1;
        e.1 += margin;
    }
    println!("group\tbest_human_hit\tn\tmean_margin_bits");
    for ((g, b), (n, m)) in summ {
        println!("{g}\t{b}\t{n}\t{:.1}", m / n as f64);
    }
}

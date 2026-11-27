//! Интерфейс двух цепей модели: остатки с тяжёлыми атомами ближе 4.5 Å к другой цепи.
//! interface <complex.pdb> [posA ...]  — дополнительно печатает, где лежат заданные остатки цепи A
//! (минимальное расстояние до цепи B).
use std::collections::{BTreeMap, BTreeSet};
use std::fs;

struct Atom {
    chain: char,
    res: usize,
    name: String,
    aa: String,
    xyz: [f64; 3],
    b: f64,
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let atoms: Vec<Atom> = fs::read_to_string(&a[1]).unwrap().lines().filter(|l| l.starts_with("ATOM")).map(|l| {
        let p = |i: usize| l[i..i + 8].trim().parse::<f64>().unwrap();
        Atom {
            chain: l.as_bytes()[21] as char,
            res: l[22..26].trim().parse().unwrap(),
            name: l[12..16].trim().to_string(),
            aa: l[17..20].to_string(),
            xyz: [p(30), p(38), p(46)],
            b: l[60..66].trim().parse().unwrap(),
        }
    }).filter(|x| !x.name.starts_with('H')).collect();
    let chains: BTreeSet<char> = atoms.iter().map(|x| x.chain).collect();
    let ch: Vec<char> = chains.into_iter().collect();
    let d2 = |x: &Atom, y: &Atom| (0..3).map(|k| (x.xyz[k] - y.xyz[k]).powi(2)).sum::<f64>();
    let ca = |c: char| atoms.iter().filter(move |x| x.chain == c && x.name == "CA");
    // минимальное расстояние от каждого остатка до другой цепи
    let mut mind: BTreeMap<(char, usize), (f64, String, f64)> = BTreeMap::new();
    for x in &atoms {
        let other = if x.chain == ch[0] { ch[1] } else { ch[0] };
        let m = atoms.iter().filter(|y| y.chain == other).map(|y| d2(x, y)).fold(f64::MAX, f64::min).sqrt();
        let e = mind.entry((x.chain, x.res)).or_insert((f64::MAX, x.aa.clone(), 0.0));
        if m < e.0 {
            e.0 = m;
        }
        if x.name == "CA" {
            e.2 = x.b;
        }
    }
    for &c in &ch {
        let n = ca(c).count();
        let iface: Vec<String> = mind.iter().filter(|((cc, _), v)| *cc == c && v.0 < 4.5).map(|((_, r), v)| format!("{}{r}", &v.1)).collect();
        let mean_pl = mind.iter().filter(|((cc, _), v)| *cc == c && v.0 < 4.5).map(|(_, v)| v.2).sum::<f64>() / iface.len().max(1) as f64;
        println!("chain {c}: {n} residues, interface {} residues (mean pLDDT {mean_pl:.0})", iface.len());
        println!("  {}", iface.join(" "));
    }
    for p in &a[2..] {
        let r: usize = p.parse().unwrap();
        if let Some(v) = mind.get(&(ch[0], r)) {
            println!("A {}{r}: {:.1} Å to chain {}", v.1, v.0, ch[1]);
        }
    }
}

//! Консервация остатков белка человека по ортологам (одна последовательность на вид).
//! conserve <orth.fa> <fam.tsv с lineage> <human_acc> <model.pdb> <out.tsv> [pos ...]
//! Для каждого вида берётся лучшая по SW-счёту последовательность (защита от паралогов),
//! выравнивания с покрытием < 60% длины человеческого белка отбрасываются.
//! Отдельно считается «далёкий» набор — не позвоночные.
use gene_hunt::align::{encode, sw_align};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::io::Write;

fn read_fasta(p: &str) -> Vec<(String, String, Vec<u8>)> {
    let mut out = Vec::new();
    for rec in fs::read_to_string(p).unwrap().split('>').skip(1) {
        let mut l = rec.lines();
        let h = l.next().unwrap();
        let acc = h.split('|').nth(1).unwrap_or("").to_string();
        let ox = h.split("OX=").nth(1).and_then(|x| x.split(' ').next()).unwrap_or("").to_string();
        out.push((acc, ox, l.collect::<String>().into_bytes()));
    }
    out
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let seqs = read_fasta(&a[1]);
    let lineage: HashMap<String, String> = fs::read_to_string(&a[2]).unwrap().lines().skip(1).filter_map(|l| {
        let v: Vec<&str> = l.split('\t').collect();
        (v.len() > 2).then(|| (v[0].to_string(), v[2].to_string()))
    }).collect();
    let vert: HashMap<String, bool> = lineage.iter().map(|(k, v)| (k.clone(), v.contains("Vertebrata"))).collect();
    // контрольная последовательность (паралог, например импортин-α) для проверки ортологии
    let control: Option<Vec<u8>> = std::env::var("CONTROL").ok().map(|p| encode(&read_fasta(&p)[0].2));
    let human = seqs.iter().find(|s| s.0 == a[3]).expect("human").2.clone();
    let hq = encode(&human);
    let n = human.len();
    // соседи Cα (экспонированность) из модели
    let pdb = fs::read_to_string(&a[4]).unwrap();
    let ca: Vec<[f64; 3]> = pdb.lines().filter(|l| l.starts_with("ATOM") && &l[12..16] == " CA ").map(|l| {
        let p = |i: usize| l[i..i + 8].trim().parse::<f64>().unwrap();
        [p(30), p(38), p(46)]
    }).collect();
    let plddt: Vec<f64> = pdb.lines().filter(|l| l.starts_with("ATOM") && &l[12..16] == " CA ").map(|l| l[60..66].trim().parse().unwrap()).collect();
    let nb: Vec<usize> = (0..n).map(|i| (0..n).filter(|&j| j != i && (0..3).map(|k| (ca[i][k] - ca[j][k]).powi(2)).sum::<f64>() < 100.0).count()).collect();

    // выравнивания
    let alns: Vec<(String, String, i32, Vec<Option<u8>>)> = seqs.par_iter().filter(|s| s.0 != a[3]).map(|(acc, ox, s)| {
        let al = sw_align(&hq, &encode(s));
        let mut col: Vec<Option<u8>> = vec![None; n];
        let mut qi = al.q.0;
        for (x, y) in al.qa.iter().zip(&al.sa) {
            if *x != b'-' {
                if *y != b'-' {
                    col[qi] = Some(*y);
                }
                qi += 1;
            }
        }
        (acc.clone(), ox.clone(), al.score, col)
    }).collect();
    let mut best: HashMap<&str, usize> = HashMap::new();
    for (i, x) in alns.iter().enumerate() {
        let cov = x.3.iter().filter(|c| c.is_some()).count() as f64 / n as f64;
        if cov < 0.6 {
            continue;
        }
        let e = best.entry(x.1.as_str()).or_insert(i);
        if alns[*e].2 < x.2 {
            *e = i;
        }
    }
    let chosen: Vec<usize> = best.values().copied().collect();
    let far: Vec<usize> = chosen.iter().copied().filter(|&i| !vert.get(&alns[i].0).copied().unwrap_or(true)).collect();
    eprintln!("species used: {} (non-vertebrate {})", chosen.len(), far.len());

    let stat = |set: &[usize], p: usize| -> (f64, f64) {
        let present: Vec<u8> = set.iter().filter_map(|&i| alns[i].3[p]).collect();
        let id = present.iter().filter(|&&c| c == human[p]).count() as f64 / present.len().max(1) as f64;
        (id, present.len() as f64 / set.len().max(1) as f64)
    };
    // по видам: группа, остатки в запрошенных позициях, счёт против человека и против контроля
    let groups = ["Mammalia", "Aves", "Testudines", "Crocodylia", "Lepidosauria", "Amphibia", "Coelacanthiformes", "Dipnomorpha", "Actinopterygii", "Chondrichthyes", "Cyclostomata", "Tunicata", "Cephalochordata", "Echinodermata", "Hemichordata", "Mollusca", "Annelida", "Brachiopoda", "Cnidaria", "Filasterea", "Fungi", "Viridiplantae"];
    let pos: Vec<usize> = a[6..].iter().map(|s| s.parse::<usize>().unwrap() - 1).collect();
    let mut fs_ = fs::File::create(format!("{}.species", a[5])).unwrap();
    writeln!(fs_, "acc\ttaxid\tgroup\tscore_vs_query\tscore_vs_control\tresidues").unwrap();
    for &i in &chosen {
        let lin = lineage.get(&alns[i].0).map(|s| s.as_str()).unwrap_or("");
        let g = groups.iter().find(|g| lin.contains(&format!("{} (", g))).unwrap_or(&"other");
        let res: String = pos.iter().map(|&p| alns[i].3[p].map(|c| c as char).unwrap_or('-')).collect();
        let sc_ctrl = control.as_ref().map(|c| {
            let s = seqs.iter().find(|s| s.0 == alns[i].0).unwrap();
            sw_align(c, &encode(&s.2)).score
        }).unwrap_or(0);
        writeln!(fs_, "{}\t{}\t{g}\t{}\t{sc_ctrl}\t{res}", alns[i].0, alns[i].1, alns[i].2).unwrap();
    }
    let mut f = fs::File::create(&a[5]).unwrap();
    writeln!(f, "pos\taa\tid_all\tid_nonvert\tcov_nonvert\tneighbors\tpLDDT").unwrap();
    let mut rows = Vec::new();
    for p in 0..n {
        let (ia, _) = stat(&chosen, p);
        let (iv, cv) = stat(&far, p);
        writeln!(f, "{}\t{}\t{ia:.2}\t{iv:.2}\t{cv:.2}\t{}\t{:.0}", p + 1, human[p] as char, nb[p], plddt[p]).unwrap();
        rows.push((p, ia, iv, cv));
    }
    for s in &a[6..] {
        let p: usize = s.parse::<usize>().unwrap() - 1;
        let r = rows[p];
        // распределение остатков в далёком наборе
        let mut cnt: HashMap<u8, usize> = HashMap::new();
        for &i in &far {
            if let Some(c) = alns[i].3[p] {
                *cnt.entry(c).or_default() += 1;
            }
        }
        let mut c: Vec<_> = cnt.into_iter().collect();
        c.sort_by(|x, y| y.1.cmp(&x.1));
        let top: Vec<String> = c.iter().take(4).map(|(k, v)| format!("{}:{v}", *k as char)).collect();
        println!("{}{}\tall {:.2}\tnonvert {:.2} (cov {:.2})\tnb {}\t{}", human[p] as char, p + 1, r.1, r.2, r.3, nb[p], top.join(" "));
    }
    // самые консервативные экспонированные остатки (кандидаты интерфейса)
    let mut exp: Vec<_> = rows.iter().filter(|r| nb[r.0] <= 14 && r.3 >= 0.7 && plddt[r.0] >= 70.0).collect();
    exp.sort_by(|x, y| y.2.partial_cmp(&x.2).unwrap());
    println!("most conserved exposed residues (non-vertebrate identity):");
    for r in exp.iter().take(30) {
        print!("{}{}:{:.2} ", human[r.0] as char, r.0 + 1, r.2);
    }
    println!();
}

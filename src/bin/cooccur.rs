//! Филогенетическое совместное присутствие двух семейств в референсных протеомах.
//! cooccur <refprot.tsv> <famA.tsv> <famB.tsv> <min_busco_C>
//! refprot: upid, taxid, lineage, protein_count, BUSCO; fam: accession, taxid, …
use std::collections::{BTreeMap, HashSet};
use std::fs;

const GROUPS: &[&str] = &[
    "Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria", "Amphibia", "Coelacanthiformes", "Dipnomorpha",
    "Teleostei", "Holostei", "Chondrostei", "Polypteriformes", "Chondrichthyes", "Cyclostomata", "Tunicata", "Cephalochordata",
    "Echinodermata", "Hemichordata", "Insecta", "Crustacea", "Chelicerata", "Myriapoda", "Mollusca", "Annelida",
    "Brachiopoda", "Platyhelminthes", "Rotifera", "Nematoda", "Tardigrada", "Cnidaria", "Placozoa", "Porifera", "Ctenophora",
    "Choanoflagellata", "Ichthyosporea", "Filasterea", "Fungi", "Amoebozoa", "Streptophyta", "Chlorophyta", "Rhodophyta",
    "Stramenopiles", "Alveolata", "Rhizaria", "Discoba", "Metamonada", "Haptista", "Cryptophyceae",
];

fn busco_c(s: &str) -> f64 {
    s.strip_prefix("C:").and_then(|r| r.split('%').next()).and_then(|x| x.parse().ok()).unwrap_or(0.0)
}

fn taxa(path: &str) -> HashSet<String> {
    fs::read_to_string(path).unwrap().lines().skip(1).filter_map(|l| l.split('\t').nth(1).map(String::from)).collect()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let min_c: f64 = a[4].parse().unwrap();
    let ta = taxa(&a[2]);
    let tb = taxa(&a[3]);
    // группа → [neither, A only, B only, both]
    let mut g: BTreeMap<usize, [usize; 4]> = BTreeMap::new();
    let mut tot = [0usize; 4];
    let mut seen = HashSet::new();
    let mut only: BTreeMap<usize, Vec<String>> = BTreeMap::new();
    for l in fs::read_to_string(&a[1]).unwrap().lines().skip(1) {
        let v: Vec<&str> = l.split('\t').collect();
        if v.len() < 5 || busco_c(v[4]) < min_c || !seen.insert(v[1].to_string()) {
            continue;
        }
        let lin: Vec<&str> = v[2].split(", ").collect();
        let Some(gi) = GROUPS.iter().position(|grp| lin.contains(grp)) else { continue };
        let (ha, hb) = (ta.contains(v[1]), tb.contains(v[1]));
        let k = (ha as usize) + 2 * (hb as usize);
        g.entry(gi).or_insert([0; 4])[k] += 1;
        tot[k] += 1;
        if k == 1 || k == 2 {
            let name = lin.last().unwrap_or(&"");
            only.entry(gi).or_default().push(format!("{}:{}", if k == 1 { "A" } else { "B" }, name));
        }
    }
    println!("group\tn\tneither\tA_only\tB_only\tboth");
    for (gi, c) in &g {
        let n: usize = c.iter().sum();
        println!("{}\t{n}\t{}\t{}\t{}\t{}", GROUPS[*gi], c[0], c[1], c[2], c[3]);
    }
    let n: f64 = tot.iter().sum::<usize>() as f64;
    let (n00, n10, n01, n11) = (tot[0] as f64, tot[1] as f64, tot[2] as f64, tot[3] as f64);
    let phi = (n11 * n00 - n10 * n01) / ((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00)).sqrt();
    println!("TOTAL n={n} neither={n00} A_only={n10} B_only={n01} both={n11} phi={phi:.3}");
    // где наборы расходятся — с названием вида
    for (gi, v) in &only {
        if v.len() <= 12 {
            println!("discordant {}: {}", GROUPS[*gi], v.join("; "));
        } else {
            println!("discordant {}: {} species", GROUPS[*gi], v.len());
        }
    }
}

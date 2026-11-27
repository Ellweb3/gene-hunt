//! Neighbor-joining (NJ) tree from pairwise alignments with bootstrap support.
//! njtree <seqs.fa> <out_prefix> [reps]
//! Distance: fraction of mismatches in aligned SW columns, Kimura correction for proteins:
//! d = −ln(1 − p − 0.2p²). Bootstrap: for each pair the alignment columns are resampled
//! with replacement, NJ is rebuilt, and bipartition support is counted.
//! Sequence labels must start with a group tag: TAG|acc|name.
use gene_hunt::align::{encode, sw_align};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::io::Write;

type Split = Vec<u64>;

fn kimura(p: f64) -> f64 {
    let x = 1.0 - p - 0.2 * p * p;
    if x <= 0.05 { 3.0 } else { -x.ln() }
}

struct Node {
    ch: Vec<(usize, f64)>,
}

/// NJ; returns nodes (leaves 0..n) and the root index (last three-way join).
fn nj(d0: &[Vec<f64>]) -> (Vec<Node>, usize) {
    let n = d0.len();
    let mut nodes: Vec<Node> = (0..n).map(|_| Node { ch: vec![] }).collect();
    let mut active: Vec<usize> = (0..n).collect();
    let mut d: HashMap<(usize, usize), f64> = HashMap::new();
    for i in 0..n {
        for j in 0..n {
            d.insert((i, j), d0[i][j]);
        }
    }
    let get = |d: &HashMap<(usize, usize), f64>, a: usize, b: usize| if a == b { 0.0 } else { d[&(a, b)] };
    while active.len() > 3 {
        let m = active.len() as f64;
        let r: Vec<f64> = active.iter().map(|&a| active.iter().map(|&b| get(&d, a, b)).sum()).collect();
        let (mut bi, mut bj, mut bq) = (0, 1, f64::MAX);
        for i in 0..active.len() {
            for j in (i + 1)..active.len() {
                let q = (m - 2.0) * get(&d, active[i], active[j]) - r[i] - r[j];
                if q < bq {
                    bq = q;
                    bi = i;
                    bj = j;
                }
            }
        }
        let (a, b) = (active[bi], active[bj]);
        let dab = get(&d, a, b);
        let la = (0.5 * dab + (r[bi] - r[bj]) / (2.0 * (m - 2.0))).max(0.0);
        let lb = (dab - la).max(0.0);
        let u = nodes.len();
        nodes.push(Node { ch: vec![(a, la), (b, lb)] });
        for &k in &active {
            if k != a && k != b {
                let v = 0.5 * (get(&d, a, k) + get(&d, b, k) - dab);
                d.insert((u, k), v);
                d.insert((k, u), v);
            }
        }
        active.retain(|&x| x != a && x != b);
        active.push(u);
    }
    let (a, b, c) = (active[0], active[1], active[2]);
    let (dab, dac, dbc) = (get(&d, a, b), get(&d, a, c), get(&d, b, c));
    let root = nodes.len();
    nodes.push(Node { ch: vec![(a, ((dab + dac - dbc) / 2.0).max(0.0)), (b, ((dab + dbc - dac) / 2.0).max(0.0)), (c, ((dac + dbc - dab) / 2.0).max(0.0))] });
    (nodes, root)
}

fn leaves(nodes: &[Node], u: usize, n: usize, out: &mut Split) {
    if u < n {
        out[u / 64] |= 1 << (u % 64);
    }
    for &(c, _) in &nodes[u].ch {
        leaves(nodes, c, n, out);
    }
}

/// Canonical bipartition representation: the side that does not contain leaf 0.
fn canon(s: &Split, n: usize) -> Split {
    if s[0] & 1 == 0 {
        return s.clone();
    }
    let mut c: Split = s.iter().map(|x| !x).collect();
    let extra = c.len() * 64 - n;
    if extra > 0 {
        let last = c.len() - 1;
        c[last] &= u64::MAX >> extra;
    }
    c
}

fn splits(nodes: &[Node], root: usize, n: usize) -> HashMap<usize, Split> {
    let w = n.div_ceil(64);
    let mut m = HashMap::new();
    for u in n..nodes.len() {
        if u == root {
            continue;
        }
        let mut s = vec![0u64; w];
        leaves(nodes, u, n, &mut s);
        m.insert(u, canon(&s, n));
    }
    m
}

fn newick(nodes: &[Node], u: usize, n: usize, names: &[String], sup: &HashMap<usize, f64>) -> String {
    if u < n {
        return names[u].clone();
    }
    let inner: Vec<String> = nodes[u].ch.iter().map(|&(c, l)| format!("{}:{l:.4}", newick(nodes, c, n, names, sup))).collect();
    let s = sup.get(&u).map(|x| format!("{:.0}", x)).unwrap_or_default();
    format!("({}){s}", inner.join(","))
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let reps: usize = a.get(3).map(|x| x.parse().unwrap()).unwrap_or(100);
    let mut names = Vec::new();
    let mut seqs = Vec::new();
    for rec in fs::read_to_string(&a[1]).unwrap().split('>').skip(1) {
        let mut l = rec.lines();
        names.push(l.next().unwrap().to_string());
        seqs.push(encode(l.collect::<String>().as_bytes()));
    }
    let n = names.len();
    // pairwise match vectors
    let pairs: Vec<(usize, usize)> = (0..n).flat_map(|i| ((i + 1)..n).map(move |j| (i, j))).collect();
    let matches: Vec<Vec<bool>> = pairs.par_iter().map(|&(i, j)| {
        let al = sw_align(&seqs[i], &seqs[j]);
        al.qa.iter().zip(&al.sa).filter(|(x, y)| **x != b'-' && **y != b'-').map(|(x, y)| x == y).collect()
    }).collect();
    let dist = |m: &[bool]| if m.len() < 80 { 3.0 } else { kimura(m.iter().filter(|x| !**x).count() as f64 / m.len() as f64) };
    let mut d = vec![vec![0.0; n]; n];
    for (k, &(i, j)) in pairs.iter().enumerate() {
        let v = dist(&matches[k]);
        d[i][j] = v;
        d[j][i] = v;
    }
    let (tree, root) = nj(&d);
    let ref_splits = splits(&tree, root, n);
    // bootstrap
    let boots: Vec<Vec<Split>> = (0..reps).into_par_iter().map(|r| {
        let mut rng = Rng(0x9E3779B97F4A7C15 ^ (r as u64 + 1).wrapping_mul(0xBF58476D1CE4E5B9));
        let mut db = vec![vec![0.0; n]; n];
        for (k, &(i, j)) in pairs.iter().enumerate() {
            let m = &matches[k];
            let res: Vec<bool> = (0..m.len()).map(|_| m[(rng.next() % m.len().max(1) as u64) as usize]).collect();
            let v = dist(&res);
            db[i][j] = v;
            db[j][i] = v;
        }
        let (t, rt) = nj(&db);
        splits(&t, rt, n).into_values().collect()
    }).collect();
    let mut count: HashMap<Split, usize> = HashMap::new();
    for b in &boots {
        for s in b {
            *count.entry(s.clone()).or_default() += 1;
        }
    }
    let sup: HashMap<usize, f64> = ref_splits.iter().map(|(&u, s)| (u, 100.0 * *count.get(s).unwrap_or(&0) as f64 / reps as f64)).collect();
    fs::write(format!("{}.nwk", a[2]), newick(&tree, root, n, &names, &sup) + ";\n").unwrap();

    // hypothesis tests: exact bipartitions by tag
    let tag = |i: usize| names[i].split('|').next().unwrap().to_string();
    let w = n.div_ceil(64);
    let set_of = |pred: &dyn Fn(usize) -> bool| {
        let mut s = vec![0u64; w];
        for i in (0..n).filter(|&i| pred(i)) {
            s[i / 64] |= 1 << (i % 64);
        }
        canon(&s, n)
    };
    let mut f = fs::File::create(format!("{}_tests.txt", a[2])).unwrap();
    let tests: Vec<(&str, Split)> = vec![
        ("vertebrate TMCO6 are monophyletic", set_of(&|i| tag(i) == "TMCO6v")),
        ("non-vertebrate PF30746 are monophyletic", set_of(&|i| tag(i) == "PF30746")),
        ("TMCO6v + PF30746 — single clade vs. all KPNA", set_of(&|i| tag(i) != "KPNA")),
    ];
    for (name, s) in &tests {
        let in_tree = ref_splits.values().any(|x| x == s);
        let b = 100.0 * *count.get(s).unwrap_or(&0) as f64 / reps as f64;
        writeln!(f, "{name}: in NJ tree {}; bootstrap {b:.0}%", if in_tree { "YES" } else { "no" }).unwrap();
    }
    // smallest clade (bipartition side) containing all TMCO6v: its composition and support
    let tm: Vec<usize> = (0..n).filter(|&i| tag(i) == "TMCO6v").collect();
    let contains_all = |s: &Split, inv: bool| tm.iter().all(|&i| ((s[i / 64] >> (i % 64)) & 1 == 1) != inv);
    let mut best: Option<(usize, Split, bool, usize)> = None;
    for (&u, s) in &ref_splits {
        for inv in [false, true] {
            if contains_all(s, inv) {
                let size = (0..n).filter(|&i| ((s[i / 64] >> (i % 64)) & 1 == 1) != inv).count();
                if size > tm.len() && best.as_ref().is_none_or(|b| size < b.3) {
                    best = Some((u, s.clone(), inv, size));
                }
            }
        }
    }
    if let Some((u, s, inv, size)) = best {
        let mem: Vec<&String> = (0..n).filter(|&i| ((s[i / 64] >> (i % 64)) & 1 == 1) != inv && tag(i) != "TMCO6v").map(|i| &names[i]).collect();
        let mut by: HashMap<String, usize> = HashMap::new();
        for m in &mem {
            *by.entry(m.split('|').next().unwrap().to_string()).or_default() += 1;
        }
        writeln!(f, "smallest clade with vertebrate TMCO6: {size} leaves, support {:.0}%, excluding TMCO6v: {by:?}", sup.get(&u).unwrap_or(&0.0)).unwrap();
        for m in mem.iter().take(40) {
            writeln!(f, "   {m}").unwrap();
        }
    }
    // nearest neighbors by distance for each group
    for t in ["PF30746", "TMCO6v"] {
        let mut near: HashMap<String, usize> = HashMap::new();
        for i in (0..n).filter(|&i| tag(i) == t) {
            let j = (0..n).filter(|&j| tag(j) != t).min_by(|&x, &y| d[i][x].partial_cmp(&d[i][y]).unwrap()).unwrap();
            *near.entry(tag(j)).or_default() += 1;
        }
        writeln!(f, "nearest out-group neighbor for {t}: {near:?}").unwrap();
    }
    // nearest neighbor overall (another species): same group or importin?
    for t in ["PF30746", "TMCO6v"] {
        let mut near: HashMap<String, usize> = HashMap::new();
        for i in (0..n).filter(|&i| tag(i) == t) {
            let j = (0..n).filter(|&j| j != i).min_by(|&x, &y| d[i][x].partial_cmp(&d[i][y]).unwrap()).unwrap();
            *near.entry(tag(j)).or_default() += 1;
        }
        writeln!(f, "nearest overall neighbor for {t}: {near:?}").unwrap();
    }
    // mean distances between groups
    let plant = ["Arabidopsis", "Oryza", "Solanum", "Viridiplantae"];
    let fungi = ["Saccharomyces", "Schizosaccharomyces", "Penicillium", "Aspergillus", "Fungi"];
    let grp = |i: usize| {
        let nm = &names[i];
        let t = tag(i);
        let k = if plant.iter().any(|p| nm.contains(p)) { "plants" } else if fungi.iter().any(|p| nm.contains(p)) { "fungi" } else if nm.contains("Dictyostelium") { "amoebae" } else { "animals" };
        if t == "TMCO6v" { "vertebrate TMCO6".to_string() } else { format!("{t} {k}") }
    };
    let gs: Vec<String> = { let mut v: Vec<String> = (0..n).map(grp).collect(); v.sort(); v.dedup(); v };
    write!(f, "\nmean Kimura distance\t").unwrap();
    for g in &gs { write!(f, "{g}\t").unwrap(); }
    writeln!(f).unwrap();
    for g1 in &gs {
        write!(f, "{g1}\t").unwrap();
        for g2 in &gs {
            let (mut s, mut c) = (0.0, 0);
            for i in (0..n).filter(|&i| &grp(i) == g1) { for j in (0..n).filter(|&j| &grp(j) == g2 && j != i) { s += d[i][j]; c += 1; } }
            write!(f, "{:.2}\t", if c > 0 { s / c as f64 } else { f64::NAN }).unwrap();
        }
        writeln!(f).unwrap();
    }
    eprintln!("n={n}, reps={reps}");
}

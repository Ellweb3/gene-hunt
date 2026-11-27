//! Связан ли мотив ZNF511 (участок интерфейса) с наличием PF30746/TMCO6 у того же вида?
//! motif_presence <cons_ZNF511.tsv.species> <ref_TMCO6.tsv> [max_gaps]
//! Мотив «есть», если в выбранных позициях пропусков ≤ max_gaps. Точный тест Фишера по группам.
use std::collections::{BTreeMap, HashSet};
use std::fs;

fn ln_fact(n: usize) -> f64 {
    (1..=n).map(|x| (x as f64).ln()).sum()
}

/// двусторонний точный тест Фишера для [[a,b],[c,d]]
fn fisher(a: usize, b: usize, c: usize, d: usize) -> f64 {
    let (r1, r2, c1, n) = (a + b, c + d, a + c, a + b + c + d);
    let p = |x: usize| (ln_fact(r1) + ln_fact(r2) + ln_fact(c1) + ln_fact(n - c1) - ln_fact(n) - ln_fact(x) - ln_fact(r1 - x) - ln_fact(c1 - x) - ln_fact(r2 + x - c1)).exp();
    let lo = c1.saturating_sub(r2);
    let hi = r1.min(c1);
    let p0 = p(a);
    (lo..=hi).map(p).filter(|&v| v <= p0 * (1.0 + 1e-9)).sum::<f64>().min(1.0)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let max_gaps: usize = a.get(3).map(|x| x.parse().unwrap()).unwrap_or(2);
    let has_t: HashSet<String> = fs::read_to_string(&a[2]).unwrap().lines().skip(1).filter_map(|l| l.split('\t').nth(1).map(String::from)).collect();
    // группа → [мотив нет & T нет, мотив нет & T есть, мотив есть & T нет, мотив есть & T есть]
    let mut g: BTreeMap<String, [usize; 4]> = BTreeMap::new();
    for l in fs::read_to_string(&a[1]).unwrap().lines().skip(1) {
        let v: Vec<&str> = l.split('\t').collect();
        let motif = v[5].chars().filter(|c| *c == '-').count() <= max_gaps;
        let t = has_t.contains(v[1]);
        let k = 2 * motif as usize + t as usize;
        g.entry(v[2].to_string()).or_insert([0; 4])[k] += 1;
        g.entry("ВСЕ".into()).or_insert([0; 4])[k] += 1;
        if v[2] != "Mammalia" && v[2] != "Aves" && v[2] != "Actinopterygii" && v[2] != "Testudines" && v[2] != "Lepidosauria" && v[2] != "Amphibia" && v[2] != "Crocodylia" && v[2] != "Chondrichthyes" && v[2] != "Cyclostomata" && v[2] != "Coelacanthiformes" {
            g.entry("не позвоночные".into()).or_insert([0; 4])[k] += 1;
        }
    }
    println!("группа\tn\tмотив|T+\tмотив|T−\tp Фишера");
    for (k, c) in &g {
        let n: usize = c.iter().sum();
        let (t_pos, t_neg) = (c[1] + c[3], c[0] + c[2]);
        let f = |x: usize, y: usize| if y == 0 { "—".to_string() } else { format!("{x}/{y} ({:.0}%)", 100.0 * x as f64 / y as f64) };
        println!("{k}\t{n}\t{}\t{}\t{:.1e}", f(c[3], t_pos), f(c[2], t_neg), fisher(c[3], c[2], c[1], c[0]));
    }
}

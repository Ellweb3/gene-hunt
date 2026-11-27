//! Map positions from a target protein onto a query via a Foldseek structural alignment.
//! map_sites <aln.json> <pos[:label]> ...   (positions in target numbering, 1-based)
use std::collections::HashMap;
use std::fs;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let v: serde_json::Value = serde_json::from_str(&fs::read_to_string(&a[1]).unwrap()).unwrap();
    let qa = v["qAln"].as_str().unwrap().as_bytes();
    let ta = v["dbAln"].as_str().unwrap().as_bytes();
    let (mut qi, mut ti) = (v["qStartPos"].as_u64().unwrap() as usize, v["dbStartPos"].as_u64().unwrap() as usize);
    let mut map: HashMap<usize, (u8, usize, u8)> = HashMap::new();
    for (&q, &t) in qa.iter().zip(ta) {
        if t != b'-' && q != b'-' {
            map.insert(ti, (t, qi, q));
        }
        if q != b'-' {
            qi += 1;
        }
        if t != b'-' {
            ti += 1;
        }
    }
    println!("target {} ({} id {})", v["target"], v["eval"], v["seqId"]);
    let (mut kept, mut total) = (0, 0);
    for p in &a[2..] {
        let (pos, label) = p.split_once(':').unwrap_or((p, ""));
        let pos: usize = pos.parse().unwrap();
        total += 1;
        match map.get(&pos) {
            Some((t, q, qq)) => {
                let same = t == qq;
                if same {
                    kept += 1;
                }
                println!("  {}{pos}\t-> {}{q}\t{}\t{label}", *t as char, *qq as char, if same { "=" } else { "≠" });
            }
            None => println!("  {pos}\t-> gap\t\t{label}"),
        }
    }
    println!("conserved {kept}/{total}");
}

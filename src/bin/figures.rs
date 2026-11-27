//! Generate the manuscript figures as SVG.
//! figures <out_dir>
//! Reads the project's own result files, so the figures regenerate from the data rather than
//! from hard-coded numbers:
//!   fig1  data/tmco6/AF-TMCO6.pdb                  burial and pLDDT along the sequence
//!   fig2  af3/results/*/..._summary_confidences_0.json   ipTM across the six AF3 jobs
//!   fig3  data/tmco6/cons_TMCO6.tsv                conservation of the NLS-groove positions
use std::collections::BTreeMap;
use std::fs;

const W: f64 = 900.0;
const PAD_L: f64 = 70.0;
const PAD_R: f64 = 20.0;

struct Svg {
    body: String,
    h: f64,
}

impl Svg {
    fn new(h: f64) -> Self {
        Svg { body: String::new(), h }
    }
    fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str, extra: &str) {
        self.body += &format!("<rect x='{x:.1}' y='{y:.1}' width='{w:.1}' height='{h:.1}' fill='{fill}' {extra}/>\n");
    }
    fn line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, stroke: &str, w: f64) {
        self.body += &format!("<line x1='{x1:.1}' y1='{y1:.1}' x2='{x2:.1}' y2='{y2:.1}' stroke='{stroke}' stroke-width='{w}'/>\n");
    }
    fn path(&mut self, d: &str, stroke: &str, w: f64) {
        self.body += &format!("<path d='{d}' fill='none' stroke='{stroke}' stroke-width='{w}' stroke-linejoin='round'/>\n");
    }
    fn text(&mut self, x: f64, y: f64, s: &str, size: f64, anchor: &str, extra: &str) {
        let esc = s.replace('&', "&amp;").replace('<', "&lt;");
        self.body += &format!("<text x='{x:.1}' y='{y:.1}' font-family='Helvetica,Arial,sans-serif' font-size='{size}' text-anchor='{anchor}' {extra}>{esc}</text>\n");
    }
    fn save(&self, path: &str) {
        let s = format!("<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{:.0}' viewBox='0 0 {W} {:.0}'>\n<rect width='100%' height='100%' fill='white'/>\n{}</svg>\n", self.h, self.h, self.body);
        fs::write(path, s).unwrap();
    }
}

fn aa1(s: &str) -> char {
    let t = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"];
    let o = "ARNDCQEGHILKMFPSTWYV";
    t.iter().position(|x| *x == s).map(|i| o.as_bytes()[i] as char).unwrap_or('X')
}

/// Figure 1: the two annotated transmembrane segments sit at the burial maximum.
fn fig1(pdb_path: &str, out: &str) {
    let txt = fs::read_to_string(pdb_path).unwrap();
    let mut ca: Vec<[f64; 3]> = Vec::new();
    let mut pl: Vec<f64> = Vec::new();
    let mut seq: Vec<char> = Vec::new();
    for l in txt.lines().filter(|l| l.starts_with("ATOM") && &l[12..16] == " CA ") {
        let p = |i: usize| l[i..i + 8].trim().parse::<f64>().unwrap();
        ca.push([p(30), p(38), p(46)]);
        pl.push(l[60..66].trim().parse().unwrap());
        seq.push(aa1(&l[17..20]));
    }
    let n = ca.len();
    let nb: Vec<f64> = (0..n).map(|i| (0..n).filter(|&j| j != i && (0..3).map(|k| (ca[i][k] - ca[j][k]).powi(2)).sum::<f64>() < 100.0).count() as f64).collect();
    // smooth over a 9-residue window for legibility
    let smooth = |v: &[f64]| -> Vec<f64> {
        (0..n).map(|i| {
            let a = i.saturating_sub(4);
            let b = (i + 5).min(n);
            v[a..b].iter().sum::<f64>() / (b - a) as f64
        }).collect()
    };
    let nbs = smooth(&nb);
    let charged: Vec<f64> = (0..n).map(|i| if "KRDE".contains(seq[i]) { 1.0 } else { 0.0 }).collect();
    let chs = smooth(&charged);

    let (h_top, h_each, gap) = (46.0, 150.0, 54.0);
    let mut s = Svg::new(h_top + 2.0 * h_each + gap + 70.0);
    let plot_w = W - PAD_L - PAD_R;
    let x_of = |i: usize| PAD_L + plot_w * (i as f64) / (n as f64 - 1.0);

    s.text(PAD_L, 24.0, "Figure 1. Burial does not distinguish the annotated TMCO6 segments from real transmembrane helices", 15.0, "start", "font-weight='bold'");
    s.text(PAD_L, 40.0, "AlphaFold model AF-Q96DC7-F1; shaded bands mark UniProt TRANSMEM 338-358 and 386-406 (ECO:0000255)", 11.5, "start", "fill='#555'");

    let bands = [(338usize, 358usize), (386, 406)];
    let panels: Vec<(&str, &Vec<f64>, f64, &str)> = vec![
        ("A  Burial: Calpha neighbours within 10 A (9-residue mean)", &nbs, 26.0, "neighbours"),
        ("B  Fraction of charged residues, K/R/D/E (9-residue mean)", &chs, 0.5, "charged fraction"),
    ];
    for (pi, (title, data, ymax, ylab)) in panels.into_iter().enumerate() {
        let y0 = h_top + pi as f64 * (h_each + gap);
        let y1 = y0 + h_each;
        s.text(PAD_L, y0 - 8.0, title, 12.5, "start", "font-weight='bold'");
        // shaded TM bands
        for (a, b) in bands {
            s.rect(x_of(a - 1), y0, x_of(b - 1) - x_of(a - 1), h_each, "#ffd9b3", "");
        }
        // axes
        s.line(PAD_L, y1, W - PAD_R, y1, "#333", 1.0);
        s.line(PAD_L, y0, PAD_L, y1, "#333", 1.0);
        let ticks = 5;
        for t in 0..=ticks {
            let v = ymax * t as f64 / ticks as f64;
            let y = y1 - h_each * (v / ymax);
            s.line(PAD_L - 4.0, y, PAD_L, y, "#333", 1.0);
            let lbl = if ymax > 10.0 { format!("{v:.0}") } else { format!("{v:.1}") };
            s.text(PAD_L - 8.0, y + 4.0, &lbl, 10.5, "end", "fill='#333'");
        }
        s.text(18.0, (y0 + y1) / 2.0, ylab, 11.0, "middle", &format!("transform='rotate(-90 18 {:.1})' fill='#333'", (y0 + y1) / 2.0));
        // curve
        let mut d = String::new();
        for i in 0..n {
            let y = y1 - h_each * (data[i] / ymax).min(1.0);
            d += &format!("{}{:.1},{:.1} ", if i == 0 { "M" } else { "L" }, x_of(i), y);
        }
        s.path(&d, "#1a4b8c", 1.6);
        // x labels
        for p in [1usize, 100, 200, 300, 400, 493] {
            let x = x_of(p - 1);
            s.line(x, y1, x, y1 + 4.0, "#333", 1.0);
            s.text(x, y1 + 17.0, &p.to_string(), 10.5, "middle", "fill='#333'");
        }
        s.text((PAD_L + W - PAD_R) / 2.0, y1 + 34.0, "residue", 11.5, "middle", "fill='#333'");
    }
    let yb = h_top + 2.0 * (h_each + gap) - gap + 42.0;
    s.text(PAD_L, yb, "Burial is uninformative here: the transmembrane helices of AQP1, GLUT1 and ABCG2 are equally or more", 11.5, "start", "fill='#333'");
    s.text(PAD_L, yb + 16.0, "buried in their own models. The reassignment of TMCO6 rests on the fold match, not on this measure.", 11.5, "start", "fill='#333'");
    s.save(out);
}

/// Figure 2: ipTM of the pair across the six AlphaFold 3 jobs.
fn fig2(results_dir: &str, out: &str) {
    // collect job -> (iptm of the whole complex, pair iptm) from the summary files
    let mut jobs: BTreeMap<String, (f64, f64)> = BTreeMap::new();
    let mut stack = vec![results_dir.to_string()];
    while let Some(d) = stack.pop() {
        let Ok(rd) = fs::read_dir(&d) else { continue };
        for e in rd.flatten() {
            let p = e.path();
            if p.is_dir() {
                stack.push(p.to_string_lossy().into());
                continue;
            }
            let name = p.file_name().unwrap().to_string_lossy().to_string();
            if let Some(key) = name.strip_suffix("_summary_confidences_0.json") {
                let v: serde_json::Value = serde_json::from_str(&fs::read_to_string(&p).unwrap()).unwrap();
                let global = v["iptm"].as_f64().unwrap_or(f64::NAN);
                // chains A and B are always TMCO6 and ZNF511 in these jobs
                let pair = v["chain_pair_iptm"][0][1].as_f64().unwrap_or(f64::NAN);
                jobs.insert(key.trim_start_matches("fold_").to_string(), (global, pair));
            }
        }
    }
    let labels: Vec<(&str, &str)> = vec![
        ("4_tmco6_znf511_wt", "pair alone\n(wild type)"),
        ("5_tmco6_e134r_y171w_x_znf511", "TMCO6\nE134R/Y171W"),
        ("6_tmco6_x_znf511_r42a_f50a_r57a", "ZNF511\nR42A/F50A/R57A"),
        ("2_tmco6_znf511_lsm2_8", "with LSm2-8\n(nuclear)"),
        ("1_tmco6_znf511_lsm2_8_u6", "with LSm2-8\n+ U6 RNA"),
        ("3_control_tmco6_znf511_lsm1_7", "with LSm1-7\n(control)"),
    ];
    let h_plot = 210.0;
    let mut s = Svg::new(h_plot + 150.0);
    s.text(PAD_L, 24.0, "Figure 2. The TMCO6-ZNF511 interface is reproducible; the predicted mutants are not", 15.0, "start", "font-weight='bold'");
    s.text(PAD_L, 40.0, "AlphaFold 3, top-ranked model of each job. Bars: ipTM of the TMCO6-ZNF511 chain pair.", 11.5, "start", "fill='#555'");
    let y0 = 60.0;
    let y1 = y0 + h_plot;
    let plot_w = W - PAD_L - PAD_R;
    let bw = plot_w / labels.len() as f64;
    s.line(PAD_L, y1, W - PAD_R, y1, "#333", 1.0);
    s.line(PAD_L, y0, PAD_L, y1, "#333", 1.0);
    for t in 0..=5 {
        let v = t as f64 * 0.2;
        let y = y1 - h_plot * v;
        s.line(PAD_L - 4.0, y, PAD_L, y, "#333", 1.0);
        s.text(PAD_L - 8.0, y + 4.0, &format!("{v:.1}"), 10.5, "end", "fill='#333'");
        if t > 0 {
            s.line(PAD_L, y, W - PAD_R, y, "#eee", 1.0);
        }
    }
    s.text(20.0, (y0 + y1) / 2.0, "ipTM of the pair", 12.0, "middle", &format!("transform='rotate(-90 20 {:.1})'", (y0 + y1) / 2.0));
    for (i, (key, lab)) in labels.iter().enumerate() {
        let Some(&(_, pair)) = jobs.get(*key) else { continue };
        let x = PAD_L + bw * i as f64 + bw * 0.22;
        let w = bw * 0.56;
        let hh = h_plot * pair;
        let fill = if i == 0 { "#1a4b8c" } else if i <= 2 { "#7aa3d0" } else { "#9bbf9b" };
        s.rect(x, y1 - hh, w, hh, fill, "");
        s.text(x + w / 2.0, y1 - hh - 7.0, &format!("{pair:.2}"), 12.0, "middle", "font-weight='bold' fill='#222'");
        for (li, line) in lab.split('\n').enumerate() {
            s.text(x + w / 2.0, y1 + 18.0 + li as f64 * 13.0, line, 10.8, "middle", "fill='#333'");
        }
    }
    let yb = y1 + 62.0;
    s.text(PAD_L, yb, "Two point substitutions in the rebuilt TMCO6 groove left the predicted interface unchanged, so they are", 11.5, "start", "fill='#333'");
    s.text(PAD_L, yb + 16.0, "not a usable experiment. The pair survives unchanged inside the LSm ring in every configuration tested.", 11.5, "start", "fill='#333'");
    s.save(out);
}

/// Figure 3: the NLS-groove positions are rebuilt in vertebrates and ancestral outside them.
fn fig3(cons_path: &str, out: &str) {
    // pos -> (aa, identity in all species, identity in non-vertebrates)
    let mut m: BTreeMap<usize, (char, f64, f64)> = BTreeMap::new();
    for l in fs::read_to_string(cons_path).unwrap().lines().skip(1) {
        let v: Vec<&str> = l.split('\t').collect();
        if v.len() < 4 {
            continue;
        }
        let p: usize = v[0].parse().unwrap();
        m.insert(p, (v[1].chars().next().unwrap(), v[2].parse().unwrap_or(0.0), v[3].parse().unwrap_or(0.0)));
    }
    let sites: Vec<(usize, &str)> = vec![
        (130, "major site"), (134, "major site"), (171, "major site"),
        (213, "major site"), (217, "major site"),
        (362, "minor site"), (366, "minor site"),
    ];
    let h_plot = 200.0;
    let mut s = Svg::new(h_plot + 170.0);
    s.text(PAD_L, 24.0, "Figure 3. The major NLS-binding site is rebuilt in vertebrates; the minor site is not", 15.0, "start", "font-weight='bold'");
    s.text(PAD_L, 40.0, "Per-position identity to human TMCO6 across 751 species with a PF30746 family member.", 11.5, "start", "fill='#555'");
    let y0 = 60.0;
    let y1 = y0 + h_plot;
    let plot_w = W - PAD_L - PAD_R;
    let gw = plot_w / sites.len() as f64;
    s.line(PAD_L, y1, W - PAD_R, y1, "#333", 1.0);
    s.line(PAD_L, y0, PAD_L, y1, "#333", 1.0);
    for t in 0..=5 {
        let v = t as f64 * 0.2;
        let y = y1 - h_plot * v;
        s.line(PAD_L - 4.0, y, PAD_L, y, "#333", 1.0);
        s.text(PAD_L - 8.0, y + 4.0, &format!("{v:.1}"), 10.5, "end", "fill='#333'");
        if t > 0 {
            s.line(PAD_L, y, W - PAD_R, y, "#eee", 1.0);
        }
    }
    s.text(20.0, (y0 + y1) / 2.0, "identity to human residue", 12.0, "middle", &format!("transform='rotate(-90 20 {:.1})'", (y0 + y1) / 2.0));
    for (i, (pos, kind)) in sites.iter().enumerate() {
        let Some(&(aa, all, nonvert)) = m.get(pos) else { continue };
        let bx = PAD_L + gw * i as f64;
        let bw = gw * 0.34;
        for (j, (v, fill)) in [(all, "#1a4b8c"), (nonvert, "#c46a2f")].iter().enumerate() {
            let x = bx + gw * 0.13 + j as f64 * bw;
            let hh = h_plot * v;
            s.rect(x, y1 - hh, bw * 0.9, hh, fill, "");
            s.text(x + bw * 0.45, y1 - hh - 6.0, &format!("{v:.2}"), 10.0, "middle", "fill='#222'");
        }
        s.text(bx + gw / 2.0, y1 + 18.0, &format!("{aa}{pos}"), 12.0, "middle", "font-weight='bold' fill='#222'");
        s.text(bx + gw / 2.0, y1 + 33.0, kind, 10.2, "middle", "fill='#666'");
    }
    // legend
    let ly = y1 + 56.0;
    s.rect(PAD_L, ly, 13.0, 13.0, "#1a4b8c", "");
    s.text(PAD_L + 19.0, ly + 11.0, "all 751 species", 11.5, "start", "fill='#333'");
    s.rect(PAD_L + 160.0, ly, 13.0, 13.0, "#c46a2f", "");
    s.text(PAD_L + 179.0, ly + 11.0, "non-vertebrate family members only (n = 65)", 11.5, "start", "fill='#333'");
    let yb = ly + 40.0;
    s.text(PAD_L, yb, "At the five major-site positions the human residue is conserved across vertebrates but absent outside them,", 11.5, "start", "fill='#333'");
    s.text(PAD_L, yb + 16.0, "where the ancestral importin W/N ladder is retained. The minor site W362/N366 is invariant in both groups.", 11.5, "start", "fill='#333'");
    s.save(out);
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let out = a.get(1).map(|s| s.as_str()).unwrap_or("preprint/figures");
    fs::create_dir_all(out).unwrap();
    fig1("data/tmco6/AF-TMCO6.pdb", &format!("{out}/fig1_burial.svg"));
    fig2("af3/results", &format!("{out}/fig2_iptm.svg"));
    fig3("data/tmco6/cons_TMCO6.tsv", &format!("{out}/fig3_conservation.svg"));
    eprintln!("wrote 3 SVG figures to {out}");
}

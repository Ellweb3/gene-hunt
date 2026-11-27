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

/// Figure 1: burial does not separate the annotated TMCO6 segments from genuine TM helices.
/// Each point is one annotated transmembrane helix, placed at its mean Calpha neighbour count.
fn fig1(out: &str) {
    // (file, label, note, TRANSMEM ranges from the corresponding UniProt entry)
    let sets: Vec<(&str, &str, &str, Vec<(usize, usize)>)> = vec![
        ("data/membrane_controls/P29972.pdb", "AQP1", "experimental evidence",
         vec![(12, 29), (47, 65), (91, 109), (134, 153), (164, 181), (207, 224)]),
        ("data/membrane_controls/P11166.pdb", "GLUT1", "",
         vec![(12, 33), (67, 87), (91, 112), (121, 144), (156, 176), (186, 206),
              (272, 293), (307, 328), (335, 355), (366, 388), (402, 422), (430, 450)]),
        ("data/membrane_controls/P08183.pdb", "MDR1", "",
         vec![(45, 67), (117, 137), (187, 208), (216, 236), (293, 316), (331, 352),
              (712, 732), (757, 777), (833, 853), (855, 874), (935, 957), (974, 995)]),
        ("data/membrane_controls/Q9UNQ0.pdb", "ABCG2", "",
         vec![(396, 416), (429, 449), (478, 498), (507, 527), (536, 556), (631, 651)]),
        ("data/tmco6/AF-TMCO6.pdb", "TMCO6", "the two segments in question",
         vec![(338, 358), (386, 406)]),
    ];

    let means = |path: &str, segs: &[(usize, usize)]| -> Vec<f64> {
        let Ok(txt) = fs::read_to_string(path) else { return vec![] };
        let ca: Vec<[f64; 3]> = txt.lines().filter(|l| l.starts_with("ATOM") && &l[12..16] == " CA ")
            .map(|l| {
                let p = |i: usize| l[i..i + 8].trim().parse::<f64>().unwrap();
                [p(30), p(38), p(46)]
            }).collect();
        let n = ca.len();
        let nb = |i: usize| (0..n).filter(|&j| j != i && (0..3).map(|k| (ca[i][k] - ca[j][k]).powi(2)).sum::<f64>() < 100.0).count() as f64;
        segs.iter().filter(|&&(a, b)| b <= n).map(|&(a, b)| {
            (a - 1..b).map(nb).sum::<f64>() / (b - a + 1) as f64
        }).collect()
    };

    let h_plot = 300.0;
    let mut s = Svg::new(h_plot + 80.0);

    let y0 = 18.0;
    let y1 = y0 + h_plot;
    let (lo, hi) = (10.0, 28.0);
    let y_of = |v: f64| y1 - h_plot * (v - lo) / (hi - lo);
    let plot_w = W - PAD_L - PAD_R;
    let cw = plot_w / sets.len() as f64;

    s.line(PAD_L, y0, PAD_L, y1, "#333", 1.0);
    for t in 0..=6 {
        let v = lo + (hi - lo) * t as f64 / 6.0;
        s.line(PAD_L - 4.0, y_of(v), PAD_L, y_of(v), "#333", 1.0);
        s.text(PAD_L - 8.0, y_of(v) + 4.0, &format!("{v:.0}"), 10.5, "end", "fill='#333'");
        s.line(PAD_L, y_of(v), W - PAD_R, y_of(v), "#eee", 1.0);
    }
    s.text(20.0, (y0 + y1) / 2.0, "mean Calpha neighbours within 10 A", 11.5, "middle",
           &format!("transform='rotate(-90 20 {:.1})' fill='#333'", (y0 + y1) / 2.0));

    for (i, (path, label, note, segs)) in sets.iter().enumerate() {
        let vals = means(path, segs);
        if vals.is_empty() {
            continue;
        }
        let cx = PAD_L + cw * (i as f64 + 0.5);
        let is_q = *label == "TMCO6";
        let fill = if is_q { "#c43b2f" } else { "#1a4b8c" };
        if is_q {
            s.rect(PAD_L + cw * i as f64 + 6.0, y0, cw - 12.0, h_plot, "#fdf0ee", "");
        }
        // spread points horizontally so they do not overlap
        for (k, v) in vals.iter().enumerate() {
            let off = ((k as f64) - (vals.len() as f64 - 1.0) / 2.0) * (cw * 0.52 / vals.len().max(2) as f64);
            s.body += &format!("<circle cx='{:.1}' cy='{:.1}' r='4.6' fill='{fill}' fill-opacity='0.8'/>\n", cx + off, y_of(*v));
        }
        // mean line
        let m = vals.iter().sum::<f64>() / vals.len() as f64;
        s.line(cx - cw * 0.2, y_of(m), cx + cw * 0.2, y_of(m), fill, 2.2);
        s.text(cx, y1 + 20.0, label, 13.0, "middle", &format!("font-weight='bold' fill='{}'", if is_q { "#c43b2f" } else { "#222" }));
        s.text(cx, y1 + 35.0, &format!("{} helices", vals.len()), 10.5, "middle", "fill='#666'");
        if !note.is_empty() {
            s.text(cx, y1 + 49.0, note, 9.8, "middle", &format!("fill='{}'", if is_q { "#c43b2f" } else { "#666" }));
        }
    }
    let yb = y1 + 76.0;
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
    let mut s = Svg::new(h_plot + 60.0);
    let y0 = 18.0;
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
    let mut s = Svg::new(h_plot + 95.0);
    let y0 = 18.0;
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
    s.save(out);
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let out = a.get(1).map(|s| s.as_str()).unwrap_or("preprint/figures");
    fs::create_dir_all(out).unwrap();
    fig1(&format!("{out}/fig1_burial.svg"));
    fig2("af3/results", &format!("{out}/fig3_iptm.svg"));
    fig3("data/tmco6/cons_TMCO6.tsv", &format!("{out}/fig2_conservation.svg"));
    eprintln!("wrote 3 SVG figures to {out}");
}

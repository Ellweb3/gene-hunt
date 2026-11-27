//! Поиск белка в ДНК: трансляция в 6 рамках + Смит–Уотерман (BLOSUM62, аффинные штрафы 11/1).

const AA: &[u8] = b"ARNDCQEGHILKMFPSTWYVBZX*";
#[rustfmt::skip]
const BLOSUM62: [[i8; 24]; 24] = [
 [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4],
 [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4],
 [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4],
 [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4],
 [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
 [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4],
 [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
 [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4],
 [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4],
 [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4],
 [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4],
 [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4],
 [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4],
 [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4],
 [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
 [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4],
 [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4],
 [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4],
 [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4],
 [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4],
 [-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4],
 [-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
 [ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4],
 [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1],
];
const GAP_OPEN: i32 = 11;
const GAP_EXT: i32 = 1;
/// Параметры Карлина–Альтшуля для BLOSUM62 11/1 (как в BLAST).
pub const LAMBDA: f64 = 0.267;
pub const K: f64 = 0.041;

pub fn encode(p: &[u8]) -> Vec<u8> {
    p.iter().map(|c| AA.iter().position(|a| a == &c.to_ascii_uppercase()).unwrap_or(22) as u8).collect()
}

fn codon(c: &[u8]) -> u8 {
    const B: &[u8] = b"TCAG";
    const T: &[u8] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let idx = |x: u8| B.iter().position(|b| *b == x.to_ascii_uppercase());
    match (idx(c[0]), idx(c[1]), idx(c[2])) {
        (Some(a), Some(b), Some(d)) => T[a * 16 + b * 4 + d],
        _ => b'X',
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|c| match c.to_ascii_uppercase() { b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b'N' }).collect()
}

/// Шесть рамок: (смещение, прямая цепь?, белок).
pub fn six_frames(dna: &[u8]) -> Vec<(usize, bool, Vec<u8>)> {
    let rc = revcomp(dna);
    let mut out = vec![];
    for (fwd, s) in [(true, dna), (false, rc.as_slice())] {
        for f in 0..3 {
            out.push((f, fwd, s[f..].chunks_exact(3).map(codon).collect()));
        }
    }
    out
}

/// Лучший локальный счёт Смита–Уотермана (Готох) и позиция конца в субъекте.
pub fn sw_best(q: &[u8], s: &[u8]) -> (i32, usize) {
    let n = s.len();
    let (mut h, mut e) = (vec![0i32; n + 1], vec![i32::MIN / 2; n + 1]);
    let (mut best, mut end) = (0, 0);
    for &qa in q {
        let row = &BLOSUM62[qa as usize];
        let (mut diag, mut f) = (0i32, i32::MIN / 2);
        h[0] = 0;
        for j in 1..=n {
            e[j] = (e[j] - GAP_EXT).max(h[j] - GAP_OPEN);
            f = (f - GAP_EXT).max(h[j - 1] - GAP_OPEN);
            let v = (diag + sc(row, s[j - 1])).max(e[j]).max(f).max(0);
            diag = h[j];
            h[j] = v;
            if v > best { best = v; end = j; }
        }
    }
    (best, end)
}

pub fn bits(score: i32) -> f64 { (LAMBDA * score as f64 - K.ln()) / std::f64::consts::LN_2 }
pub fn evalue(score: i32, m: usize, n: usize) -> f64 { K * m as f64 * n as f64 * (-LAMBDA * score as f64).exp() }

/// Выравнивание с обратным проходом: (счёт, q_start, q_end, s_start, s_end, строки q/s).
pub struct Aln { pub score: i32, pub q: (usize, usize), pub s: (usize, usize), pub qa: Vec<u8>, pub sa: Vec<u8> }

pub fn sw_align(q: &[u8], s: &[u8]) -> Aln {
    let (m, n) = (q.len(), s.len());
    let w = n + 1;
    const NEG: i32 = i32::MIN / 2;
    let mut h = vec![0i32; (m + 1) * w];
    let mut e = vec![NEG; (m + 1) * w];
    let mut f = vec![NEG; (m + 1) * w];
    let (mut best, mut bi, mut bj) = (0, 0, 0);
    for i in 1..=m {
        let row = &BLOSUM62[q[i - 1] as usize];
        for j in 1..=n {
            let k = i * w + j;
            e[k] = (e[k - w] - GAP_EXT).max(h[k - w] - GAP_OPEN);
            f[k] = (f[k - 1] - GAP_EXT).max(h[k - 1] - GAP_OPEN);
            h[k] = (h[k - w - 1] + sc(row, s[j - 1])).max(e[k]).max(f[k]).max(0);
            if h[k] > best { best = h[k]; bi = i; bj = j; }
        }
    }
    let (mut i, mut j, mut state) = (bi, bj, 0u8); // 0 = H, 1 = E (пропуск в s), 2 = F (пропуск в q)
    let (mut qa, mut sa) = (vec![], vec![]);
    while i > 0 && j > 0 {
        let k = i * w + j;
        match state {
            0 => {
                if h[k] == 0 { break; }
                if h[k] == e[k] { state = 1; continue; }
                if h[k] == f[k] { state = 2; continue; }
                qa.push(AA[q[i - 1] as usize]); sa.push(AA[s[j - 1] as usize]); i -= 1; j -= 1;
            }
            1 => { qa.push(AA[q[i - 1] as usize]); sa.push(b'-'); if e[k] == h[k - w] - GAP_OPEN { state = 0; } i -= 1; }
            _ => { qa.push(b'-'); sa.push(AA[s[j - 1] as usize]); if f[k] == h[k - 1] - GAP_OPEN { state = 0; } j -= 1; }
        }
    }
    qa.reverse(); sa.reverse();
    Aln { score: best, q: (i, bi), s: (j, bj), qa, sa }
}

/// Штраф за стоп-кодон в субъекте: выравнивание не должно проходить сквозь стоп (интрон/псевдоген).
pub const STOP_PENALTY: i32 = -25;
#[inline]
fn sc(row: &[i8; 24], s: u8) -> i32 { if s == 23 { STOP_PENALTY } else { row[s as usize] as i32 } }

/// Балл BLOSUM62 для двух закодированных аминокислот.
pub fn blosum(a: u8, b: u8) -> i32 { BLOSUM62[a as usize][b as usize] as i32 }

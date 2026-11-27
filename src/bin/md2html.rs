//! Render the manuscript Markdown as a print-ready HTML page.
//! md2html <manuscript.md> <out.html>
//!
//! Supports the subset of Markdown the manuscript uses: ATX headings, paragraphs, bold,
//! italic, inline code, links, fenced code blocks, blockquotes, ordered and unordered lists,
//! pipe tables, and `---` rules. A `<!--FIG n file caption-->` comment inserts a figure.
use std::fs;

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}

/// inline formatting, applied after escaping
fn inline(s: &str) -> String {
    let mut out = esc(s);
    // links [text](url) before emphasis, so underscores in URLs survive.
    // Walk by char boundaries: byte-wise iteration would split multi-byte UTF-8.
    let mut res = String::new();
    let mut i = 0;
    while i < out.len() {
        let ch = out[i..].chars().next().unwrap();
        if ch == '[' {
            if let Some(close) = out[i..].find("](") {
                let text_end = i + close;
                if let Some(paren) = out[text_end + 2..].find(')') {
                    let url = &out[text_end + 2..text_end + 2 + paren];
                    let text = &out[i + 1..text_end];
                    res += &format!("<a href=\"{url}\">{text}</a>");
                    i = text_end + 3 + paren;
                    continue;
                }
            }
        }
        res.push(ch);
        i += ch.len_utf8();
    }
    out = res;
    // `code`
    while let (Some(a), Some(_)) = (out.find('`'), out.rfind('`')) {
        let Some(b2) = out[a + 1..].find('`') else { break };
        let code = out[a + 1..a + 1 + b2].to_string();
        out = format!("{}<code>{}</code>{}", &out[..a], code, &out[a + b2 + 2..]);
    }
    // **bold** then *italic*
    for (mark, tag) in [("**", "strong"), ("*", "em")] {
        loop {
            let Some(a) = out.find(mark) else { break };
            let Some(b2) = out[a + mark.len()..].find(mark) else { break };
            let inner = out[a + mark.len()..a + mark.len() + b2].to_string();
            out = format!("{}<{tag}>{}</{tag}>{}", &out[..a], inner, &out[a + 2 * mark.len() + b2..]);
        }
    }
    out
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let md = fs::read_to_string(&a[1]).unwrap();
    let mut body = String::new();
    let mut lines = md.lines().peekable();
    let mut para: Vec<String> = Vec::new();
    let mut fig_n = 0;

    macro_rules! flush {
        () => {
            if !para.is_empty() {
                body += &format!("<p>{}</p>\n", inline(&para.join(" ")));
                para.clear();
            }
        };
    }

    while let Some(l) = lines.next() {
        let t = l.trim_end();
        if t.trim().is_empty() {
            flush!();
            continue;
        }
        if let Some(rest) = t.strip_prefix("```") {
            flush!();
            let lang = rest.trim().to_string();
            let mut code = String::new();
            for c in lines.by_ref() {
                if c.trim_start().starts_with("```") {
                    break;
                }
                code += &esc(c);
                code.push('\n');
            }
            body += &format!("<pre class=\"lang-{lang}\"><code>{code}</code></pre>\n");
            continue;
        }
        if t.trim() == "---" {
            flush!();
            body += "<hr/>\n";
            continue;
        }
        if let Some(h) = t.strip_prefix("##### ").map(|x| (5, x)).or(t.strip_prefix("#### ").map(|x| (4, x)))
            .or(t.strip_prefix("### ").map(|x| (3, x))).or(t.strip_prefix("## ").map(|x| (2, x)))
            .or(t.strip_prefix("# ").map(|x| (1, x)))
        {
            flush!();
            body += &format!("<h{0}>{1}</h{0}>\n", h.0, inline(h.1));
            continue;
        }
        if t.trim_start().starts_with("> ") || t.trim() == ">" {
            flush!();
            let mut q = vec![t.trim_start().trim_start_matches('>').trim().to_string()];
            while let Some(p) = lines.peek() {
                if p.trim_start().starts_with('>') {
                    q.push(lines.next().unwrap().trim_start().trim_start_matches('>').trim().to_string());
                } else {
                    break;
                }
            }
            body += &format!("<blockquote><p>{}</p></blockquote>\n", inline(&q.join(" ")));
            continue;
        }
        // tables
        if t.trim_start().starts_with('|') {
            flush!();
            let mut rows: Vec<String> = vec![t.to_string()];
            while let Some(p) = lines.peek() {
                if p.trim_start().starts_with('|') {
                    rows.push(lines.next().unwrap().to_string());
                } else {
                    break;
                }
            }
            let cells = |r: &str| -> Vec<String> {
                r.trim().trim_matches('|').split('|').map(|c| inline(c.trim())).collect()
            };
            body += "<table>\n<thead><tr>";
            for c in cells(&rows[0]) {
                body += &format!("<th>{c}</th>");
            }
            body += "</tr></thead>\n<tbody>\n";
            for r in rows.iter().skip(2) {
                body += "<tr>";
                for c in cells(r) {
                    body += &format!("<td>{c}</td>");
                }
                body += "</tr>\n";
            }
            body += "</tbody></table>\n";
            continue;
        }
        // lists
        let is_ul = t.trim_start().starts_with("- ");
        let is_ol = t.trim_start().split_once(". ").is_some_and(|(n, _)| !n.is_empty() && n.trim().chars().all(|c| c.is_ascii_digit()));
        if is_ul || is_ol {
            flush!();
            let tag = if is_ul { "ul" } else { "ol" };
            let mut items: Vec<String> = Vec::new();
            let mut cur = String::new();
            let mut push_line = |s: &str, items: &mut Vec<String>, cur: &mut String| {
                let st = s.trim_start();
                let new = if let Some(r) = st.strip_prefix("- ") { Some(r.to_string()) }
                    else { st.split_once(". ").and_then(|(n, r)| (!n.is_empty() && n.chars().all(|c| c.is_ascii_digit())).then(|| r.to_string())) };
                match new {
                    Some(r) => {
                        if !cur.is_empty() { items.push(std::mem::take(cur)); }
                        *cur = r;
                    }
                    None => { cur.push(' '); cur.push_str(st); }
                }
            };
            push_line(t, &mut items, &mut cur);
            while let Some(p) = lines.peek() {
                let pt = p.trim_start();
                if pt.is_empty() { break; }
                let cont = p.starts_with("  ") || p.starts_with('\t');
                let item = pt.starts_with("- ") || pt.split_once(". ").is_some_and(|(n, _)| !n.is_empty() && n.chars().all(|c| c.is_ascii_digit()));
                if item || cont {
                    let line = lines.next().unwrap().to_string();
                    push_line(&line, &mut items, &mut cur);
                } else { break; }
            }
            if !cur.is_empty() { items.push(cur); }
            body += &format!("<{tag}>\n");
            for it in items {
                body += &format!("<li>{}</li>\n", inline(&it));
            }
            body += &format!("</{tag}>\n");
            continue;
        }
        if let Some(fig) = t.trim().strip_prefix("<!--FIG ").and_then(|x| x.strip_suffix("-->")) {
            flush!();
            fig_n += 1;
            let (file, cap) = fig.trim().split_once(' ').unwrap_or((fig.trim(), ""));
            body += &format!("<figure><img src=\"{file}\" alt=\"Figure {fig_n}\"/><figcaption>{}</figcaption></figure>\n", inline(cap));
            continue;
        }
        para.push(t.trim().to_string());
    }
    flush!();

    let css = r#"
@page { size: A4; margin: 20mm 18mm; }
body { font-family: Georgia, 'Times New Roman', serif; font-size: 10.5pt; line-height: 1.5;
       color: #111; max-width: 170mm; margin: 0 auto; }
h1 { font-size: 17pt; line-height: 1.3; margin: 0 0 4pt; }
h2 { font-size: 13pt; margin: 20pt 0 6pt; border-bottom: 1px solid #ccc; padding-bottom: 3pt; }
h3 { font-size: 11.5pt; margin: 14pt 0 4pt; }
p { margin: 0 0 7pt; text-align: justify; }
table { border-collapse: collapse; width: 100%; margin: 8pt 0 12pt; font-size: 9pt;
        font-family: Helvetica, Arial, sans-serif; page-break-inside: avoid; }
th, td { border: 1px solid #bbb; padding: 3pt 5pt; text-align: left; vertical-align: top; }
th { background: #f0f0f0; }
pre { background: #f6f6f6; border: 1px solid #ddd; padding: 6pt; font-size: 8.5pt;
      overflow-x: auto; page-break-inside: avoid; }
code { font-family: 'SF Mono', Menlo, Consolas, monospace; font-size: 9pt; }
pre code { font-size: 8.5pt; }
blockquote { margin: 8pt 0; padding: 6pt 10pt; background: #fbf6ee; border-left: 3px solid #d9a441; }
blockquote p { margin: 0; }
figure { margin: 12pt 0; page-break-inside: avoid; text-align: center; }
figure img { width: 100%; border: 1px solid #ddd; }
figcaption { font-size: 8.5pt; color: #444; text-align: left; margin-top: 4pt;
             font-family: Helvetica, Arial, sans-serif; }
ul, ol { margin: 0 0 8pt; padding-left: 18pt; }
li { margin-bottom: 3pt; }
hr { border: none; border-top: 1px solid #ccc; margin: 14pt 0; }
a { color: #11456e; text-decoration: none; }
h2, h3 { page-break-after: avoid; }
"#;
    let html = format!("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"/>\n<title>TMCO6-ZNF511 preprint</title>\n<style>{css}</style>\n</head>\n<body>\n{body}</body></html>\n");
    fs::write(&a[2], html).unwrap();
    eprintln!("wrote {} ({} figures)", a[2], fig_n);
}

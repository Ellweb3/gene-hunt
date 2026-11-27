# Preprint

`manuscript.md` is the source; `manuscript.pdf` is what gets submitted.

## Building

No external tooling is required beyond Rust, `rsvg-convert` and Chrome.

```bash
cargo build --release

# figures: SVG from the project's own result files, then PNG for embedding
./target/release/figures preprint/figures
for f in preprint/figures/*.svg; do rsvg-convert -f png -w 1600 "$f" -o "${f%.svg}.png"; done

# manuscript: Markdown -> print-ready HTML -> PDF
./target/release/md2html preprint/manuscript.md preprint/manuscript.html
cd preprint && "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
  --headless --disable-gpu --no-pdf-header-footer \
  --print-to-pdf=manuscript.pdf manuscript.html
```

Figures are generated from the result files rather than from hard-coded numbers, so they
change if the underlying analysis changes. Figure 1 additionally needs the membrane-protein
control models in `data/membrane_controls/` (see `data/README.md`).

## Figures

| file | in the manuscript | what it shows |
|---|---|---|
| `fig1_burial.svg` | Figure 1 | burial does not separate the annotated TMCO6 segments from genuine transmembrane helices — a negative control |
| `fig2_conservation.svg` | Figure 2 | the major NLS-binding site is rebuilt in vertebrates, the minor site is not |
| `fig3_iptm.svg` | Figure 3 | the pair interface is reproducible; the predicted point mutants are not |

## Submission

Order matters: **bioRxiv first, Zenodo only as a fallback.** bioRxiv does not accept
manuscripts already posted on another preprint server, so depositing elsewhere first would
foreclose it.

bioRxiv notes to keep in mind:
- An affiliation must be supplied at registration. Whether "independent researcher" is
  accepted is not stated on their pages and has not been confirmed.
- No fee, no endorsement requirement, screening typically 24-48 hours, about 5% of
  submissions are not posted.
- Their FAQ welcomes computational papers but excludes "simple automated/computational
  analyses of public data, molecular modeling, sequence alignments and results of database
  searches". This manuscript is written to sit on the other side of that line: it reports
  controls, a negative control, and a prediction of ours that failed.
- A single PDF is the simplest accepted format.

# gene-hunt

**Finding under-studied human proteins for which enough open data already exists to reach a
substantive conclusion without a lab of one's own.**

Starting from the 5501 proteins in the Tdark class (next to no publications, no known
function), one pair — **TMCO6 and ZNF511** — was selected and analysed in depth. All
computation is in Rust; all data comes from public databases.

> ⚠️ Everything here is computational. These are hypotheses worth testing at the bench, not
> established facts. Where a prediction failed, it is marked as such — see sections 4, 7
> and 9.3 of the main analysis.

---

## The finding

**TMCO6 is not a membrane protein. It is a relative of importin-alpha.** The gene name
("transmembrane and coiled-coil domains 6") and the UniProt annotation "multi-pass membrane"
do not match the structure: Foldseek retrieves importins-alpha over the full length of the
protein (E ~1e-12), and both annotated "transmembrane" helices are buried in the hydrophobic
core of the solenoid (pLDDT 94-97, 0-1 charged residues, ~19 Calpha neighbours within 10 A).
The similarity is invisible at the sequence level, at only ~19% identity.

**TMCO6 and ZNF511 form a stable pair.** Three independent lines of evidence:

| source | what it shows |
|---|---|
| DepMap, 1208 cell lines | reciprocal best partners, r = 0.60; the next-best partner of each is only 0.24 |
| BioPlex 3.0, two cell lines | they co-precipitate with each other, and both with an LSm ring |
| AlphaFold, two independent models | ipTM 0.79 and 0.84, interface of 133 / 90 residues |

**The pair appears to be associated with the nuclear LSm2-8 ring, the core of the U6 snRNP.**
The key independent evidence is published experimental data from another kingdom: in an
AtLSM8 purification from *Arabidopsis* (Golisz et al. 2013, PMID 23620288) the ZNF511
homolog was recovered alongside SART3/Prp24 and Brr2, and was absent from the AtLSM1a
purifications.

**Possible significance.** In North Country Cheviot sheep, a loss-of-function variant in
TMCO6 causes inherited motor neuron disease (PMID 37488055). If this module does serve U6,
that phenotype would sit alongside the other snRNP-assembly disorders. That comparison is
ours, not a conclusion of the paper.

📄 **[Full analysis: TMCO6_ZNF511.md](TMCO6_ZNF511.md)** — nine sections, including the
hypotheses that were refuted.
✉️ **[Draft letter to a lab](letter_draft.md)** — who to contact and what to propose testing.

---

## How the candidates were found

1. **Starting set** — 5501 Tdark proteins from Pharos (the IDG consortium).
2. **Data layers** — publication counts (NCBI gene2pubmed), cell line dependency
   (DepMap Chronos 26Q1), physical interactions (BioPlex 3.0), subcellular localisation
   (Human Protein Atlas).
3. **Co-dependency** — Pearson correlation of dependency profiles across ~1100 lines, for
   genes measured in at least 500 of them.
4. **Removing artifacts** — this turned out to be the bulk of the work:
   - cytoband neighbours excluded (copy-number artifact);
   - multi-copy families excluded (OR, KRTAP, GOLGA6/8, AGAP, POTE, SPDYE and others): a
     guide cutting many loci at once makes a gene look falsely essential;
   - genes measured in few lines dropped, as they produce spurious correlations;
   - manually removed: subunits of known complexes that are formally classed as unstudied
     (mitochondrial ribosome, SAGA, minor spliceosome), and genes with recent functional
     papers (TTC27, GPATCH1/LENG1, both 2026).

### Finalists

| gene | DepMap | co-dependency | BioPlex | status in the literature |
|---|---|---|---|---|
| **TMCO6 + ZNF511** | mean −0.11, 7-8% of lines < −0.5 | **reciprocal best, r = 0.60** | with each other and with an LSm ring | link not described → **taken forward** |
| **SS18L2** | mean −2.59, required by 99% of lines | weak (0.32) | the entire BAF (SWI/SNF) complex | listed as a BAF subunit, but why it is indispensable is unknown |
| **C19orf25** | mean −0.13, 12% of lines | BNIP1 0.52, NBAS 0.47 | all 8 subunits of the syntaxin-18/NRZ complex | mentioned in systematic surveys, no paper of its own |

---

## Repository layout

```
TMCO6_ZNF511.md      main analysis, nine sections
letter_draft.md      draft letter to a lab
src/bin/             12 analysis programs in Rust
data/                computed results (raw downloads are not stored)
af3/                 AlphaFold Server jobs and results
```

### Programs

| program | what it does |
|---|---|
| `main.rs` | joins Tdark x DepMap x BioPlex x HPA and ranks candidates |
| `pair.rs` | DepMap profile of a pair: gene panel, tumour types, co-dependency module |
| `nullpc.rs` | **artifact control**: removes screening noise estimated from 193 non-expressed olfactory receptor genes |
| `plddt.rs` | AlphaFold model profile: pLDDT, helicity, hydrophobicity, burial |
| `map_sites.rs` | transfers functional positions from a reference protein via structural alignment |
| `conserve.rs` | per-residue conservation across orthologs, one sequence per species |
| `cooccur.rs` | co-occurrence of two families across 2087 proteomes |
| `rbh.rs` | **orthology control**: best hit against a panel of reference proteins |
| `ident.rs` | identity to two references at once |
| `njtree.rs` | neighbor-joining tree with bootstrap and explicit tests of named clades |
| `motif_presence.rs` | association between a motif and a partner's presence, Fisher's exact test |
| `interface.rs` | chain interface in a complex model |
| `af3_report.rs` | parses AlphaFold Server output: per-chain-pair ipTM and contacts |

Chains and proteins are identified by sequence rather than by file order, so a wild type
cannot be confused with a variant.

---

## Reproducing

Requires Rust (edition 2024) and `curl`.

```bash
cargo build --release

# Candidate search (large inputs go to a temporary directory)
S=$(mktemp -d)
curl -sL https://ndownloader.figshare.com/files/67214582 -o $S/gene_effect.csv
curl -s https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz -o $S/gene2pubmed.gz
cargo run --release -- $S/gene_effect.csv $S/gene2pubmed.gz

# Analysis of the pair
./target/release/pair      $S/gene_effect.csv data/tmco6/ach_disease.tsv data/tmco6/dm
./target/release/nullpc    $S/gene_effect.csv 10 TMCO6 ZNF511 LSM8 USB1 LARP7 SNUPN
./target/release/plddt     data/tmco6/AF-TMCO6.pdb 1-60 338-358 386-406
./target/release/interface data/tmco6/AF-complex-TMCO6-ZNF511.pdb 130 134 171 175 213 217
./target/release/af3_report af3/results af3/inputs.fa
```

Remaining downloads are listed in **[data/README.md](data/README.md)**.

---

## Sources and terms of use

This work rests on open data from **DepMap** (Broad Institute), **BioPlex 3.0**,
**Human Protein Atlas**, **UniProt**, **InterPro/Pfam**, **Pharos/IDG**, **NCBI Gene**,
**Cellosaurus**, **Ensembl**, **Foldseek**, **AlphaFold DB** and **Europe PMC**. Each source
carries its own terms; check them at the source.

Raw downloads from those databases are not stored here: each carries its own licence, and
together they run to some 60 MB. `data/README.md` lists the commands that produced them.

> **The models in `af3/results/` were produced on AlphaFold Server and are distributed under
> its terms: non-commercial use only.** The full text sits next to each model
> (`af3/results/*/terms_of_use.md`) and is published at
> [alphafoldserver.com/output-terms](https://alphafoldserver.com/output-terms). Only the
> top-ranked models and their confidence summaries are stored; the rest is reproduced by
> re-running `af3/alphafold_server_jobs.json`.

The code under `src/` was written for this work — use it freely.

If any of these conclusions is already known, is wrong, or rests on a flaw you can see,
please say so. That would be useful.

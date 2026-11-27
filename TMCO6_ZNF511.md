# TMCO6–ZNF511: an importin-like module at the nuclear LSm2-8 ring

Detailed analysis of finalist #1 from the [README](README.md). Date: 2026-11-27.
Data: `data/tmco6/`.

> Everything here is computational. These are hypotheses worth testing at the bench, not
> established facts. Where a prediction failed, it is marked as such — see sections 4, 7
> and 9.3.

## Summary

TMCO6 (493 aa) is **not a membrane protein**. It is a distant relative of importin-alpha,
matching it along its entire length. ZNF511 (252 aa, three C2H2 fingers) is its constant
partner. The pair associates with an LSm ring. In plants, the ZNF511 homolog was recovered
in an **LSm8** purification together with SART3/Prp24 and Brr2, which points to the nuclear
**LSm2-8** ring, the core of the U6 snRNP. Genetic dependencies in DepMap outline the same
module: USB1, LARP7, SNUPN, LSM8.

Two independent predictions (AlphaFold DB, ipTM 0.79; our own AlphaFold 3 run, ipTM 0.84)
place the N-terminus of ZNF511 in the NLS-binding groove of TMCO6. Of the LSm ring, it is
ZNF511 that makes contact; TMCO6 is held via ZNF511 (section 9).

The main site of that groove is rebuilt in vertebrates and conserved from lamprey to human.
The ZNF511 motif that occupies the groove arose in vertebrates at the same time as the
rebuild. TMCO6 itself is a fast-evolving member of an ancient family of importin-like
proteins, PF30746 (sections 7-8).

A loss-of-function variant in TMCO6 causes recessive motor neuron disease in sheep, a
picture also seen in human snRNP-assembly disorders. Working hypothesis:
**TMCO6-ZNF511 is an importin-like chaperone for the U6/LSm2-8 snRNP.**

## 1. DepMap: the link is real, not a screening artifact (`src/bin/pair.rs`, `src/bin/nullpc.rs`)

| test | result |
|---|---|
| r(TMCO6, ZNF511), 1208 cell lines | **0.601** |
| after removing 10 principal components of a null set (193 non-expressed OR genes) | **0.612** |
| OR gene pairs after the same correction: median / 99.9% / max | −0.005 / 0.20 / 0.34 |
| next-best partners of TMCO6 after correction | DSCC1 0.20, **LARP7 0.18**, **SNUPN 0.18**, PRRC1, RNF130, ARMC5, POLE3, KIF2C, **USB1 0.16** |
| LSm subunits: correlation with TMCO6 / ZNF511 | **LSM8 0.15/0.17**, LSM3 0.11/0.14, LSM6 0.11, LSM1 **0.00**/−0.01 |

The correction matters: olfactory receptor genes are not expressed in cancer lines, so any
shared structure in their dependency profiles is screening noise. Removing it left the pair
untouched while collapsing the null distribution.

The partners that survive (`data/tmco6/dm_panel.tsv`) all sit in the same process:
- **SNUPN** (snurportin) — nuclear import of snRNPs;
- **USB1** — 3'-end processing of U6;
- **LARP7** — 7SK RNP and U6 modification.

A caveat worth stating: the first version of this list returned CLEC20A, CD8B and similar
genes. Those are measured in only ~70 lines — small-sample noise, not biology. The threshold
is now ≥1000 lines.

Tumor type does not matter much: the dependency is spread across tissues, with only a weak
shift in glioblastoma (t = −2.2; `dm_lineage.tsv`). The most dependent lines
(`dm_toplines.tsv`) are melanoma, prostate carcinoma and esophageal squamous cell carcinoma,
reaching Chronos −1.4.

## 2. Physical association (BioPlex 3.0, two cell lines)

| | 293T | HCT116 |
|---|---|---|
| TMCO6 | ZNF511, LSM4, LSM5, LSM7, DPP9 | ZNF511, LSM6, LSM7, DPP9, PSME4 |
| ZNF511 | TMCO6, LSM3-7, WAVE complex (WASF1/2, CYFIP1/2, ABI1/2, NCKAP1), KLHL15 | TMCO6, LSM4/6/7, KLHL15 |

LSM4 and LSM7 baits recover both proteins reciprocally. LSM1 and LSM8 are absent from
BioPlex, so the assignment to LSm2-8 rests on the plant experiment (section 5) and on DepMap
(LSM8 > LSM1).

## 3. Structure of TMCO6 (AlphaFold + Foldseek, `src/bin/plddt.rs`, `src/bin/map_sites.rs`)

- **Foldseek** retrieves importin-alpha from human, yeast, plants and fungi: full-length
  match over residues 1-486, E ~1e-12, at only ~19% sequence identity. PDB hits include
  importin-alpha bound to NLS peptides (1bk6, 6wx8 and others). The similarity is invisible
  at the sequence level, which is presumably why it was missed.
- **The UniProt "transmembrane" segments** (338-358, 386-406) are two ordinary helices of the
  repeat stack. Both carry only computational evidence (ECO:0000255), and the subcellular
  location is inferred from them (ECO:0000305); the entry cites no experimental support for a
  membrane localisation.
  > ⚠ **Corrected 2026-11-27.** An earlier version of this section argued that the segments
  > are "buried, unlike a real transmembrane helix". That argument is wrong. Controls against
  > experimentally supported membrane proteins show their transmembrane helices are equally or
  > more buried in their own AlphaFold models (AQP1 20.8-23.6, GLUT1 20.3-25.6, ABCG2
  > 16.0-21.2 Calpha neighbours, versus 18.3 and 20.0 here), because they pack against one
  > another. A second attempt, the hydrophobicity of the solvent-exposed face, separates GLUT1
  > (+3.8) and ABCG2 (+4.4) from their own soluble regions but fails on AQP1 (-0.3), which is
  > almost entirely membrane-embedded; TMCO6 gives -0.4. Neither measure discriminates.
  > The reassignment rests on the fold match and on the entry's own cross-references, not on
  > these. See `src/bin/lipid_face.rs`.
- **The N-terminus (1-60)** is disordered (pLDDT 46) and basic: 26 of 60 residues are
  charged, with RRRRER and RKARR clusters. Structurally this is the counterpart of the IBB
  domain by which importin-alpha binds importin-beta.
- **The NLS-binding pocket**, mapped from mouse importin-alpha2 (P52293):

| importin site | residue -> TMCO6 |
|---|---|
| major, ARM2 W142/N146 | **R130/E134** (changed) |
| ARM3 W184/N188/D192 | Y171/**N175**/E179 |
| ARM4 W231/N235 | Y213/Q217 |
| ARM5 W273 | **W261** |
| minor, ARM7 W357/N361 | **W362/N366** (preserved) |
| ARM8 E396/W399/N403 | M400/T403/**N407** |

The major site is rebuilt while the minor site is preserved, suggesting TMCO6 holds
something other than a classical NLS peptide.

ZNF511 for comparison: C2H2 fingers at 80-169 (pLDDT 91-94); the C-terminus 170-252 is
disordered and carries basic, NLS-like clusters (KPKKSRS...KKKTK).

## 4. Phylogenetic distribution (`src/bin/cooccur.rs`) — corrected in sections 7-8

> **Correction.** A best-hit test against a panel of all human importins-alpha
> (`src/bin/rbh.rs`) showed that PF30746 members outside vertebrates are **not demonstrated
> orthologs of TMCO6**. They are equidistant from TMCO6 and from KPNA1-8: 28-33% identity,
> margins of only 1-9 bits. Unambiguous TMCO6 exists only in vertebrates, from lamprey
> onward (margins 195-760 bits). The table below therefore describes the
> "importin-like PF30746" family, not TMCO6. And the original conclusion that "TMCO6 occurs
> only together with ZNF511" is uninformative within vertebrates, where almost every species
> has both genes.

2087 eukaryotic UniProt reference proteomes with BUSCO completeness ≥90%. Labels:
PF30746 (Pfam) and IPR039258 (InterPro, ZNF511).

| | ZNF511 only | both | PF30746 only |
|---|---|---|---|
| vertebrates | 4 | 254 | 6 |
| molluscs, annelids, brachiopods | 6 | 10 | 0 |
| cnidarians | 0 | 6 | 0 |
| *Capsaspora* (Filasterea) | 0 | 1 | 0 |
| land plants | 191 | 15 | 0 |
| insects / nematodes / chelicerates | 133 / 6 / 11 | 0 | 0 |
| fungi | 965 | 2 | 1 |
| **total** | **1354** | **291** | **7** |

The family is ancient, present in both animals and plants, and was lost independently in
several lineages (insects, nematodes, arachnids, most fungi and plants) while ZNF511
persisted. The converse is almost never seen: the 7 cases of PF30746 without ZNF511 are
scattered birds and a marten, and look like annotation gaps.

In *Drosophila*, the ZNF511 homolog is **l(2)k10201**; its insertional mutant is lethal
(Berkeley Drosophila Genome Project, 1999).

## 5. Plants: the ZNF511 homolog sits on LSm8

Golisz et al. 2013, *Nucleic Acids Res* (PMID 23620288), FLAG purification of *Arabidopsis*
LSm proteins. **At5g14140**, the ZNF511 homolog, was recovered in 2 of 3 **AtLSM8**
purifications alongside SART3/Prp24 (EMB140), Brr2 (EMB1507) and SmG, and was absent from
the AtLSM1a purifications.

We initially treated this as the strongest independent evidence in the whole analysis:
published experimental data, from another kingdom, placing this protein on the nuclear LSm2-8
ring rather than on the cytoplasmic LSm1-7 ring.

> ⚠ **Weakened 2026-11-27, on the authority of the original author.** We wrote to Prof. Joanna
> Kufel, who replied that nothing is known about At5g14140 in *Arabidopsis*, that the
> interaction "has not been validated", and that its nuclear localisation there is predicted
> rather than confirmed. The co-purification is therefore an unvalidated AP-MS observation,
> which is normal for such datasets but is not the firm experimental anchor we treated it as.
> With this caveat, the assignment to the nuclear ring rests on the DepMap correlations
> (LSM8 yes, LSM1 no) and on structure prediction — that is, entirely on computation.

## 6. Phenotypes and literature

- **Sheep** (Genomics 2023, PMID 37488055). In North Country Cheviot lambs, a homozygous
  4 bp deletion in TMCO6 truncates the protein by 50%, causing recessive early-onset
  progressive ataxia with motor neuron degeneration. The picture resembles human disorders
  of snRNA assembly and processing (SMA — SMN; pontocerebellar hypoplasia — EXOSC3, TOE1).
  That comparison is ours, not a conclusion of the paper.
- **Mouse** (Cancer Res 2024, PMID 38381538). Describes TMCO6 as a membrane receptor for
  neutrophil extracellular trap DNA on CD8+ T cells, with DNA bound by the N-terminus.
  TMCO6−/− mice are viable. The membrane interpretation is inconsistent with the fold
  (section 3). DNA binding by the N-terminus is consistent with its basic, IBB-like
  character, but such binding may also be non-specific.
- **gnomAD**: both genes tolerate heterozygous loss (LoF o/e: TMCO6 0.76, ZNF511 1.05;
  pLI ~0), as expected for a recessive gene.

## 7. The complex and interface coevolution (second pass)

### Conservation (`src/bin/conserve.rs`)

751 species, best sequence per species, coverage ≥60%. The rebuilt major site is conserved
across all vertebrates:

| position | identity in vertebrates | lamprey / sharks | PF30746 outside vertebrates |
|---|---|---|---|
| R130 | 0.82 | R / Q,G | W (59 of 65) |
| E134 | 0.89 | R / E,K | N (63) |
| Y171 | 0.91 | Y | W (64) |
| Y213 | 0.69 | F/Y | F/W |
| Q217 | 0.91 | Q | N (63) |
| W362, N366 (minor site) | 1.00 | W, N | W, N |

The ancestral W/N ladder of classical importins is replaced in vertebrate TMCO6 by an
R/E/Y/Q set, and that substitution is fixed from lamprey to human — the signature of
selection rather than decay.

### The AlphaFold DB complex (`src/bin/interface.rs`)

AF-0000000211739682, computed by NVIDIA: ipTM 0.79, ipSAE 0.81, pDockQ 0.72.

- Interface: 120 TMCO6 residues and 86 ZNF511 residues, mean pLDDT ~88.
- ZNF511 runs along the entire concave groove of TMCO6. Distances to it from R130, E134,
  Y171, Y213, Q217 are 3.5-4.5 A; from W261, W362, N366, N407, 2.8-3.2 A.
- What occupies the groove is the **N-terminus of ZNF511** (1-74). Disordered in the
  monomer, in the complex it packs the way an NLS peptide packs into importin-alpha. Parts
  of ZF1/ZF2 and segment 168-181 also contribute.
- The basic N-terminus of TMCO6 (19-42, the IBB counterpart) also contacts ZNF511.
- AlphaFold DB contains only Sm-Sm pairs for LSm2-8; no TMCO6-LSm or ZNF511-LSm pairs,
  consistent with the ring being bound as a whole.

### Coevolution of the ZNF511 motif

Fraction of species in which segment 35-66 aligns (`cons_ZNF511.tsv.species`):

| group | PF30746 present | motif aligns | typical sequence |
|---|---|---|---|
| mammals / turtles / sharks | yes | 0.87 / 0.92 / 1.00 | RFRRFRFRYLVLT |
| teleosts / birds | yes | 0.97 / 0.65 | TFQRLRFRYLLFI |
| plants | PF30746-like only | 0.36 | RFPFRLAVAL |
| fungi, molluscs, cnidarians, amphioxus, insects | no | 0.00-0.14 | — |

> **Correction (section 8).** A per-species test (`src/bin/motif_presence.rs`) shows no
> association outside vertebrates. The motif is present in 24% of plants with PF30746 and in
> 37% of plants without it (p = 0.25); among molluscs with PF30746, 1 of 12 species has it.
> The overall p = 3e-59 is driven by species relatedness: the motif and the rebuilt site
> arose **once, in the vertebrate ancestor**. That is consistent with the two proteins being
> matched to each other, but it is a single event, not repeated independent evidence.

## 8. Family tree (`src/bin/njtree.rs`, `data/phylo/`)

167 sequences:
- 73 canonical importins-alpha from animals, plants, fungi and *Dictyostelium*, including
  importins from the very molluscs, cnidarians, *Capsaspora* and solanaceous plants that
  carry PF30746;
- 27 vertebrate TMCO6;
- 67 PF30746 members outside vertebrates.

Method: pairwise Smith-Waterman alignments, Kimura distance, neighbor-joining, 100 bootstrap
replicates (columns resampled within each pairwise alignment).

| test | result |
|---|---|
| vertebrate TMCO6 form a single clade | yes, bootstrap **100%** |
| nearest neighbour of any PF30746 member | **another PF30746 from a different species: 66 of 67**; an importin: 1 |
| nearest non-self-group neighbour of vertebrate TMCO6 | PF30746: 26 of 27; an importin: 1 |
| PF30746 outside vertebrates form a single clade | no (0%) |
| TMCO6 + PF30746 against all importins | no (0%) |

Mean Kimura distances:
- importins within one kingdom: 0.46-0.96;
- PF30746 to importins: 1.3-1.8;
- vertebrate TMCO6 to animal PF30746 1.65, to plant PF30746 1.61, to importins 1.86-1.96.

**Conclusion.** PF30746 is a **distinct ancient family of importin-like proteins**, not
mislabelled importins: every member is closer to relatives in other species than to the
importins of its own species. Vertebrate TMCO6 is a fast-evolving member on a long branch,
which is why neither the best-hit test nor NJ can place it reliably. Nearest neighbours and
mean distances favour affinity with PF30746. A maximum-likelihood tree from a proper
multiple alignment would settle it.

**Evolutionary picture:**
1. The ancestor of animals and plants already had an importin-like PF30746 protein with the
   classical W/N ladder.
2. In vertebrates it rebuilt its major site (R130/E134/Y171/Q217).
3. At the same time ZNF511 acquired an R/F-rich N-terminal motif, which in the model lies on
   exactly that site.
4. The specialised TMCO6-ZNF511 pair is therefore a vertebrate innovation.
5. PF30746 proteins of invertebrates and plants most likely do not bind ZNF511 this way.

## 9. AlphaFold 3: the pair on the LSm ring (third pass)

Six jobs on AlphaFold Server (`af3/alphafold_server_jobs.json`), parsed with
`src/bin/af3_report.rs`; full output in `af3/report.txt`. Top-ranked model throughout, no
atom clashes.

### 9.1 Quality control: the model reproduces known biology

This is the main reason to trust the remaining numbers.

- **The ring assembled in the canonical order** LSm2-LSm3-LSm6-LSm5-LSm7-LSm4-LSm8-LSm2
  (ipTM 0.77-0.85 for adjacent pairs). The order was not specified; the model derived it.
- **The ring binds RNA as it should**: each subunit grips the 3' oligo(U) of U6 (u102-u107)
  through the conserved loop 3 and loop 5 residues — LSM2 N37/R63/G64, LSM3 N51/R88/G89,
  LSM8 N35/R63/G64 and so on, pLDDT 92-96. This is the established mechanism of Sm/LSm
  rings.

### 9.2 The pair is confirmed independently

| job | pair ipTM | PAE_min | interface TMCO6 / ZNF511 | of which ZNF511 N-term (1-74) |
|---|---|---|---|---|
| 4, wild type | **0.84** | 1.1 | 133 / 90 res. | 61 |
| 5, TMCO6 E134R+Y171W | 0.83 | 1.0 | 127 / 88 | 61 |
| 6, ZNF511 R42A+F50A+R57A | 0.73 | 1.2 | 126 / 91 | 62 |
| 2, within the LSm2-8 ring | 0.75 | 1.8 | 133 / 87 | 62 |

The interface is the same across all six runs and matches the AlphaFold DB model
(section 7), produced by a different group with a different pipeline. Two thirds of the
ZNF511 contacts come from its N-terminus 1-74, lying in the TMCO6 groove. This is the most
confident result in the whole analysis.

### 9.3 The mutation prediction failed

Substitutions in TMCO6 (E134R, Y171W) **changed nothing**: ipTM 0.83 versus 0.84, same
interface. The triple substitution in ZNF511 gave a moderate drop (0.73 versus 0.84), but
the interface survived there too.

Two readings, and these data cannot distinguish them:
1. The interface is long (133 residues along the whole groove), so 2-3 point substitutions
   should not break it. In that case the bench test needs a deletion of the ZNF511
   N-terminus (1-74) or replacement of the entire cluster, not point mutants.
2. AlphaFold is poor at point substitutions in general: it relies on the evolutionary
   profile of the family rather than on the physics of individual side chains. This is a
   known limitation.

Either way: **the proposed point mutants are not a usable experiment**, and the plan now
calls for an N-terminal deletion instead.

### 9.4 Where the pair sits on the ring

| pair | ipTM | contact size | contact pLDDT |
|---|---|---|---|
| ZNF511 — LSM2 | 0.78-0.80 | 19 / 15 res. | 84-87 |
| ZNF511 — LSM6 | 0.70-0.73 | 12 / 20 | 79-90 |
| ZNF511 — LSM4, LSM5, LSM7, LSM8 | 0.59-0.71 | 5-16 res. | 61-89 |
| **TMCO6 — any subunit** | 0.61-0.80 | **1-4 res.** | 20-85 |

The picture is unambiguous: **the ring holds ZNF511, and TMCO6 hangs off ZNF511** without
touching the ring itself. The high ipTM of the TMCO6-LSM8 pair (0.80) alongside a contact of
only 3 residues is not an interface: it reflects TMCO6 being rigidly tethered to ZNF511 and
therefore confidently positioned relative to the whole complex. Contact size is the more
informative number here.

On the ZNF511 side, the ring is held by its zinc fingers (V108-R140), while the 3' end of U6
(c98-a102, u105-u107) is contacted by the basic C-terminal tail (R240, K243, K246, K248) and
the fingers (K113, R114, R176, K179, K182). ZNF511 therefore lies on the same face of the
ring that the RNA approaches.

⚠ The RNA contacts are low-confidence (pLDDT 55-64), and placing zinc fingers on nucleic
acid is something AlphaFold does readily. Treat this as a hypothesis, not a result.

### 9.5 Nuclear versus cytoplasmic ring

| | LSm2-8 (nuclear, U6) | LSm1-7 (cytoplasmic, mRNA decay) |
|---|---|---|
| ipTM of the whole complex | **0.76** | 0.73 |
| TMCO6 to the ring-specific subunit | LSM8 **0.76-0.80** | LSM1 **0.62** |
| ZNF511 to the ring-specific subunit | LSM8 0.68-0.71 | LSM1 0.58 |
| ZNF511 to LSM2 (shared) | 0.78 | 0.74 |

The preference for the nuclear ring is real but weak. What makes it interesting is where the
difference sits: the ring-specific subunits separate by 0.14-0.18 (LSM8 versus LSM1), while
the shared subunit LSM2 differs by only 0.04. This agrees with DepMap (LSM8 correlates,
LSM1 does not) and with the plant LSm8 purification, but on its own it is not decisive.

### 9.6 Bottom line

Confirmed: the TMCO6-ZNF511 pair is real and built as section 7 described. Refined: it is
ZNF511, not TMCO6, that engages the ring. Refuted: the proposed point mutants are not a
usable test. Unresolved: whether the pair binds the ring at all — ipTM 0.6-0.8 with contacts
of a few residues is too weak, and an experiment is needed.

## Hypothesis and how to test it

**TMCO6-ZNF511 is an importin-like chaperone for the nuclear LSm2-8 / U6 snRNP.** Possible
roles: nuclear import of the ring, or assembly and recycling of U6 together with SART3 and
USB1.

### Computational, done or remaining

- ~~Which site of the ARM solenoid is occupied?~~ Done (section 7): the groove is occupied
  along its length by the ZNF511 N-terminus.
- ~~Conservation of R130/E134/N175/W362/N366~~ — done (section 7).
- ~~AlphaFold 3 with the LSm2-8 ring and U6~~ — done (section 9): ZNF511 engages the ring,
  TMCO6 hangs off it.
- ~~Phylogeny of PF30746 and importins~~ — NJ done (section 8); a maximum-likelihood tree
  remains.
- ~~Point mutants E134R/Y171W and R42A/F50A/R57A~~ — tested (section 9.3), they do not work.
  Model an N-terminal deletion of ZNF511 (1-74) and a truncation to the fingers instead.

### At the bench

- **The key test of the pair** (replacing the point mutants, section 9.3): ZNF511 lacking its
  N-terminus (delta 1-74) should fail to co-precipitate TMCO6, while the full-length protein
  does.
- If the pair really sits on the ring, ZNF511 should be the one that brings it there: a
  TMCO6 pull-down should recover LSm proteins only when ZNF511 is present, and not the
  reverse.
- Knock out TMCO6 or ZNF511 in a line from `dm_toplines.tsv`, then assay LSm8
  nuclear/cytoplasmic distribution, U6 snRNA levels and its 3' end.
- Co-IP against LSM8 versus LSM1.
- Tissue from carrier sheep: splicing in motor neurons.

### Database corrections

- UniProt Q96DC7: remove "Multi-pass membrane", add the importin-alpha-like fold.
- The gene name itself is misleading: "transmembrane and coiled-coil domains 6".

## Reproducing

```bash
cargo build --release
S=<directory holding gene_effect.csv>   # DepMap 26Q1, figshare file 67214582

./target/release/pair      $S/gene_effect.csv data/tmco6/ach_disease.tsv data/tmco6/dm
./target/release/nullpc    $S/gene_effect.csv 10 TMCO6 ZNF511 LSM8 USB1 LARP7 SNUPN
./target/release/plddt     data/tmco6/AF-TMCO6.pdb 1-60 338-358 386-406
./target/release/interface data/tmco6/AF-complex-TMCO6-ZNF511.pdb 130 134 171 175 213 217 362 366 407
./target/release/map_sites data/tmco6/aln_P52293.json 142 146 184 188 231 235 357 361 399 403
./target/release/af3_report af3/results af3/inputs.fa

cd data/tmco6
CONTROL=KPNA2_human.fa ../../target/release/conserve \
  orth_TMCO6.fa fam_TMCO6.tsv Q96DC7 AF-TMCO6.pdb cons_TMCO6.tsv 130 134 171 175 213 217 362 366
../../target/release/rbh orth_TMCO6.fa human_importins.fa all_acc.txt cons_TMCO6.tsv.species
../../target/release/cooccur refprot_euk.tsv ref_TMCO6.tsv ref_ZNF511.tsv 90
cd ../..
```

Input files that are not stored in the repository are listed with download commands in
[data/README.md](data/README.md).

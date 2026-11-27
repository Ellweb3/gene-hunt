# TMCO6 is not a membrane protein but an importin-alpha-like solenoid that binds ZNF511, a candidate module at the nuclear LSm2-8 ring

**Eldar Sufiyanov**

Independent researcher, Punta del Este, Uruguay.

Correspondence: fluttereldar@gmail.com  ·  Code and data: https://github.com/Ellweb3/gene-hunt

---

## Abstract

TMCO6 ("transmembrane and coiled-coil domains 6") is annotated in UniProt as a multi-pass
membrane protein, and that annotation has been used as the premise of published work. We
show it is incorrect. Structural comparison assigns TMCO6 to the importin-alpha family over
its entire length (E ~1e-12 at ~19% sequence identity), placing the two annotated
transmembrane helices among the ordinary helices of an armadillo-repeat solenoid. The entry's
own InterPro, Gene3D and SUPFAM cross-references already describe this fold, including an
importin-beta-binding (IBB) domain signature, so the membrane assignment — which rests only
on sequence-based prediction — is internally inconsistent with the rest of the record. We
also report two measures that do *not* distinguish buried solenoid helices from genuine
transmembrane helices, having tested them against experimentally supported membrane proteins.

We then ask what this solenoid binds. TMCO6 and ZNF511, a zinc-finger protein of unknown
function, are reciprocal top co-dependencies across 1208 cancer cell lines in DepMap
(r = 0.60; next-best partner 0.24 for each), a correlation that survives subtraction of
screen-wide noise estimated from non-expressed olfactory receptor genes. They
co-precipitate with each other and with LSm3-7 in two independent BioPlex cell lines. Two
independently computed structure predictions place the disordered N-terminus of ZNF511
(residues 1-74) in the NLS-binding groove of TMCO6, forming a 133/90-residue interface. That
groove is rebuilt in vertebrates: the canonical W/N ladder of importin-alpha is replaced by
R130/E134/Y171/Q217, a substitution fixed from lamprey to human, while the minor site
(W362/N366) is preserved.

Several weak lines of evidence point to the nuclear LSm2-8 ring rather than the cytoplasmic
LSm1-7 ring: the pair correlates with LSM8, SNUPN, USB1 and LARP7 but not with LSM1; in
*Arabidopsis*, the ZNF511 homolog was recovered in AtLSM8 but not AtLSM1a purifications
(Golisz et al. 2013, though its author tells us that interaction was never validated); and
in a model of the full ring, the pair contacts the ring-specific subunit LSM8 more
confidently than LSM1. We state plainly that this last step
is not established: the predicted pair-to-ring interfaces are small, and a prediction that
point substitutions would disrupt the pair was tested and failed. We propose specific
experiments, including the one we consider decisive, and we have deposited all code and
intermediate results.

**Keywords:** TMCO6, ZNF511, importin-alpha, armadillo repeat, LSm2-8, U6 snRNP,
protein annotation, structure prediction

---

## 1. Introduction

A large fraction of the human proteome remains functionally uncharacterized. The Illuminating
the Druggable Genome consortium classifies about 5500 human proteins as "Tdark": proteins
with almost no dedicated literature and no assigned function. For many of these, the limiting
factor is no longer data. Genome-wide CRISPR dependency screens, systematic affinity-
purification interactomes, proteome-wide structure prediction and structure-based search have
accumulated enough public evidence that specific, testable hypotheses can sometimes be
formed before any experiment is performed.

Such hypotheses are only useful if they are honest about their own strength. Much of what can
be computed about an uncharacterized protein is correlational, and structure prediction
produces confident-looking models of interactions that may not occur. We therefore separate
throughout what we consider established from what we consider suggestive, and we report a
prediction of ours that failed.

This work began as a systematic survey of Tdark proteins for which several independent public
data types converge. One candidate pair, TMCO6 and ZNF511, emerged with unusually consistent
signals, and examining it revealed that a basic structural annotation of TMCO6 is wrong.

TMCO6 (UniProt Q96DC7, 493 aa) is named for predicted transmembrane and coiled-coil domains
and is annotated as a multi-pass membrane protein. It has essentially no functional
literature. A 2024 study reported TMCO6 as a cell-surface receptor for neutrophil
extracellular trap DNA on CD8+ T cells. A homozygous loss-of-function variant in ovine TMCO6
causes recessive early-onset motor neuron degeneration.

ZNF511 (UniProt Q8NB15, 252 aa) contains three C2H2 zinc fingers and is annotated only as
"may be involved in transcriptional regulation" by similarity. Its *Arabidopsis* homolog
At5g14140 appeared in a 2013 survey of LSm-associated proteins.

---

## 2. Results

### 2.1 TMCO6 is an importin-alpha-like solenoid, not a membrane protein

A structural search with the AlphaFold model of TMCO6 (AF-Q96DC7-F1, global pLDDT 87.3)
against AFDB/Swiss-Prot and PDB100 returns importin-alpha proteins from animals, plants and
fungi as the highest-scoring hits outside the TMCO6 family itself. The match spans residues
1-486 at E ~1e-12 despite only ~19% sequence identity, which is why it is invisible to
sequence-based search. PDB hits include importin-alpha in complex with NLS peptides.

Within that fold, the two annotated transmembrane segments (338-358, 386-406) are two
ordinary helices of the repeat stack. Both features carry only computational evidence
(ECO:0000255), and the subcellular location is inferred from them (ECO:0000305); no
experimental support for a membrane localisation is cited in the entry.

**A negative control worth reporting.** We first attempted to support this reassignment with
properties of the segments themselves, and both attempts failed. Burial does not
discriminate: the mean number of Calpha neighbours within 10 A is 18.3 and 20.0 for the two
TMCO6 segments, but experimentally supported transmembrane helices of AQP1 (20.8-23.6),
GLUT1 (20.3-25.6) and ABCG2 (16.0-21.2) are equally or more buried in their own AlphaFold
models, because in a multi-pass membrane protein the transmembrane helices pack against one
another. The hydrophobicity of the solvent-exposed face is better but still unreliable: the
exposed residues of the annotated segments are more hydrophobic than the rest of the chain by
+3.8 Kyte-Doolittle units in GLUT1 and +4.4 in ABCG2, as expected for a lipid-facing surface,
but by −0.3 in AQP1, which is almost entirely membrane-embedded so that any relative exposure
threshold misfires. The TMCO6 segments give −0.4, indistinguishable from AQP1.

<!--FIG figures/fig1_burial.png **Figure 1. Burial does not distinguish the annotated TMCO6 segments from genuine transmembrane helices.** Each point is one annotated transmembrane helix, placed at its mean number of Calpha atoms within 10 A in that protein's own AlphaFold model. The TMCO6 segments (red) fall inside the range of AQP1, GLUT1, MDR1 and ABCG2, whose transmembrane helices pack against one another and are therefore just as buried.-->

We therefore do not rest the reassignment on either measure. The evidence is the fold
assignment itself, together with the entry's own structural cross-references, and we report
the two failed discriminators so that others do not repeat them.

The UniProt entry is internally inconsistent on this point. Alongside the membrane
annotation it carries Gene3D 1.25.10.10 (ARM-like), InterPro IPR011989 (ARM-like), IPR016024
(ARM-type fold), IPR000225 (Armadillo), IPR002652 (importin-alpha, importin-beta-binding
domain) and SUPFAM SSF48371 (ARM repeat). Consistent with the IBB signature, residues 1-60
are disordered (pLDDT 46) and strongly basic, with 26 of 60 residues charged and RRRRER and
RKARR clusters — the architecture of an IBB domain, by which importin-alpha engages
importin-beta, rather than a cytoplasmic tail.

### 2.2 The NLS-binding groove is rebuilt in vertebrates

Mapping the NLS-binding residues of mouse importin-alpha2 (P52293) onto TMCO6 through the
structural alignment shows that the major site is altered while the minor site is preserved:

| importin-alpha site | TMCO6 |
|---|---|
| major, ARM2 W142/N146 | R130/E134 |
| ARM3 W184/N188/D192 | Y171/N175/E179 |
| ARM4 W231/N235 | Y213/Q217 |
| ARM5 W273 | W261 |
| minor, ARM7 W357/N361 | W362/N366 |
| ARM8 E396/W399/N403 | M400/T403/N407 |

Across 751 species with a family member, the substituted residues are conserved within
vertebrates (identity 0.69-0.91) while the corresponding positions in non-vertebrate members
of the family retain the ancestral tryptophan/asparagine ladder (W in 59-64 of 65 species,
N in 63 of 65). The minor site is invariant in both groups. A site that had simply decayed
would not be conserved; this pattern is more consistent with the major site having been
repurposed for a different ligand in the vertebrate lineage.

<!--FIG figures/fig2_conservation.png **Figure 2. The major NLS-binding site is rebuilt in vertebrates; the minor site is not.** Per-position identity to the human residue across 751 species carrying a PF30746 family member (blue), and across the 65 non-vertebrate members alone (orange). At the five major-site positions the human residue is absent outside vertebrates, where the ancestral importin W/N ladder is retained.-->

### 2.3 TMCO6 and ZNF511 are reciprocal top co-dependencies

Across 1208 cell lines in DepMap (Chronos, 26Q1), the dependency profiles of TMCO6 and
ZNF511 correlate at r = 0.601. For each gene the next-best partner is far behind (0.24 and
0.23). Neither gene is broadly essential (mean gene effect −0.11; 7-8% of lines below −0.5).

Because co-dependency correlations are vulnerable to screen-wide technical structure, we
constructed a null set from 193 single-copy olfactory receptor genes, which are not expressed
in cancer lines and whose shared profile structure is therefore noise. Removing the first ten
principal components of that null set does not weaken the pair (r = 0.612). The same
correction collapses the null distribution itself (median −0.005, 99.9th percentile 0.20,
maximum 0.34), placing the pair far outside it. After correction, the next-ranked partners of
TMCO6 are DSCC1 (0.20), LARP7 (0.18), SNUPN (0.18) and USB1 (0.16).

Two earlier versions of this analysis produced artifacts that we report because they are
instructive. Genes measured in few cell lines generated spurious high correlations and were
removed by requiring at least 1000 shared lines. Multi-copy gene families produced false
essentiality through multiple cutting and were excluded by family prefix.

### 2.4 Both proteins associate with an LSm ring

In BioPlex 3.0, TMCO6 recovers ZNF511, LSM4, LSM5 and LSM7 in HEK293T, and ZNF511, LSM6 and
LSM7 in HCT116. ZNF511 recovers TMCO6 and LSM3-7 in both lines. LSM4 and LSM7 baits recover
both proteins reciprocally.

The DepMap correlations distinguish the two rings that share subunits LSm2-7. The pair
correlates with LSM8 (0.15/0.17), which is specific to the nuclear LSm2-8 ring, and not at
all with LSM1 (0.00/−0.01), which is specific to the cytoplasmic LSm1-7 ring. The partners
surviving noise correction belong to the same process: SNUPN (snurportin, nuclear import of
snRNPs), USB1 (3'-end processing of U6 snRNA) and LARP7 (U6 modification).

Experimental support of a kind comes from another kingdom. In a FLAG purification of
*Arabidopsis* LSm proteins, At5g14140 — the ZNF511 homolog — was recovered in two of three
AtLSM8 purifications alongside Prp24/SART3 and Brr2, and was absent from the AtLSM1a
purifications.

We initially treated this as the strongest independent evidence here. After submitting this
manuscript we wrote to its senior author, who replied that nothing is known about At5g14140
in *Arabidopsis*, that the interaction "has not been validated", and that its nuclear
localisation there is predicted rather than confirmed (J. Kufel, personal communication,
2026). The co-purification is therefore an unvalidated AP-MS observation. That is
unremarkable for such datasets, but it is not the firm cross-kingdom anchor we took it for,
and with this caveat the assignment to the nuclear ring rests on the DepMap correlations and
on structure prediction alone.

### 2.5 ZNF511 binds the TMCO6 groove

Two structure predictions computed independently — an AlphaFold DB heterodimer
(AF-0000000211739682, ipTM 0.79, ipSAE 0.81, pDockQ 0.72) and our own AlphaFold 3 run
(ipTM 0.84, PAE_min 1.1) — agree on the same arrangement. The interface comprises 133
residues of TMCO6 and 90 of ZNF511. Two thirds of the ZNF511 contribution comes from its
N-terminal 74 residues, which are disordered in the monomer and which pack into the concave
groove the way an NLS peptide packs into importin-alpha. The rebuilt major-site residues
(R130, E134, Y171, Y213, Q217) lie 3.5-4.5 A from ZNF511 and the preserved minor-site
residues 2.8-3.2 A.

The N-terminal motif of ZNF511 that occupies the groove shows a lineage distribution
matching the rebuilt site: it aligns in 87-100% of mammals, turtles, sharks and teleosts, and
in 0-14% of fungi, molluscs, cnidarians and amphioxus. A per-species test finds no
association outside vertebrates (24% of plants with a family member carry the motif versus
37% without, p = 0.25), indicating that the motif and the rebuilt site arose once in the
vertebrate ancestor rather than repeatedly. This is one event, consistent with the two
proteins having been matched to each other, but it is not independent repeated evidence.

<!--FIG figures/fig3_iptm.png **Figure 3. The TMCO6-ZNF511 interface is reproducible; the predicted mutants are not.** ipTM of the TMCO6-ZNF511 chain pair in the top-ranked AlphaFold 3 model of each job. The pair is unchanged by two substitutions in the rebuilt groove, and survives inside the LSm ring in every configuration tested.-->

### 2.6 A model of the pair with the LSm2-8 ring, and its limits

We modelled TMCO6 + ZNF511 + LSm2-8 + U6 snRNA with AlphaFold 3. As a control on whether the
model can be trusted at all, we examined features it was not given: the ring assembled in the
canonical order LSm2-LSm3-LSm6-LSm5-LSm7-LSm4-LSm8-LSm2 (adjacent-pair ipTM 0.77-0.85), and
each subunit engaged the 3' oligo(U) of U6 through the expected loop 3 and loop 5 residues
(LSM2 N37/R63/G64, LSM3 N51/R88/G89, LSM8 N35/R63/G64; pLDDT 92-96). The model therefore
reproduces established biology it was not told about.

Within that model, the ring is engaged by ZNF511, not by TMCO6. ZNF511 contacts LSM2 over 19
residues (pLDDT 84-87) and other subunits over 5-16 residues, while TMCO6 contacts any
subunit over only 1-4 residues and is held via ZNF511. The high pair-wise ipTM of
TMCO6-LSM8 (0.80) alongside a 3-residue contact reflects TMCO6 being rigidly tethered to
ZNF511 and therefore confidently positioned, not an interface; contact size is the more
informative quantity.

A control against the cytoplasmic ring is consistent with, but does not establish, a
preference for the nuclear one. The whole-complex ipTM is 0.76 for LSm2-8 against 0.73 for
LSm1-7, and the separation sits where it should: on the ring-specific subunits (TMCO6 to
LSM8 0.76-0.80 versus to LSM1 0.62) rather than on the shared subunit LSM2 (0.78 versus
0.74).

We do not claim the pair binds the ring. The predicted pair-to-ring interfaces are small, and
the contacts between ZNF511 and U6 are low-confidence (pLDDT 55-64); placing zinc fingers on
nucleic acid is moreover something structure predictors do readily.

### 2.7 A prediction that failed

We predicted that point substitutions in the rebuilt groove would disrupt the pair, and
tested this in silico. TMCO6 E134R/Y171W changed nothing (ipTM 0.83 versus 0.84 for wild
type; interface 127 versus 133 residues). ZNF511 R42A/F50A/R57A gave a moderate reduction
(0.73) with the interface again intact.

Two readings are available and these data do not separate them. The interface may simply be
too extensive — 133 residues along the whole groove — for two or three substitutions to
disrupt. Alternatively, structure predictors are known to respond poorly to point
substitutions, relying on family-level evolutionary profiles rather than individual side-chain
physics. Either way, the practical conclusion is that these point mutants are not a usable
experiment, and a deletion of the ZNF511 N-terminus should be used instead.

### 2.8 The family is ancient; vertebrate TMCO6 is a fast-evolving member

A best-hit test against a panel of all human importins-alpha shows that family members
outside vertebrates are not demonstrated orthologs of TMCO6: they are equidistant from TMCO6
and from KPNA1-8 (28-33% identity, margins of 1-9 bits), whereas vertebrate TMCO6 is
unambiguous (margins 195-760 bits). A neighbor-joining tree of 167 sequences (73 canonical
importins-alpha, 27 vertebrate TMCO6, 67 non-vertebrate family members; Kimura distances,
100 bootstrap replicates) places vertebrate TMCO6 in a single clade with 100% support. For
66 of 67 non-vertebrate family members the nearest neighbour is another family member from a
different species rather than an importin of their own species, and for 26 of 27 vertebrate
TMCO6 the nearest non-self neighbour is a family member.

The family is therefore distinct and ancient rather than a set of mislabelled importins, and
vertebrate TMCO6 is a fast-evolving member on a long branch. Neither the best-hit test nor
neighbor-joining can place that long branch reliably; a maximum-likelihood analysis from a
curated multiple alignment would be needed to settle its position.

---

## 3. Discussion

The structural reassignment of TMCO6 is, in our view, settled: it can be reproduced in under
an hour from public resources, it is consistent with every structural cross-reference already
in the UniProt entry, and the contrary annotation rests entirely on sequence-based prediction
with no experimental support. Its practical consequence is that a protein currently treated
as a membrane receptor is a soluble nuclear-transport-related solenoid. We have submitted a
correction request to UniProt.

The pair TMCO6-ZNF511 is likewise well supported, by three data types that fail in different
ways: genetic co-dependency, affinity purification, and structure prediction. The interface
is the same in two independently produced models and is unchanged whether the pair is modelled
alone or inside the ring.

The association with the LSm2-8 ring is a hypothesis, and a weaker one than we first judged.
We had leant on the recovery of the ZNF511 homolog in an AtLSM8 purification in plants; its
senior author has since told us that this interaction was never validated (section 2.4).
What remains is concordant but circumstantial: the human genetic data distinguishing LSM8
from LSM1, and the modelling. What our modelling adds is the direction of the
interaction — ZNF511 would engage the ring, with TMCO6 carried along — and this is the part
most in need of testing.

If the module does serve the U6/LSm2-8 snRNP, the natural functions to consider are nuclear
import of the ring, given the importin-alpha architecture and the IBB-like N-terminus of
TMCO6 and the co-dependency with snurportin, or participation in U6 assembly and recycling
alongside SART3 and USB1. We note, without putting weight on it, that recessive loss of TMCO6
in sheep causes early-onset motor neuron degeneration, a class of phenotype also produced by
defects in snRNP assembly and snRNA processing in humans.

### 3.1 Limitations

- No experiments were performed. Every result here is computational.
- Whether the module binds an LSm ring is not established (section 2.6).
- The contacts between ZNF511 and U6 snRNA are low-confidence and should not be relied on.
- One of our own predictions was tested and failed (section 2.7).
- The position of vertebrate TMCO6 within its family is not resolved by the tree we computed
  (section 2.8).
- Human evidence for the LSm2-8 over LSm1-7 preference is correlational; the direct
  experimental evidence is from *Arabidopsis*.
- The co-dependency analysis is restricted to cancer cell lines and may not reflect normal
  physiology.
- The plant co-purification we relied on is an unvalidated AP-MS observation, per its own
  author (section 2.4).

### 3.2 Proposed experiments

We consider the following decisive, in order of informativeness:

1. ZNF511 lacking residues 1-74 should fail to co-precipitate TMCO6, while the full-length
   protein does. This tests the one claim we consider firm.
2. A TMCO6 pull-down should recover LSm proteins only when ZNF511 is present, and depleting
   ZNF511 should abolish it; the reverse should not hold. This tests the predicted direction.
3. Co-immunoprecipitation against LSM8 versus LSM1, to test ring specificity directly.
4. Loss of either protein should affect U6 snRNA levels, its 3' end, or LSm8 nuclear
   localisation, rather than mRNA decay.

---

## 4. Methods

All analysis code is written in Rust and is available in the repository cited below. Each tool
is a standalone binary; command lines reproducing every figure and table are given in the
repository README.

**Dependency correlation.** DepMap Chronos gene effect (26Q1 public release) was used.
Pearson correlations were computed over cell lines in which both genes were measured,
requiring at least 1000 shared lines. For the artifact control, a null set was defined as
genes whose symbol begins with OR followed by a digit and whose mean gene effect lies within
±0.15 of zero (193 genes). Profiles were z-scored with missing values set to the gene mean;
principal components of the null set were obtained by power iteration with deflation in cell-
line space, and partial correlations computed after projecting them out.

**Structure comparison.** AlphaFold models were obtained from the AlphaFold Database.
Structural search used the Foldseek web server (3Di+AA mode) against PDB100 and
AFDB/Swiss-Prot. Burial was quantified as the number of Calpha atoms within 10 A; secondary
structure was assigned from Calpha geometry (d(i,i+3) and d(i,i+4) criteria); hydrophobicity
used the Kyte-Doolittle scale over a 19-residue window. Functional positions were transferred
between proteins through the Foldseek alignment.

**Complex prediction.** Six jobs were run on AlphaFold Server (AlphaFold 3) with seed 1:
the pair alone; the pair with LSm2-8; the pair with LSm2-8 and U6 snRNA (RNU6-1,
ENSG00000206625, 107 nt); the pair with LSm1-7 as a ring control; and two mutant pairs. The
top-ranked model of each job was analysed. Interfaces were defined as residues with any heavy
atom within 4.5 A of the other chain. Chains were assigned to proteins by sequence matching
rather than by input order.

**Conservation and phylogeny.** Family members were retrieved from UniProt by Pfam PF30746
and InterPro IPR039258. One sequence per species was kept, selected by Smith-Waterman score,
requiring at least 60% coverage of the human query. Co-occurrence was assessed across 2087
eukaryotic reference proteomes with BUSCO completeness at least 90%. The neighbor-joining
tree used Kimura-corrected protein distances from pairwise Smith-Waterman alignments, with
100 bootstrap replicates resampling alignment columns within each pair; named clades were
tested as explicit bipartitions. Motif-presence association was tested with Fisher's exact
test per taxonomic group.

**Data sources.** DepMap (Broad Institute), BioPlex 3.0, Human Protein Atlas, UniProt,
InterPro, Pfam, Pharos/IDG, NCBI Gene, Cellosaurus, Ensembl, Foldseek, AlphaFold Database and
Europe PMC.

---

## 5. Data and code availability

All code, intermediate results, structure models and the commands to re-download every input
are at **https://github.com/Ellweb3/gene-hunt**. Models produced on AlphaFold Server are
redistributed there under the AlphaFold Server Output Terms of Use (non-commercial use only).

## 6. Competing interests

None.

## 7. Acknowledgements

This analysis used public data produced by many groups; we particularly note that the single
most informative piece of evidence is an experimental observation published by others.

---

## References

1. Golisz A, Sikorski PJ, Kruszka K, Kufel J. Arabidopsis thaliana LSM proteins function in mRNA splicing and degradation. *Nucleic acids research* 2013;41:6232-6249. doi:10.1093/nar/gkt296. PMID:23620288.
2. Montemayor EJ, Didychuk AL, Yake AD, Sidhu GK, Brow DA, Butcher SE. Architecture of the U6 snRNP reveals specific recognition of 3'-end processed U6 snRNA. *Nature communications* 2018;9:1749. doi:10.1038/s41467-018-04145-4. PMID:29717126.
3. Montemayor EJ, Virta JM, Hayes SM, Nomura Y, Brow DA, Butcher SE. Molecular basis for the distinct cellular functions of the Lsm1-7 and Lsm2-8 complexes. *RNA (New York, N.Y.)* 2020;26:1400-1413. doi:10.1261/rna.075879.120. PMID:32518066.
4. Huttlin EL, Bruckner RJ, Navarrete-Perea J, Cannon JR, Baltier K, Gebreab F, et al. Dual proteome-scale networks reveal cell-specific remodeling of the human interactome. *Cell* 2021;184:3022-3040.e28. doi:10.1016/j.cell.2021.04.011. PMID:33961781.
5. Dempster JM, Boyle I, Vazquez F, Root DE, Boehm JS, Hahn WC, et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome biology* 2021;22:343. doi:10.1186/s13059-021-02540-7. PMID:34930405.
6. Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 2021;596:583-589. doi:10.1038/s41586-021-03819-2. PMID:34265844.
7. Abramson J, Adler J, Dunger J, Evans R, Green T, Pritzel A, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 2024;630:493-500. doi:10.1038/s41586-024-07487-w. PMID:38718835.
8. van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, et al. Fast and accurate protein structure search with Foldseek. *Nature biotechnology* 2024;42:243-246. doi:10.1038/s41587-023-01773-0. PMID:37156916.
9. UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic acids research* 2025;53:D609-D617. doi:10.1093/nar/gkae1010. PMID:39552041.
10. Sheils T, Mathias SL, Siramshetty VB, Bocci G, Bologa CG, Yang JJ, et al. How to Illuminate the Druggable Genome Using Pharos. *Current protocols in bioinformatics* 2020;69:e92. doi:10.1002/cpbi.92. PMID:31898878.
11. Saitou N, Nei M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular biology and evolution* 1987;4:406-425. doi:10.1093/oxfordjournals.molbev.a040454. PMID:3447015.
12. Letko A, Brülisauer F, Häfliger IM, Corr E, Scholes S, Drögemüller C. Loss-of-function variant in the ovine TMCO6 gene in North Country Cheviot sheep with motor neuron disease. *Genomics* 2023;115:110689. doi:10.1016/j.ygeno.2023.110689. PMID:37488055.
13. Song M, Zhang C, Cheng S, Ouyang D, Ping Y, Yang J, et al. DNA of Neutrophil Extracellular Traps Binds TMCO6 to Impair CD8+ T-cell Immunity in Hepatocellular Carcinoma. *Cancer research* 2024;84:1613-1629. doi:10.1158/0008-5472.can-23-2986. PMID:38381538.

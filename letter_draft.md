# Draft letter to a lab

One letter, two possible recipients; only the second paragraph differs.

| recipient | why them | reference |
|---|---|---|
| **Samuel Butcher / Eric Montemayor** (UW–Madison) | structures of the U6 snRNP and of the Lsm2-8 ring bound to processed U6 (Nat Commun 2018; PDB 6ASO, 6PPP); they work on U6 assembly mechanism | [paper](https://www.nature.com/articles/s41467-018-04145-4) |
| **Joanna Kufel** (IGIB, University of Warsaw) | their own AtLSM8 purification contains the ZNF511 homolog, which is the strongest piece of evidence here; they work in both yeast and plants | [lab](https://kufel.igib.uw.edu.pl/lsm-complexes/) |

Before sending: fill in the name and return address in place of `[…]`, check the lab's
current contact details on their site, and link this repository.

---

## Subject

`An uncharacterized candidate partner of the nuclear Lsm2-8 ring: TMCO6-ZNF511`

## Body

Dear Professor [Name],

I am writing with a computational observation about two uncharacterized human proteins
that appear to form a module associated with the nuclear Lsm2-8 ring. I have no wet-lab
capacity, so I am sharing it in case it is useful to your group. All data and code are
open; I would be glad to be told if this is already known or if the reasoning is flawed.

**[variant for Butcher / Montemayor]**
Your structures of the U6 snRNP and of Lsm2-8 bound to processed U6 are what let me
interpret the prediction below, and the question it raises is a structural one.

**[variant for Kufel]**
The strongest independent evidence comes from your own data: in Golisz et al. 2013
(NAR, PMID 23620288), At5g14140 was recovered in two of three AtLSM8 purifications,
alongside Prp24/SART3 and Brr2, and was absent from the AtLSM1a purifications. At5g14140
is the Arabidopsis homolog of human ZNF511.

**The observation.** TMCO6 and ZNF511 are reciprocal top co-dependencies in DepMap
(r = 0.60 across 1208 lines; the next-best partner of each is 0.24). This survives
correction for screen-wide noise, estimated from 193 non-expressed olfactory receptor
genes. In BioPlex 3.0 they co-precipitate with each other and with LSm3-7 in two
independent cell lines, and LSM4/LSM7 baits recover both. In human DepMap the pair
correlates with LSM8, SNUPN, USB1 and LARP7, but not with LSM1.

**Structure.** Despite its name and its UniProt annotation as a multi-pass membrane
protein, TMCO6 is not a membrane protein: Foldseek matches it to importin-alpha along
its entire length (E ~1e-12 against importins from fungi, plants and animals), and both
annotated "transmembrane" helices are buried in the hydrophobic core of the solenoid
(pLDDT 94-97, 0-1 charged residues, ~19 Calpha neighbours within 10 A). Its NLS-binding
groove is rebuilt in vertebrates: the canonical W/N ladder of importin-alpha is replaced
by R130/E134/Y171/Q217, and this substitution is fixed from lamprey to human.

Two independent predictions - the AlphaFold DB heterodimer AF-0000000211739682
(ipTM 0.79) and my own AlphaFold 3 run (ipTM 0.84) - place the disordered N-terminus of
ZNF511 (residues 1-74) in that rebuilt groove, with a 133/90-residue interface.

In an AlphaFold 3 model of TMCO6 + ZNF511 + Lsm2-8 + U6 snRNA, the ring assembles in the
canonical order and binds the 3' oligo(U) through the expected loop 3/5 residues, which
I took as a quality control. In that model ZNF511 contacts the ring (19 residues with
LSm2, pLDDT 84-87) while TMCO6 contacts it only marginally (1-4 residues) and is held
via ZNF511.

**What I am not claiming.** Whether this module binds the ring at all is not established:
the predicted pair-to-ring interfaces are small, and the contacts between ZNF511 and U6
are low-confidence (pLDDT 55-64). Placing zinc fingers on RNA is also something AlphaFold
does readily. The TMCO6-ZNF511 interface itself is the only confident part.

**Testable predictions**, in case any are cheap for you:

1. A TMCO6 pull-down should recover LSm proteins only when ZNF511 is present; depleting
   ZNF511 should abolish it. The reverse should not hold.
2. ZNF511 lacking residues 1-74 should fail to co-precipitate TMCO6 while remaining
   competent to bind the ring. (I tested point mutants in silico - E134R/Y171W in TMCO6,
   R42A/F50A/R57A in ZNF511 - and they did not reduce the predicted ipTM, so a deletion
   rather than point substitutions seems the right experiment.)
3. Loss of either protein should affect U6 snRNA levels, its 3' end, or LSm8 nuclear
   localisation, rather than mRNA decay.

**Possible relevance.** A homozygous loss-of-function variant in ovine TMCO6 causes
recessive early-onset motor neuron degeneration (Genomics 2023, PMID 37488055). If this
module does serve U6/Lsm2-8, that phenotype would sit alongside the other snRNP-assembly
diseases, though this is speculation on my part. Note also that a 2024 study
(Cancer Res, PMID 38381538) treats TMCO6 as a membrane receptor for NET DNA on CD8+
T cells; the DNA-binding N-terminus is consistent with its basic IBB-like character, but
the membrane interpretation is not consistent with the fold.

Data, code and the full write-up: https://github.com/Ellweb3/gene-hunt

With thanks for your time,
[Name]
[Return address]

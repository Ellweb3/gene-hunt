# Correction request for UniProt entry Q96DC7 (TMCO6_HUMAN)

Submit via the "Entry feedback" link on https://www.uniprot.org/uniprotkb/Q96DC7/entry
(it pre-fills the accession), or via https://www.uniprot.org/contact

Entry state at the time of writing: version 171, last annotation update 2026-11-04.

---

**Subject:** Q96DC7 (TMCO6_HUMAN) — transmembrane annotation appears incorrect; the protein
is an importin-alpha-like ARM solenoid

Dear UniProt curators,

I would like to flag what looks like an incorrect annotation in Q96DC7 (TMCO6_HUMAN), and
to point out that the entry is currently internally inconsistent.

**The entry simultaneously asserts two incompatible things.**

It is annotated as a membrane protein:

- TRANSMEM 338..358 "Helical" (ECO:0000255)
- TRANSMEM 386..406 "Helical" (ECO:0000255)
- SUBCELLULAR LOCATION: Membrane (ECO:0000305)
- Keywords: Membrane, Transmembrane, Transmembrane helix

Yet every structural signature already attached to the same entry describes a soluble
alpha-solenoid of the importin-alpha type:

- Gene3D 1.25.10.10 (ARM-like / alpha-alpha superhelix)
- InterPro IPR011989 (ARM-like), IPR016024 (ARM-type fold), IPR000225 (Armadillo)
- InterPro **IPR002652 (Importin-alpha, importin-beta-binding domain)**
- SUPFAM SSF48371 (ARM repeat)
- Pfam PF30746 / InterPro IPR064072 (TMCO6, ARM-like domain)

The membrane assignment rests only on a sequence-based prediction (ECO:0000255) and on a
curator inference from it (ECO:0000305). No experimental evidence is cited.

**Structural evidence that the two helices are buried, not membrane-spanning.**

In the AlphaFold model AF-Q96DC7-F1 (global pLDDT 87.3) both annotated segments sit in the
hydrophobic core of the solenoid rather than on its surface:

| segment | mean pLDDT | Calpha neighbours within 10 A | charged residues |
|---|---|---|---|
| 338-358 | 94 | 18.3 | 1 |
| 386-406 | 97 | 20.0 | 0 |

A genuine transmembrane helix in a monomeric protein is solvent-exposed in such a model.
These two are packed against the rest of the repeat stack, which is the expected failure
mode of sequence-based TM predictors on buried hydrophobic helices of alpha-solenoids.

**Fold assignment.** A Foldseek search with AF-Q96DC7-F1 against AFDB/Swiss-Prot and PDB100
returns importin-alpha proteins from animals, plants and fungi as the top hits outside the
TMCO6 family itself, matching over residues 1-486 at E ~1e-12 despite only ~19% sequence
identity. PDB hits include importin-alpha bound to NLS peptides (for example 1bk6, 6wx8).

Consistent with IPR002652, residues 1-60 are disordered (pLDDT 46) and strongly basic
(26 of 60 residues charged, with RRRRER and RKARR clusters), which is the architecture of an
IBB domain rather than of a cytoplasmic tail.

**Suggested changes.**

1. Remove TRANSMEM 338..358 and TRANSMEM 386..406.
2. Remove SUBCELLULAR LOCATION: Membrane (ECO:0000305) and the keywords Membrane,
   Transmembrane and Transmembrane helix, unless experimental support exists that I have
   not found.
3. Consider recording the importin-alpha-like ARM-repeat architecture, which the entry's own
   InterPro, Gene3D and SUPFAM cross-references already support.

**Why this matters beyond the entry.** The membrane annotation and the gene name are being
used downstream as a premise: a 2024 study (Cancer Res, PMID 38381538) describes TMCO6 as a
cell-surface receptor for neutrophil-extracellular-trap DNA on CD8+ T cells. Its
experimental observations are of course independent of this annotation, but the structural
assignment of TMCO6 to a membrane protein is, I believe, not supportable.

I have no connection to UniProt and no commercial interest; this came out of an open
reanalysis of under-studied human proteins. The supporting computation, including the
Foldseek output and the per-residue burial figures, is at
https://github.com/Ellweb3/gene-hunt (see section 3 of TMCO6_ZNF511.md).

Happy to provide any further detail.

Eldar Sufiyanov

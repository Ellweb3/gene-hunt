# AlphaFold Server jobs

`alphafold_server_jobs.json` holds six jobs. Protein sequences come from UniProt
(`inputs.fa`); the U6 snRNA comes from Ensembl (RNU6-1, ENSG00000206625, 107 nt).

| # | composition | purpose |
|---|---|---|
| 1 | TMCO6 + ZNF511 + LSm2-8 + U6 snRNA | the main question: where the pair sits on the U6 ring |
| 2 | TMCO6 + ZNF511 + LSm2-8 | the same without RNA: does binding depend on it |
| 3 | TMCO6 + ZNF511 + LSm1-7 | control: the cytoplasmic mRNA-decay ring |
| 4 | TMCO6 + ZNF511 | wild type, the baseline for 5 and 6 |
| 5 | TMCO6 E134R/Y171W + ZNF511 | disrupting the rebuilt TMCO6 site |
| 6 | TMCO6 + ZNF511 R42A/F50A/R57A | disrupting the ZNF511 motif |

To run them: sign in at https://alphafoldserver.com, use **Upload JSON**, then submit each
job (**Continue and preview job** followed by **Confirm and submit job**). The daily quota is
30 jobs.

## Results (2026-11-27)

All six jobs completed; archives are extracted into `results/<job name>/`. To parse them:

```bash
cargo build --release --bin af3_report
./target/release/af3_report af3/results af3/inputs.fa > af3/report.txt

# ALL_PAIRS=1 also prints contacts that do not involve TMCO6 or ZNF511,
# which is how the LSm ring geometry and its RNA binding were checked
ALL_PAIRS=1 ./target/release/af3_report af3/results/fold_1_tmco6_znf511_lsm2_8_u6 af3/inputs.fa
```

Outcome, in short: the pair was confirmed (ipTM 0.84, a 133/90-residue interface); the LSm
ring is engaged by ZNF511 rather than TMCO6; and the point-mutant predictions failed. Full
discussion in section 9 of [../TMCO6_ZNF511.md](../TMCO6_ZNF511.md).

> The models here were produced on AlphaFold Server and are distributed under its terms:
> non-commercial use only. The full text sits in `results/*/terms_of_use.md` and at
> [alphafoldserver.com/output-terms](https://alphafoldserver.com/output-terms). Only the
> top-ranked model and its confidence summary are stored per job.

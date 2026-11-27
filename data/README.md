# Obtaining the input data

Raw downloads from third-party databases are not stored in this repository: each carries its
own licence, and together they run to some 60 MB. Below are the commands that produced them.
Computed results (`tdark_profile.tsv`, `tmco6/dm_*.tsv`, `tmco6/cons_*.tsv`, `phylo/`,
`../af3/`) are stored here.

Put the large files (DepMap, gene2pubmed) in a temporary directory; the rest go where the
command points.

## Candidate search

```bash
S=$(mktemp -d)

# DepMap Chronos 26Q1 dependency matrix (413 MB), figshare 31660582
curl -sL https://ndownloader.figshare.com/files/67214582 -o $S/gene_effect.csv

# NCBI: publication counts per gene (292 MB), plus symbols and cytobands
curl -s https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz -o $S/gene2pubmed.gz
curl -s https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz \
  -o data/Homo_sapiens.gene_info.gz

# BioPlex 3.0, two cell lines
curl -s https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv \
  -o data/BioPlex_293T_Network_10K_Dec_2019.tsv
curl -s https://bioplex.hms.harvard.edu/data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv \
  -o data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv

# Human Protein Atlas: subcellular localisation
curl -sL https://www.proteinatlas.org/download/subcellular_location.tsv.zip -o data/hpa.zip
unzip -o data/hpa.zip -d data/

# Pharos (IDG): the Tdark protein list -> data/pharos_tdark.json
# GraphQL, paginated via the nested targets(top, skip) field
```

## TMCO6-ZNF511 analysis (`data/tmco6/`)

```bash
# AlphaFold models (a normal User-Agent is required, otherwise the server returns 403)
curl -s -A "Mozilla/5.0" https://alphafold.ebi.ac.uk/api/prediction/Q96DC7 \
  | python3 -c "import sys,json;print(json.load(sys.stdin)[0]['pdbUrl'])" \
  | xargs curl -s -A "Mozilla/5.0" -o data/tmco6/AF-TMCO6.pdb
curl -s -A "Mozilla/5.0" https://alphafold.ebi.ac.uk/api/prediction/Q8NB15 \
  | python3 -c "import sys,json;print(json.load(sys.stdin)[0]['pdbUrl'])" \
  | xargs curl -s -A "Mozilla/5.0" -o data/tmco6/AF-ZNF511.pdb

# The predicted TMCO6-ZNF511 complex from AlphaFold DB (served gzipped)
curl -s -A "Mozilla/5.0" https://alphafold.ebi.ac.uk/api/complex/Q96DC7 \
  -o data/tmco6/afdb_complex.json
curl -s -A "Mozilla/5.0" https://alphafold.ebi.ac.uk/files/AF-0000000211739682-model_v1.pdb -o c.gz
gunzip -c c.gz > data/tmco6/AF-complex-TMCO6-ZNF511.pdb && rm c.gz

# Foldseek structural search with the TMCO6 model (the response is ~31 MB)
curl -s -F q=@data/tmco6/AF-TMCO6.pdb -F mode=3diaa \
  -F 'database[]=pdb100' -F 'database[]=afdb-swissprot' https://search.foldseek.com/api/ticket
# then, once the ticket status is COMPLETE:
# curl -s https://search.foldseek.com/api/result/<id>/0 -o data/tmco6/foldseek_tmco6.json

# Family members from UniProt: PF30746 (TMCO6) and IPR039258 (ZNF511)
curl -s "https://rest.uniprot.org/uniprotkb/stream?query=xref:pfam-PF30746&fields=accession,organism_id,lineage,length,protein_name&format=tsv" \
  -o data/tmco6/fam_TMCO6.tsv
curl -s "https://rest.uniprot.org/uniprotkb/stream?query=xref:interpro-IPR039258&fields=accession,organism_id,lineage,length,protein_name&format=tsv" \
  -o data/tmco6/fam_ZNF511.tsv
# the same queries with "+AND+keyword:KW-1185" (reference proteome) give ref_TMCO6.tsv and ref_ZNF511.tsv

# Eukaryotic reference proteomes (note: the working filter is reference:true)
curl -s "https://rest.uniprot.org/proteomes/stream?query=reference:true+AND+taxonomy_id:2759&fields=upid,organism_id,lineage,protein_count,busco&format=tsv" \
  -o data/tmco6/refprot_euk.tsv

# Sequences by accession list, in batches of up to 400
# curl -s "https://rest.uniprot.org/uniprotkb/accessions?accessions=<comma-separated>&format=fasta"

# Mapping DepMap cell lines to tumour types
# https://api.cellosaurus.org/search/cell-line?q=dr:DepMap&fields=ac,id,dr,di&format=json&rows=1000&start=<N>

# Full text of Golisz et al. 2013 (the AtLSM8 purification)
curl -s "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC3695525/unicode" \
  -o data/tmco6/golisz.json
```

## Family tree (`data/phylo/`)

```bash
# Reviewed importins-alpha from all kingdoms
curl -s "https://rest.uniprot.org/uniprotkb/stream?query=family:%22importin%20alpha%20family%22+AND+reviewed:true&format=fasta" \
  -o data/phylo/kpna_reviewed.fa

# Plus importins from the species that also carry PF30746
# (organism_id 2607531, 29159, 2652724, 595528, 4081, 283909, 45351),
# queried with protein_name:"importin subunit alpha" and length:[450 TO 620]
```

`phylo_set.fa` is assembled from these files plus `tmco6/orth_TMCO6.fa`. Sequence labels
carry the group tag as a prefix: `TMCO6v`, `PF30746` or `KPNA`.

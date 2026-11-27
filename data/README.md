# Как получить исходные данные

Сырые выгрузки из чужих баз в репозитории не хранятся: у них свои лицензии, и вместе они
весят под 60 МБ. Здесь команды, которыми они были получены. Собственные результаты
(`tdark_profile.tsv`, `tmco6/dm_*.tsv`, `tmco6/cons_*.tsv`, `phylo/`, `../af3/`) лежат
в репозитории.

Большие файлы (DepMap, gene2pubmed) кладите во временную папку, остальные — туда,
куда указывает команда.

## Поиск кандидатов

```bash
S=<временная папка>

# DepMap Chronos 26Q1, матрица зависимостей (413 МБ), figshare 31660582
curl -sL https://ndownloader.figshare.com/files/67214582 -o $S/gene_effect.csv

# NCBI: число публикаций на ген (292 МБ) и символы/цитобэнды
curl -s https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz -o $S/gene2pubmed.gz
curl -s https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz \
  -o data/Homo_sapiens.gene_info.gz

# BioPlex 3.0, две клеточные линии
curl -s https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv \
  -o data/BioPlex_293T_Network_10K_Dec_2019.tsv
curl -s https://bioplex.hms.harvard.edu/data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv \
  -o data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv

# Human Protein Atlas: локализация в клетке
curl -sL https://www.proteinatlas.org/download/subcellular_location.tsv.zip -o data/hpa.zip
unzip -o data/hpa.zip -d data/

# Pharos (IDG): список белков класса Tdark → data/pharos_tdark.json
# GraphQL, постранично по полю targets(top, skip); см. scripts/ и историю коммитов
```

## Разбор TMCO6–ZNF511 (`data/tmco6/`)

```bash
# Модели AlphaFold (нужен обычный User-Agent, иначе 403)
curl -s -A "Mozilla/5.0" https://alphafold.ebi.ac.uk/api/prediction/Q96DC7 \
  | python3 -c "import sys,json;print(json.load(sys.stdin)[0]['pdbUrl'])" \
  | xargs curl -s -A "Mozilla/5.0" -o data/tmco6/AF-TMCO6.pdb
curl -s -A "Mozilla/5.0" https://alphafold.ebi.ac.uk/api/prediction/Q8NB15 \
  | python3 -c "import sys,json;print(json.load(sys.stdin)[0]['pdbUrl'])" \
  | xargs curl -s -A "Mozilla/5.0" -o data/tmco6/AF-ZNF511.pdb

# Предсказанный комплекс TMCO6–ZNF511 из AlphaFold DB (файл отдаётся сжатым)
curl -s -A "Mozilla/5.0" https://alphafold.ebi.ac.uk/api/complex/Q96DC7 -o data/tmco6/afdb_complex.json
curl -s -A "Mozilla/5.0" https://alphafold.ebi.ac.uk/files/AF-0000000211739682-model_v1.pdb -o c.gz
gunzip -c c.gz > data/tmco6/AF-complex-TMCO6-ZNF511.pdb && rm c.gz

# Структурный поиск Foldseek по модели TMCO6 (ответ ~31 МБ)
curl -s -F q=@data/tmco6/AF-TMCO6.pdb -F mode=3diaa \
  -F 'database[]=pdb100' -F 'database[]=afdb-swissprot' https://search.foldseek.com/api/ticket
# затем по полученному id, когда статус станет COMPLETE:
# curl -s https://search.foldseek.com/api/result/<id>/0 -o data/tmco6/foldseek_tmco6.json

# Семейства по UniProt: PF30746 (TMCO6) и IPR039258 (ZNF511)
curl -s "https://rest.uniprot.org/uniprotkb/stream?query=xref:pfam-PF30746&fields=accession,organism_id,lineage,length,protein_name&format=tsv" \
  -o data/tmco6/fam_TMCO6.tsv
curl -s "https://rest.uniprot.org/uniprotkb/stream?query=xref:interpro-IPR039258&fields=accession,organism_id,lineage,length,protein_name&format=tsv" \
  -o data/tmco6/fam_ZNF511.tsv
# то же с "+AND+keyword:KW-1185" (референсный протеом) → ref_TMCO6.tsv, ref_ZNF511.tsv

# Эукариотические референсные протеомы (фильтр именно reference:true)
curl -s "https://rest.uniprot.org/proteomes/stream?query=reference:true+AND+taxonomy_id:2759&fields=upid,organism_id,lineage,protein_count,busco&format=tsv" \
  -o data/tmco6/refprot_euk.tsv

# Последовательности по спискам accession (пачками до 400 штук)
# curl -s "https://rest.uniprot.org/uniprotkb/accessions?accessions=<через запятую>&format=fasta"

# Соответствие линий DepMap и типов опухолей
# https://api.cellosaurus.org/search/cell-line?q=dr:DepMap&fields=ac,id,dr,di&format=json&rows=1000&start=<N>

# Полный текст Golisz et al. 2013 (очистка AtLSM8)
curl -s "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC3695525/unicode" \
  -o data/tmco6/golisz.json
```

## Дерево семейства (`data/phylo/`)

```bash
# Рецензированные импортины-альфа всех царств
curl -s "https://rest.uniprot.org/uniprotkb/stream?query=family:%22importin%20alpha%20family%22+AND+reviewed:true&format=fasta" \
  -o data/phylo/kpna_reviewed.fa
# плюс импортины тех видов, у которых есть PF30746 (organism_id 2607531, 29159, 2652724,
# 595528, 4081, 283909, 45351) — запрос protein_name:"importin subunit alpha" length:[450 TO 620]
```

Набор `phylo_set.fa` собирается из этих файлов плюс `tmco6/orth_TMCO6.fa`; метки групп —
`TMCO6v`, `PF30746`, `KPNA`.

# depclaims

An atlas of published claims that silencing a gene reduces the growth of a named human cancer cell line, each scored against genome-wide
CRISPR and RNA-interference screens and gene expression of the same gene in the same line in the Cancer Dependency Map (DepMap), and
the `depclaims` Python package, which scores any such claim.

This repository accompanies the manuscript *Genome-wide perturbation screens calibrate the landscape of published gene dependencies in
human cancer* (Xu, Cheng, Pan, Chen, Luo, Zhu and Guo; submitted).

## Contents

| Folder | What it holds |
|---|---|
| `atlas/` | The atlas tables: 13,610 scored claims (`claims_scored.tsv.gz`), all 17,920 extracted statements (`statements_all.tsv.gz`) and a summary of the 18,530 CRISPR-screened genes (`genes.tsv.gz`); column definitions in `atlas/README.md` |
| `src/depclaims/`, `tests/` | The depclaims package and its unit tests |
| `analysis/` | The analysis code used for the paper (see below) |

## Install

```bash
pip install .            # Python >= 3.10; pandas, numpy, pyarrow
```

## Build the local cache (once)

Download from the DepMap portal the CRISPR gene effect (Chronos) table, the DEMETER2 combined RNAi gene effect table, the protein-coding
expression table (log2 TPM + 1) and the model annotation table, plus the CEGv2 and NEGv1 gene lists, then:

```bash
depclaims build --chronos CRISPRGeneEffect.csv --demeter2 D2_combined_gene_dep_scores.csv \
  --expression OmicsExpressionProteinCodingGenesTPMLogp1.csv --models Model.csv \
  --ceg CEGv2.txt --neg NEGv1.txt --out depmap_cache
```

## Use

```bash
depclaims score --cache depmap_cache KRAS PANC-1
depclaims lineage --cache depmap_cache ERBB2 Breast
depclaims published --cache depmap_cache --claims claims_scored.tsv.gz EGFR U87MG
```

```python
from depclaims import Atlas
A = Atlas("depmap_cache", claims="claims_scored.tsv.gz")
A.score("KRAS", "PANC-1").as_dict()
A.lineage_profile("ERBB2", "Breast").head()
A.published("EGFR")
```

Cell lines can be given as DepMap model IDs, cell line names or stripped names (case and punctuation are ignored); genes as HGNC
symbols or Entrez IDs.

## Tests

```bash
python -m pytest tests
```

## Licence

MIT. DepMap data are subject to the DepMap terms of use.

## Analysis code

`analysis/` holds the scripts exactly as they were run: candidate selection from the PubMed 2026 baseline and the PubTator3 annotations
(`m5/01`, `m5/02`), claim reading (`m5/05` for the local model, `29_claude_read.py` for the hosted model), mapping and scoring (`03`–`12`),
replication classes (`m5/14`, `15`), figures (`16`, `m5/17`, `20`–`24`), the atlas tables (`25`) and the additional analyses (`28`, `30`–`35`,
`38`). They read local copies of the public inputs from paths set at the top of each script and record provenance with a small internal
wrapper (`rdata.proj`) that can be removed without changing any result. The figure scripts use the authors' plotting style package. The
`depclaims` package is the portable way to reproduce the scored values: it rebuilds every value in `atlas/claims_scored.tsv.gz` from the
public DepMap files.

## Data sources

PubMed 2026 baseline (US National Library of Medicine); PubTator3; Cellosaurus; DepMap Public 26Q1 CRISPR gene effects
(doi:10.6084/m9.figshare.31660582.v2); DEMETER2 combined RNAi scores (doi:10.6084/m9.figshare.9170975.v1); DepMap 24Q4 expression and model
annotations (doi:10.25452/figshare.plus.27993248.v1); Genomics of Drug Sensitivity in Cancer release 8.5; NIH iCite; Europe PMC.

## License

Code: MIT (`LICENSE`). Atlas tables: CC BY 4.0 (`atlas/LICENSE`).

## Citation

See `CITATION.cff`.

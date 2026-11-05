# Atlas of published cancer-cell growth claims scored against DepMap screens

Three tables. In the release folder they are tab-separated and gzipped; in the PNAS submission they are Datasets S1–S3 as CSV files.

## Dataset S1 (claims_scored) — 13,610 claims from 10,420 papers

One row per published claim that perturbing a gene reduces the growth of a named cell line (or increases its apoptosis), anchored to a
verbatim sentence of the abstract and deduplicated per paper, gene and line, for genes and lines screened by CRISPR in DepMap.

| column | meaning |
|---|---|
| pmid, year, journal | the paper |
| symbol, entrez_id | the gene (HGNC symbol, Entrez ID) |
| depmap_model, cell_line, lineage, disease | the named line (DepMap model ID, stripped name, Oncotree lineage and primary disease) |
| perturbation | RNAi (siRNA or shRNA), CRISPR, or other (antisense oligonucleotides, or knockdown, knockout, silencing or depletion without the reagent named) |
| growth_effect, apoptosis_increased | the effect as stated in the abstract |
| gene_class | common essential (dependent lines ≥ 90% of screened lines), selective (> 0%), never (none) |
| crispr_effect_scaled | Chronos gene effect (DepMap Public 26Q1) in the named line, scaled so that common-essential genes have median −1 and non-essential genes median 0 |
| crispr_supported | the named line is dependent on the gene (scaled effect ≤ −0.5) |
| crispr_rank_in_lineage, rnai_rank_in_lineage | rank of the named line among screened lines of its lineage (0 = most dependent, 0.5 = a random line); empty when the lineage has fewer than five screened lines |
| rnai_effect_scaled | DEMETER2 effect in the named line, scaled the same way; empty when the line was not screened by RNAi |
| expression_log2_tpm1, not_expressed | expression of the gene in the named line (DepMap 24Q4, log2 TPM + 1); not_expressed = TPM < 1 |
| gene_dependent_fraction_all, gene_dependent_fraction_lineage | share of all screened lines, and of lines in the named lineage, that are dependent on the gene |
| n_papers_pair, repeat_class | papers making the same gene–line claim; single, independent repeat (at least two papers share no author) or shared-author repeat |
| retracted | the paper is retracted |
| relative_citation_ratio | NIH iCite relative citation ratio |
| split | analysis set used in the paper (development or confirmatory) |

## Dataset S2 (statements_all) — 17,920 gene–line growth statements from 12,825 abstracts

Every statement read from the abstracts that perturbing a gene changes the growth of a cell line, in any direction, with the gene and
line as written, their mapping to Entrez and DepMap where one exists, the perturbation method, the direction of the effect and whether the
statement is anchored to a verbatim sentence of the abstract (grounded).

## Dataset S3 (genes) — 18,530 CRISPR-screened genes

Per gene: share of screened lines that are dependent on it, gene class, number of scored claims, papers and supported claims, and number
of lines named in claims.

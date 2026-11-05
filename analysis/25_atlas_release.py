"""Release tables of the atlas (local package under release/; upload only on the user's decision).
claims_scored: every reduced-growth claim scored in DepMap (both splits), with screen evidence and repeat class.
statements_all: every gene-line growth statement read from the abstracts, mapped or not.
genes: per-gene summary over all CRISPR-screened genes."""
import sys, os, unicodedata, itertools, numpy as np, pandas as pd, duckdb
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
from rdata.proj import Project
P = Project(__file__)
fold = lambda s: ''.join(c for c in unicodedata.normalize('NFKD', s or '') if not unicodedata.combining(c)).lower().replace('-', ' ').strip()
os.makedirs('release', exist_ok=True)
with P.run(__file__, params=dict(author_key='last name + first initial or ORCID')) as run:
    C = pd.concat([pd.read_parquet('db/parquet/claims_development.parquet'), pd.read_parquet('db/parquet/claims_confirmatory.parquet')], ignore_index=True)
    B = pd.read_parquet('work/fig17/out/claims_base.parquet')[['pmid', 'gid', 'ach', 'split', 'symbol', 'frac_all', 'frac_lin']]
    C = C.merge(B, on=['pmid', 'gid', 'ach'], how='left')
    M = pd.read_csv('raw/Model.csv', usecols=['ModelID', 'StrippedCellLineName', 'OncotreePrimaryDisease']).rename(columns={'ModelID': 'ach'})
    C = C.merge(M, on='ach', how='left')
    A = pd.read_parquet('work/dcauth/out/authors.parquet')
    A['key'] = A['last'].map(fold) + '|' + A.fore.map(lambda f: fold(f)[:1])
    au = A.groupby('pmid').apply(lambda d: set(d.key) | {'orcid:' + o for o in d.orcid.dropna()}, include_groups=False).to_dict()
    rows = []
    for (g, a), d in C.groupby(['gid', 'ach']):
        ps = sorted(d.pmid.unique())
        if len(ps) == 1: cls = 'single'
        elif any(p not in au or not au[p] for p in ps): cls = 'undetermined'
        else: cls = 'independent repeat' if any(not (au[x] & au[y]) for x, y in itertools.combinations(ps, 2)) else 'shared-author repeat'
        rows.append((g, a, len(ps), cls))
    C = C.merge(pd.DataFrame(rows, columns=['gid', 'ach', 'n_papers_pair', 'repeat_class']), on=['gid', 'ach'])
    out = C.rename(columns={'gid': 'entrez_id', 'ach': 'depmap_model', 'StrippedCellLineName': 'cell_line', 'OncotreeLineage': 'lineage',
                            'OncotreePrimaryDisease': 'disease', 'yr': 'year', 'meth': 'perturbation', 'gclass': 'gene_class',
                            'effect': 'crispr_effect_scaled', 'sup': 'crispr_supported', 'pct_CRISPR': 'crispr_rank_in_lineage',
                            'v_RNAi': 'rnai_effect_scaled', 'pct_RNAi': 'rnai_rank_in_lineage', 'expr': 'expression_log2_tpm1',
                            'frac_all': 'gene_dependent_fraction_all', 'frac_lin': 'gene_dependent_fraction_lineage'})
    cols = ['pmid', 'year', 'journal', 'symbol', 'entrez_id', 'depmap_model', 'cell_line', 'lineage', 'disease', 'perturbation', 'growth_effect',
            'apoptosis_increased', 'gene_class', 'crispr_effect_scaled', 'crispr_supported', 'crispr_rank_in_lineage', 'rnai_effect_scaled',
            'rnai_rank_in_lineage', 'expression_log2_tpm1', 'not_expressed', 'gene_dependent_fraction_all', 'gene_dependent_fraction_lineage',
            'n_papers_pair', 'repeat_class', 'retracted', 'relative_citation_ratio', 'split']
    out = out[cols].sort_values(['symbol', 'cell_line', 'year', 'pmid'])
    out.to_csv('release/claims_scored.tsv.gz', sep='\t', index=False, float_format='%.4g')
    S = pd.read_parquet('work/reader_mapped_explore.parquet').merge(M, on='ach', how='left')
    meta = duckdb.sql("select pmid, any_value(pub_year) as pub_year from read_parquet('work/cand/cand_*.parquet') group by 1").df()
    S = S.merge(meta, on='pmid', how='left').rename(columns={'pub_year': 'year', 'gene_str': 'gene_as_written', 'line_str': 'line_as_written', 'gid': 'entrez_id',
                                                               'ach': 'depmap_model', 'StrippedCellLineName': 'cell_line'})
    S['entrez_id'] = S.entrez_id.astype('Int64')
    S = S[['pmid', 'year', 'gene_as_written', 'entrez_id', 'line_as_written', 'cell_line', 'depmap_model', 'method', 'growth_effect',
           'apoptosis_increased', 'grounded']].sort_values(['pmid', 'gene_as_written'])
    S.to_csv('release/statements_all.tsv.gz', sep='\t', index=False)
    G = pd.read_parquet('work/fig17/out/genes.parquet').rename(columns={'gid': 'entrez_id', 'frac': 'dependent_fraction', 'gclass': 'gene_class',
                                                                       'n_sup': 'n_claims_supported', 'n_lines': 'n_lines_claimed'})
    G.sort_values('n_papers', ascending=False).to_csv('release/genes.tsv.gz', sep='\t', index=False, float_format='%.4g')
    run.table('release_claims_scored', out, note='atlas: scored reduced-growth claims, both splits, with repeat class')
    print('claims_scored', len(out), out.pmid.nunique(), '| repeat classes', out.repeat_class.value_counts().to_dict())
    print('statements_all', len(S), S.pmid.nunique(), '| genes', len(G))
    print('independent repeat support (all claims)', out.groupby('repeat_class').crispr_supported.mean().round(3).to_dict())

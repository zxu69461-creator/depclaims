"""Reader-model items: one row per paper with >= 1 claim triple whose DepMap cell line is named in the title/abstract.
genes = human DepMap genes with a perturbation-tagged sentence (symbol plus the forms written in the paper); lines = DepMap lines named in the abstract."""
import duckdb, pandas as pd
con = duckdb.connect()
sym = {int(c.split('(')[-1].rstrip(')')): c.split(' (')[0] for c in open('raw/gene_effect_26q1_chronos.csv').readline().strip().split(',')[1:] if '(' in c}
cn = pd.read_csv('raw/cvcl2ach.tsv', sep='\t').drop_duplicates('ach').set_index('ach').cell_name
tn = con.execute("SELECT DISTINCT pmid, gid, ach FROM read_parquet('work/triples_named.parquet') WHERE named").df()
ab = con.execute("SELECT pmid, any_value(title) title, any_value(abstract) abstract FROM read_parquet('work/cand/cand_*.parquet') GROUP BY 1").df()
g = tn.groupby('pmid').gid.apply(lambda s: sorted({sym[int(x)] for x in s})).rename('genes')
l = tn.groupby('pmid').ach.apply(lambda s: sorted({cn[x] for x in s})).rename('lines')
it = ab.merge(g, on='pmid').merge(l, on='pmid')
it.to_parquet('work/llm_items_all.parquet', index=False)
qc = pd.read_csv('validation/claim_qc_v1_key.tsv', sep='\t')
it[it.pmid.isin(qc.pmid)].to_parquet('work/llm_items_qc150.parquet', index=False)
print('items', len(it), '| QC items with a named line', it.pmid.isin(qc.pmid).sum())

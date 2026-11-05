"""Claim triples (pmid, human gene, DepMap model) and the DepMap 26Q1 Chronos effect of that gene in that model.
Scaling per model (charter ruler): effect' = (e - med_NEG) / (med_NEG - med_CEG) so that core essentials (Hart CEGv2) sit at -1 and
non-essentials (Hart NEGv1) at 0. Writes work/triples.parquet and work/depmap_long.parquet (only claimed genes, all models)."""
import numpy as np, pandas as pd, duckdb
con = duckdb.connect()
cl = con.execute("""SELECT pmid, gene_id, list(DISTINCT direction) dirs, list(DISTINCT method) meths, any_value(sentence) sentence, any_value(trigger) trig
                    FROM read_parquet('work/claims/claims_*.parquet', union_by_name=true) GROUP BY 1, 2""").df()
ge = pd.read_csv('raw/gene_effect_26q1_chronos.csv', index_col=0)
ge = ge[[c for c in ge.columns if '(' in c]]  # one unnamed column in the release file is dropped
ge.columns = [int(c.split('(')[-1].rstrip(')')) for c in ge.columns]
ceg = pd.read_csv('raw/CEGv2.txt', sep='\t').ENTREZ_ID.dropna().astype(int); neg = pd.read_csv('raw/NEGv1.txt', sep='\t').ENTREZ_ID.dropna().astype(int)
mc = ge[[g for g in ceg if g in ge.columns]].median(axis=1); mn = ge[[g for g in neg if g in ge.columns]].median(axis=1)
print('raw medians across models: CEG', round(mc.median(), 3), 'NEG', round(mn.median(), 3), '| models', ge.shape[0], 'genes', ge.shape[1])
sc = ge.sub(mn, axis=0).div((mn - mc), axis=0)
cl['gid'] = pd.to_numeric(cl.gene_id, errors='coerce')
genes = sorted(set(cl.gid.dropna().astype(int)) & set(sc.columns))
long = sc[genes].stack().rename('effect').reset_index().rename(columns={'level_0': 'ach', 'index': 'ach', 'level_1': 'gid'})
long.columns = ['ach', 'gid', 'effect']
long['pct'] = long.groupby('ach').effect.rank(pct=True)  # within-model percentile among claimed genes (descriptive only)
long.to_parquet('work/depmap_long.parquet', index=False)
cells = con.execute("SELECT DISTINCT pmid, ach, cvcl FROM read_parquet('work/claims/cell_ann.parquet') WHERE ach IS NOT NULL").df()
nl = cells.groupby('pmid').ach.nunique().rename('n_lines_abstract')
t = cl[cl.gid.isin(genes)].merge(cells, on='pmid').merge(long, left_on=['ach', 'gid'], right_on=['ach', 'gid'], how='left').merge(nl, on='pmid')
t['in_depmap'] = t.effect.notna()
t.to_parquet('work/triples.parquet', index=False)
print('pmid-gene claims', len(cl), '| human DepMap genes', cl.gid.isin(genes).sum(), '| triples', len(t), 'with effect', int(t.in_depmap.sum()),
      '| pmids', t.pmid.nunique(), '| genes', t.gid.nunique(), '| models', t.ach.nunique())

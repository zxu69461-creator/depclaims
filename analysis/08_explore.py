"""Exploration on DEVELOPMENT papers only (plans/split_2026-11-03.md). Usage: 08_explore.py READER_JSONL
Claims = frozen reader v2, grounded, growth_effect decreased or apoptosis increased, gene and line mapped. DepMap support: scaled Chronos <= -0.5.
Baselines: (a) same gene, other DepMap lines of the same NCIt disease; (b) other claimed genes (random draw from all claimed genes) in the same line."""
import sys, hashlib, subprocess, numpy as np, pandas as pd, duckdb
inp = sys.argv[1]; mp = 'work/reader_mapped_explore.parquet'
subprocess.run([sys.executable, '-I', 'scripts/07_map_reader.py', inp, mp], check=True)
M = pd.read_parquet(mp)
dev = lambda p: hashlib.md5(str(p).encode()).hexdigest()[-1] in '01'
M = M[M.pmid.map(dev)]
C = M[M.grounded & M.gid.notna() & M.ach.notna() & ((M.growth_effect == 'decreased') | M.apoptosis_increased)].drop_duplicates(['pmid', 'gid', 'ach']).copy()
C['gid'] = C.gid.astype(int)
L = pd.read_parquet('work/depmap_long.parquet')
C = C.merge(L[['ach', 'gid', 'effect']], on=['ach', 'gid'], how='inner')
meta = duckdb.sql("select pmid, any_value(pub_year) yr, any_value(country) country, any_value(retracted) retracted, any_value(journal) journal from read_parquet('work/cand/cand_*.parquet') group by 1").df()
C = C.merge(meta, on='pmid', how='left')
lin = pd.read_csv('raw/depmap_lineage.tsv', sep='\t').drop_duplicates('ach')[['ach', 'ncit', 'disease']]
C = C.merge(lin, on='ach', how='left')
C['sup'] = C.effect <= -0.5
print(f'development claims {len(C)} | papers {C.pmid.nunique()} | genes {C.gid.nunique()} | lines {C.ach.nunique()}')
print(f'supported (effect <= -0.5): {C.sup.mean():.3f} | median effect {C.effect.median():.3f}')
# baseline (a): same gene, other lines with the same disease code
Lg = L.merge(lin, on='ach')
a = []
for r in C.itertuples():
    o = Lg[(Lg.gid == r.gid) & (Lg.ncit == r.ncit) & (Lg.ach != r.ach)].effect
    a.append((o <= -0.5).mean() if len(o) else np.nan)
C['base_same_gene_lineage'] = a
# baseline (b): random claimed genes in the same line
rng = np.random.default_rng(20261103); pool = C.gid.unique()
Lw = L.pivot_table(index='ach', columns='gid', values='effect')
C['base_random_gene_line'] = [np.nanmean(Lw.loc[r.ach, rng.choice(pool, 200)] <= -0.5) if r.ach in Lw.index else np.nan for r in C.itertuples()]
print('baselines: same gene other same-disease lines', round(np.nanmean(C.base_same_gene_lineage), 3), '| random claimed genes same line', round(np.nanmean(C.base_random_gene_line), 3))
C['era'] = pd.cut(C.yr, [0, 2009, 2014, 2019, 2030], labels=['<=2009', '2010-14', '2015-19', '2020+'])
C['meth'] = C.method.str.lower().map(lambda m: 'CRISPR' if m and 'crispr' in m else 'RNAi' if m in ('sirna', 'shrna') else 'other')
for k in ['era', 'meth']: print(C.groupby(k, observed=True).agg(n=('sup', 'size'), supported=('sup', 'mean'), base=('base_same_gene_lineage', 'mean')).round(3))
C['cn'] = np.where(C.country == 'China', 'China', np.where(C.country == 'unknown', 'unknown', 'other'))
print(C.groupby('cn').agg(n=('sup', 'size'), supported=('sup', 'mean'), base=('base_same_gene_lineage', 'mean')).round(3))
C.to_parquet('work/explore_dev_claims.parquet', index=False)

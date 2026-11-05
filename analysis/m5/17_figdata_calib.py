"""Figure data for Figs 2 and 5 (M5 job; inputs in mjob IN_DIR, outputs to OUT_DIR). Rulers identical to scripts/12_confirm.py.
calib_decile: non-common confirmatory claims by expression decile of the gene in the named line (deciles within model): claim support
  (bootstrap over papers, 2,000), expression-matched background, pooled background and RNAi-hit ceiling in the same decile.
class_comp: gene-class composition of genes, RNAi-hit gene-line pairs and confirmatory claims.
claims_base: every claim (both splits) with the gene's dependent fraction over all lines and within the named line's lineage.
genes: per gene dependent fraction, class, claims and papers (all claims)."""
import os, numpy as np, pandas as pd
I, O = os.environ['IN_DIR'], os.environ['OUT_DIR']; p = lambda f: os.path.join(I, f)
rng = np.random.default_rng(20261103); B = 2000
ge = pd.read_csv(p('gene_effect_26q1_chronos.csv'), index_col=0); ge = ge[[c for c in ge.columns if '(' in c]]
sym = {int(c.split('(')[-1].rstrip(')')): c.split(' (')[0] for c in ge.columns}
ge.columns = [int(c.split('(')[-1].rstrip(')')) for c in ge.columns]
ceg = pd.read_csv(p('CEGv2.txt'), sep='\t').ENTREZ_ID.dropna().astype(int); neg = pd.read_csv(p('NEGv1.txt'), sep='\t').ENTREZ_ID.dropna().astype(int)
scale = lambda X: X.sub(X[[g for g in neg if g in X.columns]].median(axis=1), axis=0).div(
    X[[g for g in neg if g in X.columns]].median(axis=1) - X[[g for g in ceg if g in X.columns]].median(axis=1), axis=0)
sc = scale(ge)
frac = (sc <= -0.5).mean(axis=0)
gclass = pd.Series(np.where(frac >= 0.9, 'common', np.where(frac > 0, 'selective', 'never')), index=frac.index)
d2 = pd.read_csv(p('d2_gene_effect.csv'), index_col=0); d2 = d2[[c for c in d2.columns if '(' in c and '&' not in c]]
d2.columns = [int(c.split('(')[-1].rstrip(')')) for c in d2.columns]; d2 = d2[d2.index.str.startswith('ACH-')]; d2 = scale(d2[~d2.index.duplicated()])
cr = sc.stack().rename('c').rename_axis(['ach', 'gid']).reset_index(); rn = d2.stack().rename('r').rename_axis(['ach', 'gid']).reset_index()
BR = cr.merge(rn, on=['ach', 'gid']); BR['gclass'] = BR.gid.map(gclass)
Mo = pd.read_csv(p('Model.csv'), usecols=['ModelID', 'OncotreeLineage']).rename(columns={'ModelID': 'ach'})
X = pd.read_csv(p('OmicsExpressionProteinCodingGenesTPMLogp1.csv'), low_memory=False)
if 'IsDefaultEntryForModel' in X.columns: X = X[X.IsDefaultEntryForModel.astype(str).str.lower().isin(['yes', 'true'])]
X = X.set_index('ModelID' if 'ModelID' in X.columns else X.columns[0]); X = X[[c for c in X.columns if '(' in c]]
ids = [int(c.split('(')[-1].rstrip(')')) for c in X.columns]; dup = pd.Series(ids).duplicated(keep=False).values
X = X[[c for c, i, d in zip(X.columns, ids, dup) if not d or c.split(' (')[0] == sym.get(i, c.split(' (')[0])]]
X.columns = [int(c.split('(')[-1].rstrip(')')) for c in X.columns]; X = X.loc[:, ~X.columns.duplicated()]; X = X[~X.index.duplicated()]
EX = X.stack().rename('expr').rename_axis(['ach', 'gid']).reset_index().merge(cr, on=['ach', 'gid'])
EX = EX[EX.gid.map(gclass) != 'common']
EX['dec'] = EX.groupby('ach').expr.transform(lambda v: pd.qcut(v.rank(method='first'), 10, labels=False))
EX['s'] = EX.c <= -0.5
bg = EX.groupby(['ach', 'dec']).s.mean().rename('bg').reset_index()

# calibration by decile (confirmatory, non-common)
C = pd.read_parquet(p('claims_confirmatory.parquet'))
N = C[C.gclass != 'common'].merge(EX[['ach', 'gid', 'dec']], on=['ach', 'gid'], how='inner').merge(bg, on=['ach', 'dec'], how='left')
N['dec'] = N.dec.astype(int)
pp = {k: i for i, k in enumerate(N.pmid.unique())}; pi = N.pmid.map(pp).values
K = np.zeros((len(pp), 10)); M = np.zeros((len(pp), 10)); G = np.zeros((len(pp), 10))
np.add.at(K, (pi, N.dec.values), N.sup.values.astype(float)); np.add.at(M, (pi, N.dec.values), 1.0); np.add.at(G, (pi, N.dec.values), N.bg.fillna(0).values)
bs = np.empty((B, 10)); rb = np.empty(B)
for b in range(B):
    w = np.bincount(rng.integers(0, len(pp), len(pp)), minlength=len(pp))
    k, m, g = w @ K, w @ M, w @ G; bs[b] = k / m; rb[b] = k.sum() / g.sum()
hit = BR[(BR.gclass != 'common') & (BR.r <= -0.5)].merge(EX[['ach', 'gid', 'dec']], on=['ach', 'gid'])
D = pd.DataFrame({'dec': range(10), 'n': M.sum(0).astype(int), 'claim_sup': K.sum(0) / M.sum(0),
                  'lo': np.percentile(bs, 2.5, 0), 'hi': np.percentile(bs, 97.5, 0), 'bg_matched': G.sum(0) / M.sum(0),
                  'bg_pooled': EX.groupby('dec').s.mean().values,
                  'ceiling': hit.assign(s=hit.c <= -0.5).groupby('dec').s.mean().reindex(range(10)).values,
                  'n_hits': hit.groupby('dec').size().reindex(range(10)).values,
                  'expr_med': EX.groupby('dec').expr.median().values})
D.to_parquet(os.path.join(O, 'calib_decile.parquet'), index=False); print(D.round(4).to_string(index=False))
print('ratio to matched background', round(K.sum() / G.sum(), 4), np.percentile(rb, [2.5, 97.5]).round(4))

# class composition
hits = BR[BR.r <= -0.5]
comp = pd.concat([gclass.value_counts(normalize=True).rename('share').to_frame().assign(population='Genes'),
                  hits.gclass.value_counts(normalize=True).rename('share').to_frame().assign(population='RNAi hits'),
                  C.gclass.value_counts(normalize=True).rename('share').to_frame().assign(population='Claims')]).rename_axis('gclass').reset_index()
comp['n'] = comp.population.map({'Genes': len(gclass), 'RNAi hits': len(hits), 'Claims': len(C)})
comp.to_parquet(os.path.join(O, 'class_comp.parquet'), index=False); print(comp.round(4).to_string(index=False))

# per claim gene base rates (both splits)
A = pd.concat([pd.read_parquet(p('claims_development.parquet')).assign(split='development'), C.assign(split='confirmatory')])
L = cr.merge(Mo, on='ach'); L['d'] = L.c <= -0.5
fl = L[L.gid.isin(A.gid.unique())].groupby(['gid', 'OncotreeLineage']).d.agg(['mean', 'size']).rename(columns={'mean': 'frac_lin', 'size': 'n_lin'}).reset_index()
A = A.merge(fl, on=['gid', 'OncotreeLineage'], how='left'); A['frac_all'] = A.gid.map(frac); A['symbol'] = A.gid.map(sym)
A[['pmid', 'gid', 'symbol', 'ach', 'OncotreeLineage', 'split', 'gclass', 'sup', 'v_CRISPR', 'pct_CRISPR', 'expr', 'not_expressed', 'yr', 'meth', 'retracted',
   'frac_all', 'frac_lin', 'n_lin']].to_parquet(os.path.join(O, 'claims_base.parquet'), index=False)
g = A.groupby('gid').agg(n_claims=('pmid', 'size'), n_papers=('pmid', 'nunique'), n_sup=('sup', 'sum'), n_lines=('ach', 'nunique')).reset_index()
Gt = pd.DataFrame({'gid': frac.index, 'symbol': frac.index.map(sym), 'frac': frac.values, 'gclass': gclass.values}).merge(g, on='gid', how='left').fillna({'n_claims': 0, 'n_papers': 0, 'n_sup': 0, 'n_lines': 0})
Gt.to_parquet(os.path.join(O, 'genes.parquet'), index=False)
print('genes', len(Gt), 'claimed', int((Gt.n_claims > 0).sum()), '| claims', len(A))
nc = A[A.gclass != 'common']
print('non-common claims: mean frac_lin', round(nc.frac_lin.mean(), 4), '| support', round(nc.sup.mean(), 4))

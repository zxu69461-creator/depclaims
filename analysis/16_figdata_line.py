"""Figure data for line specificity (Fig. 3): per-claim claimed-line and RNAi-top control percentiles at r_x >= 0.3/0.4/0.5, lineage means,
and the specimen gene strip. Confirmatory claims; rulers and percentile definition identical to 12_confirm.py / 09_specificity.py."""
import numpy as np, pandas as pd
C = pd.read_parquet('db/parquet/claims_confirmatory.parquet')
ge = pd.read_csv('raw/gene_effect_26q1_chronos.csv', index_col=0); ge = ge[[c for c in ge.columns if '(' in c]]
sym = {int(c.split('(')[-1].rstrip(')')): c.split(' (')[0] for c in ge.columns}
ge.columns = [int(c.split('(')[-1].rstrip(')')) for c in ge.columns]
ceg = pd.read_csv('raw/CEGv2.txt', sep='\t').ENTREZ_ID.dropna().astype(int); neg = pd.read_csv('raw/NEGv1.txt', sep='\t').ENTREZ_ID.dropna().astype(int)
scale = lambda X: X.sub(X[[g for g in neg if g in X.columns]].median(axis=1), axis=0).div(
    X[[g for g in neg if g in X.columns]].median(axis=1) - X[[g for g in ceg if g in X.columns]].median(axis=1), axis=0)
genes = C[C.r_x >= 0.3].gid.unique()
sc = scale(ge)[[g for g in genes if g in ge.columns]]
d2 = pd.read_csv('raw/d2_gene_effect.csv', index_col=0); d2 = d2[[c for c in d2.columns if '(' in c and '&' not in c]]
d2.columns = [int(c.split('(')[-1].rstrip(')')) for c in d2.columns]; d2 = d2[d2.index.str.startswith('ACH-')]
d2 = scale(d2[~d2.index.duplicated()]); d2 = d2[[g for g in genes if g in d2.columns]]
Mo = pd.read_csv('raw/Model.csv', usecols=['ModelID', 'OncotreeLineage', 'StrippedCellLineName']).rename(columns={'ModelID': 'ach'})
BL = sc.stack().rename('c').rename_axis(['ach', 'gid']).reset_index().merge(
     d2.stack().rename('r').rename_axis(['ach', 'gid']).reset_index(), on=['ach', 'gid']).merge(Mo, on='ach')
S = C[(C.r_x >= 0.3) & C.pct_CRISPR.notna()].copy()
ctl = {}
for (g, lin), o in BL.merge(S[['gid', 'OncotreeLineage']].drop_duplicates(), on=['gid', 'OncotreeLineage']).groupby(['gid', 'OncotreeLineage']):
    if len(o) < 5: continue
    t = o.loc[o.r.idxmin()]; oth = o[o.ach != t.ach].c
    ctl[(g, lin)] = (oth < t.c).mean() + 0.5 * (oth == t.c).mean()
S['ctl'] = [ctl.get((g, l), np.nan) for g, l in zip(S.gid, S.OncotreeLineage)]
S['symbol'] = S.gid.map(sym)
S[['pmid', 'gid', 'symbol', 'ach', 'OncotreeLineage', 'r_x', 'gclass', 'pct_CRISPR', 'ctl', 'v_CRISPR']].to_parquet('work/fig/line_claims.parquet', index=False)
print('claims', len(S), 'with control', S.ctl.notna().sum(), '| mean', S.pct_CRISPR.mean().round(3), S.ctl.mean().round(3))
for t in (0.3, 0.4, 0.5):
    s = S[S.r_x >= t]; print(t, len(s), s.pct_CRISPR.mean().round(3), s.ctl.mean().round(3))
# specimen candidates: gene-lineage pairs with many distinct named lines and a clear dependent tail
cand = S.groupby(['symbol', 'OncotreeLineage']).agg(n=('pmid', 'nunique'), lines=('ach', 'nunique'), pct=('pct_CRISPR', 'mean'), rx=('r_x', 'first')).reset_index()
print(cand[cand.lines >= 4].sort_values(['lines', 'n'], ascending=False).head(25).to_string(index=False))
BL.merge(S[['gid', 'OncotreeLineage']].drop_duplicates(), on=['gid', 'OncotreeLineage']).assign(symbol=lambda d: d.gid.map(sym)).to_parquet('work/fig/line_strips.parquet', index=False)

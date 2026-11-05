"""Line specificity of literature growth claims (exploration; DEVELOPMENT claims only). Usage: 09_specificity.py CLAIMS_PARQUET
For each claim (gene g, line y): where does y rank among DepMap lines of the same cancer lineage for g?
  pct = within-lineage percentile of y's effect (0 = most dependent line, 0.5 = no line information).
Rulers: CRISPR Chronos 26Q1 (scaled, charter), RNAi DEMETER2 (Achilles+DRIVE+Marcotte), RNA expression 24Q4 (log2 TPM+1).
Gene classes: common essential (scaled effect <= -0.5 in >= 90% of lines), selective (<= -0.5 in 1-89%), never (<= -0.5 in no line).
Positive control: claims on canonical lineage oncogene addictions should rank y near 0."""
import sys, numpy as np, pandas as pd
inp = sys.argv[1]
C = pd.read_parquet(inp)
L = pd.read_parquet('work/depmap_long.parquet')
M = pd.read_csv('raw/Model.csv', usecols=['ModelID', 'CCLEName', 'OncotreeLineage', 'OncotreePrimaryDisease']).rename(columns={'ModelID': 'ach'})
C = C.merge(M[['ach', 'OncotreeLineage']], on='ach', how='left')
# gene classes from CRISPR across all lines
fr = L.assign(d=L.effect <= -0.5).groupby('gid').d.mean()
cls = pd.cut(fr, [-1, 0, 0.8999, 1.01], labels=['never', 'selective', 'common']).rename('gclass')
C = C.merge(cls, left_on='gid', right_index=True, how='left')

def ruler_long():
    out = {'CRISPR': L[['ach', 'gid', 'effect']].rename(columns={'effect': 'v'})}
    d2 = pd.read_csv('raw/d2_gene_effect.csv', index_col=0)            # rows CCLE names, columns 'SYMBOL (entrez)'
    d2 = d2[[c for c in d2.columns if '(' in c and '&' not in c]]       # multi-gene shRNA collapses dropped
    d2.columns = [int(c.split('(')[-1].rstrip(')')) for c in d2.columns]
    d2 = d2[d2.index.str.startswith('ACH-')]                            # rows are already DepMap model ids
    d2 = d2[~d2.index.duplicated()]
    ceg = pd.read_csv('raw/CEGv2.txt', sep='\t').ENTREZ_ID.dropna().astype(int); neg = pd.read_csv('raw/NEGv1.txt', sep='\t').ENTREZ_ID.dropna().astype(int)
    mc = d2[[g for g in ceg if g in d2.columns]].median(axis=1); mn = d2[[g for g in neg if g in d2.columns]].median(axis=1)
    d2 = d2.sub(mn, axis=0).div(mn - mc, axis=0)                        # same per-line scaling as the CRISPR ruler
    out['RNAi'] = d2.stack().rename('v').rename_axis(['ach', 'gid']).reset_index()
    return out

def pct_table(R, C, lab):
    R = R.merge(M[['ach', 'OncotreeLineage']], on='ach')
    R = R[R.gid.isin(C.gid.unique())]
    g = {k: v for k, v in R.groupby(['gid', 'OncotreeLineage'])}
    rows = []
    for r in C.itertuples():
        o = g.get((r.gid, r.OncotreeLineage))
        if o is None or r.ach not in set(o.ach): rows.append((np.nan, np.nan, np.nan)); continue
        y = o.loc[o.ach == r.ach, 'v'].iloc[0]; oth = o[o.ach != r.ach].v
        if len(o) < 5: rows.append((np.nan, y, np.nan)); continue          # effect kept; rank only within lineages of >= 5 lines
        rows.append(((oth < y).mean() + 0.5 * (oth == y).mean(), y, y - oth.median()))
    return pd.DataFrame(rows, columns=[f'pct_{lab}', f'v_{lab}', f'd_{lab}'], index=C.index)

RL = ruler_long()
for lab, R in RL.items(): C = C.join(pct_table(R, C, lab))
# expression in the claimed line
import os
if os.path.exists('raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv'):
    E = pd.read_csv('raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0, nrows=0).columns
    meta = [c for c in E if '(' not in c]
    X = pd.read_csv('raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv', usecols=lambda c: c in ('Unnamed: 0',) or c in meta or '(' in c, low_memory=False)
    if 'ModelID' in X.columns:
        if 'IsDefaultEntryForModel' in X.columns: X = X[X.IsDefaultEntryForModel.astype(str).str.lower().isin(['yes', 'true'])]
        X = X.set_index('ModelID')
    else:
        X = X.set_index(X.columns[0])
    X = X[[c for c in X.columns if '(' in c]]
    # 17 Entrez ids label two expression columns (readthrough/antisense loci): keep the one whose symbol matches the CRISPR release
    sym = {int(c.split('(')[-1].rstrip(')')): c.split(' (')[0] for c in open('raw/gene_effect_26q1_chronos.csv').readline().strip().split(',')[1:] if '(' in c}
    ids = [int(c.split('(')[-1].rstrip(')')) for c in X.columns]
    dup = pd.Series(ids).duplicated(keep=False).values
    X = X[[c for c, i, d in zip(X.columns, ids, dup) if not d or c.split(' (')[0] == sym.get(i, c.split(' (')[0])]]
    X.columns = [int(c.split('(')[-1].rstrip(')')) for c in X.columns]; X = X.loc[:, ~X.columns.duplicated()]
    X = X[~X.index.duplicated()]
    C['expr'] = [X.at[a, g] if (a in X.index and g in X.columns) else np.nan for a, g in zip(C.ach, C.gid)]
    C['not_expressed'] = C.expr < 1.0   # log2(TPM+1) < 1, i.e. TPM < 1
else:
    C['expr'] = np.nan; C['not_expressed'] = np.nan
C.to_parquet(inp.replace('.parquet', '_spec.parquet'), index=False)

def summ(df, lab):
    s = df[f'pct_{lab}'].dropna(); n = len(s)
    return f'{lab}: n {n} | mean pct {s.mean():.3f} +- {s.std() / np.sqrt(max(n, 1)):.3f} | in most-dependent quartile {(s < .25).mean():.3f}'
print('claims', len(C), '| lineage mapped', C.OncotreeLineage.notna().mean().round(3))
for lab in RL: print(summ(C, lab))
print('gene class:'); print(C.groupby('gclass', observed=True).agg(n=('sup', 'size'), supported=('sup', 'mean'), pct_CRISPR=('pct_CRISPR', 'mean'), pct_RNAi=('pct_RNAi', 'mean')).round(3))
print('expression in claimed line: n', C.expr.notna().sum(), '| TPM < 1:', round(float(C.not_expressed[C.expr.notna()].mean()), 3), '| median log2(TPM+1)', round(float(C.expr.median()), 2))
print('RNAi support (D2 <= -0.5):', (C.v_RNAi <= -0.5)[C.v_RNAi.notna()].mean().round(3), 'n', C.v_RNAi.notna().sum(),
      '| either ruler', ((C.v_RNAi <= -0.5) | C.sup)[C.v_RNAi.notna()].mean().round(3))
print('by method:'); print(C.groupby('meth').agg(n=('sup', 'size'), pct_CRISPR=('pct_CRISPR', 'mean'), pct_RNAi=('pct_RNAi', 'mean')).round(3))
# positive control for the metric: the line an independent screen (RNAi) ranks most dependent in the lineage, ranked under CRISPR
cr = RL['CRISPR'].merge(M[['ach', 'OncotreeLineage']], on='ach'); rn = RL['RNAi'].merge(M[['ach', 'OncotreeLineage']], on='ach')
both = cr.merge(rn, on=['ach', 'gid', 'OncotreeLineage'], suffixes=('_c', '_r'))
pairs = C[['gid', 'OncotreeLineage', 'gclass']].drop_duplicates(['gid', 'OncotreeLineage'])
pc = []
for r in pairs.itertuples():
    o = both[(both.gid == r.gid) & (both.OncotreeLineage == r.OncotreeLineage)]
    if len(o) < 5: continue
    top = o.loc[o.v_r.idxmin()]; oth = o[o.ach != top.ach].v_c
    pc.append((r.gclass, (oth < top.v_c).mean() + 0.5 * (oth == top.v_c).mean()))
pc = pd.DataFrame(pc, columns=['gclass', 'pct'])
print('positive control (RNAi-top line ranked by CRISPR, same gene-lineage pairs): n', len(pc), '| mean pct', round(pc.pct.mean(), 3), '| by class', pc.groupby('gclass', observed=True).pct.mean().round(3).to_dict())
s = C[['pct_CRISPR', 'pct_RNAi']].dropna(); print('claimed line: CRISPR vs RNAi percentile Spearman', round(s.corr('spearman').iloc[0, 1], 3), 'n', len(s))

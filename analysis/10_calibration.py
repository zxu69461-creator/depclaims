"""Ruler calibration (exploration; DEVELOPMENT claims). Usage: 10_calibration.py  (reads work/explore_dev_claims_spec.parquet)
Ceiling for agreement: how often a DepMap RNAi (DEMETER2) hit is a CRISPR hit in the same gene-line, overall and outside common essentials;
gene-level CRISPR-RNAi agreement across lines (r_x); power of the within-lineage percentile where the screens agree (RNAi-top line ranked by CRISPR)."""
import pandas as pd, numpy as np, runpy, sys
C = pd.read_parquet('work/explore_dev_claims_spec.parquet')
L = pd.read_parquet('work/depmap_long.parquet')[['ach','gid','effect']]
d2 = pd.read_csv('raw/d2_gene_effect.csv', index_col=0); d2 = d2[[c for c in d2.columns if '(' in c and '&' not in c]]
d2.columns = [int(c.split('(')[-1].rstrip(')')) for c in d2.columns]; d2 = d2[d2.index.str.startswith('ACH-')]
ceg = pd.read_csv('raw/CEGv2.txt', sep='\t').ENTREZ_ID.dropna().astype(int); neg = pd.read_csv('raw/NEGv1.txt', sep='\t').ENTREZ_ID.dropna().astype(int)
mc = d2[[g for g in ceg if g in d2.columns]].median(axis=1); mn = d2[[g for g in neg if g in d2.columns]].median(axis=1)
d2 = d2.sub(mn, axis=0).div(mn-mc, axis=0)
R = d2.stack().rename('rnai').rename_axis(['ach','gid']).reset_index()
B = L.merge(R, on=['ach','gid'])
print('paired gene-line cells', len(B), 'genes', B.gid.nunique(), 'lines', B.ach.nunique())
# ceiling: RNAi hit -> CRISPR hit (all claimed genes), and excluding common essentials
fr = L.assign(d=L.effect<=-.5).groupby('gid').d.mean()
B['common'] = B.gid.map(fr) >= .9
for lab, b in [('all claimed genes', B), ('non-common', B[~B.common])]:
    h = b[b.rnai<=-.5]; print(f'{lab}: RNAi hits {len(h)} -> CRISPR <= -0.5 {(h.effect<=-.5).mean():.3f} | CRISPR hits -> RNAi {(b[b.effect<=-.5].rnai<=-.5).mean():.3f} | base CRISPR rate {(b.effect<=-.5).mean():.3f}')
# literature claims, same restriction
c = C.dropna(subset=['v_RNAi']); c = c.assign(common=c.gclass=='common')
for lab, x in [('claims all', c), ('claims non-common', c[~c.common])]:
    print(f'{lab}: n {len(x)} CRISPR {x.sup.mean():.3f} RNAi {(x.v_RNAi<=-.5).mean():.3f}')
# gene-level cross-screen agreement
r = B.groupby('gid').apply(lambda g: g.effect.corr(g.rnai) if len(g)>50 else np.nan).rename('r_x')
C2 = C.merge(r, left_on='gid', right_index=True, how='left')
for t in [.2,.3,.4]:
    x = C2[C2.r_x>=t]; print(f'genes with CRISPR-RNAi r>={t}: claims {len(x)} pct_CRISPR {x.pct_CRISPR.mean():.3f} pct_RNAi {x.pct_RNAi.mean():.3f} support {x.sup.mean():.3f}')
M = pd.read_csv('raw/Model.csv', usecols=['ModelID','OncotreeLineage']).rename(columns={'ModelID':'ach'})
Bl = B.merge(M, on='ach')
rg = r.dropna()
for t in [0, .3, .5]:
    gs = set(rg[rg>=t].index)
    pcs=[]
    for (g, lin), o in Bl[Bl.gid.isin(gs)].groupby(['gid','OncotreeLineage']):
        if len(o)<5: continue
        top = o.loc[o.rnai.idxmin()]; oth = o[o.ach!=top.ach].effect
        pcs.append((oth<top.effect).mean())
    print(f'all gene-lineage pairs, genes r>={t}: n {len(pcs)} RNAi-top line CRISPR pct {np.mean(pcs):.3f} | genes {len(gs)}')
print('claimed genes r distribution', C2.drop_duplicates('gid').r_x.describe().round(2).to_dict())

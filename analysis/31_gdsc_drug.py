"""Additional analysis C of notes/analysis_plan_review2_2026-11-04.md (not prespecified): pharmacological inhibition (GDSC release 8.5, GDSC1 and
GDSC2 fitted dose response, Sanger) as an external test of screen-scored claims. Compounds whose putative targets name the claimed gene (HGNC
symbol or alias) and list at most two targets; lines mapped to DepMap by COSMIC identifier. Per compound and lineage, a line's sensitivity
rank among the lineage's lines (AUC; 0 = most sensitive), averaged over the gene's compounds; >= 5 lineage lines with data.
C1 named-line drug rank, CRISPR-supported vs not supported claims. C2 the line the CRISPR screen ranks most dependent in the lineage (other
than the named line, with drug data) vs the named line. Confirmatory claims; bootstrap over genes."""
import sys, re, numpy as np, pandas as pd
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
from rdata.proj import Project
P = Project(__file__); rng = np.random.default_rng(20261104); B = 2000
G = '/Volumes/ResearchData/Claude/store/raw/GDSC_rel8.5/'
out = []
add = lambda a, k, n, e, lo=np.nan, hi=np.nan: out.append(dict(analysis=a, quantity=k, n=n, est=e, lo=lo, hi=hi))
with P.run(__file__, params=dict(B=B, seed=20261104, max_targets=2, min_lines=5)) as run:
    for f in ('GDSC1', 'GDSC2'):
        run.input(f, path=G + f'{f}_fitted_dose_response_27Oct23.xlsx', note=f'{f} fitted dose response, release 8.5')
    D = pd.concat([pd.read_excel(G + f'{f}_fitted_dose_response_27Oct23.xlsx', usecols=['DATASET', 'COSMIC_ID', 'DRUG_ID', 'DRUG_NAME',
                                                                                     'PUTATIVE_TARGET', 'AUC']) for f in ('GDSC1', 'GDSC2')])
    M = pd.read_csv('raw/Model.csv', usecols=['ModelID', 'OncotreeLineage', 'COSMICID']).dropna(subset=['COSMICID'])
    M['COSMICID'] = M.COSMICID.astype(int)
    D = D.merge(M.rename(columns={'COSMICID': 'COSMIC_ID'}), on='COSMIC_ID')
    H = pd.read_csv('/Volumes/ResearchData/Claude/ref/raw/hgnc_complete_set.txt', sep='\t', usecols=['symbol', 'alias_symbol', 'prev_symbol', 'entrez_id'],
                    low_memory=False).dropna(subset=['entrez_id'])
    name2id = {}
    for r in H.itertuples():
        name2id.setdefault(r.symbol.upper(), int(r.entrez_id))
    for r in H.itertuples():
        for col in (r.alias_symbol, r.prev_symbol):
            if isinstance(col, str):
                for a in col.split('|'): name2id.setdefault(a.upper(), int(r.entrez_id))
    tg = D[['DATASET', 'DRUG_ID', 'PUTATIVE_TARGET']].drop_duplicates()
    rows = []
    for r in tg.itertuples():
        ts = [t.strip().upper() for t in re.split(r'[,;/]', str(r.PUTATIVE_TARGET)) if t.strip()]
        if not 1 <= len(ts) <= 2: continue
        for t in ts:
            if t in name2id: rows.append((r.DATASET, r.DRUG_ID, name2id[t]))
    T = pd.DataFrame(rows, columns=['DATASET', 'DRUG_ID', 'entrez_id']).drop_duplicates()
    D['cpd'] = D.DATASET + ':' + D.DRUG_ID.astype(str)
    D = D.groupby(['cpd', 'DATASET', 'DRUG_ID', 'ModelID', 'OncotreeLineage'], as_index=False).AUC.mean()
    D['n_lin'] = D.groupby(['cpd', 'OncotreeLineage']).ModelID.transform('size')
    D = D[D.n_lin >= 5].copy()
    D['drug_rank'] = D.groupby(['cpd', 'OncotreeLineage']).AUC.transform(lambda a: (a.rank(method='average') - 1) / (len(a) - 1))
    R = D.merge(T, on=['DATASET', 'DRUG_ID']).groupby(['entrez_id', 'ModelID'], as_index=False).agg(drug_rank=('drug_rank', 'mean'), n_cpd=('cpd', 'nunique'))
    add('C coverage', 'genes with targeted compounds', int(R.entrez_id.nunique()), np.nan)

    A = pd.read_csv('release/claims_scored.tsv.gz', sep='\t')
    C = A[(A.split == 'confirmatory') & A.crispr_effect_scaled.notna()].drop_duplicates(['entrez_id', 'depmap_model']).copy()
    C = C.merge(R.rename(columns={'ModelID': 'depmap_model'}), on=['entrez_id', 'depmap_model']).reset_index(drop=True)
    add('C coverage', 'claims (gene-line pairs) with drug data', len(C), np.nan)
    add('C coverage', 'genes', int(C.entrez_id.nunique()), np.nan)
    def gboot(d, f):
        g = {k: v.index.values for k, v in d.groupby('entrez_id')}; kk = np.array(list(g))
        est = f(d); res = [f(d.loc[np.concatenate([g[k] for k in rng.choice(kk, len(kk))])]) for _ in range(B)]
        lo, hi = np.nanpercentile(res, [2.5, 97.5]); return est, lo, hi
    for k, m in (('supported', C.crispr_supported), ('not supported', ~C.crispr_supported.astype(bool))):
        d = C[m].reset_index(drop=True); add('C1 named-line drug-sensitivity rank', k, len(d), *gboot(d, lambda x: x.drug_rank.mean()))
    add('C1 named-line drug-sensitivity rank', 'supported minus not supported', len(C),
        *gboot(C, lambda x: x.loc[x.crispr_supported.astype(bool), 'drug_rank'].mean() - x.loc[~x.crispr_supported.astype(bool), 'drug_rank'].mean()))
    rho = C[['crispr_effect_scaled', 'drug_rank']].corr('spearman').iloc[0, 1]
    add('C1 Spearman rho, CRISPR effect vs drug rank (named lines)', 'all', len(C), *gboot(C, lambda x: x[['crispr_effect_scaled', 'drug_rank']].corr('spearman').iloc[0, 1]))

    cr = pd.read_parquet('work/tool_cache/crispr.parquet')
    lin = pd.read_parquet('work/tool_cache/models.parquet').set_index('ModelID').OncotreeLineage
    pick = []
    for r in C.itertuples():
        g = str(r.entrez_id)
        if g not in cr.columns: pick.append(None); continue
        cand = R[(R.entrez_id == r.entrez_id) & (R.ModelID != r.depmap_model)].ModelID
        cand = [m for m in cand if lin.get(m) == r.lineage and m in cr.index and pd.notna(cr.at[m, g])]
        pick.append(min(cand, key=lambda m: cr.at[m, g]) if cand else None)
    C['screen_line'] = pick
    E = C[C.screen_line.notna()].merge(R.rename(columns={'ModelID': 'screen_line', 'drug_rank': 'drug_rank_screen', 'n_cpd': 'n2'}),
                                       on=['entrez_id', 'screen_line']).reset_index(drop=True)
    add('C2 drug-sensitivity rank', 'named line', len(E), *gboot(E, lambda x: x.drug_rank.mean()))
    add('C2 drug-sensitivity rank', 'screen-chosen line', len(E), *gboot(E, lambda x: x.drug_rank_screen.mean()))
    add('C2 drug-sensitivity rank', 'screen-chosen minus named', len(E), *gboot(E, lambda x: (x.drug_rank_screen - x.drug_rank).mean()))
    add('C2 share where screen-chosen line is more sensitive', 'all', len(E), *gboot(E, lambda x: (x.drug_rank_screen < x.drug_rank).mean()))
    O = pd.DataFrame(out); print(O.round(4).to_string(index=False))
    run.table('gdsc_drug_additional', O, note='GDSC drug response as an external test (additional analysis, plan 2026-11-04)')
    E.to_parquet('work/gdsc_pairs.parquet', index=False)

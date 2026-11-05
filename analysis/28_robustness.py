"""Additional analyses (not prespecified; added after an internal pre-submission review). Confirmatory claims unless stated.
R1 line rank: named line vs (a) a popularity-weighted null (lines drawn in proportion to how often the literature names them, other genes'
   claims only), (b) the line with the highest expression of the gene in the lineage (information available to authors), split by era
   (before 2017 / 2017 on, DepMap Achilles release); genes with CRISPR-RNAi r >= 0.3; bootstrap over genes.
R2 replication within gene: Mantel-Haenszel odds ratio, independent repeat vs single claims, stratified by gene and by gene x lineage.
R3 author matching on full forename + surname (stricter than initial): repeat classes and support.
R4 dependency cutoff -0.3 / -0.7 / -1.0: support by repeat class and the expression-matched calibrated ratio."""
import sys, unicodedata, itertools, numpy as np, pandas as pd
from statsmodels.stats.contingency_tables import StratifiedTable
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
from rdata.proj import Project
P = Project(__file__); rng = np.random.default_rng(20261104); B = 2000
fold = lambda s: ''.join(c for c in unicodedata.normalize('NFKD', s or '') if not unicodedata.combining(c)).lower().replace('-', ' ').strip()
out = []
add = lambda a, k, n, e, lo=np.nan, hi=np.nan: out.append(dict(analysis=a, quantity=k, n=n, est=e, lo=lo, hi=hi))
with P.run(__file__, params=dict(B=B, seed=20261104)) as run:
    run.input('atlas_claims', path='release/claims_scored.tsv.gz', note='atlas release table')
    A = pd.read_csv('release/claims_scored.tsv.gz', sep='\t')
    Cf = pd.read_parquet('db/parquet/claims_confirmatory.parquet')[['pmid', 'gid', 'ach', 'r_x']].rename(columns={'gid': 'entrez_id', 'ach': 'depmap_model'})
    C = A[A.split == 'confirmatory'].merge(Cf, on=['pmid', 'entrez_id', 'depmap_model'])
    cr = pd.read_parquet('work/tool_cache/crispr.parquet'); ex = pd.read_parquet('work/tool_cache/expression.parquet')
    M = pd.read_parquet('work/tool_cache/models.parquet').set_index('ModelID')
    lin = M.OncotreeLineage

    # R1 line rank
    S = C[(C.r_x >= 0.3) & C.crispr_rank_in_lineage.notna()].copy()
    use = A.groupby(['depmap_model', 'entrez_id']).size().rename('k').reset_index()
    tot = use.groupby('depmap_model').k.sum()
    def ranks(v):
        v = v.dropna(); r = {}
        for m, y in v.items():
            oth = v.drop(m); r[m] = (oth < y).mean() + 0.5 * (oth == y).mean()
        return pd.Series(r)
    cache = {}; pop, exr = [], []
    for row in S.itertuples():
        key = (row.entrez_id, row.lineage)
        if key not in cache:
            same = lin[lin == row.lineage].index
            v = cr[str(row.entrez_id)].reindex(same) if str(row.entrez_id) in cr.columns else pd.Series(dtype=float)
            rk = ranks(v)
            w = (tot.reindex(rk.index).fillna(0) - use[use.entrez_id == row.entrez_id].set_index('depmap_model').k.reindex(rk.index).fillna(0)).clip(lower=0)
            p_null = float((rk * w).sum() / w.sum()) if w.sum() > 0 else np.nan
            e = ex[str(row.entrez_id)].reindex(rk.index).dropna() if str(row.entrez_id) in ex.columns else pd.Series(dtype=float)
            e_top = float(rk[e.idxmax()]) if len(e) else np.nan
            cache[key] = (p_null, e_top)
        pop.append(cache[key][0]); exr.append(cache[key][1])
    S['pop_null'] = pop; S['expr_top'] = exr
    S = S.reset_index(drop=True)
    g = {k: v.index.values for k, v in S.groupby('entrez_id')}; keys = np.array(list(g))
    def bs(cols, sub=None):
        d = S if sub is None else S[sub].reset_index(drop=True)
        gg = {k: v.index.values for k, v in d.groupby('entrez_id')}; kk = np.array(list(gg))
        est = [d[c].mean() for c in cols]
        res = np.array([[d.loc[ix, c].mean() for c in cols] for ix in
                        (np.concatenate([gg[k] for k in rng.choice(kk, len(kk))]) for _ in range(B))])
        return est, np.nanpercentile(res, [2.5, 97.5], axis=0)
    est, ci = bs(['crispr_rank_in_lineage', 'pop_null', 'expr_top'])
    for i, k in enumerate(['named line', 'popularity-weighted null', 'highest-expression line']):
        add('R1 line rank (r_x >= 0.3)', k, len(S), est[i], *ci[:, i])
    S['d_pop'] = S.crispr_rank_in_lineage - S.pop_null
    e2, c2 = bs(['d_pop']); add('R1 line rank (r_x >= 0.3)', 'named minus popularity null', len(S), e2[0], *c2[:, 0])
    for lab, m in (('before 2017', S.year < 2017), ('2017 on', S.year >= 2017)):
        e3, c3 = bs(['crispr_rank_in_lineage'], m); add('R1 named-line rank by era', lab, int(m.sum()), e3[0], *c3[:, 0])

    # R2 replication within gene (non-common genes; independent repeat vs single)
    R = C[(C.gene_class != 'common') & C.repeat_class.isin(['independent repeat', 'single'])].copy()
    R['ind'] = R.repeat_class == 'independent repeat'
    for lab, strata in (('gene', ['entrez_id']), ('gene x lineage', ['entrez_id', 'lineage'])):
        tabs = []
        for _, d in R.groupby(strata):
            if d.ind.nunique() < 2 or d.crispr_supported.nunique() < 2: continue
            t = np.array([[((d.ind) & d.crispr_supported).sum(), ((d.ind) & ~d.crispr_supported).sum()],
                          [((~d.ind) & d.crispr_supported).sum(), ((~d.ind) & ~d.crispr_supported).sum()]], float)
            tabs.append(t)
        st = StratifiedTable(tabs); lo, hi = st.oddsratio_pooled_confint()
        add(f'R2 MH odds ratio within {lab}', 'independent repeat vs single', int(sum(t.sum() for t in tabs)), st.oddsratio_pooled, lo, hi)
        add(f'R2 MH odds ratio within {lab}', 'strata', len(tabs), np.nan)

    # R3 stricter author matching (full forename)
    au = pd.read_parquet('work/dcauth/out/authors.parquet')
    au['key'] = au['last'].map(fold) + '|' + au.fore.map(fold)
    sets = au.groupby('pmid').apply(lambda d: set(d.key) | {'orcid:' + o for o in d.orcid.dropna()}, include_groups=False).to_dict()
    cls = {}
    for (gid, m), d in C.groupby(['entrez_id', 'depmap_model']):
        ps = sorted(d.pmid.unique())
        cls[(gid, m)] = 'single' if len(ps) == 1 else ('independent repeat' if any(not (sets.get(x, set()) & sets.get(y, set())) for x, y in itertools.combinations(ps, 2)) else 'shared-author repeat')
    C['rc_full'] = [cls[(a, b)] for a, b in zip(C.entrez_id, C.depmap_model)]
    for k, d in C.groupby('rc_full'):
        add('R3 repeat class, full-forename matching', k, len(d), d.crispr_supported.mean())
    add('R3 repeat class, full-forename matching', 'claims changing class', len(C), (C.rc_full != C.repeat_class).mean())

    # R4 cutoffs
    for t in (-0.3, -0.5, -0.7, -1.0):
        s = C.crispr_effect_scaled <= t
        for k in ('independent repeat', 'single', 'shared-author repeat'):
            m = C.repeat_class == k; add(f'R4 support at cutoff {t}', k, int(m.sum()), s[m].mean())
    # expression-matched background at each cutoff (non-common genes, same line, same within-line expression decile)
    ge = cr.loc[:, cr.columns.isin(ex.columns)]; ee = ex.reindex(index=ge.index, columns=ge.columns)
    frac = (cr <= -0.5).sum() / len(cr); common = set(frac[frac >= 0.9].index)
    keep = [c for c in ge.columns if c not in common]
    ge, ee = ge[keep], ee[keep]
    dec = ee.rank(axis=1, method='first').apply(lambda r: pd.qcut(r, 10, labels=False) if r.notna().sum() >= 10 else r * np.nan, axis=1)
    N = C[C.gene_class != 'common']
    for t in (-0.3, -0.5, -0.7, -1.0):
        hit = (ge <= t)
        bgs, sup = [], []
        for row in N.itertuples():
            m, col = row.depmap_model, str(row.entrez_id)
            if m not in dec.index or col not in dec.columns or pd.isna(dec.at[m, col]): continue
            dd = dec.loc[m]; sel = dd == dec.at[m, col]
            bgs.append(hit.loc[m, sel.index[sel]].mean()); sup.append(row.crispr_effect_scaled <= t)
        add('R4 calibrated ratio (claims / expression-matched background)', f'cutoff {t}', len(sup), np.mean(sup) / np.mean(bgs))
    O = pd.DataFrame(out)
    print(O.round(4).to_string(index=False))
    run.table('robustness_additional', O, note='additional (not prespecified) robustness analyses after internal review')

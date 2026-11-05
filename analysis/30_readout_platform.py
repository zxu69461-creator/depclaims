"""Additional analyses A and B of notes/analysis_plan_review2_2026-11-04.md (not prespecified). Confirmatory claims.
A readout type of the verbatim growth sentence (clonal or in vivo > proliferation or viability > apoptosis or cell cycle only): CRISPR support
  among non-common genes and its ratio to the expression-matched background, named-line rank (r_x >= 0.3), support of not-expressed claims.
B CRISPR vs RNAi in the same claims: support on each platform, both, either, kappa, effect correlation; reagent-matched ruler; RNAi-ruler
  calibrated ratio, named-line rank and not-expressed support.
Usage: 30_readout_platform.py [READER_JSONL MAPPED_PARQUET] (defaults: the production reader files)."""
import sys, re, json, numpy as np, pandas as pd
from sklearn.metrics import cohen_kappa_score
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
from rdata.proj import Project
P = Project(__file__); rng = np.random.default_rng(20261104); B = 2000
rj, mp = (sys.argv[1:3] if len(sys.argv) > 2 else ('work/reader_all.jsonl', 'work/reader_mapped_explore.parquet'))
CLON = re.compile(r'colon(y|ies)|clonogen|soft agar|anchorage|xenograft|tumou?r growth|tumorigen', re.I)
PROL = re.compile(r'prolifer|viab|growth|\bMTT\b|CCK-?8|\bEdU\b|BrdU|cell number|surviv', re.I)
APCC = re.compile(r'apoptos|cell death|arrest|cell[- ]cycle|\bG[012](/[GSM]1?)?\b|\bS[- ]phase', re.I)
def readout(t):
    return 'clonal or in vivo' if CLON.search(t) else 'proliferation or viability' if PROL.search(t) else 'apoptosis or cell cycle only' if APCC.search(t) else 'other'
out = []
add = lambda a, k, n, e, lo=np.nan, hi=np.nan: out.append(dict(analysis=a, quantity=k, n=n, est=e, lo=lo, hi=hi))
with P.run(__file__, params=dict(B=B, seed=20261104, reader=rj)) as run:
    run.input('atlas_claims', path='release/claims_scored.tsv.gz', note='atlas release table')
    A = pd.read_csv('release/claims_scored.tsv.gz', sep='\t')
    Cf = pd.read_parquet('db/parquet/claims_confirmatory.parquet')[['pmid', 'gid', 'ach', 'r_x']].rename(columns={'gid': 'entrez_id', 'ach': 'depmap_model'})
    C = A[A.split == 'confirmatory'].merge(Cf, on=['pmid', 'entrez_id', 'depmap_model']).reset_index(drop=True)
    # verbatim growth sentences, aligned row by row with the mapped reader output
    M = pd.read_parquet(mp); qs = []
    for l in open(rj):
        d = json.loads(l)
        for c in d['claims'] or []:
            if isinstance(c, dict): qs.append(c.get('growth_quote') or '')
    assert len(qs) == len(M); M['q'] = qs
    M = M[M.gid.notna() & M.ach.notna()].assign(entrez_id=lambda m: m.gid.astype(int), depmap_model=lambda m: m.ach)
    Qt = M.groupby(['pmid', 'entrez_id', 'depmap_model']).q.agg(' | '.join).rename('quote').reset_index()
    C = C.merge(Qt, on=['pmid', 'entrez_id', 'depmap_model'], how='left'); C['readout'] = C.quote.fillna('').map(readout)
    for k, d in C.groupby('readout'): add('A readout share', k, len(d), len(d) / len(C))

    cr = pd.read_parquet('work/tool_cache/crispr.parquet'); rn = pd.read_parquet('work/tool_cache/rnai.parquet')
    ex = pd.read_parquet('work/tool_cache/expression.parquet')
    def background(mat, rows, t=-0.5):
        """expression-matched background per claim: share of non-common genes in the same line and within-line expression decile with
        effect <= t (common = dependent in >= 90% of the platform's lines)."""
        ge = mat.loc[:, mat.columns.isin(ex.columns)]
        frac = (mat <= -0.5).sum() / len(mat); ge = ge[[c for c in ge.columns if frac.get(c, 0) < 0.9]]
        ee = ex.reindex(index=ge.index, columns=ge.columns)
        dec = ee.rank(axis=1, method='first').apply(lambda r: pd.qcut(r, 10, labels=False) if r.notna().sum() >= 10 else r * np.nan, axis=1)
        hit = ge <= t; bg = []
        for m, g in zip(rows.depmap_model, rows.entrez_id.astype(str)):
            if m not in dec.index or g not in dec.columns or pd.isna(dec.at[m, g]): bg.append(np.nan); continue
            sel = dec.loc[m] == dec.at[m, g]; bg.append(hit.loc[m, sel.index[sel]].mean())
        return np.array(bg)
    N = C[C.gene_class != 'common'].reset_index(drop=True)
    N['bg_c'] = background(cr, N); N['sup_r'] = N.rnai_effect_scaled <= -0.5
    N['bg_r'] = np.where(N.rnai_effect_scaled.notna(), background(rn, N), np.nan)
    def gboot(d, f):
        g = {k: v.index.values for k, v in d.groupby('entrez_id')}; kk = np.array(list(g))
        est = f(d); res = [f(d.loc[np.concatenate([g[k] for k in rng.choice(kk, len(kk))])]) for _ in range(B)]
        lo, hi = np.nanpercentile(res, [2.5, 97.5]); return est, lo, hi
    ratio = lambda s, b: (lambda d: d[s].mean() / d[b].mean())
    for k, d in N.groupby('readout'):
        d = d[d.bg_c.notna()].reset_index(drop=True)
        add('A CRISPR support, non-common genes', k, len(d), d.crispr_supported.mean())
        add('A calibrated ratio (CRISPR, expression-matched)', k, len(d), *gboot(d, ratio('crispr_supported', 'bg_c')))
    S = C[(C.r_x >= 0.3) & C.crispr_rank_in_lineage.notna()].reset_index(drop=True)
    for k, d in S.groupby('readout'):
        add('A named-line rank (CRISPR, r_x >= 0.3)', k, len(d), *gboot(d.reset_index(drop=True), lambda x: x.crispr_rank_in_lineage.mean()))
    for k, d in C[C.not_expressed == True].groupby('readout'):
        add('A support of claims on not-expressed genes', k, len(d), d.crispr_supported.mean())

    # B cross-platform
    D = C[C.crispr_effect_scaled.notna() & C.rnai_effect_scaled.notna()].copy()
    c, r = D.crispr_effect_scaled <= -0.5, D.rnai_effect_scaled <= -0.5
    for k, v in (('CRISPR', c), ('RNAi', r), ('both', c & r), ('either', c | r)): add('B support, claims with both screens', k, len(D), v.mean())
    add('B platform agreement', 'Cohen kappa', len(D), cohen_kappa_score(c, r))
    add('B platform agreement', 'raw agreement', len(D), (c == r).mean())
    add('B platform agreement', 'Pearson r of scaled effects', len(D), np.corrcoef(D.crispr_effect_scaled, D.rnai_effect_scaled)[0, 1])
    for k, v in (('RNAi reagent, RNAi screen', (C.perturbation == 'RNAi', 'rnai')), ('RNAi reagent, CRISPR screen', (C.perturbation == 'RNAi', 'crispr')),
                 ('CRISPR reagent, CRISPR screen', (C.perturbation == 'CRISPR', 'crispr')), ('CRISPR reagent, RNAi screen', (C.perturbation == 'CRISPR', 'rnai'))):
        m, pl = v; e = C.loc[m, f'{pl}_effect_scaled'].dropna(); add('B reagent-matched ruler', k, len(e), (e <= -0.5).mean())
    Nr = N[N.bg_r.notna()].reset_index(drop=True)
    add('B calibrated ratio (RNAi, expression-matched)', 'all readouts', len(Nr), *gboot(Nr, ratio('sup_r', 'bg_r')))
    Nc = N[N.bg_c.notna()].reset_index(drop=True)
    add('B calibrated ratio (CRISPR, expression-matched)', 'all readouts', len(Nc), *gboot(Nc, ratio('crispr_supported', 'bg_c')))
    Sr = C[(C.r_x >= 0.3) & C.rnai_rank_in_lineage.notna()].reset_index(drop=True)
    add('B named-line rank (RNAi, r_x >= 0.3)', 'all readouts', len(Sr), *gboot(Sr, lambda x: x.rnai_rank_in_lineage.mean()))
    ne = C[C.not_expressed == True]
    add('B support of claims on not-expressed genes', 'CRISPR', int(ne.crispr_effect_scaled.notna().sum()), (ne.crispr_effect_scaled.dropna() <= -0.5).mean())
    add('B support of claims on not-expressed genes', 'RNAi', int(ne.rnai_effect_scaled.notna().sum()), (ne.rnai_effect_scaled.dropna() <= -0.5).mean())
    ex_ = C[C.not_expressed == False]
    add('B support of claims on expressed genes', 'RNAi', int(ex_.rnai_effect_scaled.notna().sum()), (ex_.rnai_effect_scaled.dropna() <= -0.5).mean())
    O = pd.DataFrame(out); print(O.round(4).to_string(index=False))
    run.table('readout_platform_additional', O, note='readout-type strata and CRISPR-RNAi consistency (additional analyses, plan 2026-11-04)')
    C[['pmid', 'entrez_id', 'depmap_model', 'readout']].to_parquet('work/claims_readout.parquet', index=False)

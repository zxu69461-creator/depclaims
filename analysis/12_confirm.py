"""Confirmatory analysis (plans/confirmatory_2026-11-03.md; confirmatory papers = md5(str(pmid)) last hex digit not in {0, 1}).
P1 calibrated agreement (expression-matched background, RNAi-hit ceiling), P2 line specificity in genes with CRISPR-RNAi r_x >= 0.3 (vs 0.5 and
the RNAi-top positive control), P3 gene-class profile, P4 expression; secondary strata. Bootstrap 2,000 over papers (P1, P3, P4) or genes (P2).
Per-claim lineage percentiles and expression come from scripts/09_specificity.py, unchanged. Tables go to db/depclaims.duckdb."""
import sys, hashlib, subprocess, numpy as np, pandas as pd, duckdb
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
from rdata.proj import Project
P = Project(__file__); B = 2000; rng = np.random.default_rng(20261103)
dev = lambda p: hashlib.md5(str(p).encode()).hexdigest()[-1] in '01'

def claims(split):
    M = pd.read_parquet('work/reader_mapped_explore.parquet')            # 07_map_reader.py on work/reader_all.jsonl (all papers)
    M = M[M.pmid.map(dev) == (split == 'development')]
    C = M[M.grounded & M.gid.notna() & M.ach.notna() & ((M.growth_effect == 'decreased') | M.apoptosis_increased)].drop_duplicates(['pmid', 'gid', 'ach']).copy()
    C['gid'] = C.gid.astype(int)
    C = C.merge(pd.read_parquet('work/depmap_long.parquet')[['ach', 'gid', 'effect']], on=['ach', 'gid'])
    meta = duckdb.sql("select pmid, any_value(pub_year) yr, any_value(country) country, any_value(retracted) retracted, any_value(journal) journal from read_parquet('work/cand/cand_*.parquet') group by 1").df()
    C = C.merge(meta, on='pmid', how='left').merge(pd.read_parquet('work/icite.parquet')[['pmid', 'relative_citation_ratio']], on='pmid', how='left')
    C['sup'] = C.effect <= -0.5
    C['era'] = pd.cut(C.yr, [0, 2009, 2014, 2019, 2030], labels=['<=2009', '2010-14', '2015-19', '2020+']).astype(str)
    C['meth'] = C.method.str.lower().map(lambda m: 'CRISPR' if m and 'crispr' in m else 'RNAi' if m in ('sirna', 'shrna') else 'other')
    C['cn'] = np.where(C.country == 'China', 'China', 'other')
    f = f'work/claims_{split}.parquet'; C.to_parquet(f, index=False)
    subprocess.run([sys.executable, '-I', 'scripts/09_specificity.py', f], check=True, stdout=subprocess.DEVNULL)
    return pd.read_parquet(f.replace('.parquet', '_spec.parquet'))

# DepMap-only rulers: full 26Q1 CRISPR scaled per model, DEMETER2 scaled the same way, expression deciles within model
ge = pd.read_csv('raw/gene_effect_26q1_chronos.csv', index_col=0); ge = ge[[c for c in ge.columns if '(' in c]]
ge.columns = [int(c.split('(')[-1].rstrip(')')) for c in ge.columns]
ceg = pd.read_csv('raw/CEGv2.txt', sep='\t').ENTREZ_ID.dropna().astype(int); neg = pd.read_csv('raw/NEGv1.txt', sep='\t').ENTREZ_ID.dropna().astype(int)
scale = lambda X: X.sub(X[[g for g in neg if g in X.columns]].median(axis=1), axis=0).div(
    X[[g for g in neg if g in X.columns]].median(axis=1) - X[[g for g in ceg if g in X.columns]].median(axis=1), axis=0)
sc = scale(ge)
frac = (sc <= -0.5).mean(axis=0)
gclass = pd.Series(np.where(frac >= 0.9, 'common', np.where(frac > 0, 'selective', 'never')), index=frac.index)
d2 = pd.read_csv('raw/d2_gene_effect.csv', index_col=0); d2 = d2[[c for c in d2.columns if '(' in c and '&' not in c]]
d2.columns = [int(c.split('(')[-1].rstrip(')')) for c in d2.columns]; d2 = d2[d2.index.str.startswith('ACH-')]; d2 = scale(d2[~d2.index.duplicated()])
cr = sc.stack().rename('c').rename_axis(['ach', 'gid']).reset_index(); rn = d2.stack().rename('r').rename_axis(['ach', 'gid']).reset_index()
BR = cr.merge(rn, on=['ach', 'gid'])
BR['nc'] = BR.gid.map(gclass) != 'common'
h = BR[BR.nc & (BR.r <= -0.5)]; ceiling = (h.c <= -0.5).mean()
never_rnai = (BR[BR.r <= -0.5].gid.map(gclass) == 'never').mean()
rx = BR.groupby('gid').apply(lambda g: g.c.corr(g.r) if len(g) >= 50 else np.nan, include_groups=False).rename('r_x')
Mo = pd.read_csv('raw/Model.csv', usecols=['ModelID', 'OncotreeLineage']).rename(columns={'ModelID': 'ach'})
BL = BR.merge(Mo, on='ach')
X = pd.read_csv('raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv', low_memory=False)
if 'IsDefaultEntryForModel' in X.columns: X = X[X.IsDefaultEntryForModel.astype(str).str.lower().isin(['yes', 'true'])]
X = X.set_index('ModelID' if 'ModelID' in X.columns else X.columns[0]); X = X[[c for c in X.columns if '(' in c]]
sym = {int(c.split('(')[-1].rstrip(')')): c.split(' (')[0] for c in open('raw/gene_effect_26q1_chronos.csv').readline().strip().split(',')[1:] if '(' in c}
ids = [int(c.split('(')[-1].rstrip(')')) for c in X.columns]; dup = pd.Series(ids).duplicated(keep=False).values
X = X[[c for c, i, d in zip(X.columns, ids, dup) if not d or c.split(' (')[0] == sym.get(i, c.split(' (')[0])]]   # same rule as 09
X.columns = [int(c.split('(')[-1].rstrip(')')) for c in X.columns]; X = X.loc[:, ~X.columns.duplicated()]; X = X[~X.index.duplicated()]
EX = X.stack().rename('expr').rename_axis(['ach', 'gid']).reset_index().merge(cr, on=['ach', 'gid'])
EX = EX[EX.gid.map(gclass) != 'common']
EX['dec'] = EX.groupby('ach').expr.transform(lambda v: pd.qcut(v.rank(method='first'), 10, labels=False))
bg = EX.assign(s=EX.c <= -0.5).groupby(['ach', 'dec']).s.mean().rename('bg').reset_index()

def boot(df, unit, fn):
    g = {k: v.index.values for k, v in df.groupby(unit)}; keys = np.array(list(g))
    est = fn(df); bs = [fn(df.loc[np.concatenate([g[k] for k in rng.choice(keys, len(keys))])]) for _ in range(B)]
    return est, *np.nanpercentile(bs, [2.5, 97.5])

def analyse(C, split):
    C = C.copy(); C['gclass'] = C.gid.map(gclass); C['r_x'] = C.gid.map(rx)
    rows = []
    add = lambda k, n, e, lo, hi, ref=np.nan: rows.append(dict(split=split, quantity=k, n=n, est=e, lo=lo, hi=hi, ref=ref))
    # P1
    N = C[C.gclass != 'common'].merge(EX[['ach', 'gid', 'dec']], on=['ach', 'gid'], how='left').merge(bg, on=['ach', 'dec'], how='left').reset_index(drop=True)
    add('P1 non-common support', len(N), *boot(N, 'pmid', lambda d: d.sup.mean()))
    Ne = N.dropna(subset=['bg']).reset_index(drop=True)
    add('P1 expression-matched background', len(Ne), *boot(Ne, 'pmid', lambda d: d.bg.mean()))
    add('P1 ratio to background', len(Ne), *boot(Ne, 'pmid', lambda d: d.sup.mean() / d.bg.mean()), ref=1)
    add('P1 ratio to RNAi-hit ceiling', len(N), *boot(N, 'pmid', lambda d: d.sup.mean() / ceiling), ref=ceiling)
    # P2
    S = C[(C.r_x >= 0.3) & C.pct_CRISPR.notna()].copy()
    pairs = S[['gid', 'OncotreeLineage']].drop_duplicates()
    ctl = {}
    for (g, lin), o in BL.merge(pairs, on=['gid', 'OncotreeLineage']).groupby(['gid', 'OncotreeLineage']):
        if len(o) < 5: continue
        t = o.loc[o.r.idxmin()]; oth = o[o.ach != t.ach].c
        ctl[(g, lin)] = (oth < t.c).mean() + 0.5 * (oth == t.c).mean()
    S['ctl'] = [ctl.get((g, l), np.nan) for g, l in zip(S.gid, S.OncotreeLineage)]
    S = S.reset_index(drop=True)
    add('P2 claimed-line percentile (r_x >= 0.3)', len(S), *boot(S, 'gid', lambda d: d.pct_CRISPR.mean()), ref=0.5)
    Sc = S.dropna(subset=['ctl']).reset_index(drop=True)
    add('P2 positive control percentile', len(Sc), *boot(Sc, 'gid', lambda d: d.ctl.mean()))
    add('P2 claims minus control', len(Sc), *boot(Sc, 'gid', lambda d: (d.pct_CRISPR - d.ctl).mean()), ref=0)
    for t in (0.4, 0.5):
        s = C[(C.r_x >= t) & C.pct_CRISPR.notna()].reset_index(drop=True)
        if len(s) >= 20: add(f'S P2 claimed-line percentile (r_x >= {t})', len(s), *boot(s, 'gid', lambda d: d.pct_CRISPR.mean()), ref=0.5)
    # P3, P4
    Cr = C.reset_index(drop=True)
    for k in ('never', 'selective', 'common'):
        add(f'P3 share {k}', len(Cr), *boot(Cr, 'pmid', lambda d: (d.gclass == k).mean()), ref=(gclass == k).mean())
    rows[-3]['ref2'] = never_rnai
    Ce = C[C.expr.notna()].reset_index(drop=True)
    add('P4 TPM < 1', len(Ce), *boot(Ce, 'pmid', lambda d: d.not_expressed.astype(bool).mean()))
    # secondary strata (support and claimed-line percentile)
    for k in ('era', 'meth', 'cn', 'retracted', 'gclass'):
        for v, d in C.groupby(k, dropna=False):
            d = d.reset_index(drop=True)
            if len(d) >= 20: add(f'S support | {k}={v}', len(d), *boot(d, 'pmid', lambda x: x.sup.mean()))
    rep = C.groupby(['gid', 'ach']).pmid.nunique().rename('n_papers'); C['repeat'] = C.set_index(['gid', 'ach']).index.map(rep) >= 2
    for v, d in C.groupby('repeat'):
        d = d.reset_index(drop=True); add(f'S support | repeated in >= 2 papers={v}', len(d), *boot(d, 'pmid', lambda x: x.sup.mean()))
    for v, d in C.dropna(subset=['relative_citation_ratio']).groupby('sup'):
        d = d.reset_index(drop=True); add(f'S median RCR | supported={v}', len(d), *boot(d, 'pmid', lambda x: x.relative_citation_ratio.median()))
    return pd.DataFrame(rows), C

with P.run(__file__, params=dict(bootstrap=B, seed=20261103, plan='plans/confirmatory_2026-11-03.md')) as run:
    run.input('depmap_26q1_chronos', path='raw/gene_effect_26q1_chronos.csv', note='figshare 31660582, scaled per model to CEGv2 -1 / NEGv1 0')
    run.input('demeter2_combined', path='raw/d2_gene_effect.csv', note='figshare 9170975, same scaling')
    run.input('depmap_24q4_expression', path='raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv', note='figshare 27993248')
    import os; rd = os.environ.get('DEPCLAIMS_READER', 'Qwen3.5-9B reader v2 (e4f3f87)')
    run.input('depclaims_reader', path='work/reader_all.jsonl', note=f'{rd} on 32,939 candidate abstracts')
    out, keep = [], {}
    for split in ('confirmatory', 'development'):
        r, C = analyse(claims(split), split); out.append(r); keep[split] = C
    R = pd.concat(out); R['ceiling'] = ceiling
    run.table('confirm_results', R, note='P1-P4 and secondary estimates with bootstrap 95% CI, by split (plan 2026-11-03)')
    run.parquet('claims_confirmatory', keep['confirmatory'].drop(columns=['gene_str', 'line_str'], errors='ignore'), note='confirmatory claim set with rulers')
    run.parquet('claims_development', keep['development'].drop(columns=['gene_str', 'line_str'], errors='ignore'), note='development claim set with rulers')
    pd.set_option('display.width', 200); pd.set_option('display.max_rows', 200)
    print(f'RNAi-hit ceiling {ceiling:.3f} | never share of RNAi hits {never_rnai:.3f} | genes with r_x {rx.notna().sum()}')
    print(R.round(3).to_string(index=False))

"""Pre-specified secondary analysis (plans/confirmatory_2026-11-03.md): independent repeat claims. Confirmatory claims only.
A (gene, model) pair is an independent repeat when at least two of its papers share no author (author = accent-folded last name + first
initial, or the same ORCID); a shared-author repeat when every pair of its papers shares an author; single when one paper makes it.
Support (scaled Chronos <= -0.5) per class with 2,000 bootstrap resamples over papers (seed 20261103, as in 12_confirm.py)."""
import sys, unicodedata, itertools, numpy as np, pandas as pd
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
from rdata.proj import Project
P = Project(__file__); B = 2000; rng = np.random.default_rng(20261103)
fold = lambda s: ''.join(c for c in unicodedata.normalize('NFKD', s or '') if not unicodedata.combining(c)).lower().replace('-', ' ').strip()
with P.run(__file__, params=dict(B=B, seed=20261103, author_key='last name + first initial or ORCID')) as run:
    run.input('claims_confirmatory', path='db/parquet/claims_confirmatory.parquet', note='confirmatory claim set from 12_confirm.py')
    C = pd.read_parquet('db/parquet/claims_confirmatory.parquet')
    run.input('pubmed_2026_authors', path='work/dcauth/out/authors.parquet', note='author lists from the PubMed 2026 baseline (scripts/m5/14_authors.py)')
    A = pd.read_parquet('work/dcauth/out/authors.parquet')
    A['key'] = A['last'].map(fold) + '|' + A.fore.map(lambda f: fold(f)[:1])
    au = A.groupby('pmid').apply(lambda d: set(d.key) | {'orcid:' + o for o in d.orcid.dropna()}, include_groups=False).to_dict()
    rows = []
    for (g, a), d in C.groupby(['gid', 'ach']):
        ps = sorted(d.pmid.unique())
        if len(ps) == 1: cls = 'single'
        elif any(p not in au or not au[p] for p in ps): cls = 'undetermined'
        else: cls = 'independent repeat' if any(not (au[x] & au[y]) for x, y in itertools.combinations(ps, 2)) else 'shared-author repeat'
        rows.append((g, a, len(ps), cls))
    R = pd.DataFrame(rows, columns=['gid', 'ach', 'n_papers', 'repeat_class'])
    C = C.merge(R, on=['gid', 'ach'])
    def boot(df, fn):
        u = df.pmid.unique(); idx = {p: i for p, i in df.groupby('pmid').indices.items()}
        est = fn(df); bs = []
        for _ in range(B):
            s = rng.choice(u, len(u)); bs.append(fn(df.iloc[np.concatenate([idx[p] for p in s])]))
        return est, *np.nanpercentile(bs, [2.5, 97.5])
    out = []
    for k, d in C.groupby('repeat_class'):
        d = d.reset_index(drop=True); e, lo, hi = boot(d, lambda x: x.sup.mean())
        out.append(dict(analysis='S support | repeat class', level=k, n_claims=len(d), n_pairs=d[['gid', 'ach']].drop_duplicates().shape[0],
                        n_papers=d.pmid.nunique(), est=e, lo=lo, hi=hi))
    # same contrast outside common-essential genes, where support is not near-automatic
    for k, d in C[C.gclass != 'common'].groupby('repeat_class'):
        d = d.reset_index(drop=True); e, lo, hi = boot(d, lambda x: x.sup.mean())
        out.append(dict(analysis='S support | repeat class, non-common genes', level=k, n_claims=len(d),
                        n_pairs=d[['gid', 'ach']].drop_duplicates().shape[0], n_papers=d.pmid.nunique(), est=e, lo=lo, hi=hi))
    O = pd.DataFrame(out); print(O.round(3).to_string(index=False))
    print('papers with author lists', round(np.mean([p in au for p in C.pmid.unique()]), 4))
    run.table('independent_repeats', O, note='support by repeat class (pre-specified secondary)'); run.parquet('claims_confirmatory_repeat_class', C[['pmid', 'gid', 'ach', 'n_papers', 'repeat_class']], note='repeat class per confirmatory claim')

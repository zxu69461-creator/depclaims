"""Claim QC v3 scoring (validation/claim_qc_protocol_v3.md). Merges reader B parts, reports per-question agreement (raw and Cohen kappa),
claim correctness (Q1, Q2, Q3, Q6, Q7 all yes) per reader and per stratum with Wilson 95% CI, and writes the disagreement list for a third reading."""
import numpy as np, pandas as pd
from statsmodels.stats.proportion import proportion_confint
from sklearn.metrics import cohen_kappa_score
V = 'validation/claim_qc_v3_'
rd = lambda f: pd.read_csv(f, sep='\t', dtype=str).assign(item=lambda d: d.item.astype(int))
A = rd(V + 'reader_A.tsv')
Bp = pd.concat([rd(V + f'reader_B{s}.tsv') for s in ('', '_part2', '_part3')]).drop_duplicates('item').sort_values('item')
Bp.to_csv(V + 'reader_B_merged.tsv', sep='\t', index=False)
assert sorted(A.item) == sorted(Bp.item) == list(range(1, 201))
K = pd.read_csv(V + 'key.tsv', sep='\t')
Q = ['Q1', 'Q2', 'Q3', 'Q6', 'Q7']
D = A.merge(Bp, on='item', suffixes=('_A', '_B')).merge(K[['item', 'stratum', 'era']], on='item')
for q in Q + ['Q4', 'Q5']:
    print(f'{q}: agreement {(D[q + "_A"] == D[q + "_B"]).mean():.3f} kappa {cohen_kappa_score(D[q + "_A"], D[q + "_B"]):.2f}')
for r in 'AB': D[f'ok_{r}'] = (D[[f'{q}_{r}' for q in Q]] == 'yes').all(axis=1)
D['ok_both'] = D.ok_A & D.ok_B; D['disagree'] = D.ok_A != D.ok_B
def w(x): k, n = int(x.sum()), len(x); lo, hi = proportion_confint(k, n, method='wilson'); return f'{k}/{n} = {k / n:.3f} ({lo:.3f}-{hi:.3f})'
for s, d in D.groupby('stratum'):
    print(f'{s}: correct A {w(d.ok_A)} | B {w(d.ok_B)} | both {w(d.ok_both)} | disagree on correctness {int(d.disagree.sum())}')
print('main by era (both readers correct):', {e: w(d.ok_both) for e, d in D[D.stratum == 'main'].groupby('era')})
fail = D[~D.ok_both]
for q in Q:
    print(f'{q} failed by >= 1 reader:', {s: int(((d[q + "_A"] != "yes") | (d[q + "_B"] != "yes")).sum()) for s, d in fail.groupby('stratum')})
D[D.disagree][['item'] + [f'{q}_{r}' for r in 'AB' for q in Q] + ['note_A', 'note_B']].to_csv(V + 'disagreements.tsv', sep='\t', index=False)
D.to_csv(V + 'scored_two_readers.tsv', sep='\t', index=False)
print('disagreements written:', int(D.disagree.sum()))

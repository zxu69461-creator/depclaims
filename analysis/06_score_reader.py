"""Score reader-model output against blinded human-protocol readers. Usage: 06_score_reader.py OUT_JSONL v1|v2
Model claims are mapped to ids with scripts/07_map_reader.py. A QC triple is a true claim when both readers answered Q1, Q2, Q3 = yes (v1: limited to triples whose line is named in the abstract).
The model calls the triple a claim when it returns that gene and cell line with growth_effect 'decreased' or apoptosis_increased true."""
import sys, json, pandas as pd
out, ver = sys.argv[1:3]
R = {}
for l in open(out):
    d = json.loads(l); R[d['pmid']] = d['claims']
A = pd.read_csv(f'validation/claim_qc_{ver}_reader_A.tsv', sep='\t'); B = pd.read_csv(f'validation/claim_qc_{ver}_reader_B.tsv', sep='\t')
K = pd.read_csv(f'validation/claim_qc_{ver}_key.tsv', sep='\t')
P = json.load(open(f'validation/claim_qc_{ver}_packet.json')); P = pd.DataFrame(P)[['item', 'gene', 'cell_line']]
m = A.merge(B, on='item', suffixes=('_a', '_b')).merge(K, on='item').merge(P, on='item')
m = m[m.pmid.isin(R.keys())]
ok = lambda r, s: r[f'Q1_{s}'] == 'yes' and r[f'Q2_{s}'] == 'yes' and r[f'Q3_{s}'] == 'yes'
m['truth'] = m.apply(lambda r: ok(r, 'a') and ok(r, 'b'), axis=1); m['disagree'] = m.apply(lambda r: ok(r, 'a') != ok(r, 'b'), axis=1)
import subprocess; mp = out.replace('.jsonl', '_mapped.parquet')
subprocess.run([sys.executable, '-I', 'scripts/07_map_reader.py', out, mp], check=True)
M = pd.read_parquet(mp); M = M[((M.growth_effect == 'decreased') | M.apoptosis_increased) & (M.grounded if 'grounded' in M and len(sys.argv) < 4 else True)]
pos = set(zip(M.pmid, M.gid, M.ach))
m['pred'] = [(p, float(g), a) in {(x, float(y) if pd.notna(y) else None, z) for x, y, z in pos} for p, g, a in zip(m.pmid, m.gid, m.ach)]
e = m[~m.disagree]
tp = (e.pred & e.truth).sum(); fp = (e.pred & ~e.truth).sum(); fn = (~e.pred & e.truth).sum(); tn = (~e.pred & ~e.truth).sum()
print(f'{ver}: items {len(m)} (reader disagreement excluded {m.disagree.sum()}); TP {tp} FP {fp} FN {fn} TN {tn}; '
      f'precision {tp / max(tp + fp, 1):.3f} recall {tp / max(tp + fn, 1):.3f}; unparsed {sum(v is None for v in R.values())}')
m[['item', 'pmid', 'gene', 'cell_line', 'truth', 'disagree', 'pred']].to_csv(out.replace('.jsonl', f'_score_{ver}.tsv'), sep='\t', index=False)

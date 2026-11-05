"""Additional analysis (2026-11-05, after the PNAS pre-submission review; not prespecified): two sources of claims that screens do not
support. Confirmatory claims are split into (i) genes not expressed in the named line (TPM < 1), (ii) expressed genes whose product is
secreted (UniProt reviewed human entries with keyword KW-0964 'Secreted' and without the keywords Membrane or Transmembrane, so that
receptors with a secreted isoform such as EGFR are not counted; mapped to Entrez through HGNC), whose growth effects may
act through paracrine or microenvironmental routes that a monolayer screen does not measure, and (iii) all other expressed genes. Per class:
claims, share of all claims and of unsupported claims, CRISPR support (bootstrap over genes), and the most claimed genes.
Writes db table secreted_additional. Usage: 38_secreted.py"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
from rdata.proj import Project
P = Project(__file__); rng = np.random.default_rng(20261105); B = 2000
with P.run(__file__, params=dict(B=B, seed=20261105, keyword='KW-0964 without Membrane/Transmembrane')) as run:
    run.input('uniprot_secreted', path='raw/uniprot_secreted_human_2026-11-05.tsv', note='UniProt reviewed human, keyword Secreted, 2,146 entries, 2026-11-05')
    U = pd.read_csv('raw/uniprot_secreted_human_2026-11-05.tsv', sep='\t')
    U = U[~U.Keywords.fillna('').str.contains(r'\bMembrane\b|Transmembrane', regex=True)]
    H = pd.read_csv('/Volumes/ResearchData/Claude/ref/raw/hgnc_complete_set.txt', sep='\t', usecols=['symbol', 'entrez_id', 'uniprot_ids'],
                    low_memory=False).dropna(subset=['entrez_id'])
    acc2id = {a: int(r.entrez_id) for r in H.itertuples() if isinstance(r.uniprot_ids, str) for a in r.uniprot_ids.split('|')}
    sym2id = {s: int(e) for s, e in zip(H.symbol, H.entrez_id)}
    sec = {acc2id.get(a) or sym2id.get(g) for a, g in zip(U['Entry'], U['Gene Names (primary)'])} - {None}
    A = pd.read_csv('release/claims_scored.tsv.gz', sep='\t')
    C = A[(A.split == 'confirmatory') & A.expression_log2_tpm1.notna()].reset_index(drop=True)
    C['cls'] = np.where(C.not_expressed.astype(bool), 'not expressed (TPM < 1)', np.where(C.entrez_id.isin(sec), 'expressed, secreted product', 'expressed, other'))
    C['sup'] = C.crispr_supported.astype(bool)
    def gboot(d):
        g = {k: v.index.values for k, v in d.groupby('entrez_id')}; kk = np.array(list(g))
        res = [d.loc[np.concatenate([g[k] for k in rng.choice(kk, len(kk))])].sup.mean() for _ in range(B)]
        return d.sup.mean(), *np.percentile(res, [2.5, 97.5])
    out = []
    uns = (~C.sup).sum()
    for k, d in C.groupby('cls'):
        est, lo, hi = gboot(d.reset_index(drop=True))
        top = d.symbol.value_counts().head(6).index.tolist()
        out.append(dict(cls=k, claims=len(d), genes=d.entrez_id.nunique(), share_claims=len(d) / len(C), share_unsupported=(~d.sup).sum() / uns,
                        support=est, lo=lo, hi=hi, top_genes=', '.join(top)))
    O = pd.DataFrame(out); print(O.round(4).to_string(index=False)); print('secreted genes in the list', len(sec))
    run.table('secreted_additional', O, note='unsupported claims by class: unexpressed, secreted, other expressed (additional analysis 2026-11-05)')

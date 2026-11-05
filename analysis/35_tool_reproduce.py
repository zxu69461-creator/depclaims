"""Check that the depclaims package reproduces the scored values of every claim in release/claims_scored.tsv.gz (CRISPR and RNAi effects,
support, gene class, within-lineage ranks, dependent fractions), at the precision written to the release (3-4 decimals). Writes work/tool_check/mismatches.tsv."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, 'tool/src')
from depclaims.core import Atlas
A = pd.read_csv('release/claims_scored.tsv.gz', sep='\t'); T = Atlas('work/tool_cache')
cmp = [('crispr_effect_scaled', 'crispr_effect_scaled'), ('crispr_supported', 'needed_in_line'), ('gene_class', 'gene_class'),
       ('crispr_rank_in_lineage', 'crispr_rank_in_lineage'), ('rnai_effect_scaled', 'rnai_effect_scaled'), ('rnai_rank_in_lineage', 'rnai_rank_in_lineage'),
       ('gene_dependent_fraction_all', 'dependent_fraction_all'), ('gene_dependent_fraction_lineage', 'dependent_fraction_lineage')]
cache, bad = {}, []
for r in A.itertuples():
    k = (int(r.entrez_id), r.depmap_model)
    if k not in cache: cache[k] = T.score(*k).__dict__
    s = cache[k]
    for a, b in cmp:
        x, y = getattr(r, a), s[b]
        same = (pd.isna(x) and y is None) or (y is not None and not pd.isna(x) and (x == y if isinstance(y, (str, bool)) else np.isclose(float(x), float(y), rtol=0, atol=5.01e-4)))
        if not same: bad.append(dict(pmid=r.pmid, entrez_id=k[0], model=k[1], field=a, atlas=x, package=y))
B = pd.DataFrame(bad); B.to_csv('work/tool_check/mismatches.tsv', sep='\t', index=False)
print('claims', len(A), 'pairs', len(cache), 'mismatching claims', B.pmid.nunique() if len(B) else 0)
if len(B): print(B.field.value_counts())

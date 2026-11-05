"""Between-reader comparison (notes/analysis_plan_review2_2026-11-04.md, section E): Claude Sonnet 5.5 (primary) vs the local Qwen3.5-9B reader,
same prompt, same mapping and filters. Claim-level overlap by split, support of shared and reader-specific claims, and the confirmatory
estimates (confirm_results) under each reader side by side."""
import sys, pandas as pd
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
from rdata.proj import Project
P = Project(__file__)
k = ['pmid', 'gid', 'ach']
with P.run(__file__) as run:
    out = []
    for split in ('confirmatory', 'development'):
        S = pd.read_parquet(f'db/parquet/claims_{split}.parquet'); Q = pd.read_parquet(f'work/qwen/db_parquet/claims_{split}.parquet')
        m = S[k + ['sup']].merge(Q[k + ['sup']], on=k, how='outer', suffixes=('_s', '_q'), indicator=True)
        both, so, qo = (m._merge == 'both'), (m._merge == 'left_only'), (m._merge == 'right_only')
        for lab, v in (('claims, Sonnet', len(S)), ('claims, Qwen', len(Q)), ('shared', both.sum()), ('Sonnet only', so.sum()), ('Qwen only', qo.sum()),
                       ('Jaccard', both.sum() / len(m)), ('share of Sonnet claims also read by Qwen', both.sum() / len(S)),
                       ('support, shared', m.loc[both, 'sup_s'].mean()), ('support, Sonnet only', m.loc[so, 'sup_s'].mean()),
                       ('support, Qwen only', m.loc[qo, 'sup_q'].mean())):
            out.append(dict(split=split, quantity=lab, value=float(v)))
    O = pd.DataFrame(out); print(O.round(4).to_string(index=False))
    Rs = pd.read_parquet('work/confirm_results_sonnet.parquet')        # exported right after the re-run of the confirmatory script
    Rq = pd.read_parquet('work/qwen/db_tables/confirm_results.parquet')
    key = ['split', 'quantity']
    J = Rs.merge(Rq, on=key, suffixes=('_sonnet', '_qwen'), how='outer'); print(J.round(3).to_string(index=False))
    run.table('reader_agreement', O, note='claim-level agreement between Claude Sonnet 5.5 and Qwen3.5-9B readers')
    run.table('confirm_results_by_reader', J, note='confirmatory estimates under each reader')

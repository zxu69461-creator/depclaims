"""Claim QC v3 sample (reader precision on production claims; readers never see DepMap values). Usage: 11_qc_v3_sample.py
Main sample: 150 model-positive claims (frozen reader v2 + mapping 07, grounded, growth decreased or apoptosis increased, gene and model in
DepMap 26Q1), all papers, stratified by era (<=2009, 2010-14, 2015-19, 2020+; 37-38 each). Selection uses reader output and DepMap membership only.
Over-sample: 50 development claims whose gene is not expressed in the claimed model (TPM < 1), to split mapping errors from genuine cases;
drawn from development papers only so no confirmatory expression value is looked at. Items are shuffled together; the key holds the stratum."""
import json, numpy as np, pandas as pd, duckdb
rng = np.random.default_rng(20261104)
M = pd.read_parquet('work/reader_mapped_explore.parquet')
C = M[M.grounded & M.gid.notna() & M.ach.notna() & ((M.growth_effect == 'decreased') | M.apoptosis_increased)].drop_duplicates(['pmid', 'gid', 'ach']).copy()
C['gid'] = C.gid.astype(int)
L = duckdb.sql("select distinct ach, gid from read_parquet('work/depmap_long.parquet')").df()
C = C.merge(L, on=['ach', 'gid'])
meta = duckdb.sql("select pmid, any_value(pub_year) yr, any_value(title) title, any_value(abstract) abstract from read_parquet('work/cand/cand_*.parquet') group by 1").df()
C = C.merge(meta, on='pmid', how='left')
C['era'] = pd.cut(C.yr, [0, 2009, 2014, 2019, 2030], labels=['<=2009', '2010-14', '2015-19', '2020+']).astype(str)
main = pd.concat([g.sample(n, random_state=int(rng.integers(10**9))) for (e, g), n in zip(C.groupby('era'), [38, 37, 38, 37])])
main['stratum'] = 'main'
D = pd.read_parquet('work/explore_dev_claims_spec.parquet')
D = D[D.not_expressed.fillna(False).astype(bool)][['pmid', 'gid', 'ach']].merge(C, on=['pmid', 'gid', 'ach'])
D = D.merge(main[['pmid', 'gid', 'ach']], how='left', indicator=True).query('_merge == "left_only"').drop(columns='_merge')
over = D.sample(min(50, len(D)), random_state=int(rng.integers(10**9))); over['stratum'] = 'not_expressed'
S = pd.concat([main, over]).sample(frac=1, random_state=int(rng.integers(10**9))).reset_index(drop=True)
S['item'] = np.arange(1, len(S) + 1)
H = pd.read_csv('/Volumes/ResearchData/Claude/ref/raw/hgnc_complete_set.txt', sep='\t', usecols=['symbol', 'name', 'alias_symbol', 'entrez_id'], low_memory=False)
H = H.dropna(subset=['entrez_id']).assign(gid=lambda h: h.entrez_id.astype(int)).drop_duplicates('gid')
S = S.merge(H[['gid', 'symbol', 'name', 'alias_symbol']], on='gid', how='left')
mod = pd.read_csv('raw/Model.csv', usecols=['ModelID', 'CellLineName', 'OncotreePrimaryDisease']).rename(columns={'ModelID': 'ach'})
S = S.merge(mod, on='ach', how='left').sort_values('item')
pk = [dict(item=int(r.item), pmid=int(r.pmid), gene_symbol=r.symbol, gene_name=r.name, gene_aliases=r.alias_symbol if isinstance(r.alias_symbol, str) else '',
           cell_line=r.CellLineName, cell_line_disease=r.OncotreePrimaryDisease, title=r.title, abstract=r.abstract) for r in S.itertuples()]
json.dump(pk, open('validation/claim_qc_v3_packet.json', 'w'), indent=1, ensure_ascii=False)
S[['item', 'pmid', 'gid', 'ach', 'stratum', 'era', 'method', 'gene_str', 'line_str']].to_csv('validation/claim_qc_v3_key.tsv', sep='\t', index=False)
print('production claims', len(C), '| by era', C.era.value_counts().to_dict(), '| over-sample pool', len(D), '| packet', len(S), S.stratum.value_counts().to_dict())

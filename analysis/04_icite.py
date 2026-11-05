"""iCite citation metrics for every claim paper (citation count, RCR, clinical citations). Writes work/icite.parquet."""
import json, subprocess, duckdb, pandas as pd, concurrent.futures as cf, sys
sys.path.insert(0, '/Volumes/ResearchData/Claude/netroute')
p = duckdb.sql("select distinct pmid from read_parquet('work/triples.parquet')").df().pmid.tolist()
def get(ch):
    u = 'https://icite.od.nih.gov/api/pubs?fl=pmid,year,citation_count,relative_citation_ratio,is_clinical,cited_by_clin&pmids=' + ','.join(map(str, ch))
    for k in range(6):
        r = subprocess.run(['ncurl', '-s', '-m', '120', u], capture_output=True, text=True)
        try: return json.loads(r.stdout)['data']
        except Exception: __import__("time").sleep(5 * (k + 1))
    raise RuntimeError('iCite failed for chunk starting ' + str(ch[0]))
chunks = [p[i:i + 200] for i in range(0, len(p), 200)]
with cf.ThreadPoolExecutor(3) as ex: rows = [x for d in ex.map(get, chunks) for x in d]
df = pd.DataFrame(rows); df['n_clin_cites'] = df.cited_by_clin.apply(lambda v: len(v) if isinstance(v, list) else 0)
df = df.drop(columns=['_id', 'cited_by_clin'], errors='ignore'); df.to_parquet('work/icite.parquet', index=False)
assert len(df) > 0.95 * len(p), f'only {len(df)} of {len(p)} pmids returned'
print('pmids', len(p), 'returned', len(df))

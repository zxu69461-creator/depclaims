"""Map reader-model claims to (pmid, Entrez gene id, DepMap model) using each paper's own PubTator3 mention forms, the DepMap symbol and the
Cellosaurus name. Usage: 07_map_reader.py IN_JSONL OUT_PARQUET. Rows: pmid, gid, ach, method, growth_effect, apoptosis_increased, gene_str, line_str."""
import sys, json, re, duckdb, pandas as pd
inp, outp = sys.argv[1:3]
norm = lambda s: re.sub(r'[^A-Z0-9]', '', str(s).upper())
R = [json.loads(l) for l in open(inp)]
pm = sorted({d['pmid'] for d in R}); con = duckdb.connect(); con.register('p', pd.DataFrame({'pmid': pm}))
gm = con.execute("SELECT g.pmid, try_cast(g.gene_id AS BIGINT) gid, g.mention FROM read_parquet('work/claims/gene_ann.parquet') g JOIN p USING (pmid)").df()
cm = con.execute("SELECT c.pmid, c.ach, c.mention FROM read_parquet('work/claims/cell_ann.parquet') c JOIN p USING (pmid) WHERE c.ach IS NOT NULL").df()
sym = {int(c.split('(')[-1].rstrip(')')): c.split(' (')[0] for c in open('raw/gene_effect_26q1_chronos.csv').readline().strip().split(',')[1:] if '(' in c}
cn = pd.read_csv('raw/cvcl2ach.tsv', sep='\t')
G, C = {}, {}
for r in gm.itertuples():
    if pd.isna(r.gid) or int(r.gid) not in sym: continue
    for f in [sym[int(r.gid)]] + str(r.mention).split('|'): G.setdefault(r.pmid, {}).setdefault(norm(f), int(r.gid))
for r in cm.itertuples():
    for f in str(r.mention).split('|') + cn[cn.ach == r.ach].cell_name.tolist(): C.setdefault(r.pmid, {}).setdefault(norm(f), r.ach)
ab = con.execute("SELECT c.pmid, any_value(c.title) || ' ' || coalesce(any_value(c.abstract), '') tx FROM read_parquet('work/cand/cand_*.parquet') c JOIN p USING (pmid) GROUP BY 1").df().set_index('pmid').tx
sq = lambda s: re.sub(r'\s+', ' ', str(s or '')).strip().lower()
GW = re.compile(r'prolifer|growth|viab|colon|clonogen|tumou?r|xenograft|apoptos|surviv|cell cycle|arrest|cell death|g[012]/?[ms]? ?(phase)?')
def grounded(pmid, c, ach):
    tx = sq(ab.get(pmid, '')); gq, lq = sq(c.get('growth_quote')), sq(c.get('line_quote'))
    if len(gq) < 15 or gq[:200] not in tx or not GW.search(gq): return False
    forms = {f for f, a in C.get(pmid, {}).items() if a == ach} if ach else set()
    ctx = [norm(gq)] + ([norm(lq)] if len(lq) >= 30 and lq[:200] in tx else [])   # line named in the growth sentence or in a full quoted sentence
    return any(f and f in x for f in forms for x in ctx)
rows = []
for d in R:
    for c in d['claims'] or []:
        if not isinstance(c, dict): continue
        g = G.get(d['pmid'], {}).get(norm(c.get('gene'))); a = C.get(d['pmid'], {}).get(norm(c.get('cell_line')))
        rows.append((d['pmid'], g, a, c.get('method'), c.get('growth_effect'), bool(c.get('apoptosis_increased')), c.get('gene'), c.get('cell_line'), grounded(d['pmid'], c, a)))
df = pd.DataFrame(rows, columns=['pmid', 'gid', 'ach', 'method', 'growth_effect', 'apoptosis_increased', 'gene_str', 'line_str', 'grounded'])
df.to_parquet(outp, index=False)
print('papers', len(R), 'unparsed', sum(d['claims'] is None for d in R), '| claims', len(df), 'gene mapped', df.gid.notna().mean().round(3), 'line mapped', df.ach.notna().mean().round(3), 'grounded', df.grounded.mean().round(3))

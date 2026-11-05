"""Claude reader for gene-perturbation growth claims: the prompt of scripts/m5/05_llm_read.py, sent through the Claude Code CLI (no tools,
no settings, no MCP). Usage: 29_claude_read.py MODEL ITEMS_PARQUET OUT_JSONL [WORKERS] [EFFORT]
Output: one JSON line per pmid with the parsed claim list (or the raw text when parsing fails), as the local reader writes it. Resumes."""
import sys, os, json, re, subprocess, threading, time, importlib.util
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
model, items_p, out_p = sys.argv[1:4]
W = int(sys.argv[4]) if len(sys.argv) > 4 else 6; EFF = sys.argv[5] if len(sys.argv) > 5 else 'medium'
src = open(os.path.join(os.path.dirname(__file__), 'm5', '05_llm_read.py')).read()
SYS = re.search(r'SYS = \((.*?)\)\n', src, re.S).group(1); SYS = eval('(' + SYS + ')')
TPL = eval(re.search(r'TPL = (""".*?""")', src, re.S).group(1))
done = set()
if os.path.exists(out_p):
    for l in open(out_p):
        try: done.add(json.loads(l)['pmid'])
        except Exception: pass
it = pd.read_parquet(items_p); it = it[~it.pmid.isin(done)]
CWD = os.environ.get('READER_CWD', '/private/tmp')
CMD = ['claude', '-p', '--model', model, '--system-prompt', SYS, '--tools', '', '--strict-mcp-config', '--setting-sources', '',
       '--exclude-dynamic-system-prompt-sections', '--no-session-persistence', '--output-format', 'json', '--effort', EFF]
def parse(txt):
    # first JSON object carrying 'claims' (the model sometimes adds a code fence or a sentence after the object)
    dec = json.JSONDecoder(); txt = txt or ''
    for m in re.finditer(r'\{', txt):
        try: o, _ = dec.raw_decode(txt, m.start())
        except Exception: continue
        if isinstance(o, dict) and 'claims' in o: return o['claims']
    return None
lock = threading.Lock(); n = [0]; cost = [0.0]; t0 = time.time()
def read(r):
    prompt = TPL.format(title=r.title, abstract=r.abstract or '', genes=', '.join(r.genes), lines=', '.join(r.lines))
    for attempt in range(4):
        p = subprocess.run(CMD, input=prompt, capture_output=True, text=True, cwd=CWD, timeout=600)
        t = p.stdout; i = t.find('{"duration')
        try: d = json.loads(t[i:]) if i >= 0 else {}
        except Exception: d = {}
        if d and not d.get('is_error'): break
        time.sleep(20 * (attempt + 1))
    txt = d.get('result', '') if d else ''
    claims = parse(txt)
    with lock:
        out.write(json.dumps({'pmid': int(r.pmid), 'claims': claims, 'raw': None if claims is not None else (txt or t)[:800]}) + '\n'); out.flush()
        n[0] += 1; cost[0] += d.get('total_cost_usd', 0) or 0
        if n[0] % 25 == 0: print(n[0], f'{(time.time() - t0) / n[0]:.2f}s/item', f'notional ${cost[0]:.2f}', flush=True)
with open(out_p, 'a') as out, ThreadPoolExecutor(W) as ex:
    list(ex.map(read, it.itertuples()))
print('done', len(it), f'{time.time() - t0:.0f}s', f'notional ${cost[0]:.2f}')

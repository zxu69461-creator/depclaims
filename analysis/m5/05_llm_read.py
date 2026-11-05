"""Local reader model (MLX on M5) for gene-perturbation growth claims. Usage: 05_llm_read.py MODEL ITEMS_PARQUET OUT_JSONL
ITEMS: pmid, title, abstract, genes (list of symbols), lines (list of cell-line names). One prompt per abstract; the model may only use the
listed genes and lines. Output: one JSON line per pmid with the parsed claim list (or the raw text when parsing fails). Resumes."""
import sys, os, json, re, time
import pandas as pd
from mlx_lm import load, generate, batch_generate
from mlx_lm.sample_utils import make_sampler
model_id, items_p, out_p = sys.argv[1:4]; BS = int(sys.argv[4]) if len(sys.argv) > 4 else 16
SYS = ("You extract experimental claims from biomedical abstracts. Report only what the authors themselves did and found, as stated in the "
       "title or abstract. Answer with JSON only.")
TPL = """Title: {title}
Abstract: {abstract}

Candidate genes: {genes}
Candidate cell lines: {lines}

List every claim in which the authors experimentally REDUCED or ABOLISHED one candidate gene (siRNA, shRNA, knockdown, CRISPR/knockout, silencing,
antisense, depletion) in one candidate cell line and report an effect on growth-type phenotypes (proliferation, viability, growth, colony formation,
survival, tumorigenicity/xenograft growth, apoptosis, cell-cycle arrest). Do not count genes that were only measured or overexpressed, genes reduced
indirectly through another molecule (e.g. a microRNA, lncRNA or drug acting on them), or cell lines not stated to carry that experiment. If the
abstract states the experiment in a short list of cell lines, give one claim per line. Count a claim only when the abstract reports the growth effect of reducing that gene itself; not when the reduction only serves to reverse or modify the effect of a drug, another gene or another treatment.
For each claim copy, word for word from the title or abstract, the sentence that states the growth-type effect of reducing the gene
("growth_quote") and the sentence that names the cell line as used in that experiment ("line_quote"; may be the same sentence).
Return JSON: {{"claims": [{{"gene": "...", "cell_line": "...", "method": "siRNA|shRNA|CRISPR|antisense|other", "growth_effect": "decreased|increased|no_change",
"apoptosis_increased": true|false, "growth_quote": "...", "line_quote": "..."}}]}}. Return {{"claims": []}} if there is none."""
done = set()
if os.path.exists(out_p):
    for l in open(out_p):
        try: done.add(json.loads(l)['pmid'])
        except Exception: pass
it = pd.read_parquet(items_p); it = it[~it.pmid.isin(done)]
model, tok = load(model_id); samp = make_sampler(temp=0.0)
def prompt_of(r):
    msg = [{'role': 'system', 'content': SYS}, {'role': 'user', 'content': TPL.format(title=r.title, abstract=r.abstract or '',
           genes=', '.join(r.genes), lines=', '.join(r.lines))}]
    try: return tok.apply_chat_template(msg, add_generation_prompt=True, tokenize=False, enable_thinking=False)
    except TypeError: return tok.apply_chat_template(msg, add_generation_prompt=True, tokenize=False)
def parse(txt):
    txt = re.sub(r'(?s)<think>.*?</think>', '', txt).strip(); m = re.search(r'(?s)\{.*\}', txt)
    try: return (json.loads(m.group(0))['claims'] if m else None), txt
    except Exception: return None, txt
t0 = time.time(); rows = list(it.itertuples())
with open(out_p, 'a') as out:
    for k in range(0, len(rows), BS):
        chunk = rows[k:k + BS]
        toks = [tok.encode(prompt_of(r), add_special_tokens=False) for r in chunk]
        res = batch_generate(model, tok, toks, max_tokens=900, sampler=samp).texts if BS > 1 else [generate(model, tok, prompt=prompt_of(chunk[0]), max_tokens=900, sampler=samp)]
        for r, txt in zip(chunk, res):
            claims, txt = parse(txt)
            out.write(json.dumps({'pmid': int(r.pmid), 'claims': claims, 'raw': None if claims is not None else txt[:800]}) + '\n')
        out.flush(); print(k + len(chunk), f'{(time.time() - t0) / (k + len(chunk)):.2f}s/item', flush=True)
print('done', len(it), f'{time.time() - t0:.0f}s')

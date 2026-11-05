"""Additional analysis D of notes/analysis_plan_review2_2026-11-04.md (not prespecified): do open-access full texts describe the experiment the
abstract claim states? Confirmatory papers with PMC full text (Europe PMC JATS); 50 papers drawn from claims on genes not expressed in the named
line and 50 from all other claims (one claim per paper, seed 20261104, drawn in random order until 50 full texts per stratum are retrieved).
Per claim, text checks: (1) gene named in a sentence with a loss-of-function reagent; (2) named line in the Methods; (3) growth or viability
assay named; (gene forms: HGNC symbol, aliases, previous symbols and the paper's own annotated mentions, matched ignoring hyphens, and ignoring case for forms with a digit or of five or more characters) (4) knockdown or knockout verification (immunoblot, western blot, qPCR) named in a sentence with the gene; (5) gene, line and
reagent in one sentence or figure legend. Writes validation/fulltext_items.tsv (passages for the human questionnaire) and a summary table.
Usage: 32_fulltext.py [OUT_SUFFIX]"""
import sys, os, re, json, subprocess, html, numpy as np, pandas as pd
import xml.etree.ElementTree as ET
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
from rdata.proj import Project
P = Project(__file__); rng = np.random.default_rng(20261104)
suf = sys.argv[1] if len(sys.argv) > 1 else ''
XD = 'work/fulltext/xml'; os.makedirs(XD, exist_ok=True)
pmc = json.load(open('work/fulltext/pmid2pmc.json')); pmc = {int(k): v for k, v in pmc.items()}
REAG = r'siRNA|shRNA|sgRNA|CRISPR|knock[- ]?down|knock[- ]?out|silenc|antisense|depleti|deplet|si-|sh-|RNAi|Cas9'
ASSAY = re.compile(r'CCK-?8|\bMTT\b|\bMTS\b|CellTiter|colony[- ]formation|clonogenic|\bEdU\b|BrdU|proliferation assay|cell counting|xenograft|soft agar|'
                   r'trypan|growth curve|SRB|real-time cell|IncuCyte', re.I)
VER = r'western blot|immunoblot|qRT-PCR|RT-qPCR|real-time PCR|qPCR|RT-PCR|immunofluorescen'
H = pd.read_csv('/Volumes/ResearchData/Claude/ref/raw/hgnc_complete_set.txt', sep='\t', usecols=['symbol', 'alias_symbol', 'prev_symbol', 'entrez_id'],
                low_memory=False).dropna(subset=['entrez_id']).assign(entrez_id=lambda h: h.entrez_id.astype(int)).drop_duplicates('entrez_id').set_index('entrez_id')
cn = pd.read_csv('raw/cvcl2ach.tsv', sep='\t')
GA = pd.read_parquet('work/claims/gene_ann.parquet', columns=['pmid', 'gene_id', 'mention'])
GA['gene_id'] = pd.to_numeric(GA.gene_id, errors='coerce')
def gene_forms(g, sym, pm):
    f = {sym} | {m for ms in GA.loc[(GA.pmid == pm) & (GA.gene_id == g), 'mention'] for m in str(ms).split('|') if len(m) >= 3}
    if g in H.index:
        for c in ('alias_symbol', 'prev_symbol'):
            v = H.at[g, c]
            if isinstance(v, str): f |= {a for a in v.split('|') if len(a) >= 3}
    return f
def line_forms(ach, name):
    f = {name} | set(cn[cn.ach == ach].cell_name.astype(str))
    return {x for x in f if len(re.sub(r'[^A-Za-z0-9]', '', x)) >= 3}
norm = lambda s: re.sub(r'[^a-z0-9]', '', s.lower())
def fetch(pc):
    fp = f'{XD}/{pc}.xml'
    if not os.path.exists(fp):
        subprocess.run(['ncurl', '-sS', '-L', '-o', fp, f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pc}/fullTextXML'], timeout=120)
    try:
        root = ET.parse(fp).getroot(); return root if root.find('.//body') is not None else None
    except Exception: return None
def texts(root):
    body = root.find('.//body'); secs = []
    for s in body.iter('sec'):
        t = (s.findtext('title') or '').lower(); st = (s.get('sec-type') or '').lower()
        secs.append(('methods' if 'method' in st or 'method' in t or 'material' in t else 'other', ' '.join(s.itertext())))
    legends = [' '.join(c.itertext()) for c in root.iter('caption')]
    full = ' '.join(body.itertext()) + ' ' + ' '.join(legends)
    meth = ' '.join(t for k, t in secs if k == 'methods')
    sents = re.split(r'(?<=[.!?])\s+(?=[A-Z(])', re.sub(r'\s+', ' ', full))
    return full, meth, sents + legends
def pat(f):
    # case-insensitive; optional hyphen or space between letter and digit runs (WNT5A ~ Wnt-5a)
    parts = re.findall(r'[A-Za-z]+|[0-9]+|[^A-Za-z0-9]+', f)
    return r'(?<![A-Za-z0-9])' + r'[-\s]?'.join(re.escape(p) for p in parts if p.strip('-  ')) + r'(?![A-Za-z0-9])'
def has(s, forms):
    # ignore case only for forms that cannot be ordinary words (a digit, or five or more characters)
    return any(re.search(pat(f), s, re.I if (re.search(r'\d', f) or len(f) >= 5) else 0) for f in forms)
def hasl(s, forms):
    ns = norm(s); return any(norm(f) in ns for f in forms)
with P.run(__file__, params=dict(n_per_stratum=50, seed=20261104)) as run:
    A = pd.read_csv('release/claims_scored.tsv.gz', sep='\t')
    C = A[(A.split == 'confirmatory') & A.pmid.isin(pmc.keys())].copy()
    C['stratum'] = np.where(C.not_expressed == True, 'not expressed', 'other')
    rows = []
    for s, d in C.groupby('stratum'):
        papers = d.pmid.unique(); rng.shuffle(papers); got = 0
        for pm in papers:
            if got >= 50: break
            root = fetch(pmc[pm])
            if root is None: continue
            r = d[d.pmid == pm].sample(1, random_state=int(rng.integers(10**9))).iloc[0]
            full, meth, sents = texts(root)
            gf, lf = gene_forms(int(r.entrez_id), r.symbol, int(pm)), line_forms(r.depmap_model, r.cell_line)
            g_re = [x for x in sents if has(x, gf)]
            c1 = any(re.search(REAG, x, re.I) for x in g_re)
            c2 = hasl(meth, lf) if meth else hasl(full, lf)
            c3 = bool(ASSAY.search(full))
            c4 = any(re.search(VER, x, re.I) for x in g_re)
            c5 = [x for x in g_re if hasl(x, lf) and re.search(REAG, x, re.I)]
            rows.append(dict(stratum=s, pmid=int(pm), pmcid=pmc[pm], gene_symbol=r.symbol, entrez_id=int(r.entrez_id), cell_line=r.cell_line,
                             depmap_model=r.depmap_model, crispr_supported=bool(r.crispr_supported), methods_section=bool(meth),
                             reagent_with_gene=c1, line_in_methods=c2, growth_assay=c3, verification_with_gene=c4, gene_line_reagent_sentence=bool(c5),
                             passages=' || '.join((c5 or [x for x in g_re if re.search(REAG, x, re.I)])[:3])[:3000]))
            got += 1
    F = pd.DataFrame(rows)
    chk = ['reagent_with_gene', 'line_in_methods', 'growth_assay', 'verification_with_gene', 'gene_line_reagent_sentence']
    S = F.groupby('stratum')[chk].mean().T; S['all'] = F[chk].mean(); print(F.groupby('stratum').size()); print(S.round(3))
    F['all_checks'] = F[chk].all(axis=1); print(F.groupby('stratum').all_checks.mean())
    F.to_csv(f'validation/fulltext_items{suf}.tsv', sep='\t', index=False)
    O = S.reset_index().melt(id_vars='index', var_name='stratum', value_name='share').rename(columns={'index': 'check'})
    run.table('fulltext_checks_additional', O, note='open-access full-text checks of abstract claims (additional analysis, plan 2026-11-04)')

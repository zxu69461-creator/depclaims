"""Rule-based tagging of gene-perturbation growth claims (M5 job). Usage: 02_tag_claims.py CAND_DIR GENE_GZ CELL_GZ CVCL2ACH_TSV OUTDIR NPROC
For each candidate abstract: PubTator3 gene mentions (any species id kept; human filter later) and Cellosaurus cell lines; a gene is the
perturbed gene of a sentence if a perturbation pattern wraps one of its mention strings; the claim direction comes from growth-phenotype verbs in
the same sentence (or the next one when the perturbation sentence has none). Output one row per (pmid, gene, sentence) claim with method and
direction, plus the abstract's DepMap-mapped cell lines."""
import sys, os, re, gzip, glob, duckdb, multiprocessing as mp
import pandas as pd
cand_dir, gene_gz, cell_gz, map_tsv, outdir, nproc = sys.argv[1:7]; nproc = int(nproc); os.makedirs(outdir, exist_ok=True)
PV = r"(?:knock(?:ing|ed)?[- ]?down|knock(?:ing|ed)?[- ]?out|silenc(?:ing|ed|e)|deplet(?:ion|ing|ed)|ablation|deletion|loss|down-?regulation|inhibition|suppression|disruption|depleting|knocking down|knocking out)"
PN = r"(?:knock-?down|knock-?out|KD|KO|silencing|silenced|depletion|depleted|deficien(?:cy|t)|-?null|siRNAs?|shRNAs?|sgRNAs?|ASOs?)"
GROW = r"(?:prolifer\w*|growth|viabilit\w*|colon(?:y|ies)|clonogen\w*|tumou?rigen\w*|xenograft\w*|survival|cell cycle|apoptos\w*|DNA synthesis)"
NEG = re.compile(r"(?i)\b(inhibit\w*|suppress\w*|reduc\w*|decreas\w*|impair\w*|attenuat\w*|abrogat\w*|block\w*|abolish\w*|diminish\w*|repress\w*|"
                 r"retard\w*|restrain\w*|arrest\w*|lower\w*|slow\w*|compromis\w*|halt\w*|prevent\w*|limit\w*)\b[^.;]{0,80}?\b" + GROW + r"|"
                 + GROW + r"[^.;]{0,40}?\b(was|were|is|are)?\s*(markedly |significantly |dramatically |strongly )?(inhibit\w*|suppress\w*|reduc\w*|decreas\w*|impair\w*|attenuat\w*|abrogat\w*|block\w*|arrest\w*)|"
                 r"\b(induc\w*|trigger\w*|promot\w*|increas\w*|enhanc\w*)\b[^.;]{0,40}?\bapoptos\w*|\bG[012]/?[MS]?(?:-phase)? (?:phase )?arrest")
POS = re.compile(r"(?i)\b(promot\w*|enhanc\w*|increas\w*|accelerat\w*|stimulat\w*|augment\w*|facilitat\w*|boost\w*)\b[^.;]{0,60}?\b(prolifer\w*|growth|viabilit\w*|colon(?:y|ies)|clonogen\w*|tumou?rigen\w*)")
NOEFF = re.compile(r"(?i)\b(no|not|did not|does not|failed to|without)\b[^.;]{0,40}?\b(affect|alter|change|influence|impact|inhibit|reduce|effect)\w*[^.;]{0,40}?\b(prolifer\w*|growth|viabilit\w*|colon(?:y|ies)|survival)")
METH = [('CRISPR', re.compile(r'(?i)crispr|cas9|sgrna|knock-?out|\bKO\b|gene editing|deletion')), ('RNAi', re.compile(r'(?i)sirna|shrna|\bsi-|\bsh-|rna interference|\brnai\b|knock-?down|silenc')),
        ('ASO', re.compile(r'(?i)antisense|\bASOs?\b|morpholino'))]
SENT = re.compile(r'(?<=[.;!?])\s+(?=[A-Z(\[])')
def pats(m):
    e = re.escape(m)
    return re.compile(r'(?:\b' + PV + r'\s+(?:of\s+)?(?:the\s+)?(?:endogenous\s+)?' + e + r'(?![\w-])|(?<![\w-])' + e + r'(?:[- ]\s*(?:gene\s+|expression\s+|mRNA\s+)?)' + PN +
                      r'\b|(?<![\w])(?:si|sh|sg)-?' + e + r'(?![\w])|(?:siRNAs?|shRNAs?|sgRNAs?|antisense\s+oligonucleotides?)\s+(?:against|targeting|specific\s+(?:for|to))\s+' + e + r'(?![\w])|'
                      + e + r'-(?:specific|targeting|targeted)\s+(?:siRNAs?|shRNAs?|sgRNAs?))', re.I if len(m) > 4 else 0)
def work(seq_path):
    seq, path = seq_path
    out = os.path.join(outdir, f'claims_{seq:04d}.parquet')
    if os.path.exists(out): return seq, 'skip', 0
    con = duckdb.connect()
    ab = con.execute(f"SELECT pmid, title, abstract FROM read_parquet('{path}')").df()
    if not len(ab): pd.DataFrame().to_parquet(out); return seq, 'ok', 0
    con.register('ab', ab[['pmid']])
    gm = con.execute("SELECT pmid, gene_id, mention FROM read_parquet('" + os.path.join(outdir, 'gene_ann.parquet') + "') JOIN ab USING (pmid)").df()
    rows = []
    gby = {p: g for p, g in gm.groupby('pmid')}
    for r in ab.itertuples():
        g = gby.get(r.pmid)
        if g is None: continue
        text = (r.title or '') + ' ' + (r.abstract or ''); sents = SENT.split(text)
        for gid, ms in zip(g.gene_id, g.mention):
            names = sorted({m.strip() for m in str(ms).split('|') if len(m.strip()) >= 2}, key=len, reverse=True)
            if not names: continue
            P = [pats(m) for m in names]
            for i, s in enumerate(sents):
                hit = next((p.search(s) for p in P if p.search(s)), None)
                if not hit: continue
                s2 = s if (NEG.search(s) or POS.search(s) or NOEFF.search(s)) else (s + ' ' + sents[i + 1] if i + 1 < len(sents) else s)
                d = 'none' if NOEFF.search(s2) else 'reduces' if NEG.search(s2) else 'increases' if POS.search(s2) else 'unclear'
                meth = next((k for k, rx in METH if rx.search(hit.group(0))), next((k for k, rx in METH if rx.search(s2)), 'other'))
                rows.append((r.pmid, gid, i, d, meth, hit.group(0)[:80], s2[:600]))
    df = pd.DataFrame(rows, columns=['pmid', 'gene_id', 'sent_idx', 'direction', 'method', 'trigger', 'sentence'])
    df.to_parquet(out, index=False); return seq, 'ok', len(df)
if __name__ == '__main__':
    con = duckdb.connect()
    con.execute(f"CREATE TABLE p AS SELECT DISTINCT pmid FROM read_parquet('{cand_dir}/cand_*.parquet')")
    con.execute(f"""COPY (SELECT g.column0::BIGINT pmid, g.column2 gene_id, g.column3 mention FROM read_csv('{gene_gz}', delim='\t', header=false, quote='', all_varchar=true) g
                  JOIN p ON p.pmid = g.column0::BIGINT) TO '{outdir}/gene_ann.parquet' (FORMAT parquet)""")
    con.execute(f"""COPY (SELECT c.column0::BIGINT pmid, c.column2 cvcl, c.column3 mention, m.ach FROM read_csv('{cell_gz}', delim='\t', header=false, quote='', all_varchar=true) c
                  JOIN p ON p.pmid = c.column0::BIGINT LEFT JOIN read_csv('{map_tsv}', delim='\t', header=true) m ON m.cvcl = c.column2) TO '{outdir}/cell_ann.parquet' (FORMAT parquet)""")
    print('annotations filtered', con.execute(f"select count(*) from '{outdir}/gene_ann.parquet'").fetchone(), con.execute(f"select count(*) from '{outdir}/cell_ann.parquet'").fetchone(), flush=True)
    files = sorted(glob.glob(os.path.join(cand_dir, 'cand_*.parquet')))
    jobs = [(int(re.search(r'cand_(\d+)', f).group(1)), f) for f in files]; n = 0
    with mp.Pool(nproc) as pool:
        for seq, st, k in pool.imap_unordered(work, jobs): n += k
    print('claims rows', n)

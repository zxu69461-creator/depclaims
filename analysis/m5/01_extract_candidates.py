"""Candidate abstracts for gene-perturbation growth claims in cell lines, from the PubMed 2026 baseline (M5 job).
Usage: 01_extract_candidates.py BASELINE_DIR OUTDIR NPROC. Recall-first: a record is kept if title+abstract name a perturbation and a
growth/viability phenotype and mention cells; claim reading happens later. Writes OUTDIR/cand_<seq>.parquet (+ all-record metadata
OUTDIR/meta_<seq>.parquet: pmid, year, journal, issn, country, ptypes, retracted, for denominators)."""
import sys, os, gzip, re, glob, time, multiprocessing as mp
import xml.etree.ElementTree as ET
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pmparse as L
PERT = re.compile(r"(?i)\b(knock(?:ed|ing)?[- ]?down|knock(?:ed|ing)?[- ]?out|si-?rnas?|sh-?rnas?|sgrnas?|silenc\w*|deplet\w*|crispr|"
                  r"loss[- ]of[- ]function|antisense|morpholino|downregulation of|down-regulation of|inhibition of \w+ expression)")
PHEN = re.compile(r"(?i)\b(prolifer\w*|viabilit\w*|cell growth|growth of|colony|colonies|clonogenic|tumou?rigen\w*|xenograft|apoptos\w*|"
                  r"cell cycle arrest|growth inhibition|cell survival|survival of)")
CELL = re.compile(r"(?i)\bcells?\b|cell lines?")
CTY = None
def work(job):
    seq, path, outdir = job
    out = os.path.join(outdir, f"cand_{seq:04d}.parquet")
    if os.path.exists(out): return seq, 'skip', 0
    import pyarrow as pa, pyarrow.parquet as pq
    cty = L.load_countries(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'countries_v1.tsv'))
    C = {k: [] for k in ('pmid', 'version', 'pub_year', 'journal', 'issn', 'country', 'ptypes', 'retracted', 'title', 'abstract')}
    M = {k: [] for k in ('pmid', 'version', 'pub_year', 'issn', 'country', 'ptypes', 'retracted')}
    with gzip.open(path) as fh:
        for ev, el in ET.iterparse(fh, events=('end',)):
            if el.tag != 'PubmedArticle': continue
            mc = el.find('MedlineCitation'); a = mc.find('Article'); pm = mc.find('PMID')
            title = L.txt(a.find('ArticleTitle')); ab = ' '.join(L.txt(t) for t in a.findall('Abstract/AbstractText'))
            ptypes = '|'.join(t.text for t in a.findall('PublicationTypeList/PublicationType') if t.text)
            auths = a.findall('AuthorList/Author'); aff = L.txt(auths[0].find('AffiliationInfo/Affiliation')) if auths else ''
            ret = int(any(c.get('RefType') == 'RetractionIn' for c in mc.findall('CommentsCorrectionsList/CommentsCorrections')))
            j = a.find('Journal'); issn = j.findtext('ISSN') or ''; yr = L.pub_year(a); co = L.country(aff, cty)
            for k, v in zip(M, (int(pm.text), int(pm.get('Version') or 1), yr, issn, co, ptypes, ret)): M[k].append(v)
            text = title + ' ' + ab
            if PERT.search(text) and PHEN.search(text) and CELL.search(text):
                for k, v in zip(C, (int(pm.text), int(pm.get('Version') or 1), yr, j.findtext('Title') or '', issn, co, ptypes, ret, title, ab)):
                    C[k].append(v)
            el.clear()
    pq.write_table(pa.table(M), os.path.join(outdir, f"meta_{seq:04d}.parquet"))
    pq.write_table(pa.table(C), out + '.tmp'); os.replace(out + '.tmp', out)
    return seq, 'ok', len(C['pmid'])
if __name__ == '__main__':
    base, outdir, nproc = sys.argv[1], sys.argv[2], int(sys.argv[3]); os.makedirs(outdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(base, 'pubmed26n*.xml.gz')))
    jobs = [(int(re.search(r'n(\d{4})', f).group(1)), f, outdir) for f in files]
    t0 = time.time(); n = 0
    with mp.Pool(nproc) as p:
        for seq, st, k in p.imap_unordered(work, jobs):
            n += k; print(seq, st, k, f'{time.time() - t0:.0f}s', flush=True)
    print('done candidates', n)

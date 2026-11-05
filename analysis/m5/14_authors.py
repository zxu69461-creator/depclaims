"""Author lists for claim papers, from the PubMed 2026 baseline (M5 job). Usage: 14_authors.py BASELINE_DIR PMIDS_TXT OUT_PARQUET NPROC (defaults: mjob IN_DIR/OUT_DIR).
One row per (pmid, author position): last name, fore name, initials, ORCID when deposited. Last record version wins."""
import sys, os, gzip, glob, multiprocessing as mp
import xml.etree.ElementTree as ET
WANT = None
def init(f):
    global WANT; WANT = {int(x) for x in open(f) if x.strip()}
def work(path):
    rows = []; seq = int(os.path.basename(path)[9:13])
    with gzip.open(path) as fh:
        for ev, el in ET.iterparse(fh, events=('end',)):
            if el.tag != 'PubmedArticle': continue
            mc = el.find('MedlineCitation'); pm = int(mc.findtext('PMID'))
            if pm in WANT:
                for i, a in enumerate(mc.findall('Article/AuthorList/Author')):
                    orcid = next((x.text for x in a.findall('Identifier') if x.get('Source') == 'ORCID' and x.text), None)
                    rows.append((pm, seq, i, a.findtext('LastName') or a.findtext('CollectiveName') or '', a.findtext('ForeName') or '',
                                 a.findtext('Initials') or '', orcid))
            el.clear()
    return rows
if __name__ == '__main__':
    a = sys.argv[1:] or [os.path.expanduser('~/echo/raw/baseline'), os.path.join(os.environ['IN_DIR'], 'claim_pmids.txt'),
                         os.path.join(os.environ['OUT_DIR'], 'authors.parquet'), '14']
    base, pm, out, n = a[0], a[1], a[2], int(a[3])
    files = sorted(glob.glob(os.path.join(base, 'pubmed26n*.xml.gz')))
    with mp.Pool(n, initializer=init, initargs=(pm,)) as p: R = [r for rs in p.imap_unordered(work, files) for r in rs]
    import pandas as pd
    D = pd.DataFrame(R, columns=['pmid', 'file_seq', 'pos', 'last', 'fore', 'initials', 'orcid'])
    D = D[D.file_seq == D.groupby('pmid').file_seq.transform('max')].drop(columns='file_seq')
    D.to_parquet(out, index=False); print('papers', D.pmid.nunique(), 'author rows', len(D), 'with ORCID', D.orcid.notna().mean().round(3))

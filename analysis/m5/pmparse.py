"""PubMed parse helpers copied verbatim from echo/scripts/lib_echo_parse.py (country, txt, date_of, pub_year)."""
import re, csv
def load_countries(tsv):
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    exact_cs = {"US", "UK", "USA", "PRC", "DRC", "UAE", "KSA", "DPRK", "HKSAR", "U.S.", "U.K.", "U.S.A.", "U.S.A", "U.A.E."}
    names = {}
    for r in rows:
        if r["kind"] == "us_state_abbr_after_comma":
            continue
        names[r["alias"] if r["alias"] in exact_cs else r["alias"].lower()] = r["country"]
    abbr = {r["alias"] for r in rows if r["kind"] == "us_state_abbr_after_comma"}
    low = sorted((w for w in names if w.islower() or not w.isupper()), key=len, reverse=True)
    low = [w for w in low if w not in exact_cs and len(w) > 3]
    low_re = re.compile("(?<![\\w])(" + "|".join(re.escape(w) for w in low) + ")(?![\\w])")
    cs = sorted((w for w in exact_cs if w in names), key=len, reverse=True)
    cs_re = re.compile("(?<![\\w.])(" + "|".join(re.escape(w) for w in cs) + ")(?![\\w])")
    return names, abbr, low_re, cs_re
EMAIL = re.compile(r"\S+@\S+|electronic address.*$", re.I)
ZIP = re.compile(r"^([A-Z]{2})\s+\d{5}(?:-\d{4})?$")
def country(aff, cdict):
    names, abbr, low_re, cs_re = cdict
    if not aff: return "unknown"
    a = EMAIL.sub("", aff)
    if re.search(r"democratic people'?s republic of korea", a, re.I): return "North Korea"
    segs = [x.strip(" .") for x in re.split(r"[,;]", a) if x.strip(" .")]
    for seg in reversed(segs):
        if seg in names: return names[seg]
        lowseg = seg.lower()
        if lowseg in names: return names[lowseg]
        m = ZIP.match(seg)
        if (m and m.group(1) in abbr) or seg in abbr: return "United States"
        m = cs_re.search(seg)
        if m: return names[m.group(1)]
        m = low_re.search(lowseg)
        if m: return names[m.group(1)]
    return "unknown"

def txt(e):
    return re.sub(r"\s+", " ", "".join(e.itertext())).strip() if e is not None else ""
MON = {"jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6, "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12}
def date_of(e):
    if e is None: return None
    y, m, d = e.findtext("Year"), e.findtext("Month") or "1", e.findtext("Day") or "1"
    if not y: return None
    m = MON.get(m[:3].lower(), m) if not str(m).isdigit() else m
    try: return f"{int(y):04d}-{int(m):02d}-{int(d):02d}"
    except Exception: return None
def pub_year(a):
    pd = a.find("Journal/JournalIssue/PubDate")
    if pd is not None:
        y = pd.findtext("Year")
        if y and y[:4].isdigit(): return int(y[:4])
        md = pd.findtext("MedlineDate") or ""
        m = re.search(r"(\d{4})", md)
        if m: return int(m.group(1))
    ad = a.find("ArticleDate")
    if ad is not None and (ad.findtext("Year") or "").isdigit(): return int(ad.findtext("Year"))
    return None

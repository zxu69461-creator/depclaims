"""Score a published gene-line growth claim against DepMap screens of the same cell line.

The rulers are the ones used for the atlas: CRISPR (Chronos) and RNAi (DEMETER2) gene effects are scaled per cell line so that the
median of common-essential genes (CEGv2) is -1 and the median of non-essential genes (NEGv1) is 0; a gene is needed in a line when its
scaled effect is <= -0.5. Gene classes follow the share of screened lines that need the gene: common (>= 90%), selective (> 0) or never.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd

THRESHOLD = -0.5
COMMON = 0.9


def _entrez_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Keep 'SYMBOL (ENTREZ)' columns and index them by Entrez id; remember symbols."""
    keep = [c for c in df.columns if '(' in c and '&' not in c]
    df = df[keep]
    df.columns = [int(c.split('(')[-1].rstrip(')')) for c in keep]
    return df.loc[:, ~df.columns.duplicated()]


def scale(X: pd.DataFrame, ceg: list[int], neg: list[int]) -> pd.DataFrame:
    """Per-line scaling: non-essential median -> 0, common-essential median -> -1."""
    mn = X[[g for g in neg if g in X.columns]].median(axis=1)
    mc = X[[g for g in ceg if g in X.columns]].median(axis=1)
    return X.sub(mn, axis=0).div(mn - mc, axis=0)


def _key(s: str) -> str:
    return re.sub(r'[^A-Z0-9]', '', str(s).upper())


def build(chronos: str, demeter2: str, expression: str, models: str, ceg: str, neg: str, out: str) -> Path:
    """Build a compact cache from the public DepMap downloads (paths to the CSV/TSV files as released)."""
    out_dir = Path(out); out_dir.mkdir(parents=True, exist_ok=True)
    head = open(chronos).readline().strip().split(',')[1:]
    sym = {int(c.split('(')[-1].rstrip(')')): c.split(' (')[0] for c in head if '(' in c}
    ceg_ids = pd.read_csv(ceg, sep='\t').ENTREZ_ID.dropna().astype(int).tolist()
    neg_ids = pd.read_csv(neg, sep='\t').ENTREZ_ID.dropna().astype(int).tolist()
    cr = scale(_entrez_columns(pd.read_csv(chronos, index_col=0)), ceg_ids, neg_ids).astype('float32')
    d2 = _entrez_columns(pd.read_csv(demeter2, index_col=0))
    d2 = d2[d2.index.str.startswith('ACH-')]
    d2 = scale(d2[~d2.index.duplicated()], ceg_ids, neg_ids).astype('float32')
    X = pd.read_csv(expression, low_memory=False)
    if 'IsDefaultEntryForModel' in X.columns:
        X = X[X.IsDefaultEntryForModel.astype(str).str.lower().isin(['yes', 'true'])]
    X = X.set_index('ModelID' if 'ModelID' in X.columns else X.columns[0])
    cols = [c for c in X.columns if '(' in c]
    ids = [int(c.split('(')[-1].rstrip(')')) for c in cols]
    dup = pd.Series(ids).duplicated(keep=False).values
    cols = [c for c, i, d in zip(cols, ids, dup) if not d or c.split(' (')[0] == sym.get(i, c.split(' (')[0])]
    X = _entrez_columns(X[cols]); X = X[~X.index.duplicated()].astype('float32')
    M = pd.read_csv(models)
    M = M[['ModelID', 'StrippedCellLineName', 'CellLineName', 'OncotreeLineage', 'OncotreePrimaryDisease']].copy()
    for name, df in (('crispr', cr), ('rnai', d2), ('expression', X)):
        df.columns = df.columns.astype(str); df.to_parquet(out_dir / f'{name}.parquet')
    M.to_parquet(out_dir / 'models.parquet')
    pd.Series(sym, name='symbol').rename_axis('entrez_id').reset_index().to_parquet(out_dir / 'genes.parquet')
    return out_dir


@dataclass
class Score:
    gene: str
    entrez_id: int
    cell_line: str
    depmap_model: str
    lineage: str
    crispr_effect_scaled: float | None
    needed_in_line: bool | None
    gene_class: str
    dependent_fraction_all: float
    dependent_fraction_lineage: float | None
    crispr_rank_in_lineage: float | None
    rnai_effect_scaled: float | None
    rnai_rank_in_lineage: float | None
    rnai_top_line_in_lineage: str | None
    expression_tpm: float | None
    lines_in_lineage: int

    def as_dict(self) -> dict:
        return asdict(self)


class Atlas:
    """Look-ups against a cache made by build(); optionally joined to the published-claims table (claims_scored.tsv.gz)."""

    def __init__(self, cache: str, claims: str | None = None):
        c = Path(cache)
        self.cr = pd.read_parquet(c / 'crispr.parquet'); self.rn = pd.read_parquet(c / 'rnai.parquet')
        self.ex = pd.read_parquet(c / 'expression.parquet'); self.models = pd.read_parquet(c / 'models.parquet')
        g = pd.read_parquet(c / 'genes.parquet')
        self.sym2id = {s.upper(): int(i) for i, s in zip(g.entrez_id, g.symbol)}
        self.id2sym = dict(zip(g.entrez_id.astype(int), g.symbol))
        self.line2id = {}
        for r in self.models.itertuples():
            for k in (r.ModelID, r.StrippedCellLineName, r.CellLineName):
                if isinstance(k, str): self.line2id.setdefault(_key(k), r.ModelID)
        self.lineage = self.models.set_index('ModelID').OncotreeLineage
        self.names = self.models.set_index('ModelID').StrippedCellLineName
        self.claims = pd.read_csv(claims, sep='\t') if claims else None

    def gene_id(self, gene) -> int:
        if str(gene).isdigit(): return int(gene)
        try: return self.sym2id[str(gene).upper()]
        except KeyError: raise KeyError(f'gene {gene!r} is not in the screened gene set') from None

    def model_id(self, line) -> str:
        try: return self.line2id[_key(line)]
        except KeyError: raise KeyError(f'cell line {line!r} is not a DepMap model') from None

    @staticmethod
    def _rank(v: pd.Series, model: str) -> float | None:
        v = v.dropna()
        if model not in v.index or len(v) < 5: return None
        y, oth = v[model], v.drop(model)
        return float((oth < y).mean() + 0.5 * (oth == y).mean())

    def score(self, gene, line) -> Score:
        g, m = self.gene_id(gene), self.model_id(line); col = str(g)
        lin = self.lineage.get(m)
        same = self.lineage[self.lineage == lin].index
        c = self.cr[col] if col in self.cr.columns else pd.Series(dtype=float)
        r = self.rn[col] if col in self.rn.columns else pd.Series(dtype=float)
        frac = float((c <= THRESHOLD).sum() / len(c)) if len(c) else float('nan')   # over all lines in the release, as in the atlas
        cls = 'common' if frac >= COMMON else 'selective' if frac > 0 else 'never'
        cl = c.reindex(same).dropna(); rl = r.reindex(same).dropna()
        ce = float(c[m]) if m in c.index and pd.notna(c[m]) else None
        re_ = float(r[m]) if m in r.index and pd.notna(r[m]) else None
        ex = self.ex[col].get(m) if col in self.ex.columns else None
        return Score(gene=self.id2sym.get(g, str(g)), entrez_id=g, cell_line=self.names.get(m, m), depmap_model=m, lineage=lin,
                     crispr_effect_scaled=ce, needed_in_line=None if ce is None else ce <= THRESHOLD, gene_class=cls,
                     dependent_fraction_all=frac, dependent_fraction_lineage=float((cl <= THRESHOLD).mean()) if len(cl) else None,
                     crispr_rank_in_lineage=self._rank(cl, m), rnai_effect_scaled=re_, rnai_rank_in_lineage=self._rank(rl, m),
                     rnai_top_line_in_lineage=self.names.get(rl.idxmin(), rl.idxmin()) if len(rl) >= 5 else None,
                     expression_tpm=None if ex is None or pd.isna(ex) else float(2 ** ex - 1), lines_in_lineage=len(cl))

    def lineage_profile(self, gene, lineage: str) -> pd.DataFrame:
        """Scaled CRISPR and RNAi effects of the gene in every screened line of a lineage, most dependent first."""
        col = str(self.gene_id(gene)); same = self.lineage[self.lineage == lineage].index
        df = pd.DataFrame({'crispr_effect_scaled': self.cr[col].reindex(same) if col in self.cr.columns else np.nan,
                           'rnai_effect_scaled': self.rn[col].reindex(same) if col in self.rn.columns else np.nan})
        df.insert(0, 'cell_line', [self.names.get(m, m) for m in df.index])
        return df.dropna(subset=['crispr_effect_scaled']).sort_values('crispr_effect_scaled')

    def published(self, gene, line=None) -> pd.DataFrame:
        """Published claims in the atlas for a gene (and line)."""
        if self.claims is None: raise ValueError('pass claims=claims_scored.tsv.gz to Atlas() to look up published claims')
        d = self.claims[self.claims.entrez_id == self.gene_id(gene)]
        return d if line is None else d[d.depmap_model == self.model_id(line)]

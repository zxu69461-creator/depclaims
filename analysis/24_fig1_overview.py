"""Figure 1. A screen-calibrated atlas of published cancer-cell growth claims. a, design: claims that perturbing gene X reduces growth of
cell line Y are read from PubMed abstracts, mapped to DepMap genes and models, and scored against genome-wide CRISPR and RNAi screens of the
same line; b, the corpus from PubMed records to scored claims; c, claims per year by perturbation method; d, a typical case, KRAS in
pancreatic cancer lines."""
import os
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import figkit as fk
ROOT = Path(__file__).resolve().parents[1]
rp = lambda f: duckdb.sql(f"select * from '{ROOT / f}'").df()
A = rp('work/fig17/out/claims_base.parquet'); L = rp('work/fig/line_strips.parquet')
VER, BLU, GD, GR, GL = fk.c.vermilion, fk.c.blue, fk.c.grey_dark, fk.c.grey, fk.c.grey_light
FB, FT = fk.pt('base'), fk.pt('tick')
names = pd.read_csv(ROOT / 'raw/Model.csv', usecols=['ModelID', 'StrippedCellLineName']).set_index('ModelID').StrippedCellLineName

# a design strip (real-data tiles)
W = fk.art.Art(166, 44)
xs = [10, 48, 84, 119, 153]
st = ['Published claim', 'Gene and line', 'Same line in DepMap', 'Genome-wide screens', 'Calibrated atlas']
W.track([(xs[0], 6), (xs[-1], 6)], color='grey_light', width=1.1)
for k, (x, nm) in enumerate(zip(xs, st)):
    W.station(x, 6, color='vermilion' if k == 4 else 'blue', r=1.25)
    W.text(x, 1.8, nm, role='legend', anchor='middle', bold=True)
W.stack(0.6, 13, 17, 21, n=3, dx=1.6, dy=-1.6)
W.doc(0.6, 13, 17, 21, color='grey', lines=6, mark=3)
W.text(10, 39, 'PubMed abstract', role='annot', anchor='middle')
for k, (lab, c) in enumerate((('gene X', 'blue'), ('line Y', 'vermilion'), ('growth ↓', 'grey_dark'))):
    W.rect(38, 14 + k * 7, 20, 5.2, r=2.6, fill=f'{c}@t0.8' if c != 'grey_dark' else 'grey_light')
    W.text(48, 16.6 + k * 7, lab, role='annot', anchor='middle', color=c)
W.text(48, 39, 'HGNC, Cellosaurus', role='annot', anchor='middle')
for k, (cx, cy, r) in enumerate(((77, 20, 4.2), (87, 17, 3.6), (93, 24, 3.2), (82, 28, 3.6), (91, 31, 2.8))):
    W.cell(cx, cy, r, color='vermilion' if k == 0 else 'grey')
W.text(84, 39, 'DepMap model', role='annot', anchor='middle')
br = L[L.OncotreeLineage == 'Breast'].pivot_table(index='symbol', columns='ach', values='c').dropna(axis=1)
br = br.loc[br.var(axis=1).sort_values(ascending=False).index[:10]].iloc[:, :16]
W.heat(108, 12, 22, 19, np.digitize(br.values, [-1.0, -0.5]).astype(float), cmap='diverging', vmin=-1.2, vmax=3.2)
W.text(119, 39, 'CRISPR and RNAi effects', role='annot', anchor='middle')
sup = A.sup.mean()
W.waffle(143.5, 12, 19, 19, [int(round(sup * 100)), 100 - int(round(sup * 100))], ['vermilion', 'grey_light'], ncol=10)
W.text(153, 39, 'gene, line rank, expression', role='annot', anchor='middle')
for x0, x1 in zip(xs[:-1], xs[1:]):
    W.block_arrow(x0 + 12, 22, x1 - 12, 22, width=1.3, color='grey_light')

PROFILE = __import__('os').environ.get('FIG_PROFILE', 'natcomms')     # journal profile; FIG_PROFILE=pnas writes figures/pnas/
fig = fk.figure('double', height=135, profile=PROFILE)
gs = fig.add_gridspec(2, 3, height_ratios=[0.62, 1], width_ratios=[0.95, 1.1, 1.2])
ax = fig.add_subplot(gs[0, :]); fk.art.place(ax, W); ax.set_anchor('NW')
fk.panel(ax, 'a', 'Each published growth claim is scored against genome-wide screens of the same cell line')

# b claims per year by method
ax = fig.add_subplot(gs[1, 1]); fk.panel(ax, 'c', 'Claims accumulate from 2005')
yr = A.groupby(['yr', 'meth']).size().unstack(fill_value=0).reindex(range(1995, 2026), fill_value=0)[['RNAi', 'other', 'CRISPR']]
bottom = np.zeros(len(yr))
cols = {'RNAi': GD, 'other': GL, 'CRISPR': BLU}
for m in yr.columns:
    ax.bar(yr.index, yr[m], bottom=bottom, width=0.82, color=cols[m], lw=0); bottom += yr[m].values
ax.text(1996, 1110, 'siRNA, shRNA', color=GD, fontsize=FT, ha='left'); ax.text(1996, 1035, 'knockdown, unspecified', color=GR, fontsize=FT, ha='left')
ax.text(1996, 960, 'CRISPR', color=BLU, fontsize=FT, ha='left')
ax.set_xlim(1994.3, 2025.7); ax.set_xticks([2000, 2010, 2020]); ax.set_ylim(0, 1250); ax.set_yticks([0, 500, 1000])
fk.chart(ax, grid='y', unit='Claims per year')

# c corpus
ax = fig.add_subplot(gs[1, 0]); fk.panel(ax, 'b', 'From PubMed to scored claims')
NST, NRED = duckdb.sql(f"""select count(*), sum(case when grounded and (growth_effect = 'decreased' or apoptosis_increased) then 1 else 0 end)
                       from read_csv('{ROOT}/release/statements_all.tsv.gz')""").fetchone()
NSC = duckdb.sql(f"select count(*) from read_csv('{ROOT}/release/claims_scored.tsv.gz')").fetchone()[0]
steps = [('PubMed records', 40.0e6), ('abstracts reporting a perturbation and growth', 264080), ('with a candidate gene and DepMap line', 32939),
         ('gene–line growth statements', NST), ('claims of reduced growth', NRED), ('scored in DepMap', NSC)]
sub = {'scored in DepMap': f'scored in DepMap (n = {A.gid.nunique():,} genes, {A.ach.nunique()} lines)'}
y = np.arange(len(steps))[::-1]
for yi, (lab, v) in zip(y, steps):
    ax.barh(yi, np.log10(v), height=0.5, color=VER if lab == 'scored in DepMap' else GL, lw=0)
    ax.text(0.1, yi + 0.42, sub.get(lab, lab), ha='left', va='bottom', fontsize=FT, color=GD)
    ax.text(np.log10(v) + 0.1, yi, f'{v / 1e6:.1f}M' if v >= 1e6 else f'{v:,.0f}', ha='left', va='center', fontsize=FT, color=VER if lab == 'scored in DepMap' else GD)
ax.set_xlim(0, 9.4); ax.set_ylim(-0.5, len(steps) - 0.05); ax.set_yticks([]); ax.set_xticks([0, 2, 4, 6, 8]); ax.set_xticklabels(['1', '10$^{2}$', '10$^{4}$', '10$^{6}$', '10$^{8}$'])
fk.chart(ax, grid=None); ax.spines['left'].set_visible(False); ax.spines['bottom'].set_bounds(0, 8)

# d typical case
ax = fig.add_subplot(gs[1, 2]); fk.panel(ax, 'd', 'A typical case: KRAS in pancreatic lines')
MO = duckdb.sql(f"select ModelID as ach, OncotreeLineage from '{ROOT}/work/tool_cache/models.parquet'").df()
CR = duckdb.sql(f"""select "__index_level_0__" as ach, "3845" as c from '{ROOT}/work/tool_cache/crispr.parquet'""").df()
o = CR.merge(MO, on='ach').query("OncotreeLineage == 'Pancreas'").dropna(subset=['c']).sort_values('c').reset_index(drop=True)   # all CRISPR-screened lines
named = A[(A.symbol == 'KRAS') & (A.OncotreeLineage == 'Pancreas')].groupby('ach').pmid.nunique()
assert named.index.isin(o.ach).all()
yj = np.random.default_rng(5).uniform(-0.18, 0.18, len(o))

ax.scatter(o.c, yj, s=14, color=GR, lw=0, zorder=2)
nm = o[o.ach.isin(named.index)]
ax.scatter(nm.c, yj[nm.index], s=30 + 8 * named.reindex(nm.ach).values, color=VER, lw=0, zorder=3)
ax.axvline(-0.5, color=GD, lw=0.6, ls=(0, (3, 2)))
ax.text(-0.6, -0.36, '← dependent', ha='right', va='center', fontsize=FT, color=VER)
ax.text(-3.4, 0.5, f'{len(o)} screened lines; named in claims:\n' + ', '.join(names.get(a, a) for a in nm.sort_values('c').ach[:3]) + ',\n' + ', '.join(names.get(a, a) for a in nm.sort_values('c').ach[3:6]),
        ha='left', va='center', fontsize=FT, color=GD)
ax.set_xlim(-3.5, 0.6); ax.set_ylim(-0.45, 0.72); ax.set_yticks([]); ax.set_xticks([-3, -2, -1, -0.5, 0])
fk.chart(ax, grid=None); ax.spines['left'].set_visible(False)
ax.set_xlabel('Scaled CRISPR effect', fontsize=FB, color=GD)

out = os.environ.get('FIG_DRAFT') or ROOT / ('figures/pnas' if PROFILE == 'pnas' else 'figures/main') / 'fig1_overview'
saved = fk.save(fig, out, strict=not os.environ.get('FIG_DRAFT')); print(saved.report)

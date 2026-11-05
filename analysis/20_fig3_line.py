"""Figure 3. The lines named in papers carry little information about which lines depend on the gene (confirmatory claims on genes whose
CRISPR and RNAi profiles agree, r_x >= 0.3). a, within-lineage rank of the named line vs the line an RNAi screen ranks most dependent;
b, stricter gene agreement thresholds; c, by lineage; d, two specimen genes."""
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import figkit as fk
ROOT = Path(__file__).resolve().parents[1]
rp = lambda f: duckdb.sql(f"select * from '{ROOT / f}'").df()
S = rp('work/fig/line_claims.parquet'); L = rp('work/fig/line_strips.parquet')
R = duckdb.connect(str(ROOT / 'db/depclaims.duckdb'), read_only=True).execute(
    "select quantity, n, est, lo, hi from confirm_results where split = 'confirmatory' and quantity like '%P2%'").df().set_index('quantity')
VER, BLU, GD, GR, GL = fk.c.vermilion, fk.c.blue, fk.c.grey_dark, fk.c.grey, fk.c.grey_light
FB, FT = fk.pt('base'), fk.pt('tick')
PROFILE = __import__('os').environ.get('FIG_PROFILE', 'natcomms')     # journal profile; FIG_PROFILE=pnas writes figures/pnas/
fig = fk.figure('double', height=125, profile=PROFILE)
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.9], width_ratios=[1, 1, 1.05])

# a hero
ax = fig.add_subplot(gs[0:2, 0:2]); fk.panel(ax, 'a', 'Named lines rank near random selection, whereas screens pinpoint dependent models')
bins = np.linspace(0, 1, 21); w = bins[1] - bins[0]
hc, _ = np.histogram(S.pct_CRISPR, bins, density=True); hk, _ = np.histogram(S.ctl.dropna(), bins, density=True)
ax.fill_between(np.repeat(bins, 2)[1:-1], 0, np.repeat(hc, 2), color=VER, alpha=0.35, lw=0, zorder=2)
ax.step(bins, np.r_[hc, hc[-1]], where='post', color=VER, lw=1.8, zorder=4)
ax.step(bins, np.r_[hk, hk[-1]], where='post', color=BLU, lw=1.2, zorder=3)
ax.axhline(1, color=GR, lw=0.7, ls=(0, (3, 2)), zorder=1)
ax.text(1.015, 1, 'a random\nline', ha='left', va='center', fontsize=FB, color=GR, linespacing=1.05)
ax.vlines(S.pct_CRISPR, -0.30, -0.08, color=VER, lw=0.25, alpha=0.18, rasterized=True)
ax.vlines(S.ctl.dropna(), -0.58, -0.36, color=BLU, lw=0.25, alpha=0.18, rasterized=True)
ax.text(-0.012, -0.19, 'claims', ha='right', va='center', fontsize=FB, color=VER)
ax.text(-0.012, -0.47, 'screen', ha='right', va='center', fontsize=FB, color=BLU)
top = max(hk.max(), hc.max())
for q, col, y, lab in (('P2 claimed-line percentile (r_x >= 0.3)', VER, top * 0.62, 'line named in the paper'),
                       ('P2 positive control percentile', BLU, top * 0.74, 'line an RNAi screen ranks first')):
    r = R.loc[q]; ax.plot([r.lo, r.hi], [y, y], color=col, lw=2.2, solid_capstyle='butt', zorder=5)
    ax.scatter([r.est], [y], s=30, color=col, zorder=6, lw=0)
    ax.text(r.hi + 0.02, y, f'{lab}: mean {r.est:.2f} ({r.lo:.2f}–{r.hi:.2f})', ha='left', va='center', fontsize=FB, color=col)
d = R.loc['P2 claims minus control']
fk.numbox(ax, 0.62, top * 0.50, f'gap {d.est:.2f} ({d.lo:.2f}–{d.hi:.2f})', color=VER)
ax.set_xlim(-0.13, 1.12); ax.set_ylim(-0.66, top * 1.04); ax.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax.set_yticks(range(int(top) + 1))
fk.chart(ax, grid='y', unit='Density of claims'); ax.spines['bottom'].set_bounds(0, 1)
ax.set_xlabel('Rank of the line among CRISPR-screened lines of its lineage (0 = most dependent)', fontsize=FB, color=GD)
fk.pin_inside(ax)

# b thresholds
ax = fig.add_subplot(gs[0, 2]); fk.panel(ax, 'b', 'Stricter gene agreement sharpens the screen, not the papers')
rows = []
for t in (0.3, 0.4, 0.5):
    s = S[S.r_x >= t]; rows.append((f'r ≥ {t}  (n = {len(s):,})', s.ctl.mean(), s.pct_CRISPR.mean()))
for i, (lab, a_, b_) in enumerate(rows):
    y = 2 - i
    ax.plot([a_, b_], [y, y], color=GL, lw=2.4, solid_capstyle='butt', zorder=1)
    ax.scatter([a_], [y], s=26, color=BLU, marker='D', zorder=3, lw=0); ax.scatter([b_], [y], s=30, color=VER, zorder=3, lw=0)
    ax.text(0.0, y + 0.32, lab, ha='left', va='center', fontsize=FB, color=GD)
ax.axvline(0.5, color=GR, lw=0.6, ls=(0, (3, 2)))
ax.set_xlim(0, 0.62); ax.set_ylim(-0.5, 2.7); ax.set_yticks([]); ax.set_xticks([0, 0.25, 0.5])
fk.chart(ax, grid=None); ax.spines['left'].set_visible(False)
ax.set_xlabel('Mean within-lineage rank', fontsize=FB, color=GD)

# c lineages
ax = fig.add_subplot(gs[1, 2]); fk.panel(ax, 'c', 'The pattern holds across lineages')
g = S.groupby('OncotreeLineage').agg(n=('pmid', 'size'), m=('pct_CRISPR', 'mean'), k=('ctl', 'mean')).query('n >= 80').sort_values('m')
y = np.arange(len(g))
ax.barh(y, g.m, height=0.62, color=GL, lw=0, zorder=2)
ax.scatter(g.k, y, s=18, marker='D', color=BLU, zorder=3, lw=0)
ax.scatter(g.m, y, s=14, color=VER, zorder=3, lw=0)
for yi, (lin, r) in zip(y, g.iterrows()):
    ax.text(0.69, yi, f'{int(r.n)}', ha='right', va='center', fontsize=FT, color=GR)
ax.set_yticks(y); ax.set_yticklabels([l.replace('Esophagus/Stomach', 'Oesophagus/stomach').replace('CNS/Brain', 'CNS/brain').replace('Ovary/Fallopian Tube', 'Ovary') for l in g.index], fontsize=FT, color=GD)
ax.text(0.69, len(g) - 0.2, 'claims', ha='right', va='bottom', fontsize=FT, color=GR)
ax.axvline(0.5, color=GR, lw=0.6, ls=(0, (3, 2)))
ax.set_xlim(0, 0.70); ax.set_ylim(-0.6, len(g) + 0.6); ax.set_xticks([0, 0.25, 0.5]); ax.tick_params(axis='y', length=0)
fk.chart(ax, grid=None); ax.spines['left'].set_visible(False)
ax.set_xlabel('Mean within-lineage rank', fontsize=FB, color=GD)

# d specimens
ax = fig.add_subplot(gs[2, :]); fk.panel(ax, 'd', 'Two specimens: HER2 claims name the dependent breast lines; EGFR claims name glioma lines that do not depend on it')
names = pd.read_csv(ROOT / 'raw/Model.csv', usecols=['ModelID', 'StrippedCellLineName']).set_index('ModelID').StrippedCellLineName
for row, (gene, lin, lab) in enumerate((('ERBB2', 'Breast', 'ERBB2 (HER2), breast lines'), ('EGFR', 'CNS/Brain', 'EGFR, CNS/brain lines'))):
    y0 = 1 - row
    o = L[(L.symbol == gene) & (L.OncotreeLineage == lin)]
    named = S[(S.symbol == gene) & (S.OncotreeLineage == lin)].ach.unique()
    top = o.loc[o.r.idxmin(), 'ach']
    ax.scatter(o.c, np.full(len(o), y0), s=14, color=GL, lw=0, zorder=2)
    nm = o[o.ach.isin(named)].sort_values('c')
    ax.scatter(nm.c, np.full(len(nm), y0), s=34, color=VER, lw=0, zorder=4)
    t = o[o.ach == top]; ax.scatter(t.c, [y0], s=40, marker='D', facecolor='none', edgecolor=BLU, lw=1.1, zorder=5)
    ax.text(-3.62, y0 + 0.3, f'{lab} ({len(o)} lines)', ha='left', va='center', fontsize=FB, color=GD)
    ax.text(-3.62, y0 - 0.3, 'named: ' + ', '.join(names.get(a, a) for a in nm.ach), ha='left', va='center', fontsize=FT, color=VER)
ax.axvline(-0.5, color=GR, lw=0.6, ls=(0, (3, 2))); ax.text(-0.52, 1.42, 'dependent  ←', ha='right', va='center', fontsize=FB, color=GR)
ax.set_xlim(-3.7, 0.8); ax.set_ylim(-0.5, 1.6); ax.set_yticks([]); 
fk.chart(ax, grid=None); ax.spines['left'].set_visible(False); ax.spines['bottom'].set_bounds(-3.5, 0.8); ax.set_xticks([-3, -2, -1, -0.5, 0, 0.5])
ax.set_xlabel('CRISPR effect, scaled (0 = non-essential, −1 = common essential); open diamond: line an RNAi screen ranks first', fontsize=FB, color=GD)
fk.pin_inside(ax)
saved = fk.save(fig, ROOT / ('figures/pnas' if PROFILE == 'pnas' else 'figures/main') / 'fig3_line_specificity', strict=True); print(saved.report)

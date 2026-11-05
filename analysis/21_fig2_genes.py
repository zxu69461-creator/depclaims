"""Figure 2. Published growth claims are drawn towards genes that cancer cells need. a, claims per gene across genes ordered by the share of
CRISPR-screened lines that depend on them (bar width = share of genes; all claims, 95% CI by bootstrap over papers); b, gene classes
(confirmatory claims); c, support against the gene's dependent fraction in the named line's lineage; d, support of claims on non-common genes by
expression decile against genes at the same expression in the same line (pre-specified calibrated agreement, confirmatory claims)."""
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import figkit as fk
ROOT = Path(__file__).resolve().parents[1]
rp = lambda f: duckdb.sql(f"select * from '{ROOT / f}'").df()
D = rp('work/fig17/out/calib_decile.parquet'); CC = rp('work/fig17/out/class_comp.parquet')
A = rp('work/fig17/out/claims_base.parquet'); G = rp('work/fig17/out/genes.parquet')
R = duckdb.connect(str(ROOT / 'db/depclaims.duckdb'), read_only=True).execute(
    "select quantity, n, est, lo, hi, ref from confirm_results where split = 'confirmatory'").df().set_index('quantity')
VER, BLU, GD, GR, GL = fk.c.vermilion, fk.c.blue, fk.c.grey_dark, fk.c.grey, fk.c.grey_light
FB, FT = fk.pt('base'), fk.pt('tick')
rng = np.random.default_rng(20261103)
PROFILE = __import__('os').environ.get('FIG_PROFILE', 'natcomms')     # journal profile; FIG_PROFILE=pnas writes figures/pnas/
fig = fk.figure('double', height=140, profile=PROFILE)
gs = fig.add_gridspec(3, 3, width_ratios=[1, 1, 1.05], height_ratios=[1, 1, 0.78])

# a hero: claims per gene by the share of lines that need the gene (Marimekko: width = gene share, height = claim share / gene share)
ax = fig.add_subplot(gs[0:2, 0:2]); fk.panel(ax, 'a', 'Genes with more dependent lines attract about twice their share of claims')
edges = [0, 1e-9, 0.01, 0.05, 0.2, 0.5, 0.9, 1.01]; labs = ['no line', '<1%', '1–5%', '5–20%', '20–50%', '50–90%', '≥90%']
bin_of = lambda f: np.digitize(f, edges[1:-1], right=False)
gb = np.bincount(bin_of(G.frac.values), minlength=7); gsh = gb / gb.sum()
cb = bin_of(A.frac_all.values)
pidx = pd.factorize(A.pmid)[0]; Kp = np.zeros((pidx.max() + 1, 7)); np.add.at(Kp, (pidx, cb), 1.0)
csh = Kp.sum(0) / Kp.sum(); ratio = csh / gsh
bs = np.array([(lambda k: k / k.sum())(np.bincount(rng.integers(0, len(Kp), len(Kp)), minlength=len(Kp)) @ Kp) / gsh for _ in range(2000)])
lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
xs = np.arange(7)
for i in range(7):
    up = ratio[i] > 1
    ax.bar(xs[i], ratio[i], width=0.72, color=fk.c['vermilion@t0.6'] if up else GL, edgecolor=VER if up else GL, lw=0.9 if up else 0, zorder=2)
    ax.vlines(xs[i], lo[i], hi[i], color=GD, lw=0.8, zorder=3)
    ax.text(xs[i], hi[i] + 0.05 if up else ratio[i] - 0.08, f'{ratio[i]:.1f}×', ha='center', va='bottom' if up else 'top', fontsize=FB, color=VER if up else GD,
            weight='bold' if up else None)
    ax.text(xs[i], -0.30, f'{gsh[i] * 100:.1f}', ha='center', va='center', fontsize=FT, color=GR)
    ax.text(xs[i], -0.47, f'{csh[i] * 100:.1f}', ha='center', va='center', fontsize=FT, color=VER if ratio[i] > 1 else GD)
ax.text(-0.5, -0.30, '% of genes', ha='right', va='center', fontsize=FT, color=GR)
ax.text(-0.5, -0.47, '% of claims', ha='right', va='center', fontsize=FT, color=VER)
ax.axhline(1, color=GD, lw=0.7, ls=(0, (3, 2)), zorder=4)
ax.text(-0.45, 1.04, 'even spread', ha='left', va='bottom', fontsize=FB, color=GD)
need = cb >= 2; gneed = (G.frac >= 0.01)
fk.numbox(ax, -0.42, 2.25, f'{need.mean() * 100:.0f}% of claims', color=VER)
ax.text(-0.42, 1.85, f'go to the {gneed.mean() * 100:.0f}% of\ngenes essential in\n≥ 1% of lines', ha='left', va='center', fontsize=FB, color=GD)
ax.set_xlim(-0.6, 6.6); ax.set_ylim(-0.74, 2.8); ax.set_xticks(xs); ax.set_xticklabels(labs); ax.set_yticks([0, 0.5, 1, 1.5, 2])
fk.chart(ax, grid='y', unit='Claims per gene, relative to an even spread'); ax.spines['bottom'].set_position(('data', 0))
ax.tick_params(axis='x', pad=1)
ax.text(3, -0.66, 'Share of CRISPR-screened lines that depend on the gene', ha='center', va='center', fontsize=FB, color=GD)

# b gene classes
ax = fig.add_subplot(gs[0, 2]); fk.panel(ax, 'b', 'Gene classes')
s = CC.pivot(index='gclass', columns='population', values='share').reindex(['common', 'selective', 'never'])
sup = [R.loc[f'S support | gclass={k}', 'est'] * 100 for k in s.index]
fk.bartable(ax, ['Common essential', 'Selective', 'Never essential'], [
    {'name': 'Genes', 'values': s.Genes.values * 100, 'fmt': '{:.1f}', 'bar': 'right', 'color': 'grey_light'},
    {'name': 'Claims', 'values': s.Claims.values * 100, 'fmt': '{:.1f}', 'bar': 'right', 'color': 'grey'},
    {'name': 'Supported', 'values': sup, 'fmt': '{:.0f}', 'bar': 'right', 'color': 'grey_dark'}], row_header='Gene class (%)')

# c lineage base rate
ax = fig.add_subplot(gs[1, 2]); fk.panel(ax, 'c', 'Support follows how widely the lineage depends on the gene')
nc = A[(A.gclass != 'common') & A.frac_lin.notna()].copy()
nc['bin'] = pd.cut(nc.frac_lin, [-0.001, 0, 0.05, 0.1, 0.2, 0.35, 0.5, 1.0])
b = nc.groupby('bin', observed=True).agg(n=('sup', 'size'), obs=('sup', 'mean'), exp=('frac_lin', 'mean')).query('n >= 30')
ax.plot([0, 0.62], [0, 0.62], color=GR, lw=0.7, ls=(0, (3, 2)), zorder=1)
ax.text(0.64, 0.50, 'random\nline', ha='right', va='center', fontsize=FT, color=GR)
se = np.sqrt(b.obs * (1 - b.obs) / b.n)
ax.vlines(b.exp, b.obs - 1.96 * se, b.obs + 1.96 * se, color=GD, lw=0.8, zorder=2)
ax.scatter(b.exp, b.obs, s=np.clip(b.n / 30, 8, 50), color=GD, lw=0, zorder=3)
rho = nc[['frac_lin']].assign(s=nc.sup.astype(float)).corr('spearman').iloc[0, 1]
ax.text(0.64, 0.04, f'Spearman ρ = {rho:.2f}\n{len(nc):,} claims', ha='right', va='bottom', fontsize=FT, color=GD)
ax.set_xlim(-0.02, 0.66); ax.set_ylim(-0.02, 0.66); ax.set_xticks([0, 0.25, 0.5]); ax.set_yticks([0, 0.25, 0.5])
ax.set_xticklabels(['0', '25', '50%']); ax.set_yticklabels(['0', '25', '50%'])
fk.chart(ax, grid='y', unit='Claims supported')
ax.set_xlabel('Lineage lines dependent on the gene', fontsize=FB, color=GD)

# d expression-matched calibration
ax = fig.add_subplot(gs[2, :]); fk.panel(ax, 'd', 'Claimed genes are supported more often than other genes expressed at the same level in the same line')
x = D.dec.values + 1
ax.fill_between(x, 0, D.bg_matched * 100, color=GL, lw=0, zorder=1); ax.plot(x, D.bg_matched * 100, color=GR, lw=0.9, zorder=2)
ax.fill_between(x, D.lo * 100, D.hi * 100, color=GR, alpha=0.25, lw=0, zorder=3)
ax.plot(x, D.claim_sup * 100, color=GD, lw=1.3, zorder=4); ax.scatter(x, D.claim_sup * 100, s=12, color=GD, lw=0, zorder=5)
fk.endlabels(ax, [(D.claim_sup.iloc[-1] * 100, 'claimed genes', GD), (D.bg_matched.iloc[-1] * 100 - 1.2, 'other genes, same line', GR)])
r1 = R.loc['P1 ratio to background']
fk.numbox(ax, 1.2, 13, f'{r1.est:.2f}× ({r1.lo:.2f}–{r1.hi:.2f})', color=GD)
ax.set_xlim(0.7, 10.3); ax.set_ylim(0, 19); ax.set_xticks(range(1, 11)); ax.set_yticks([0, 5, 10, 15])
fk.chart(ax, grid='y', unit='% supported')
ax.set_xlabel('Expression decile of the gene in the named line (non-common genes; 1 = lowest)', fontsize=FB, color=GD)
import os
out = os.environ.get('FIG_DRAFT') or ROOT / ('figures/pnas' if PROFILE == 'pnas' else 'figures/main') / 'fig2_gene_concentration'
saved = fk.save(fig, out, strict=not os.environ.get('FIG_DRAFT')); print(saved.report)

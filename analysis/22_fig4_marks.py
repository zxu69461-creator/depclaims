"""Figure 4. Independent replication and expression in the named line mark the claims that screens support (confirmatory claims).
a, support by the number of papers making the same gene-line claim, for claims repeated by groups sharing no author and by the same groups
(95% CI by bootstrap over papers); b, the same contrast outside common-essential genes; c, expression of the gene in the named line;
d, support by retraction status, era and method."""
import os
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import figkit as fk
ROOT = Path(__file__).resolve().parents[1]
C = duckdb.sql(f"""select r.repeat_class, r.n_papers, c.* from '{ROOT}/db/parquet/claims_confirmatory_repeat_class.parquet' r
                   join '{ROOT}/db/parquet/claims_confirmatory.parquet' c using (pmid, gid, ach)""").df()
db = duckdb.connect(str(ROOT / 'db/depclaims.duckdb'), read_only=True)
IR = db.execute("select * from independent_repeats").df()
R = db.execute("select quantity, n, est, lo, hi from confirm_results where split = 'confirmatory'").df().set_index('quantity')
VER, BLU, GD, GR, GL = fk.c.vermilion, fk.c.blue, fk.c.grey_dark, fk.c.grey, fk.c.grey_light
e = fk.enc('depclaims'); RC = {k: e.kw('repeat', k)['color'] for k in ('Single', 'Shared-author repeat', 'Independent repeat')}
FB, FT = fk.pt('base'), fk.pt('tick')
rng = np.random.default_rng(20261103)

def boot(d, B=2000):
    k = d.groupby('pmid').sup.agg(['sum', 'size']).values
    bs = [(lambda w: (w @ k[:, 0]) / (w @ k[:, 1]))(np.bincount(rng.integers(0, len(k), len(k)), minlength=len(k))) for _ in range(B)]
    return d.sup.mean(), *np.percentile(bs, [2.5, 97.5])

PROFILE = __import__('os').environ.get('FIG_PROFILE', 'natcomms')     # journal profile; FIG_PROFILE=pnas writes figures/pnas/
fig = fk.figure('double', height=150, profile=PROFILE)
gs = fig.add_gridspec(3, 3, width_ratios=[1, 1, 1.05], height_ratios=[1, 1, 0.72])

# a hero: dose-response in the number of papers
ax = fig.add_subplot(gs[0:2, 0:2]); fk.panel(ax, 'a', 'Claims repeated by independent groups are supported more often, rising with the number of papers')
C['nb'] = np.select([C.n_papers == 1, C.n_papers == 2, C.n_papers <= 4], [1, 2, 3], 4)
single = boot(C[C.repeat_class == 'single'])
ax.axhspan(single[1] * 100, single[2] * 100, color=GL, lw=0, zorder=0)
ax.axhline(single[0] * 100, color=GR, lw=0.8, zorder=1)
ax.vlines(1, single[1] * 100, single[2] * 100, color=GR, lw=1.0, zorder=3); ax.scatter([1], [single[0] * 100], s=60, color=GR, lw=0, zorder=4)
ax.text(1, single[2] * 100 + 1.2, f'single paper\n{single[0] * 100:.1f}%', ha='center', va='bottom', fontsize=FB, color=GR)
for cls, col, dx, mk in (('independent repeat', VER, 0.06, 'o'), ('shared-author repeat', GD, -0.06, 's')):
    d = C[C.repeat_class == cls]
    pts = [(b, *boot(g), len(g), g[['gid', 'ach']].drop_duplicates().shape[0]) for b, g in d.groupby('nb') if len(g) >= 30]
    P = pd.DataFrame(pts, columns=['nb', 'est', 'lo', 'hi', 'n', 'pairs'])
    xs = P.nb + dx
    ax.vlines(xs, P.lo * 100, P.hi * 100, color=col, lw=1.2 if col == VER else 0.9, zorder=3)
    ax.plot(xs, P.est * 100, color=col, lw=1.8 if col == VER else 1.0, zorder=4)
    ax.scatter(xs, P.est * 100, s=np.clip(P.n / 6, 14, 90), color=col, marker=mk, lw=0, zorder=5)
    for x, r in zip(xs, P.itertuples()):
        if col == VER: ax.text(x + 0.09, r.est * 100, f'{r.pairs:,} pairs', ha='left', va='center', fontsize=FT, color=col)
        else: ax.text(x, r.lo * 100 - 0.6, f'{r.pairs:,} pairs', ha='center', va='top', fontsize=FT, color=col)
    lab = 'independent groups' if col == VER else 'the same group'
    r = P.iloc[-2] if col == VER else P.iloc[0]
    ax.text(r.nb + dx - 0.1, r.est * 100 + (2.6 if col == VER else 0), lab, ha='right',
            va='center', fontsize=FB, color=col, weight='bold' if col == VER else None)
ir = IR[(IR.analysis == 'S support | repeat class') & (IR.level == 'independent repeat')].iloc[0]
fk.numbox(ax, 0.75, 41, f'independent repeats {ir.est * 100:.1f}% ({ir.lo * 100:.1f}–{ir.hi * 100:.1f}) vs single {single[0] * 100:.1f}%', color=VER)
ax.set_xlim(0.6, 4.7); ax.set_ylim(0, 44); ax.set_yticks([0, 10, 20, 30, 40]); ax.set_xticks([1, 2, 3, 4]); ax.set_xticklabels(['1', '2', '3–4', '≥5'])
fk.chart(ax, grid='y', unit='% of claims supported by CRISPR')
ax.set_xlabel('Papers making the same gene–line claim', fontsize=FB, color=GD)
fk.pin_inside(ax)

# b non-common genes
ax = fig.add_subplot(gs[0, 2]); fk.panel(ax, 'b', 'The gap widens outside common-essential genes')
T = IR.pivot(index='level', columns='analysis', values='est') * 100
T.columns = ['all', 'nc'] if T.columns[0].endswith('class') else ['nc', 'all']
labs = {'independent repeat': 'independent groups', 'single': 'single paper', 'shared-author repeat': 'same group'}
def spread(v, gap=1.6):
    # label positions at least `gap` apart, in the order of the values
    o = np.argsort(v); p = np.array(v, float)[o]
    for i in range(1, len(p)): p[i] = max(p[i], p[i - 1] + gap)
    out = np.empty(len(p)); out[o] = p - (p - np.array(v, float)[o]).mean(); return out
ks = ['single', 'shared-author repeat', 'independent repeat']
ya, yb = spread(T.loc[ks, 'all'].values), spread(T.loc[ks, 'nc'].values)
for (k, c), la, lb in zip((('single', GR), ('shared-author repeat', GD), ('independent repeat', VER)), ya, yb):
    a_, b_ = T.loc[k, 'all'], T.loc[k, 'nc']
    ax.plot([0, 1], [a_, b_], color=c, lw=1.8 if c == VER else 1.0, zorder=3); ax.scatter([0, 1], [a_, b_], s=18, color=c, lw=0, zorder=4)
    ax.text(-0.06, la, f'{a_:.1f}', ha='right', va='center', fontsize=FT, color=c)
    ax.text(1.06, lb, f'{b_:.1f}  {labs[k]}', ha='left', va='center', fontsize=FT, color=c, weight='bold' if c == VER else None)
ri, rn = T.loc['independent repeat'] / T.loc['single']
ax.text(0.5, 3, f'independent / single: {ri:.2f}× → {rn:.2f}×', ha='center', va='center', fontsize=FT, color=VER)
ax.set_xlim(-0.35, 2.2); ax.set_ylim(0, 28); ax.set_xticks([0, 1]); ax.set_xticklabels(['all genes', 'non-common\ngenes'], fontsize=FT); ax.set_yticks([0, 10, 20])
fk.chart(ax, grid='y', unit='% supported'); ax.spines['bottom'].set_bounds(0, 1)

# c expression
ax = fig.add_subplot(gs[1, 2]); fk.panel(ax, 'c', 'Support requires expression in the named line')
Ce = C[C.expr.notna()].copy()
Ce['eb'] = pd.cut(Ce.expr, [-0.01, 1, 2, 3, 4, 5, 6, 7, 20], labels=['<1', '1–2', '2–3', '3–4', '4–5', '5–6', '6–7', '≥7'])
b = Ce.groupby('eb', observed=True).agg(n=('sup', 'size'), k=('sup', 'sum')).reset_index()
xs = np.arange(len(b))
ax.bar(xs, b.n, width=0.72, color=GL, lw=0, zorder=2); ax.bar(xs, b.k, width=0.72, color=GD, lw=0, zorder=3)
for x, r in zip(xs, b.itertuples()):
    ax.text(x, r.n + 60, f'{r.k / r.n * 100:.0f}%' if r.k / r.n >= 0.01 else f'{r.k / r.n * 100:.1f}%', ha='center', va='bottom', fontsize=FT, color=GD)
ax.text(len(b) - 0.6, 2900, 'dark: supported', ha='right', va='center', fontsize=FT, color=GD)
ax.set_xticks(xs); ax.set_xticklabels(b.eb.astype(str), fontsize=FT); ax.set_xlim(-0.6, len(b) - 0.4); ax.set_ylim(0, 3100)
fk.chart(ax, grid='y', unit='Claims')
ax.set_xlabel('Expression in the named line, log$_{2}$(TPM + 1)', fontsize=FB, color=GD)

# d strata
ax = fig.add_subplot(gs[2, :]); fk.panel(ax, 'd', 'Claims from retracted papers are supported less often than claims from other papers')
groups = [('Retraction', [('retracted=0', 'not retracted'), ('retracted=1', 'retracted')]),
          ('Era', [('era=<=2009', '≤2009'), ('era=2010-14', '2010–14'), ('era=2015-19', '2015–19'), ('era=2020+', '2020–26')]),
          ('Method', [('meth=RNAi', 'siRNA, shRNA'), ('meth=CRISPR', 'CRISPR'), ('meth=other', 'unspecified')])]
x = 0; ticks = []; tl = []
for gname, items in groups:
    x0 = x
    for key, lab in items:
        r = R.loc[f'S support | {key}']
        hl = key == 'retracted=1'
        ax.bar(x, r.est * 100, width=0.66, color=fk.c['vermilion@t0.6'] if hl else GL, lw=0, zorder=2)
        ax.vlines(x, r.lo * 100, r.hi * 100, color=VER if hl else GD, lw=0.8, zorder=3)
        ax.text(x, r.hi * 100 + 0.8, f'{r.est * 100:.1f}', ha='center', va='bottom', fontsize=FT, color=VER if hl else GD)
        ticks.append(x); tl.append(f'{lab}\n({int(r.n):,})'); x += 1
    ax.text((x0 + x - 1) / 2, 30, gname, ha='center', va='bottom', fontsize=FB, color=GD, weight='bold')
    x += 0.7
ax.set_xticks(ticks); ax.set_xticklabels(tl, fontsize=FT); ax.set_xlim(-1.2, x - 0.6); ax.set_ylim(0, 34); ax.set_yticks([0, 10, 20])
fk.chart(ax, grid='y', unit='% supported')

out = os.environ.get('FIG_DRAFT') or ROOT / ('figures/pnas' if PROFILE == 'pnas' else 'figures/main') / 'fig4_replication_expression'
saved = fk.save(fig, out, strict=not os.environ.get('FIG_DRAFT')); print(saved.report)

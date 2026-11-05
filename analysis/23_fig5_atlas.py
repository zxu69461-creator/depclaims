"""Figure 5. An attention-dependency map of 18,530 CRISPR-screened genes (all claims). a, papers with a growth claim on the gene against the
share of screened lines that depend on it; b, the most-claimed genes; c, share of genes with any claim by dependent fraction; d, claims per
lineage by gene class, with the supported share."""
import os
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import figkit as fk
ROOT = Path(__file__).resolve().parents[1]
rp = lambda f: duckdb.sql(f"select * from '{ROOT / f}'").df()
G = rp('work/fig17/out/genes.parquet'); A = rp('work/fig17/out/claims_base.parquet')
VER, BLU, GD, GR, GL = fk.c.vermilion, fk.c.blue, fk.c.grey_dark, fk.c.grey, fk.c.grey_light
e = fk.enc('depclaims'); KC = {k: e.kw('gclass', k)['color'] for k in ('Never', 'Selective', 'Common')}
FB, FT = fk.pt('base'), fk.pt('tick')
rng = np.random.default_rng(7)
PROFILE = __import__('os').environ.get('FIG_PROFILE', 'natcomms')     # journal profile; FIG_PROFILE=pnas writes figures/pnas/
fig = fk.figure('double', height=140, profile=PROFILE)
gs = fig.add_gridspec(3, 3, width_ratios=[1, 1, 1.1], height_ratios=[1.0, 0.62, 1.05])

# a hero: attention vs dependency
ax = fig.add_subplot(gs[0:3, 0:2]); fk.panel(ax, 'a', 'Gene by gene, claims on broadly essential genes hold up')
fx = lambda f: np.where(f > 0, np.log10(np.clip(f, 1e-3, 1)), -3.55)            # genes no line needs sit in a strip left of the axis
G['x'] = fx(G.frac) + np.where(G.frac > 0, 0, rng.uniform(-0.12, 0.12, len(G)))
G['y'] = np.log10(G.n_papers + 1) + rng.uniform(-0.035, 0.035, len(G)) * (G.n_papers <= 3)
G['s'] = np.where(G.n_claims > 0, G.n_sup / G.n_claims.clip(lower=1), np.nan)
un = G[G.n_claims == 0]; cl = G[G.n_claims > 0]
ax.scatter(un.x, np.full(len(un), 0.06) + rng.uniform(-0.05, 0.05, len(un)), s=1.2, color=GL, lw=0, rasterized=True, zorder=1)
ax.scatter(cl.x, cl.y, s=3 + 0.8 * cl.n_papers.clip(upper=60), c=np.where(cl.s >= 0.5, VER, GR), alpha=np.where(cl.s >= 0.5, 0.8, 0.45), lw=0, rasterized=True, zorder=2)
ax.axvline(-3.3, color=GL, lw=0.6); ax.axhline(0.17, color=GL, lw=0.6)
nl = ax.text(-3.55, 1.72, 'no line', ha='center', va='bottom', fontsize=FT, color=GD)
ax.text(-0.15, 0.06, f'{len(un):,} genes without a claim', ha='right', va='center', fontsize=FT, color=GR)
box = fk.numbox(ax, -2.05, 2.24, f'{len(cl):,} genes · {len(A):,} claims · {A.ach.nunique()} lines', color=VER)
ax.set_xlim(-3.8, 0.12); ax.set_ylim(-0.03, 2.34)
# Labels of the named genes: each label sits next to its point on a short straight leader. Positions come from a small search in print
# millimetres (32 directions x 5 distances per label, then passes of coordinate descent) that forbids label-label and label-leader
# overlaps and crossing leaders, keeps labels clear of dots, the boxed number and the axes edges, and keeps leaders short and clear of
# other labelled points.
def place_labels(ax, labels, X, Y, sizes, avoid, vlines=(), hlines=(), soft_hlines=()):
    fig = ax.figure; fig.canvas.draw(); r = fig.canvas.get_renderer()
    mm = lambda xy: ax.transData.transform(xy) / fig.dpi * 25.4
    P = mm(np.c_[X, Y]); R = np.sqrt(sizes) / 2 / 72 * 25.4 + 0.25          # every dot: centre and radius (mm)
    ab = ax.get_window_extent(r); ax0, ay0, ax1, ay1 = (np.array([ab.x0, ab.y0, ab.x1, ab.y1]) / fig.dpi * 25.4)
    av = [np.array([b.x0, b.y0, b.x1, b.y1]) / fig.dpi * 25.4 + np.array([-0.8, -0.8, 0.8, 0.8]) for b in (t.get_window_extent(r) for t in avoid)]
    VL = [mm((v, 0))[0] for v in vlines]; HL = [mm((0, h))[1] for h in hlines]; SH = [mm((0, h))[1] for h in soft_hlines]
    T = []
    for sym, x, y, s, col in labels:
        t = ax.text(0, 0, sym, fontsize=FT); e = t.get_window_extent(r); t.remove()
        T.append((sym, mm((x, y)), np.sqrt(s) / 2 / 72 * 25.4, e.width / fig.dpi * 25.4 + 0.6, e.height / fig.dpi * 25.4 + 0.4, col))
    def cands(c, rad, w, h):
        out = []
        for d in (2.5, 4.5, 7.0, 10.0, 13.0):
            for th in np.deg2rad(np.arange(0, 360, 11.25)):
                ux, uy = np.cos(th), np.sin(th); ax_, ay_ = c[0] + (rad + d) * ux, c[1] + (rad + d) * uy
                ha = 'left' if ux > 0.38 else 'right' if ux < -0.38 else 'center'; va = 'bottom' if uy > 0.38 else 'top' if uy < -0.38 else 'center'
                x0 = ax_ if ha == 'left' else ax_ - w if ha == 'right' else ax_ - w / 2
                y0 = ay_ if va == 'bottom' else ay_ - h if va == 'top' else ay_ - h / 2
                out.append(dict(box=np.array([x0, y0, x0 + w, y0 + h]), a=np.array([ax_, ay_]), b=c + rad * np.array([ux, uy]), ha=ha, va=va, d=d))
        return out
    def hit_box(bx, P, R):
        cx = np.clip(P[:, 0], bx[0], bx[2]); cy = np.clip(P[:, 1], bx[1], bx[3])
        return (P[:, 0] - cx) ** 2 + (P[:, 1] - cy) ** 2 < R ** 2
    def seg_box(p, q, bx):
        for t in np.linspace(0, 1, 25):
            z = p + (q - p) * t
            if bx[0] <= z[0] <= bx[2] and bx[1] <= z[1] <= bx[3]: return True
        return False
    def cross(p1, p2, q1, q2):
        o = lambda a, b, c: np.sign((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        return o(p1, p2, q1) != o(p1, p2, q2) and o(q1, q2, p1) != o(q1, q2, p2)
    def seg_pts(p, q, P, R):
        v = q - p; L2 = max(v @ v, 1e-9); t = np.clip(((P - p) @ v) / L2, 0, 1); z = p + t[:, None] * v
        return ((P - z) ** 2).sum(1) < R ** 2
    C = [cands(c, rad, w, h) for _, c, rad, w, h, _ in T]
    centres = np.array([t[1] for t in T]); radii = np.array([t[2] for t in T])
    def own(i, k):
        o = C[i][k]; bx = o['box']
        if bx[0] < ax0 + 0.5 or bx[2] > ax1 - 0.5 or bx[1] < ay0 + 0.5 or bx[3] > ay1 - 0.5: return 1e6
        cost = 60 * hit_box(bx, P, R).sum() + 4 * seg_pts(o['b'], o['a'], P, R * 0.6).sum() + 3 * o['d']
        cost += 1e6 * sum(bx[0] - 0.4 < v < bx[2] + 0.4 for v in VL) + 1e6 * sum(bx[1] - 0.4 < h < bx[3] + 0.4 for h in HL)
        cost += 400 * sum(bx[1] - 0.3 < h < bx[3] + 0.3 for h in SH)
        cost += 5000 * sum(abs(o['a'][1] - h) < 0.6 and abs(o['b'][1] - h) < 0.6 for h in SH)          # leader lying along a gridline
        cost += sum(1e6 for z in av if not (bx[2] < z[0] or bx[0] > z[2] or bx[3] < z[1] or bx[1] > z[3]) or seg_box(o['b'], o['a'], z))
        oth = np.ones(len(T), bool); oth[i] = False
        cost += 2e3 * seg_pts(o['b'], o['a'], centres[oth], radii[oth] + 0.6).sum()
        cost += 2e3 * hit_box(bx, centres[oth], radii[oth] + 0.6).sum()
        return cost
    base = [np.array([own(i, k) for k in range(len(C[i]))]) for i in range(len(T))]
    def pair(i, k, j, l):
        a, b = C[i][k], C[j][l]; A_, B_ = a['box'], b['box']
        c = 1e6 if not (A_[2] < B_[0] or A_[0] > B_[2] or A_[3] < B_[1] or A_[1] > B_[3]) else 0
        c += 1e6 * (cross(a['b'], a['a'], b['b'], b['a']) or seg_box(a['b'], a['a'], B_) or seg_box(b['b'], b['a'], A_))
        return c
    sel = [int(np.argmin(c)) for c in base]
    for _ in range(12):
        changed = False
        for i in range(len(T)):
            tot = base[i] + np.array([sum(pair(i, k, j, sel[j]) for j in range(len(T)) if j != i) for k in range(len(C[i]))])
            k = int(np.argmin(tot))
            if k != sel[i]: sel[i], changed = k, True
        if not changed: break
    inv = ax.transData.inverted()
    for (sym, c, rad, w, h, col), i in zip(T, range(len(T))):
        o = C[i][sel[i]]; xy = inv.transform(o['a'] / 25.4 * fig.dpi); pt_ = inv.transform(c / 25.4 * fig.dpi)
        fk.callout(ax, tuple(pt_), sym, tuple(xy), ha=o['ha'], va=o['va'], color=col, fontsize=FT, dot=False)
    return sel
names = ['STAT3', 'TP53', 'BCL2', 'VEGFA', 'MYC', 'PLK1', 'EGFR', 'AURKA', 'SPP1', 'HIF1A', 'CXCR4', 'KRAS', 'ERBB2', 'YAP1']
L = G[G.symbol.isin(names)].set_index('symbol')
size = lambda n: 3 + 0.8 * min(n, 60)
axA = ax
ax.set_xticks([-3, -2, -1, 0]); ax.set_xticklabels(['0.1', '1', '10', '100%'])
ax.set_yticks(np.log10([2, 6, 11, 31, 101])); ax.set_yticklabels(['1', '5', '10', '30', '100'])
fk.chart(ax, grid='y', unit='Papers with a growth claim on the gene')
ax.set_xlabel('Share of CRISPR-screened lines that depend on the gene', fontsize=FB, color=GD)
fk.pin_inside(ax)

# b most-claimed genes
ax = fig.add_subplot(gs[0, 2]); fk.panel(ax, 'b', 'The most-claimed genes')
T = G.sort_values('n_papers', ascending=False).head(10)
fk.bartable(ax, list(T.symbol), [
    {'name': 'Papers', 'values': T.n_papers.values, 'fmt': '{:.0f}', 'bar': 'right', 'color': 'grey_light'},
    {'name': 'Lines dep. (%)', 'values': T.frac.values * 100, 'fmt': '{:.0f}', 'bar': 'right', 'color': 'grey'},
    {'name': 'Supported (%)', 'values': (T.n_sup / T.n_claims).values * 100, 'fmt': '{:.0f}', 'bar': 'right', 'color': 'vermilion@t0.4'}], row_header='Gene')

# c lineages by gene class
ax = fig.add_subplot(gs[2, 2]); fk.panel(ax, 'd', 'Claims by lineage and gene class')
lt = A.groupby('OncotreeLineage').agg(n=('sup', 'size'), s=('sup', 'mean'))
cmp_ = A.groupby(['OncotreeLineage', 'gclass']).size().unstack(fill_value=0)
lt = lt.join(cmp_.div(cmp_.sum(1), axis=0)).sort_values('n', ascending=False).head(8).iloc[::-1]
y = np.arange(len(lt)); nm = lambda l: {'Esophagus/Stomach': 'Oesophagus/stomach', 'CNS/Brain': 'CNS/brain', 'Ovary/Fallopian Tube': 'Ovary',
                                        'Bladder/Urinary Tract': 'Bladder', 'Bowel': 'Bowel', 'Head and Neck': 'Head and neck'}.get(l, l)
left = np.zeros(len(lt))
for k in ('never', 'selective', 'common'):
    ax.barh(y, lt[k], left=left, height=0.62, color=KC[k.capitalize()], lw=0); left += lt[k].values
for yi, r in zip(y, lt.itertuples()):
    ax.text(1.03, yi, f'{r.n:,}', ha='left', va='center', fontsize=FT, color=GR)
    ax.text(1.29, yi, f'{r.s * 100:.0f}%', ha='left', va='center', fontsize=FT, color=GD)
ax.text(1.03, len(lt) - 0.3, 'claims', ha='left', va='bottom', fontsize=FT, color=GR); ax.text(1.29, len(lt) - 0.3, 'supp.', ha='left', va='bottom', fontsize=FT, color=GD)
ax.set_yticks(y); ax.set_yticklabels([nm(l) for l in lt.index], fontsize=FT, color=GD); ax.tick_params(axis='y', length=0)
ax.set_xlim(0, 1.45); ax.set_ylim(-0.6, len(lt) + 0.5); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(['0', '50', '100%'])
for k, xx in (('never', 0.0), ('selective', 0.36), ('common', 0.72)):
    ax.text(xx, len(lt) - 0.3, k, ha='left', va='bottom', fontsize=FT, color=GD)
fk.chart(ax, grid=None); ax.spines['left'].set_visible(False); ax.spines['bottom'].set_bounds(0, 1)

# d open dependencies
ax = fig.add_subplot(gs[1, 2]); fk.panel(ax, 'c', 'Open map: selective dependencies of a fifth of lines')
sel = G[(G.gclass == 'selective') & (G.frac >= 0.2)]
k = int(round((sel.n_claims > 0).mean() * 100))
fk.waffle(ax, [k, 100 - k], [fk.c.grey_dark, fk.c.grey_light], ncol=20)
ax.text(1.04, 0.85, f'{(sel.n_claims > 0).sum():,} with a\ngrowth claim ({k}%)', transform=ax.transAxes, ha='left', va='center', fontsize=FT, color=GD)
ax.text(1.04, 0.2, f'{(sel.n_claims == 0).sum():,} without\none ({100 - k}%)', transform=ax.transAxes, ha='left', va='center', fontsize=FT, color=GR)

# label panel A last, once every panel exists and the layout is settled
place_labels(axA, [(sym, r.x, r.y, size(r.n_papers), VER if r.s >= 0.5 else GD) for sym, r in L.iterrows()],
             cl.x.values, cl.y.values, 3 + 0.8 * cl.n_papers.clip(upper=60).values, avoid=[box, nl],
             vlines=[-3.3], hlines=[0.17], soft_hlines=list(np.log10([2, 6, 11, 31, 101])))
out = os.environ.get('FIG_DRAFT') or ROOT / ('figures/pnas' if PROFILE == 'pnas' else 'figures/main') / 'fig5_attention_map'
saved = fk.save(fig, out, strict=not os.environ.get('FIG_DRAFT')); print(saved.report)

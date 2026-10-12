# =============================================================================
# Cohort fate map library (Fig 10.3.6 method, generalized)
# =============================================================================
# Description : Pooling + drawing for CARTA-based consensus fate maps, refactored from
#               Npp53_tracing/.../datadriven_fatemap/fate_map_datadriven_cohort.py (Fig 10.3.6) and
#               fate_map_datadriven_3panel.py (Fig 10.3.5) without changing their rules:
#               - terminal progenitor = the single smallest CARTA potency set per tumor (2026-09-07 rule)
#               - node size = fraction of tumors whose terminal set holds the state
#               - edge width = fraction of tumors whose terminal set holds both states (min 0.12)
#               - direction = pooled Mann-Whitney on within-tumor pseudotime rank + cross-tumor agreement
#                 (solid: p<0.01, |AUC-0.5|>=0.10, tumors agree; double: strong but tumors disagree;
#                 dashed: p<0.05, effect>=0.05; else headless)
#               - fade = <20 cells in every tumor; one dashed non-mesenchymal->mesenchymal bridge if none drawn
#               States, colours, positions, candidate edges and fallbacks come from a config dict, so the
#               original 12-state Npp53 configuration reproduces Fig 10.3.6 exactly (test_reproduce_fig10_3_6.py).
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/scvi/bin/python
# Validated   : 2026-10-11
# =============================================================================
import ast, itertools, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
from matplotlib.legend_handler import HandlerPatch
from scipy.stats import mannwhitneyu

MINC = 20   # fade threshold and per-tumor state threshold, as the original

# ── CARTA parsing and pooling ────────────────────────────────────────────────
def parse_progenitors(path, rename=None):
    rename = rename or {}
    out = []
    for ln in Path(path).read_text().splitlines():
        ln = ln.strip()
        if ln: out.append(frozenset(rename.get(x, x) for x in ast.literal_eval(ln)))
    return out

def tumor_summary(sets, order):
    """Original 2026-09-07 rule: terminal = min(sorted(sets, len desc), key=len); direct = root states in no other set."""
    srt = sorted(sets, key=len, reverse=True)
    root, inter = srt[0], srt[1:]
    term = min(srt, key=len)
    n_ties = sum(1 for s in srt if len(s) == len(term))
    direct = {x for x in root if not any(x in q for q in inter)}
    pairs = {tuple(sorted(p, key=order.index)) for p in itertools.combinations(term, 2)}
    comem = set()
    for q in inter:
        comem |= {tuple(sorted(p, key=order.index)) for p in itertools.combinations(q, 2)}
    return dict(terminal=set(term), direct=direct, cooccur=pairs, comember=comem, n_tied_smallest=n_ties)

def pool(summaries, arm_of, order):
    """summaries: {tumor: tumor_summary}; arm_of: {tumor: arm}. Returns the consensus JSON structure."""
    arms = sorted(set(arm_of.values()), key=['Intact', 'Cas', 'Cas+Enz'].index)
    out = dict(states=list(order), n_by_arm={a: 0 for a in arms},
               in_terminal_by_arm={a: {s: 0 for s in order} for a in arms},
               direct_from_root_by_arm={a: {s: 0 for s in order} for a in arms},
               cooccur_by_arm={a: {} for a in arms}, comember_by_arm={a: {} for a in arms})
    for t, sm in summaries.items():
        a = arm_of[t]; out['n_by_arm'][a] += 1
        for s in sm['terminal']: out['in_terminal_by_arm'][a][s] += 1
        for s in sm['direct']: out['direct_from_root_by_arm'][a][s] += 1
        for key, src in (('cooccur_by_arm', 'cooccur'), ('comember_by_arm', 'comember')):
            for x, y in sm[src]:
                k = f'{x}|{y}'; out[key][a][k] = out[key][a].get(k, 0) + 1
    return out

# ── k* selection (compute_best_k_v2.py logic) ────────────────────────────────
def objective_curve(res_dir, tumor, kmax=8):
    rows = []
    for k in range(1, kmax + 1):
        f = Path(res_dir) / f'{tumor}_results_{k}.txt'
        if f.exists():
            r = pd.read_csv(f, sep='\t').iloc[0]; rows.append((int(r['k']), float(r['objective_score'])))
    d = pd.DataFrame(rows, columns=['k', 'obj']).drop_duplicates('k').sort_values('k').reset_index(drop=True)
    d['is_plateau'] = False
    for i in range(1, len(d)):
        if d['obj'].iloc[i] >= d['obj'].iloc[i - 1] - 1e-9:
            d.loc[i:, 'is_plateau'] = True; break
    return d

def chord_elbow(d):
    u = d[~d['is_plateau']]
    xs, ys = u['k'].astype(float).values, u['obj'].astype(float).values
    if len(xs) < 3: return int(xs[0]) if len(xs) else np.nan
    xn = (xs - xs.min()) / max(1e-12, xs.max() - xs.min()); yn = (ys - ys.min()) / max(1e-12, ys.max() - ys.min())
    return int(xs[int(np.argmax(np.abs(xn + yn - 1) / np.sqrt(2)))])

# ── direction statistics (verbatim logic of pooled()) ─────────────────────────
class Directions:
    def __init__(self, pdf, states):
        """pdf columns: sample, state, pt_rank (within-tumor rank of the cells used)."""
        g = pdf.groupby(['sample', 'state']).agg(n=('pt_rank', 'size'), m=('pt_rank', 'mean')).reset_index()
        self.mpiv = g.pivot(index='sample', columns='state', values='m'); self.npiv = g.pivot(index='sample', columns='state', values='n')
        self.cells = {s: pdf[pdf.state == s]['pt_rank'].values for s in states}
        self.nsamp = {s: (pdf[pdf.state == s].groupby('sample').size() > 20).sum() for s in states}
    def both(self, a, b):
        if a not in self.npiv or b not in self.npiv: return []
        return [s for s in self.npiv.index if self.npiv.loc[s, a] > 20 and self.npiv.loc[s, b] > 20]
    def pooled(self, a, b):
        da, db = self.cells[a], self.cells[b]
        if len(da) < 5 or len(db) < 5 or self.nsamp[a] < 2 or self.nsamp[b] < 2: return None
        U, pv = mannwhitneyu(da, db, alternative='two-sided'); auc = U / (len(da) * len(db)); eff = abs(auc - 0.5)
        direc = (a, b) if np.median(da) <= np.median(db) else (b, a)
        bo = self.both(a, b); npair = len(bo)
        agree = (sum(1 for s in bo if ((self.mpiv.loc[s, a] < self.mpiv.loc[s, b]) == (direc == (a, b)))) / npair) if npair else None
        strong = pv < 0.01 and eff >= 0.10
        if strong and npair >= 2 and agree <= 0.5: return (direc, 'bidir')
        if strong and (npair < 2 or agree > 0.5): return (direc, 'solid')
        if pv < 0.05 and eff >= 0.05 and (npair < 2 or agree >= 0.5): return (direc, 'lean')
        return (direc, 'headless')
    def stats(self, a, b):
        da, db = self.cells[a], self.cells[b]
        if len(da) < 5 or len(db) < 5: return dict(p=np.nan, effect=np.nan, n_tumors_both=len(self.both(a, b)))
        U, pv = mannwhitneyu(da, db, alternative='two-sided')
        return dict(p=pv, effect=abs(U / (len(da) * len(db)) - 0.5), n_tumors_both=len(self.both(a, b)))

# ── drawing: cohort panel (Fig 10.3.6 sequence of calls, parameterized) ──────
class _HArrow(HandlerPatch):
    def create_artists(self, legend, oh, xd, yd, w, h, fs, trans):
        y = h / 2. - yd
        p = FancyArrowPatch((-xd, y), (-xd + w, y), arrowstyle=oh._astyle, color=oh._acol, lw=oh._alw, linestyle=oh._als, mutation_scale=9)
        p.set_transform(trans); return [p]

def _legarrow(color, ls, astyle, lw):
    h = FancyArrowPatch((0, 0), (1, 0), arrowstyle=astyle, color=color, lw=lw, linestyle=ls, mutation_scale=9)
    h._astyle = astyle; h._acol = color; h._als = ls; h._alw = lw; return h

def draw_cohort(cfg, data, faded, dirs, out_png, title, header, footer, legend_labels=None, edges_log=None):
    """cfg keys: colors (ordered dict state->hex), labels, index, pos, mes, canon (list of pairs), fallback (dict),
    init_state, bridge_order. data = pooled JSON (all arms summed). dirs = Directions."""
    SC, LB, SI, POS, MES = cfg['colors'], cfg['labels'], cfg['index'], cfg['pos'], cfg['mes']
    order = list(SC)
    N_TOT = sum(data['n_by_arm'].values())
    in_term = {s: sum(data['in_terminal_by_arm'][a].get(s, 0) for a in data['n_by_arm']) for s in SC}
    co = {}
    for a in data['n_by_arm']:
        for k, v in data['cooccur_by_arm'][a].items(): co[k] = co.get(k, 0) + v
    def cf(a, b):
        return co.get('|'.join(sorted([a, b], key=order.index)), 0) / max(1, N_TOT)
    def nrad(s): return 0.032 + 0.052 * (in_term.get(s, 0) / max(1, N_TOT))
    def circle(ax, s, r, fd):
        x, y = POS[s]; c = SC[s]; ec = c if not fd else '#bbb'; lw = 2.8 if not fd else 1.4; al = 1 if not fd else 0.35
        ax.add_patch(Circle((x, y), r, facecolor='white', edgecolor=ec, linewidth=lw, zorder=5, alpha=al))
        tx = ax.text(x, y + r * 0.14, LB[s], ha='center', va='center', fontsize=8.4, color='#111' if not fd else '#aaa', fontweight='bold', zorder=6, linespacing=0.9, alpha=1 if not fd else 0.55)
        if cfg.get('fit_labels'):   # shrink the label until it sits inside the circle (off for the original config)
            rend = fig.canvas.get_renderer(); fs = 8.4
            while fs > 4.0:
                bb = tx.get_window_extent(renderer=rend).transformed(ax.transData.inverted())
                if bb.width <= 1.70 * r and bb.height <= 1.45 * r: break
                fs -= 0.2; tx.set_fontsize(fs)
        by = y - r - 0.016
        ax.add_patch(FancyBboxPatch((x - 0.019, by - 0.015), 0.038, 0.028, boxstyle='round,pad=0.001,rounding_size=0.004', facecolor='white', edgecolor='#666' if not fd else '#ccc', linewidth=0.9, zorder=6, alpha=1 if not fd else 0.5))
        ax.text(x, by - 0.001, f'{SI[s]}', ha='center', va='center', fontsize=10, fontweight='bold', color='#111' if not fd else '#bbb', zorder=7)
    def edge(ax, src, dst, w, color, ls, rad, rs, rd, double=False, headless=False):
        if edges_log is not None:
            edges_log.append(dict(src=src, dst=dst, w=round(w, 3), color=color, double=double, headless=headless, dashed=ls != '-'))
        x0, y0 = POS[src]; x1, y1 = POS[dst]; rs += 0.004; rd += 0.009; dx, dy = x1 - x0, y1 - y0; L = (dx ** 2 + dy ** 2) ** .5
        if L < 1e-6: return
        style = '<|-|>' if double else ('-' if headless else '-|>')
        ax.add_patch(FancyArrowPatch((x0 + dx / L * rs, y0 + dy / L * rs), (x1 - dx / L * rd, y1 - dy / L * rd), arrowstyle=style, color=color, linewidth=1.0 + 3.2 * w, linestyle=ls, mutation_scale=18, connectionstyle=f'arc3,rad={rad}', zorder=3))
    node_r = {s: nrad(s) for s in SC}; avoid = list(SC)
    def cfree(src, dst):
        def bez(rad, n=15):
            x0, y0 = POS[src]; x1, y1 = POS[dst]; dx, dy = x1 - x0, y1 - y0; mx_, my_ = (x0 + x1) / 2, (y0 + y1) / 2; cx, cy = mx_ + rad * dy, my_ - rad * dx
            return [((1 - t / n) ** 2 * x0 + 2 * (t / n) * (1 - t / n) * cx + (t / n) ** 2 * x1, (1 - t / n) ** 2 * y0 + 2 * (t / n) * (1 - t / n) * cy + (t / n) ** 2 * y1) for t in range(1, n)]
        for rad in [0.0, 0.22, -0.22, 0.40, -0.40, 0.55, -0.55]:
            ok = True
            for nm in avoid:
                if nm in (src, dst): continue
                nx, ny = POS[nm]
                if any(((px - nx) ** 2 + (py - ny) ** 2) ** .5 < 0.08 for px, py in bez(rad)): ok = False; break
            if ok: return rad
        return 0.55
    def cfree2(src, dst):
        # new maps: route around every other circle by its actual radius + margin (original keeps cfree / rad 0.10)
        def bez(rad, n=40):
            x0, y0 = POS[src]; x1, y1 = POS[dst]; dx, dy = x1 - x0, y1 - y0; mx_, my_ = (x0 + x1) / 2, (y0 + y1) / 2; cx, cy = mx_ + rad * dy, my_ - rad * dx
            return [((1 - t / n) ** 2 * x0 + 2 * (t / n) * (1 - t / n) * cx + (t / n) ** 2 * x1, (1 - t / n) ** 2 * y0 + 2 * (t / n) * (1 - t / n) * cy + (t / n) ** 2 * y1) for t in range(1, n)]
        best, best_gap = 0.10, -1
        for rad in [0.10, -0.10, 0.0, 0.22, -0.22, 0.32, -0.32, 0.42, -0.42, 0.55, -0.55, 0.70, -0.70]:
            gap = min((((px - POS[nm][0]) ** 2 + (py - POS[nm][1]) ** 2) ** .5 - node_r[nm] - 0.035 for nm in avoid if nm not in (src, dst) for px, py in bez(rad)), default=1)
            if gap > 0: return rad
            if gap > best_gap: best, best_gap = rad, gap
        return best
    route = cfg.get('route_all_edges', False)
    init = cfg['init_state']
    fig, ax = plt.subplots(figsize=(11, 11)); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal'); ax.axis('off')
    ax.text(0.50, 0.955, 'Tumor initiation', ha='center', va='center', fontsize=11, fontweight='bold', color='#111', zorder=6)
    r1 = nrad(init); xi, yi = POS[init]
    ax.add_patch(FancyArrowPatch((xi, yi + r1 + 0.028), (xi, yi + r1 + 0.006), arrowstyle='-|>', color='#111', lw=2.0, mutation_scale=14, zorder=3))
    ax.add_patch(FancyBboxPatch((0.02, 0.965), 0.96, 0.045, boxstyle='round,pad=0.002,rounding_size=0.01', facecolor='#0f2747', edgecolor='none', alpha=0.12, zorder=1))
    ax.text(0.50, 0.9875, header, ha='center', va='center', fontsize=15, fontweight='bold', color='#0f2747')
    bridge_drawn = False; has_incoming = set(); connected = set()
    for a_, b_ in cfg['canon']:
        if a_ in faded or b_ in faded: continue
        has_incoming.add(b_); connected |= {a_, b_}
        pr = dirs.pooled(a_, b_); w = max(cf(a_, b_), 0.12); db = False
        if (a_, b_) == cfg.get('boost_pair'): w = max(w, 0.30)
        if pr is None or pr[1] == 'headless':
            edge(ax, a_, b_, min(w, 1), '#9a9a9a', (0, (4, 3)), cfree2(a_, b_) if route else cfree(a_, b_), node_r[a_], node_r[b_], headless=True); continue
        direc, lvl = pr
        if lvl == 'solid': col, ls = '#111', '-'
        elif lvl == 'bidir': col, ls, db = '#111', '-', True
        else: col, ls = '#7a86a0', (0, (4, 3))
        edge(ax, direc[0], direc[1], min(w, 1), col, ls, cfree2(direc[0], direc[1]) if route else 0.10, node_r[direc[0]], node_r[direc[1]], double=db)
        if (direc[0] not in MES and direc[1] in MES) or (direc[1] not in MES and direc[0] in MES): bridge_drawn = True
    present_mes = [s for s in MES if s not in faded]
    if not bridge_drawn and present_mes:
        target = max(present_mes, key=lambda s: len(dirs.cells[s])); up = init
        for c in cfg['bridge_order']:
            if c not in faded and c not in MES: up = c; break
        edge(ax, up, target, 0.20, '#7a86a0', (0, (4, 3)), cfree2(up, target) if route else cfree(up, target), node_r[up], node_r[target]); has_incoming.add(target); connected |= {up, target}
    for s in SC:
        if s in faded or s == init or s in has_incoming: continue
        if cfg.get('fallback_only_isolated') and s in connected: continue   # new maps: connectivity edge only for states with no edge at all
        up = init
        for c in cfg['fallback'].get(s, [init]):
            if c not in faded: up = c; break
        edge(ax, up, s, 0.10, '#c0c0c0', (0, (4, 3)), cfree2(up, s) if route else cfree(up, s), node_r[up], node_r[s], headless=True)
    axx = 0.985
    ax.add_patch(FancyArrowPatch((axx, 0.87), (axx, 0.05), arrowstyle='-|>', color='black', lw=1.6, mutation_scale=12, zorder=4))
    ax.text(axx + 0.012, 0.46, 'Trajectory', ha='center', va='center', fontsize=10.5, color='black', fontweight='bold', rotation=270)
    for s in SC: circle(ax, s, node_r[s], fd=(s in faded))
    fig.suptitle(title, fontsize=13.5, fontweight='bold', color='#0f2747', y=0.965)
    leg = [_legarrow('#111', '-', '-|>', 2.6), _legarrow('#111', '-', '<|-|>', 2.6),
           _legarrow('#7a86a0', (0, (4, 3)), '-|>', 2.2), _legarrow('#9a9a9a', (0, (4, 3)), '-', 1.8)]
    labels = legend_labels or ['confident direction — single arrow (pooled MWU p<0.01, effect≥0.10, tumors agree)',
                               'strong link, direction unresolved — double arrow ↔ (tumors disagree)',
                               'weak lean — dashed single arrow (p<0.05) / guaranteed mesenchymal bridge',
                               'headless dashed - - - (below threshold / unresolved / connectivity; no arrowhead)']
    ax.legend(handles=leg, labels=labels, handler_map={FancyArrowPatch: _HArrow()}, loc='lower center', bbox_to_anchor=(0.5, -0.10), ncol=2, fontsize=8.5, frameon=False, handlelength=3.4, handletextpad=0.7, columnspacing=2.4)
    fig.text(0.5, 0.005, footer, ha='center', fontsize=8.5, color='#333', style='italic')
    plt.subplots_adjust(top=0.92, bottom=0.10, left=0.02, right=0.98)
    out_png = Path(out_png)
    for fmt in ['png', 'pdf']:
        fig.savefig(out_png.with_suffix(f'.{fmt}'), dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return dict(node_size_fraction={s: in_term[s] / max(1, N_TOT) for s in SC}, n_tumors=N_TOT)


# ── drawing v2: same rules as draw_cohort, layout fixes for the new maps ─────
# Edge choice, classes, widths, bridge and connectivity rules are copied from draw_cohort unchanged.
# Layout only: (1) every edge is routed around all circles AND number boxes; (2) arrow ends that meet the
# same circle are kept >= MIN_SEP degrees apart (no stacked arrowheads); (3) ends sit where the curve
# actually leaves/enters the circle; (4) header band not clipped; "Tumor initiation" placed above the
# initiation node; (5) circle labels are shrunk to fit after the final layout.
def draw_cohort_v2(cfg, data, faded, dirs, out_png, title, header, footer, legend_labels, edges_log=None):
    import math
    SC, LB, SI, POS, MES = cfg['colors'], cfg['labels'], cfg['index'], cfg['pos'], cfg['mes']
    order = list(SC)
    N_TOT = sum(data['n_by_arm'].values())
    in_term = {s: sum(data['in_terminal_by_arm'][a].get(s, 0) for a in data['n_by_arm']) for s in SC}
    co = {}
    for a in data['n_by_arm']:
        for k, v in data['cooccur_by_arm'][a].items(): co[k] = co.get(k, 0) + v
    def cf(a, b): return co.get('|'.join(sorted([a, b], key=order.index)), 0) / max(1, N_TOT)
    K = cfg.get('node_scale', 1.0)   # uniform display scale of the original radius formula (relative sizes unchanged)
    def nrad(s): return K * (0.032 + 0.052 * (in_term.get(s, 0) / max(1, N_TOT)))
    node_r = {s: nrad(s) for s in SC}
    BOX_DY, BOX_R, MARGIN, MIN_SEP = 0.016, 0.027, 0.012, math.radians(16)
    obst = [(POS[s][0], POS[s][1], node_r[s], s, 'circle') for s in SC] + \
           [(POS[s][0], POS[s][1] - node_r[s] - BOX_DY, BOX_R, s, 'box') for s in SC]
    init = cfg['init_state']
    _xi, _yi = POS[init]; _ty = min(_yi + node_r[init] + 0.045, 0.932)
    obst += [(_xi + dx_, _ty, 0.020, '_init', 'label') for dx_ in (-0.06, -0.03, 0.0, 0.03, 0.06)] + \
            [(_xi, (_ty + _yi + node_r[init]) / 2, 0.012, '_init', 'label')]
    def bound(px, py): return min(0.935 - py, px - 0.012, 0.972 - px, py - 0.010)   # stay inside the drawing area
    # 1) decide edges exactly as draw_cohort
    specs = []   # (src, dst, w, color, ls, double, headless)
    bridge_drawn = False; has_incoming = set(); connected = set()
    for a_, b_ in cfg['canon']:
        if a_ in faded or b_ in faded: continue
        has_incoming.add(b_); connected |= {a_, b_}
        pr = dirs.pooled(a_, b_); w = max(cf(a_, b_), 0.12)
        if pr is None or pr[1] == 'headless':
            specs.append((a_, b_, min(w, 1), '#9a9a9a', (0, (4, 3)), False, True)); continue
        direc, lvl = pr
        if lvl == 'solid': specs.append((direc[0], direc[1], min(w, 1), '#111', '-', False, False))
        elif lvl == 'bidir': specs.append((direc[0], direc[1], min(w, 1), '#111', '-', True, False))
        else: specs.append((direc[0], direc[1], min(w, 1), '#7a86a0', (0, (4, 3)), False, False))
        if (direc[0] not in MES and direc[1] in MES) or (direc[1] not in MES and direc[0] in MES): bridge_drawn = True
    present_mes = [s for s in MES if s not in faded]
    if not bridge_drawn and present_mes:
        target = max(present_mes, key=lambda s: len(dirs.cells[s])); up = init
        for c in cfg['bridge_order']:
            if c not in faded and c not in MES: up = c; break
        specs.append((up, target, 0.20, '#7a86a0', (0, (4, 3)), False, False)); has_incoming.add(target); connected |= {up, target}
    for s in SC:
        if s in faded or s == init or s in has_incoming or s in connected: continue
        up = init
        for c in cfg['fallback'].get(s, [init]):
            if c not in faded: up = c; break
        specs.append((up, s, 0.10, '#c0c0c0', (0, (4, 3)), False, True))
    # 2) route each edge: clearance from circles/boxes + separation of arrow ends at shared circles
    used = {s: [] for s in SC}
    def geom(src, dst, rad, os_=0.0, od_=0.0):
        (x0, y0), (x1, y1) = POS[src], POS[dst]; dx, dy = x1 - x0, y1 - y0
        ts = math.atan2(dy / 2 - rad * dx, dx / 2 + rad * dy) + os_
        td = math.atan2(-(dy / 2 + rad * dx), -(dx / 2 - rad * dy)) + od_
        rs, rd = node_r[src] + 0.004, node_r[dst] + 0.009
        p0 = (x0 + rs * math.cos(ts), y0 + rs * math.sin(ts)); p1 = (x1 + rd * math.cos(td), y1 + rd * math.sin(td))
        ddx, ddy = p1[0] - p0[0], p1[1] - p0[1]; mx_, my_ = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
        c = (mx_ + rad * ddy, my_ - rad * ddx)
        pts = [((1 - t) ** 2 * p0[0] + 2 * t * (1 - t) * c[0] + t ** 2 * p1[0], (1 - t) ** 2 * p0[1] + 2 * t * (1 - t) * c[1] + t ** 2 * p1[1]) for t in np.linspace(0.02, 0.98, 60)]
        return ts, td, p0, p1, pts
    def angdist(a, b): d = abs(a - b) % (2 * math.pi); return min(d, 2 * math.pi - d)
    drawn, paths, heads, HEAD_CLEAR, problems = [], [], [], 0.012, []
    # route short edges first (fewest clean options); long arcs then bend around them
    specs = sorted(specs, key=lambda e: math.hypot(POS[e[1]][0] - POS[e[0]][0], POS[e[1]][1] - POS[e[0]][1]))
    for (src, dst, w, color, ls, double, headless) in specs:
        best = None
        RADS = [0.10, -0.10, 0.0, 0.18, -0.18, 0.28, -0.28, 0.38, -0.38, 0.50, -0.50, 0.65, -0.65, 0.80, -0.80, 0.95, -0.95]
        cands = [(r_, 0.0, 0.0) for r_ in RADS] + [(r_, o1, o2) for o1 in (0.0, 0.26, -0.26) for o2 in (0.0, 0.26, -0.26) if (o1, o2) != (0.0, 0.0) for r_ in RADS]
        for rad, os_, od_ in cands:
            ts, td, p0, p1, pts = geom(src, dst, rad, os_, od_)
            clear = min((math.hypot(px - ox, py - oy) - orad for px, py in pts for ox, oy, orad, s, kind in obst
                         if not (kind == 'circle' and s in (src, dst))), default=1.0)
            clear = min(clear, min(bound(px, py) for px, py in pts))
            sep = min([angdist(ts, a) for a in used[src]] + [angdist(td, a) for a in used[dst]] + [math.pi])
            ends = [] if headless else pts[-8:] + (pts[:8] if double else [])   # arrowhead region(s)
            head = min([math.hypot(px - qx, py - qy) for px, py in ends for q in paths for qx, qy in q] +      # my heads off other lines
                       [math.hypot(px - hx, py - hy) for px, py in pts for hx, hy in heads] + [1.0])        # my line off other heads
            ok_c, ok_h, ok_s = clear > MARGIN, head > HEAD_CLEAR, sep >= MIN_SEP
            par = _parallel_run(pts, paths)
            hard = (clear > 0.0, sep >= math.radians(5), head > 0.011, par <= 0.02)   # no touching, no overlapping ends, no arrowhead on a line, no side-by-side runs
            score = (hard, ok_c, ok_h, ok_s, 0 if ok_c else clear, 0 if ok_h else head, 0 if ok_s else sep, -abs(rad))
            if best is None or score > best[0]: best = (score, rad, ts, td, p0, p1, pts)
            if ok_c and ok_h and ok_s: break
        sc, rad, ts, td, p0, p1, pts = best
        if not all(sc[1:4]) or not all(sc[0]): problems.append((src, dst, dict(hard=sc[0], clear=sc[1], head=sc[2], sep=sc[3])))
        used[src].append(ts); used[dst].append(td); paths.append(pts)
        if not headless: heads.extend(pts[-8:] + (pts[:8] if double else []))
        drawn.append((src, dst, w, color, ls, double, headless, rad, p0, p1))
        if edges_log is not None:
            edges_log.append(dict(src=src, dst=dst, w=round(w, 3), color=color, double=double, headless=headless, dashed=ls != '-', rad=rad))
    # 3) draw
    fig, ax = plt.subplots(figsize=(11, 11)); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal'); ax.axis('off')
    r1 = node_r[init]; xi, yi = POS[init]
    ty = min(yi + r1 + 0.045, 0.932)   # stay below the header band
    ax.text(xi, ty, 'Tumor initiation', ha='center', va='center', fontsize=11, fontweight='bold', color='#111', zorder=6)
    if ty - 0.014 > yi + r1 + 0.010:
        ax.add_patch(FancyArrowPatch((xi, ty - 0.014), (xi, yi + r1 + 0.006), arrowstyle='-|>', color='#111', lw=2.0, mutation_scale=14, zorder=3))
    ax.add_patch(FancyBboxPatch((0.02, 0.955), 0.96, 0.042, boxstyle='round,pad=0.002,rounding_size=0.01', facecolor='#0f2747', edgecolor='none', alpha=0.12, zorder=1, clip_on=False))
    ax.text(0.50, 0.976, header, ha='center', va='center', fontsize=15, fontweight='bold', color='#0f2747')
    for (src, dst, w, color, ls, double, headless, rad, p0, p1) in drawn:
        style = '<|-|>' if double else ('-' if headless else '-|>')
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, color=color, linewidth=1.0 + 3.2 * w, linestyle=ls, mutation_scale=18, connectionstyle=f'arc3,rad={rad}', zorder=3))
    axx = 0.985
    ax.add_patch(FancyArrowPatch((axx, 0.87), (axx, 0.05), arrowstyle='-|>', color='black', lw=1.6, mutation_scale=12, zorder=4))
    ax.text(axx + 0.012, 0.46, 'Trajectory', ha='center', va='center', fontsize=10.5, color='black', fontweight='bold', rotation=270)
    labels_to_fit = []
    for s in SC:
        x, y = POS[s]; r = node_r[s]; fd = s in faded; c = SC[s]
        ax.add_patch(Circle((x, y), r, facecolor='white', edgecolor=c if not fd else '#bbb', linewidth=2.8 if not fd else 1.4, zorder=5, alpha=1 if not fd else 0.35))
        tx = ax.text(x, y + r * 0.10, LB[s], ha='center', va='center', fontsize=8.4, color='#111' if not fd else '#8a8a8a', fontweight='bold', zorder=6, linespacing=0.9)
        labels_to_fit.append((tx, r))
        by = y - r - BOX_DY
        ax.add_patch(FancyBboxPatch((x - 0.019, by - 0.015), 0.038, 0.028, boxstyle='round,pad=0.001,rounding_size=0.004', facecolor='white', edgecolor='#666' if not fd else '#ccc', linewidth=0.9, zorder=6, alpha=1 if not fd else 0.5))
        ax.text(x, by - 0.001, f'{SI[s]}', ha='center', va='center', fontsize=10, fontweight='bold', color='#111' if not fd else '#999', zorder=7)
    fig.suptitle(title, fontsize=13.5, fontweight='bold', color='#0f2747', y=0.965)
    leg = [_legarrow('#111', '-', '-|>', 2.6), _legarrow('#111', '-', '<|-|>', 2.6),
           _legarrow('#7a86a0', (0, (4, 3)), '-|>', 2.2), _legarrow('#9a9a9a', (0, (4, 3)), '-', 1.8)]
    ax.legend(handles=leg, labels=legend_labels, handler_map={FancyArrowPatch: _HArrow()}, loc='lower center', bbox_to_anchor=(0.5, -0.10), ncol=2, fontsize=8.5, frameon=False, handlelength=3.4, handletextpad=0.7, columnspacing=2.4)
    fig.text(0.5, 0.005, footer, ha='center', fontsize=8.5, color='#333', style='italic')
    plt.subplots_adjust(top=0.92, bottom=0.10, left=0.02, right=0.98)
    fig.canvas.draw()   # apply the equal-aspect squeeze before measuring labels
    rend = fig.canvas.get_renderer()
    lw_data = (2.8 / 72) / (ax.get_window_extent(renderer=rend).height / fig.dpi)   # ring line width in data units
    for tx, r in labels_to_fit:   # fit after the final layout: every corner of the text box inside the ring's inner edge
        fs = 8.4; r_in = r - lw_data - 0.003
        ccx, ccy = tx.get_position()[0], tx.get_position()[1] - r * 0.10   # circle centre
        while fs > 4.0:
            bb = tx.get_window_extent(renderer=rend).transformed(ax.transData.inverted())
            far = max(((px - ccx) ** 2 + (py - ccy) ** 2) ** .5 for px in (bb.x0, bb.x1) for py in (bb.y0, bb.y1))
            if far <= r_in: break
            fs -= 0.2; tx.set_fontsize(fs)
    out_png = Path(out_png)
    for fmt in ['png', 'pdf']:
        fig.savefig(out_png.with_suffix(f'.{fmt}'), dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return dict(node_size_fraction={s: in_term[s] / max(1, N_TOT) for s in SC}, n_tumors=N_TOT, routing_problems=problems,
                min_label_pt={tx.get_text().replace(chr(10), ' '): round(tx.get_fontsize(), 1) for tx, r in labels_to_fit})



def _parallel_run(pts, paths, near=0.006):
    """Longest stretch (plot units) of the interior of pts that stays within `near` of any already-drawn path."""
    import math
    if not paths: return 0.0
    inner = pts[6:-6]
    close = [min(math.hypot(px - qx, py - qy) for q in paths for qx, qy in q) < near for px, py in inner]
    best = run = 0.0
    for i in range(1, len(inner)):
        seg = math.hypot(inner[i][0] - inner[i - 1][0], inner[i][1] - inner[i - 1][1])
        run = run + seg if (close[i] and close[i - 1]) else 0.0
        best = max(best, run)
    return best

# ── per-treatment panels v2 (Fig 10.3.5 rules, layout as draw_cohort_v2) ─────
# Rules copied from fate_map_datadriven_3panel.py: per-arm node size (0.032 + 0.048 x terminal fraction),
# per-arm edge width (co-occurrence / n_arm, min 0.12), per-arm direction statistics, states with <20 cells
# in every tumor of the arm are faded and get no links, below-threshold (headless) edges are NOT drawn, and
# one dashed non-mesenchymal -> mesenchymal bridge is added if none is drawn. Same candidate scaffold for
# every arm (the cohort map's list), as the original used one list for all arms.
def _route(specs, POS, node_r, SC, MARGIN=0.012, MIN_SEP_DEG=16, HEAD_CLEAR=0.012, BOX_DY=0.016, BOX_R=0.027, extra=()):
    import math
    MIN_SEP = math.radians(MIN_SEP_DEG)
    obst = [(POS[s][0], POS[s][1], node_r[s], s, 'circle') for s in SC] + \
           [(POS[s][0], POS[s][1] - node_r[s] - BOX_DY, BOX_R, s, 'box') for s in SC] + list(extra)
    def bound(px, py): return min(0.935 - py, px - 0.012, 0.988 - px, py - 0.010)   # stay inside the panel
    used = {s: [] for s in SC}
    def geom(src, dst, rad, os_=0.0, od_=0.0):
        (x0, y0), (x1, y1) = POS[src], POS[dst]; dx, dy = x1 - x0, y1 - y0
        ts = math.atan2(dy / 2 - rad * dx, dx / 2 + rad * dy) + os_
        td = math.atan2(-(dy / 2 + rad * dx), -(dx / 2 - rad * dy)) + od_
        rs, rd = node_r[src] + 0.004, node_r[dst] + 0.009
        p0 = (x0 + rs * math.cos(ts), y0 + rs * math.sin(ts)); p1 = (x1 + rd * math.cos(td), y1 + rd * math.sin(td))
        ddx, ddy = p1[0] - p0[0], p1[1] - p0[1]; mx_, my_ = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
        c = (mx_ + rad * ddy, my_ - rad * ddx)
        pts = [((1 - t) ** 2 * p0[0] + 2 * t * (1 - t) * c[0] + t ** 2 * p1[0], (1 - t) ** 2 * p0[1] + 2 * t * (1 - t) * c[1] + t ** 2 * p1[1]) for t in np.linspace(0.02, 0.98, 60)]
        return ts, td, p0, p1, pts
    def angdist(a, b): d = abs(a - b) % (2 * math.pi); return min(d, 2 * math.pi - d)
    drawn, paths, heads, problems = [], [], [], []
    specs = sorted(specs, key=lambda e: math.hypot(POS[e['dst']][0] - POS[e['src']][0], POS[e['dst']][1] - POS[e['src']][1]))   # short edges first
    for sp in specs:
        src, dst, double = sp['src'], sp['dst'], sp['double']
        best = None
        RADS = [0.10, -0.10, 0.0, 0.18, -0.18, 0.28, -0.28, 0.38, -0.38, 0.50, -0.50, 0.65, -0.65, 0.80, -0.80, 0.95, -0.95]
        cands = [(r_, 0.0, 0.0) for r_ in RADS] + [(r_, o1, o2) for o1 in (0.0, 0.26, -0.26) for o2 in (0.0, 0.26, -0.26) if (o1, o2) != (0.0, 0.0) for r_ in RADS]
        for rad, os_, od_ in cands:
            ts, td, p0, p1, pts = geom(src, dst, rad, os_, od_)
            clear = min((math.hypot(px - ox, py - oy) - orad for px, py in pts for ox, oy, orad, s, kind in obst
                         if not (kind == 'circle' and s in (src, dst))), default=1.0)
            clear = min(clear, min(bound(px, py) for px, py in pts))
            sep = min([angdist(ts, a) for a in used[src]] + [angdist(td, a) for a in used[dst]] + [math.pi])
            ends = pts[-8:] + (pts[:8] if double else [])
            head = min([math.hypot(px - qx, py - qy) for px, py in ends for q in paths for qx, qy in q] +
                       [math.hypot(px - hx, py - hy) for px, py in pts for hx, hy in heads] + [1.0])
            ok_c, ok_h, ok_s = clear > MARGIN, head > HEAD_CLEAR, sep >= MIN_SEP
            par = _parallel_run(pts, paths)
            hard = (clear > 0.0, sep >= math.radians(5), head > 0.011, par <= 0.02)   # no touching, no overlapping ends, no arrowhead on a line, no side-by-side runs
            score = (hard, ok_c, ok_h, ok_s, 0 if ok_c else clear, 0 if ok_h else head, 0 if ok_s else sep, -abs(rad))
            if best is None or score > best[0]: best = (score, rad, ts, td, p0, p1, pts)
            if ok_c and ok_h and ok_s: break
        sc, rad, ts, td, p0, p1, pts = best
        if not all(sc[0]) or not all(sc[1:4]): problems.append((src, dst, dict(hard=sc[0], clear=sc[1], head=sc[2], sep=sc[3])))
        used[src].append(ts); used[dst].append(td); paths.append(pts); heads.extend(pts[-8:] + (pts[:8] if double else []))
        drawn.append(dict(sp, rad=rad, p0=p0, p1=p1))
    return drawn, problems

def _fit_labels(fig, ax, labels_to_fit, start=8.4, lw_pt=2.6):
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    lw_data = (lw_pt / 72) / (ax.get_window_extent(renderer=rend).height / fig.dpi)
    out = []
    for tx, r, ccx, ccy in labels_to_fit:
        fs = start; r_in = r - lw_data - 0.003; tx.set_fontsize(fs)
        while fs > 4.0:
            bb = tx.get_window_extent(renderer=rend).transformed(ax.transData.inverted())
            if max(((px - ccx) ** 2 + (py - ccy) ** 2) ** .5 for px in (bb.x0, bb.x1) for py in (bb.y0, bb.y1)) <= r_in: break
            fs -= 0.2; tx.set_fontsize(fs)
        out.append(fs)
    return out

def draw_arms_v2(cfg, data, faded_by_arm, dirs_by_arm, out_png, title, footer, legend_labels, arm_display, arm_accent, edges_log=None, arm_text=None, hidden_by_arm=None):
    arm_text = arm_text or arm_accent
    hidden_by_arm = hidden_by_arm or {}
    SC, LB, SI, POS, MES = cfg['colors'], cfg['labels'], cfg['index'], cfg['pos'], cfg['mes']
    order = list(SC); K = cfg.get('node_scale', 1.0); init = cfg['init_state']
    arms = [a for a in ['Intact', 'Cas', 'Cas+Enz'] if data['n_by_arm'].get(a, 0) > 0]
    fig, axes = plt.subplots(1, len(arms), figsize=(8 * len(arms), 9.2), squeeze=False); axes = axes[0]
    all_labels, min_pt, problems = [], {}, {}
    for ax, arm in zip(axes, arms):
        n_arm = data['n_by_arm'][arm]; in_term = data['in_terminal_by_arm'][arm]
        co = {k: v for k, v in data['cooccur_by_arm'][arm].items()}
        def cf(a, b): return co.get('|'.join(sorted([a, b], key=order.index)), 0) / max(1, n_arm)
        node_r = {s: K * (0.032 + 0.048 * (in_term.get(s, 0) / max(1, n_arm))) for s in SC}
        hidden = set(hidden_by_arm.get(arm, set())); faded = set(faded_by_arm[arm]) | hidden; dirs = dirs_by_arm[arm]
        specs = []; bridge_drawn = False
        for a_, b_ in cfg['canon']:
            if a_ in faded or b_ in faded: continue
            pr = dirs.pooled(a_, b_)
            if pr is None or pr[1] == 'headless': continue          # below threshold: not shown (Fig 10.3.5 rule)
            direc, lvl = pr; w = max(cf(a_, b_), 0.12)
            col, ls, db = ('#111', '-', False) if lvl == 'solid' else (('#111', '-', True) if lvl == 'bidir' else ('#7a86a0', (0, (4, 3)), False))
            specs.append(dict(src=direc[0], dst=direc[1], w=min(w, 1), color=col, ls=ls, double=db, cls=lvl))
            if (direc[0] not in MES and direc[1] in MES) or (direc[1] not in MES and direc[0] in MES): bridge_drawn = True
        present_mes = [s for s in MES if s not in faded]
        if not bridge_drawn and present_mes:
            target = max(present_mes, key=lambda s: len(dirs.cells[s])); up = init
            for c in cfg['bridge_order']:
                if c not in faded and c not in MES: up = c; break
            specs.append(dict(src=up, dst=target, w=0.20, color='#7a86a0', ls=(0, (4, 3)), double=False, cls='bridge'))
        _xi, _yi = POS[init]; _ty = min(_yi + node_r[init] + 0.045, 0.932)
        extra = [(_xi + dx_, _ty, 0.020, '_init', 'label') for dx_ in (-0.06, -0.03, 0.0, 0.03, 0.06)] + [(_xi, (_ty + _yi + node_r[init]) / 2, 0.012, '_init', 'label')]
        SCa = {k: v for k, v in SC.items() if k not in hidden}   # hidden states are neither drawn nor obstacles
        drawn, prob = _route(specs, POS, node_r, SCa, extra=extra); problems[arm] = prob
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal'); ax.axis('off')
        ax.add_patch(FancyBboxPatch((0.02, 0.955), 0.96, 0.042, boxstyle='round,pad=0.002,rounding_size=0.01', facecolor=arm_accent[arm], edgecolor='none', alpha=0.14, zorder=1, clip_on=False))
        ax.text(0.50, 0.976, f'{arm_display[arm]} (n = {n_arm} tumors)', ha='center', va='center', fontsize=14, fontweight='bold', color=arm_text[arm])
        r1 = node_r[init]; xi, yi = POS[init]; ty = min(yi + r1 + 0.045, 0.932)
        ax.text(xi, ty, 'Tumor initiation', ha='center', va='center', fontsize=10, fontweight='bold', color='#111', zorder=6)
        if ty - 0.014 > yi + r1 + 0.010:
            ax.add_patch(FancyArrowPatch((xi, ty - 0.014), (xi, yi + r1 + 0.006), arrowstyle='-|>', color='#111', lw=1.8, mutation_scale=12, zorder=3))
        for d in drawn:
            style = '<|-|>' if d['double'] else '-|>'
            ax.add_patch(FancyArrowPatch(d['p0'], d['p1'], arrowstyle=style, color=d['color'], linewidth=0.8 + 3.0 * d['w'], linestyle=d['ls'], mutation_scale=15, connectionstyle=f"arc3,rad={d['rad']}", zorder=3))
            if edges_log is not None:
                edges_log.append(dict(arm=arm, src=d['src'], dst=d['dst'], cls=d['cls'], width=round(d['w'], 3)))
        for s in SC:
            if s in hidden: continue
            x, y = POS[s]; r = node_r[s]; fd = s in faded; c = SC[s]
            ax.add_patch(Circle((x, y), r, facecolor='white', edgecolor=c if not fd else '#bbb', linewidth=2.6 if not fd else 1.4, zorder=5, alpha=1 if not fd else 0.35))
            tx = ax.text(x, y + r * 0.10, LB[s], ha='center', va='center', fontsize=8.4, color='#111' if not fd else '#8a8a8a', fontweight='bold', zorder=6, linespacing=0.9)
            all_labels.append((ax, tx, r, x, y))
            by = y - r - 0.016
            ax.add_patch(FancyBboxPatch((x - 0.019, by - 0.015), 0.038, 0.028, boxstyle='round,pad=0.001,rounding_size=0.004', facecolor='white', edgecolor='#666' if not fd else '#ccc', linewidth=0.8, zorder=6, alpha=1 if not fd else 0.5))
            ax.text(x, by - 0.001, f'{SI[s]}', ha='center', va='center', fontsize=9, fontweight='bold', color='#111' if not fd else '#999', zorder=7)
    fig.suptitle(title, fontsize=15, fontweight='bold', color='#0f2747', y=0.995)
    leg = [_legarrow('#111', '-', '-|>', 2.4), _legarrow('#111', '-', '<|-|>', 2.4), _legarrow('#7a86a0', (0, (4, 3)), '-|>', 2.0)]
    fig.legend(handles=leg, labels=legend_labels, handler_map={FancyArrowPatch: _HArrow()}, loc='lower center', ncol=3, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 0.025), handlelength=3.4)
    fig.text(0.5, 0.005, footer, ha='center', fontsize=10, color='#333', style='italic')
    plt.subplots_adjust(top=0.93, bottom=0.10, left=0.01, right=0.99, wspace=0.04)
    for ax in axes:
        labs = [(tx, r, x, y) for a, tx, r, x, y in all_labels if a is ax]
        pts = _fit_labels(fig, ax, labs, lw_pt=2.6)
        min_pt[ax.get_title() or len(min_pt)] = min(pts)
    out_png = Path(out_png)
    for fmt in ['png', 'pdf']:
        fig.savefig(out_png.with_suffix(f'.{fmt}'), dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return dict(min_label_pt=min(min_pt.values()), routing_problems=problems, arms={a: data['n_by_arm'][a] for a in arms})

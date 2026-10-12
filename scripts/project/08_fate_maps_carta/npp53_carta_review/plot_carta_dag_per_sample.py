#!/usr/bin/env python
"""Per-sample CARTA differentiation maps drawn in the paper's DAG convention
(Sashittal 2025 Nat Methods, Fig 1c/2d).

Layout per subplot:
  Top:    "Multipotent Progenitor" (root, {all observed states})
  Middle: intermediate progenitor boxes — each box is SPLIT into colored
          segments, one segment per cell type in that progenitor's potency set
  Bottom: individual observed cell-type leaves (single-color boxes labeled c#)
  Edges:  Hasse arrows from each parent progenitor to its direct children
          (child = smaller potency set strictly contained in parent)
          + from each terminal progenitor to the leaves for each state it contains

Grid: Intact (14, 4 rows × 4 cols), Cas (8, 2 rows × 4 cols), Cas+Enz (3, 1 row × 4 cols)
"""
from pathlib import Path
import ast
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.gridspec import GridSpec

M = '/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
RES = Path(f'{M}/analysis/evolution_features_v2/carta_fine12/results_fine12_greedy')
OUT = Path(f'{M}/analysis/evolution_features_v2/carta_fine12/biology_review')
BASE = Path(f'{M}/analysis/evolution_features_v2/carta_fine12')
VERD = pd.read_csv(OUT/'per_sample_verdict.tsv', sep='\t')

# FINAL k* consumer (chord-distance + manual overrides): read from best_k_per_sample_v2_chord.csv
BEST_K = pd.read_csv(BASE/'best_k_per_sample_v2_chord.csv')
K_FROM_CONSENSUS = BEST_K.set_index('sample')['k_star_final'].to_dict()
import json as _json
# Suppress cell states with <=20 cells (or absent) in a given sample: >20-cell states only.
KEPT = {k: set(v) for k, v in _json.load(open(BASE/'per_sample_kept_states_gt20.json')).items()}
# Fresh verdicts at final k* (recomputed against TarCA + Fitch checks at new terminal)
try:
    VERD_V2 = pd.read_csv(OUT/'per_sample_verdict_v2.tsv', sep='\t')
    VERDICT_AT_FINAL_K = VERD_V2.set_index('sample')['verdict_score'].to_dict()
except FileNotFoundError:
    VERDICT_AT_FINAL_K = {}

STATE_ORDER = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
               'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
               'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
SHORT = {s: f'c{i+1}' for i, s in enumerate(STATE_ORDER)}
STATE_COLORS = {
    'Luminal-Distal':'#1f5fa8','Luminal-Proximal':'#9467bd','Luminal-ARdep':'#d23030',
    'Luminal-CR':'#f1903d','Luminal-Persistent':'#117c70','Luminal-Squamous':'#f78ec7',
    'Basal-like':'#8c564b','Partial-EMT':'#bcbd22','Mes-Quiescent':'#c49c94',
    'Mes-Collagen':'#5b5b5b','Mes-Stemness':'#8c1a1a','Mes-myCAF':'#4daf4a',
}
ARM_COLORS = {'Intact':'#4a7fbf','Cas':'#e0955f','Cas+Enz':'#c44e52'}

K_OVERRIDES = {}   # historical manual overrides — now empty because k=2..8 consensus
                   # already picks k=4 for JZ202 and JZ206 (previously the only overrides).

# Per-sample per-leaf rad OVERRIDES for arcs that global _clearing_rad cannot
# tune well. Each leaf has ≤1 incoming arrow (from its minimal container OR
# root if dashed) so keying by (sample, leaf_c#) is unique.
#   Positive rad → arc bulges LEFT of arrow direction (matplotlib arc3).
#   Negative rad → arc bulges RIGHT.
# Leave empty by default; parallel per-panel review appends here.
ARC_RAD_OVERRIDES = {}

def _override_rad(sample, leaf_short):
    return ARC_RAD_OVERRIDES.get(sample, {}).get(leaf_short)

def read_progen(sample, k):
    """Return list of frozensets, one per progenitor node in the CARTA map."""
    f = RES/sample/f'{sample}_progenitors_{k}.txt'
    if not f.exists(): return []
    out = []
    for ln in f.read_text().splitlines():
        ln = ln.strip()
        if not ln: continue
        try:
            out.append(frozenset(ast.literal_eval(ln)))
        except Exception:
            pass
    return out

def all_edges(nodes):
    """Return full subset-containment DAG edges (parent ⊃ child, any depth)."""
    edges = []
    for i, a in enumerate(nodes):
        for j, b in enumerate(nodes):
            if i == j: continue
            if b < a: edges.append((i, j))
    return edges

def hasse_edges(nodes):
    """Transitive-reduction (Hasse) edges: parent ⊃ child with NO intermediate."""
    edges = []
    for i, a in enumerate(nodes):
        for j, b in enumerate(nodes):
            if i == j: continue
            if b < a:
                intermediate = any(k != i and k != j and b < nodes[k] < a
                                   for k in range(len(nodes)))
                if not intermediate:
                    edges.append((i, j))
    return edges

def layout_dag(nodes, edges):
    """Assign (x, y) to each node. y = topological depth from Root, x = spread."""
    n = len(nodes)
    depth = [0]*n
    changed = True
    while changed:
        changed = False
        for p, c in edges:
            if depth[c] < depth[p] + 1:
                depth[c] = depth[p] + 1
                changed = True
    max_d = max(depth) if depth else 0
    # group by depth, spread x
    by_d = {}
    for i, d in enumerate(depth): by_d.setdefault(d, []).append(i)
    xy = {}
    for d, idxs in by_d.items():
        # size-ordered from left (larger set on left) for readability
        idxs_sorted = sorted(idxs, key=lambda k: (-len(nodes[k]), min(STATE_ORDER.index(s) for s in nodes[k])))
        m = len(idxs_sorted)
        for k, i in enumerate(idxs_sorted):
            x = (k + 1) / (m + 1)   # 0..1
            xy[i] = (x, d)
    return xy, max_d

def _dark(hex_color):
    r,g,b = [int(hex_color.lstrip('#')[i:i+2],16)/255 for i in (0,2,4)]
    return (0.299*r + 0.587*g + 0.114*b) < 0.55

def draw_split_box(ax, cx, cy, states, box_w=0.14, box_h=0.06, edge='#333', lw=0.6):
    """Rectangle split into colored segments (one per state in potency),
    with c# label INSIDE each segment (white on dark, black on light)."""
    n = len(states)
    if n == 0: return
    seg_w = box_w / n
    x0 = cx - box_w/2
    y0 = cy - box_h/2
    states_ord = sorted(states, key=lambda s: STATE_ORDER.index(s))
    # Scale label font to segment width (segment width in figure units → point size)
    seg_pts = seg_w * 1000  # rough scale — 1 axes-unit ≈ 900 pt at default fig size
    label_fs = min(9.0, max(4.8, seg_pts * 1.05))
    for k, s in enumerate(states_ord):
        ax.add_patch(Rectangle((x0 + k*seg_w, y0), seg_w, box_h,
                                facecolor=STATE_COLORS[s], edgecolor=edge, lw=lw))
        text_c = '#fff' if _dark(STATE_COLORS[s]) else '#111'
        ax.text(x0 + k*seg_w + seg_w/2, cy, SHORT[s],
                ha='center', va='center', fontsize=label_fs, color=text_c,
                fontweight='bold')


def _seg_cx(box_cx, box_w, states, s):
    """x of the segment center for state s within a split box (arrows leave
    from under their own segment, not the box center — avoids sweeping across
    the box bottom)."""
    ordered = sorted(states, key=lambda x: STATE_ORDER.index(x))
    seg_w = box_w / max(1, len(ordered))
    return (box_cx - box_w/2) + (ordered.index(s) + 0.5) * seg_w

def draw_leaf(ax, cx, cy, state, box_w=0.052, box_h=0.062):
    ax.add_patch(Rectangle((cx-box_w/2, cy-box_h/2), box_w, box_h,
                            facecolor=STATE_COLORS[state], edgecolor='#333', lw=0.5))
    ax.text(cx, cy-box_h*0.9, SHORT[state], ha='center', va='top', fontsize=8.5, color='#111')

def draw_arrow(ax, x0, y0, x1, y1, color='#333', lw=0.7, curved=False,
               rad_sign=None, explicit_rad=None, dashed=False):
    """Straight arrow if curved=False; arc-3 curved arrow if True.
    rad_sign: force curve direction (+1 = right, -1 = left); default = sign(dx).
    explicit_rad: bypass formula and use this rad value directly.
    dashed: dashed line style (used for root→leaf 'direct-from-progenitor' arcs)."""
    linestyle = (0, (2.5, 2.0)) if dashed else '-'
    if curved:
        if explicit_rad is not None:
            rad = explicit_rad
        else:
            dx = x1 - x0
            s = rad_sign if rad_sign is not None else (1 if dx >= 0 else -1)
            rad = 0.38 * s * min(1.1, abs(dx) * 3.5 + 0.18)
        arrow = FancyArrowPatch((x0, y0), (x1, y1),
                                connectionstyle=f'arc3,rad={rad}',
                                arrowstyle='-|>', color=color, linewidth=lw,
                                linestyle=linestyle,
                                mutation_scale=6, shrinkA=1, shrinkB=1)
        ax.add_patch(arrow)
    else:
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle='-|>', color=color, lw=lw,
                                    linestyle=linestyle,
                                    shrinkA=2, shrinkB=2, mutation_scale=7))

def _outward_rad(lx, panel_center=0.5):
    """Base outward rad for a DASHED root→leaf arc (used as the starting
    magnitude before collision-checking widens it further if needed)."""
    away = abs(lx - panel_center)
    mag = 0.95 - 1.5 * away
    mag = max(0.40, min(0.92, mag))
    sign = +1 if lx <= panel_center else -1
    return sign * mag

def _dashed_arc_rad(x0, y0, x1, y1, boxes, exclude_indices=(),
                    max_rad=1.20, step=0.04, tol=0.028):
    """Rad for a DASHED root→leaf arc: start with `_outward_rad`'s magnitude
    (sign fixed by leaf side) and iteratively grow |rad| until the arc clears
    every intermediate box in its way. Guarantees the arc always curves
    outward AND clears every box, even the wide top-level intermediates."""
    base = _outward_rad(x1)
    side = +1 if base >= 0 else -1
    dx, dy = x1 - x0, y1 - y0
    r_mag = max(abs(base), 0.30)
    excl = set(exclude_indices)
    while r_mag <= max_rad + 1e-9:
        rad = side * r_mag
        cx = (x0 + x1)/2 + rad*dy
        cy = (y0 + y1)/2 - rad*dx
        all_clear = True
        y_top = max(y0, y1); y_bot = min(y0, y1)
        for xb, yb, hw, hh, idx in boxes:
            if idx in excl: continue
            if yb + hh <= y_bot or yb - hh >= y_top: continue
            for f in (-0.95, -0.5, 0.0, 0.5, 0.95):
                y_probe = yb + f * hh
                for arc_x in _arc_x_at(y_probe, x0, y0, x1, y1, cx, cy):
                    if xb - hw - tol < arc_x < xb + hw + tol:
                        all_clear = False; break
                if not all_clear: break
            if not all_clear: break
        if all_clear: return rad
        r_mag += step
    return side * max_rad

def _box_extents(nodes, node_xy_pts, root_idx, box_h=0.078):
    """Return list of (xb, yb, half_w, half_h, node_index) for each non-root
    progenitor box. Widths match draw_split_box's computed width."""
    out = []
    for i, (cx, cy) in node_xy_pts.items():
        if i == root_idx: continue
        n_st = len(nodes[i])
        box_w = min(0.62, max(0.16, 0.052 * n_st)) if n_st >= 1 else 0.02
        out.append((cx, cy, box_w/2, box_h/2, i))
    return out

def _straight_blockers(x0, y0, x1, y1, boxes, exclude_indices=(), tol=0.003):
    """List of (xb, yb, hw, hh) for every box the straight segment passes
    through (excluding the endpoint boxes)."""
    if y1 == y0: return []
    excl = set(exclude_indices)
    out = []
    y_top = max(y0, y1); y_bot = min(y0, y1)
    for xb, yb, hw, hh, idx in boxes:
        if idx in excl: continue
        if yb + hh - tol <= y_bot or yb - hh + tol >= y_top: continue
        # Sample line at 3 y's within box vertical span
        for y_probe in (yb - hh*0.9, yb, yb + hh*0.9):
            if not (y_bot < y_probe < y_top): continue
            t = (y_probe - y0) / (y1 - y0)
            if not (0 < t < 1): continue
            xline = x0 + t * (x1 - x0)
            if xb - hw - tol <= xline <= xb + hw + tol:
                out.append((xb, yb, hw, hh))
                break
    return out

def _arc_x_at(y_probe, x0, y0, x1, y1, cx, cy):
    """Return list of arc x-values where the bezier reaches y=y_probe (0-2 hits)."""
    a = y0 - 2*cy + y1
    b = 2*(cy - y0)
    c = y0 - y_probe
    ts = []
    if abs(a) < 1e-12:
        if abs(b) > 1e-12: ts.append(-c / b)
    else:
        disc = b*b - 4*a*c
        if disc >= 0:
            sd = np.sqrt(disc)
            ts.extend([(-b - sd)/(2*a), (-b + sd)/(2*a)])
    xs = []
    for t in ts:
        if 0 < t < 1:
            xs.append((1-t)**2*x0 + 2*(1-t)*t*cx + t**2*x1)
    return xs

def _clearing_rad(x0, y0, x1, y1, boxes, exclude_indices=(),
                  max_rad=1.15, step=0.04, tol=0.028):
    """Find the smallest |rad| with correct sign such that the matplotlib
    arc3 curve from (x0,y0)->(x1,y1) clears every intermediate box its
    straight-line path would hit. Returns None if no collision.
    Iteratively grows rad until every blocker's x-range is missed at all
    sampled y within the box's vertical span, with a visible tol margin
    (~2% axes ≈ 4 pixels at 170 dpi) so cleared arcs read as *not touching*
    the box, not just mathematically outside."""
    intersecting = _straight_blockers(x0, y0, x1, y1, boxes, exclude_indices, tol)
    if not intersecting: return None
    # Side sign: use first blocker's straight-line hit position
    xb0, yb0, hw0, _ = intersecting[0]
    t0 = (yb0 - y0) / (y1 - y0)
    xline0 = x0 + t0 * (x1 - x0)
    side = +1 if xline0 < xb0 else -1
    dx, dy = x1 - x0, y1 - y0
    # Start at r_mag=0.30 (never draw a barely-visible arc) and step 0.04.
    r_mag = 0.30
    while r_mag <= max_rad + 1e-9:
        rad = side * r_mag
        cx = (x0 + x1)/2 + rad*dy
        cy = (y0 + y1)/2 - rad*dx
        all_clear = True
        for xb, yb, hw, hh in intersecting:
            # Sample 5 y within box span for a stricter clearance check
            for f in (-0.95, -0.5, 0.0, 0.5, 0.95):
                y_probe = yb + f * hh
                for arc_x in _arc_x_at(y_probe, x0, y0, x1, y1, cx, cy):
                    if xb - hw - tol < arc_x < xb + hw + tol:
                        all_clear = False; break
                if not all_clear: break
            if not all_clear: break
        if all_clear: return rad
        r_mag += step
    return side * max_rad

def _crosses_any_box(x0, y0, x1, y1, boxes, exclude_indices=(), tol=0.003):
    """Back-compat wrapper: returns side_sign for callers that only need
    'does the straight line collide?' — see _clearing_rad for the geometry
    that actually clears the blocker."""
    blockers = _straight_blockers(x0, y0, x1, y1, boxes, exclude_indices, tol)
    if not blockers: return None
    xb, yb, _, _ = blockers[0]
    t = (yb - y0) / (y1 - y0)
    xline = x0 + t * (x1 - x0)
    return +1 if xline < xb else -1

def draw_carta_map(ax, sample, arm, k, chain, verdict, n_leaves, n_states):
    # Suppress states with <=20 cells (or absent) in THIS sample: intersect every
    # progenitor potency set with the sample's >20-cell states, drop empties.
    kept = KEPT.get(sample)
    if kept is not None:
        chain = [frozenset(pset) & kept for pset in chain]
        chain = [pset for pset in chain if len(pset) > 0]
    # Deduplicate progenitor nodes and layout
    nodes = []
    for s in chain:
        if s not in nodes: nodes.append(s)
    if not nodes:
        ax.text(0.5, 0.5, '(no CARTA output)', transform=ax.transAxes, ha='center')
        ax.axis('off'); return
    h_edges = hasse_edges(nodes)                   # BLACK: transitive-reduction (direct parent-child)
    a_edges = all_edges(nodes)                     # ALL subset-containment edges
    r_edges = [e for e in a_edges if e not in h_edges]   # RED: shortcut / convergent-diff edges
    xy, max_d = layout_dag(nodes, h_edges)

    # Discover observed cell-type leaves = union of all sets present at any level
    # (since Root is {all present states}, we take those)
    root_idx = max(range(len(nodes)), key=lambda i: len(nodes[i]))
    all_states = sorted(nodes[root_idx], key=lambda s: STATE_ORDER.index(s))
    n_leaf_states = len(all_states)

    # Y layout: pack progenitor rows near the TOP so the arrows to leaves are SHORT
    # (matches CARTA paper convention — final leaf-connection lines stay low).
    # Progenitor rows compress into the upper ~60% of the axes; leaves sit at y=0.08
    # with a small gap above.
    n_prog_rows = max_d + 1
    top_y  = 0.94
    prog_bottom_y = 0.35   # last progenitor row anchor — leaves 0.08 → arrows short
    leaf_y = 0.08
    if n_prog_rows == 1:
        y_of_depth = {0: (top_y + prog_bottom_y) / 2}
    else:
        y_of_depth = {d: top_y - (top_y - prog_bottom_y) * d / (n_prog_rows - 1)
                      for d in range(n_prog_rows)}

    # Panel title — plain black, no verdict tag (per request — subjective verdicts removed)
    ax.text(0.02, 1.02, f'{sample} · {arm}',
            transform=ax.transAxes, ha='left', va='bottom',
            fontsize=12, fontweight='bold', color='#111')
    ax.text(0.98, 1.02, f'k*={k} · {int(n_leaves)} cells · {len(all_states)} states (>20 cells)',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=9.5, color='#111')

    # Draw progenitor nodes
    node_xy_pts = {}
    for i, (x, d) in xy.items():
        y = y_of_depth[d]
        cx = 0.08 + 0.84 * x
        cy = y
        node_xy_pts[i] = (cx, cy)
        if d == 0:
            # Root: plain "Progenitor" text label at top — no split-color box (per user spec)
            ax.text(cx, cy, 'Progenitor',
                    ha='center', va='center', fontsize=11.5, fontweight='bold', color='#111')
        else:
            # Intermediate progenitor: split-color box with c# labels inside segments.
            # Wider per-state so c# labels have breathing room (was 0.024, now 0.045; caps up to 0.55).
            w = min(0.62, max(0.16, 0.052 * len(nodes[i])))
            draw_split_box(ax, cx, cy, nodes[i], box_w=w, box_h=0.078)

    # Compute intermediate-progenitor box extents for collision-avoidance curving.
    boxes = _box_extents(nodes, node_xy_pts, root_idx)

    # BLACK arrows between progenitors: Hasse (direct parent-child) only.
    # Non-Hasse subset-containment edges are TRANSITIVELY IMPLIED by the Hasse
    # chain (root ⊃ L1 ⊃ … ⊃ terminal), so drawing them adds visual noise but
    # zero information — drop them entirely (matches CARTA-paper convention).
    for p, c in h_edges:
        (x0, y0) = node_xy_pts[p]; (x1, y1) = node_xy_pts[c]
        y0a, y1a = y0 - 0.038, y1 + 0.033
        rad_val = _clearing_rad(x0, y0a, x1, y1a, boxes, exclude_indices=(p, c))
        if rad_val is None:
            draw_arrow(ax, x0, y0a, x1, y1a, color='#111', lw=0.7, curved=False)
        else:
            draw_arrow(ax, x0, y0a, x1, y1a, color='#111', lw=0.7,
                       curved=True, explicit_rad=rad_val)

    # Draw observed cell-type leaves at the bottom row
    for k_leaf, s in enumerate(all_states):
        cx = 0.08 + 0.84 * (k_leaf + 0.5) / max(1, n_leaf_states)
        cy = leaf_y
        draw_leaf(ax, cx, cy, s)

    # Progenitor → observed cell-type arrows — two rules:
    #   (a) SOLID arrow from MINIMAL intermediate container(s) for states that
    #       appear in at least one intermediate progenitor.
    #   (b) DASHED outward arc from ROOT for states that appear ONLY in root
    #       (no intermediate contains them → biologically "direct-from-
    #       progenitor" differentiation). The dashed arc curves away from
    #       panel center, bypassing all intermediate boxes.
    for s in all_states:
        inter_containers = [i for i in range(len(nodes))
                            if s in nodes[i] and i != root_idx]
        k_leaf = all_states.index(s)
        lx = 0.08 + 0.84 * (k_leaf + 0.5) / max(1, n_leaf_states)

        override = _override_rad(sample, SHORT[s])
        if inter_containers:
            # (a) SOLID arrows from minimal-among-intermediate containers only
            minimals = [i for i in inter_containers
                        if not any(j != i and nodes[j] < nodes[i]
                                   and s in nodes[j]
                                   for j in inter_containers)]
            for i in minimals:
                (px, py) = node_xy_pts[i]
                w_i = min(0.62, max(0.16, 0.052 * len(nodes[i])))
                sx = _seg_cx(px, w_i, nodes[i], s)     # leave from under the segment
                y0a = py - 0.041; y1a = leaf_y + 0.033
                if override is not None:
                    rad_val = override
                else:
                    rad_val = _clearing_rad(sx, y0a, lx, y1a, boxes,
                                            exclude_indices=(i,))
                if rad_val is None:
                    draw_arrow(ax, sx, y0a, lx, y1a, color='#111', lw=0.55,
                               curved=False)
                else:
                    draw_arrow(ax, sx, y0a, lx, y1a, color='#111', lw=0.55,
                               curved=True, explicit_rad=rad_val)
        else:
            # (b) DASHED outward arc from Progenitor — direct differentiation.
            (px, py) = node_xy_pts[root_idx]
            if override is not None:
                rad_val = override
            else:
                rad_val = _dashed_arc_rad(px, py - 0.020, lx, leaf_y + 0.033,
                                           boxes, exclude_indices=())
            draw_arrow(ax, px, py - 0.020, lx, leaf_y + 0.033,
                       color='#555', lw=0.5,
                       curved=True, explicit_rad=rad_val, dashed=True)

    # arm-tinted background
    ax.set_facecolor((*[int(ARM_COLORS[arm][i:i+2],16)/255 for i in (1,3,5)], 0.05))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.05)
    ax.axis('off')

# --- Build grid figure ---
intact = sorted(VERD[VERD['arm']=='Intact']['sample'].tolist())
cas    = sorted(VERD[VERD['arm']=='Cas']['sample'].tolist())
casenz = sorted(VERD[VERD['arm']=='Cas+Enz']['sample'].tolist())

NCOL = 4
n_int_r = int(np.ceil(len(intact)/NCOL))
n_cas_r = int(np.ceil(len(cas)/NCOL))
n_ce_r  = int(np.ceil(len(casenz)/NCOL))
total_r = n_int_r + n_cas_r + n_ce_r

fig = plt.figure(figsize=(20, total_r*3.6 + 1.4))
gs = GridSpec(total_r, NCOL, figure=fig, hspace=0.50, wspace=0.15)

def fill(samples, row_off):
    for i, samp in enumerate(samples):
        r = i // NCOL + row_off; c = i % NCOL
        ax = fig.add_subplot(gs[r, c])
        rr = VERD[VERD['sample']==samp].iloc[0]
        # Precedence: K_OVERRIDES > new k=2..8 consensus > legacy verdict k_star
        if samp in K_OVERRIDES:
            k = K_OVERRIDES[samp]
        elif samp in K_FROM_CONSENSUS and pd.notna(K_FROM_CONSENSUS[samp]):
            k = int(K_FROM_CONSENSUS[samp])
        else:
            k = int(rr['k_star'])
        chain = read_progen(samp, k)
        # verdict at NEW k* (fresh scoring) — falls back to legacy per_sample_verdict.tsv if missing
        verdict_at_new_k = int(VERDICT_AT_FINAL_K.get(samp, rr['verdict_score']))
        draw_carta_map(ax, samp, rr['arm'], k, chain,
                       verdict_at_new_k, int(rr['n_leaves_tree']), int(rr['n_states_kept']))
fill(intact, 0)
fill(cas,    n_int_r)
fill(casenz, n_int_r + n_cas_r)

# Legend — 12 states color→c#
from matplotlib.patches import Patch
handles = [Patch(facecolor=STATE_COLORS[s], edgecolor='#333', label=f'{SHORT[s]} = {s}') for s in STATE_ORDER]
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.998),
           ncol=6, fontsize=11, frameon=True, columnspacing=1.2, handlelength=1.4)
fig.suptitle('Per-sample CARTA differentiation maps at k*  ·  DAG convention (Sashittal 2025 Nat Methods style)',
             y=1.005, fontsize=15, fontweight='bold', color='#0f2747')
fig.subplots_adjust(top=0.94, bottom=0.02, left=0.02, right=0.99)

for fmt in ['pdf','png']:
    fig.savefig(OUT/f'carta_dag_per_sample.{fmt}', dpi=170, bbox_inches='tight')
plt.close(fig)
print(f'saved: {OUT}/carta_dag_per_sample.png')

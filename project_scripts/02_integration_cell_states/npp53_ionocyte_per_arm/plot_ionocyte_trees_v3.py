#!/usr/bin/env python
# =============================================================================
# Radial Greedy-tree plots — v3
# Layout: CIRCLE CLADOGRAM (same as Npp53_tree_viewer_interactive.html default).
#   - leaves at outer radius, evenly spaced angles
#   - internal-node radius = (H_max - hgt(node)) / H_max * rMax  (root at center)
#   - hgt(leaf)=0, hgt(parent-of-leaves)=1, ... hgt(root)=tree height
# Two annotation rings:
#   - inner: detailed subtype (putative_state, 14 states)
#   - outer: major cell type (6 groups: Basal / Luminal / Mesenchymal / EMT / Cycling / Ionocyte)
# Ionocyte highlighting preserved: big red leaves + purple Ion-ancestor nodes + bold red Ion-lineage edges.
# =============================================================================
import numpy as np, pandas as pd
from pathlib import Path
from ete3 import Tree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Wedge
from matplotlib.backends.backend_pdf import PdfPages

M   = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing')
OUT = M/'1-mouse_based/Npp53_tracing/analysis/evolution_features_v2'
GREEDY_DIR = M/'cassiopeia_trees_output'
SAMPLES = [('JZ204', 'Cas',    78),
           ('JZ205', 'Cas',    62),
           ('JZ201', 'Cas+Enz', 39),
           ('JZ203', 'Cas+Enz', 36),
           ('JZ214', 'Intact', 13),
           ('JZ202', 'Cas+Enz', 12)]

STATE_COLOR = {
    'Basal-like':          '#4a7fbf',
    'Cycling':             '#9e9e9e',
    'Ionocyte':            '#c44e52',
    'Luminal-ARdep':        '#7cb342',
    'Luminal-CR':          '#ffb300',
    'Luminal-Distal':      '#26a69a',
    'Luminal-Persistent':  '#5c6bc0',
    'Luminal-Proximal':    '#ec407a',
    'Luminal-Squamous':    '#8d6e63',
    'Mes-Collagen':        '#78909c',
    'Mes-Quiescent':       '#bcaaa4',
    'Mes-Stemness':        '#ce93d8',
    'Mes-myCAF':           '#a1887f',
    'Partial-EMT':         '#ffca28',
    'None':                '#e0e0e0',
}
MAJOR_COLOR = {
    'Basal':       '#4a7fbf',
    'Luminal':     '#26a69a',
    'Mesenchymal': '#78909c',
    'EMT':         '#ffca28',
    'Cycling':     '#9e9e9e',
    'Ionocyte':    '#c44e52',
    'Other':       '#e0e0e0',
}
def state_to_major(s):
    if s is None or (isinstance(s,float) and np.isnan(s)): return 'Other'
    if s == 'Basal-like':   return 'Basal'
    if s == 'Ionocyte':     return 'Ionocyte'
    if s == 'Cycling':      return 'Cycling'
    if s == 'Partial-EMT':  return 'EMT'
    if s.startswith('Luminal-'): return 'Luminal'
    if s.startswith('Mes-'):     return 'Mesenchymal'
    return 'Other'

# ---- CLADOGRAM LAYOUT (ported from Npp53_tree_viewer_interactive.html: layoutCladogram) ----
def layout_cladogram(tree, R_leaf=1.0):
    """Circle cladogram: leaves at r=R_leaf on evenly-spaced angles.
       Internal-node r = (H_max - hgt(n)) / H_max * R_leaf.  Root at r=0."""
    leaves = tree.get_leaves()
    n_leaves = len(leaves)
    angle = {}
    for i, l in enumerate(leaves):
        angle[l] = 2*np.pi * i / n_leaves
    # internal-node angle = mean of children's angles (postorder)
    for n in tree.traverse('postorder'):
        if n.is_leaf(): continue
        angle[n] = np.mean([angle[c] for c in n.children])
    # graph height from leaves
    hgt = {}
    def compute_hgt(n):
        if n.is_leaf():
            hgt[n] = 0
            return 0
        ch_hgt = [compute_hgt(c) for c in n.children]
        hgt[n] = 1 + max(ch_hgt)
        return hgt[n]
    H_max = compute_hgt(tree)
    if H_max == 0: H_max = 1
    pos = {}
    for n in tree.traverse():
        r = ((H_max - hgt[n]) / H_max) * R_leaf
        t = angle[n]
        pos[n] = (r, t)
    return pos, H_max

def polar_to_cart(r, theta):
    return r*np.cos(theta), r*np.sin(theta)

def draw_tree(ax, tree, pos, state_map, sample_name, arm, n_ion,
              R_leaf=1.0, ring1_gap=0.02, ring_width=0.055,
              ion_leaf_size=45, ion_node_size=18, leaf_size=3.0,
              edge_lw_base=0.30, edge_lw_ion=1.5, arc_n=24):
    # Ionocyte-descendant nodes
    ion_desc_nodes = set()
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            if state_map.get(node.name) == 'Ionocyte':
                ion_desc_nodes.add(node)
        else:
            if any(c in ion_desc_nodes for c in node.children):
                ion_desc_nodes.add(node)
    # -----------------------------------------------------------------------
    # CLADOGRAM branch drawing (same as Npp53_tree_viewer_interactive.html):
    # for each parent→child edge, draw TWO parts:
    #  1) ARC at parent's radius sweeping from θ_parent to θ_child (approximated by segments)
    #  2) STRAIGHT RADIAL spoke from (r_parent, θ_child) to (r_child, θ_child)
    # -----------------------------------------------------------------------
    edge_segments, edge_widths, edge_colors = [], [], []
    # arcs at parent's radius: for each internal node with children, draw one arc from min-child-angle to max-child-angle at r_parent
    for parent in tree.traverse():
        if parent.is_leaf() or not parent.children: continue
        r_p = pos[parent][0]
        child_angles = sorted([pos[c][1] for c in parent.children])
        if len(child_angles) < 2: continue
        t_min, t_max = child_angles[0], child_angles[-1]
        # detect if arc goes across 0/2π boundary — for cladogram with mean-angle interior nodes this normally doesn't happen
        # arc as segments
        thetas = np.linspace(t_min, t_max, arc_n)
        pts = [(r_p*np.cos(t), r_p*np.sin(t)) for t in thetas]
        # color: red if the parent is an Ion-ancestor AND every child in this arc range is also Ion-desc (whole parent is Ion-lineage)
        parent_is_ion = parent in ion_desc_nodes
        all_children_ion = all(c in ion_desc_nodes for c in parent.children) if parent_is_ion else False
        arc_color = '#c44e52' if all_children_ion else '#888'
        arc_lw    = edge_lw_ion if all_children_ion else edge_lw_base
        for i in range(len(pts)-1):
            edge_segments.append([pts[i], pts[i+1]])
            edge_widths.append(arc_lw)
            edge_colors.append(arc_color)
    # radial spokes: for each child, draw straight line from (r_parent, θ_child) to (r_child, θ_child)
    for node in tree.traverse('preorder'):
        if node.up is None: continue
        r_p = pos[node.up][0]; t_c = pos[node][1]; r_c = pos[node][0]
        x1, y1 = polar_to_cart(r_p, t_c)
        x2, y2 = polar_to_cart(r_c, t_c)
        edge_segments.append([(x1,y1), (x2,y2)])
        if node in ion_desc_nodes and node.up in ion_desc_nodes:
            edge_widths.append(edge_lw_ion); edge_colors.append('#c44e52')
        else:
            edge_widths.append(edge_lw_base); edge_colors.append('#888')
    ax.add_collection(LineCollection(edge_segments, colors=edge_colors, linewidths=edge_widths, alpha=0.85))
    # non-Ion leaves as tiny dots colored by state
    xs, ys, cs = [], [], []
    ion_leaf_x, ion_leaf_y = [], []
    for l in tree.get_leaves():
        r, t = pos[l]
        x, y = polar_to_cart(r, t)
        st = state_map.get(l.name, 'None')
        if st == 'Ionocyte':
            ion_leaf_x.append(x); ion_leaf_y.append(y)
        else:
            xs.append(x); ys.append(y); cs.append(STATE_COLOR.get(st, '#e0e0e0'))
    ax.scatter(xs, ys, c=cs, s=leaf_size, edgecolors='none', zorder=3)
    ax.scatter(ion_leaf_x, ion_leaf_y, c='#c44e52', s=ion_leaf_size,
               edgecolors='k', linewidths=0.6, zorder=5)
    # Ion-ancestor internal nodes
    ion_int_x, ion_int_y = [], []
    for n in ion_desc_nodes:
        if n.is_leaf(): continue
        r, t = pos[n]
        x, y = polar_to_cart(r, t)
        ion_int_x.append(x); ion_int_y.append(y)
    ax.scatter(ion_int_x, ion_int_y, s=ion_node_size, c='#8e24aa', marker='o',
               edgecolors='k', linewidths=0.4, zorder=4)
    # ---- TWO ANNOTATION RINGS ----
    # Ring 1 (inner) — detailed subtype (14 states)
    # Ring 2 (outer) — major cell type (6 groups)
    R1_in  = R_leaf + ring1_gap
    R1_out = R1_in + ring_width
    R2_in  = R1_out + 0.008
    R2_out = R2_in + ring_width
    n_leaves = len(tree.get_leaves())
    half_w   = np.degrees(np.pi / n_leaves)
    wedges1, wedges2 = [], []
    for l in tree.get_leaves():
        r, t = pos[l]
        st = state_map.get(l.name, 'None')
        col1 = STATE_COLOR.get(st, '#e0e0e0')
        col2 = MAJOR_COLOR.get(state_to_major(st), '#e0e0e0')
        theta_deg = np.degrees(t)
        wedges1.append(Wedge((0,0), R1_out, theta_deg - half_w, theta_deg + half_w,
                              width=ring_width, facecolor=col1, edgecolor='none'))
        wedges2.append(Wedge((0,0), R2_out, theta_deg - half_w, theta_deg + half_w,
                              width=ring_width, facecolor=col2, edgecolor='none'))
    ax.add_collection(PatchCollection(wedges1, match_original=True))
    ax.add_collection(PatchCollection(wedges2, match_original=True))
    # ring labels — small text at top (12 o'clock) — "subtype" and "major"
    ax.text(0, (R1_in + R1_out)/2, 'inner: subtype', fontsize=7, ha='center', va='center', color='#666',
            rotation=0, bbox=dict(facecolor='white', edgecolor='none', pad=1, alpha=0.85))
    ax.text(0, (R2_in + R2_out)/2, 'outer: major', fontsize=7, ha='center', va='center', color='#666',
            rotation=0, bbox=dict(facecolor='white', edgecolor='none', pad=1, alpha=0.85))
    # cosmetics
    ax.set_aspect('equal')
    max_r = R2_out * 1.06
    ax.set_xlim(-max_r, max_r); ax.set_ylim(-max_r, max_r)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    n_leaves_matched = len([l for l in tree.get_leaves() if l.name in state_map])
    n_ion_leaves = sum(1 for l in tree.get_leaves() if state_map.get(l.name) == 'Ionocyte')
    n_ion_anc = len([n for n in ion_desc_nodes if not n.is_leaf()])
    ax.set_title(f'{sample_name} ({arm}) — Greedy tree · cladogram layout\n'
                 f'n_leaves = {len(tree.get_leaves())} (matched = {n_leaves_matched}) · '
                 f'Ionocyte leaves = {n_ion_leaves} · Ion-ancestor internal nodes = {n_ion_anc}',
                 fontsize=13)
    # legend (upper-left)
    states_present = sorted({state_map.get(l.name) for l in tree.get_leaves() if state_map.get(l.name)})
    handles = []
    # subtype legend
    handles.append(plt.Line2D([0],[0], marker='s', color='w', markerfacecolor='none', markeredgecolor='none',
                              label='─ Subtype (inner ring)', markersize=0))
    for st in states_present:
        col = STATE_COLOR.get(st, '#e0e0e0')
        handles.append(plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=col,
                                    markersize=11 if st == 'Ionocyte' else 7, label=st,
                                    markeredgecolor='k' if st == 'Ionocyte' else 'none', markeredgewidth=0.5))
    handles.append(plt.Line2D([0],[0], marker='s', color='w', markerfacecolor='none', markeredgecolor='none',
                              label='─ Major type (outer ring)', markersize=0))
    majors_present = sorted({state_to_major(st) for st in states_present})
    for mj in majors_present:
        col = MAJOR_COLOR.get(mj, '#e0e0e0')
        handles.append(plt.Line2D([0],[0], marker='s', color='w', markerfacecolor=col,
                                    markersize=8, label=mj, markeredgecolor='none'))
    handles.append(plt.Line2D([0],[0], marker='o', color='w', markerfacecolor='#8e24aa',
                              markeredgecolor='k', markeredgewidth=0.5, markersize=8, label='Ion-ancestor node'))
    ax.legend(handles=handles, loc='upper left', fontsize=8, ncol=1, frameon=False, bbox_to_anchor=(-0.08, 1.02))

def main():
    sc_fp = pd.read_csv(M/'1-mouse_based/Npp53_tracing/analysis/evolution_fitness_plasticity/sc_fitness_plasticity.csv')
    pdf_path = OUT/'Fig_7_2_5_greedy_cladogram_per_sample.pdf'
    with PdfPages(pdf_path) as pdf:
        for s, arm, n_ion in SAMPLES:
            treef = GREEDY_DIR/f'{s}_tree.nwk'
            if not treef.exists():
                print(f'  {s}: no tree, skip'); continue
            tr = Tree(str(treef), format=1)
            sub = sc_fp[sc_fp['sample']==s].drop_duplicates('cellBC').set_index('cellBC')
            state_map = sub['putative_state'].to_dict()
            pos, H_max = layout_cladogram(tr, R_leaf=1.0)
            fig, ax = plt.subplots(1, 1, figsize=(13, 13))
            draw_tree(ax, tr, pos, state_map, s, arm, n_ion)
            fig.tight_layout()
            png_path = OUT/f'Fig_7_2_5_greedy_clado_{s}.png'
            fig.savefig(png_path, dpi=180, bbox_inches='tight')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            print(f'saved: {png_path.name}  (n_leaves={len(tr.get_leaves())}  H_max={H_max})')
    print(f'\nsaved combined PDF: {pdf_path.name}')

if __name__ == '__main__':
    main()

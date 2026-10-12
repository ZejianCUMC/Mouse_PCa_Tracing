#!/usr/bin/env python3
"""
Update the tree viewer's CARTA panel with new v2 DAG data.
1. Replace "carta":{...} in the DATA JSON with new DAG-structured data
2. Replace the renderCarta() function with SVG DAG renderer
"""
import ast, csv, json, re
from pathlib import Path

M = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing')
RES = M / 'analysis/evolution_features_v2/carta_fine12/results_fine12_greedy'
BASE = M / 'analysis/evolution_features_v2/carta_fine12'
INPUTS = BASE / 'inputs_fine12'
VIEWER = Path('/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/Npp53_tree_viewer_interactive.html')

STATE_ORDER = ['Luminal-Distal','Luminal-Proximal','Luminal-ARdep','Luminal-CR',
               'Luminal-Persistent','Luminal-Squamous','Basal-like','Partial-EMT',
               'Mes-Quiescent','Mes-Collagen','Mes-Stemness','Mes-myCAF']
SHORT = {s: f'c{i+1}' for i, s in enumerate(STATE_ORDER)}
VIEWER_NAME = {
    'Luminal-Distal':'c1_Lum-Distal','Luminal-Proximal':'c2_Lum-Proximal',
    'Luminal-ARdep':'c3_Lum-ARdep','Luminal-CR':'c4_Lum-CR',
    'Luminal-Persistent':'c5_Lum-Persistent','Luminal-Squamous':'c6_Lum-Squamous',
    'Basal-like':'c7_Basal-like','Partial-EMT':'c8_Partial-EMT',
    'Mes-Quiescent':'c9_Mes-Quiescent','Mes-Collagen':'c10_Mes-Collagen',
    'Mes-Stemness':'c11_Mes-Stemness','Mes-myCAF':'c12_Mes-myCAF',
}

# ── Read k* data ─────────────────────────────────────────────────────────────
kstar = {}
with open(BASE / 'best_k_per_sample_v2_chord.csv') as f:
    for row in csv.DictReader(f):
        kstar[row['sample']] = {
            'k': int(float(row['k_star_final'])),
            'arm': row['arm'],
            'n_states': int(row['n_states']),
            'n_leaves': int(row['n_leaves']),
        }

# ── Read observed states per sample from labels files ────────────────────────
def get_observed_states(sample):
    f = INPUTS / f'{sample}_labels_fine12.tsv'
    if not f.exists():
        return set()
    states = set()
    for ln in f.read_text().splitlines()[1:]:  # skip header
        parts = ln.strip().split('\t')
        if len(parts) >= 2 and parts[1] in STATE_ORDER:
            states.add(parts[1])
    return states

# ── Read progenitor sets (skip universal root = first line) ──────────────────
def read_progen(sample, k):
    f = RES / sample / f'{sample}_progenitors_{k}.txt'
    if not f.exists():
        return []
    out = []
    for ln in f.read_text().splitlines():
        ln = ln.strip()
        if not ln:
            continue
        try:
            out.append(frozenset(ast.literal_eval(ln)))
        except:
            pass
    # Skip line 1 (universal root = all 12 states)
    if len(out) > 1:
        return out[1:]  # intermediate progenitors only
    return out

# ── Hasse edges (transitive reduction of containment) ────────────────────────
def hasse_edges(nodes):
    edges = []
    for i, a in enumerate(nodes):
        for j, b in enumerate(nodes):
            if i == j:
                continue
            if b < a:  # b is proper subset of a
                # Check no intermediate: no node c where b < c < a
                if not any(k != i and k != j and b < nodes[k] < a
                           for k in range(len(nodes))):
                    edges.append([i, j])
    return edges

# ── Build CARTA DAG data ─────────────────────────────────────────────────────
carta_data = {}
for sample, info in sorted(kstar.items()):
    k = info['k']
    observed = get_observed_states(sample)
    if not observed:
        print(f"  SKIP {sample}: no observed states")
        continue

    chain = read_progen(sample, k)
    if not chain:
        print(f"  SKIP {sample}: no progenitor sets")
        continue

    # Deduplicate and filter to observed states
    seen = set()
    nodes = []
    for s in chain:
        s_obs = frozenset(s & observed)  # filter to observed only
        if len(s_obs) > 1 and s_obs not in seen:
            seen.add(s_obs)
            nodes.append(s_obs)

    if not nodes:
        # k=1: no intermediate progenitors, just root → leaves
        nodes = [frozenset(observed)]

    # Sort: largest first (root), then by size desc
    nodes.sort(key=lambda s: -len(s))

    # Determine root
    root_idx = 0  # largest set
    all_states = sorted(nodes[root_idx],
                        key=lambda s: STATE_ORDER.index(s) if s in STATE_ORDER else 99)

    # Compute Hasse edges between multi-state nodes
    edges = hasse_edges(nodes)

    # Depths via BFS from root
    depth = [0] * len(nodes)
    changed = True
    while changed:
        changed = False
        for p, c in edges:
            if depth[c] < depth[p] + 1:
                depth[c] = depth[p] + 1
                changed = True

    # Build node list for JSON
    node_list = []
    for i, ns in enumerate(nodes):
        ss = sorted(ns, key=lambda s: STATE_ORDER.index(s) if s in STATE_ORDER else 99)
        node_list.append({
            'id': i,
            'st': [VIEWER_NAME.get(s, s) for s in ss],     # viewer state names
            'cs': [SHORT.get(s, s) for s in ss],            # short codes
            'd': depth[i],
            'r': i == root_idx,
        })

    # Leaf states = observed states (in canonical order)
    leaves = [VIEWER_NAME.get(s, s) for s in all_states]

    # Leaf-path map: for each leaf, which node IDs are ancestors (contain it)
    lp = {}
    for s in all_states:
        vn = VIEWER_NAME.get(s, s)
        lp[vn] = [i for i in range(len(nodes)) if s in nodes[i]]

    carta_data[sample] = {
        'k': k, 'arm': info['arm'], 'nl': info['n_leaves'], 'ns': info['n_states'],
        'nodes': node_list, 'edges': edges, 'leaves': leaves, 'lp': lp,
    }
    print(f"  {sample}: k={k}, {len(all_states)} leaves, {len(nodes)} nodes, {len(edges)} edges")

carta_json = json.dumps(carta_data, separators=(',', ':'))
print(f"\nNew CARTA data: {len(carta_data)} samples, {len(carta_json):,} chars")

# ── Patch the HTML ───────────────────────────────────────────────────────────
html = VIEWER.read_text(encoding='utf-8')

# 1. Replace "carta":{...} in the DATA JSON
carta_start = html.find('"carta":{')
if carta_start < 0:
    print("ERROR: Could not find '\"carta\":{' in HTML")
    exit(1)

# Brace-match to find end
brace_count = 0
i = carta_start + len('"carta":')
while i < len(html):
    if html[i] == '{':
        brace_count += 1
    elif html[i] == '}':
        brace_count -= 1
        if brace_count == 0:
            carta_end = i + 1
            break
    i += 1
else:
    print("ERROR: Could not find end of carta data")
    exit(1)

old_carta = html[carta_start:carta_end]
new_carta = f'"carta":{carta_json}'
html = html[:carta_start] + new_carta + html[carta_end:]
print(f"Replaced carta data: {len(old_carta):,} → {len(new_carta):,} chars")

# 2. Replace renderCarta function block
new_render_fn = r'''// ==================== CARTA DAG PANEL ====================
function cartaLum(hex) {
  hex = hex.replace('#','');
  var r = parseInt(hex.substr(0,2),16)/255, g = parseInt(hex.substr(2,2),16)/255, b = parseInt(hex.substr(4,2),16)/255;
  return 0.299*r + 0.587*g + 0.114*b;
}
function renderCarta() {
  var el = document.getElementById('cartaContent');
  if (!el || !DATA.carta) return;
  var sample = document.getElementById('selSample').value;
  var info = DATA.carta[sample];
  if (!info) { el.innerHTML = '<i style="color:#889;font-size:10px">No CARTA data for ' + sample + '.</i>'; return; }
  var nodes = info.nodes, edges = info.edges, leaves = info.leaves, lp = info.lp;
  var pal = DATA.palettes.state;
  // Which leaves are highlighted?
  var hlLeaves = new Set();
  highlightSet.forEach(function(s) { if (leaves.indexOf(s) >= 0) hlLeaves.add(s); });
  var anyHl = hlLeaves.size > 0;
  // Which nodes are on a highlighted leaf's path?
  var hlNodes = new Set();
  hlLeaves.forEach(function(lf) { (lp[lf]||[]).forEach(function(nid) { hlNodes.add(nid); }); });
  // Which edges connect two highlighted nodes?
  var hlEdges = new Set();
  edges.forEach(function(e, ei) { if (hlNodes.has(e[0]) && hlNodes.has(e[1])) hlEdges.add(ei); });
  // Layout: compute rows by depth
  var maxD = 0;
  nodes.forEach(function(n) { if (n.d > maxD) maxD = n.d; });
  var byD = {};
  nodes.forEach(function(n) { if (!byD[n.d]) byD[n.d] = []; byD[n.d].push(n); });
  // Sort each row: largest first, then by id
  Object.keys(byD).forEach(function(d) {
    byD[d].sort(function(a,b) { return b.st.length - a.st.length || a.id - b.id; });
  });
  var W = 236, rowH = 26, leafH = 22, gapY = 18, padTop = 12;
  var svgH = padTop + (maxD + 1) * (rowH + gapY) + gapY + leafH + 8;
  var h = '<div class="carta-info"><b>k* = ' + info.k + '</b> &middot; ' + info.arm + ' &middot; ' + info.ns + ' states &middot; ' + info.nl + ' cells</div>';
  h += '<svg width="' + W + '" height="' + svgH + '" style="display:block;margin:0 auto">';
  // Arrow markers
  h += '<defs>';
  h += '<marker id="cah" markerWidth="5" markerHeight="4" refX="4" refY="2" orient="auto"><path d="M0,0 L5,2 L0,4Z" fill="#999"/></marker>';
  h += '<marker id="cahH" markerWidth="5" markerHeight="4" refX="4" refY="2" orient="auto"><path d="M0,0 L5,2 L0,4Z" fill="#0f2747"/></marker>';
  h += '</defs>';
  // Position each node
  var np = {};
  for (var d = 0; d <= maxD; d++) {
    var row = byD[d] || [];
    var y = padTop + d * (rowH + gapY);
    // Calculate total width of this row
    var tw = 0;
    row.forEach(function(n) {
      var bw = n.r ? 100 : Math.min(n.st.length * 18 + 8, W * 0.9);
      tw += bw;
    });
    var gap = row.length > 1 ? Math.min(12, (W - tw) / (row.length + 1)) : 0;
    var x = (W - tw - gap * (row.length - 1)) / 2;
    row.forEach(function(n) {
      var bw = n.r ? 100 : Math.min(n.st.length * 18 + 8, W * 0.9);
      np[n.id] = { x: x, y: y, w: bw, h: rowH };
      x += bw + gap;
    });
  }
  // Draw edges between progenitor nodes
  edges.forEach(function(e, ei) {
    var p = np[e[0]], c = np[e[1]];
    if (!p || !c) return;
    var x1 = p.x + p.w/2, y1 = p.y + p.h;
    var x2 = c.x + c.w/2, y2 = c.y;
    var isH = hlEdges.has(ei);
    var col = anyHl ? (isH ? '#0f2747' : '#ddd') : '#888';
    var sw = isH ? 2.2 : 0.8;
    var mk = isH ? 'cahH' : 'cah';
    if (Math.abs(x2 - x1) > 15) {
      var my = (y1 + y2) / 2;
      h += '<path d="M'+x1+','+y1+' C'+x1+','+my+' '+x2+','+my+' '+x2+','+y2+'" fill="none" stroke="'+col+'" stroke-width="'+sw+'" marker-end="url(#'+mk+')"/>';
    } else {
      h += '<line x1="'+x1+'" y1="'+y1+'" x2="'+x2+'" y2="'+y2+'" stroke="'+col+'" stroke-width="'+sw+'" marker-end="url(#'+mk+')"/>';
    }
  });
  // Draw progenitor node boxes
  nodes.forEach(function(n) {
    var p = np[n.id]; if (!p) return;
    var nH = anyHl && hlNodes.has(n.id);
    var dim = anyHl && !nH;
    if (n.r) {
      // Root: text label
      h += '<rect x="'+p.x+'" y="'+p.y+'" width="'+p.w+'" height="'+p.h+'" fill="none" stroke="'+(dim?'#ddd':'#555')+'" stroke-width="0.8" rx="3" stroke-dasharray="3,2"/>';
      h += '<text x="'+(p.x+p.w/2)+'" y="'+(p.y+p.h/2+4)+'" text-anchor="middle" font-size="9" font-weight="bold" fill="'+(dim?'#ccc':'#333')+'" opacity="'+(dim?0.4:1)+'">Progenitor ('+n.st.length+')</text>';
    } else {
      // Intermediate: split-color box
      var segW = Math.min(18, (p.w - 4) / n.st.length);
      var sx = p.x + (p.w - segW * n.st.length) / 2;
      var sy = p.y + 2, sh = p.h - 4;
      // Background rect for grouping
      h += '<rect x="'+(sx-1)+'" y="'+(sy-1)+'" width="'+(segW*n.st.length+2)+'" height="'+(sh+2)+'" fill="none" stroke="'+(dim?'#eee':'#777')+'" stroke-width="0.5" rx="2"/>';
      n.st.forEach(function(st, si) {
        var c = pal[st] || '#ccc';
        var op = dim ? 0.15 : 1;
        h += '<rect x="'+(sx+si*segW)+'" y="'+sy+'" width="'+segW+'" height="'+sh+'" fill="'+c+'" stroke="#555" stroke-width="0.3" opacity="'+op+'" rx="1"/>';
        var lm = cartaLum(c);
        var tc = lm > 0.45 ? '#111' : '#fff';
        var fs = segW < 13 ? 5.5 : 7;
        h += '<text x="'+(sx+si*segW+segW/2)+'" y="'+(sy+sh/2+3)+'" text-anchor="middle" font-size="'+fs+'" fill="'+tc+'" font-weight="bold" opacity="'+op+'" style="pointer-events:none">'+n.cs[si]+'</text>';
      });
    }
  });
  // Leaf row
  var leafY = padTop + (maxD + 1) * (rowH + gapY) + 4;
  var lW = Math.min(18, (W - 8) / leaves.length);
  var lX = (W - lW * leaves.length) / 2;
  // Edges from deepest ancestor node to each leaf
  leaves.forEach(function(lf, li) {
    var lx = lX + li * lW + lW / 2;
    var ly = leafY;
    var path = lp[lf] || [];
    // Find deepest non-root ancestor
    var deepest = -1, deepD = -1;
    path.forEach(function(nid) {
      var nd = nodes.find(function(nn) { return nn.id === nid; });
      if (nd && nd.d > deepD) { deepD = nd.d; deepest = nid; }
    });
    var isLH = hlLeaves.has(lf);
    if (deepest >= 0 && np[deepest]) {
      var pp = np[deepest];
      var isRoot = nodes.find(function(nn) { return nn.id === deepest; });
      var dashed = isRoot && isRoot.r;
      var c2 = anyHl ? (isLH ? '#0f2747' : '#ddd') : (dashed ? '#aaa' : '#777');
      var s2 = isLH ? 1.8 : 0.5;
      var dashAttr = dashed ? ' stroke-dasharray="3,2"' : '';
      h += '<line x1="'+(pp.x+pp.w/2)+'" y1="'+(pp.y+pp.h)+'" x2="'+lx+'" y2="'+ly+'" stroke="'+c2+'" stroke-width="'+s2+'"'+dashAttr+'/>';
    }
  });
  // Leaf boxes
  leaves.forEach(function(lf, li) {
    var lx = lX + li * lW;
    var ly = leafY;
    var c = pal[lf] || '#ccc';
    var isLH = hlLeaves.has(lf);
    var op = anyHl ? (isLH ? 1 : 0.15) : 1;
    h += '<rect x="'+lx+'" y="'+ly+'" width="'+(lW-1)+'" height="'+leafH+'" fill="'+c+'" stroke="'+(isLH?'#0f2747':'#555')+'" stroke-width="'+(isLH?2:0.4)+'" opacity="'+op+'" rx="2" class="carta-leaf" data-state="'+lf+'" style="cursor:pointer"/>';
    var lm = cartaLum(c);
    var tc = lm > 0.45 ? '#111' : '#fff';
    var fs = lW < 14 ? 5 : 6.5;
    h += '<text x="'+(lx+(lW-1)/2)+'" y="'+(ly+leafH/2+3)+'" text-anchor="middle" font-size="'+fs+'" fill="'+tc+'" font-weight="bold" opacity="'+op+'" style="pointer-events:none">'+lf.split('_')[0]+'</text>';
  });
  h += '</svg>';
  h += '<div class="carta-hint">Click leaf to highlight path &middot; multi-select</div>';
  el.innerHTML = h;
  // Attach click handlers to leaf boxes
  el.querySelectorAll('.carta-leaf').forEach(function(rect) {
    rect.addEventListener('click', function() {
      var st = this.dataset.state;
      var selColor = document.getElementById('selColor');
      if (selColor.value !== 'state') { selColor.value = 'state'; }
      if (highlightSet.has(st)) highlightSet.delete(st); else highlightSet.add(st);
      renderLegend(); applyLegendHighlight(); renderCarta();
    });
  });
}
'''

old_start = html.find('// ==================== CARTA')
if old_start < 0:
    old_start = html.find('function cartaLum')
    if old_start < 0:
        print("ERROR: Could not find CARTA function block")
        exit(1)

old_end = html.find('// ==================== MODE MANAGEMENT', old_start)
if old_end < 0:
    print("ERROR: Could not find MODE MANAGEMENT marker after CARTA")
    exit(1)

old_block = html[old_start:old_end]
html = html[:old_start] + new_render_fn.strip() + '\n\n' + html[old_end:]
print(f"Replaced renderCarta: {len(old_block):,} → {len(new_render_fn):,} chars")

# 3. Also update the section header in HTML
old_h3 = '<h3>CARTA progenitor chain</h3>'
new_h3 = '<h3>CARTA DAG</h3>'
# also strip any pre-existing "(v*)" suffix (e.g. from a prior "CARTA DAG (v2)" run)
import re as _re
html = _re.sub(r'<h3>CARTA DAG \(v[0-9]+\)</h3>', new_h3, html)
if old_h3 in html:
    html = html.replace(old_h3, new_h3)
    print(f"Updated panel header")

# ── Write ────────────────────────────────────────────────────────────────────
VIEWER.write_bytes(html.encode('utf-8'))
print(f"\nUpdated viewer: {VIEWER}")
print(f"File size: {VIEWER.stat().st_size:,} bytes")

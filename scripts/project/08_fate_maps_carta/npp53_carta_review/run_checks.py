import json

carta = json.load(open('carta_extracted.json'))
meta = json.load(open('state_meta.json'))
state_order = set(meta['stateOrder'])
palette = set(meta['palette'].keys())

results = {}
all_issues = []

for sample, d in sorted(carta.items()):
    issues = []
    nodes = d.get('nodes', [])
    edges = d.get('edges', [])
    leaves = d.get('leaves', [])
    lp = d.get('lp', {})

    node_by_id = {n['id']: n for n in nodes}
    node_ids = set(node_by_id.keys())

    # ---- Check 1: state name consistency ----
    bad_states = set()
    for lf in leaves:
        if lf not in state_order:
            bad_states.add(lf)
    for n in nodes:
        for s in n.get('st', []):
            if s not in state_order:
                bad_states.add(s)
    for k in lp.keys():
        if k not in state_order:
            bad_states.add(k)
    if bad_states:
        issues.append(f"State names not in stateOrder: {sorted(bad_states)}")

    # ---- Check 2: color coverage ----
    no_color = [lf for lf in leaves if lf not in palette]
    if no_color:
        issues.append(f"Leaves missing from palette (would render gray): {sorted(no_color)}")

    # ---- Check 3: lp coverage ----
    missing_lp_keys = [lf for lf in leaves if lf not in lp]
    if missing_lp_keys:
        issues.append(f"Leaves missing lp entry: {sorted(missing_lp_keys)}")
    extra_lp_keys = [k for k in lp.keys() if k not in leaves]
    if extra_lp_keys:
        issues.append(f"lp keys not in leaves list: {sorted(extra_lp_keys)}")
    bad_lp_node_refs = {}
    for lf, path in lp.items():
        bad_ids = [nid for nid in path if nid not in node_ids]
        if bad_ids:
            bad_lp_node_refs[lf] = bad_ids
    if bad_lp_node_refs:
        issues.append(f"lp paths reference invalid node ids: {bad_lp_node_refs}")

    # ---- Check 4: edge validity + Hasse containment ----
    bad_edge_refs = []
    non_superset_edges = []
    for e in edges:
        if len(e) != 2:
            bad_edge_refs.append(e)
            continue
        p, c = e
        if p not in node_ids or c not in node_ids:
            bad_edge_refs.append(e)
            continue
        pst = set(node_by_id[p]['st'])
        cst = set(node_by_id[c]['st'])
        if not cst.issubset(pst):
            non_superset_edges.append((p, c, sorted(cst - pst)))
        elif cst == pst:
            non_superset_edges.append((p, c, "EQUAL sets (not proper superset)"))
    if bad_edge_refs:
        issues.append(f"Edges referencing invalid node ids: {bad_edge_refs}")
    if non_superset_edges:
        issues.append(f"Edges violating Hasse containment (parent st must be proper superset of child st): {non_superset_edges}")

    # ---- Check 5: root consistency ----
    roots = [n for n in nodes if n.get('r') is True]
    if len(roots) != 1:
        issues.append(f"Expected exactly 1 root node, found {len(roots)}: ids={[n['id'] for n in roots]}")
    else:
        root = roots[0]
        root_st = set(root['st'])
        leaf_set = set(leaves)
        if root_st != leaf_set:
            missing_in_root = leaf_set - root_st
            extra_in_root = root_st - leaf_set
            msg = []
            if missing_in_root:
                msg.append(f"leaves not in root st: {sorted(missing_in_root)}")
            if extra_in_root:
                msg.append(f"root st has states not in leaves: {sorted(extra_in_root)}")
            issues.append(f"Root state set != union of leaf states ({'; '.join(msg)})")

    # ---- Check 6: leaf coverage across non-root nodes ----
    non_root_nodes = [n for n in nodes if not n.get('r')]
    covered = set()
    for n in non_root_nodes:
        covered |= set(n['st'])
    orphaned = [lf for lf in leaves if lf not in covered]
    if orphaned:
        issues.append(f"Orphaned leaves not contained in any non-root node: {sorted(orphaned)}")

    # ---- Check 7: node depth consistency ----
    depth_issues = []
    if roots:
        root_id = roots[0]['id']
        if node_by_id[root_id].get('d', None) != 0:
            depth_issues.append(f"Root node depth != 0 (d={node_by_id[root_id].get('d')})")
    # build parent->children from edges
    children_of = {}
    for e in edges:
        if len(e) == 2 and e[0] in node_ids and e[1] in node_ids:
            children_of.setdefault(e[0], []).append(e[1])
    for p, kids in children_of.items():
        pd = node_by_id[p].get('d')
        for c in kids:
            cd = node_by_id[c].get('d')
            if pd is None or cd is None or cd != pd + 1:
                depth_issues.append(f"Edge {p}->{c}: depth {pd} -> {cd} (expected {pd}+1)")
    if depth_issues:
        issues.append(f"Depth inconsistencies: {depth_issues}")

    results[sample] = issues
    if issues:
        all_issues.append((sample, issues))

# ---- Print report ----
passed = [s for s, iss in results.items() if not iss]
failed = [s for s, iss in results.items() if iss]

print("="*80)
print("PER-SAMPLE RESULTS")
print("="*80)
for sample in sorted(carta.keys()):
    status = "PASS" if not results[sample] else "FAIL"
    print(f"[{status}] {sample}  (arm={carta[sample].get('arm')}, nl={carta[sample].get('nl')}, ns={carta[sample].get('ns')}, k={carta[sample].get('k')}, nodes={len(carta[sample].get('nodes',[]))}, edges={len(carta[sample].get('edges',[]))}, leaves={len(carta[sample].get('leaves',[]))})")
    for issue in results[sample]:
        print(f"    - {issue}")

print()
print("="*80)
print("SUMMARY")
print("="*80)
print(f"Total samples: {len(carta)}")
print(f"PASSED: {len(passed)}  -> {sorted(passed)}")
print(f"FAILED: {len(failed)}  -> {sorted(failed)}")
print()
if all_issues:
    print("All issues found:")
    for sample, issues in all_issues:
        for issue in issues:
            print(f"  {sample}: {issue}")
else:
    print("No issues found across all 25 samples.")

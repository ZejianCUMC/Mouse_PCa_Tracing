#!/usr/bin/env python
# =============================================================================
# Step 5: k* review sheet, pooling and the six cohort fate maps (Npp53 / Npp53Rb1 / both x c10 in / out)
# =============================================================================
# Description : Uses fatemap_lib (Fig 10.3.6 rules, unit-tested pixel-identical on the original) with the
#               merged 10-state labels. Two stages:
#                 --stage kstar : chord-elbow k* per tumor and version + objective-curve sheets for ZJ's review
#                 --stage maps  : pool at the final k (kstar_final.tsv: k_star_chord, or k_review where filled)
#                                 and draw the six maps + inc-vs-exc comparison tables
#               Differences from the original, all agreed in PLAN.md: merged labels and current trees;
#               candidate edges = state pairs that share the terminal (smallest) CARTA potency set in >=2 tumors
#               of the map, i.e. the same statistic that sets edge width (replaces the hand-written CANON list;
#               the looser 'any non-root set' rule drew 21-36 edges per map and is kept as edge_rule='any');
#               fallback edge = strongest terminal / co-membership partner; labels shrunk to fit their circles;
#               no manual width boost on one edge; titles without em-dashes.
# Input       : ../step2_carta/results/{inc,exc}/{T}/ (CARTA outputs copied from HPC), ../step1_inputs/,
#               ../step4_pseudotime/pseudotime_percell.csv, reviewed10 cell_inputs_PRIVATE.csv (fade rule)
# Output      : kstar_auto.tsv, kstar_curves_{Npp53,Npp53Rb1,exc}.png, kstar_final.tsv (copy to edit),
#               maps/{map}_{inc,exc}.{png,pdf}, pooled/{map}_{ver}.json, edges/{map}_{ver}.tsv,
#               compare_inc_exc_edges.tsv, compare_inc_exc_tumors.tsv
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/scvi/bin/python
# Validated   : 2026-10-11 (library); this driver runs once CARTA outputs exist
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fatemap_lib as F

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
HERE = Path(__file__).resolve().parent; TASK = HERE.parent
RV10 = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp_tumor_merged/analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42/approved_merges/reviewed10')
DEFAULTS = {
    'results': TASK / 'step2_carta/results',          # CARTA outputs copied back from HPC
    'inputs': TASK / 'step1_inputs',                   # labels, trees, carta_tasks.tsv, build_summary.tsv
    'pseudotime': TASK / 'step4_pseudotime/pseudotime_percell.csv',
    'cells': RV10 / 'fsu_m3_merged_20261008/cell_inputs_PRIVATE.csv',   # all cells, for the fade rule (original used all h5ad cells)
    'c10': 'c10_Mes-Proliferation',                    # cycling state (ZJ 2026-10-12)
    'min_comember_tumors': 2,                          # candidate edge if the pair qualifies in >= this many tumors
    'edge_rule': 'terminal',                           # 'terminal' = shares the smallest potency set; 'any' = shares any non-root set
    'min_tumors_shown': 2,                             # ZJ 2026-10-12: show a state only if >= 2 tumors have >= 20 of its cells
    'out_suffix': '',                                  # appended to map file names (e.g. '_anyrule' for the sensitivity version)
}
ORDER = ['c1_Lum-general', 'c2_Lum-Proximal-like', 'c3_Lum-basal-hybrid', 'c4_Basal-general', 'c5_Basal-Squamous',
         'c6_Partial-EMT(Lum)', 'c7_Partial-EMT(Mes)', 'c8_Mes-ECMremodeling', 'c9_Mes-Fibroblast-like', 'c10_Mes-Proliferation']
CFG = dict(   # colours = reviewed10/relabel_20261008/STATE_MAP.csv; positions follow the Npp53 counterparts (step5_maps/crosswalk)
    colors=dict(zip(ORDER, ['#2166ac', '#4393c3', '#35978f', '#1b7837', '#7fbc41', '#E69F00', '#d95f02', '#9970ab', '#b4232a', '#762a83'])),
    labels=dict(zip(ORDER, ['Lum-\ngeneral', 'Lum-\nProximal-\nlike', 'Lum-\nbasal-\nhybrid', 'Basal-\ngeneral', 'Basal-\nSquamous',
                            'Partial-\nEMT\n(Lum)', 'Partial-\nEMT\n(Mes)', 'Mes-\nECM-\nremodel-\ning', 'Mes-\nFibro-\nblast-\nlike', 'Mes-\nProlif-\neration'])),
    pos={s: (x, round(0.10 + (y - 0.12) * 0.70 / 0.72, 4)) for s, (x, y) in zip(ORDER, [   # y squeezed to 0.10-0.80 so the
         (0.50, 0.84), (0.25, 0.74), (0.76, 0.66), (0.50, 0.55), (0.11, 0.49),             # largest circles clear the header; c2/c5 spaced so c2->c5 is visible
         (0.32, 0.32), (0.68, 0.32), (0.50, 0.12), (0.78, 0.14), (0.22, 0.14)])},
    mes={'c6_Partial-EMT(Lum)', 'c7_Partial-EMT(Mes)', 'c8_Mes-ECMremodeling', 'c9_Mes-Fibroblast-like', 'c10_Mes-Proliferation'},
    init_state='c1_Lum-general', boost_pair=None, fit_labels=True, fallback_only_isolated=True, route_all_edges=True, node_scale=1.2,
    bridge_order=['c4_Basal-general', 'c3_Lum-basal-hybrid', 'c6_Partial-EMT(Lum)', 'c1_Lum-general'])
CFG['index'] = {s: i + 1 for i, s in enumerate(ORDER)}
ARM = {'Intact': 'Intact', 'Castration': 'Cas', 'Castration+ENZA': 'Cas+Enz', 'Castration + Enzalutamide': 'Cas+Enz'}
DISPLAY = {'Intact': 'Intact', 'Cas': 'Cas', 'Cas+Enz': 'Cas+Enza'}   # display only, as the Part 2 report
ARM_ACCENT = {'Intact': '#1f5fa8', 'Cas': '#e0955f', 'Cas+Enz': '#c44e52'}   # as fate_map_datadriven_3panel.py (bands)
ARM_TEXT = {'Intact': '#1f5fa8', 'Cas': '#a85a1c', 'Cas+Enz': '#a83a3e'}       # darker header text for contrast on the bands
MAPS = {'Npp53': ['Npp53'], 'Npp53Rb1': ['Npp53Rb1'], 'Npp53_Npp53Rb1': ['Npp53', 'Npp53Rb1']}
LEGEND = ['confident direction: single arrow (pooled MWU p<0.01, effect≥0.10, tumors agree)',
          'strong link, direction unresolved: double arrow ↔ (tumors disagree)',
          'weak lean: dashed single arrow (p<0.05) / guaranteed mesenchymal bridge',
          'headless dashed - - - (below threshold / unresolved / connectivity; no arrowhead)']
LEGEND3 = LEGEND[:3]   # per-treatment panels omit below-threshold links (Fig 10.3.5 rule)

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def load_meta(P):
    bs = pd.read_csv(P['inputs'] / 'build_summary.tsv', sep='\t')
    tasks = pd.read_csv(P['inputs'] / 'carta_tasks.tsv', sep='\t')
    xw = pd.read_csv(P['inputs'] / 'state_crosswalk.csv')
    rename = dict(zip(xw.carta_state, xw.new_state))
    return bs, set(map(tuple, tasks[['version', 'tumor']].values)), rename

def res_dir(P, tasks, ver, t):
    """exc tumors without their own run reuse the inc run (inputs byte-identical)."""
    return P['results'] / (ver if (ver, t) in tasks else 'inc') / t

def kstar_stage(P):
    bs, tasks, _ = load_meta(P)
    rows, curves = [], {}
    for r in bs.itertuples():
        if r.version == 'exc' and ('exc', r.tumor) not in tasks: continue
        d = F.objective_curve(res_dir(P, tasks, r.version, r.tumor), r.tumor)
        if d.empty: rows.append(dict(tumor=r.tumor, model=r.model, treatment=r.treatment, version=r.version, status='missing')); continue
        d['obj_norm'] = d['obj'] / r.cells_kept
        curves[(r.version, r.tumor)] = d
        rows.append(dict(tumor=r.tumor, model=r.model, treatment=r.treatment, version=r.version, n_states=r.n_states,
                         n_leaves=r.cells_kept, k_run_max=int(d.k.max()), k_max_informative=int(d[~d.is_plateau].k.max()),
                         k_star_chord=F.chord_elbow(d), status='ok'))
    tab = pd.DataFrame(rows); tab.to_csv(HERE / 'kstar_auto.tsv', sep='\t', index=False)
    fin = tab.copy(); fin['k_review'] = ''; fin.to_csv(HERE / 'kstar_final.tsv', sep='\t', index=False)
    for name, sel in (('Npp53', (tab.model == 'Npp53') & (tab.version == 'inc')), ('Npp53Rb1', (tab.model == 'Npp53Rb1') & (tab.version == 'inc')),
                      ('exc', tab.version == 'exc')):
        sub = tab[sel & (tab.status == 'ok')].sort_values(['treatment', 'tumor'])
        n = len(sub); nc = 5; nr = int(np.ceil(n / nc))
        fig, axes = plt.subplots(nr, nc, figsize=(16, 2.7 * nr), constrained_layout=True, squeeze=False)
        for ax, rr in zip(axes.flat, sub.itertuples()):
            d = curves[(rr.version, rr.tumor)]; inf = d[~d.is_plateau]; pl = d[d.is_plateau]
            ax.plot(inf.k, inf.obj_norm, 'o-', color='#2166ac', lw=1.6, ms=4)
            if len(pl): ax.plot(list(inf.k.tail(1)) + list(pl.k), list(inf.obj_norm.tail(1)) + list(pl.obj_norm), ':', color='#888', lw=1.0)
            if len(inf) >= 2: ax.plot([inf.k.iloc[0], inf.k.iloc[-1]], [inf.obj_norm.iloc[0], inf.obj_norm.iloc[-1]], '-', color='#bbb', lw=0.7)
            y = d.loc[d.k == rr.k_star_chord, 'obj_norm'].values
            if len(y): ax.scatter([rr.k_star_chord], [y[0]], s=120, marker='*', facecolor='#ffde3a', edgecolor='#111', zorder=5)
            ax.set_title(f'{rr.tumor} ({DISPLAY.get(ARM.get(rr.treatment), rr.treatment)})  states={rr.n_states}  k*={rr.k_star_chord}', fontsize=9)
            ax.set_xticks(range(1, int(d.k.max()) + 1))
            ax.set_xlabel('k', fontsize=8); ax.set_ylabel('objective / leaves', fontsize=8); ax.tick_params(labelsize=7); ax.grid(alpha=0.3)
        for ax in list(axes.flat)[n:]: ax.axis('off')
        ttl = {'Npp53': 'Npp53, c10 included', 'Npp53Rb1': 'Npp53Rb1, c10 included', 'exc': 'c10 excluded, the 9 tumors whose inputs change'}[name]
        fig.suptitle(f'CARTA objective vs k, {ttl}. ★ = chord-elbow k* (k=1 is the default when a tumor has only 3 states)', fontsize=12, fontweight='bold')
        fig.savefig(HERE / f'kstar_curves_{name}.png', dpi=150, facecolor='white'); plt.close(fig)
    print(tab.groupby(['version', 'model']).k_star_chord.value_counts().to_string())

def pt_ranks(P, ver, tumors):
    """Within-tumor rank over the cells CARTA saw in this version (as build_combined_pseudotime.py)."""
    pt = pd.read_csv(P['pseudotime']); out = []
    for t in tumors:
        lab = pd.read_csv(P['inputs'] / ver / f'{t}_labels.tsv', sep='\t')
        m = pt[pt.tumor == t].merge(lab, left_on='leaf', right_on='cellBC')
        m = m.assign(sample=t, pt_rank=rankdata(m.pt_raw) / len(m))
        out.append(m[['sample', 'state', 'pt_rank']])
    return pd.concat(out)

def maps_stage(P):
    bs, tasks, rename = load_meta(P)
    kf = pd.read_csv(HERE / 'kstar_final.tsv', sep='\t', dtype={'k_review': str})
    kf['k_final'] = [int(float(r.k_review)) if str(r.k_review).strip() not in ('', 'nan') else int(r.k_star_chord) for r in kf.itertuples()]
    kmap = {(r.version, r.tumor): r.k_final for r in kf.itertuples()}
    cells = pd.read_csv(P['cells'], usecols=['sample_id', 'model', 'new_state'])
    (HERE / 'maps').mkdir(exist_ok=True); (HERE / 'pooled').mkdir(exist_ok=True); (HERE / 'edges').mkdir(exist_ok=True)
    summaries, edge_rows = {}, []
    for ver in ('inc', 'exc'):
        for r in bs[bs.version == ver].itertuples():
            k = kmap.get((ver, r.tumor), kmap.get(('inc', r.tumor)))
            sets = F.parse_progenitors(res_dir(P, tasks, ver, r.tumor) / f'{r.tumor}_progenitors_{k}.txt', rename)
            summaries[(ver, r.tumor)] = dict(F.tumor_summary(sets, ORDER), k=k, sets=sets, model=r.model, arm=ARM[r.treatment])
    for mname, models in MAPS.items():
        for ver in ('inc', 'exc'):
            tums = sorted(t for (v, t), s in summaries.items() if v == ver and s['model'] in models)
            sm = {t: summaries[(ver, t)] for t in tums}
            data = F.pool(sm, {t: sm[t]['arm'] for t in tums}, ORDER)
            (HERE / 'pooled' / f'{mname}_{ver}.json').write_text(json.dumps(data, indent=1))
            states = [s for s in ORDER if not (ver == 'exc' and s == P['c10'])]
            sub = cells[cells.sample_id.isin(tums)]
            if ver == 'exc': sub = sub[sub.new_state != P['c10']]
            n20 = (sub.groupby(['new_state', 'sample_id']).size() >= F.MINC).groupby('new_state').sum().to_dict()
            states = [s for s in states if n20.get(s, 0) >= P['min_tumors_shown']]   # hidden: >=20 cells in < 2 tumors
            cfg = dict(CFG, colors={s: CFG['colors'][s] for s in states}, mes=CFG['mes'] & set(states))
            mx = sub.groupby(['sample_id', 'new_state']).size().groupby('new_state').max().to_dict()
            faded = {s for s in states if mx.get(s, 0) < F.MINC}
            pdf_all = pt_ranks(P, ver, tums)
            dirs = F.Directions(pdf_all, states)
            com, trm = {}, {}
            for a in data['comember_by_arm']:
                for kk, v in data['comember_by_arm'][a].items(): com[kk] = com.get(kk, 0) + v
                for kk, v in data['cooccur_by_arm'][a].items(): trm[kk] = trm.get(kk, 0) + v
            basis = trm if P['edge_rule'] == 'terminal' else com
            cand = []
            for kk, v in sorted(basis.items(), key=lambda x: -x[1]):
                a_, b_ = kk.split('|')
                if v < P['min_comember_tumors'] or a_ not in states or b_ not in states: continue
                ma, mb = np.median(dirs.cells[a_]) if len(dirs.cells[a_]) else 1, np.median(dirs.cells[b_]) if len(dirs.cells[b_]) else 1
                cand.append((a_, b_) if ma <= mb else (b_, a_))
            fb = {}
            for s in states:
                part = sorted(((trm.get(kk, 0), v, kk) for kk, v in com.items() if s in kk.split('|')), reverse=True)
                part = [(v, kk) for _, v, kk in part]
                fb[s] = [x for _, kk in part for x in kk.split('|') if x != s and x in states] + [CFG['init_state']]
            cfg.update(canon=cand, fallback=fb)
            log = []
            label = {'Npp53': 'Npp53', 'Npp53Rb1': 'Npp53Rb1', 'Npp53_Npp53Rb1': 'Npp53 + Npp53Rb1'}[mname]
            vlab = 'c10 (cycling) included' if ver == 'inc' else 'c10 (cycling) excluded'
            idx_txt = 'c1–c10' if ver == 'inc' else 'c1–c9; c10 removed'
            res = F.draw_cohort_v2(cfg, data, faded, dirs, HERE / 'maps' / f'{mname}_{ver}{P["out_suffix"]}.png',
                          title=f'{label} cohort-wide consensus fate map, {vlab} (CARTA + PhyloVelo pseudotime)',
                          header=f'Cohort-wide (all {len(tums)} tumors pooled)',
                          footer=f'Circle size ∝ CARTA terminal-progenitor frequency (cohort) · number in box = state index ({idx_txt}) · states with ≥20 cells in fewer than 2 tumors are not shown.',
                          legend_labels=LEGEND, edges_log=log)
            print(mname, ver, 'cohort shown:', [x.split('_')[0] for x in states], '| smallest label pt:', min(res['min_label_pt'].values()), '| routing compromises:', res['routing_problems'])
            for e in log:
                st = dirs.stats(e['src'], e['dst'])
                cls = 'headless' if e['headless'] else ('double' if e['double'] else ('dashed' if e['dashed'] else 'solid'))
                edge_rows.append(dict(map=mname, version=ver, src=e['src'], dst=e['dst'], cls=cls, width=e['w'], **st))
            pd.DataFrame([r for r in edge_rows if r['map'] == mname and r['version'] == ver]).to_csv(HERE / 'edges' / f'{mname}_{ver}{P["out_suffix"]}.tsv', sep='\t', index=False)
            # per-treatment panels (Fig 10.3.5 rules): per-arm sizes, widths, directions and fade; same scaffold
            arm_of = {t: sm[t]['arm'] for t in tums}; faded_by_arm, dirs_by_arm = {}, {}
            hidden_by_arm = {}
            for arm in ('Intact', 'Cas', 'Cas+Enz'):
                ta = [t for t in tums if arm_of[t] == arm]
                if not ta: continue
                ca = sub[sub.sample_id.isin(ta)].groupby(['sample_id', 'new_state']).size()
                mxa = ca.groupby('new_state').max().to_dict(); n20a = (ca >= F.MINC).groupby('new_state').sum().to_dict()
                hidden_by_arm[arm] = {s for s in states if n20a.get(s, 0) < P['min_tumors_shown']}
                faded_by_arm[arm] = {s for s in states if mxa.get(s, 0) < F.MINC}
                dirs_by_arm[arm] = F.Directions(pdf_all[pdf_all['sample'].isin(ta)], states)
            print(mname, ver, 'hidden by treatment:', {a: sorted(x.split('_')[0] for x in h) for a, h in hidden_by_arm.items()})
            alog = []
            ra = F.draw_arms_v2(cfg, data, faded_by_arm, dirs_by_arm, HERE / 'maps' / f'{mname}_{ver}_by_treatment{P["out_suffix"]}.png',
                                title=f'{label} consensus fate map by treatment, {vlab} (CARTA + PhyloVelo pseudotime)',
                                footer=f'Circle size ∝ per-treatment CARTA terminal-progenitor frequency · number in box = state index ({idx_txt}) · states with ≥20 cells in fewer than 2 tumors of that treatment are not shown · below-threshold links omitted.',
                                legend_labels=LEGEND3, arm_display=DISPLAY, arm_accent=ARM_ACCENT, edges_log=alog, arm_text=ARM_TEXT, hidden_by_arm=hidden_by_arm)
            print(mname, ver, 'by treatment:', ra)
            pd.DataFrame(alog).assign(map=mname, version=ver).to_csv(HERE / 'edges' / f'{mname}_{ver}_by_treatment{P["out_suffix"]}.tsv', sep='\t', index=False)
    # inc vs exc: edges among non-c10 states, and per-tumor potency-set Jaccard (9 tumors)
    ed = pd.DataFrame(edge_rows); ed = ed[(ed.src != P['c10']) & (ed.dst != P['c10'])]
    ed['pair'] = ed.apply(lambda r: '|'.join(sorted([r.src, r.dst], key=ORDER.index)), axis=1)
    ed['arrow'] = ed.src.str.split('_').str[0] + '→' + ed.dst.str.split('_').str[0]
    cmp_e = ed.pivot_table(index=['map', 'pair'], columns='version', values=['cls', 'arrow'], aggfunc='first')
    cmp_e.columns = [f'{a}_{b}' for a, b in cmp_e.columns]
    cmp_e = cmp_e.reset_index(); cmp_e['changed'] = (cmp_e.get('cls_inc') != cmp_e.get('cls_exc')) | (cmp_e.get('arrow_inc') != cmp_e.get('arrow_exc'))
    cmp_e.to_csv(HERE / f'compare_inc_exc_edges{P["out_suffix"]}.tsv', sep='\t', index=False)
    rows = []
    for (v, t), s in summaries.items():
        if v != 'exc' or ('exc', t) not in tasks: continue
        inc = summaries[('inc', t)]
        strip = lambda sets: {frozenset(x - {P['c10']}) for x in sets} - {frozenset()}
        A, B = strip(inc['sets']), strip(s['sets'])
        rows.append(dict(tumor=t, model=s['model'], k_inc=inc['k'], k_exc=s['k'], jaccard_nonc10_sets=len(A & B) / len(A | B),
                         terminal_inc='|'.join(sorted(inc['terminal'], key=ORDER.index)), terminal_exc='|'.join(sorted(s['terminal'], key=ORDER.index))))
    pd.DataFrame(rows).to_csv(HERE / 'compare_inc_exc_tumors.tsv', sep='\t', index=False)
    print(pd.DataFrame(rows).to_string(index=False)); print(cmp_e.groupby('map').changed.agg(['sum', 'size']))

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(stage, **kwargs):
    P = {**DEFAULTS, **kwargs}
    (kstar_stage if stage == 'kstar' else maps_stage)(P)
    return {'output_path': str(HERE), 'params_used': {k: str(v) for k, v in P.items()}, 'summary': stage}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--stage', choices=['kstar', 'maps'], required=True)
    ap.add_argument('--edge_rule', choices=['terminal', 'any'], default='terminal')
    a = ap.parse_args()
    run_pipeline(a.stage, edge_rule=a.edge_rule, out_suffix='' if a.edge_rule == 'terminal' else '_anyrule')

#!/usr/bin/env python
"""Convert Zejian's canonical 75-pathway human gene-set JSON to mouse symbols
using the pyviper canonical human↔mouse ortholog table (28,795 pairs, based
on Ensembl/Entrez), NOT the case-conversion heuristic used by the AUCell
scoring script. Save as both JSON (for programmatic use) and GMT (for
clusterProfiler / fgsea)."""
import json, csv
from pathlib import Path

SRC_JSON   = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/LnCaP_MyCCaP/post_scvi/scripts/pathway_genesets_human.json')
ORTHO_CSV  = Path('/media/jrchristin/Backup10TB1/agent/2-tools/pyviper/pyviper/data/translate/human_mouse_gene_translation.csv')
OUT_DIR    = Path('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/aucell_arhi')
OUT_JSON   = OUT_DIR / 'pathway_genesets_mouse_via_ortholog.json'
OUT_GMT    = OUT_DIR / 'pathway_genesets_mouse_via_ortholog.gmt'
FALLBACK_JSON = OUT_DIR / 'pathway_genesets_mouse.json'   # case-conversion version (already saved)

# Load ortholog table → dict[human_symbol] = mouse_symbol
h2m = {}
with ORTHO_CSV.open() as f:
    reader = csv.DictReader(f)
    for row in reader:
        h, m = row['human_symbol'], row['mouse_symbol']
        if h and m: h2m[h] = m
print(f'ortholog pairs loaded: {len(h2m):,}')

def mouse_case(gene):
    return gene[0] + gene[1:].lower() if gene else gene

# Convert every human gene set → mouse
gs_h = json.loads(SRC_JSON.read_text())
gs_m_ortho    = {}
gs_m_fallback = json.loads(FALLBACK_JSON.read_text()) if FALLBACK_JSON.exists() else {}

n_total_hits = n_total_orthomiss = 0
n_recovered_by_case = 0
for pathway, genes in gs_h.items():
    mouse_genes = []
    orth_miss   = 0
    case_recover= 0
    for g in genes:
        if g in h2m:
            mouse_genes.append(h2m[g])
            n_total_hits += 1
        else:
            # Fallback: try mouse_case heuristic — some human symbols may not be in
            # the ortholog table (e.g. pseudogenes, non-canonical names) but their
            # case-converted form matches a mouse gene symbol.
            cased = mouse_case(g)
            mouse_genes.append(cased)
            orth_miss += 1
            case_recover += 1
            n_total_orthomiss += 1
            n_recovered_by_case += 1
    gs_m_ortho[pathway] = mouse_genes
    if orth_miss > 0:
        print(f'  {pathway:35s}: {len(genes)-orth_miss:3d}/{len(genes):3d} ortholog + {orth_miss:2d} case-fallback')

print(f'\ntotal genes: {n_total_hits} ortholog + {n_total_orthomiss} case-fallback = {n_total_hits+n_total_orthomiss}')
print(f'ortholog hit rate: {n_total_hits/(n_total_hits+n_total_orthomiss)*100:.1f}%')

OUT_JSON.write_text(json.dumps(gs_m_ortho, indent=2))
print(f'\nsaved JSON: {OUT_JSON}')

# Also write GMT (tab-separated: pathway\tdescription\tgene1\tgene2\t...)
with OUT_GMT.open('w') as f:
    for pathway, genes in gs_m_ortho.items():
        # Deduplicate + drop empty
        uniq = sorted(set(g for g in genes if g))
        f.write('\t'.join([pathway, 'na'] + uniq) + '\n')
print(f'saved GMT:  {OUT_GMT}')

# Diff: how different is proper-ortholog vs simple case-conversion?
if gs_m_fallback:
    print('\ndiff proper-ortholog vs case-conversion (per pathway):')
    for k in gs_m_ortho:
        s_ortho = set(gs_m_ortho[k])
        s_case  = set(gs_m_fallback.get(k, []))
        only_ortho = s_ortho - s_case
        only_case  = s_case  - s_ortho
        if only_ortho or only_case:
            print(f'  {k:35s}: +{len(only_ortho)} ortho-only, -{len(only_case)} case-only')

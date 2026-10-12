"""Evidence-reviewed, one-to-one display names; not validated malignant states."""
from pathlib import Path
import json, hashlib, shutil
import pandas as pd

out=Path(__file__).resolve().parents[1]/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42'
dest=out/'cluster_state_inference'
assert (dest/'ENRICHMENT_COMPLETE.json').exists()
snapshot=dest/'before_reporter_reorder_20260922'
if not snapshot.exists():
 snapshot.mkdir()
 for f in ['state_summary.csv','state_marker_evidence.csv','cell_state_labels_PRIVATE.csv','STATE_NAMES_COMPLETE.json','ANNOTATION_VERIFIED.json']:
  shutil.copy2(dest/f,snapshot/f)
 for f in ['candidate_r1.1_embeddings.png','candidate_r1.1_markers.png','PUBLICATION_VERIFIED.json','CANDIDATE_REVIEW_PLOTS.json']:
  shutil.copy2(out/f,snapshot/f)
 shutil.copy2(Path(__file__).resolve().parents[1]/'Npp53Rb1_tracing_report_part1_preprocessing_tumor_states.html',snapshot/'report.html')
# Original r1.1 IDs are immutable; display order follows the user's requested
# luminal-to-mesenchymal arrangement, NOT an inferred evolutionary trajectory.
spec=[
 (18,'Lum-Lcn2-Clu','Luminal-like','Lcn2 Clu Clca3a2 Slpi','Secretory-associated luminal program; strong within-sample support.'),
 (2,'Lum-Krt19-Nupr1','Luminal-like','Krt19 Nupr1 Wfdc2 Aldh1a3','Luminal epithelial/stress-associated program; not evidence of treatment resistance.'),
 (9,'Basal-Krt5-Col17a1','Basal-like','Krt5 Krt14 Col17a1 Krt17','Broad basal identity; marker enrichment supported across samples.'),
 (12,'Basal-Mecom-Ccnd2','Basal-like','Mecom Ccnd2 Trp63 Robo2','Basal identity supported; independent substate less certain; no stemness claim.'),
 (15,'Basal-lowComplexity','Basal-like, QC flagged','Grip1 Nrg1 Mecom','Median 1,257 genes / 2,141 UMIs; basal markers retained. QC qualifier, not a proven biological state.'),
 (13,'Basal-Squamous','Basal/squamous-like','Krt13 Krt10 Krt6a Cstdc5','Coherent keratinizing program; Krt13 alone is not sufficient for this name.'),
 (17,'Basal-Squamous-Hypoxia','Basal/squamous-like','Bnip3 Egln1 Ndrg1 Vegfa','Hypoxia-associated transcriptional program, not measured oxygen; distinct from cl13 in ten eligible samples.'),
 (1,'Cycling-Basal','Cycling, basal-like','Mki67 Top2a Pclaf Birc5','Strong cycling program on a basal background; retained in the cohort.'),
 (16,'Cycling-Stress','Cycling, lineage unresolved','Top2a Hells Smc2 Txnip','Only 36 cells; cycling/stress program supported, lineage and stable substate unresolved.'),
 (14,'Myeloid-MesMixed','Mixed / unresolved','C1qa Lyz2 Tyrobp Fcer1g','Myeloid/macrophage-like plus mesenchymal signals; possible admixture/doublets or unusual biology. Not confirmed malignant.'),
 (10,'Transitional-Spp1-Igfbp2','Transitional / partial-EMT-like','Igfbp2 Spp1 Tspan8 Selenop','82.2% JZ138; epithelial and mesenchymal signals coexist. Partial-EMT-like is descriptive, not a demonstrated transition.'),
 (11,'Transitional-Krt7','Transitional / partial-EMT-like','Krt7 Krt14 Dcun1d5','91.9% JZ138; hybrid keratin/mesenchymal profile. Pooled Igfbp2 enrichment reverses within samples, so omitted from the name.'),
 (6,'Mes-lowComplexity','Mesenchymal-like, QC flagged','Arhgap22 Frmd4a','Median 1,004 genes / 1,760 UMIs; weakly resolved low-complexity mesenchymal identity; not a stemness/quiescence claim.'),
 (0,'Mes-Crabp1-Id3','Mesenchymal-like','Crabp1 Id3 Aspn','89.4% JZ135; Crabp1/Id3 enriched versus mesenchymal peers. Aspn is matrix support, with inconsistent within-sample peer enrichment.'),
 (3,'Mes-Chemokine','Mesenchymal-like','Cxcl12 Ccl2 Ccl7 Cxcl1','91.0% JZ135; chemokine-associated program, limited cross-sample support.'),
 (4,'Mes-Acta2-Tagln','Mesenchymal-like','Acta2 Tagln Stmn2','Contractile/myofibroblast-like program; does not establish host fibroblast versus malignant origin.'),
 (7,'Mes-Spp1-Il33','Mesenchymal-like','Spp1 Il33 Cxcl2 Mt1','94.2% JZ135; inflammatory/stress-associated mesenchymal program.'),
 (8,'Mes-Dcn-Inhba','Mesenchymal-like','Dcn Inhba Cdh13 Gas6','65.6% JZ266; matrix/remodelling-associated program; Cdh13 within-sample support is weaker.'),
 (5,'Mes-ECM-Collagen','Mesenchymal-like','Col1a1 Col3a1 Col5a2 Ecrg4','Strong extracellular-matrix/collagen program; malignant origin is not proven by this signature.')]
# Explicit original-cluster order makes this operation idempotent. The user's
# c8/c9/c13 refer to the first naming revision, not the renumbered result.
original_order=[2,18,9,12,15,13,17,14,10,11,0,3,4,7,8,5,1,16,6]
by_original={row[0]:row for row in spec}
spec=[by_original[c] for c in original_order]
qc=pd.read_csv(out/'expression_review/r1.1_QC.csv').set_index('cluster')
labels=pd.read_csv(out/'resolution_assignments_PRIVATE.csv',index_col=0)['leiden_r1.1_seed42'].astype(int)
rows=[];evidence=[]
for rank,(cl,name,comp,genes,note) in enumerate(spec,1):
 state=f'c{rank}_{name}'
 rows.append(dict(rank=rank,state=state,original_cluster=f'cl{cl}',compartment=comp,n_cells=int(qc.loc[cl,'n']),representative_genes=genes.replace(' ', ', '),interpretation_and_limits=note))
 for comparison in ['peer','nearest']:
  tab=pd.read_csv(dest/f'cl{cl}_{comparison}_all_genes.csv.gz').set_index('gene')
  for gene in genes.split():
   r=tab.loc[gene].to_dict();r.update(state=state,gene=gene);evidence.append(r)
summary=pd.DataFrame(rows)
assert summary.n_cells.sum()==42806 and len(summary)==19 and len({r[0] for r in spec})==19
summary.to_csv(dest/'state_summary.csv',index=False)
old=pd.read_csv(snapshot/'state_summary.csv').set_index('original_cluster')
crosswalk=summary[['original_cluster','state','n_cells']].copy()
crosswalk.insert(1,'previous_state',crosswalk.original_cluster.map(old.state))
crosswalk.to_csv(dest/'state_reorder_crosswalk.csv',index=False)
pd.DataFrame(evidence).to_csv(dest/'state_marker_evidence.csv',index=False)
mapping={cl:rows[i]['state'] for i,(cl,*_) in enumerate(spec)}
cell=pd.DataFrame({'original_cluster':labels,'inferred_state':labels.map(mapping)})
assert cell.inferred_state.notna().all() and cell.inferred_state.value_counts().to_dict()==dict(zip(summary.state,summary.n_cells))
cell.to_csv(dest/'cell_state_labels_PRIVATE.csv')
for name in ['candidate_r1.1_embeddings.png','candidate_r1.1_markers.png','interpretation_summary.json']:
 target=dest/('before_naming_'+name)
 if not target.exists():shutil.copy2(out/name,target)
manifest=dict(cells=42806,clusters=19,labels='putative gene-enrichment-informed states; not independently validated',
 ordering='Swap prior c1/c2; preserve other relative order except prior c8/c9/c13 become final c17/c18/c19. Display order, not trajectory.',
 merges=0,additional_exclusions=0,seed=42,
 protected_input_sha256={f:hashlib.sha256((out/f).read_bytes()).hexdigest() for f in ['metadata_PRIVATE.csv','resolution_assignments_PRIVATE.csv','scvi_latent_PRIVATE.csv','umap_PRIVATE.csv','tsne_PRIVATE.csv']})
(dest/'STATE_NAMES_COMPLETE.json').write_text(json.dumps(manifest,indent=2))
print(summary[['state','original_cluster','n_cells']].to_string(index=False))
print('PASS: 19 one-to-one states; 42,806 cells; no merges, exclusions or re-embedding.')

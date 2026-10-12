"""Select statistically supported, nonredundant signatures without forcing state labels."""
from pathlib import Path
import json,hashlib,textwrap
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parents[1];parent=p/'analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42/approved_merges/reviewed10';d=parent/'pathway_statistics'
assert json.loads((d/'SCORING_COMPLETE.json').read_text())['status']=='PASS'
assert json.loads((d/'R_AUCELL_VERIFIED.json').read_text())['status']=='PASS'
stats={seed:pd.read_csv(d/f'statistics_seed{seed}.csv').set_index(['cluster','signature']) for seed in [42]}
s=stats[42].copy();s['q_worst_seed']=np.maximum.reduce([v.q_value.to_numpy() for v in stats.values()]);s['delta_min_seed']=np.minimum.reduce([v.median_delta.to_numpy() for v in stats.values()]);s['positive_fraction_min_seed']=np.minimum.reduce([v.positive_fraction.to_numpy() for v in stats.values()])
s['robust_enriched']=(s.q_worst_seed<.05)&(s.delta_min_seed>=.005)&(s.positive_fraction_min_seed>=.7)&s.testable
s['robust_significant']=(s.q_worst_seed<.05)&s.testable
s.reset_index().to_csv(d/'all_pathway_state_statistics.csv',index=False)
sets=json.loads((d/'genesets_mouse_scored.json').read_text());clusters=[f'c{i}' for i in range(1,11)]
collections={'Hallmark':9,'ProstateCancerRelated':5,'Reactome':7,'KEGG_Legacy':5}
# Curated display pool: biological relevance and complementary programs, not pathway-size rank.
# Inferential family remains ALL 16,190 tests; no recalculation on this display subset.
pool={
 'Hallmark':['HALLMARK_ANDROGEN_RESPONSE','HALLMARK_PROTEIN_SECRETION','HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_MYC_TARGETS_V1','HALLMARK_E2F_TARGETS','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_ANGIOGENESIS'],
 'ProstateCancerRelated':['Dehm_luminal','Hieronymus_AR_signaling','Tang_CRPC_SCL','MetC','CCP_31'],
 'Reactome':['REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY','REACTOME_KERATINIZATION','REACTOME_ECM_PROTEOGLYCANS','REACTOME_COLLAGEN_FORMATION','REACTOME_COLLAGEN_DEGRADATION','REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS','REACTOME_FOLDING_OF_ACTIN_BY_CCT_TRIC'],
 'KEGG_Legacy':['KEGG_FOCAL_ADHESION','KEGG_MAPK_SIGNALING_PATHWAY','KEGG_GLUTATHIONE_METABOLISM','KEGG_PYRUVATE_METABOLISM','KEGG_CITRATE_CYCLE_TCA_CYCLE']}
selected=[];reasons=[]
for collection,cap in collections.items():
 chosen=[]
 for name in pool[collection]:
  assert name in sets
  candidates=s.reset_index();candidates=candidates[candidates.robust_enriched&candidates.signature.eq(name)].sort_values('delta_min_seed',ascending=False)
  if candidates.empty:continue
  gs=set(sets[name])
  if any(len(gs&set(sets[v]))/len(gs|set(sets[v]))>=.5 for v in chosen):continue
  row=candidates.iloc[0];chosen.append(name)
  reasons.append(dict(signature=name,collection=collection,selected_for_state=row.cluster,q_seed42=row.q_worst_seed,median_delta_seed42=row.median_delta,positive_fraction_seed42=row.positive_fraction))
 selected.extend(chosen)
assert selected and len(selected)==len(set(selected))
chosen=pd.DataFrame(reasons);chosen.to_csv(d/'selected_pathways.csv',index=False)
means=pd.read_csv(d/'sample_balanced_auc_seed42.csv',index_col=0).loc[clusters,selected].T
z=means.sub(means.mean(axis=1),axis=0).div(means.std(axis=1,ddof=0).replace(0,1),axis=0)
means.to_csv(d/'selected_mean_aucell.csv');z.to_csv(d/'selected_heatmap_zscores.csv')
s.loc[pd.MultiIndex.from_product([clusters,selected])].reset_index().to_csv(d/'selected_pathway_state_statistics.csv',index=False)
summary=pd.read_csv(parent/'state_summary.csv').set_index('cluster')
ns=[int(s.loc[c].eligible_samples.iloc[0]) for c in clusters]
labels=clusters
short=[]
readable={'HALLMARK_TNFA_SIGNALING_VIA_NFKB':'TNFα signaling via NF-κB','HALLMARK_E2F_TARGETS':'E2F targets','HALLMARK_MYC_TARGETS_V1':'MYC targets V1','Dehm_luminal':'Dehm luminal','Hieronymus_AR_signaling':'Hieronymus AR signaling','Tang_CRPC_SCL':'Tang CRPC-SCL','MetC':'MetC','CCP_31':'CCP-31 (cell-cycle progression)','REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY':'Type I hemidesmosome assembly','REACTOME_ECM_PROTEOGLYCANS':'ECM proteoglycans','REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS':'TGF-β receptor signaling → SMADs','REACTOME_FOLDING_OF_ACTIN_BY_CCT_TRIC':'Actin folding by CCT/TRiC','KEGG_MAPK_SIGNALING_PATHWAY':'MAPK signaling','KEGG_CITRATE_CYCLE_TCA_CYCLE':'Citrate cycle (TCA cycle)'}
for n in selected:
 short.append(readable.get(n,n.replace('HALLMARK_','').replace('REACTOME_','').replace('KEGG_','').replace('_',' ').lower().capitalize()))
assert len(selected)<30
fig,ax=plt.subplots(figsize=(12,max(6,len(selected)*.29+1.8)))
im=ax.imshow(z.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
ax.set_xticks(range(10),labels,fontsize=11)
ax.set_yticks(range(len(selected)),[textwrap.fill(v,48) for v in short],fontsize=10)
ax.tick_params(length=0);ax.set_xticks(np.arange(-.5,10,1),minor=True);ax.set_yticks(np.arange(-.5,len(selected),1),minor=True);ax.grid(which='minor',color='white',linewidth=.6);ax.tick_params(which='minor',length=0)
ann=[]
for i,name in enumerate(selected):
 for j,c in enumerate(clusters):
  row=s.loc[(c,name)];mark=''
  if not row.testable:mark='NT'
  elif row.robust_significant:
   q=row.q_worst_seed;mark=('***' if q<.001 else '**' if q<.01 else '*')+('↑' if row.median_delta>0 else '↓')
  if mark:ax.text(j,i,mark,ha='center',va='center',fontsize=8 if mark=='NT' else 9,color='white' if abs(z.iloc[i,j])>1.15 else '#222222')
  ann.append(dict(cluster=c,signature=name,annotation=mark))
offset=0
for coll in collections:
 count=int(chosen.collection.eq(coll).sum())
 if count:
  if offset:ax.axhline(offset-.5,color='black',linewidth=1.3)
  title={'ProstateCancerRelated':'Prostate cancer','KEGG_Legacy':'KEGG legacy'}.get(coll,coll)
  ax.text(10.05,offset+(count-1)/2,title,va='center',ha='center',rotation=90,fontsize=10,weight='bold')
  offset+=count
pd.DataFrame(ann).to_csv(d/'heatmap_annotations.csv',index=False)
ax.set_title('Representative AUCell programs · seed 42',fontsize=13,pad=12)
fig.subplots_adjust(left=.40,right=.86,top=.94,bottom=.16)
cax=fig.add_axes([.93,.70,.014,.19]);cb=fig.colorbar(im,cax=cax);cb.set_label('Row z-score\nsample-balanced mean AUCell',fontsize=9)
fig.text(.40,.025,'* q<0.05   ** q<0.01   *** q<0.001 (BH FDR)\n↑ / ↓: higher / lower than other states within sample.\nNT: insufficient samples; color is descriptive.\nPaired samples per c1–c10: '+', '.join(map(str,ns)),fontsize=9,va='bottom')
fig.savefig(d/'representative_aucell_heatmap.png',dpi=180,bbox_inches='tight');plt.close(fig)
record=dict(status='PASS',pathways=len(selected),collections=chosen.collection.value_counts().to_dict(),eligible_samples=dict(zip(clusters,ns)),selection='Curated complementary biological-program display pool after testing all pathways; seed42 q<0.05, delta>=0.005 and positive fraction>=0.70 for seed42; Jaccard<0.50 within collection; at most 26 pathways; full results retained',stars='Seed42 BH-adjusted q-values; no additional seeds used',q_family=len(s),current_states=summary.putative_name.to_dict())
(d/'SELECTION_COMPLETE.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))

from pathlib import Path
import json,os,shutil
import anndata as ad
import pandas as pd
p=Path(__file__).resolve().parents[1];d=p/'analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42/approved_merges';n=d/'reviewed10';n.mkdir(exist_ok=True)
assert not (n/'VERIFIED.json').exists()
a=ad.read_h5ad(d/'named_graph_PRIVATE.h5ad')
a.obs['review12_cluster']=a.obs.final_cluster.astype(str);a.obs['review12_state']=a.obs.putative_state.astype(str)
groups=[(['c1'],'Lum-general'),(['c2','c4'],'Lum-Proximal-like'),(['c3'],'Lum-basal-hybrid-Inflam'),(['c5'],'Basal-general'),(['c6'],'Basal-Squamous'),(['c7'],'Partial-EMT'),(['c8'],'Mes-CAF-like'),(['c9'],'Mes-ECMremodeling'),(['c10'],'Mes-GlycogenMetabolism'),(['c11','c12'],'Mes-Proliferation')]
spec=[];cross=[]
for i,(old,name) in enumerate(groups):
 label=f'c{i+1}';mask=a.obs.review12_cluster.isin(old)
 a.obs.loc[mask,'raw_cluster']=str(i)
 spec.append([i,name,' + '.join(old),'Reviewed annotation union; not an independent lineage or functional claim'])
 for k in old:cross.append(dict(previous_cluster=k,new_cluster=label,putative_name=name,cells=int((a.obs.review12_cluster==k).sum())))
a.obs.to_csv(n/'raw_metadata_PRIVATE.csv');a.write_h5ad(n/'clustered_graph_PRIVATE.h5ad',compression='gzip')
pd.DataFrame(cross).to_csv(n/'merge_crosswalk.csv',index=False)
(n/'NAME_SPEC.json').write_text(json.dumps(spec,indent=2))
(n/'STARTED.json').write_text(json.dumps(dict(pid=os.getpid(),estimated_ram_gb=16,estimated_minutes=10,cells=len(a),groups=10,geometry='unchanged'),indent=2))
print(pd.DataFrame(cross).to_string(index=False))

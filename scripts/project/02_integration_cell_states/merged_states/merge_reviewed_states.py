"""Apply exactly the four user-approved unions, retaining premerge results."""
from pathlib import Path
import json,os
import anndata as ad
import pandas as pd
p=Path(__file__).resolve().parents[1];d=p/'analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42';n=d/'approved_merges';n.mkdir(exist_ok=True)
assert not (n/'VERIFIED.json').exists()
a=ad.read_h5ad(d/'named_graph_PRIVATE.h5ad');a.obs['premerge_cluster']=a.obs.final_cluster.astype(str);a.obs['premerge_state']=a.obs.putative_state.astype(str);a.obs['premerge_raw_cluster']=a.obs.raw_cluster.astype(str)
groups=[(['c1','c4'],'Lum-general'),(['c2'],'Lum-IFN'),(['c5','c3'],'Lum-basal-hybrid-Inflam'),(['c6'],'Lum-Proximal'),(['c7'],'Basal-general'),(['c9','c8'],'Basal-Squamous'),(['c10','c11'],'Partial-EMT'),(['c12'],'Mes-Collagen'),(['c13'],'Mes-ECMremodeling'),(['c14'],'Mes-GlycogenMetabolism'),(['c15'],'Mes-Chemotactic'),(['c16'],'Mes-Actomyosin')]
mapping={old:str(i) for i,(olds,_) in enumerate(groups) for old in olds};a.obs['raw_cluster']=a.obs.premerge_cluster.map(mapping);assert not a.obs.raw_cluster.isna().any()
spec=[];cross=[]
for i,(olds,name) in enumerate(groups):
 spec.append([i,name,' + '.join(olds),'User-approved merge' if len(olds)>1 else ('Proposed process name; unchanged membership' if i>=10 else 'Reviewed family; unchanged membership')])
 for old in olds:cross.append(dict(previous_cluster=old,new_cluster=f'c{i+1}',putative_name=name,cells=int(a.obs.premerge_cluster.eq(old).sum())))
a.obs.to_csv(n/'raw_metadata_PRIVATE.csv');a.write_h5ad(n/'clustered_graph_PRIVATE.h5ad',compression='gzip');pd.DataFrame(cross).to_csv(n/'merge_crosswalk.csv',index=False)
(n/'NAME_SPEC.json').write_text(json.dumps(spec,indent=2));(n/'STARTED.json').write_text(json.dumps(dict(pid=os.getpid(),cells=55133,groups=12,merges={'c4':'c1','c8':'c9','c3':'c5','c11':'c10'},geometry='unchanged',estimated_ram_gb=16,estimated_minutes=20),indent=2))
print(pd.DataFrame(cross).to_string(index=False))

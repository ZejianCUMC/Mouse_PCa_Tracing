# Are the calls recovered by the re-run (UMI read cutoff 1 vs 6) real?  Tests per mouse (v8 tree cells only):
#  A concordance: cell x intBC called in both runs -> same allele?
#  B read support of recovered vs original calls
#  C allele sharing: fraction of calls whose allele is carried by >=2 other cells
#  D independent phylogenetic check: on the ORIGINAL v8 tree (built without the recovered calls), does a call's allele
#    match the majority allele of that intBC among the cell's nearest relatives (>=5 cells with a call)? Compared between
#    calls present originally and calls only in the re-run.
import pandas as pd, numpy as np
from ete3 import Tree
from pathlib import Path
N=Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing")
OUT=Path(__file__).resolve().parents[1]
def calls(a):
    a=a[(a.allele!="")]
    a=a.assign(UMI=a.UMI.astype(int),readCount=a.readCount.astype(int)).sort_values("UMI",ascending=False).drop_duplicates(["cellBC","intBC"])
    return a.set_index(["cellBC","intBC"])
rows=[]
for s in ["JZ201","JZ202","JZ203","JZ204","JZ205"]:
    t=Tree(open(N/f"lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy/{s}_tree_greedy.nwk").read(),format=1); leaves=set(t.get_leaf_names())
    o=pd.read_csv(N/f"lineage_tracing/allele_tables/per_sample/{s}_allele_table.csv",dtype=str,keep_default_na=False); o=o[o.cellBC.isin(leaves)]
    n=pd.read_csv(OUT/f"inputs/{s}_allele_table.csv",dtype=str,keep_default_na=False); n=n[n.cellBC.isin(leaves)]
    co,cn=calls(o),calls(n)
    both=co.index.intersection(cn.index); rec=cn.index.difference(co.index)
    conc=(co.loc[both,"allele"]==cn.loc[both,"allele"]).mean()
    al=cn.allele; cnt=cn.reset_index().groupby(["intBC","allele"]).cellBC.nunique()
    shared=lambda idx: np.mean([cnt[(i,al[(c,i)])]>=3 for c,i in idx])
    # D: nearest-relative agreement on the original tree
    node={l.name:l for l in t.get_leaves()}; by_int={}
    for (c,i),a_ in al.items(): by_int.setdefault(i,{})[c]=a_
    def agree(idx, maxn=4000):
        rng=np.random.default_rng(1); idx=list(idx); 
        if len(idx)>maxn: idx=[idx[k] for k in rng.choice(len(idx),maxn,replace=False)]
        ok=[]
        for c,i in idx:
            obs=by_int[i]; nd=node[c]; nb=[]
            while nd.up is not None:
                nd=nd.up; nb=[x for x in nd.get_leaf_names() if x!=c and x in obs]
                if len(nb)>=5: break
            if len(nb)<5: continue
            vc=pd.Series([obs[x] for x in nb]).value_counts()
            ok.append(vc.index[0]==obs[c])
        return np.mean(ok), len(ok)
    a_orig,n_o=agree(both); a_rec,n_r=agree(rec)
    rows.append(dict(sample=s, calls_original=len(co), calls_rerun=len(cn), recovered_calls=len(rec), pct_more=round(100*len(rec)/len(co),1),
        concordance_shared_calls=round(100*conc,2),
        UMI_median_original_calls=float(cn.loc[both,"UMI"].median()), UMI_median_recovered=float(cn.loc[rec,"UMI"].median()),
        reads_per_UMI_recovered=round(float((cn.loc[rec,"readCount"]/cn.loc[rec,"UMI"]).median()),1),
        allele_shared_ge3_original=round(100*shared(both),1), allele_shared_ge3_recovered=round(100*shared(rec),1),
        tree_agree_original=round(100*a_orig,1), tree_agree_recovered=round(100*a_rec,1), n_tested=f"{n_o}/{n_r}"))
    print(rows[-1],flush=True)
d=pd.DataFrame(rows); d.to_csv(OUT/"THRESHOLD_VALIDATION.csv",index=False)
pd.set_option("display.width",250); print(d.drop(columns=["calls_original","calls_rerun"]).to_string(index=False))

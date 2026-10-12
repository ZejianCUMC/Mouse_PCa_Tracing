# For cell x intBC called in both runs but with a DIFFERENT allele: which allele agrees with the cell's nearest relatives on
# the ORIGINAL tree (neighbour alleles from the original calls, i.e. the data that built that tree, and from the re-run calls)?
# Also: is the difference in one site only, and how do the original vs re-run supports compare?
import pandas as pd, numpy as np
from ete3 import Tree
from pathlib import Path
N=Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"); OUT=Path(__file__).resolve().parents[1]
def calls(a):
    a=a[a.allele!=""].assign(UMI=lambda x:x.UMI.astype(int)).sort_values("UMI",ascending=False).drop_duplicates(["cellBC","intBC"])
    return a.set_index(["cellBC","intBC"])
rows=[]
for s in ["JZ201","JZ202","JZ203","JZ204","JZ205"]:
    t=Tree(open(N/f"lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy/{s}_tree_greedy.nwk").read(),format=1); L=set(t.get_leaf_names())
    co=calls(pd.read_csv(N/f"lineage_tracing/allele_tables/per_sample/{s}_allele_table.csv",dtype=str,keep_default_na=False).query("cellBC in @L"))
    cn=calls(pd.read_csv(OUT/f"inputs/{s}_allele_table.csv",dtype=str,keep_default_na=False).query("cellBC in @L"))
    both=co.index.intersection(cn.index); ch=[k for k in both if co.at[k,"allele"]!=cn.at[k,"allele"]]
    nsite=[sum(co.at[k,r]!=cn.at[k,r] for r in ["r1","r2","r3"]) for k in ch]
    blank_old=np.mean([any(co.at[k,r]=="" for r in ["r1","r2","r3"]) for k in ch])
    node={l.name:l for l in t.get_leaves()}
    def maj(src,c,i):
        obs={x for (x,j) in src.index if j==i} if False else None
    idx_o={}; idx_n={}
    for (c,i),a in co.allele.items(): idx_o.setdefault(i,{})[c]=a
    for (c,i),a in cn.allele.items(): idx_n.setdefault(i,{})[c]=a
    rng=np.random.default_rng(1); sel=[ch[k] for k in rng.choice(len(ch),min(3000,len(ch)),replace=False)]
    res={"o_vs_onb":[],"n_vs_onb":[],"o_vs_nnb":[],"n_vs_nnb":[]}
    for c,i in sel:
        for nbsrc,tag in [(idx_o,"onb"),(idx_n,"nnb")]:
            obs=nbsrc.get(i,{}); nd=node[c]; nb=[]
            while nd.up is not None:
                nd=nd.up; nb=[x for x in nd.get_leaf_names() if x!=c and x in obs]
                if len(nb)>=5: break
            if len(nb)<5: continue
            top=pd.Series([obs[x] for x in nb]).value_counts().index[0]
            res[f"o_vs_{tag}"].append(co.at[(c,i),"allele"]==top); res[f"n_vs_{tag}"].append(cn.at[(c,i),"allele"]==top)
    rows.append(dict(sample=s, shared_calls=len(both), changed=len(ch), pct_changed=round(100*len(ch)/len(both),1),
        changed_1site=round(100*np.mean(np.array(nsite)==1),1), changed_old_had_blank_site=round(100*blank_old,1),
        UMI_old_median=float(co.loc[ch,"UMI"].median()), UMI_new_median=float(cn.loc[ch,"UMI"].median()),
        original_allele_matches_relatives_origdata=round(100*np.mean(res["o_vs_onb"]),1), rerun_allele_matches_relatives_origdata=round(100*np.mean(res["n_vs_onb"]),1),
        original_allele_matches_relatives_rerundata=round(100*np.mean(res["o_vs_nnb"]),1), rerun_allele_matches_relatives_rerundata=round(100*np.mean(res["n_vs_nnb"]),1)))
    print(rows[-1],flush=True)
d=pd.DataFrame(rows); d.to_csv(OUT/"CHANGED_CALLS_CHECK.csv",index=False); pd.set_option("display.width",250); print(d.to_string(index=False))

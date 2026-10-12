# Pilot: are "uncut" calls that contradict their tree neighbourhood (sisters share an edit) weaker-supported?
import pandas as pd, numpy as np, sys
from ete3 import Tree
B="/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based"
CASES={"JZ201_original":(f"{B}/Npp53_tracing/lineage_tracing/allele_tables/per_sample/JZ201_allele_table.csv",f"{B}/Npp53_tracing/lineage_tracing/cassiopeia_trees_v8_rebuilt/greedy/JZ201_tree_greedy.nwk"),
       "JZ201_rerun":(f"{B}/Npp53_tracing/lineage_tracing/MJZ019_rerun_trees_20261005/inputs/JZ201_allele_table.csv",f"{B}/Npp53_tracing/lineage_tracing/MJZ019_rerun_trees_20261005/greedy/JZ201_tree_greedy.nwk"),
       "JZ204_rerun":(f"{B}/Npp53_tracing/lineage_tracing/MJZ019_rerun_trees_20261005/inputs/JZ204_allele_table.csv",f"{B}/Npp53_tracing/lineage_tracing/MJZ019_rerun_trees_20261005/greedy/JZ204_tree_greedy.nwk")}
MINN, AGREE = 5, 0.8
out=[]
for s,(af,tf) in CASES.items():
    a=pd.read_csv(af,dtype={"cellBC":str,"intBC":str},keep_default_na=False)
    t=Tree(open(tf).read(),format=1); leaves=set(t.get_leaf_names()); a=a[a.cellBC.isin(leaves)]
    # one call per cell x intBC x site; keep cell-intBC with a single allele row (ambiguous ones analysed separately)
    n_al=a.groupby(["cellBC","intBC"]).allele.transform("nunique"); amb=a[n_al>1]; a1=a[n_al==1]
    long=a1.melt(id_vars=["cellBC","intBC","UMI","readCount"],value_vars=["r1","r2","r3"],var_name="site",value_name="call")
    long=long[long.call!=""]
    long["uncut"]=long.call.str.contains(r"\[None\]")
    # neighbourhood: walk up from each leaf until >= MINN other leaves have an observed call at (intBC,site)
    call={(r.cellBC,r.intBC,r.site):r.call for r in long.itertuples()}
    keyset={}
    for (c,i,st),v in call.items(): keyset.setdefault((i,st),{})[c]=v
    anc={l.name:l for l in t.get_leaves()}
    # recurrence of an indel = number of distinct intBC-site positions where the same indel string occurs in this mouse
    E0=long[~long.uncut]; rec=E0.groupby("call").apply(lambda d: d[["intBC","site"]].drop_duplicates().shape[0]).to_dict()
    res=[]
    for r in long[long.uncut].itertuples():
        obs=keyset[(r.intBC,r.site)]; node=anc[r.cellBC]; others=[]
        while node.up is not None:
            node=node.up
            others=[l for l in node.get_leaf_names() if l!=r.cellBC and l in obs]
            if len(others)>=MINN: break
        if len(others)<MINN: continue
        vals=pd.Series([obs[o] for o in others]); top=vals.value_counts()
        ed=top[~top.index.str.contains(r"\[None\]")]
        frac_same_edit=ed.iloc[0]/len(vals) if len(ed) else 0.0
        frac_uncut=vals.str.contains(r"\[None\]").mean()
        cls="discordant" if frac_same_edit>=AGREE else ("concordant" if frac_uncut>=AGREE else "mixed")
        res.append(dict(sample=s,cls=cls,UMI=r.UMI,reads=r.readCount,clade=len(others),edit=(ed.index[0] if len(ed) else None),rec=(rec.get(ed.index[0],0) if len(ed) else None)))
    R=pd.DataFrame(res); pass
    E=long[~long.uncut]
    D=R[R.cls=="discordant"].copy(); D["recurrent"]=D.rec>=2
    for lab,g in D.groupby("recurrent"):
        u=pd.to_numeric(g.UMI); print(s,"discordant, neighbour edit", "RECURRENT (>=2 positions)" if lab else "UNIQUE (1 position)", "n=",len(g),"UMI median",u.median(),"single-UMI",round((u==1).mean(),3))
    for name,df in [("edited (all)",E),("uncut concordant",R[R.cls=="concordant"]),("uncut mixed",R[R.cls=="mixed"]),("uncut DISCORDANT",R[R.cls=="discordant"])]:
        u=pd.to_numeric(df.UMI); rd=pd.to_numeric(df["readCount" if "readCount" in df else "reads"])
        out.append(dict(sample=s,group=name,n=len(df),UMI_median=u.median(),UMI_eq1=round((u==1).mean(),3),reads_median=rd.median(),reads_per_UMI=round((rd/u).median(),1)))
    out.append(dict(sample=s,group=f"cell-intBC with >1 allele (ambiguous)",n=int(amb.groupby(['cellBC','intBC']).ngroups)))
pd.DataFrame(out).to_csv("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/scripts/barcode_error_report_20261004/data/uncut_pilot_rerun.csv",index=False)

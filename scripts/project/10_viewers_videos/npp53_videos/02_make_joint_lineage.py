#!/usr/bin/env python
"""Joint tree <-> UMAP animation, LINEAGE-SEQUENCED (same reveal pattern as {t}_evolution_lineage.mp4).
LEFT  = circular lineage tree revealed one STATE-COMMITTED LINEAGE (clone/trunk) at a time, in per-cell-type
        median PhyloVelo pseudotime order (ancestral scaffold first). Each lineage grows outward.
RIGHT = sample UMAP; each cell's dot is PAIRED to its tree leaf (same barcode) and lights up in the SAME step its
        leaf appears — so as each clone reveals on the tree, its cells light up in the UMAP.
Internal branches grey; terminal branches + leaf dots + UMAP dots coloured by major cell type.

Data: report/anim/{t}_treecoords.csv + data/bench/{t}/phylovelo/sidecar.h5ad (obsm X_umap).
Usage: make_joint_tree_umap_clade.py TUMOR [--stills] [MINTIP] [PURITY]
Output: report/anim/{t}_joint_lineage.mp4 + .gif (+ {t}_jointlin_frame_*.png)
"""
import sys, numpy as np, pandas as pd
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba_array

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Nimbus Sans","Helvetica","Arial"],
                     "pdf.fonttype":42,"ps.fonttype":42})
ROOT="/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/videos_out"; ANIM=f"{ROOT}/anim"
PAL={"Basal":"#4a7fbf","Luminal":"#26a69a","Mesenchymal":"#78909c","EMT":"#ffca28"}
GREY="#c7ccd2"; BG="#e6e8eb"; ANCHOR="__ancestral__"

def load(t, mintip=40, purity=0.60):
    d=pd.read_csv(f"{ANIM}/{t}_treecoords.csv")
    r=d.set_index("node")["x"].to_dict()
    theta=(d.set_index("node")["angle"]*np.pi/180.0).to_dict()
    par=d.set_index("node")["parent"].to_dict(); istip=d.set_index("node")["isTip"].to_dict()
    ct=d.set_index("node")["CellType"].to_dict(); lab=d.set_index("node")["label"].to_dict()
    kids=defaultdict(list)
    for n in d["node"]:
        p=par[n]
        if p!=n: kids[p].append(n)
    root=[n for n in d["node"] if par[n]==n][0]
    def dtips(s):
        o=[];st=[s]
        while st:
            u=st.pop()
            if istip[u]:o.append(u)
            else:st.extend(kids[u])
        return o
    def alln(s):
        o=[];st=[s]
        while st:
            u=st.pop(); o.append(u); st.extend(kids[u])
        return o
    chap={}; committed=[]
    def visit(n):
        tps=dtips(n)
        if len(tps)>=mintip:
            f=pd.Series([ct[x] for x in tps]).value_counts(normalize=True)
            if f.iloc[0]>=purity:
                dom=f.index[0]
                for u in alln(n): chap[u]=dom
                committed.append((n,dom,len(tps))); return
        for c in kids[n]:
            if not istip[c]: visit(c)
    visit(root)
    for n in d["node"]: chap.setdefault(n, ANCHOR)
    # order chapters by per-cell-type median PhyloVelo pseudotime (the real state pseudotime)
    tipdf=d[d["isTip"]]; type_medpt=tipdf.groupby("CellType")["Pseudotime"].median().to_dict()
    pt=d.set_index("node")["Pseudotime"].to_dict()
    anc_tips=[n for n in d["node"] if istip[n] and chap[n]==ANCHOR]
    labels=[{"label":ANCHOR,"dom":None,"med":(float(np.median([pt[n] for n in anc_tips])) if anc_tips else 0),
             "ntips":len(anc_tips),"anc":True}]
    by=defaultdict(int)
    for (n,dom,nt) in committed: by[dom]+=nt
    for dom,nt in by.items():
        labels.append({"label":dom,"dom":dom,"med":float(type_medpt.get(dom,9)),"ntips":nt,"anc":False})
    chapters=[labels[0]]+sorted([l for l in labels if not l["anc"]], key=lambda l:l["med"])
    # UMAP coords by barcode
    a=ad.read_h5ad(f"{ROOT}/sidecar/{t}_sidecar.h5ad"); a.obs_names_make_unique()
    bc=[str(x).split('-')[0] for x in a.obs_names]; xy=np.asarray(a.obsm["X_umap"],float)
    umap=pd.DataFrame({"bc":bc,"ux":xy[:,0],"uy":xy[:,1]}).drop_duplicates("bc").set_index("bc")
    return d,r,theta,par,kids,istip,ct,lab,chap,root,chapters,umap

def chap_geom(lbl, d,r,theta,par,kids,istip,ct,lab,chap,umap):
    nodes=[n for n in d["node"] if chap[n]==lbl]
    RI=[];RT=[];RTc=[]
    for n in nodes:
        p=par[n]
        if p==n: continue
        row=(r[p],r[n],theta[n])
        (RT.append(row) or RTc.append(PAL.get(ct[n],"#999999"))) if istip[n] else RI.append(row)
    ARCS=[]
    for p in nodes:
        ch=kids[p]
        if not ch: continue
        ths=[theta[c] for c in ch]; t0,t1=min(ths),max(ths); rp=r[p]
        k=max(2,int(abs(t1-t0)/(2*np.pi)*400)); tt=np.linspace(t0,t1,k)
        ARCS.append((rp,np.c_[rp*np.cos(tt),rp*np.sin(tt)]))
    tips=[n for n in nodes if istip[n]]
    lr=np.array([r[n] for n in tips]); lx=np.array([r[n]*np.cos(theta[n]) for n in tips])
    ly=np.array([r[n]*np.sin(theta[n]) for n in tips])
    cols=[PAL.get(ct[n],"#999999") for n in tips]; rgba=to_rgba_array(cols) if cols else np.empty((0,4))
    # umap coords aligned to tips
    ux=np.array([umap["ux"].get(lab[n],np.nan) for n in tips]); uy=np.array([umap["uy"].get(lab[n],np.nan) for n in tips])
    return (np.array(RI,float),np.array(RT,float),np.array(RTc,dtype=object),ARCS,lx,ly,lr,rgba,ux,uy,lbl==ANCHOR)

def radsegs(R,tt):
    if len(R)==0: return np.empty((0,2,2)),np.zeros(0,bool)
    m=R[:,0]<tt; rp=R[m,0]; rn=np.minimum(R[m,1],tt); th=R[m,2]
    p0=np.stack([rp*np.cos(th),rp*np.sin(th)],1); p1=np.stack([rn*np.cos(th),rn*np.sin(th)],1)
    return np.stack([p0,p1],1),m

def animate(t, make_stills=False, mintip=40, purity=0.60):
    d,r,theta,par,kids,istip,ct,lab,chap,root,chapters,umap=load(t,mintip,purity)
    rmax=float(d["x"].max()); ntip=int(d["isTip"].sum())
    G=[dict(meta=c, geom=chap_geom(c["label"],d,r,theta,par,kids,istip,ct,lab,chap,umap)) for c in chapters]
    allux=umap["ux"].values; alluy=umap["uy"].values
    ms_leaf=5 if ntip>1600 else (9 if ntip>=700 else 20)
    ms_dot=9 if ntip>1600 else (14 if ntip>=700 else 26)
    tlw=0.6 if ntip>1600 else (0.9 if ntip>=700 else 1.4)
    SWEEP=40; HOLD=14
    schedule=[]
    for gi in range(len(G)):
        for f in range(SWEEP): schedule.append((gi, rmax*1.03*(f+1)/SWEEP))
        for _ in range(HOLD): schedule.append((gi, rmax*1.03))
    for _ in range(20): schedule.append((len(G)-1, rmax*1.03))

    fig=plt.figure(figsize=(17,8.8)); fig.patch.set_facecolor("white")
    axL=fig.add_axes([0.005,0.06,0.47,0.85]); axL.set_aspect("equal"); axL.axis("off")
    axR=fig.add_axes([0.545,0.09,0.43,0.80]); axR.set_aspect("equal")
    ilc=LineCollection([],colors=GREY,linewidths=0.55,zorder=1); axL.add_collection(ilc)
    alc=LineCollection([],colors=GREY,linewidths=0.55,zorder=1); axL.add_collection(alc)
    tlc=LineCollection([],linewidths=tlw,zorder=2); axL.add_collection(tlc)
    leafsc=axL.scatter([],[],s=ms_leaf,zorder=3,linewidths=0)
    lim=rmax*1.08; axL.set_xlim(-lim,lim); axL.set_ylim(-lim,lim)
    axL.plot([0],[0],marker="o",ms=5,color="0.35",zorder=5)
    axL.annotate("common\nancestor",xy=(0,0),xytext=(6,-8),textcoords="offset points",fontsize=9,color="0.45",va="top")
    axL.set_title("Lineage tree  (revealed clone by clone, in pseudotime order)",fontsize=14,fontweight="bold")
    axR.scatter(allux,alluy,s=ms_dot*0.7,c=BG,zorder=1,linewidths=0)
    umsc=axR.scatter([],[],s=ms_dot,zorder=2,linewidths=0)
    padx=(np.nanmax(allux)-np.nanmin(allux))*0.05; pady=(np.nanmax(alluy)-np.nanmin(alluy))*0.05
    axR.set_xlim(np.nanmin(allux)-padx,np.nanmax(allux)+padx); axR.set_ylim(np.nanmin(alluy)-pady,np.nanmax(alluy)+pady)
    axR.set_xticks([]); axR.set_yticks([])
    for sp in ("top","right"): axR.spines[sp].set_visible(False)
    axR.spines["left"].set_color("0.6"); axR.spines["bottom"].set_color("0.6")
    axR.set_xlabel("UMAP_1",fontsize=12); axR.set_ylabel("UMAP_2",fontsize=12)
    axR.set_title("Sample UMAP  (paired cells light up with each clone)",fontsize=14,fontweight="bold")
    fig.suptitle(f"Tumour evolution — paired lineage tree & UMAP, clone by clone — {t}",fontsize=19,fontweight="bold",y=0.985)
    ptxt=fig.text(0.985,0.02,"",ha="right",va="bottom",fontsize=13,color="#111",fontweight="bold")
    cap=fig.text(0.5,0.02,"",ha="center",va="bottom",fontsize=15,fontweight="bold")
    present=[k for k in PAL if (d.loc[d['isTip'],'CellType']==k).any()]
    handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=PAL[k],markersize=9,label=k) for k in present]
    fig.legend(handles=handles,title="Major cell type / state",loc="lower left",bbox_to_anchor=(0.01,0.02),
               fontsize=10,title_fontsize=11,frameon=False,handletextpad=0.4,labelspacing=0.3)

    def draw(geom, front):
        RI,RT,RTc,ARCS,lx,ly,lr,rgba,ux,uy,is_anc=geom
        s_i,_=radsegs(RI,front); s_t,mt=radsegs(RT,front)
        arcs=[ARCS[i][1] for i in range(len(ARCS)) if ARCS[i][0]<=front]
        lm=lr<=front
        return s_i,(s_t,list(RTc[mt]) if mt.any() else []),arcs,(lx[lm],ly[lm],rgba[lm],ux[lm],uy[lm])

    def upd(fr):
        gi,front=schedule[fr]
        SI=[];ST=[];STc=[];ARC=[];LX=[];LY=[];LRG=[];UX=[];UY=[]
        for j in range(len(G)):
            f=rmax*1.03 if j<gi else (front if j==gi else -1)
            if f<0: continue
            s_i,(s_t,stc),arcs,(lxx,lyy,rg,uxx,uyy)=draw(G[j]["geom"], f)
            if len(s_i): SI.append(s_i)
            if len(s_t): ST.append(s_t); STc+=stc
            ARC+=arcs
            if len(lxx): LX.append(lxx); LY.append(lyy); LRG.append(rg); UX.append(uxx); UY.append(uyy)
        ilc.set_segments(np.concatenate(SI) if SI else np.empty((0,2,2)))
        tlc.set_segments(np.concatenate(ST) if ST else np.empty((0,2,2))); tlc.set_color(STc if STc else [])
        alc.set_segments(ARC if ARC else [])
        if LX:
            rg=np.concatenate(LRG)
            leafsc.set_offsets(np.c_[np.concatenate(LX),np.concatenate(LY)]); leafsc.set_facecolors(rg)
            umsc.set_offsets(np.c_[np.concatenate(UX),np.concatenate(UY)]); umsc.set_facecolors(rg)
            shown=sum(len(a) for a in LX)
        else:
            leafsc.set_offsets(np.empty((0,2))); umsc.set_offsets(np.empty((0,2))); shown=0
        m=G[gi]["meta"]
        if m["anc"]:
            cap.set_text(f"Ancestral scaffold + uncommitted cells  (n={m['ntips']})"); cap.set_color("#555")
        else:
            cap.set_text(f"{m['label']} clone emerges   (n={m['ntips']} cells,  PhyloVelo pseudotime {m['med']:.2f})")
            cap.set_color(PAL.get(m["dom"],"#333"))
        ptxt.set_text(f"cells shown: {shown}/{ntip}")
        return ilc,alc,tlc,leafsc,umsc,ptxt,cap

    anim=FuncAnimation(fig,upd,frames=len(schedule),interval=55,blit=False)
    mp4=f"{ANIM}/{t}_joint_lineage.mp4"; gif=f"{ANIM}/{t}_joint_lineage.gif"
    anim.save(mp4,writer=FFMpegWriter(fps=20,bitrate=3200),dpi=130)
    anim.save(gif,writer=PillowWriter(fps=16),dpi=72)
    order=" -> ".join(("ancestral" if c["anc"] else c["label"]) for c in chapters)
    print(f"wrote {mp4} and {gif}\n  clone order: {order}")
    if make_stills:
        acc=0; marks=[]
        for gi in range(len(G)):
            acc+=SWEEP+HOLD; marks.append(min(acc-2,len(schedule)-1))
        for i,fr in enumerate(marks):
            upd(fr); fig.savefig(f"{ANIM}/{t}_jointlin_frame_{i:02d}.png",dpi=110)
        print(f"wrote {len(marks)} stills")
    plt.close(fig)

if __name__=="__main__":
    a=[x for x in sys.argv[1:] if not x.startswith("--")]
    t=a[0] if a else "JZ138"; mintip=int(a[1]) if len(a)>1 else 40; pur=float(a[2]) if len(a)>2 else 0.60
    animate(t, make_stills=("--stills" in sys.argv), mintip=mintip, purity=pur)

#!/usr/bin/env python
"""Lineage-SEQUENCED circular tumour-evolution animation.
The tree is partitioned into STATE-COMMITTED LINEAGES = maximal clades with >=MINTIP tips that are >=PURITY
dominated by one cell type. Uncommitted/ancestral nodes form an 'ancestral scaffold' chapter. Chapters are
revealed one at a time in order of median pseudotime (ancestral scaffold first), each grown outward by a local
radius front; already-shown chapters stay drawn. Terminal branches + leaf dots coloured by major cell type.
This tells the real story: which cell-state lineages differentiate first vs later.  Reuses {t}_treecoords.csv.

Usage: make_tree_animation_clade.py TUMOR [--stills] [MINTIP] [PURITY]
Output: report/anim/{t}_evolution_lineage.mp4 + .gif (+ {t}_lineage_frame_*.png with --stills)
"""
import sys, numpy as np, pandas as pd
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.lines import Line2D

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Nimbus Sans","Helvetica","Arial"],
                     "pdf.fonttype":42,"ps.fonttype":42})
ANIM="/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/videos_out/anim"
PAL={"Basal":"#4a7fbf","Luminal":"#26a69a","Mesenchymal":"#78909c","EMT":"#ffca28"}
GREY="#c7ccd2"; ANC="#9aa0a6"; GAP=0.06
ANCHOR="__ancestral__"

def load(t, mintip=40, purity=0.60):
    d=pd.read_csv(f"{ANIM}/{t}_treecoords.csv")
    r=d.set_index("node")["x"].to_dict()
    theta=(d.set_index("node")["angle"]*np.pi/180.0).to_dict()  # ggtree fan angle -> matches poster orientation
    par=d.set_index("node")["parent"].to_dict(); istip=d.set_index("node")["isTip"].to_dict()
    ct=d.set_index("node")["CellType"].to_dict()
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
    # maximal committed clades
    chap={}                    # node -> chapter label (cell type) ; unset => ancestral
    committed=[]               # (root_node, domtype, ntips, purity, medpt)
    def visit(n):
        tps=dtips(n)
        if len(tps)>=mintip:
            f=pd.Series([ct[x] for x in tps]).value_counts(normalize=True)
            if f.iloc[0]>=purity:
                dom=f.index[0]
                for u in [n]+ [c for c in _all_nodes(n,kids)]: chap[u]=dom
                committed.append((n,dom,len(tps),round(f.iloc[0],2),float(np.median([r[x] for x in tps]))))
                return
        for c in kids[n]:
            if not istip[c]: visit(c)
    def _all_nodes(s,kids):
        o=[];st=[s]
        while st:
            u=st.pop(); o.append(u); st.extend(kids[u])
        return o
    visit(root)
    for n in d["node"]:
        chap.setdefault(n, ANCHOR)
    # Order chapters by PER-CELL-TYPE median PhyloVelo pseudotime (the real state pseudotime), NOT by clade
    # median (which lumps off-type cells + uses tree-branch radius and can invert e.g. Early-EMT/Partial-EMT).
    pt=d.set_index("node")["Pseudotime"].to_dict()
    tipdf=d[d["isTip"]]
    type_medpt=tipdf.groupby("CellType")["Pseudotime"].median().to_dict()
    labels=[]
    anc_tips=[n for n in d["node"] if istip[n] and chap[n]==ANCHOR]
    labels.append({"label":ANCHOR,"dom":None,
                   "med":(float(np.median([pt[n] for n in anc_tips])) if anc_tips else 0),
                   "ntips":len(anc_tips),"anc":True})
    by=defaultdict(int)
    for (n,dom,nt,pu,md) in committed: by[dom]+=nt
    for dom,nt in by.items():
        labels.append({"label":dom,"dom":dom,"med":float(type_medpt.get(dom,9)),"ntips":nt,"anc":False})
    committed_sorted=sorted([l for l in labels if not l["anc"]], key=lambda l:l["med"])
    chapters=[labels[0]]+committed_sorted            # ancestral scaffold first, then states by PhyloVelo pseudotime
    return d,r,theta,par,kids,istip,ct,chap,root,chapters

def chap_geom(lbl, d,r,theta,par,kids,istip,ct,chap):
    nodes=[n for n in d["node"] if chap[n]==lbl]
    RI=[];RT=[];RTc=[]
    nodeset=set(nodes)
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
    ly=np.array([r[n]*np.sin(theta[n]) for n in tips]); lc=np.array([PAL.get(ct[n],"#999999") for n in tips],dtype=object)
    is_anc = (lbl==ANCHOR)
    return (np.array(RI,float),np.array(RT,float),np.array(RTc,dtype=object),ARCS,lx,ly,lr,lc,is_anc)

def radsegs(R,tt):
    if len(R)==0: return np.empty((0,2,2)),np.zeros(0,bool)
    m=R[:,0]<tt; rp=R[m,0]; rn=np.minimum(R[m,1],tt); th=R[m,2]
    p0=np.stack([rp*np.cos(th),rp*np.sin(th)],1); p1=np.stack([rn*np.cos(th),rn*np.sin(th)],1)
    return np.stack([p0,p1],1),m

def animate(t, make_stills=False, mintip=40, purity=0.60):
    d,r,theta,par,kids,istip,ct,chap,root,chapters=load(t,mintip,purity)
    rmax=float(d["x"].max()); ntip=int(d["isTip"].sum())
    G=[dict(meta=c, geom=chap_geom(c["label"],d,r,theta,par,kids,istip,ct,chap)) for c in chapters]
    ms = 26 if ntip<700 else (9 if ntip<1600 else 5)
    tlw = 1.4 if ntip<700 else (0.9 if ntip<1600 else 0.6)
    SWEEP=46; HOLD=12
    schedule=[]
    for gi in range(len(G)):
        for f in range(SWEEP): schedule.append((gi, rmax*1.03*(f+1)/SWEEP))
        for _ in range(HOLD): schedule.append((gi, rmax*1.03))
    for _ in range(20): schedule.append((len(G)-1, rmax*1.03))

    fig,ax=plt.subplots(figsize=(9.6,9.9)); ax.set_aspect("equal")
    fig.subplots_adjust(top=0.9,bottom=0.075,left=0.02,right=0.98)
    fig.patch.set_facecolor("white"); ax.set_facecolor("white"); ax.axis("off")
    ilc=LineCollection([],colors=GREY,linewidths=0.55,zorder=1); ax.add_collection(ilc)
    alc=LineCollection([],colors=GREY,linewidths=0.55,zorder=1); ax.add_collection(alc)
    tlc=LineCollection([],linewidths=tlw,zorder=2); ax.add_collection(tlc)
    scat=ax.scatter([],[],s=ms,c=[],zorder=3,linewidths=0)
    lim=rmax*1.08; ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
    ax.plot([0],[0],marker="o",ms=5,color="0.35",zorder=5)
    ax.annotate("common\nancestor",xy=(0,0),xytext=(6,-8),textcoords="offset points",fontsize=9,color="0.45",va="top")
    fig.suptitle(f"Cell-state lineages emerge in pseudotime order — {t}",fontsize=17,fontweight="bold",y=0.965)
    cap=fig.text(0.5,0.032,"",ha="center",va="bottom",fontsize=14.5,fontweight="bold")
    fig.text(0.5,0.006,"State-committed lineages revealed in order of per-cell-type median PhyloVelo pseudotime",
             ha="center",va="bottom",fontsize=9.5,color="0.45")
    present=[k for k in PAL if (d.loc[d['isTip'],'CellType']==k).any()]
    handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=PAL[k],markersize=9,label=k) for k in present]
    ax.legend(handles=handles,title="Major cell type / state",loc="upper left",bbox_to_anchor=(-0.02,1.0),
              fontsize=9.5,title_fontsize=10.5,frameon=False,handletextpad=0.4,labelspacing=0.32)

    def draw(geom, front):
        RI,RT,RTc,ARCS,lx,ly,lr,lc,is_anc=geom
        s_i,_=radsegs(RI,front); s_t,mt=radsegs(RT,front)
        arcs=[ARCS[i][1] for i in range(len(ARCS)) if ARCS[i][0]<=front]
        lm=lr<=front
        return s_i,is_anc,(s_t,list(RTc[mt]) if mt.any() else []),arcs,(lx[lm],ly[lm],list(lc[lm]) if lm.any() else [])

    def upd(fr):
        gi,front=schedule[fr]
        SIg=[];SIa=[];ST=[];STc=[];ARC=[];LX=[];LY=[];LC=[]
        for j in range(len(G)):
            f = rmax*1.03 if j<gi else (front if j==gi else -1)
            if f<0: continue
            si,is_anc,(st,stc),arcs,(lxx,lyy,lcc)=draw(G[j]["geom"], f)
            (SIa if is_anc else SIg).append(si) if len(si) else None
            if len(st): ST.append(st); STc+=stc
            ARC+=arcs
            if len(lxx): LX.append(lxx); LY.append(lyy); LC+=lcc
        ilc.set_segments(np.concatenate(SIa+SIg) if (SIa or SIg) else np.empty((0,2,2)))
        tlc.set_segments(np.concatenate(ST) if ST else np.empty((0,2,2))); tlc.set_color(STc if STc else [])
        alc.set_segments(ARC if ARC else [])
        if LX:
            scat.set_offsets(np.c_[np.concatenate(LX),np.concatenate(LY)]); scat.set_color(LC)
        else:
            scat.set_offsets(np.empty((0,2))); scat.set_color([])
        m=G[gi]["meta"]
        if m["anc"]:
            cap.set_text(f"Ancestral scaffold + uncommitted cells  (n={m['ntips']})"); cap.set_color("#555")
        else:
            cap.set_text(f"{m['label']} lineage emerges   (n={m['ntips']} cells,  PhyloVelo pseudotime {m['med']:.2f})")
            cap.set_color(PAL.get(m["dom"],"#333"))
        return ilc,alc,tlc,scat,cap

    anim=FuncAnimation(fig,upd,frames=len(schedule),interval=55,blit=False)
    mp4=f"{ANIM}/{t}_evolution_lineage.mp4"; gif=f"{ANIM}/{t}_evolution_lineage.gif"
    anim.save(mp4,writer=FFMpegWriter(fps=20,bitrate=2800),dpi=140)
    anim.save(gif,writer=PillowWriter(fps=16),dpi=80)
    order=" -> ".join(("ancestral" if c["anc"] else c["label"]) for c in chapters)
    print(f"wrote {mp4} and {gif}\n  chapter order (by median PT): {order}")
    if make_stills:
        acc=0; marks=[]
        for gi in range(len(G)):
            acc+=SWEEP+HOLD; marks.append(min(acc-2,len(schedule)-1))
        for i,fr in enumerate(marks):
            upd(fr); fig.savefig(f"{ANIM}/{t}_lineage_frame_{i:02d}.png",dpi=110)
        print(f"wrote {len(marks)} QC stills")
    plt.close(fig)

if __name__=="__main__":
    a=[x for x in sys.argv[1:] if not x.startswith("--")]
    t=a[0] if a else "JZ135"; mintip=int(a[1]) if len(a)>1 else 40; pur=float(a[2]) if len(a)>2 else 0.60
    animate(t, make_stills=("--stills" in sys.argv), mintip=mintip, purity=pur)

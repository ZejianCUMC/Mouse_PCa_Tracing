import h5py, scipy.sparse as sp, scipy.io as sio, numpy as np, pandas as pd, os, gzip, sys
lib=sys.argv[1]
BASE="/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing"; RD=f"{BASE}/seurat_objs/redemux"; PP=f"{BASE}/preprocessing/cellranger_prep"; C=f"{BASE}/counts"
FO=f"{BASE}/final_objs"; OUT=f"{FO}/{lib}_mtx"; os.makedirs(OUT,exist_ok=True)
# demuxmix singlets
d=pd.read_csv(f"{RD}/{lib}_demuxmix.csv"); d=d[d.Type=="singlet"][["cell","HTO"]]
# scDblFinder singlets (intersection)
s=pd.read_csv(f"{RD}/{lib}_scdblfinder.csv")
s_sing=set(s.loc[s.scdbl_class=="singlet","cell"])
n_dmm=len(d); d=d[d.cell.isin(s_sing)]
final=list(d.cell)
# cellbender GEX
with h5py.File(f"{C}/{lib}/outs/cellbender_output_file_filtered.h5","r") as f:
    M=sp.csc_matrix((f["matrix/data"][:],f["matrix/indices"][:],f["matrix/indptr"][:]),
                    shape=tuple(int(x) for x in f["matrix/shape"][:]))
    name=[x.decode() for x in f["matrix/features/name"][:]]
    fid =[x.decode() for x in f["matrix/features/id"][:]] if "id" in f["matrix/features"] else name
    ftype=[x.decode() for x in f["matrix/features/feature_type"][:]]
    bc=[x.decode() for x in f["matrix/barcodes"][:]]
gex=np.array(ftype)=="Gene Expression"
bidx={b:i for i,b in enumerate(bc)}; cols=[bidx[c] for c in final]
Msub=M[gex,:][:,cols].tocoo(); Msub.data=np.round(Msub.data).astype(int)
gname=np.array(name)[gex]; gid=np.array(fid)[gex]
sio.mmwrite(f"{OUT}/matrix.mtx", Msub); os.system(f"gzip -f {OUT}/matrix.mtx")
with gzip.open(f"{OUT}/features.tsv.gz","wt") as g:
    for i in range(len(gname)): g.write(f"{gid[i]}\t{gname[i]}\tGene Expression\n")
with gzip.open(f"{OUT}/barcodes.tsv.gz","wt") as g:
    for c in final: g.write(c+"\n")
# meta
key=pd.read_csv(f"{PP}/{lib}_hashtag_key.csv")
smap=dict(zip(key.hashtag,key.jz_sample)); tmap=dict(zip(key.hashtag,key.treatment))
mmap=dict(zip(key.hashtag,key.model)); rmap=dict(zip(key.hashtag,key.route_project))
scr=dict(zip(s.cell,s.scdbl_score))
meta=pd.DataFrame({"cell":final}); meta["library"]=lib
hmap=dict(zip(d.cell,d.HTO)); meta["hashtag"]=meta.cell.map(hmap)
meta["sample_id"]=meta.hashtag.map(smap); meta["treatment"]=meta.hashtag.map(tmap)
meta["model"]=meta.hashtag.map(mmap); meta["route_project"]=meta.hashtag.map(rmap)
meta["scdbl_score"]=meta.cell.map(scr)
meta.to_csv(f"{FO}/{lib}_final_meta.csv",index=False)
print(f"{lib}: demuxmix singlets={n_dmm} -> after scDblFinder AND = {len(final)} final ({100*len(final)/n_dmm:.0f}% kept)")

import h5py, numpy as np, scipy.sparse as sp, sys, csv, os
lib=sys.argv[1]
C="/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/counts"; out="/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/seurat_objs/redemux"; os.makedirs(out,exist_ok=True)
fn=f"{C}/{lib}/outs/cellbender_output_file_filtered.h5"
with h5py.File(fn,"r") as f:
    data=f["matrix/data"][:]; indices=f["matrix/indices"][:]; indptr=f["matrix/indptr"][:]
    shape=tuple(f["matrix/shape"][:])
    names=[x.decode() if isinstance(x,bytes) else x for x in f["matrix/features/name"][:]]
    ftype=[x.decode() if isinstance(x,bytes) else x for x in f["matrix/features/feature_type"][:]]
    bcs=[x.decode() if isinstance(x,bytes) else x for x in f["matrix/barcodes"][:]]
M=sp.csc_matrix((data,indices,indptr),shape=shape)          # features x cells
hto_idx=[i for i,t in enumerate(ftype) if "Antibody" in t]
H=np.asarray(M[hto_idx,:].todense()).astype(int)            # 6 x ncells
hto_names=[names[i] for i in hto_idx]
with open(f"{out}/{lib}_cbHTO.csv","w",newline="") as g:
    w=csv.writer(g); w.writerow(["barcode"]+hto_names)
    for j,bc in enumerate(bcs):
        w.writerow([bc]+H[:,j].tolist())
print(lib,"wrote cbHTO",H.shape,"tags",hto_names,"e.g.bc",bcs[0])

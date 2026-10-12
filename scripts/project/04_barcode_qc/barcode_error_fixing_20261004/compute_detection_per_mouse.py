import pandas as pd, json, glob, os, numpy as np
B="/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based"
rows=[]
for proj,d in [("Npp53",f"{B}/Npp53_tracing/figures/barcode_heatmaps_v2_20261004"),("Np53Rb1",f"{B}/Npp53Rb1_Tracing/figures/barcode_heatmaps_v2_20261004")]:
    J=pd.read_csv(f"{d}/JOINS.csv")
    for f in sorted(glob.glob(f"{d}/_display_matrices/*_display_matrix.csv")):
        s=os.path.basename(f).replace("_display_matrix.csv","")
        cm=pd.read_csv(f,index_col=0,dtype=str); keys=json.loads(open(f.replace("_display_matrix.csv","_columns.json")).read())
        miss=(cm=="-1")
        grp={}
        for j,k in enumerate(keys): grp.setdefault(k[:-2],[]).append(j)
        cov=pd.Series({g:(~miss.iloc[:,js]).any(axis=1).mean() for g,js in grp.items()})
        top=cov.sort_values(ascending=False).head(20)          # 20 best-detected intBCs per mouse (fair across mice)
        cols=[j for g in top.index for j in grp[g]]
        rows.append(dict(project=proj,sample=s,cells=len(cm),intBCs=len(grp),median_intBC_detect=round(cov.median(),3),
            top20_mean_detect=round(top.mean(),3),missing_frac_top20=round(miss.iloc[:,cols].values.mean(),3),
            intBCs_per_cell_median=int((~miss.iloc[:,[grp[g][0] for g in grp]]).sum(axis=1).median()),
            treatment=J.set_index("sample").file.get(s,"").split("_")[1]))
d=pd.DataFrame(rows); pd.set_option("display.width",200)
d.to_csv("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/scripts/barcode_error_report_20261004/data/detection_per_mouse.csv",index=False); print(len(d))


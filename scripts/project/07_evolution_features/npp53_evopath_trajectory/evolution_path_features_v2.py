#!/usr/bin/env python
"""Genes / pathways / regulons / FSU modules following the REAL tumor evolution
trajectory (corrected excluded PhyloVelo pseudotime, best-run-per-sample).

Uses the ACTUAL Npp53 pyscenic/AUCell results:
  regulons = analysis/pyscenic_out/auc_mtx.csv        (172 regulons, per Npp53 cell)
  pathways = analysis/aucell_arhi/aucell_scores.csv   (75 Hallmark, per Npp53 cell)
  FSU      = Hotspot module gene lists scored on Npp53 expression (sc.tl.score_genes)
  genes    = Npp53 log-norm expression
(Barcodes in these matrices are the Npp53 cells' garbled names 'MJZ0xx_<bc>-1'
 with sample_id=JZ*; mapped to (sample,bc16) via the h5ad, then to pseudotime.)

Spearman(feature, pseudotime): cohort-wide + per treatment arm + per fate-BRANCH."""
import warnings; warnings.filterwarnings('ignore')
import os, re, numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr, rankdata

ROOT='/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing'
AN=f'{ROOT}/analysis'; OUT=os.environ['OUTDIR']; PT=os.environ['PTCSV']
BRANCH={'Luminal-Distal':'LumAR','Luminal-Proximal':'LumAR','Luminal-ARdep':'LumAR','Luminal-CR':'LumAR','Luminal-Squamous':'LumAR',
        'Luminal-Persistent':'BasalPersist','Basal-like':'BasalPersist',
        'Partial-EMT':'Mes','Mes-Quiescent':'Mes','Mes-Collagen':'Mes','Mes-Stemness':'Mes','Mes-myCAF':'Mes'}

pt=pd.read_csv(PT); pt['k']=pt['sample']+'|'+pt['cellBC']
k2pt=dict(zip(pt['k'],pt['pt_rank'])); k2arm=dict(zip(pt['k'],pt['arm'])); k2st=dict(zip(pt['k'],pt['state']))
print(f'pt cells: {len(pt)}',flush=True)

a=sc.read_h5ad(f'{ROOT}/data/npp53_merge/npp53_edited.h5ad')
a.obs['bc']=[re.search(r'[ACGTN]{14,18}',x).group(0) for x in a.obs_names]
a.obs['k']=a.obs['sample_id'].astype(str)+'|'+a.obs['bc']
a.obs['pt']=a.obs['k'].map(k2pt); a.obs['arm']=a.obs['k'].map(k2arm); a.obs['st']=a.obs['k'].map(k2st)
a=a[a.obs['pt'].notna()].copy()
a.obs['branch']=a.obs['st'].map(BRANCH)
obs2pt=dict(zip(a.obs_names,a.obs['pt'])); obs2arm=dict(zip(a.obs_names,a.obs['arm'])); obs2br=dict(zip(a.obs_names,a.obs['branch']))
print(f'matched cells: {a.n_obs}',flush=True)

# ---- real Npp53 AUCell matrices ----
def load_mtx(fp,prefix):
    d=pd.read_csv(fp,index_col=0)
    d=d[d.index.isin(obs2pt)]
    d.columns=[f'{prefix}:{c}' for c in d.columns]
    return d
reg=load_mtx(f'{AN}/pyscenic_out/auc_mtx.csv','REG')
pth=load_mtx(f'{AN}/aucell_arhi/aucell_scores.csv','PATH')
print(f'regulons {reg.shape}, pathways {pth.shape}',flush=True)
# FSU modules scored on expression
fsu_df=pd.read_csv(f'{AN}/hotspot_canonical/cohort_genes_by_module_Intact.csv')
present=set(a.var_names)
for m in sorted(fsu_df.Module.unique()):
    g=[x for x in fsu_df[fsu_df.Module==m]['Gene'] if x in present]
    sc.tl.score_genes(a,g,score_name=f'FSU:M{m}',use_raw=False)
fsu=a.obs[[c for c in a.obs.columns if c.startswith('FSU:')]].copy(); fsu.index=a.obs_names
sets=reg.join(pth,how='outer').join(fsu,how='outer').reindex(a.obs_names)
print(f'combined feature-set matrix: {sets.shape}',flush=True)

ptv=np.array([obs2pt[o] for o in a.obs_names],float)
armv=np.array([obs2arm[o] for o in a.obs_names]); brv=np.array([obs2br[o] for o in a.obs_names])
def corr_rows(df):
    rows=[]
    for c in df.columns:
        v=df[c].values.astype(float); ok=np.isfinite(v)
        if ok.sum()<50: continue
        rho,p=spearmanr(v[ok],ptv[ok]); rec=dict(feature=c,type=c.split(':')[0],rho_cohort=round(rho,3),p_cohort=p,n=int(ok.sum()))
        for arm in ['Intact','Cas','Cas+Enz']:
            mm=ok&(armv==arm); rec[f'rho_{arm}']=round(spearmanr(v[mm],ptv[mm])[0],3) if mm.sum()>=50 else np.nan
        for br in ['LumAR','BasalPersist','Mes']:
            mm=ok&(brv==br); rec[f'rho_{br}']=round(spearmanr(v[mm],ptv[mm])[0],3) if mm.sum()>=50 else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)
rs=corr_rows(sets).sort_values('rho_cohort',key=lambda s:-s.abs())
rs.to_csv(f'{OUT}/evopath_sets_ranked.csv',index=False)
print('\n=== TOP feature-sets following trajectory (|rho| cohort) ===',flush=True)
print(rs.head(30).to_string(index=False),flush=True)
for t in ['REG','PATH','FSU']:
    print(f'\n--- top {t} UP(late) / DOWN(early) ---',flush=True)
    sub=rs[rs.type==t]
    print('UP  :',list(sub.sort_values('rho_cohort',ascending=False).head(8)['feature']),flush=True)
    print('DOWN:',list(sub.sort_values('rho_cohort').head(8)['feature']),flush=True)

# ---- genes ----
X=a.X; X=X.toarray() if hasattr(X,'toarray') else np.asarray(X)
frac=(X>0).mean(0); keep=frac>=0.05
Xk=X[:,keep].astype(np.float32); gk=np.array(a.var_names)[keep]
print(f'\ngenes after filter: {Xk.shape[1]}',flush=True)
Xr=np.apply_along_axis(rankdata,0,Xk).astype(np.float32); ptr=rankdata(ptv).astype(np.float32)
Xr-=Xr.mean(0); ptc=ptr-ptr.mean()
rho=np.where((np.sqrt((Xr**2).sum(0))*np.sqrt((ptc**2).sum()))>0, (Xr.T@ptc)/(np.sqrt((Xr**2).sum(0))*np.sqrt((ptc**2).sum())+1e-12),0)
# per-branch gene corr
def gcorr(mask):
    xr=np.apply_along_axis(rankdata,0,Xk[mask]).astype(np.float32); pr=rankdata(ptv[mask]).astype(np.float32)
    xr-=xr.mean(0); pc=pr-pr.mean()
    return np.where((np.sqrt((xr**2).sum(0))*np.sqrt((pc**2).sum()))>0,(xr.T@pc)/(np.sqrt((xr**2).sum(0))*np.sqrt((pc**2).sum())+1e-12),0)
g=pd.DataFrame({'gene':gk,'rho_cohort':np.round(rho,3),'frac':np.round(frac[keep],3)})
for br in ['LumAR','BasalPersist','Mes']:
    mask=brv==br
    if mask.sum()>=100: g[f'rho_{br}']=np.round(gcorr(mask),3)
g=g.sort_values('rho_cohort',key=lambda s:-s.abs())
g.to_csv(f'{OUT}/evopath_genes_ranked.csv',index=False)
print('\n=== genes UP along trajectory ===',flush=True); print(g.sort_values('rho_cohort',ascending=False).head(25).to_string(index=False),flush=True)
print('\n=== genes DOWN along trajectory ===',flush=True); print(g.sort_values('rho_cohort').head(25).to_string(index=False),flush=True)
print('\nDONE',flush=True)

"""Exclude historical r0.3/seed19 cl5 and cl8; refit projections and validators.

No scVI training, HVG reselection, marker-driven ranking, or automatic annotation.
Run in a new output directory; never overwrite the historical analysis.
"""
import os
import tempfile
os.environ.setdefault('NUMBA_CACHE_DIR', tempfile.mkdtemp(prefix='rb1-filtered-numba-'))
os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='rb1-filtered-mpl-'))
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ[key] = '4'
import gzip, hashlib, json, shutil, subprocess, sys, time
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread, mmwrite
from sklearn.manifold import TSNE
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

project = Path(__file__).resolve().parents[1]
source = project / 'analysis/tumor_states_20260922'
out = source / 'exclude_cl5_cl8_seed42'
rscript = '/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript'
start = time.time()
manifest = {}

def save(name, value):
    (out / name).write_text(json.dumps(value, indent=2))

def status(stage, **kwargs):
    value = {**manifest, 'stage': stage, 'updated': time.strftime('%Y-%m-%dT%H:%M:%S'), **kwargs}
    save('STATUS.json', value)
    print(stage, kwargs, flush=True)

def run_r(script, logname, timeout, *args):
    command = [rscript, str(project / 'scripts' / script), str(out), *args]
    status(script, command=command, log=logname)
    with (out / logname).open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
    if result.returncode:
        raise RuntimeError(f'{script} exited {result.returncode}; no automatic retry')

def main():
    global manifest
    out.mkdir(exist_ok=True)
    assert not (out / 'RUNNING.json').exists(), 'Existing run: refuse overwrite/retry'
    manifest = {'cells': 42806, 'original_cells': 43057, 'excluded_cells': 251,
        'exclusion_clusters': {'5': 203, '8': 48}, 'screening_resolution': 0.3,
        'screening_seed': 19, 'primary_seed': 42, 'sensitivity_seeds': [19, 888],
        'new_scvi_fit': False, 'new_hvg_selection': False, 'new_choir_fit': True,
        'new_umap': True, 'new_tsne': True, 'neighbours': 15, 'top_n': 5,
        'started': time.strftime('%Y-%m-%dT%H:%M:%S'), 'pid': os.getpid(),
        'source_sha256': {n: hashlib.sha256((source/n).read_bytes()).hexdigest() for n in
            ['scvi_latent_PRIVATE.csv', 'reviewed_labels_PRIVATE.csv', 'hvg_counts_PRIVATE.mtx.gz', 'features.txt']}}
    save('RUNNING.json', manifest)
    status('preparing filtered cohort')
    latent = pd.read_csv(source/'scvi_latent_PRIVATE.csv', index_col=0)
    meta = pd.read_csv(source/'metadata_PRIVATE.csv', index_col=0)
    labels = pd.read_csv(source/'reviewed_labels_PRIVATE.csv', index_col=0, dtype=str)
    assert latent.shape == (43057, 30) and latent.index.is_unique
    assert latent.index.equals(meta.index) and latent.index.equals(labels.index)
    keep = ~labels.cluster.isin(['5', '8'])
    assert labels.loc[~keep, 'cluster'].value_counts().to_dict() == {'5': 203, '8': 48}
    audit = meta[['sample_id', 'run_id', 'treatment']].copy()
    audit['screening_cluster_seed19_r0.3'] = labels.cluster
    audit['retained'] = keep
    audit['reason'] = np.where(labels.cluster.eq('5'), 'Excluded: ribo-high QC-suspect (user decision)',
        np.where(labels.cluster.eq('8'), 'Excluded: ionocyte-like, outside analysis scope', 'Retained'))
    audit.to_csv(out/'cell_retention_PRIVATE.csv')
    audit.groupby(['sample_id','run_id','treatment','reason']).size().rename('cells').reset_index().to_csv(out/'retention_by_sample.csv',index=False)
    counts = mmread(gzip.open(source/'hvg_counts_PRIVATE.mtx.gz', 'rb')).tocsc()
    barcodes = (source/'barcodes_PRIVATE.txt').read_text().splitlines()
    assert barcodes == latent.index.tolist() and counts.shape == (3000, 43057)
    counts = counts[:, keep.to_numpy()]
    meta = meta.loc[keep].copy(); latent = latent.loc[keep].copy()
    assert len(meta) == 42806 and meta.sample_id.nunique() == 22 and meta.run_id.nunique() == 6
    assert (labels.loc[keep, 'cluster'].eq('2')).sum() == 3138
    with gzip.open(out/'hvg_counts_PRIVATE.mtx.gz','wb') as f:
        mmwrite(f, counts)
    shutil.copy2(source/'features.txt', out/'features.txt')
    (out/'barcodes_PRIVATE.txt').write_text('\n'.join(latent.index)+'\n')
    meta.to_csv(out/'metadata_PRIVATE.csv'); latent.to_csv(out/'scvi_latent_PRIVATE.csv')
    labels.loc[keep].to_csv(out/'screening_labels_PRIVATE.csv')
    save('PREPARED.json', manifest)

    status('rebuilding neighbours, UMAP and tSNE')
    a = ad.AnnData(np.zeros((len(meta),1), dtype=np.float32), obs=meta.copy())
    a.obsm['X_scVI'] = latent.to_numpy(dtype=np.float32)
    sc.pp.neighbors(a, use_rep='X_scVI', n_neighbors=15, metric='euclidean', random_state=42)
    sc.tl.umap(a, random_state=42)
    umap = a.obsm['X_umap']
    pd.DataFrame(umap,index=meta.index,columns=['UMAP_1','UMAP_2']).to_csv(out/'umap_PRIVATE.csv')
    tsne = TSNE(n_components=2,perplexity=30,init='pca',learning_rate='auto',max_iter=1000,
                random_state=42,n_jobs=4,verbose=1).fit_transform(latent.to_numpy())
    pd.DataFrame(tsne,index=meta.index,columns=['tSNE_1','tSNE_2']).to_csv(out/'tsne_PRIVATE.csv')
    assert np.isfinite(umap).all() and np.isfinite(tsne).all()
    treatment = meta.treatment.replace({'Castration+Enzalutamide':'Castration+ENZA'})
    order = ['Intact','Castration','Castration+ENZA']
    assert set(treatment) == set(order)
    colours = {'Intact':'#278c74','Castration':'#e99b27','Castration+ENZA':'#9569b7'}
    plt.rcParams.update({'font.size':15,'axes.titlesize':19})
    fig,axes = plt.subplots(1,2,figsize=(16,7),layout='constrained')
    plot_order = np.random.default_rng(42).permutation(len(meta))
    for ax,xy,title,prefix in zip(axes,[umap,tsne],['Rebuilt UMAP','Rebuilt tSNE'],['UMAP','tSNE']):
        ax.scatter(xy[plot_order,0],xy[plot_order,1],c=treatment.map(colours).to_numpy()[plot_order],s=7,linewidths=0,rasterized=True)
        for t in order:
            ax.scatter([],[],c=colours[t],s=35,label=f'{t} (n={int(treatment.eq(t).sum()):,})')
        ax.set(title=title,xlabel=prefix+'_1',ylabel=prefix+'_2',xticks=[],yticks=[])
    axes[1].legend(loc='upper center',bbox_to_anchor=(.5,-.1),ncol=1,frameon=False,fontsize=12)
    fig.suptitle('42,806 retained cells • fixed scVI latent space • seed 42',fontsize=21)
    fig.savefig(out/'filtered_treatment_embeddings.png',dpi=180);plt.close(fig)
    save('EMBEDDINGS_COMPLETE.json',{**manifest,'tsne_perplexity':30,'tsne_max_iter':1000,
        'umap_min_dist':0.5,'treatment_counts':treatment.value_counts().to_dict(),'elapsed_seconds':time.time()-start})

    status('Leiden sweep on filtered graph')
    assign = pd.DataFrame(index=meta.index); sweep=[]
    for resolution in np.round(np.arange(.1,1.21,.1),1):
        solutions=[]
        for seed in [42,19,888]:
            key=f'leiden_r{resolution:.1f}_seed{seed}'
            sc.tl.leiden(a,resolution=float(resolution),key_added=key,flavor='igraph',
                         n_iterations=-1,directed=False,random_state=seed)
            assign[key]=a.obs[key].astype(str).to_numpy();solutions.append(assign[key])
        sweep.append({'resolution':float(resolution),'ARI_42_19':adjusted_rand_score(solutions[0],solutions[1]),
            'ARI_42_888':adjusted_rand_score(solutions[0],solutions[2]),'ARI_19_888':adjusted_rand_score(solutions[1],solutions[2])})
        print('Leiden',resolution,[x.nunique() for x in solutions],flush=True)
    assign.to_csv(out/'resolution_assignments_PRIVATE.csv')
    pd.DataFrame(sweep).to_csv(out/'resolution_sweep.csv',index=False)
    run_r('tumor_states_choir.R','choir_execution.log',21600)
    choir = pd.read_csv(out/'choir_labels_PRIVATE.csv',index_col=0,dtype=str).loc[meta.index,'CHOIR_cluster']
    assert choir.notna().all()
    ranks=[]
    for resolution in np.round(np.arange(.1,1.21,.1),1):
        values=[assign[f'leiden_r{resolution:.1f}_seed{seed}'] for seed in [42,19,888]]
        amis=[adjusted_mutual_info_score(choir,x,average_method='arithmetic') for x in values]
        ranks.append({'resolution':float(resolution),'clusters':values[0].nunique(),
                      'AMI':amis[0],'AMI_seed19':amis[1],'AMI_seed888':amis[2]})
    rank=pd.DataFrame(ranks).sort_values(['AMI','resolution'],ascending=[False,True])
    rank['rank']=range(1,13)
    rank.to_csv(out/'AMI_ranking.csv',index=False);rank.head(5).to_csv(out/'scSHC_candidates.csv',index=False)
    save('RANKED.json',{**manifest,'candidates':rank.head(5).resolution.tolist()})
    run_r('tumor_states_scshc.R','scshc_execution.log',21600,'42')
    summary=pd.read_csv(out/'scSHC_summary.csv')
    assert len(summary)==5 and np.allclose(summary.resolution,rank.head(5).resolution)
    table=rank.head(5)[['rank','resolution','AMI','clusters']].merge(summary,on='resolution',validate='one_to_one')
    assert (table.clusters==table.input_clusters).all()
    table[['rank','resolution','AMI','input_clusters','defended_clusters','retained_fraction']].to_csv(out/'top5_evidence.csv',index=False)
    save('COMPLETE.json',{**manifest,'elapsed_seconds':time.time()-start,'status':'completed'})
    subprocess.run([sys.executable,str(project/'scripts/plot_top5_evidence.py'),str(out)],check=True,timeout=120)
    status('computation complete; awaiting report review')
    print(table.to_string(index=False),flush=True)

if __name__=='__main__':
    try:
        main()
    except Exception as error:
        if out.exists():
            status('FAILED; no automatic retry',error=repr(error))
        raise

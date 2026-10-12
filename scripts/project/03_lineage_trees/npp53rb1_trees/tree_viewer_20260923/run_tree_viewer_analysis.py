"""State-independent tree variants, LBI or PhyloVelo, one bounded task."""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import sys,json,time,random,argparse,hashlib,platform,resource,traceback
from pathlib import Path
import numpy as np,pandas as pd

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False))

def variant(casdir,sample,method,out):
    import cassiopeia as cas
    import networkx as nx
    from ete3 import Tree
    src=casdir/'results'/sample
    data=json.loads((src/'character_matrix.json').read_text())
    cm=pd.DataFrame([[tuple(v) if isinstance(v,list) else v for v in row] for row in data['states']],index=data['cells'],columns=data['characters'])
    ambiguous=any(isinstance(v,tuple) for row in cm.values for v in row)
    assert not ambiguous or method=='NJ', 'Ambiguous input requires a supported solver; never silently resolve'
    solvers=dict(NJ=lambda:cas.solver.NeighborJoiningSolver(add_root=True),Spectral=cas.solver.SpectralSolver,
        MaxCut=cas.solver.MaxCutGreedySolver,Percolation=lambda:cas.solver.PercolationSolver(joining_solver=cas.solver.VanillaGreedySolver()))
    tree=cas.data.CassiopeiaTree(character_matrix=cm)
    solvers[method]().solve(tree,collapse_mutationless_edges=False,logfile=str(out/'solver.log'))
    graph=tree.get_tree_topology();assert nx.is_arborescence(graph)
    assert set(tree.leaves)==set(cm.index) and len(tree.leaves)==len(cm)
    target=out/'tree.nwk';target.write_text(tree.get_newick(record_branch_lengths=True)+'\n')
    assert set(Tree(target.read_text(),format=1).get_leaf_names())==set(cm.index)
    return dict(cells=len(cm),characters=cm.shape[1],tree_sha256=sha(target),character_sha256=sha(src/'character_matrix.json'),cassiopeia=cas.__version__,no_cell_exclusions=True)

def lbi(casdir,sample,out):
    import cassiopeia as cas
    from cassiopeia.tools.fitness_estimator import LBIJungle
    src=casdir/'results'/sample
    qc=json.loads((src/'COMPLETE.json').read_text());nwk=src/f'{sample}_tree_greedy.nwk'
    assert sha(nwk)==qc['tree_sha256']
    tree=cas.data.CassiopeiaTree(tree=str(nwk))
    LBIJungle(random_seed=42).estimate_fitness(tree)
    rows=[]
    for node in tree.nodes:
        if node==tree.root:continue
        value=float(tree.get_attribute(node,'fitness'));assert np.isfinite(value)
        rows.append(dict(node=node,is_leaf=tree.is_leaf(node),fitness=value))
    df=pd.DataFrame(rows);leaves=df[df.is_leaf]
    assert set(leaves.node)==set(tree.leaves)
    target=out/'lbi_fitness.csv';df.to_csv(target,index=False)
    return dict(cells=len(leaves),tree_sha256=sha(nwk),fitness_sha256=sha(target),cassiopeia=cas.__version__,
        fitness_min=float(leaves.fitness.min()),fitness_max=float(leaves.fitness.max()),
        interpretation='LBIJungle mean_fitness, Greedy topological lengths; relative proxy, not observed growth rate')

def phylovelo(root,casdir,sample,out):
    import phylovelo as pv
    from phylovelo import pseudo_time as pp
    import inspect,heapq
    from Bio import Phylo
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components
    from sklearn.neighbors import kneighbors_graph
    inp=root/'inputs'/sample
    manifest=json.loads((root/'PREPARED.json').read_text())
    rec=next(r for r in manifest['inputs'] if r['sample']==sample)
    assert all(sha(inp/f)==h for f,h in rec['files'].items())
    x=sparse.load_npz(inp/'counts.npz').astype(float)
    cells=pd.read_csv(inp/'cells.csv',index_col=0).index
    genes=pd.read_csv(inp/'genes.csv',index_col=0).index
    emb=pd.read_csv(inp/'umap.csv',index_col=0).loc[cells]
    assert ((x>0).sum(axis=1)>=100).all(), 'Expression QC would exclude cells: needs review'
    keep=np.asarray((x>0).sum(axis=0)).ravel()>=10;x=x[:,keep];genes=genes[keep]
    totals=np.asarray(x.sum(axis=1)).ravel();assert (totals>0).all()
    norm_target=float(np.median(totals))
    x=sparse.diags(norm_target/totals)@x;x.data=np.log1p(x.data)
    nwk=casdir/'results'/sample/f'{sample}_tree_greedy.nwk'
    qc=json.loads((nwk.parent/'COMPLETE.json').read_text());assert sha(nwk)==qc['tree_sha256']
    tree=Phylo.read(str(nwk),'newick');tips={t.name for t in tree.get_terminals()};assert tips==set(cells)
    depths=tree.depths(unit_branch_lengths=True);tipdepth={t.name:depths[t] for t in tree.get_terminals()}
    sd=pv.scData(x_normed=pd.DataFrame(x.toarray(),index=cells,columns=genes),phylo_tree=tree)
    sd.Xdr=emb.copy()
    cutoff=.97;pool=min(100,len(cells)//3)
    print('inference',sample,len(cells),len(genes),'cutoff',cutoff,'pool',pool,flush=True)
    pv.velocity_inference(sd,[tipdepth[c] for c in cells],cutoff=cutoff,target='x_normed')
    assert np.isfinite(sd.velocity).all() and np.any(np.asarray(sd.velocity)!=0),'No finite nonzero velocity'
    pd.DataFrame({'gene':genes,'velocity':sd.velocity}).to_csv(out/'gene_velocity.csv',index=False)
    pv.velocity_embedding(sd,target='x_normed',n_neigh=pool)
    assert np.isfinite(sd.velocity_embeded).all(),'Non-finite embedding velocity'
    # Preserve Prim's root/edge rule but use a heap; explicitly reject disconnected
    # graphs instead of accepting the upstream 1e99 sentinel as pseudotime.
    graph_review={}
    def prim(graph,start):
        visited={start};heap=[];path=[]
        def push(u):
            for v,w in graph[u].items():
                if v!=u and v not in visited and np.isfinite(w) and w<1e98:heapq.heappush(heap,(float(w),u,v))
        push(start)
        while heap and len(visited)<len(graph):
            _,u,v=heapq.heappop(heap)
            if v in visited:continue
            visited.add(v);path.append((u,v));push(v)
        graph_review['reachable_cells']=len(visited)
        assert len(visited)==len(graph), 'Disconnected finite velocity graph; pseudotime withheld'
        return path
    pp.prim=prim
    # Predeclared per-sample k selection solely on UMAP graph connectivity.
    candidates=sorted(set(min(k,len(cells)-1) for k in [100,200,400,800]))
    k=None;components=[]
    for trial in candidates:
        g=kneighbors_graph(emb,n_neighbors=trial-1,mode='connectivity')
        n=connected_components(g,directed=True,connection='strong',return_labels=False);components.append(dict(k=trial,components=int(n)))
        if n==1:k=trial;break
    assert k is not None
    pv.calc_phylo_pseudotime(sd,n_neighbors=k,r_sample=1)
    pt=np.asarray(sd.phylo_pseudotime);assert len(pt)==len(cells) and np.isfinite(pt).all()
    nuniq=int(len(np.unique(np.round(pt,4))))
    pd.DataFrame({'pseudotime':pt},index=cells).to_csv(out/'pseudotime_candidate.csv',index_label='cell')
    pd.DataFrame(sd.velocity_embeded,index=cells,columns=['velocity_UMAP1','velocity_UMAP2']).to_csv(out/'velocity_umap.csv',index_label='cell')
    valid=nuniq>=20 and float(np.ptp(pt))>0
    if valid:pd.DataFrame({'pseudotime':pt},index=cells).to_csv(out/'pseudotime.csv',index_label='cell')
    return dict(cells=len(cells),genes=len(genes),tree_sha256=sha(nwk),cutoff=cutoff,alpha=.05,target='x_normed',
        normalization_target=norm_target,gene_min_cells=10,cell_min_genes=100,velocity_pool=pool,n_neighbors=k,
        r_sample=1,seed=42,embedding='existing final shared UMAP; no new state labels',graph_checks=components,
        **graph_review,unique_pseudotime_rounded4=nuniq,pseudotime_valid=valid,
        inference_source_sha256=sha(inspect.getfile(pv.velocity_inference)),
        prim_implementation='heap-based finite Prim; root index0; full reachability required',
        selection='No state-based parameter choice or orientation; no automatic failed-job retry')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--sample',required=True);ap.add_argument('--method',required=True);args=ap.parse_args()
    random.seed(42);np.random.seed(42);sys.setrecursionlimit(20000);assert os.environ.get('PYTHONHASHSEED')=='42'
    root=args.root;casdir=root/'subsample_cassiopeia' if args.sample.startswith('JZ135_n') else root.parent/'cassiopeia_eight_20260923_seed42'
    assert args.sample in json.loads((root/'PREPARED.json').read_text())['samples']
    out=root/'results'/args.sample/args.method;out.mkdir(parents=True,exist_ok=False)
    start=time.time();base=dict(sample=args.sample,method=args.method,seed=42,host=platform.node(),job_id=os.environ.get('SLURM_JOB_ID'),started_unix=start)
    put(out/'STARTED.json',base)
    try:
        if args.method=='LBI':result=lbi(casdir,args.sample,out)
        elif args.method=='PhyloVelo':result=phylovelo(root,casdir,args.sample,out)
        else:result=variant(casdir,args.sample,args.method,out)
        put(out/'COMPLETE.json',dict(**base,**{k:v for k,v in result.items() if k!='seed'},status='COMPLETE',elapsed_seconds=time.time()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        print('COMPLETE',args.sample,args.method,result,flush=True)
    except Exception as e:
        put(out/'FAILED.json',dict(**base,error=str(e),traceback=traceback.format_exc(),elapsed_seconds=time.time()-start));raise

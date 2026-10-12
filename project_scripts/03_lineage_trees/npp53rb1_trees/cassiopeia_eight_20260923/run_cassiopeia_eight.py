"""One primary Npp53-style VanillaGreedy tree; preserve ambiguous states."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
import argparse,sys,time,json,hashlib,random,platform,resource
from pathlib import Path
import numpy as np,pandas as pd,networkx as nx
import cassiopeia as cas

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def self_test():
    # Verify native ambiguous-state support, retaining both conflicting edits.
    a=pd.DataFrame([['a','bc','A','NONE','NONE'],['a','bc','B','NONE','NONE'],
        ['b','bc','A','NONE','NONE'],['c','bc','C','NONE','NONE'],['d','bc','NONE','NONE','NONE']],
        columns=['cellBC','intBC','r1','r2','r3'])
    cm,priors,_=cas.pp.convert_alleletable_to_character_matrix(a,allele_rep_thresh=.99,collapse_duplicates=True)
    assert isinstance(cm.loc['a'].iloc[0],tuple) and len(cm.loc['a'].iloc[0])==2
    tree=cas.data.CassiopeiaTree(character_matrix=cm,priors=priors)
    solver=cas.solver.VanillaGreedySolver();assert solver.allow_ambiguous
    solver.solve(tree,logfile='/dev/null')
    assert set(tree.leaves)=={'a','b','c','d'}
    print('SELF-TEST PASS: native conflicting-state tuples and all leaves preserved; Cassiopeia',cas.__version__,flush=True)

def run(root,sample):
    started=time.time();random.seed(42);np.random.seed(42);sys.setrecursionlimit(20000)
    assert os.environ.get('PYTHONHASHSEED')=='42'
    manifest=json.loads((root/'PREPARED.json').read_text());assert sample in manifest['samples']
    source=root/'inputs'/sample;qc=json.loads((source/'INPUT_QC.json').read_text())
    assert sha(source/'allele_table.csv.gz')==qc['allele_table_sha256']
    assert sha(source/'cell_metadata.csv')==qc['metadata_sha256']
    output=root/'results'/sample;output.mkdir(parents=True,exist_ok=False)
    (output/'STARTED.json').write_text(json.dumps(dict(sample=sample,host=platform.node(),pid=os.getpid(),
        job_id=os.environ.get('SLURM_JOB_ID'),array_job_id=os.environ.get('SLURM_ARRAY_JOB_ID'),
        array_task_id=os.environ.get('SLURM_ARRAY_TASK_ID'),seed=42,started_unix=started),indent=2))
    a=pd.read_csv(source/'allele_table.csv.gz',dtype={'cellBC':str,'intBC':str})
    meta=pd.read_csv(source/'cell_metadata.csv',index_col=0)
    assert a.cellBC.nunique()==len(meta)==qc['matched_cells'] and set(a.cellBC)==set(meta.index)
    print(f'[{sample}] convert {len(a)} allele rows / {len(meta)} cells',flush=True)
    cm,priors,state_map=cas.pp.convert_alleletable_to_character_matrix(a,allele_rep_thresh=.99,
        collapse_duplicates=True,cut_sites=['r1','r2','r3'])
    assert len(cm)==len(meta) and set(cm.index)==set(meta.index) and cm.shape[1]>0
    cm=cm.loc[meta.index]
    observed=cm.apply(lambda col:col.map(lambda s:any(v!=-1 for v in s) if isinstance(s,tuple) else s!=-1)).sum(axis=1)
    assert (observed>0).all(), 'Cells with no observed characters after filtering need review'
    ambiguous=cm.apply(lambda col:col.map(lambda s:isinstance(s,tuple) and len(s)>1))
    # JSON preserves tuples as lists; CSV is a readable companion, not a lossless loader.
    cm.to_csv(output/'character_matrix.csv.gz',compression={'method':'gzip','mtime':0})
    (output/'character_matrix.json').write_text(json.dumps(dict(cells=cm.index.tolist(),characters=cm.columns.tolist(),states=cm.values.tolist())))
    (output/'state_to_indel.json').write_text(json.dumps(state_map,indent=2))
    (output/'priors.json').write_text(json.dumps(priors,indent=2))
    meta.assign(observed_characters=observed,ambiguous_characters=ambiguous.sum(axis=1)).to_csv(output/'leaf_metadata.csv',index_label='cell')
    cm_qc=dict(cells=len(cm),characters=cm.shape[1],ambiguous_entries=int(ambiguous.values.sum()),
        cells_with_ambiguity=int(ambiguous.any(axis=1).sum()),missing_fraction=float(1-observed.sum()/cm.size),
        unique_profiles=int(cm.drop_duplicates().shape[0]),mutation_priors_used=bool(priors))
    (output/'CHARACTERS_COMPLETE.json').write_text(json.dumps(cm_qc,indent=2))
    print(f'[{sample}] solving {cm.shape}; ambiguous entries={cm_qc["ambiguous_entries"]}',flush=True)
    t=time.time();tree=cas.data.CassiopeiaTree(character_matrix=cm,priors=priors,cell_meta=meta)
    solver=cas.solver.VanillaGreedySolver();assert solver.allow_ambiguous
    solver.solve(tree,collapse_mutationless_edges=False,logfile=str(output/'solver.log'))
    graph=tree.get_tree_topology()
    assert nx.is_arborescence(graph) and set(tree.leaves)==set(meta.index) and len(tree.leaves)==len(meta)
    text=tree.get_newick(record_branch_lengths=True)
    from ete3 import Tree
    parsed=Tree(text,format=1)
    assert len(parsed.get_leaf_names())==len(meta) and set(parsed.get_leaf_names())==set(meta.index)
    target=output/f'{sample}_tree_greedy.nwk';target.write_text(text+'\n')
    pd.DataFrame(graph.edges,columns=['parent','child']).to_csv(output/'edges.csv',index=False)
    result=dict(status='COMPLETE',sample=sample,seed=42,solver='VanillaGreedySolver',cassiopeia_version=cas.__version__,
        numpy_version=np.__version__,pandas_version=pd.__version__,**cm_qc,
        leaves=len(tree.leaves),nodes=len(graph),polytomies=sum(d>2 for _,d in graph.out_degree()),
        maximum_depth=max(nx.single_source_shortest_path_length(graph,tree.root).values()),
        solve_seconds=time.time()-t,elapsed_seconds=time.time()-started,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        tree_sha256=sha(target),allele_input_sha256=qc['allele_table_sha256'],metadata_input_sha256=qc['metadata_sha256'],
        branch_lengths='default topological lengths; not calibrated time',allele_rep_thresh=.99,
        ambiguity_policy=qc['ambiguity_policy'],no_new_cell_exclusions=True,
        job_id=os.environ.get('SLURM_JOB_ID'),host=platform.node())
    (output/'COMPLETE.json').write_text(json.dumps(result,indent=2))
    print(f'PASS [{sample}]: {len(meta)} exact RNA/TSL leaves; Newick round-trip; {result["elapsed_seconds"]:.1f}s',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run-dir',type=Path);parser.add_argument('--sample');parser.add_argument('--self-test',action='store_true');args=parser.parse_args()
    if args.self_test:self_test()
    else:
        assert args.run_dir and args.sample
        run(args.run_dir,args.sample)

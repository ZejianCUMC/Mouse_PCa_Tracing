"""Independent leaf/topology/metadata audit; incomplete jobs never count as success."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,json,hashlib,csv
from pathlib import Path
from ete3 import Tree

def verify(root):
    manifest=json.loads((root/'PREPARED.json').read_text());rows=[]
    for sample in manifest['samples']:
        folder=root/'results'/sample
        if not (folder/'COMPLETE.json').exists():
            rows.append(dict(sample=sample,status='INCOMPLETE',cells=None));continue
        result=json.loads((folder/'COMPLETE.json').read_text())
        target=folder/f'{sample}_tree_greedy.nwk'
        assert hashlib.sha256(target.read_bytes()).hexdigest()==result['tree_sha256']
        with (root/'inputs'/sample/'cell_metadata.csv').open() as f:meta=list(csv.DictReader(f))
        cells=[r['cell'] for r in meta]
        parsed=Tree(target.read_text(),format=1);leaves=parsed.get_leaf_names()
        assert len(leaves)==len(set(leaves))==len(cells)==result['leaves']
        assert set(leaves)==set(cells) and all(r['sample_id']==sample.split('_n')[0] for r in meta)   # intbcfix: JZ135_n* subsets carry sample_id JZ135
        assert result['seed']==42 and result['allele_rep_thresh']==.99
        assert result['no_new_cell_exclusions'] and result['solver']=='VanillaGreedySolver'
        with (folder/'leaf_metadata.csv').open() as f:outmeta=list(csv.DictReader(f))
        assert [{k:r[k] for k in meta[0]} for r in outmeta]==meta
        chars=json.loads((folder/'character_matrix.json').read_text())
        assert chars['cells']==cells and len(chars['characters'])==result['characters']
        assert len(chars['states'])==len(cells)
        ambiguous=sum(isinstance(v,list) and len(v)>1 for row in chars['states'] for v in row)
        assert ambiguous==result['ambiguous_entries']
        rows.append(dict(sample=sample,status='VERIFIED',cells=len(leaves),characters=result['characters'],
            cells_with_ambiguity=result['cells_with_ambiguity'],polytomies=result['polytomies'],
            elapsed_seconds=result['elapsed_seconds'],tree_sha256=result['tree_sha256']))
    complete=all(r['status']=='VERIFIED' for r in rows)
    result=dict(all_verified=complete,samples=rows,held_runs=manifest['held_runs'],
        total_verified_leaves=sum(r['cells'] or 0 for r in rows))
    (root/'TREE_REVIEW.json').write_text(json.dumps(result,indent=2))
    with (root/'tree_summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['sample','status','cells','characters','cells_with_ambiguity','polytomies','elapsed_seconds','tree_sha256']);writer.writeheader();writer.writerows(rows)
    print(json.dumps(result,indent=2),flush=True)
    if not complete:raise SystemExit('Not all eight trees completed: review logs; no automatic retry')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('run_dir',type=Path);args=parser.parse_args();verify(args.run_dir)

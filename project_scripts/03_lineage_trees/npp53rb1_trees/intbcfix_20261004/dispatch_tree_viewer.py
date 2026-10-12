import json,os,sys
from pathlib import Path
root=Path(sys.argv[1]);tasks=json.loads(Path(sys.argv[2]).read_text());task=tasks[int(os.environ['SLURM_ARRAY_TASK_ID'])]
env='PhyloVelo' if task['method']=='PhyloVelo' else 'cassiopeia_env'
python='/groups/ms3625_gp/zw2994/mamba/envs/'+env+'/bin/python'
if task['method']=='Greedy':
    os.execv(python,[python,str(root/'scripts/run_cassiopeia_eight.py'),'--run-dir',str(root/'cassiopeia'),'--sample',task['sample']])
os.execv(python,[python,str(root/'scripts/run_tree_viewer_analysis.py'),'--root',str(root),'--sample',task['sample'],'--method',task['method']])

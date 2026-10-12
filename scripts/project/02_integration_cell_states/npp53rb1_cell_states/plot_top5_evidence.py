"""Two-panel primary resolution evidence; seed sensitivity is supplementary."""
from pathlib import Path
import sys
import json
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'analysis/tumor_states_20260922/top5_seed42'
assert (out/'COMPLETE.json').exists()
rank=pd.read_csv(out/'AMI_ranking.csv').sort_values('resolution')
tested=pd.read_csv(out/'top5_evidence.csv').sort_values('rank')
sweep=pd.read_csv(out/'resolution_sweep.csv')
assert len(tested)==5
plt.rcParams.update({'font.size':14,'axes.titlesize':18})
fig,axes=plt.subplots(1,2,figsize=(16,6),layout='constrained')
ax=axes[0]
ax.plot(rank.resolution,rank.AMI,'o-',color='#416c94',label='Primary seed 42')
selected=rank[rank.resolution.isin(tested.resolution)]
ax.scatter(selected.resolution,selected.AMI,s=100,facecolors='none',edgecolors='#ce6025',linewidths=2,label='Top five → sc-SHC')
for row in selected.itertuples():ax.annotate(f'#{row.rank}',(row.resolution,row.AMI),xytext=(0,10),textcoords='offset points',ha='center',fontsize=11)
ax.set(xlabel='Leiden resolution',ylabel='Adjusted mutual information',title='Agreement with raw CHOIR reference')
ax.set_xticks(rank.resolution);ax.margins(y=.25);ax.legend(loc='lower right',fontsize=11)
ax=axes[1];x=np.arange(5)
ax.bar(x-.18,tested.input_clusters,.36,label='Input Leiden clusters',color='#5b8db8')
ax.bar(x+.18,tested.defended_clusters,.36,label='After sc-SHC',color='#dd9856')
for i,row in enumerate(tested.itertuples()):
    ax.text(i-.18,row.input_clusters+.2,str(row.input_clusters),ha='center')
    ax.text(i+.18,row.defended_clusters+.2,str(row.defended_clusters),ha='center')
ax.set(xticks=x,xticklabels=[f'{r:.1f}\nAMI #{int(k)}' for r,k in zip(tested.resolution,tested['rank'])],
       xlabel='Candidate resolution (AMI rank order)',ylabel='Number of clusters',title='Run-aware sc-SHC: top five candidates')
ax.set_ylim(0,max(tested.input_clusters)*1.3);ax.legend(fontsize=11)
fig.suptitle('Resolution review • primary seed 42 • fixed tumour cohort and scVI latent space',fontsize=19)
fig.savefig(out/'state_resolution_review.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(15,5.5),layout='constrained')
for col,label in [('AMI','42 (primary)'),('AMI_seed19','19'),('AMI_seed888','888')]:
    axes[0].plot(rank.resolution,rank[col],'o-',label=label)
axes[0].set(xlabel='Leiden resolution',ylabel='AMI against CHOIR',title='Sensitivity of CHOIR agreement');axes[0].legend(title='Seed')
for col in ['ARI_42_19','ARI_42_888','ARI_19_888']:
    axes[1].plot(sweep.resolution,sweep[col],'o-',label=col.replace('ARI_','').replace('_',' vs '))
axes[1].set(xlabel='Leiden resolution',ylabel='Adjusted Rand index',ylim=(0,1.03),title='Agreement between Leiden seeds');axes[1].legend()
fig.savefig(out/'state_seed_sensitivity.png',dpi=180);plt.close(fig)
(out/'PLOTS_COMPLETE.json').write_text(json.dumps({'main_panels':2,'seeds':[42,19,888],'candidates':tested.resolution.tolist()},indent=2))
print(tested.to_string(index=False))

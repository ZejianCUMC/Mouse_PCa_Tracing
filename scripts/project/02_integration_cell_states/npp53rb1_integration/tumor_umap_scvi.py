import os,sys,json
from pathlib import Path
import pandas as pd,scanpy as sc,scvi,torch
from lightning.pytorch.callbacks import Callback

out=Path(os.environ['RB1_RUN'])/'tumor_umap_20260921';mode=sys.argv[1]
assert mode in ['run','sample','run_sample']
name='scVI_'+mode;scvi.settings.seed=19;torch.set_num_threads(8)
a=sc.read_h5ad(out/'tumor_hvg_counts.h5ad')
batch='sample_id' if mode=='sample' else 'run_id'
kwargs={'categorical_covariate_keys':['sample_id']} if mode=='run_sample' else {}
scvi.model.SCVI.setup_anndata(a,layer='counts',batch_key=batch,**kwargs)
model=scvi.model.SCVI(a,n_latent=30,n_layers=2,gene_likelihood='nb')
class Progress(Callback):
    def on_train_epoch_end(self,trainer,pl_module):
        if (trainer.current_epoch+1)%10==0:
            print(name,'epoch',trainer.current_epoch+1,flush=True)
model.train(max_epochs=200,early_stopping=True,early_stopping_patience=20,
            accelerator='gpu',devices=1,batch_size=256,enable_progress_bar=False,
            callbacks=[Progress()])
pd.DataFrame(model.get_latent_representation(),index=a.obs_names).to_csv(out/f'latent_{name}.csv')
model.save(out/f'model_{name}',overwrite=False)
for key,value in model.history.items():
    if hasattr(value,'to_csv'):value.to_csv(out/f'{name}_{key}.csv')
(out/f'{name}_complete.json').write_text(json.dumps({'method':name,'batch_key':batch,
 'additional_covariates':kwargs,'likelihood':'nb','max_epochs':200,'batch_size':256,'seed':19,
 'scvi_version':scvi.__version__,'torch_version':torch.__version__,
 'note':'run+sample is joint nuisance conditioning, NOT weighted or hierarchical correction'},indent=2))
print('COMPLETED',name,flush=True)

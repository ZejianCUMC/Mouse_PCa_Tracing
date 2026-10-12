import sys,json,time,platform,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scvi,torch
from lightning.pytorch.callbacks import Callback
out=Path(sys.argv[1]);mode=sys.argv[2];seed=int(sys.argv[3]) if len(sys.argv)>3 else 42
assert mode in ['run','run_sample','sample'] and seed in [42,19,888]
variant=sys.argv[4] if len(sys.argv)>4 else '';assert variant in ['', 'endogenous']
suffix='_endogenous' if variant else ''
dest=out/f'scvi_{mode}_seed{seed}{suffix}';dest.mkdir(exist_ok=True)
assert not (dest/'COMPLETE.json').exists() and not (dest/'model').exists()
scvi.settings.seed=seed;torch.set_num_threads(8)
input_file=out/f'hvg_counts{suffix}.h5ad'
a=sc.read_h5ad(input_file);assert a.shape==(58382,3000)
batch='sample_id' if mode=='sample' else 'run_id'
kw={'categorical_covariate_keys':['sample_id']} if mode=='run_sample' else {}
scvi.model.SCVI.setup_anndata(a,batch_key=batch,**kw)
model=scvi.model.SCVI(a,n_latent=30,n_layers=2,n_hidden=128,dropout_rate=.1,gene_likelihood='nb')
class Progress(Callback):
 def on_train_epoch_end(self,trainer,pl_module):
  if (trainer.current_epoch+1)%5==0:
   print(mode,'seed',seed,'epoch',trainer.current_epoch+1,'elapsed_s',round(time.time()-start),flush=True)
start=time.time()
model.train(max_epochs=200,early_stopping=True,early_stopping_patience=20,accelerator='gpu',devices=1,batch_size=256,enable_progress_bar=False,callbacks=[Progress()])
z=model.get_latent_representation();assert z.shape==(a.n_obs,30) and np.isfinite(z).all()
pd.DataFrame(z,index=a.obs_names).to_csv(dest/'latent_PRIVATE.csv')
model.save(dest/'model',overwrite=False)
for key,value in model.history.items():
 if hasattr(value,'to_csv'):value.to_csv(dest/f'{key}.csv')
(dest/'COMPLETE.json').write_text(json.dumps({'mode':mode,'seed':seed,'batch':batch,'covariates':kw,'variant':variant,'input_file':input_file.name,'cells':a.n_obs,'latent':30,'hidden':128,'layers':2,'dropout':.1,'likelihood':'nb','max_epochs':200,'epochs':len(model.history['elbo_train']),'early_stopping_patience':20,'elapsed_s':time.time()-start,'scvi':scvi.__version__,'torch':torch.__version__,'host':platform.node(),'gpu':torch.cuda.get_device_name(),'input_sha256':hashlib.sha256(input_file.read_bytes()).hexdigest(),'note':'run+sample is joint nuisance conditioning, not weighted correction'},indent=2))
print('COMPLETE',mode,seed,flush=True)

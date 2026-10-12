import os, json
from pathlib import Path
import scanpy as sc, scvi, pandas as pd, torch
out=Path(os.environ['RB1_RUN']); torch.set_num_threads(8)
a=sc.read_h5ad(out/'merged_counts.h5ad'); a=a[:,a.var.highly_variable].copy()
for likelihood in ['nb','zinb']:
 scvi.settings.seed=19
 scvi.model.SCVI.setup_anndata(a,layer='counts',batch_key='batch')
 model=scvi.model.SCVI(a,n_latent=30,n_layers=2,gene_likelihood=likelihood)
 model.train(max_epochs=200,early_stopping=True,accelerator='gpu',devices=1)
 name='scVI_'+likelihood
 pd.DataFrame(model.get_latent_representation(),index=a.obs_names).to_csv(out/f'latent_{name}.csv')
 model.save(out/f'model_{name}',overwrite=False)
 for key,value in model.history.items():
  if hasattr(value,'to_csv'):value.to_csv(out/f'{name}_{key}.csv')
 del model
 torch.cuda.empty_cache()
(out/'scvi_versions.json').write_text(json.dumps({'scvi':scvi.__version__,'torch':torch.__version__,'seed':19,'batch_key':'library, not treatment or mouse'},indent=2))

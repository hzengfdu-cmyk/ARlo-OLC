# scTour trajectory analysis

from pathlib import Path
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import sctour as sct
from scipy.io import mmread
ROOT = Path(os.environ.get("OC_PROJECT", Path.cwd()))
OUT = Path(os.environ.get("OC_RESULTS", ROOT / "results")) / "08_scTour"
meta = pd.read_csv(OUT / "metadata.csv", index_col="cell_id")
for column in meta.select_dtypes(include=["object"]).columns:
    value = meta[column]
    meta[column] = pd.Categorical(value.astype(str).where(value.notna()))
genes = pd.read_csv(OUT / "genes.csv")["gene"].astype(str)
counts = mmread(OUT / "counts.mtx").T.tocsr().astype(np.float32)
adata = ad.AnnData(X=counts,obs=meta,var=pd.DataFrame(index=genes))
sc.pp.filter_genes(adata,min_cells=20)
sc.pp.calculate_qc_metrics(adata,percent_top=None,log1p=False,inplace=True)
trainer = sct.train.Trainer(adata,loss_mode="nb",alpha_recon_lec=.5,
                            alpha_recon_lode=.5,random_state=123)
trainer.train()
adata.obs["ptime"] = trainer.get_time()
mix_zs, zs, pred_zs = trainer.get_latentsp(alpha_z=.5,alpha_predz=.5)
adata.obsm["X_TNODE"] = mix_zs
adata.obsm["X_VF"] = trainer.get_vector_field(adata.obs.ptime.to_numpy(),mix_zs)
if (OUT / "tsne.csv").exists():
    coords = pd.read_csv(OUT / "tsne.csv",index_col=0)
    adata.obsm["X_tsne"] = coords.loc[adata.obs_names].to_numpy()
trainer.save_model(str(OUT),"scTour_model")
adata.obs.to_csv(OUT / "pseudotime_source.csv")
adata.write_h5ad(OUT / "sctour_result.h5ad")

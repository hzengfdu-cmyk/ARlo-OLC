# Per-library Scrublet analysis

from pathlib import Path
import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
from scipy import sparse
ROOT = Path(os.environ.get("OC_PROJECT", Path.cwd()))
OUT = Path(os.environ.get("OC_RESULTS", ROOT / "results")) / "01_scrublet"
OUT.mkdir(parents=True, exist_ok=True)
manifest = pd.read_csv(ROOT / "data/scrublet_libraries.csv")
required = {"library_id", "h5ad", "expected_doublet_rate"}
rows = []
for row in manifest.itertuples(index=False):
    adata = sc.read_h5ad(row.h5ad)
    counts = adata.X
    values = counts.data if sparse.issparse(counts) else np.asarray(counts)
    scrub = scr.Scrublet(counts, expected_doublet_rate=row.expected_doublet_rate,
                        random_state=123)
    scores, predicted = scrub.scrub_doublets(n_prin_comps=30)
    frame = pd.DataFrame({"cell_id": adata.obs_names, "library_id": row.library_id,
                          "doublet_score": scores, "predicted_doublet": predicted})
    rows.append(frame)
    pd.DataFrame({"simulated_score": scrub.doublet_scores_sim_}).to_csv(
        OUT / f"{row.library_id}_simulated_scores.csv", index=False)
    (OUT / f"{row.library_id}_parameters.json").write_text(json.dumps({
        "expected_doublet_rate": row.expected_doublet_rate,
        "threshold": float(scrub.threshold_), "seed": 123}, indent=2))
calls = pd.concat(rows, ignore_index=True)
calls.to_csv(OUT / "doublet_calls.csv", index=False)

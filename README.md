# Osteoclast lineage cells and androgen-response analysis code

This repository provides the main computational scripts and parameters described in the study Methods.

## Modules

- Human and mouse single-cell QC, Harmony integration, and clustering
- AUCell androgen-response scoring and ARhi/ARlo OLC states
- Mfuzz expression modules and enrichment analysis
- Monocle2 and scTour pseudotime analysis
- CellChat communication and ligand-receptor filtering
- Cytokine-array and single-cell integration
- ARlo signature scoring and survival analysis

Scripts are stored in `code/`, and the AUCell gene sets are stored in `gene_sets/`. Edit `code/00_settings.R` before running the analyses.

CUT&Tag, ChIP-seq, ImageJ, and Prism processing are outside this code package.

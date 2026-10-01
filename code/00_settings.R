# Shared paths and Methods parameters
PROJECT <- normalizePath(getwd())
INPUT <- file.path(PROJECT, "data")
RESULTS <- file.path(PROJECT, "results")
CONFIG <- file.path(PROJECT, "gene_sets")

SEED <- 123L
HUMAN_RAW <- file.path(INPUT, "human_raw.rds")
HUMAN_ANNOTATED <- file.path(INPUT, "human_annotated.rds")
HUMAN_SCRUBLET_FILE <- file.path(INPUT, "adata_scrublet.csv")
HUMAN_GENE_REMOVE_FILE <- file.path(INPUT, "genes_to_remove.txt")
HUMAN_PSEUDOGENE_FILE <- file.path(INPUT, "pseudogene.txt")
MOUSE_RAW <- file.path(INPUT, "mouse_raw.rds")

HUMAN_CLUSTER_RES <- 0.5
HUMAN_GRAPH_DIMS <- 1:30
MOUSE_GRAPH_DIMS <- 1:20
CELLTYPE_COLUMN <- "celltype1"
SAMPLE_COLUMN <- "orig.ident"
AUCMAX_FRACTION <- 0.05
MFUZZ_MARROW_TYPES <- c("benign_marrow", "distal_marrow", "involved_marrow")
MFUZZ_SEED <- 5201314L
MONOCLE_ROOT_STATE <- NA_integer_
CELLCHAT_MIN_CELLS <- 10L
CELLCHAT_POPULATION_SIZE <- TRUE
GENEFU_LIB <- Sys.getenv("OC_GENEFU_LIB", unset = "")
DE_ENRICH_LOGFC <- 0.2
DE_ENRICH_PADJ <- 0.05

dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

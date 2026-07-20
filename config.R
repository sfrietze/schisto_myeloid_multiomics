REPO_DIR <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

DATA_DIR <- file.path(REPO_DIR, "data")
RESULTS_DIR <- file.path(REPO_DIR, "results")
FIGURES_DIR <- file.path(REPO_DIR, "figures")
LOG_DIR <- file.path(REPO_DIR, "logs")

ATACSEQ_DATA_DIR <- file.path(DATA_DIR, "atacseq")
ATACSEQ_RESULTS_DIR <- file.path(RESULTS_DIR, "atacseq")
ATACSEQ_FIGURES_DIR <- file.path(FIGURES_DIR, "atacseq")

BIGWIG_DIR <- file.path(
  ATACSEQ_DATA_DIR,
  "bigwig"
)

TOBIAS_DIR <- file.path(
  ATACSEQ_DATA_DIR,
  "tobias"
)

FIGURE1_RESULTS_DIR <- file.path(
  ATACSEQ_RESULTS_DIR,
  "figure1"
)

FIGURE7A_RESULTS_DIR <- file.path(
  ATACSEQ_RESULTS_DIR,
  "figure7A"
)

FIGURE7C_RESULTS_DIR <- file.path(
  ATACSEQ_RESULTS_DIR,
  "figure7C"
)

FIGURE7D_RESULTS_DIR <- file.path(
  ATACSEQ_RESULTS_DIR,
  "figure7D"
)

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(
  FIGURE1_RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE7A_RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE7C_RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE7D_RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

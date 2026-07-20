# Repository root
REPO_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Original analysis directory (one level above the repository)
PROJECT_DIR <- normalizePath(file.path(REPO_DIR, ".."), winslash = "/", mustWork = TRUE)

# Input data
DAR_DIR        <- file.path(PROJECT_DIR, "DAR")
FEMALE_INF_DIR <- file.path(PROJECT_DIR, "female_inf")
FEMALE_SEA_DIR <- file.path(PROJECT_DIR, "female_sea")
MALE_INF_DIR   <- file.path(PROJECT_DIR, "male_inf")

# Output directories
RESULTS_DIR <- file.path(REPO_DIR, "results")
LOG_DIR     <- file.path(REPO_DIR, "logs")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "atacseq", "figure1"),
           recursive = TRUE,
           showWarnings = FALSE)

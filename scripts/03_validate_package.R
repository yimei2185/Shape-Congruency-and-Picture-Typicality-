# Validate the packaged data and saved model objects without refitting models.

required <- c("readr", "dplyr", "lme4")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

library(readr)
library(dplyr)
library(lme4)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript so the package root can be located.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
package_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")

prt_srt <- read_csv(file.path(package_dir, "data", "prt_srt_analysis.csv"), show_col_types = FALSE)
accuracy <- read_csv(file.path(package_dir, "data", "accuracy_analysis.csv"), show_col_types = FALSE)
l2 <- read_csv(file.path(package_dir, "data", "l2_proficiency_analysis.csv"), show_col_types = FALSE)
l1 <- read_csv(file.path(package_dir, "data", "l1_only_prt_analysis.csv"), show_col_types = FALSE)

stopifnot(
  nrow(prt_srt) == 3434,
  nrow(accuracy) == 3492,
  nrow(l2) == 1572,
  nrow(l1) == 1862,
  sum(accuracy$accuracy == 0) == 58,
  !anyNA(prt_srt[c("participant", "log_picture_rt", "log_sentence_rt")]),
  !anyNA(accuracy[c("participant", "accuracy")]),
  !anyNA(l2[c("participant", "proficiency_z", "log_picture_rt")]),
  !anyNA(l1[c("participant", "log_picture_rt")])
)

expected <- c(
  picture_rt = 9L, accuracy = 9L, sentence_rt = 6L,
  l2_proficiency = 9L, l1_only_prt = 9L
)
expected_n <- c(
  picture_rt = 3434L, accuracy = 3492L, sentence_rt = 3434L,
  l2_proficiency = 1572L, l1_only_prt = 1862L
)
lines <- c(
  "OSF mixed-model package validation: PASS",
  paste0("Validated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Data:"
)

lines <- c(
  lines,
  sprintf("- prt_srt_analysis.csv: %d rows", nrow(prt_srt)),
  sprintf("- accuracy_analysis.csv: %d rows (%d errors)", nrow(accuracy), sum(accuracy$accuracy == 0)),
  sprintf("- l2_proficiency_analysis.csv: %d rows", nrow(l2)),
  sprintf("- l1_only_prt_analysis.csv: %d rows", nrow(l1)),
  "",
  "Saved models:"
)

for (analysis in names(expected)) {
  candidates <- readRDS(file.path(package_dir, "models", paste0(analysis, "_all_candidates.rds")))
  final <- readRDS(file.path(package_dir, "models", paste0(analysis, "_selected_model.rds")))
  stopifnot(
    length(candidates) == expected[[analysis]],
    inherits(final, "merMod"),
    nobs(final) == expected_n[[analysis]],
    !isSingular(final, tol = 1e-4)
  )
  lines <- c(
    lines,
    sprintf(
      "- %s: %d candidates; selected model %d observations; non-singular",
      analysis, length(candidates), nobs(final)
    )
  )
}

writeLines(lines, file.path(package_dir, "VALIDATION.txt"), useBytes = TRUE)
cat(paste(lines, collapse = "\n"), "\n")

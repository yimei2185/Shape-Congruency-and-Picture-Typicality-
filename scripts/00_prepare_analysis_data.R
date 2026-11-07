# Prepare de-identified, analysis-ready data for the OSF mixed-model package.
# This script reads the frozen preprocessing outputs and does not alter them.

required <- c("readr", "dplyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

library(readr)
library(dplyr)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript so paths can be located.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
package_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
project_default <- normalizePath(file.path(package_dir, "../.."), winslash = "/", mustWork = FALSE)
project_dir <- normalizePath(
  Sys.getenv("TYPICALITY_PROJECT_DIR", unset = project_default),
  winslash = "/", mustWork = FALSE
)
data_dir <- file.path(package_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

write_osf_csv <- function(x, name) {
  write_excel_csv(x, file.path(data_dir, name), na = "")
}

# -----------------------------------------------------------------------------
# PRT and SRT data: already frozen to correct retained critical trials.
# -----------------------------------------------------------------------------

prt_srt_source <- file.path(
  project_dir,
  "outputs/corrected_mixed_models/00_picture_rt_and_srt_model_data.csv"
)

prt_srt <- read_csv(prt_srt_source, show_col_types = FALSE) |>
  transmute(
    participant = as.character(participant),
    participant_num = as.integer(participant_num),
    language = as.character(language),
    language_c = if_else(language == "L1", -0.5, 0.5),
    condition = as.character(condition),
    condition_c = if_else(condition == "congruent", -0.5, 0.5),
    target_item = as.character(target_item),
    picture_file = as.character(picture_file),
    sentence_id = as.character(sentence_id),
    picture_typicality = as.numeric(picture_typicality),
    picture_typicality_c = as.numeric(picture_typicality_c),
    sentence_implied_typicality = as.numeric(sentence_implied_typicality),
    sentence_implied_typicality_c = as.numeric(sentence_implied_typicality_c),
    sentence_length = as.numeric(sentence_length),
    sentence_length_c = as.numeric(sentence_length_c),
    trial_order_c = as.numeric(trial_order_c),
    picture_rt_ms = as.numeric(picture_rt_ms),
    sentence_rt_ms = as.numeric(sentence_rt_ms),
    log_picture_rt = as.numeric(log_picture_rt),
    log_sentence_rt = as.numeric(log_sentence_rt)
  )

stopifnot(
  nrow(prt_srt) == 3434,
  n_distinct(prt_srt$participant) == 157,
  n_distinct(prt_srt$target_item) == 24,
  n_distinct(prt_srt$picture_file) == 48,
  n_distinct(prt_srt$sentence_id) == 48
)

write_osf_csv(prt_srt, "prt_srt_analysis.csv")

# -----------------------------------------------------------------------------
# Accuracy data: critical trials retained by the same picture-RT screen, before
# incorrect trials were removed.
# -----------------------------------------------------------------------------

pre_accuracy_file <- file.path(
  project_dir,
  "outputs/corrected_preprocessing/R_outputs/06_critical_trials_after_target_item_filter.csv"
)
rt_flags_file <- file.path(
  project_dir,
  "outputs/corrected_preprocessing/R_outputs/07_picture_rt_flags_all_formal_trials.csv"
)

pre_accuracy <- read_csv(pre_accuracy_file, show_col_types = FALSE)
rt_flags <- read_csv(rt_flags_file, show_col_types = FALSE) |>
  select(participant, trial_total, picture_rt_retained)

picture_reference_mean <- prt_srt |>
  distinct(picture_file, picture_typicality) |>
  summarise(value = mean(picture_typicality), .groups = "drop") |>
  pull(value)

accuracy <- pre_accuracy |>
  left_join(
    rt_flags,
    by = c("participant", "trial_total"),
    relationship = "one-to-one"
  ) |>
  filter(picture_rt_retained) |>
  transmute(
    participant = as.character(participant),
    participant_num = as.integer(participant_num),
    language = as.character(language),
    language_c = if_else(language == "L1", -0.5, 0.5),
    condition = as.character(condition),
    condition_c = if_else(condition == "congruent", -0.5, 0.5),
    target_item = as.character(target_item),
    picture_file = as.character(picture_file),
    accuracy = as.integer(accuracy),
    picture_typicality = as.numeric(picture_typicality),
    picture_typicality_c = picture_typicality - picture_reference_mean
  )

stopifnot(
  nrow(accuracy) == 3492,
  sum(accuracy$accuracy == 1) == 3434,
  sum(accuracy$accuracy == 0) == 58
)

write_osf_csv(accuracy, "accuracy_analysis.csv")

# -----------------------------------------------------------------------------
# L2-only proficiency data: frozen output from the R three-way preparation.
# -----------------------------------------------------------------------------

l2_source <- file.path(
  project_dir,
  "outputs/l2_proficiency_congruency_r_threeway/00_l2_threeway_model_data.csv"
)

l2 <- read_csv(l2_source, show_col_types = FALSE) |>
  transmute(
    participant = as.character(participant),
    participant_num = as.integer(participant_num),
    condition = as.character(condition),
    condition_c = as.numeric(cong_c),
    target_item = as.character(target_item),
    picture_file = as.character(picture_file),
    sentence_id = as.character(sentence_id),
    vocabulary_score = as.numeric(vocabularysize_score),
    communication_score = as.numeric(communicationtest_score),
    proficiency_z = as.numeric(proficiency_z),
    picture_typicality = as.numeric(picture_typicality),
    picture_typicality_c = as.numeric(picture_typicality_c_r),
    trial_order_z = as.numeric(trial_order_z),
    picture_rt_ms = as.numeric(picture_rt_ms),
    log_picture_rt = as.numeric(log_picture_rt)
  )

stopifnot(
  nrow(l2) == 1572,
  n_distinct(l2$participant) == 77,
  n_distinct(l2$target_item) == 24,
  n_distinct(l2$picture_file) == 48
)

write_osf_csv(l2, "l2_proficiency_analysis.csv")

# -----------------------------------------------------------------------------
# Data dictionary and provenance.
# -----------------------------------------------------------------------------

dictionary <- tribble(
  ~file, ~variable, ~description,
  "all", "participant", "Pseudonymous participant identifier",
  "all", "participant_num", "Numeric pseudonymous participant identifier",
  "all", "condition", "Sentence-picture condition: congruent or incongruent",
  "all", "condition_c", "Effect coding: congruent = -0.5, incongruent = 0.5",
  "prt_srt_analysis.csv; accuracy_analysis.csv", "language", "Language group: L1 or L2",
  "prt_srt_analysis.csv; accuracy_analysis.csv", "language_c", "Effect coding: L1 = -0.5, L2 = 0.5",
  "all", "target_item", "Target noun/item identifier",
  "all", "picture_file", "Picture stimulus identifier",
  "prt_srt_analysis.csv; l2_proficiency_analysis.csv", "sentence_id", "Sentence stimulus identifier",
  "all", "picture_typicality", "Continuous picture typicality rating",
  "all", "picture_typicality_c", "Picture typicality centered on the 48-picture mean",
  "prt_srt_analysis.csv", "sentence_implied_typicality", "Typicality of the sentence's congruent picture",
  "prt_srt_analysis.csv", "sentence_implied_typicality_c", "Centered sentence-implied typicality",
  "prt_srt_analysis.csv", "sentence_length_c", "Mean-centered sentence length",
  "prt_srt_analysis.csv", "trial_order_c", "Mean-centered trial order",
  "prt_srt_analysis.csv", "log_picture_rt", "Natural log of picture RT in milliseconds",
  "prt_srt_analysis.csv", "log_sentence_rt", "Natural log of sentence RT in milliseconds",
  "accuracy_analysis.csv", "accuracy", "Picture verification accuracy: 0 or 1",
  "l2_proficiency_analysis.csv", "proficiency_z", "Standardized composite L2 proficiency",
  "l2_proficiency_analysis.csv", "vocabulary_score", "Chinese vocabulary test score",
  "l2_proficiency_analysis.csv", "communication_score", "Communicative ability score"
)

write_osf_csv(dictionary, "data_dictionary.csv")

provenance <- tibble(
  packaged_file = c(
    "prt_srt_analysis.csv",
    "accuracy_analysis.csv",
    "l2_proficiency_analysis.csv"
  ),
  source = c(
    "outputs/corrected_mixed_models/00_picture_rt_and_srt_model_data.csv",
    paste(
      "outputs/corrected_preprocessing/R_outputs/06_critical_trials_after_target_item_filter.csv",
      "+ outputs/corrected_preprocessing/R_outputs/07_picture_rt_flags_all_formal_trials.csv"
    ),
    "outputs/l2_proficiency_congruency_r_threeway/00_l2_threeway_model_data.csv"
  ),
  rows = c(nrow(prt_srt), nrow(accuracy), nrow(l2))
)

write_osf_csv(provenance, "data_provenance.csv")

message("OSF analysis data prepared in: ", data_dir)

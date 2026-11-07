# Recompute Nakagawa R2 for every selected model and refresh overview tables.
# For LMMs, MuMIn values are independently verified from mean(diag(Z G Z')).

required <- c("readr", "dplyr", "lme4", "MuMIn", "performance")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lme4)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript so the package root can be located.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
package_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
models_dir <- file.path(package_dir, "models")
results_dir <- file.path(package_dir, "results")

selected_names <- c(
  picture_rt = "M6_drop_item_condition",
  accuracy = "M9_random_intercepts",
  sentence_rt = "M6_random_intercepts",
  l2_proficiency = "M8_drop_participant_typicality",
  l1_only_prt = "M9_random_intercepts"
)

convergence_message <- function(model) {
  msg <- model@optinfo$conv$lme4$messages
  if (is.null(msg)) "" else paste(msg, collapse = "; ")
}

direct_lmm_r2 <- function(model) {
  X <- getME(model, "X")
  fixed_variance <- var(as.numeric(X %*% fixef(model)))
  Z <- t(getME(model, "Zt"))
  Lambda <- getME(model, "Lambda")
  sigma2 <- sigma(model)^2
  observation_random_variance <- sigma2 * rowSums(as.matrix(Z %*% Lambda)^2)
  random_variance <- mean(observation_random_variance)
  total <- fixed_variance + random_variance + sigma2
  tibble(
    method = "direct mean diagonal of Z G Z-transpose",
    fixed_variance = fixed_variance,
    mean_random_variance = random_variance,
    residual_variance = sigma2,
    marginal_r2 = fixed_variance / total,
    conditional_r2 = (fixed_variance + random_variance) / total
  )
}

rows <- vector("list", length(selected_names))
names(rows) <- names(selected_names)

for (analysis in names(selected_names)) {
  model <- readRDS(file.path(models_dir, paste0(analysis, "_selected_model.rds")))
  if (inherits(model, "lmerMod")) {
    mumin <- unname(MuMIn::r.squaredGLMM(model)[1, c("R2m", "R2c")])
    direct <- direct_lmm_r2(model)
    stopifnot(
      isTRUE(all.equal(mumin[[1]], direct$marginal_r2[[1]], tolerance = 1e-10)),
      isTRUE(all.equal(mumin[[2]], direct$conditional_r2[[1]], tolerance = 1e-10))
    )
    r2 <- mumin
    write_excel_csv(
      direct,
      file.path(results_dir, analysis, "selected_model_r2_verified.csv")
    )
  } else {
    performance_r2 <- performance::r2(model)
    r2 <- c(
      unname(performance_r2$R2_marginal),
      unname(performance_r2$R2_conditional)
    )
  }

  rows[[analysis]] <- tibble(
    analysis = analysis,
    selected_candidate = unname(selected_names[[analysis]]),
    final_formula = paste(deparse(formula(model)), collapse = " "),
    observations = nobs(model),
    singular = isSingular(model, tol = 1e-4),
    convergence_message = convergence_message(model),
    marginal_r2 = r2[[1]],
    conditional_r2 = r2[[2]]
  )

  write_excel_csv(
    rows[[analysis]],
    file.path(results_dir, analysis, "selected_model_overview.csv"),
    na = ""
  )
}

write_excel_csv(bind_rows(rows), file.path(results_dir, "model_overview.csv"), na = "")
message("Refreshed R2 and overview tables for ", length(rows), " selected models.")

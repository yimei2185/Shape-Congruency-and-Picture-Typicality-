# Han-and-Chen-style candidate random-effects comparisons for the current data.
# All candidate models within an analysis use the same observations and fixed
# effects. Random structures are selected without reference to fixed-effect p
# values. LMM candidates use ML; selected LMMs are refitted with REML.

required <- c(
  "readr", "dplyr", "tidyr", "lme4", "lmerTest", "broom.mixed",
  "performance", "emmeans", "MuMIn"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

library(readr)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(performance)
library(emmeans)
library(MuMIn)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript so the package root can be located.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
package_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
data_dir <- file.path(package_dir, "data")
results_dir <- file.path(package_dir, "results")
models_dir <- file.path(package_dir, "models")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

options(contrasts = c("contr.sum", "contr.poly"))
boundary_tol <- 1e-7
singular_tol <- 1e-4

lmm_control <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 300000),
  check.conv.singular = "ignore"
)
glmm_control <- glmerControl(
  optimizer = "bobyqa",
  # Accuracy has only 58 errors. A finite evaluation budget prevents a clearly
  # over-parameterized candidate from monopolizing the reproducible workflow.
  optCtrl = list(maxfun = 50000),
  check.conv.singular = "ignore"
)

write_result <- function(x, path) {
  write_excel_csv(x, path, na = "")
}

fit_error <- function(message) structure(list(error = message), class = "fit_error")
is_fit_error <- function(x) inherits(x, "fit_error")

safe_lmer <- function(formula, data, reml = FALSE) {
  tryCatch(
    lmer(formula, data = data, REML = reml, control = lmm_control),
    error = function(e) fit_error(conditionMessage(e))
  )
}

safe_glmer <- function(formula, data, nAGQ = 0) {
  tryCatch(
    glmer(
      formula,
      data = data,
      family = binomial(link = "logit"),
      control = glmm_control,
      nAGQ = nAGQ
    ),
    error = function(e) fit_error(conditionMessage(e))
  )
}

convergence_message <- function(model) {
  if (is_fit_error(model)) return(model$error)
  msg <- model@optinfo$conv$lme4$messages
  if (is.null(msg)) "" else paste(msg, collapse = "; ")
}

variance_diagnostics <- function(model) {
  if (is_fit_error(model)) {
    return(list(min_variance = NA_real_, boundary_terms = ""))
  }
  vc <- as.data.frame(VarCorr(model))
  variance_rows <- vc |> filter(is.na(var2))
  boundary <- variance_rows |> filter(vcov < boundary_tol)
  list(
    min_variance = if (nrow(variance_rows)) min(variance_rows$vcov) else NA_real_,
    boundary_terms = if (nrow(boundary)) {
      paste(paste0(boundary$grp, ":", boundary$var1), collapse = "; ")
    } else ""
  )
}

diagnostic_row <- function(fit, name, formula_text, order) {
  if (is_fit_error(fit)) {
    return(tibble(
      order = order, model = name, formula = formula_text, fitted = FALSE,
      converged = FALSE, singular = NA, npar = NA_integer_, nobs = NA_integer_,
      logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, min_random_variance = NA_real_,
      boundary_terms_below_1e_7 = "", convergence_message = fit$error
    ))
  }
  vd <- variance_diagnostics(fit)
  conv <- convergence_message(fit)
  singular_value <- isSingular(fit, tol = singular_tol)
  log_likelihood <- logLik(fit)
  tibble(
    order = order,
    model = name,
    formula = formula_text,
    fitted = TRUE,
    converged = identical(conv, ""),
    singular = singular_value,
    npar = attr(log_likelihood, "df"),
    nobs = nobs(fit),
    logLik = as.numeric(log_likelihood),
    AIC = AIC(fit),
    BIC = BIC(fit),
    min_random_variance = vd$min_variance,
    boundary_terms_below_1e_7 = vd$boundary_terms,
    convergence_message = conv
  )
}

adjacent_lrt <- function(models, diagnostics) {
  out <- vector("list", max(0, nrow(diagnostics) - 1))
  if (!length(out)) return(tibble())
  for (i in seq_len(nrow(diagnostics) - 1)) {
    complex_name <- diagnostics$model[[i]]
    simple_name <- diagnostics$model[[i + 1]]
    complex <- models[[complex_name]]
    simple <- models[[simple_name]]
    if (is_fit_error(complex) || is_fit_error(simple)) {
      out[[i]] <- tibble(
        complex_model = complex_name, simpler_model = simple_name,
        chisq = NA_real_, df = NA_real_, p = NA_real_, note = "model fit error"
      )
      next
    }
    test <- tryCatch(
      suppressMessages(anova(simple, complex, refit = FALSE)),
      error = function(e) e
    )
    if (inherits(test, "error")) {
      out[[i]] <- tibble(
        complex_model = complex_name, simpler_model = simple_name,
        chisq = NA_real_, df = NA_real_, p = NA_real_, note = conditionMessage(test)
      )
    } else {
      out[[i]] <- tibble(
        complex_model = complex_name,
        simpler_model = simple_name,
        chisq = test$Chisq[[2]],
        df = test$`Chi Df`[[2]],
        p = test$`Pr(>Chisq)`[[2]],
        note = "ordinary boundary LRT; interpret as approximate"
      )
    }
  }
  bind_rows(out)
}

select_supported <- function(models, diagnostics) {
  supported <- diagnostics |>
    filter(fitted, converged, !singular, boundary_terms_below_1e_7 == "", !is.na(AIC))
  if (!nrow(supported)) {
    supported <- diagnostics |>
      filter(fitted, converged, !singular, !is.na(AIC))
  }
  if (!nrow(supported)) stop("No converged non-singular candidate model.")
  best_aic <- min(supported$AIC)
  selected_row <- supported |>
    filter(AIC <= best_aic + 2) |>
    arrange(npar, AIC, order) |>
    slice(1)
  list(name = selected_row$model[[1]], model = models[[selected_row$model[[1]]]])
}

fit_candidate_set <- function(
    analysis, data, response, fixed, random_candidates,
    model_type = c("lmm", "glmm")) {
  model_type <- match.arg(model_type)
  analysis_result_dir <- file.path(results_dir, analysis)
  dir.create(analysis_result_dir, recursive = TRUE, showWarnings = FALSE)

  models <- list()
  rows <- list()
  for (i in seq_along(random_candidates)) {
    name <- names(random_candidates)[[i]]
    formula_text <- paste(response, "~", fixed, "+", random_candidates[[i]])
    message("[", analysis, "] fitting ", name)
    model <- if (model_type == "lmm") {
      safe_lmer(as.formula(formula_text), data, reml = FALSE)
    } else {
      # nAGQ = 0 is used only to screen the complete random-structure set.
      # The selected GLMM is refitted below with the standard nAGQ = 1.
      safe_glmer(as.formula(formula_text), data, nAGQ = 0)
    }
    models[[name]] <- model
    rows[[i]] <- diagnostic_row(model, name, formula_text, i)
    # Checkpoint after every candidate so completed fits survive a later
    # optimizer failure or an interrupted long-running candidate.
    partial_diagnostics <- bind_rows(rows)
    write_result(
      partial_diagnostics,
      file.path(analysis_result_dir, "candidate_model_comparison_partial.csv")
    )
    saveRDS(
      models,
      file.path(models_dir, paste0(analysis, "_all_candidates_partial.rds"))
    )
  }

  diagnostics <- bind_rows(rows) |>
    mutate(delta_AIC = AIC - min(AIC[fitted & converged], na.rm = TRUE))
  lrt <- adjacent_lrt(models, diagnostics)
  selected <- select_supported(models, diagnostics)

  write_result(
    diagnostics,
    file.path(analysis_result_dir, "candidate_model_comparison.csv")
  )
  write_result(
    lrt,
    file.path(analysis_result_dir, "adjacent_likelihood_ratio_tests.csv")
  )
  saveRDS(models, file.path(models_dir, paste0(analysis, "_all_candidates.rds")))

  final_model <- selected$model
  if (model_type == "lmm") {
    final_model <- safe_lmer(formula(selected$model), data, reml = TRUE)
    if (is_fit_error(final_model)) stop("REML refit failed for ", analysis)
  } else {
    final_model <- safe_glmer(formula(selected$model), data, nAGQ = 1)
    if (is_fit_error(final_model)) stop("nAGQ=1 refit failed for ", analysis)
  }
  saveRDS(final_model, file.path(models_dir, paste0(analysis, "_selected_model.rds")))

  fixed_table <- tidy(
    final_model,
    effects = "fixed",
    conf.int = TRUE,
    conf.method = if (model_type == "glmm") "Wald" else "Wald"
  )
  random_table <- as.data.frame(VarCorr(final_model))
  # MuMIn correctly includes all variance components in models expressed with
  # independent (`||`) random slopes. performance/insight 0.15.2/1.4.2 can
  # omit the corresponding group intercept variance in this configuration.
  r2 <- if (model_type == "lmm") {
    unname(MuMIn::r.squaredGLMM(final_model)[1, c("R2m", "R2c")])
  } else {
    performance_r2 <- performance::r2(final_model)
    c(
      unname(performance_r2$R2_marginal),
      unname(performance_r2$R2_conditional)
    )
  }
  overview <- tibble(
    analysis = analysis,
    selected_candidate = selected$name,
    final_formula = paste(deparse(formula(final_model)), collapse = " "),
    observations = nobs(final_model),
    singular = isSingular(final_model, tol = singular_tol),
    convergence_message = convergence_message(final_model),
    marginal_r2 = r2[[1]],
    conditional_r2 = r2[[2]]
  )

  write_result(fixed_table, file.path(analysis_result_dir, "selected_fixed_effects.csv"))
  write_result(random_table, file.path(analysis_result_dir, "selected_random_effects.csv"))
  write_result(overview, file.path(analysis_result_dir, "selected_model_overview.csv"))
  capture.output(
    summary(final_model),
    file = file.path(analysis_result_dir, "selected_model_summary.txt")
  )

  list(
    models = models,
    diagnostics = diagnostics,
    lrt = lrt,
    selected_name = selected$name,
    final = final_model,
    overview = overview
  )
}

factorize <- function(data, columns) {
  for (column in columns) data[[column]] <- factor(data[[column]])
  data
}

prt_srt <- read_csv(file.path(data_dir, "prt_srt_analysis.csv"), show_col_types = FALSE) |>
  factorize(c("participant", "target_item", "picture_file", "sentence_id"))
accuracy <- read_csv(file.path(data_dir, "accuracy_analysis.csv"), show_col_types = FALSE) |>
  factorize(c("participant", "target_item", "picture_file"))
l2 <- read_csv(file.path(data_dir, "l2_proficiency_analysis.csv"), show_col_types = FALSE) |>
  factorize(c("participant", "target_item", "picture_file", "sentence_id"))

stopifnot(nrow(prt_srt) == 3434, nrow(accuracy) == 3492, nrow(l2) == 1572)

# -----------------------------------------------------------------------------
# Picture RT. Candidate path begins with article-style correlated main-effect
# slopes, then removes correlations and unsupported slopes in a nested sequence.
# -----------------------------------------------------------------------------

prt_random <- c(
  M1_correlated_full = paste0(
    "(1 + condition_c + picture_typicality_c | participant) + ",
    "(1 + language_c + condition_c + picture_typicality_c | target_item) + ",
    "(1 + language_c + condition_c | picture_file)"
  ),
  M2_uncorrelated_full = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + language_c + condition_c + picture_typicality_c || target_item) + ",
    "(1 + language_c + condition_c || picture_file)"
  ),
  M3_drop_picture_condition = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + language_c + condition_c + picture_typicality_c || target_item) + ",
    "(1 + language_c || picture_file)"
  ),
  M4_picture_intercept = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + language_c + condition_c + picture_typicality_c || target_item) + ",
    "(1 | picture_file)"
  ),
  M5_drop_item_typicality = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + language_c + condition_c || target_item) + (1 | picture_file)"
  ),
  M6_drop_item_condition = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + language_c || target_item) + (1 | picture_file)"
  ),
  M7_item_intercept = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 | target_item) + (1 | picture_file)"
  ),
  M8_drop_participant_typicality = paste0(
    "(1 + condition_c || participant) + (1 | target_item) + (1 | picture_file)"
  ),
  M9_random_intercepts =
    "(1 | participant) + (1 | target_item) + (1 | picture_file)"
)

prt_fit <- fit_candidate_set(
  "picture_rt", prt_srt, "log_picture_rt",
  "language_c * condition_c * picture_typicality_c",
  prt_random, "lmm"
)

# -----------------------------------------------------------------------------
# Accuracy GLMM, using the same candidate path as picture RT.
# -----------------------------------------------------------------------------

accuracy_fit <- fit_candidate_set(
  "accuracy", accuracy, "accuracy",
  "language_c * condition_c * picture_typicality_c",
  prt_random, "glmm"
)

# -----------------------------------------------------------------------------
# Sentence RT. Sentence-implied typicality is constant within sentence, so the
# sentence unit receives a language slope but not a typicality slope.
# -----------------------------------------------------------------------------

srt_random <- c(
  M1_correlated_full = paste0(
    "(1 + sentence_implied_typicality_c | participant) + ",
    "(1 + language_c + sentence_implied_typicality_c | target_item) + ",
    "(1 + language_c | sentence_id)"
  ),
  M2_uncorrelated_full = paste0(
    "(1 + sentence_implied_typicality_c || participant) + ",
    "(1 + language_c + sentence_implied_typicality_c || target_item) + ",
    "(1 + language_c || sentence_id)"
  ),
  M3_sentence_intercept = paste0(
    "(1 + sentence_implied_typicality_c || participant) + ",
    "(1 + language_c + sentence_implied_typicality_c || target_item) + ",
    "(1 | sentence_id)"
  ),
  M4_drop_item_typicality = paste0(
    "(1 + sentence_implied_typicality_c || participant) + ",
    "(1 + language_c || target_item) + (1 | sentence_id)"
  ),
  M5_item_intercept = paste0(
    "(1 + sentence_implied_typicality_c || participant) + ",
    "(1 | target_item) + (1 | sentence_id)"
  ),
  M6_random_intercepts =
    "(1 | participant) + (1 | target_item) + (1 | sentence_id)"
)

srt_fit <- fit_candidate_set(
  "sentence_rt", prt_srt, "log_sentence_rt",
  paste0(
    "language_c * sentence_implied_typicality_c + ",
    "sentence_length_c + trial_order_c"
  ),
  srt_random, "lmm"
)

# -----------------------------------------------------------------------------
# L2-only proficiency model. Proficiency is constant within participant and is
# therefore not included as a by-participant random slope.
# -----------------------------------------------------------------------------

l2_random <- c(
  M1_correlated_full = paste0(
    "(1 + condition_c + picture_typicality_c | participant) + ",
    "(1 + proficiency_z + condition_c + picture_typicality_c | target_item) + ",
    "(1 + proficiency_z + condition_c | picture_file)"
  ),
  M2_uncorrelated_full = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + proficiency_z + condition_c + picture_typicality_c || target_item) + ",
    "(1 + proficiency_z + condition_c || picture_file)"
  ),
  M3_drop_picture_condition = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + proficiency_z + condition_c + picture_typicality_c || target_item) + ",
    "(1 + proficiency_z || picture_file)"
  ),
  M4_picture_intercept = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + proficiency_z + condition_c + picture_typicality_c || target_item) + ",
    "(1 | picture_file)"
  ),
  M5_drop_item_typicality = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + proficiency_z + condition_c || target_item) + (1 | picture_file)"
  ),
  M6_drop_item_condition = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + proficiency_z || target_item) + (1 | picture_file)"
  ),
  M7_item_intercept = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 | target_item) + (1 | picture_file)"
  ),
  M8_drop_participant_typicality = paste0(
    "(1 + condition_c || participant) + (1 | target_item) + (1 | picture_file)"
  ),
  M9_random_intercepts =
    "(1 | participant) + (1 | target_item) + (1 | picture_file)"
)

l2_fit <- fit_candidate_set(
  "l2_proficiency", l2, "log_picture_rt",
  "proficiency_z * condition_c * picture_typicality_c",
  l2_random, "lmm"
)

# Combined overview and software record.
overview <- bind_rows(
  prt_fit$overview,
  accuracy_fit$overview,
  srt_fit$overview,
  l2_fit$overview
)
write_result(overview, file.path(results_dir, "model_overview.csv"))

capture.output(sessionInfo(), file = file.path(package_dir, "sessionInfo.txt"))

message("All candidate model sets completed. Package: ", package_dir)

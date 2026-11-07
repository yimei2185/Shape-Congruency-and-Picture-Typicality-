# L1-only picture-RT analysis.
# Fixed effects remain constant while theoretically estimable random slopes are
# simplified in a pre-specified nested path. Candidate LMMs use ML; the selected
# model is refitted with REML. Nakagawa R2 is verified independently because
# insight/performance 1.4.2/0.15.2 can omit intercept variance for `||` models.

required <- c(
  "readr", "dplyr", "lme4", "lmerTest", "broom.mixed", "emmeans", "MuMIn"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(emmeans)
  library(MuMIn)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript so the package root can be located.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
package_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
data_dir <- file.path(package_dir, "data")
models_dir <- file.path(package_dir, "models")
results_dir <- file.path(package_dir, "results", "l1_only_prt")
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

options(contrasts = c("contr.sum", "contr.poly"))
emm_options(lmerTest.limit = 5000, pbkrtest.limit = 5000)
singular_tol <- 1e-4
boundary_tol <- 1e-7

control <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 300000),
  check.conv.singular = "ignore"
)

write_result <- function(x, name) {
  write_excel_csv(as.data.frame(x), file.path(results_dir, name), na = "")
}

fit_error <- function(message) structure(list(error = message), class = "fit_error")
is_fit_error <- function(x) inherits(x, "fit_error")

safe_lmer <- function(formula, data, reml = FALSE) {
  tryCatch(
    lmer(formula, data = data, REML = reml, control = control),
    error = function(e) fit_error(conditionMessage(e))
  )
}

convergence_message <- function(model) {
  if (is_fit_error(model)) return(model$error)
  msg <- model@optinfo$conv$lme4$messages
  if (is.null(msg)) "" else paste(msg, collapse = "; ")
}

diagnostic_row <- function(fit, name, formula_text, order) {
  if (is_fit_error(fit)) {
    return(tibble(
      order = order, model = name, formula = formula_text, fitted = FALSE,
      converged = FALSE, singular = NA, npar = NA_integer_, nobs = NA_integer_,
      logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
      min_random_variance = NA_real_, boundary_terms_below_1e_7 = "",
      convergence_message = fit$error
    ))
  }
  vc <- as.data.frame(VarCorr(fit)) |> filter(is.na(var2))
  boundary <- vc |> filter(vcov < boundary_tol)
  ll <- logLik(fit)
  conv <- convergence_message(fit)
  tibble(
    order = order,
    model = name,
    formula = formula_text,
    fitted = TRUE,
    converged = identical(conv, ""),
    singular = isSingular(fit, tol = singular_tol),
    npar = attr(ll, "df"),
    nobs = nobs(fit),
    logLik = as.numeric(ll),
    AIC = AIC(fit),
    BIC = BIC(fit),
    min_random_variance = min(vc$vcov),
    boundary_terms_below_1e_7 = if (nrow(boundary)) {
      paste(paste0(boundary$grp, ":", boundary$var1), collapse = "; ")
    } else "",
    convergence_message = conv
  )
}

adjacent_lrt <- function(models, diagnostics) {
  out <- vector("list", nrow(diagnostics) - 1)
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

direct_nakagawa_r2 <- function(model) {
  fixed_variance <- var(as.numeric(predict(model, re.form = NA)))
  Z <- t(getME(model, "Zt"))
  Lambda <- getME(model, "Lambda")
  ZLambda <- Z %*% Lambda
  observation_random_variance <- rowSums(as.matrix(ZLambda)^2) * sigma(model)^2
  random_variance <- mean(observation_random_variance)
  residual_variance <- sigma(model)^2
  denominator <- fixed_variance + random_variance + residual_variance
  tibble(
    method = "mean diagonal of Z G Z-transpose",
    fixed_variance = fixed_variance,
    mean_random_variance = random_variance,
    residual_variance = residual_variance,
    marginal_R2 = fixed_variance / denominator,
    conditional_R2 = (fixed_variance + random_variance) / denominator,
    min_observation_random_variance = min(observation_random_variance),
    max_observation_random_variance = max(observation_random_variance)
  )
}

# -----------------------------------------------------------------------------
# Data and estimability checks
# -----------------------------------------------------------------------------

l1 <- read_csv(file.path(data_dir, "prt_srt_analysis.csv"), show_col_types = FALSE) |>
  filter(language == "L1") |>
  mutate(
    participant = factor(participant),
    target_item = factor(target_item),
    picture_file = factor(picture_file)
  )

stopifnot(
  nrow(l1) == 1862,
  n_distinct(l1$participant) == 80,
  n_distinct(l1$target_item) == 24,
  n_distinct(l1$picture_file) == 48,
  !anyNA(l1[c("log_picture_rt", "condition_c", "picture_typicality_c")])
)

write_excel_csv(l1, file.path(data_dir, "l1_only_prt_analysis.csv"), na = "")

estimability <- bind_rows(lapply(
  c("participant", "target_item", "picture_file"),
  function(group) {
    l1 |>
      group_by(.data[[group]]) |>
      summarise(
        observations = n(),
        condition_levels = n_distinct(condition_c),
        typicality_levels = n_distinct(picture_typicality_c),
        condition_typicality_combinations =
          n_distinct(interaction(condition_c, picture_typicality_c)),
        .groups = "drop"
      ) |>
      summarise(
        grouping_factor = group,
        groups = n(),
        min_observations = min(observations),
        max_observations = max(observations),
        min_condition_levels = min(condition_levels),
        max_condition_levels = max(condition_levels),
        min_typicality_levels = min(typicality_levels),
        max_typicality_levels = max(typicality_levels),
        min_combinations = min(condition_typicality_combinations),
        max_combinations = max(condition_typicality_combinations)
      )
  }
))
write_result(estimability, "00_random_slope_estimability.csv")

# -----------------------------------------------------------------------------
# Candidate random-effects structures
# -----------------------------------------------------------------------------

fixed <- "condition_c * picture_typicality_c"
random_candidates <- c(
  M1_correlated_maximal = paste0(
    "(1 + condition_c * picture_typicality_c | participant) + ",
    "(1 + condition_c * picture_typicality_c | target_item) + ",
    "(1 + condition_c | picture_file)"
  ),
  M2_uncorrelated_maximal = paste0(
    "(1 + condition_c * picture_typicality_c || participant) + ",
    "(1 + condition_c * picture_typicality_c || target_item) + ",
    "(1 + condition_c || picture_file)"
  ),
  M3_picture_intercept = paste0(
    "(1 + condition_c * picture_typicality_c || participant) + ",
    "(1 + condition_c * picture_typicality_c || target_item) + ",
    "(1 | picture_file)"
  ),
  M4_drop_item_interaction = paste0(
    "(1 + condition_c * picture_typicality_c || participant) + ",
    "(1 + condition_c + picture_typicality_c || target_item) + ",
    "(1 | picture_file)"
  ),
  M5_drop_participant_interaction = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + condition_c + picture_typicality_c || target_item) + ",
    "(1 | picture_file)"
  ),
  M6_drop_item_typicality = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 + condition_c || target_item) + (1 | picture_file)"
  ),
  M7_item_intercept = paste0(
    "(1 + condition_c + picture_typicality_c || participant) + ",
    "(1 | target_item) + (1 | picture_file)"
  ),
  M8_participant_condition_only = paste0(
    "(1 + condition_c || participant) + ",
    "(1 | target_item) + (1 | picture_file)"
  ),
  M9_random_intercepts =
    "(1 | participant) + (1 | target_item) + (1 | picture_file)"
)

models <- list()
rows <- list()
for (i in seq_along(random_candidates)) {
  name <- names(random_candidates)[[i]]
  formula_text <- paste("log_picture_rt ~", fixed, "+", random_candidates[[i]])
  message("[L1-only PRT] fitting ", name)
  model <- safe_lmer(as.formula(formula_text), l1, reml = FALSE)
  models[[name]] <- model
  rows[[i]] <- diagnostic_row(model, name, formula_text, i)
  write_result(bind_rows(rows), "candidate_model_comparison_partial.csv")
  saveRDS(models, file.path(models_dir, "l1_only_prt_all_candidates_partial.rds"))
}

diagnostics <- bind_rows(rows)
valid_aic <- diagnostics$AIC[diagnostics$fitted & diagnostics$converged]
diagnostics <- diagnostics |>
  mutate(delta_AIC = AIC - min(valid_aic, na.rm = TRUE))
write_result(diagnostics, "candidate_model_comparison.csv")
write_result(adjacent_lrt(models, diagnostics), "adjacent_likelihood_ratio_tests.csv")
saveRDS(models, file.path(models_dir, "l1_only_prt_all_candidates.rds"))

supported <- diagnostics |>
  filter(
    fitted, converged, !singular,
    boundary_terms_below_1e_7 == "", !is.na(AIC)
  )
if (!nrow(supported)) {
  supported <- diagnostics |>
    filter(fitted, converged, !singular, !is.na(AIC))
}
if (!nrow(supported)) stop("No converged non-singular L1-only candidate model.")

best_supported_aic <- min(supported$AIC)
selected_row <- supported |>
  filter(AIC <= best_supported_aic + 2) |>
  arrange(npar, AIC, order) |>
  slice(1)
selected_name <- selected_row$model[[1]]
selected_ml <- models[[selected_name]]
selected_reml <- safe_lmer(formula(selected_ml), l1, reml = TRUE)
if (is_fit_error(selected_reml)) stop("REML refit failed: ", selected_reml$error)

saveRDS(selected_ml, file.path(models_dir, "l1_only_prt_selected_model_ml.rds"))
saveRDS(selected_reml, file.path(models_dir, "l1_only_prt_selected_model.rds"))

# -----------------------------------------------------------------------------
# Final model tables and independently verified R2
# -----------------------------------------------------------------------------

fixed_table <- tidy(
  selected_reml, effects = "fixed", conf.int = TRUE, conf.method = "Wald"
) |>
  mutate(
    exp_beta = exp(estimate),
    percent_rt_change = 100 * (exp_beta - 1)
  )
random_table <- as.data.frame(VarCorr(selected_reml))
write_result(fixed_table, "selected_fixed_effects.csv")
write_result(random_table, "selected_random_effects.csv")

r2_direct <- direct_nakagawa_r2(selected_reml)
r2_mumin_matrix <- suppressWarnings(MuMIn::r.squaredGLMM(selected_reml))
r2_mumin <- tibble(
  method = "MuMIn::r.squaredGLMM",
  marginal_R2 = unname(r2_mumin_matrix[1, "R2m"]),
  conditional_R2 = unname(r2_mumin_matrix[1, "R2c"])
)
if (
  abs(r2_direct$marginal_R2 - r2_mumin$marginal_R2) > 1e-7 ||
  abs(r2_direct$conditional_R2 - r2_mumin$conditional_R2) > 1e-7
) stop("Direct and MuMIn R2 calculations do not agree.")

write_result(r2_direct, "selected_model_r2_direct.csv")
write_result(r2_mumin, "selected_model_r2_mumin_check.csv")

overview <- tibble(
  analysis = "l1_only_prt",
  selected_candidate = selected_name,
  final_formula = paste(deparse(formula(selected_reml)), collapse = " "),
  observations = nobs(selected_reml),
  participants = n_distinct(l1$participant),
  target_items = n_distinct(l1$target_item),
  pictures = n_distinct(l1$picture_file),
  singular = isSingular(selected_reml, tol = singular_tol),
  convergence_message = convergence_message(selected_reml),
  AIC_REML = AIC(selected_reml),
  marginal_R2 = r2_direct$marginal_R2,
  conditional_R2 = r2_direct$conditional_R2
)
write_result(overview, "selected_model_overview.csv")
capture.output(summary(selected_reml), file = file.path(results_dir, "selected_model_summary.txt"))

# -----------------------------------------------------------------------------
# Planned decomposition of Condition x Picture Typicality
# -----------------------------------------------------------------------------

typicality_sd <- l1 |>
  distinct(picture_file, picture_typicality) |>
  summarise(value = sd(picture_typicality)) |>
  pull(value)
probe_values <- c(low = -typicality_sd, mean = 0, high = typicality_sd)

condition_emm <- emmeans(
  selected_reml,
  ~ condition_c | picture_typicality_c,
  at = list(
    condition_c = c(-0.5, 0.5),
    picture_typicality_c = unname(probe_values)
  ),
  lmer.df = "satterthwaite"
)
condition_contrasts <- contrast(
  condition_emm,
  method = list("congruent - incongruent" = c(1, -1)),
  by = "picture_typicality_c",
  adjust = "none"
) |>
  summary(infer = c(TRUE, TRUE), adjust = "none") |>
  as.data.frame() |>
  mutate(
    typicality_level = case_when(
      picture_typicality_c < -1e-8 ~ "low_minus_1SD",
      picture_typicality_c > 1e-8 ~ "high_plus_1SD",
      TRUE ~ "mean"
    ),
    picture_typicality_raw = picture_typicality_c +
      mean((l1 |> distinct(picture_file, picture_typicality))$picture_typicality),
    p_holm = p.adjust(p.value, method = "holm"),
    rt_ratio_congruent_over_incongruent = exp(estimate),
    percent_difference = 100 * (exp(estimate) - 1),
    cohens_d = estimate / sigma(selected_reml)
  )
write_result(condition_emm, "emm_condition_at_typicality_levels.csv")
write_result(condition_contrasts, "condition_effects_at_typicality_levels.csv")

typicality_trends <- emtrends(
  selected_reml,
  ~ condition_c,
  var = "picture_typicality_c",
  at = list(condition_c = c(-0.5, 0.5)),
  lmer.df = "satterthwaite"
)
typicality_slopes <- summary(typicality_trends, infer = c(TRUE, TRUE)) |>
  as.data.frame() |>
  mutate(
    condition = if_else(condition_c < 0, "congruent", "incongruent"),
    percent_rt_change_per_typicality_point =
      100 * (exp(picture_typicality_c.trend) - 1)
  )
typicality_slope_difference <- pairs(typicality_trends, adjust = "none") |>
  summary(infer = c(TRUE, TRUE), adjust = "none") |>
  as.data.frame()
write_result(typicality_slopes, "picture_typicality_slopes_by_condition.csv")
write_result(typicality_slope_difference, "picture_typicality_slope_difference.csv")

software <- tibble(
  component = c("R", "lme4", "lmerTest", "emmeans", "MuMIn"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("lme4")),
    as.character(packageVersion("lmerTest")),
    as.character(packageVersion("emmeans")),
    as.character(packageVersion("MuMIn"))
  )
)
write_result(software, "software_versions.csv")

capture.output(sessionInfo(), file = file.path(results_dir, "sessionInfo.txt"))
message("L1-only PRT analysis completed. Selected: ", selected_name)

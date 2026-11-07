# Reproducible follow-up contrasts for interactions in the selected models.
# These files are descriptive/model-based follow-ups; the candidate random-
# effects selection is performed independently in 01_fit_candidate_models.R.

required <- c("readr", "dplyr", "lme4", "lmerTest", "emmeans", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

library(readr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)

emm_options(lmerTest.limit = 4000, pbkrtest.limit = 4000)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript so the package root can be located.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
package_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
models_dir <- file.path(package_dir, "models")
results_dir <- file.path(package_dir, "results")

write_result <- function(x, path) write_excel_csv(as.data.frame(x), path, na = "")

prt <- readRDS(file.path(models_dir, "picture_rt_selected_model.rds"))
accuracy <- readRDS(file.path(models_dir, "accuracy_selected_model.rds"))
prt_data <- read_csv(
  file.path(package_dir, "data", "prt_srt_analysis.csv"),
  show_col_types = FALSE
)

# The selected picture-RT model contains a significant three-way interaction.
# Report both theoretically relevant families of simple effects. Cohen's d is
# the estimated log-RT contrast divided by the selected model residual SD.
prt_condition <- emmeans(
  prt, ~ condition_c | language_c * picture_typicality_c,
  at = list(language_c = c(-0.5, 0.5),
            condition_c = c(-0.5, 0.5),
            picture_typicality_c = c(-0.5, 0.5)),
  lmer.df = "satterthwaite"
)
prt_condition_pairs <- pairs(prt_condition, adjust = "tukey")
prt_condition_d <- eff_size(
  prt_condition, sigma = sigma(prt), edf = df.residual(prt), method = "pairwise"
)

prt_typicality <- emmeans(
  prt, ~ picture_typicality_c | language_c * condition_c,
  at = list(language_c = c(-0.5, 0.5),
            condition_c = c(-0.5, 0.5),
            picture_typicality_c = c(-0.5, 0.5)),
  lmer.df = "satterthwaite"
)
prt_typicality_pairs <- pairs(prt_typicality, adjust = "tukey")
prt_typicality_d <- eff_size(
  prt_typicality, sigma = sigma(prt), edf = df.residual(prt), method = "pairwise"
)

write_result(prt_condition, file.path(results_dir, "picture_rt", "emm_condition_cells.csv"))
write_result(prt_condition_pairs, file.path(results_dir, "picture_rt", "simple_condition_contrasts.csv"))
write_result(prt_condition_d, file.path(results_dir, "picture_rt", "simple_condition_cohens_d.csv"))
write_result(prt_typicality, file.path(results_dir, "picture_rt", "emm_typicality_cells.csv"))
write_result(prt_typicality_pairs, file.path(results_dir, "picture_rt", "simple_typicality_contrasts.csv"))
write_result(prt_typicality_d, file.path(results_dir, "picture_rt", "simple_typicality_cohens_d.csv"))

# Publication figure for the significant Language x Congruency x Picture
# Typicality interaction. The three probe values match the manuscript table:
# the mean and one picture-level SD below and above the mean. Estimates and
# confidence limits are bias-adjusted back-transformations from log ms.
typicality_reference <- prt_data |>
  distinct(picture_file, picture_typicality) |>
  summarise(
    typicality_mean = mean(picture_typicality),
    typicality_sd = sd(picture_typicality)
  )
probe_values <- c(
  low = -typicality_reference$typicality_sd,
  mean = 0,
  high = typicality_reference$typicality_sd
)

prt_probe_emm <- emmeans(
  prt, ~ condition_c | language_c * picture_typicality_c,
  at = list(
    language_c = c(-0.5, 0.5),
    condition_c = c(-0.5, 0.5),
    picture_typicality_c = unname(probe_values)
  ),
  lmer.df = "satterthwaite"
)

prt_probe_plot_data <- summary(prt_probe_emm, infer = c(TRUE, TRUE)) |>
  as.data.frame() |>
  mutate(
    language = factor(
      if_else(language_c < 0, "L1", "L2"),
      levels = c("L1", "L2")
    ),
    condition = factor(
      if_else(condition_c < 0, "Congruent", "Incongruent"),
      levels = c("Congruent", "Incongruent")
    ),
    typicality_level = factor(
      case_when(
        picture_typicality_c < -1e-8 ~ "Low (-1 SD)",
        picture_typicality_c > 1e-8 ~ "High (+1 SD)",
        TRUE ~ "Mean"
      ),
      levels = c("Low (-1 SD)", "Mean", "High (+1 SD)")
    ),
    picture_typicality_raw = picture_typicality_c + typicality_reference$typicality_mean,
    estimate_ms = exp(emmean + sigma(prt)^2 / 2),
    lower_ms = exp(lower.CL + sigma(prt)^2 / 2),
    upper_ms = exp(upper.CL + sigma(prt)^2 / 2)
  )

prt_probe_contrasts <- contrast(
  prt_probe_emm,
  method = list("congruent - incongruent" = c(1, -1)),
  by = c("language_c", "picture_typicality_c"),
  adjust = "none"
) |>
  summary(infer = c(TRUE, TRUE), adjust = "none") |>
  as.data.frame() |>
  mutate(
    language = if_else(language_c < 0, "L1", "L2"),
    typicality_level = case_when(
      picture_typicality_c < -1e-8 ~ "Low (-1 SD)",
      picture_typicality_c > 1e-8 ~ "High (+1 SD)",
      TRUE ~ "Mean"
    ),
    p_holm = ave(p.value, language, FUN = function(x) p.adjust(x, method = "holm")),
    rt_ratio_congruent_over_incongruent = exp(estimate),
    percent_difference = 100 * (exp(estimate) - 1),
    cohens_d = estimate / sigma(prt)
  )

write_result(
  prt_probe_plot_data,
  file.path(results_dir, "picture_rt", "figure_prt_three_way_plot_data.csv")
)
write_result(
  prt_probe_contrasts,
  file.path(results_dir, "picture_rt", "condition_effects_at_mean_and_sd.csv")
)
write_result(
  typicality_reference,
  file.path(results_dir, "picture_rt", "picture_typicality_reference.csv")
)

prt_three_way_figure <- ggplot(
  prt_probe_plot_data,
  aes(
    x = typicality_level,
    y = estimate_ms,
    colour = condition,
    shape = condition,
    group = condition
  )
) +
  geom_line(linewidth = 0.75) +
  geom_errorbar(
    aes(ymin = lower_ms, ymax = upper_ms),
    width = 0.09,
    linewidth = 0.55
  ) +
  geom_point(size = 2.8, fill = "white", stroke = 0.8) +
  facet_wrap(~ language, nrow = 1) +
  scale_colour_manual(values = c("Congruent" = "#0072B2", "Incongruent" = "#D55E00")) +
  scale_shape_manual(values = c("Congruent" = 16, "Incongruent" = 17)) +
  labs(
    x = "Picture Typicality",
    y = "Model-Estimated PRTs (ms)",
    colour = "Condition",
    shape = "Condition"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.spacing.x = grid::unit(1.5, "lines"),
    axis.title = element_text(face = "bold")
  )

ggsave(
  file.path(results_dir, "picture_rt", "figure_prt_three_way_interaction.pdf"),
  prt_three_way_figure, width = 7.2, height = 4.3, units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(results_dir, "picture_rt", "figure_prt_three_way_interaction.png"),
  prt_three_way_figure, width = 7.2, height = 4.3, units = "in",
  dpi = 600, bg = "white"
)

# Accuracy has a significant Language x Picture Typicality interaction.
# For the binomial model, response-scale EMMs are probabilities and contrasts
# are reported as odds ratios with Tukey adjustment.
acc_typicality <- emmeans(
  accuracy, ~ picture_typicality_c | language_c,
  at = list(language_c = c(-0.5, 0.5),
            condition_c = 0,
            picture_typicality_c = c(-0.5, 0.5)),
  type = "response"
)
acc_typicality_pairs <- pairs(acc_typicality, adjust = "tukey", type = "response")

acc_language <- emmeans(
  accuracy, ~ language_c | picture_typicality_c,
  at = list(language_c = c(-0.5, 0.5),
            condition_c = 0,
            picture_typicality_c = c(-0.5, 0.5)),
  type = "response"
)
acc_language_pairs <- pairs(acc_language, adjust = "tukey", type = "response")

write_result(acc_typicality, file.path(results_dir, "accuracy", "emm_typicality_probabilities.csv"))
write_result(acc_typicality_pairs, file.path(results_dir, "accuracy", "simple_typicality_odds_ratios.csv"))
write_result(acc_language, file.path(results_dir, "accuracy", "emm_language_probabilities.csv"))
write_result(acc_language_pairs, file.path(results_dir, "accuracy", "simple_language_odds_ratios.csv"))

message("Follow-up contrasts completed.")

##############################################################################
## Seasonal step-selection functions for adult female white-tailed deer
## Clevinger et al. -- Bath County, Virginia, USA, 2019-2021
##
## Run from the repository root:   Rscript R/fit_ssf_models.R
##
## Input:   data/ssf_steps.csv.gz   (used + available steps with covariates)
## Outputs: results/ssf_results_<season>.csv   full model output per season
##          results/table1_ssf_coefficients.csv manuscript Table 1
##          results/covariate_scaling.csv       means/SDs used for scaling
##          results/figure_ssf_coefficients.png coefficient plot
##          results/sessionInfo.txt
##
## Model, per biological season (conditional Poisson SSF; Muff et al. 2020):
##   y ~ -1 + d2emsuc + d2edge + d2herb + elev
##       + uncorrelated random slopes for each covariate by animal
##       + (1 | stratum), stratum-intercept SD fixed at 1000
##
## Changes from the original analysis script (original/):
##   1. Strata are unique to animal x year x season x step. The original
##      stratum ID (collar_step) repeated across years, so up to three
##      steps from different years shared one stratum.
##   2. Covariates are standardised across the whole dataset (all animals,
##      seasons and years). The original standardised within each
##      animal-season-year.
##   3. Distance to mature forest (d2forest) is excluded: it equals distance
##      to forested edge at 98.6% of steps (r = 0.97), so both cannot be
##      estimated in the same model.
##############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(glmmTMB)
  library(broom.mixed)
  library(ggplot2)
})

dir.create("results", showWarnings = FALSE)

seasons  <- c("Fawning", "Breeding", "Early Gestation", "Late Gestation")
raw_covs <- c("d2emsuc", "d2edge", "d2herb", "elev")
covs     <- paste0(raw_covs, "_sc")
labels   <- c(d2emsuc_sc = "Disturbed forest", d2edge_sc = "Forested edge",
              d2herb_sc  = "Open",             elev_sc   = "Elevation")

#### Data ####################################################################

steps <- read.csv("data/ssf_steps.csv.gz") %>%
  filter(if_all(all_of(raw_covs), ~ !is.na(.x))) %>%
  mutate(
    y = case_,
    step_an_id = paste(CollarSerialNumber, year, bio_season, step_id_, sep = "_"),
    across(all_of(raw_covs), ~ (.x - mean(.x)) / sd(.x), .names = "{.col}_sc")
  ) %>%
  # keep complete strata: exactly one used step plus its available steps
  group_by(step_an_id) %>%
  filter(sum(y) == 1, n() > 1) %>%
  ungroup()

scaling <- steps %>%
  summarise(across(all_of(raw_covs), list(mean = mean, sd = sd))) %>%
  pivot_longer(everything(), names_to = c("covariate", ".value"),
               names_pattern = "(.*)_(mean|sd)")
write.csv(scaling, "results/covariate_scaling.csv", row.names = FALSE)

#### Models ##################################################################

fit_ssf <- function(dat) {
  form <- as.formula(paste(
    "y ~ -1 +", paste(covs, collapse = " + "), "+",
    paste0("(0 + ", covs, " | CollarSerialNumber)", collapse = " + "),
    "+ (1 | step_an_id)"))
  k <- length(covs)
  str <- glmmTMB(form, family = poisson, data = dat,
                 map   = list(theta = factor(c(1:k, NA))),   # fix stratum SD
                 start = list(theta = c(rep(0, k), log(1e3))),
                 doFit = FALSE)
  glmmTMB:::fitTMB(str)
}

results <- list()
for (s in seasons) {
  dat <- filter(steps, bio_season == s)
  message(sprintf("%s: %d deer, %d used steps", s,
                  n_distinct(dat$CollarSerialNumber), sum(dat$y)))
  m <- fit_ssf(dat)
  out <- tidy(m) %>% mutate(season = s, .before = 1)
  write.csv(out, sprintf("results/ssf_results_%s.csv", s), row.names = FALSE)
  results[[s]] <- out %>%
    filter(effect == "fixed") %>%
    mutate(n_deer = n_distinct(dat$CollarSerialNumber), n_used_steps = sum(dat$y))
  rm(m, dat); invisible(gc())
}

#### Table 1 #################################################################

table1 <- bind_rows(results) %>%
  left_join(scaling %>% mutate(term = paste0(covariate, "_sc")), by = "term") %>%
  mutate(
    resource = labels[term],
    interpretation = case_when(
      p.value >= 0.05 ~ "Not Selected",
      term == "elev_sc" & estimate > 0 ~ "Selected",
      term == "elev_sc" ~ "Avoided",
      estimate < 0 ~ "Selected",       # closer = negative coefficient
      TRUE ~ "Avoided"),
    # % change in relative selection strength per +100 m
    pct_change_per_100m = 100 * (exp(estimate * 100 / sd) - 1)
  ) %>%
  transmute(season, resource, beta = round(estimate, 3), SE = round(std.error, 3),
            z = round(statistic, 3),
            p = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
            interpretation, pct_change_per_100m = round(pct_change_per_100m, 1),
            n_deer, n_used_steps)
write.csv(table1, "results/table1_ssf_coefficients.csv", row.names = FALSE)
print(table1, n = Inf)

#### Figure ##################################################################

fig <- bind_rows(results) %>%
  mutate(resource = factor(labels[term], levels = rev(labels)),
         season = factor(season, levels = seasons)) %>%
  ggplot(aes(x = estimate, y = resource)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error), height = 0.2) +
  geom_point(size = 2) +
  facet_wrap(~ season, nrow = 1) +
  labs(x = "Selection coefficient (estimate and 95% CI)", y = NULL) +
  theme_bw()
ggsave("results/figure_ssf_coefficients.png", fig, width = 10, height = 3.5, dpi = 300)

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("Done.")

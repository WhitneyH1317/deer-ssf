library(glmmTMB)
library(tidyverse)
library(raster)
library(survival)
library(TwoStepCLogit)
library(amt)

# Tibble with nested, animal-season-year specific data
cumu.tib <- readRDS("./glmmtmb_data4ssf.RDS")

cumu.tib

# Unnest data and renest by season for seasonal models
tmbData <- cumu.tib %>%
  dplyr::select(1,2,3,6) %>%
  mutate(bio_season = factor(bio_season, levels = c("Fawning", "Early Gestation", "Late Gestation", "Breeding"))) %>%
  unnest(data) %>%
  mutate(y = ifelse(case_ == T, 1, 0),
         CollarSerialNumber = as.numeric(factor(CollarSerialNumber))) %>%
  nest(-bio_season)

# Sampe Size Table
cumu.tib %>%
  dplyr::select(1,2,3,6) %>%
  mutate(bio_season = factor(bio_season, levels = c("Fawning", "Early Gestation", "Late Gestation", "Breeding"))) %>%
  unnest(data) %>%
  mutate(y = ifelse(case_ == T, 1, 0),
         CollarSerialNumber = as.numeric(factor(CollarSerialNumber))) %>%
  filter(case_ == 1) %>%
  nest(-CollarSerialNumber, -year, -bio_season) %>%
  mutate(N = map_dbl(data, ~nrow(.))) %>%
  dplyr::select(year, bio_season, CollarSerialNumber, N) %>%
  arrange(year, bio_season, CollarSerialNumber) %>%
  write.csv("./Results/SSF_sample_size_table.csv")

##### Summary Figures #####
##Turn Angle Distributions
windows()
tmbData %>%
  unnest(data) %>%
  filter(case_ == 1) %>%
  ggplot(., aes(x = ta_)) +
  geom_density() +
  xlab("Turn Angle (radians)") +
  ylab("Density") +
  # xlim(0,750) +
  xlim(-3.3, 3.3) +
  facet_wrap(vars(bio_season)) +
  theme(panel.background = element_rect(colour = "black", fill = "transparent"),
        panel.grid = element_line(colour = "transparent"),
        axis.title.x.bottom = element_text(margin = margin(12,0,0,0)),
        axis.title.y.left = element_text(margin = margin(0,12,0,0)),
        axis.text.y.left = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8))

## Step Length Distributions
windows()
tmbData %>%
  unnest(data) %>%
  filter(case_ == 1) %>%
  ggplot(., aes(x = sl_)) +
  geom_density() +
  xlab("Step Length (m)") +
  ylab("Density") +
  xlim(0,750) +
  # xlim(-3.3, 3.3) +
  facet_wrap(vars(bio_season)) +
  theme(panel.background = element_rect(colour = "black", fill = "transparent"),
        panel.grid = element_line(colour = "transparent"),
        axis.title.x.bottom = element_text(margin = margin(12,0,0,0)),
        axis.title.y.left = element_text(margin = margin(0,12,0,0)),
        axis.text.y.left = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8))

# -CollarSerialNumber, -year) %>%
#   arrange(year, bio_season, CollarSerialNumber)

resultsList <- list()

# Loop through season-specific data and execute models
for(i in 1:4){

  TMBStruc <- glmmTMB(y ~ -1 + d2emsuc_sc + d2edge_sc + d2herb_sc + d2forest_sc + elev_sc +
                        (0 + d2emsuc_sc | CollarSerialNumber) +
                        (0 + d2edge_sc | CollarSerialNumber) +
                        (0 + d2herb_sc | CollarSerialNumber) +
                        (0 + d2forest_sc | CollarSerialNumber) +
                        (0 + elev_sc | CollarSerialNumber) +
                        (1 | step_an_id),
                      family=poisson,
                      data = tmbData$data[[i]],
                      map = list(theta = factor(c(1:5,NA))), # Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                      start = list(theta = c(rep(0, 5),log(1e3))), #' Set the value of the standard deviation of the first random effect (here (1|step_id)):
                      doFit=FALSE)
  runtime <- system.time({
    glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
  })
  print(runtime)
  saveRDS(glmm.TMB.random, paste("./Results/random_slopes_ssf_",tmbData$bio_season[i],".RDS", sep = ""))
  write.csv(tidy(glmm.TMB.random), paste("./Results/random_slopes_ssf_results_",tmbData$bio_season[[i]],".csv",sep=""))
  resultsList[[i]] <- tidy(glmm.TMB.random)
}

# Make list of seasonal model result files
results_files <- list.files("./Results", pattern = "random_slopes_ssf", full.names = T)[grep(".csv",list.files("./Results", pattern = "random_slopes_ssf"))]

# Loop through files and clean them for table
for(i in 1:4){
  result.df <- read.csv(results_files[i])
  result.df$season <- ifelse(i == 1, "Breeding",
                             ifelse(i == 2, "Early Gestation",
                                    ifelse(i == 3, "Fawning","Late Gestation")))
  result.df <- result.df %>%
    filter(effect == "fixed")
  if(i == 1){
    ssf_results <- result.df
  } else {
    ssf_results <- rbind(ssf_results, result.df)
  }
}


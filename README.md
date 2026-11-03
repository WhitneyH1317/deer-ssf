# Seasonal step selection by female white-tailed deer

Data and code for the step-selection analysis in:

> Clevinger et al. *Seasonal resource selection by white-tailed deer reveals …* (in prep.)

Adult female white-tailed deer (*Odocoileus virginianus*) were GPS-collared in Bath County, Virginia, USA (2019–2021). Hourly fixes were used to fit season-specific step-selection functions (SSFs) for four biological seasons: Early Gestation (Dec 1–Feb 28), Late Gestation (Mar 1–May 31), Fawning (Jun 1–Aug 31) and Breeding (Sep 1–Nov 30).

## Repository contents

```
data/ssf_steps.csv.gz              used and available steps with covariate values
R/fit_ssf_models.R                 fits the seasonal SSFs and writes all results
results/                           model output, Table 1, figure, session info
original/                          original analysis script, kept for provenance
```

## How to run

From the repository root, in R ≥ 4.2:

```r
install.packages(c("dplyr", "tidyr", "glmmTMB", "broom.mixed", "ggplot2"))
```
```bash
Rscript R/fit_ssf_models.R
```

It takes about 10 minutes and about 6 GB of RAM.

## Data: `data/ssf_steps.csv.gz`

There is one row per step. Each used step (`case_ = 1`) has 5 available steps (`case_ = 0`) drawn from step-length and turn-angle distributions estimated from all collared deer. Coordinates are UTM zone 17N (NAD83, EPSG:26917). Times are in EST.

| Column | Description |
|---|---|
| `CollarSerialNumber` | Animal ID. Collar 41467 was deployed on two deer (`41467a`, `41467b`). |
| `year`, `bio_season` | Biological year and season. Early Gestation is assigned to the year in which it starts (December). |
| `burst_`, `step_id_` | Burst and step number within each animal-season-year. |
| `case_` | 1 = used step, 0 = available step. |
| `t1_`, `t2_` | Start and end time of the step. |
| `x1_`, `y1_`, `x2_`, `y2_` | Start and end coordinates (m). |
| `sl_`, `ta_` | Step length (m) and turn angle (radians). |
| `d2emsuc` | Distance (m) to disturbed forest (early/mid-successional: timber harvest or fire-created gaps within 40 years). |
| `d2edge` | Distance (m) to forested edge (40 m buffer into mature forest next to open or disturbed patches). |
| `d2herb` | Distance (m) to open areas (wildlife openings, clearings, pasture). |
| `d2forest` | Distance (m) to mature forest. **Not used in the models; see below.** |
| `elev` | Elevation (m), from the USGS DEM. |

Covariates were extracted at step end points from 10 m rasters. The source rasters (≈800 MB) are not included here and are available from the authors.

Deer-season-years were included only if the collar recorded >80% of scheduled fixes. The first 14 days after capture were removed.

## Model

For each season, `R/fit_ssf_models.R` fits a conditional Poisson model with `glmmTMB` (Muff, Signer & Fieberg 2020, *J. Anim. Ecol.* 89:80–92):

- **Fixed effects:** distance to disturbed forest, forested edge and open areas, and elevation.
- **Random slopes:** an uncorrelated random slope for every covariate by animal.
- **Strata:** a random intercept for each stratum (one used step plus its available steps), with its SD fixed at 1000.

Covariates are standardised across the full dataset; the means and SDs are written to `results/covariate_scaling.csv`. For the distance covariates, a negative coefficient means selection: use is more likely closer to that cover type.

`results/table1_ssf_coefficients.csv` also gives `pct_change_per_100m`, the percentage change in relative selection strength per 100 m increase in each covariate.

## Changes from the original analysis

The original script (`original/`) and its input reproduce the coefficients first reported in the manuscript exactly. Three problems were then corrected:

1. **Strata were not unique.** The original stratum ID (`collar_step`) repeated across years, so in about two-thirds of strata, up to three used steps from different years shared one stratum. Strata are now unique to animal × year × season × step.
2. **Covariates were standardised within each animal-season-year.** Coefficients therefore had no common unit across animals and could not be converted to distances. Covariates are now standardised across the whole dataset.
3. **Distance to mature forest was collinear with distance to forested edge.** The two were identical at 98.6% of steps (r = 0.97), so their effects could not be separated. Mature forest was removed from the models.

## Citation and license

Please cite the paper when using these data or code. *(Add a license, e.g. CC-BY 4.0 for data and MIT for code, and a DOI once archived.)*

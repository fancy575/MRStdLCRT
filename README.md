# MRStdLCRT

Model-Robust Standardization for Longitudinal Cluster-Randomized Trials (LCRTs), with a focus on multi-period crossover CRTs, parallel CRT and stepped-wedge–type treatment patterns.

This package implements **unadjusted** and **augmented (model-robust)** estimators for four longitudinal CRT estimands, with uncertainty quantified via **delete-1-cluster jackknife**.

## Key features

- Supports **multi-period** longitudinal CRT data with a **cluster** and **period** structure.
- Computes four estimands:
  - **h-iATE**: horizontal individual ATE (individual-weighted across clusters)
  - **h-cATE**: horizontal cluster ATE (cluster-weighted / cluster-period mean–based)
  - **v-iATE**: vertical individual ATE (period-weighted across periods)
  - **v-cATE**: vertical cluster-period ATE (cluster-period–weighted across periods)
- Two estimator types:
  - **Unadjusted**: design-based mean contrasts using cluster-period means.
  - **Adjusted / augmented (MRS)**: working-model predictions standardized and combined with an augmentation term to improve robustness.
- Working model options:
  - `method = "gee"` (via **gee**), with user-specified working correlation `corstr`
  - `method = "lmer"` (via **lme4**) for Gaussian outcomes
  - `method = "glmer"` (via **lme4**) for Binomial outcomes
- For longitudinal trials where **some periods are not mixed** (e.g., stepped-wedge edges), aggregation is performed over **kept periods** where both treatment arms are present (mixed assignment).

## Installation

### From GitHub

```r
# install.packages("remotes")
remotes::install_github("fancy575/MRStdLCRT")

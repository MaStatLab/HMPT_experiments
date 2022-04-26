# HMPT_experiments

This repository provides program codes to run the SMC algorithm proposed in Awaya and Ma (2022) under the settings of Section 5 (numerical experiments). 

The detailed information is provided in

Awaya, N., & Ma, L. (2020). Hidden Markov P\'olya trees for high-dimensional distributions. arXiv preprint arXiv:2011.03121.

To run the codes, the R package `SMCMP` needs to be installed  using `devtools` as follows:

```
library(devtools)
install_github("nawaya040/SMCMP")
```

The following R files are included:

1. `Density_estimation_2d.R`: 2D density estimation (Section 5.1.1)
2. `Density_estimation_multi.R`: higher dimensional density estimation (Section 5.1.2)
3. `Two_sample_comparison.R`: two-sample comparison (Section 5.2)
4. `models.R`: functions to simulate data sets

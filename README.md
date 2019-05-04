# mgp
Multivariate generalized Pareto distributions

[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) 



This package is under development: it is not ready for use and will be merged in mev in due time.

R-package to fit multivariate generalized Pareto processes, including functions to calculate the (un)censored likelihood of the process, 
conditional Gaussian densities, spatial dependence models and a pseudo-marginal algorithm to simulate from the posterior. Some routines
for latent Gaussian modelling of generalized Pareto are included.

The package currently lacks examples and references.


To install from Github, use 

```R
devtools::install_github("lbelzile/TruncatedNormal")
devtools::install_github("lbelzile/mev")
devtools::install_github("lbelzile/mgp")
```

after installing `devtools`.

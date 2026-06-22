# CA-BIBT: Covariate-Assisted Bayesian Intransitive Bradley-Terry  

[![R](https://img.shields.io/badge/R-%23276DC3.svg?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![Rcpp](https://img.shields.io/badge/Rcpp-%23999999.svg?style=flat&logo=c%2B%2B&logoColor=white)](https://rcpp.org/)

This repository provides R/Rcpp code for the **Covariate-Assisted Bayesian Intransitive Bradley-Terry (CA-BIBT)** model. 
The CA-BIBT model is a Bayesian framework for binary pairwise comparison data that may contain intransitive structures.
It extends the classical Bradley-Terry model by representing latent pairwise match-up effects as edge flows and decomposing them through combinatorial Hodge theory into covariate-induced and residual components.

This repository accompanies the following paper:
> Okahara, H., Nakagawa, T., and Sugasawa, S. (2026). *The Covariate-Assisted Bayesian Intransitive Bradley-Terry Model via Combinatorial Hodge Theory*. arXiv:2601.07158.

## OVERVIEW
Classical Bradley-Terry models represent pairwise comparisons through scalar latent strengths and therefore impose a transitive structure.  
This can be restrictive when pairwise outcomes contain cyclic patterns, such as $A \succ B \succ C \succ A$.  
The CA-BIBT model addresses this limitation by decomposing the latent match-up flow into three components:
```math
  \boldsymbol{M} = G\boldsymbol{s} + C^\top \boldsymbol{\Phi} + X_E^\top \boldsymbol{\beta}.
```

Here,
- **Covariate flow ($X_E^\top \boldsymbol{\beta}$):** explains pair-specific comparison effects using observed edge-level covariates.
- **Residual gradient flow ($G\boldsymbol{s}$):** represents the transitive strength not explained by the covariates.
- **Residual curl flow ($C^\top \boldsymbol{\Phi}$):** captures cycle-induced intransitive structure not explained by the covariates.  
The covariate flow may contain both gradient and curl components.  

This allows the model to distinguish whether observed intransitivity is explained by covariates or remains as residual cycle-induced structure.

## Contents
- `main.R`        : The primary script to fit the CA-BIBT model and reproduce results from the manuscript.
- `functions.cpp` : C++ implementation (Rcpp) of the Gibbs sampling loop for the CA-BIBT model for faster computation.
- `functions.R`   : Core utility functions, including the Gibbs sampler, data generation, and visualization tools.
- `libraries.R`   : Loads the required R packages for the project.
- `database.R`    : Example dataset used in the manuscript.
- `RJMCMC alg`    : External implementation of the **Intransitive Clustered Bradley-Terry (ICBT)** model (Spearing et al., 2023), sourced from [JessSpearing/ICBT_model](https://github.com/JessSpearing/ICBT_model).

## Getting Started
1. Ensure you have a C++ compiler installed for `Rcpp` compatibility.
2. Execute `main.R` to run a demonstration using the provided animal dominance or synthetic datasets.

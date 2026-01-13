# BIBT: Bayesian Intransitive Bradley-Terry  

[![R](https://img.shields.io/badge/R-%23276DC3.svg?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![Rcpp](https://img.shields.io/badge/Rcpp-%23999999.svg?style=flat&logo=c%2B%2B&logoColor=white)](https://rcpp.org/)

This repository provides an implementation of the **Bayesian Intransitive Bradley-Terry (BIBT)** model. 
The BIBT model is a principled Bayesian framework designed to disentangle transitive and intransitive structures in pairwise comparison data.

This repository accompanies the following paper:
> Okahara, H., Nakagawa, T., and Sugasawa, S. (2026). *The Bayesian Intransitive Bradley-Terry Model via Combinatorial Hodge Theory*. arXiv:2601.07158.

## OVERVIEW
The BIBT model extends the classical Bradley-Terry framework by embedding **combinatorial Hodge theory**. 
It decomposes paired relationships into:
- **Gradient Flow:** Represents the transitive strength of entities.
- **Curl Flow:** Captures cycle-induced patterns (e.g., $A \succ B \succ C \succ A$).

## Contents
- `main.R`        : The primary script to fit the BIBT model and reproduce results from the manuscript.
- `BIBT.cpp`      : C++ implementation (Rcpp) of the Gibbs sampling loop for the BIBT model for faster computation.
- `functions.R`   : Core utility functions, including the Gibbs sampler, data generation, and visualization tools.
- `libraries.R`   : Loads the required R packages for the project.
- `database.R`    : Example dataset used in the manuscript.
- `baseball/`     : Contains MLB game outcomes (2020–2025) sourced from [Retrosheet](https://www.retrosheet.org/).
- `BBT.stan`      : Stan implementation of the **Bayesian Bradley-Terry (BBT)** model. The original framework is adapted from [jwainer/bbtcomp](https://github.com/jwainer/bbtcomp) for baseline comparison.
- `RJMCMC alg`    : External implementation of the **Intransitive Clustered Bradley-Terry (ICBT)** model (Spearing et al., 2023), sourced from [JessSpearing/ICBT_model](https://github.com/JessSpearing/ICBT_model).

## Getting Started
1. Ensure you have a C++ compiler installed for `Rcpp` compatibility.
2. Run `libraries.R` to set up the environment.
3. Execute `main.R` to run a demonstration using the provided MLB or synthetic datasets.

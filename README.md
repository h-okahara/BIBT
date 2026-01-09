# BIBT: Bayesian Intransitive Bradley-Terry  

## OVERVIEW
This repository provides an implementation of the **Bayesian Intransitive Bradley-Terry (BIBT)** model for ranking estimation.  
The BIBT model extends the classic Bradley-Terry framework by introducing **curl flow** that explicitly capture *intransitive* patterns (e.g., rock–paper–scissors structures) in pairwise comparison data.  
Efficient Bayesian inference is carried out using **Pólya-Gamma data augmentation**, together with *Horseshoe shrinkage prior* to regularize the cycle-induced effects.

## CONTENTS
1. `libraries.R`    : Loads the required R packages for the project.
2. `database.R`     : Provides an example dataset used in the manuscript.
3. `functions.R`    : Contains core utility functions (e.g., Gibbs sampler, data generation, visualization).
4. `main.R`         : Main script to fit the BIBT model to data and reproduce key results.
5. `BBT.stan`       : Stan code for the Bayesian Bradley-Terry (BBT) model, used as a baseline comparison.
6. `BIBT.cpp`        : C++ implementation (Rcpp) of the Gibbs sampling loop for the BIBT model for faster computation.

# IBT: Intransitive Bradley-Terry

## OVERVIEW
This repository provides an implementation of the **Intransitive Bradley-Terry (IBT)** model for ranking estimation.  
The IBT model extends the classic Bradley-Terry framework by introducing **triangular parameters** that explicitly capture *intransitive* patterns (e.g., rock–paper–scissors structures) in pairwise comparison data.  
Efficient Bayesian inference is carried out using **Pólya-Gamma data augmentation**, together with *Horseshoe shrinkage prior* to regularize the cyclic effects.

## CONTENTS
1. `libraries.R`    : Loads the required R packages for the project.
2. `database.R`     : Provides an example dataset used in the manuscript.
3. `functions.R`    : Contains core utility functions (e.g., Gibbs sampler, data generation, visualization).
4. `main.R`         : Main script to fit the CBT model to data and reproduce key results.
5. `BBT.stan`       : Stan code for the Bayesian Bradley-Terry (BBT) model, used as a baseline comparison.
6. `IBT.cpp`        : C++ implementation (Rcpp) of the Gibbs sampling loop for the IBT model for faster computation.

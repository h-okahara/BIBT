# TDBT: Trans‑Dimensional Bradley‑Terry

## Overview
This repository provides an implementation of the Trans‑Dimensional Bradley‑Terry (TDBT) model for ranking estimation. It extends the classic Bradley‑Terry framework to allow multi‑dimensional latent “worth” scores and automatic selection of the appropriate dimensionality. Efficient Gibbs sampling is achieved via Polya‑Gamma data augmentation, enabling accurate ranking from pairwise comparison data.

## Key Features
- **Multi‑Dimensional Worth**: Generalizes the univariate Bradley‑Terry model to estimate latent worth in an arbitrary number of dimensions.  
- **Automatic Dimension Selection**: Employs a multiplicative gamma process to sparsify extra dimensions, letting the data determine the number of relevant latent factors.  
- **Polya‑Gamma Sampling**: Leverages data augmentation for the binomial likelihood to implement a concise and fast Gibbs sampler.  
- **Visualization & Diagnostics**: Includes functions for trace plots, violin plots, ESS/R̂ diagnostics, and credible‑interval computation.

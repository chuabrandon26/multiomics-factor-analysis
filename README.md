# Multiomics-factor-analysis
Unsupervised multi-omics integration of CLL data using MOFA+ and MEFISTO to identify latent drivers of disease variation.

# Multi-Modal Factor Analysis of CLL
This repository implements a systems biology pipeline for the unsupervised integration of multi-omics data from Chronic Lymphocytic Leukemia (CLL) patients. Using Multi-Omics Factor Analysis (MOFA+), the analysis identifies latent biological drivers that are consistent across different molecular layers, including mRNA expression, DNA methylation, somatic mutations, and drug response.
​
A key highlight of this project is the use of the MEFISTO framework, which I applied to align multi-modal data and interpolate missing values based on structured covariates. The results demonstrate that the top latent factor (Factor 1) effectively captures the clinical heterogeneity associated with IGHV mutation status, explaining over 20% of the variance in the mutation and methylation views.​

# Project Overview
This repository contains a completed notebook pipeline for integrating multiple biological modalities into a shared latent factor space and interpreting those factors with variance decomposition, sample embeddings, and feature weights.
​It includes (1) MOFA on a CLL multi-omics dataset and (2) MEFISTO on an EvoDevo-style longitudinal multi-species dataset with smooth time covariates and optional dynamic time warping alignment.

# Pipeline consists of:
1. Loading and validating per-modality matrices into AnnData objects (e.g., mRNA, methylation, mutations, drugs) and reporting shapes (200 samples across all views).
​2. Building a unified MuData container and joining clinical metadata (Gender, age, TTT, TTD, treatedAfter, died, IGHV, trisomy12).
​3. Training MOFA models via mu.tl.mofa and saving .hdf5 outputs (e.g., modelsCLL.hdf5, modelsCLL5factors.hdf5).
​4. Interpreting factors using MOFAX: variance explained (plot_r2), factor scatter plots (e.g., Factor1 vs Factor2), and weights/loadings across views.
5. ​Comparing an original MOFA model (15 factors) against a constrained MOFA model (5 factors) by printing per-view R² tables for both.
​6. Running MEFISTO with smooth_covariate=time, groups_label=species, and smooth_warping=True (DTW alignment), then visualizing embeddings in factor space.
​7. Implementing a compact Bayesian factor model in Pyro (SVI with AutoDiagonalNormal + TraceELBO) and evaluating reconstruction with MSE.

# What I did & the purpose
As a course-style personal project, I used this notebook to practice multi-modal integration and learn how latent factor models can connect heterogeneous omics signals to biological/clinical annotations.
​Key learning goals were: selecting sensible likelihoods for different data types (e.g., Bernoulli for binary mutations), checking factor usefulness via variance explained, and linking sample embeddings to metadata such as IGHV and trisomy12.
​

# Tech stack
Multi-omics containers + modeling: scanpy, muon (MuData, mu.tl.mofa)
​Factor interpretation: mofax (plot_r2, factor plots, weight plots)
​Probabilistic modeling: torch, pyro (SVI / ELBO, AutoDiagonalNormal, TraceELBO)
​General: numpy, pandas, seaborn, matplotlib
​
​
## Core Implementation
The primary workflow used in this analysis is from multi-modal object creation to factor visualization.

```python
import muon as mu
import scanpy as sc
import mofax as mofa
import matplotlib.pyplot as plt

# 1. Initialize Multi-modal data object (MuData)
mdata = mu.read_h5mu("data/cll_data.h5mu")

# 2. Train MOFA+ Model
# We extract 15 latent factors to capture clinical and molecular heterogeneity
mu.tl.mofa(mdata, 
           n_factors=15, 
           convergence_mode="medium", 
           outfile="models/CLL_trained_model.hdf5")

# 3. Downstream Analysis with mofax
model = mofa.mofamodel("models/CLL_trained_model.hdf5")

# 4. Visualize Variance Explained
# This identifies which omic layer contributes most to each latent factor
mofa.plot_r2(model, vmax=15)
plt.savefig("plots/variance_explained.png")

# 5. Factor Latent Space
# Mapping samples by Factor 1 (IGHV-linked) and Factor 2
mofa.plot_factors(model, x="Factor1", y="Factor2", color="IGHV")
plt.show()

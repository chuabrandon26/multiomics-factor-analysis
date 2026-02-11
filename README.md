# Multiomics-factor-analysis
Unsupervised multi-omics integration of CLL data using MOFA+ and MEFISTO to identify latent drivers of disease variation.

# Multi-Modal Factor Analysis of CLL
## Description
This repository implements a systems biology pipeline for the unsupervised integration of multi-omics data from Chronic Lymphocytic Leukemia (CLL) patients. Using Multi-Omics Factor Analysis (MOFA+), the analysis identifies latent biological drivers that are consistent across different molecular layers, including mRNA expression, DNA methylation, somatic mutations, and drug response.
​
A key highlight of this project is the use of the MEFISTO framework, which I applied to align multi-modal data and interpolate missing values based on structured covariates. The results demonstrate that the top latent factor (Factor 1) effectively captures the clinical heterogeneity associated with IGHV mutation status, explaining over 20% of the variance in the mutation and methylation views.
​
## Core Implementation
The primary workflow used in this analysis, from multi-modal object creation to factor visualization.

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

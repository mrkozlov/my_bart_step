# Bayesian Additive Regression Trees (BART) for Species Distribution Modeling

## Description
This R script implements an advanced BART (Bayesian Additive Regression Trees) workflow for ecological modeling with comprehensive diagnostics, including:

- Three BART implementations (`embarcadero`, `BART`, `dbarts`)
- Spatial autocorrelation analysis using Moran's I
- Overfitting detection with multiple metrics
- Automated threshold selection for FNR control
- Variable importance assessment
- Extensive model comparison and validation

## Key Features
- **Multi-model comparison**: Compares predictions from different BART packages
- **Spatial diagnostics**: Moran's I analysis for residuals at optimal neighborhood sizes
- **Overfitting detection**: AUC/TSS differences, KS-tests, spatial pattern analysis
- **Threshold optimization**: Finds optimal thresholds for target false negative rates
- **Comprehensive outputs**: ROC curves, residual plots, probability distributions

## Installation
```r
install.packages(c("embarcadero", "dbarts", "pROC", "ggplot2", "spdep", "gridExtra", "reshape2"))


Usage
Prepare your data:

r
train_data <- read.csv("train_data.csv")  # Should contain response (nf) and predictors
test_data <- read.csv("test_data.csv")    # Should have same structure as train_data
Run the analysis:

r
source("mbsMorFac.R")
results <- my_bart_step(
  x.data = train_data,
  y.data = train_data$nf,
  test_data = test_data,
  selected_vars = c("bio1", "bio12", "elevation"), # Your predictors
  target_fnr = 0.04,  # Target false negative rate
  num_trees = 200,    # Number of trees
  do_spatial_analysis = TRUE  # Set FALSE if no coordinates available
)
Examine outputs:

bart_report.txt - Detailed textual output

*.pdf files - Diagnostic plots

Returned results object contains all metrics and models

Output Files
bart_step_diagnostic_plots.pdf - Model convergence diagnostics

residuals_histograms.pdf - Residual distributions

prediction_histograms.pdf - Probability distributions with thresholds

bart_report.txt - Comprehensive results summary

Parameters
Parameter	Description	Default
target_fnr	Target false negative rate	0.04
num_trees	Number of trees in ensemble	200
num_burn_in	Burn-in iterations	250
num_iterations_after_burn_in	Post-burn-in iterations	1000
tree.step	Tree pruning step	5
auc_weight	Weight for AUC in model selection	0.5
ks_weight	Weight for KS-test in model selection	0.5
Spatial Analysis
When coordinates (x,y) are provided, the script:

Calculates optimal neighborhood size (k) for Moran's I

Tests spatial autocorrelation in residuals

Performs cross-validation spatial analysis

Generates spatial overfitting warnings

License
MIT License - See LICENSE file for details



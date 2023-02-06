# LowRankGraphQuilting

This package implements the Low Rank Graphing Quilting procedures. 

# Requirements
- MATLAB 9.6+
- R 4.1.0+
  - `huge` package, 1.3.0+

# Usage

The primary analysis pipeline can be run using the code found in LRGQ_master.m. This script is divided in to two parts:
- Generation of synthetic data used for simulation studies.
- LRGQ estimation procedure, with covariance imputation in Matlab and graphical model selection in R. 
    - Requires: K x o table of features in each block, observed covariance matrix.

Links to sample calcium imaging data sets are found in LinksToPublicDataSets.txt.

# Behind the Matrix: Predicting Numerical Error in Generalized Eigenvalue Solvers

## Overview
This repository contains the code and data used in the paper 
"Behind the Matrix: Predicting Numerical Error in Generalized Eigenvalue Solvers" 
by Anna Diaz Farias.

## Repository Contents
- `matrix_generator.py` — Python script to generate synthetic SPD 
  matrix pairs and compute relative forward error. Produces 
  `matrix_dataset.csv`.
- `matrix_dataset.csv` — Full factorial dataset of 600 synthetic SPD 
  matrix pairs with computed relative forward error across combinations 
  of dimension n, scaling factor s, and condition number κ.
- `CC_proj.R` — R script for all statistical modeling, stepwise 
  selection, and figure generation.

## Requirements
### Python
- Python 3.9+
- NumPy
- SciPy
- pandas

### R
- R 4.2+
- tidyverse
- ggplot2
- scales
- emmeans

## Usage
1. Run `matrix_generator.py` to generate the dataset
2. Run `CC_proj.R` to reproduce all models and figures

## Contact
For questions, contact Anna Diaz at anna.diaz0702@gmail.com

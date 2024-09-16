# SPAR_GLM_Paper_Code
Reproducible Code to generate all results reported in 'Data-Driven Random Projection for Generalized Linear Models with High-Dimensional Data' by [Parzer, Filzmoser and Vana-Guer 2024](https://doi.org/10.48550/arXiv.2312.00130).

This repository consists of the following folders with described contents.

- data: .csv files for the data used in the paper (darwin, lymphoma and tribology)
- data_application: R-scripts for data applications applying all methods on all datasets (with train/test splits) multiple times and saving the resulting .rds file to the folder 'saved_results'; and another R-script applying SPAR to the full dataset once and visualizing results
- functions: 3 R-scripts, glm:data_generation.R for defining a function generating data from a HD generalized linear model, glm_methods.R defining consistent wrapper functions for all considered methods, multi_assign.R to define an operator assigning multiple variables at once (by Daniel Kapla, TU Wien)
- generate_plots: R-scipts reading in .rds files from 'saved_results' and generating the plots and tables for the simulation study and the data applications and saving the plots as pdfs in 'plots'
- plots: all pdf Figures
- saved_results: .rds files produced from 'simulation' or 'data_application' folders
- simulations: R-script for simulation study (both Src+RP and Benchmark) applying all methods on all simulation settings multiple times and saving the resulting .rds file to the folder 'saved_results'
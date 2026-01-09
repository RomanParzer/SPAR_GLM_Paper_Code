# SPAR_GLM_Paper_Code
Reproducible Code to generate all results reported in 'Data-Driven Random Projection and Screening for High-Dimensional Generalized Linear Models' by [Parzer, Filzmoser and Vana-Guer 2024](https://arxiv.org/abs/2410.00971).
(published version [Parzer, Filzmoser and Vana-Guer 2025](https://doi.org/10.1177/1471082X251392705))

This repository consists of the following folders with described contents.

- data: .csv files for the data used in the paper (darwin, lymphoma and tribology)
- data_application: R-scripts for data applications applying all methods on all datasets (with train/test splits) multiple times and saving the resulting .rds file to the folder 'saved_results'; and another R-script applying SPAR to the full dataset once and visualizing results
- functions: 3 R-scripts, glm:data_generation.R for defining a function generating data from a HD generalized linear model, glm_methods.R defining consistent wrapper functions for all considered methods, multi_assign.R to define an operator assigning multiple variables at once (by Daniel Kapla, TU Wien)
- generate_plots: R-scipts reading in .rds files from 'saved_results' and generating the plots and tables for the simulation study and the data applications and saving the plots as pdfs in 'plots'
- plots: all pdf Figures
- saved_results: .rds files produced from 'simulation' or 'data_application' folders
- simulations: R-script for simulation studies (RP with different m, Src+RP, and Benchmark) applying all methods on all simulation settings multiple times and saving the resulting .rds file to the folder 'saved_results'

The simulation studies and data applications were run on a local server with 160 available cores and around 400GB RAM with the following sessionInfo()-output:

#### R version 4.3.1 (2023-06-16 ucrt)

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows Server 2022 x64 (build 20348)

Matrix products: default


#### locale:

[1] LC_COLLATE=German_Austria.1252  LC_CTYPE=German_Austria.1252    LC_MONETARY=German_Austria.1252 LC_NUMERIC=C                    LC_TIME=German_Austria.1252    

time zone: Europe/Vienna

tzcode source: internal

#### attached base packages:

[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

#### other attached packages:

 [1] randomForest_4.7-1.1 e1071_1.7-14         metrica_2.0.3        ROCR_1.0-11          SPAR_3.2.0           MASS_7.3-60          SIS_0.8-8            glmnet_4.1-7        
 [9] Matrix_1.5-4.1       pls_2.8-2            dplyr_1.1.2          tidyr_1.3.1          foreach_1.5.2       

#### loaded via a namespace (and not attached):

 [1] minerva_1.5.10    gtable_0.3.6      shape_1.4.6.1     ggplot2_3.5.1     lattice_0.21-8    vctrs_0.6.3       tools_4.3.1       generics_0.1.3    stats4_4.3.1     
[10] proxy_0.4-27      tibble_3.2.1      fansi_1.0.4       RSQLite_2.3.5     DEoptimR_1.0-14   pacman_0.5.1      blob_1.2.4        pkgconfig_2.0.3   lifecycle_1.0.4  
[19] compiler_4.3.1    stringr_1.5.1     ncvreg_3.14.1     munsell_0.5.1     codetools_0.2-19  class_7.3-22      cellWise_2.5.3    pillar_1.9.0      rrcov_1.7-5      
[28] cachem_1.0.8      iterators_1.0.14  boot_1.3-28.1     robustbase_0.99-0 svd_0.5.5         tidyselect_1.2.1  mvtnorm_1.2-4     stringi_1.7.12    reshape2_1.4.4   
[37] purrr_1.0.2       splines_4.3.1     gsl_2.1-8         pcaPP_2.0-4       fastmap_1.1.1     grid_4.3.1        colorspace_2.1-0  cli_3.6.1         magrittr_2.0.3   
[46] survival_3.5-5    utf8_1.2.3        scales_1.3.0      bit64_4.0.5       energy_1.7-11     matrixStats_1.3.0 bit_4.0.5         gridExtra_2.3     memoise_2.0.1    
[55] doParallel_1.0.17 rlang_1.1.1       Rcpp_1.0.10       glue_1.6.2        DBI_1.2.1         R6_2.5.1          plyr_1.8.8 

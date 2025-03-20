# Scripts and results for "ConNIS..." by Hanke, Harten and Foraita (2025).

- [Scripts and results for "ConNIS..." by Hanke, Harten and Foraita (2025).](#scripts-and-results-for--connis--by-hanke--harten-and-foraita--2025-)
  * [Prerequisites](#prerequisites)
  * [Structure of the repository](#structure-of-the-repository)
    + [Subdirectories](#subdirectories)
    + [R scripts](#r-scripts)
  * [Run simulation study based on synthetic data](#run-simulation-study-based-on-synthetic-data)
  * [Run real world analyses](#run-real-world-analyses)
  * [Run semi-synthetic study](#run-semi-synthetic-study)
  * [Run instability approach](#run-instability-approach)
  * [Generate plots](#generate-plots)

`R` scripts and performance results based on

 * 160 synthetic data settings (for method comparison)
 * 4 semi-synthetic data settings (for method comparison)
 * 3 real world data settings (for method comparison)
 * 3 real world and 3 synthetic data examples (for evaluation of tuning parameter selection)

 Seeds have been used for all simulations and subsample drawings.

An interactive web app of all results is available under https://connis.bips.eu.

## Prerequisites

All simulations and analyses were run on a 64 core work station. `R >= 4.3.0` and the following packages are required:
* `parallel` (base `R`)
* `tidyverse` (CRAN)
* `MASS` (CRAN)
* `insdens` (https://github.com/Kevin-walters/insdens)
* `gmp` (CRAN)
* `ggpubr` (CRAN)
* `cowplot` (CRAN) 
* `ggridges` (CRAN)
* `ggh4x` (CRAN)

## Structure of the repository

```
ConNIS_results/
    ├──14028s_data/
    ├──MG1655_data/
    ├──bw25113_data/
    ├──data_for_synthetic_data_generation/
    ├──performance/
    ├──plots/
    ├──results/
    ├──simulatedData/
    ├──tmpData/
    dataSimulation.R
    example_<simu>.R
    functions.R
    generate_stability_tables.R
    <method>Analysis.R
    performanceAnalysis<Method>.R
    plot_<...>.R
    realworld_<strain>.R
    semi_synthetic.R
    simulations.R
    stabilities_<real_world_or_synthetic>.R
```

### Subdirectories

* `14028s_data/`, `MG1655_data/` and `bw25113_data/` contain real world data for the analyses; all data have been cloned by their publicly available original sources (see Hanke et al., 2025, for references)
* `data_for_synthetic_data_generation/` contains _E. coli_ data as reference for the synthetic data generation in `dataSimulation.R`
* `performance/` contains performances of all methods for synthetic, semi-synthetic and real world data analyses based on `performanceAnalysis<Method>.R`
* `plots/` contains plots generated with `plots_<...>.R` based on performances in `performance/`
* `results/` is used for analyses results of `<method>Anaysis.R`; is used by `performanceAnalysis<Method>.R`
* `simulatedData/` is used for generated synthetic data by `dataSimulation.R`; is used by `<method>Anaysis.R`
* `tmpData` saves intermediate results by the MCMC part of the _InsDens_ method

### R scripts

* `dataSimulation.R` generates synthetic data based on parameters in `simulation.R`
* `example_<simu>.R` generates results and performances for three synthetic data sets; is needed for the evaluation of the performance of the (in)stability approach by `stabilities_<real_world_or_synthetic>.R`
* `functions.R` contains all methods and functions (with exception of _InsDens_)
* `generate_stability_tables.R` generates a table with the performance of the (in)stability approach; uses results of `example_<simu>.R`, `realworld_<strain>.R` and `stabilities_<real_world_or_synthetic>.R`
* `<method>Anaysis.R`: applies one of the following methods to analyze the synthetic data: _Binomial_ (`binomial`, implementation of `TSAS 2.0`, Burger, 2017), _ConNIS_ (`connis`, Hanke et al., 2025), _Exp. vs. Gamma_ (`expvsgamma`, implementation of `Bio-TraDIS`, Barquist, 2016), _Geometric_ (`geometric`, Goodall, 2023), _InsDens_ (`insdens`, Nlebedim, 2021) and _Tn5Gaps_ (`tn5gaps`, implementation of `TRANSIT`, DeJesus, 2015); `<method>` can be set to 
* `performanceAnalysis<Method>.R`: gives the performances for `Binomial`, `ConNIS`, `ExpVsGamma`, `Geometric`, `InsDens` and `Tn5Gaps` based on results of `<method>Anaysis.R`
* `plot_<...>.R` generates plots for (semi-)synthetic and real world data analysis
* `realworld_<strain>.R` for either `14028S`, `bw25113` or `mg1655` analyzes real world data with all 6 methods; performances are saved in `performance/`
* `semi_synthetic.R` runs the full semi-synthetic data analysis on data from Goodall, 2018
* `simulations.R` main script to run the simulation study for synthetic data
* `stabilities_<real_world_or_synthetic>.R` runs the (in)stability approach for all 6 methods



## Run simulation study based on synthetic data

The default values for the simulation study are described in Hanke et al., 2025. To (re-)run the simulation study `simulations.R` needs to be called. It sets the parameters of the simulation study and the number of workers for the parallel computation using `parLapply` and than calls three types of scripts within loops over the different parameters:

* `dataSimulation.R` for generating the synthetic data
* `<method>Analysis.R` for the analysis of the data 
* `performance<method>.R` for the performance of the chosen `<method>`

:exclamation: While the number of workers is set in `dataSimulation.R`, the cluster type is set to `PSOCK` in all `<method>Analysis.R`. If you want to use a fork approach for parallelization (e.g. by `mclapply` or `makeForkCluster` scripts have to be modified individually).

## Run real world analyses
Run the scripts `realworld_<strain>.R`. Performances will be saved in `performance/`.

## Run semi-synthetic study 
Run the script `semi_synthetic.R`. Performances will be saved in `performance/`. 

:exclamation: Uses `mclapply` and at most `detectCores()-1` workers.

## Run instability approach
To calculate the instability values for each parameter/weight/threshold of the six methods run `stabilities_<real_world_or_synthetic>.R` for the real world data and the three synthetic datasets. Next, run `realworld_<strain>.R` (if you haven't done before) and `example_<simu>.R`. These scripts will give the results of three real world data and the three examples of synthetic data which will be used to evaluate the performance of the stability approach. Finally, run `generate_stability_tables.R` to generate an CSV file under `R` with the performances of the instability approach for all six methods and six datasets. 

## Generate plots
Simply run the script `plots_for_paper.R` to generate all plots under `R/plots`.
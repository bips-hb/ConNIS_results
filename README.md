# Scripts and results for "ConNIS..." 
(by Hanke, M., Harten, T. and Foraita, R. 2025)

- [Scripts and results for "ConNIS..."](#scripts-and-results-for--connis-)
  * [Dependencies](#dependencies)
  * [Structure of the repository](#structure-of-the-repository)
    + [Subdirectories](#subdirectories)
    + [R scripts](#r-scripts)
  * [Run simulation study based on synthetic data](#run-simulation-study-based-on-synthetic-data)
  * [Run real world analyses](#run-real-world-analyses)
  * [Run semi-synthetic study](#run-semi-synthetic-study)
  * [Run instability approach](#run-instability-approach)
  * [Generate plots](#generate-plots)
  * [References](#references)

`R` scripts and performance results are based on

 * 160 synthetic data settings (for method comparison)
 * 4 semi-synthetic data settings (for method comparison)
 * 3 real world data settings (for method comparison)
 * 3 real world and 3 synthetic data examples (for evaluation of tuning parameter selection)

 Seeds have been used for all simulations and subsample drawings.

An interactive web app with all results is available at https://connis.bips.eu.

For an implementation of ConNIS and the instability approach as an `R` package see https://github.com/bips-hb/ConNIS.

## Dependencies 

All simulations and analyses were run on a 64 core workstation. `R >= 4.3.0` and the following packages are required:
* `parallel` (base `R`)
* `tidyverse` (CRAN)
* `MASS` (CRAN)
* `insdens` (https://github.com/Kevin-walters/insdens)
* `gmp` (CRAN)
* `ggpubr` (CRAN)
* `cowplot` (CRAN) 
* `ggridges` (CRAN)
* `ggh4x` (CRAN)
* `readr` (CRAN)

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

The default values for the simulation study are described in Hanke et al., 2025. To (re-)run the simulation study, call `simulations.R`. It sets the parameters of the simulation study and the number of workers for the parallel computation using `parLapply`. It then calls three types of scripts via a loop structure over the different parameters:

* `dataSimulation.R` for generating the synthetic data
* `<method>Analysis.R` for the analysis of the data 
* `performance<method>.R` for the performance of the chosen `<method>`

(the heavy work is done in `<method>Analysis.R` by parallel computation)

:exclamation: While the number of workers is set in `dataSimulation.R`, the cluster type is set to `PSOCK` in all `<method>Analysis.R` scripts. If you want to use a fork approach for parallelization (e.g. `mclapply` or `makeForkCluster`), scripts have to be modified individually.

## Run real world analyses
Run the scripts `realworld_<strain>.R`. Performances will be saved in `performance/`.

## Run semi-synthetic study 
Run the script `semi_synthetic.R`. Performances will be saved in `performance/`. 

:exclamation: Uses `mclapply` and at most `detectCores()-1` workers.

## Run instability approach
First, calculate the instability values for each parameter/weight/threshold of the six methods based on `stabilities_<real_world_or_synthetic>.R` (for three real world data and three synthetic datasets). Next, run `realworld_<strain>.R` (if you haven't done before) and `example_<simu>.R` for results of the three real world data and the three examples of synthetic data. These are used to evaluate the performance of the stability approach. Finally, run `generate_stability_tables.R` to generate an CSV file under `R/` with the performances of the instability approach. 

:exclamation: Uses `PSOCK` for parallelization.

## Generate plots
Simply run the script `plots_for_paper.R` to generate all plots under `R/plots`.

## References
* `Barquist, L. et al. The tradis toolkit: sequencing and analysis for dense transposon mutant libraries. Bioinformatics 32, 1109–1111 (2016). URL http://dx.doi.org/ 10.1093/bioinformatics/btw022`
* `Burger, B. T., Imam, S., Scarborough, M. J., Noguera, D. R. & Donohue, T. J. Combining genome-scale experimental and computational methods to identify essential genes in rhodobacter sphaeroides. mSystems 2 (2017). URL http://dx. doi.org/10.1128/msystems.00015-17`
* `DeJesus, M. A., Ambadipudi, C., Baker, R., Sassetti, C. & Ioerger, T. R. Transit - a software tool for himar1 tnseq analysis. PLOS Computational Biology 11, e1004401 (2015). URL http://dx.doi.org/10.1371/journal.pcbi.1004401`
* `Goodall, E. C. A. et al. The essential genome of escherichia coli k-12. mBio 9 (2018). URL http://dx.doi.org/10.1128/mBio.02096-17.`
* `Goodall, E. C. A. et al. A multiomic approach to defining the essential genome of the globally important pathogen corynebacterium diphtheriae. PLOS Genetics 19, e1010737 (2023). URL http://dx.doi.org/10.1371/journal.pgen.1010737`
* `Nlebedim, V. U., Chaudhuri, R. R. & Walters, K. Probabilistic identification of bacterial essential genes via insertion density using tradis data with tn5 libraries. Bioinformatics 37, 4343–4349 (2021). URL http://dx.doi.org/10.1093/ bioinformatics/btab508`

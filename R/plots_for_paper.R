# Plots for paper

library(tidyverse)
library(ggh4x)
library(cowplot)
library(reshape2)
library(ggpubr)
library(ggvenn)
library(venn)
library(ggridges)

setwd("./R")

# Plot real word data Ecoli, Salmonella enterica Serovar Typhimurium,
rm(list=ls())
source("plot_real_data.R")

# Plot real word data Synpcc (Venn diagramm)
# rm(list=ls())
# source("plot_real_data_venn_synpcc7942.R")

# Plot examples of gene-wise insertion densities for real and simulated data
rm(list=ls())
source("plot_real_vs_simulated_data.R")

# plot for synthetic settings
rm(list=ls())
source("plot_syn.R")

# plot for effect of weighting
rm(list=ls())
source("plot_effect_weighting.R")




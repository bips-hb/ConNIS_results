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

# Plot performances real word data Ecoli, Salmonella enterica Serovar Typhimurium,
rm(list=ls())
source("plot_real_data.R")

# plot performances for synthetic settings
rm(list=ls())
source("plot_syn.R")

# plot for effect of weighting
rm(list=ls())
source("plot_effect_weighting.R")

# Plot examples of gene-wise insertion densities for real and simulated data
rm(list=ls())
source("plot_real_vs_simulated_data.R")


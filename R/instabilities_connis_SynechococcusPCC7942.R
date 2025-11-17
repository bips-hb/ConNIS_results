library(tidyverse)
library(gmp)
library(ConNIS)

# Data from https://www.pnas.org/doi/full/10.1073/pnas.1519220112

setwd("./R")

data_t0 <- readxl::read_xlsx("./PCC7942/pnas.1519220112.sd01.xlsx")
IS_pos_t0 <- sort(unique(data_t0$position))
data_t7 <- readxl::read_xlsx("./PCC7942/pnas.1519220112.sd02.xlsx")
IS_pos_t7 <- sort(unique(data_t7$position))

gene_list <- readxl::read_xlsx("./PCC7942/pnas.1519220112.sd03.xlsx")
gene_list <- gene_list %>% filter(chromosome == "chromosome")


# connis_instabilities <- ConNIS::instabilities(
#   ins.positions = IS_pos_t0, 
#   gene.names = gene_list$`7942_ID`, 
#   method = "ConNIS", 
#   gene.starts = gene_list$start, 
#   gene.stops = gene_list$end, 
#   genome.length = max(gene_list$end), 
#   d = 0.5,
#   m = 500,
#   weights = seq(0.05, 1, 0.05),  
#   use.parallelization = TRUE, 
#   parallelization.type = "mclapply", 
#   numCores = 50,
#   seed = 1,
#   set.rng = "L'Ecuyer-CMRG"
#   )
# 
# saveRDS(connis_instabilities, file = "results/connis_instabilities_t0_SynechococcusPCC7942.RDS")


# load stability values from t0 and determine w
connis_instabilities <- readRDS(
  "results/connis_instabilities_SynechococcusPCC7942.RDS"
)

selected_w <- connis_instabilities$weight_value[
  which.min(connis_instabilities$instability)
]

results_t0 <- ConNIS(
  ins.positions = IS_pos_t0, 
  gene.names    = gene_list$`7942_ID`, 
  gene.starts   = gene_list$start, 
  gene.stops    = gene_list$end, 
  genome.length = max(gene_list$end), 
  weight        = selected_w
)

alpha_bonf_t0 <- 0.05 / nrow(gene_list)

essential_t0 <- gene_list$`7942_ID`[results_t0$p_value <= alpha_bonf_t0]

gene_list_no_t0 <- gene_list[!(gene_list$`7942_ID` %in% essential_t0), ]

connis_instabilities <- ConNIS::instabilities(
  ins.positions = IS_pos_t7, 
  gene.names = gene_list_no_t0$`7942_ID`, 
  method = "ConNIS", 
  gene.starts = gene_list_no_t0$start, 
  gene.stops = gene_list_no_t0$end, 
  genome.length = max(gene_list_no_t0$end), 
  d = 0.5,
  m = 500,
  weights = seq(0.05, 1, 0.05),  
  use.parallelization = TRUE, 
  parallelization.type = "mclapply", 
  numCores = 50,
  seed = 1,
  set.rng = "L'Ecuyer-CMRG"
)

saveRDS(connis_instabilities, file = "results/connis_instabilities_t7_SynechococcusPCC7942.RDS")




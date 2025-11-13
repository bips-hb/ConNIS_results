library(tidyverse)
library(gmp)
library(ConNIS)

# Data from https://www.pnas.org/doi/full/10.1073/pnas.1519220112

setwd("./R")

data_t0 <- readxl::read_xlsx("./PCC7942/pnas.1519220112.sd01.xlsx")
IS_pos_t0 <- sort(unique(data_t0$position))
data_t6 <- readxl::read_xlsx("./PCC7942/pnas.1519220112.sd02.xlsx")
IS_pos_t6 <- sort(unique(data_t6$position))

gene_list <- readxl::read_xlsx("./PCC7942/pnas.1519220112.sd03.xlsx")
gene_list <- gene_list %>% filter(chromosome == "chromosome")



library(ConNIS)

ConNIS(ins.positions = data_t0$position, 
       gene.names = gene_list$`7942_ID`, 
       gene.starts = gene_list$start, 
       gene.stops = gene_list$end, 
       genome.length = max(gene_list$end)+1, num.ins.per.gene = rep(1,nrow(gene_list)))

connis_instabilities <- ConNIS::instabilities(
  ins.positions = IS_pos_t0, 
  gene.names = gene_list$`7942_ID`, 
  method = "ConNIS", 
  gene.starts = gene_list$start, 
  gene.stops = gene_list$end, 
  genome.length = max(gene_list$end), 
  d = 0.5,
  m = 500,
  weights = seq(0.05, 1, 0.05),  
  use.parallelization = TRUE, 
  parallelization.type = "mclapply", 
  numCores = 50,
  seed = 1,
  set.rng = "L'Ecuyer-CMRG"
  )

saveRDS(connis_instabilities, file = "results/connis_instabilities_SynechococcusPCC7942.RDS")



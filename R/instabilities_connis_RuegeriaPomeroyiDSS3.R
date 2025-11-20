library(tidyverse)
library(gmp)
library(ConNIS)

# Ruegeria pomeroyi DSS-3
# Data fromhttps://www.pnas.org/doi/10.1073/pnas.2217200120#sec-1

setwd("./R")

IS_raw <- readxl::read_xlsx("./DSS3/pnas.2217200120.sd02.xlsx", skip = 2)

IS_pos_t0 <- IS_raw %>% mutate(readCount = InitialLibrary_1 + 
                                 InitialLibrary_2 +
                                 InitialLibrary_3 +
                                 InitialLibrary_4) %>%
  filter(readCount>0) %>% 
  select(genomePosition)

IS_pos_t8 <- IS_raw %>% mutate(readCount = Rp_1 + 
                                 Rp_2 +
                                 Rp_3 +
                                 Rp_4) %>%
  filter(readCount>0) %>% 
  select(genomePosition)


gene_list <- read_tsv("./DSS3//R_pomeroyi_genes.tsv")


connis_instabilities <- ConNIS::instabilities(
  ins.positions = IS_pos_t0$genomePosition,
  gene.names = gene_list$locus_tag,
  method = "ConNIS",
  gene.starts = gene_list$start,
  gene.stops = gene_list$end,
  genome.length = 4109437,
  d = 0.5,
  m = 500,
  weights = seq(0.05, 1, 0.05),
  use.parallelization = TRUE,
  parallelization.type = "mclapply",
  numCores = 50,
  seed = 2025,
  set.rng = "L'Ecuyer-CMRG"
  )

saveRDS(connis_instabilities, file = "results/connis_instabilities_t0_RuegeriaPomeroyiDSS3.RDS")


# load stability values from t0 and determine w
connis_instabilities <- readRDS(
  "results/connis_instabilities_t0_RuegeriaPomeroyiDSS3.RDS"
)

selected_w <- connis_instabilities$weight_value[
  which.min(connis_instabilities$instability)
]

results_t0 <- ConNIS(
  ins.positions = IS_pos_t0$genomePosition, 
  gene.names    = gene_list$locus_tag, 
  gene.starts   = gene_list$start, 
  gene.stops    = gene_list$end, 
  genome.length = 4109437, 
  weight        = selected_w
)

alpha_bonf_t0 <- 0.05 / nrow(gene_list)

essential_t0 <- gene_list$locus_tag[results_t0$p_value <= alpha_bonf_t0]

gene_list_no_t0 <- gene_list[!(gene_list$locus_tag %in% essential_t0), ]


connis_instabilities <- ConNIS::instabilities(
  ins.positions = IS_pos_t8$genomePosition, 
  gene.names = gene_list_no_t0$locus_tag, 
  method = "ConNIS", 
  gene.starts = gene_list_no_t0$start, 
  gene.stops = gene_list_no_t0$end, 
  genome.length = 4109437, 
  d = 0.5,
  m = 500,
  weights = seq(0.05, 1, 0.05),  
  use.parallelization = TRUE, 
  parallelization.type = "mclapply", 
  numCores = 50,
  seed = 1,
  set.rng = "L'Ecuyer-CMRG"
)

saveRDS(connis_instabilities, file = "results/connis_instabilities_t8_RuegeriaPomeroyiDSS3.RDS")




library(tidyverse)
library(gmp)
library(ConNIS)

# Ruegeria pomeroyi DSS-3
# Data fromhttps://www.pnas.org/doi/10.1073/pnas.2217200120#sec-1

setwd("./R")

IS_raw <- readxl::read_xlsx("./DSS3/pnas.2217200120.sd02.xlsx", skip = 2)

IS_pos_initial <- IS_raw %>% mutate(readCount = InitialLibrary_1 + 
                                      InitialLibrary_2 +
                                      InitialLibrary_3 +
                                      InitialLibrary_4) %>%
  filter(readCount > 0) %>% 
  select(genomePosition)

IS_pos_RP <- IS_raw %>% mutate(readCount = Rp_1 + 
                                 Rp_2 +
                                 Rp_3 +
                                 Rp_4) %>%
  filter(readCount > 0 &
           InitialLibrary_1 + 
           InitialLibrary_2 +
           InitialLibrary_3 +
           InitialLibrary_4 > 0) %>% 
  select(genomePosition)

IS_pos_RPplusV <- IS_raw %>% mutate(readCount = `Rp+V_1` + 
                                      `Rp+V_2` +
                                      `Rp+V_3` +
                                      `Rp+V_4`) %>%
  filter(readCount > 0 &
           InitialLibrary_1 + 
           InitialLibrary_2 +
           InitialLibrary_3 +
           InitialLibrary_4 > 0) %>% 
  select(genomePosition)

IS_pos_RPplusM <- IS_raw %>% mutate(readCount = `Rp+M_1` + 
                                      `Rp+M_2` +
                                      `Rp+M_3` +
                                      `Rp+M_4`) %>%
  filter(readCount > 0 &
           InitialLibrary_1 + 
           InitialLibrary_2 +
           InitialLibrary_3 +
           InitialLibrary_4 > 0) %>% 
  select(genomePosition)

IS_pos_RPplusVandM <- IS_raw %>% mutate(readCount = `Rp+V+M_1` + 
                                          `Rp+V+M_2` +
                                          `Rp+V+M_3` +
                                          `Rp+V+M_4`) %>%
  filter(readCount > 0 &
           InitialLibrary_1 + 
           InitialLibrary_2 +
           InitialLibrary_3 +
           InitialLibrary_4 > 0) %>% 
  select(genomePosition)


gene_list <- read_tsv("./DSS3//R_pomeroyi_genes.tsv")

fitness_data_original <- readxl::read_xlsx("./DSS3/pnas.2217200120.sd03.xlsx", skip = 3 )
for(i in 5:21){
  fitness_data_original[,i] <- 
    as.numeric(unlist(fitness_data_original[,i]))
}
rm(i)

# connis_instabilities <- ConNIS::instabilities(
#   ins.positions = IS_pos_initial$genomePosition,
#   gene.names = gene_list$locus_tag,
#   method = "ConNIS",
#   gene.starts = gene_list$start,
#   gene.stops = gene_list$end,
#   genome.length = 4109437,
#   d = 0.5,
#   m = 500,
#   weights = seq(0.05, 1, 0.05),
#   use.parallelization = TRUE,
#   parallelization.type = "mclapply",
#   numCores = 50,
#   seed = 2025,
#   set.rng = "L'Ecuyer-CMRG"
# )
# 
# saveRDS(connis_instabilities, file = "results/connis_instabilities_initial_RuegeriaPomeroyiDSS3.RDS")


# load stability values from initial and determine w
# connis_instabilities_initial <- readRDS(
#   "results/connis_instabilities_initial_RuegeriaPomeroyiDSS3.RDS"
# )
# 
# selected_w_initial <- connis_instabilities_initial$weight_value[
#   which.min(connis_instabilities_initial$instability)
# ]
# 
# results_initial <- ConNIS(
#   ins.positions = IS_pos_initial$genomePosition, 
#   gene.names    = gene_list$locus_tag, 
#   gene.starts   = gene_list$start, 
#   gene.stops    = gene_list$end, 
#   genome.length = 4109437, 
#   weight        = selected_w_initial
# )
# 
# alpha_bonf_initial <- 0.05 / nrow(gene_list)
# 
# essential_initial <- gene_list$locus_tag[results_initial$p_value <= alpha_bonf_initial]



# connis_instabilities <- ConNIS::instabilities(
#   ins.positions = IS_pos_RP$genomePosition,
#   gene.names = gene_list$locus_tag,
#   method = "ConNIS",
#   gene.starts = gene_list$start,
#   gene.stops = gene_list$end,
#   genome.length = 4109437,
#   d = 0.5,
#   m = 500,
#   weights = seq(0.05, 1, 0.05),
#   use.parallelization = TRUE,
#   parallelization.type = "mclapply",
#   numCores = 50,
#   seed = 2025,
#   set.rng = "L'Ecuyer-CMRG"
# )
# 
# saveRDS(connis_instabilities, file = "results/connis_instabilities_RP_RuegeriaPomeroyiDSS3.RDS")


# load stability values from initial and determine w
# connis_instabilities_RP <- readRDS(
#   "results/connis_instabilities_RP_RuegeriaPomeroyiDSS3.RDS"
# )
# 
# selected_w_RP <- connis_instabilities_RP$weight_value[
#   which.min(connis_instabilities_RP$instability)
# ]
# 
# results_RP <- ConNIS(
#   ins.positions = IS_pos_RP$genomePosition, 
#   gene.names    = gene_list$locus_tag, 
#   gene.starts   = gene_list$start, 
#   gene.stops    = gene_list$end, 
#   genome.length = 4109437, 
#   weight        = selected_w_RP
# )
# 
# alpha_bonf_RP <- 0.05 / nrow(gene_list)
# 
# essential_RP <- gene_list$locus_tag[results_RP$p_value <= alpha_bonf_RP]



connis_instabilities <- ConNIS::instabilities(
  ins.positions = IS_pos_RPplusV$genomePosition,
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

saveRDS(connis_instabilities, file = "results/connis_instabilities_RPplusV_RuegeriaPomeroyiDSS3.RDS")


# load stability values from initial and determine w
# connis_instabilities_RPplusV <- readRDS(
#   "results/connis_instabilities_RPplusV_RuegeriaPomeroyiDSS3.RDS"
# )
# 
# selected_w_RPplusV <- connis_instabilities_RPplusV$weight_value[
#   which.min(connis_instabilities_RPplusV$instability)
# ]
# 
# results_RPplusV <- ConNIS(
#   ins.positions = IS_pos_RPplusV$genomePosition, 
#   gene.names    = gene_list$locus_tag, 
#   gene.starts   = gene_list$start, 
#   gene.stops    = gene_list$end, 
#   genome.length = 4109437, 
#   weight        = selected_w_RPplusV
# )
# 
# alpha_bonf_RPplusV <- 0.05 / nrow(gene_list)
# 
# essential_RPplusV <- gene_list$locus_tag[results_RPplusV$p_value <= alpha_bonf_RPplusV]


connis_instabilities <- ConNIS::instabilities(
  ins.positions = IS_pos_RPplusM$genomePosition,
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

saveRDS(connis_instabilities, file = "results/connis_instabilities_RPplusM_RuegeriaPomeroyiDSS3.RDS")


# load stability values from initial and determine w
# connis_instabilities_RPplusM <- readRDS(
#   "results/connis_instabilities_RPplusM_RuegeriaPomeroyiDSS3.RDS"
# )
# 
# selected_w_RPplusM <- connis_instabilities_RPplusM$weight_value[
#   which.min(connis_instabilities_RPplusM$instability)
# ]
# 
# results_RPplusM<- ConNIS(
#   ins.positions = IS_pos_RPplusM$genomePosition, 
#   gene.names    = gene_list$locus_tag, 
#   gene.starts   = gene_list$start, 
#   gene.stops    = gene_list$end, 
#   genome.length = 4109437, 
#   weight        = selected_w_RPplusM
# )
# 
# alpha_bonf_RPplusM <- 0.05 / nrow(gene_list)
# 
# essential_RPplusM <- gene_list$locus_tag[results_RPplusM$p_value <= alpha_bonf_RPplusM]


connis_instabilities <- ConNIS::instabilities(
  ins.positions = IS_pos_RPplusVandM$genomePosition,
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

saveRDS(connis_instabilities, file = "results/connis_instabilities_RPplusVandM_RuegeriaPomeroyiDSS3.RDS")


# load stability values from initial and determine w
# connis_instabilities_RPplusVandM <- readRDS(
#   "results/connis_instabilities_RPplusVandM_RuegeriaPomeroyiDSS3.RDS"
# )
# 
# selected_w_RPplusVandM <- connis_instabilities_RPplusVandM$weight_value[
#   which.min(connis_instabilities_RPplusVandM$instability)
# ]
# 
# results_RPplusVandM <- ConNIS(
#   ins.positions = IS_pos_RPplusVandM$genomePosition, 
#   gene.names    = gene_list$locus_tag, 
#   gene.starts   = gene_list$start, 
#   gene.stops    = gene_list$end, 
#   genome.length = 4109437, 
#   weight        = selected_w
# )
# 
# alpha_bonf_RPplusVandM <- 0.05 / nrow(gene_list)
# 
# essential_RPplusVandM <- gene_list$locus_tag[results_RPplusVandM$p_value <= alpha_bonf_RPplusVandM]

# 
# fitness_data_connis <-tibble(
#   Gene = gene_list$locus_tag,
#   p_initial = results_initial$p_value,
#   p_RP = results_RP$p_value,
#   w_RP = log(p_RP+10^(-100),10) - log(p_initial+10^(-100), 10),
#   p_RPplusV = results_RPplusV$p_value,
#   w_RPplusV = log(p_RPplusV+10^(-100), 10) - log(p_initial+10^(-100), 10),
#   p_RPplusM = results_RPplusM$p_value,
#   w_RPplusM = log(p_RPplusM+10^(-100), 10) - log(p_initial+10^(-100), 10),
#   p_RPplusVandM = results_RPplusVandM$p_value,
#   w_RPplusVandM = log(p_RPplusVandM+10^(-100), 10) - log(p_initial+10^(-100), 10),
# ) 
# 
# ggplot(fitness_data_connis,
#        aes(x=log(-log(p_initial,10)), y=log(-log(p_RP,10)))) +
#   geom_point() +
#   geom_hline(yintercept = log(-log(0.05/nrow(fitness_data_connis),10)))+
#   geom_vline(xintercept = log(-log(0.05/nrow(fitness_data_connis),10)))
# 
# 
# ggplot(fitness_data_connis,
#        aes(x=-log(p_initial,10), y=-log(p_RP,10))) +
#   geom_point() +
#   geom_hline(yintercept = -log(0.05/nrow(fitness_data_connis),10))+
#   geom_vline(xintercept = -log(0.05/nrow(fitness_data_connis),10))
# 
# 
# genes_for_comaprison <- 
#   intersect(fitness_data_connis$Gene, 
#             fitness_data_original$Gene)
# 
# fitness_data_connis_reduced <- 
#   fitness_data_connis %>% filter(Gene %in% genes_for_comaprison)
# 
# fitness_data_original_reduced <- 
#   fitness_data_original %>% filter(Gene %in% genes_for_comaprison)
# 
# data_comparison <- 
#   right_join(x = fitness_data_connis_reduced, y = fitness_data_original_reduced)
# 
# data_comparison <- data_comparison[-which(sapply(1:nrow(data_comparison), function(i){
#   any(is.na(data_comparison[i,]))
# })),]
# 
# cor(data_comparison$w_RP, 
#     data_comparison$mean_W_Rp, 
#          use = "complete.obs")
# 
# cor(data_comparison$w_RPplusV, 
#    data_comparison$`mean_W_Rp+V`, 
#     use = "complete.obs")
# 
# cor(data_comparison$w_RPplusM, 
#    data_comparison$`mean_W_Rp+M`, 
#     use = "complete.obs")
# 
# cor(data_comparison$w_RPplusVandM, 
#     data_comparison$`mean_W_Rp+V+M`, 
#     use = "complete.obs")
# 
# 
# 
# cor(data_comparison$w_RPplusV-data_comparison$w_RP, 
#     data_comparison$`log2_foldchange_Rp+V`, 
#     use = "complete.obs")
# 
# cor(data_comparison$w_RPplusM-data_comparison$w_RP, 
#     data_comparison$`log2_foldchange_Rp+M`, 
#     use = "complete.obs")
# 
# cor(data_comparison$w_RPplusVandM-data_comparison$w_RP, 
#     data_comparison$`log2_foldchange_Rp+V+M`, 
#     use = "complete.obs")
# 
# sum(
# which(rank(data_comparison$w_RPplusV-data_comparison$w_RP) <= 20) %in%
# which(rank(data_comparison$`log2_foldchange_Rp+M`) <= 20)
# )
# plot(sort(data_comparison$w_RPplusM/data_comparison$w_RP))

library(insdens)
library(tidyverse)
library(MASS)
library(gmp)

setwd("./R")

# Load functions
source("./functions.R")

# Load gene Data of the of the BW25113 substrain; Retrieved on BioCyc
# https://biocyc.org/group?id=:ALL-GENES&orgid=GCF_000750555
# Created: 02-Mar-2025 14:09:41
gene_data <- read_tsv("./bw25113_data/All-genes-of-E.-coli-K-12-substr.-BW25113.txt")

# rename, sort by start and drop Procut and Accession-1
gene_data <- gene_data[, c(1,3,4)]
names(gene_data) <- c("gene", "start", "end")
gene_data <- gene_data %>% arrange(start)

# Load final IS from Goodall et al. (2018),
is_count_neg <-
  read.table("./bw25113_data/Ecoli_BW25113_chlor-tn_position_depth_count_neg.txt")
is_count_pos <-
  read.table("./bw25113_data/Ecoli_BW25113_chlor-tn_position_depth_count_pos.txt")

# Combine + & - (in total 884,675 IS)
is_pos <- which(is_count_neg+is_count_pos > 0)

# get number of IS per gene to gene_data
num_ins_per_gene_full <- sapply(1:nrow(gene_data), function(i) {
  sum(is_pos >= gene_data$start[i] &
        is_pos <= gene_data$end[i])
  })

# add number of Is per gene based on all IS to data
gene_data$num_IS_FULL <- num_ins_per_gene_full

# Similar to Goodad (2018) use the the Exp vs Gamma based on ALL IS as reference,
ExpVsGamma_results <-
  ExpVsGamma(gene.names = gene_data$gene,
             gene.starts = gene_data$start,
             gene.stops = gene_data$end,
             num.ins.per.gene = gene_data$num_IS_FULL,
             log2threshold = 12)

# Based on log2 threshold =12 determine the "essential" genes (based on the
# high satrated library of Goddall, 2018). This will be the "truth" to evaluate
# how the other methods perform if a subsample of the IS are drawn
ref_ess_gene <-
  ExpVsGamma_results$gene[which(ExpVsGamma_results$essential==1)]

# numbe rof "true" essential genes
length(ref_ess_gene)
# save reference essential genes
saveRDS(ref_ess_gene,"./bw25113_data/ref_ess_semisyn.RDS")

# Draw a subsamples: 200,000 and 100,000 of the total number of IS
subsample_sizes <- c(50000, 100000, 200000, 400000)

# Run Binomial, ConNIS, InsDens Geometric and Tn5Gaps on the subsample IS
# use different weights and psoterior probability thresholds r
weights <- seq(0.1, 1, 0.1)
post.prob.thresholds <- c(seq(0.1, 0.9, 0.1), 0.99)

# start lapply loop over subsample sizes
subsample_results <- lapply(subsample_sizes, function(subsample_size){
  print(paste("Running setting with subsample of size", subsample_size))

  set.seed(subsample_size)
  subsample_is_pos <- sort(sample(is_pos, subsample_size, replace = FALSE))

  # get number of IS per gene based on the subasmple
  num_ins_per_gene_subsample <- sapply(1:nrow(gene_data), function(i) {
    sum(subsample_is_pos >= gene_data$start[i] &
          subsample_is_pos <= gene_data$end[i])
  })

  # add number of Is per gene based on subsample IS to data
  gene_data$num_IS_SUBSAMPLE <- num_ins_per_gene_subsample


  print(paste("Start Binomial at", Sys.time() ))
  results_Binomial <- mclapply(X = weights,mc.cores = detectCores(), FUN = function(w){

    out <- Binomial(ins.positions = subsample_is_pos,
                    gene.names = gene_data$gene,
                    gene.starts = gene_data$start,
                    gene.stops = gene_data$end,
                    num.ins.per.gene = gene_data$num_IS_SUBSAMPLE,
                    genome.length = 4631469,
                    weight = w)


    out
  })

  print(paste("Start ConNIS at", Sys.time() ))
  results_ConNIS <- mclapply(X = weights,mc.cores = detectCores(), FUN = function(w){

    out <- ConNIS(ins.positions = subsample_is_pos,
                  gene.names = gene_data$gene,
                  gene.starts = gene_data$start,
                  gene.stops = gene_data$end,
                  genome.length = 4631469,
                  num.ins.per.gene = gene_data$num_IS_SUBSAMPLE,
                  weight = w)


    out
  })

  print(paste("Start Geometric at", Sys.time() ))
  results_Geometric <- mclapply(X = weights,mc.cores = detectCores(), FUN = function(w){

    out <- Geometric(ins.positions = subsample_is_pos,
                     gene.names = gene_data$gene,
                     gene.starts = gene_data$start,
                     gene.stops = gene_data$end,
                     num.ins.per.gene = gene_data$num_IS_SUBSAMPLE,
                     genome.length = 4631469,
                     weight = w)


    out
  })



  print(paste("Start Tn5Gaps at", Sys.time() ))
  results_Tn5Gaps <- mclapply(X = weights,mc.cores = detectCores(), FUN = function(w){

    out <- Tn5Gaps(ins.positions = subsample_is_pos,
                   gene.names = gene_data$gene,
                   gene.starts = gene_data$start,
                   gene.stops = gene_data$end,
                   genome.length = 4631469,
                   weight = w)


    out
  })

  list(results_Binomial = results_Binomial,
       results_ConNIS = results_ConNIS,
       results_Geometric = results_Geometric,
       results_Tn5Gaps = results_Tn5Gaps)

})
# save results
saveRDS(subsample_results, file="./results/Results_semisyn_BW25113.RDS")


# extract for each subsample size(i), each method (j) and each tuning
# value (k) the classification performance
semi_syn_performances <- lapply(seq_along(subsample_sizes), function(i){

  data_subsample_size <- subsample_results[[i]]
  subsample_size <- subsample_sizes[i]

  performances_subsample <- lapply(seq_along(data_subsample_size), function(j){

    data_subsample_size_method <- data_subsample_size[[j]]

    if(names(data_subsample_size[j]) %in% c("results_Binomial",
                                            "results_ConNIS",
                                            "results_Geometric",
                                            "results_Tn5Gaps")){
      method <- str_match(names(data_subsample_size[j]), pattern = "_(.*)")[2]

      out <- lapply(data_subsample_size_method, function(k){

        est_ess_genes <- k$gene[
          p.adjust(k$p_value, method = "hochberg") <= 0.05]

        out <- classification_performance(
          as.numeric(gene_data$gene %in% est_ess_genes),
          as.numeric(gene_data$gene %in% ref_ess_gene))

        out$tuning <- as.numeric(k[1,3])
        out$method<- method
        out$subsample_size <- subsample_size
        out

      })

      out <- do.call(rbind, out)
      out



    }else if(names(data_subsample_size[j]) == "results_InsDens"){

      method <- "InsDens"

      out <- lapply(data_subsample_size_method, function(k){

        est_ess_genes <- k$gene[k$essential == 1 ]

        out <- classification_performance(
          as.numeric(gene_data$gene %in% est_ess_genes),
          as.numeric(gene_data$gene %in% ref_ess_gene))

        out$tuning <- as.numeric(k[1,3])
        out$method<- method
        out$subsample_size <- subsample_size
        out

      })

      out <- do.call(rbind, out)
      out

    }

  })

  do.call(rbind, performances_subsample)

})


semi_syn_performances <- do.call(rbind, semi_syn_performances)

saveRDS(semi_syn_performances, file="./performance/Performance_semisyn_BW25113.RDS")



library(tidyverse)
library(insdens)
library(MASS)
library(gmp)

setwd("./R")

source("functions.R")

### Simulation data: 400,000 IS, 2% noise, lambda=0.8, 0 hot & cold spots 
simulated_data <- readRDS("./simulatedData/SimuData_Tnseq_GeneClusters_ess_ORF_uniform_numLoci_400000_lambda_0.8_num_hot_spots_0_hot_spot_size_10000_num_cold_spots_0_cold_spot_size_10000_NBnc_100_NBdis_1_NBp_0.3_tecnoise_0.02.RDS")

set.seed(5294)
random_sim_run <- sample(1:100, 1)
tnseq_data <- simulated_data[[random_sim_run]]$tnseqData[, c(1,4,5)]
names(tnseq_data)<- c("gene", "gene_start", "gene_end")

all_IS <- tibble(pos=simulated_data[[random_sim_run]]$observed_all_unique_loci)


original_data <- tnseq_data

true_ess <- as.numeric(simulated_data[[random_sim_run]]$tnseqData$true_essentiality == "essential")


#### trimming
list_trimming_start <- c(0.05)
list_trimming_end <- c(0.05)

alpha_value <- 0.05
ws <- c(0.01, seq(0.05,1,0.05))
ts <- c(seq(0.1,2, 0.1),3:12)
rs <- c(seq(0.01,0.1, 0.01), seq(0.2,0.9,0.1), seq(0.91,0.99,0.01))

performances <- lapply(seq_along(list_trimming_start), function(n){
  
  trimming_start <- list_trimming_start[n]
  trimming_end <- list_trimming_end[n]
  
  gene_data <- lapply(1:nrow(original_data), function(i){
    
    gene_i_start <- original_data$gene_start[i] + 
      floor(trimming_start * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    gene_i_end <- original_data$gene_end[i] - 
      floor(trimming_end * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    IS_gene_i <- 
      all_IS$pos[which(gene_i_start <= all_IS$pos &
                         gene_i_end >= all_IS$pos)]
    
    num_ins <- length(IS_gene_i)
    
    gene <- original_data$gene[i]
    if(is.null(gene)){
      gene <- NA
    }
    
    
    if(gene_i_start  < min(IS_gene_i) |
       length(IS_gene_i)==0){
      start_is <- gene_i_start-1
    }else{
      start_is <- min(IS_gene_i)
    }
    
    if(gene_i_end > max(IS_gene_i) |
       length(IS_gene_i)==0){
      end_is <- gene_i_end[i]+1
    }else{
      end_is <- max(IS_gene_i)
    }
    
    tibble(entry = original_data$locusId[i],
           gene = gene,
           gene_start = gene_i_start,
           gene_end = gene_i_end,
           length = gene_i_end - gene_i_start +1,
           num_ins = num_ins,
           nondis = max(diff(unique(c(start_is, IS_gene_i, end_is)))-1)/
             (gene_i_end - gene_i_start + 1),
           insdens = num_ins/(gene_i_end - gene_i_start + 1)
    )
    
  })
  
  
  gene_data <- do.call(rbind, gene_data)
  
  gene_data$nondis[is.na(gene_data$nondis)] <- 1
  gene_data$insdens[is.na(gene_data$insdens)] <- 0
  
  gene_data$num_ins[is.na(gene_data$num_ins)] <- 0
  
  gene_data <- 
    gene_data[!duplicated(gene_data),]
  
  gene_data$gene[is.na(gene_data$gene)] <- 
    paste("unnamed_gene_", seq_along(gene_data$gene[is.na(gene_data$gene)]),
          sep="")
  
  
  # number of insertions
  ins_sites <- 
    sort(unique(all_IS$pos))
  
  ## insdens
  ins.densities.per.gene <- gene_data %>% dplyr::select(gene, insdens)
  ins.densities.per.gene$insdens[ins.densities.per.gene$insdens == 0] <- 0.0000001
  
  
  filenname <- "./tmpData/example_simu1_insdens.csv"
  write.csv(x = ins.densities.per.gene, 
            file= filenname, 
            row.names = F)
  
  mcmc_out <- NULL
  
  mcmc_out <- 
    post_prob(
      file_path = filenname,
      print_prop_sd = T,
      print_it = T)
  
  results_InsDens <- lapply(rs, function(r){
    
    est_ess_insdens <- as.numeric(mcmc_out[[1]][,2] >= r)
    
    tibble(gene = gene_data$gene,
           essential = est_ess_insdens,
           posterior_prob = r)
    
  })
  
  
  
  ## expvsgamma
  
  results_ExpVsGamma <- lapply(ts, function(t){
    
    ExpVsGamma(gene.names = gene_data$gene, 
               gene.starts = gene_data$gene_start, 
               gene.stops = gene_data$gene_end, 
               num.ins.per.gene = gene_data$num_ins, 
               log2threshold = t
    )
    
  })
  
  
  
  
  ### binomial 
  
  results_Binomial <- lapply(ws, function(w){
    
    Binomial(ins.positions = ins_sites, 
             gene.names = gene_data$gene, 
             gene.starts = gene_data$gene_start, 
             gene.stops = gene_data$gene_end, 
             num.ins.per.gene = gene_data$num_ins, 
             genome.length =  4646332, 
             weighting = w)
    
  })
  
  
  ### ConNIS
  
  results_ConNIS <- lapply(ws, function(w){
    
    ConNIS(ins.positions = ins_sites, 
           gene.names = gene_data$gene, 
           gene.starts = gene_data$gene_start, 
           gene.stops = gene_data$gene_end, 
           num.ins.per.gene = gene_data$num_ins, 
           genome.length =  4646332, 
           weighting = w)
    
  })
  
  
  ### Geometric
  
  results_Geometric <- lapply(ws, function(w){
    
    Geometric(ins.positions = ins_sites, 
              gene.names = gene_data$gene, 
              gene.starts = gene_data$gene_start, 
              gene.stops = gene_data$gene_end, 
              num.ins.per.gene = gene_data$num_ins, 
              genome.length =  4646332, 
              weighting = w)
    
  })
  
  # Tn5Gaps
  
  ins_sites <- 
    sort(unique(
      c(1, 
        unique(all_IS$pos), 
        4646332)))
  
  
  results_Tn5Gaps <- lapply(ws, function(w){
    
    Tn5Gaps(ins.positions = ins_sites, 
            gene.names = gene_data$gene, 
            gene.starts = gene_data$gene_start, 
            gene.stops = gene_data$gene_end, 
            genome.length =  4646332, 
            weighting = w)
    
  })
  
  alpha_value <- 0.05
  
  performance_binomial <- lapply(seq_along(ws), function(i){
    
    pvalues <- results_Binomial[[i]]$p_value
    
    names(pvalues) <- seq_along(pvalues)
    sorted_pvalues <- sort(pvalues)
    
    est_ess_ben_hoch <- rep(0, length(sorted_pvalues))
    est_ess_ben_hoch[
      as.numeric(
        names(
          pvalues[
            pvalues <= max(
              sorted_pvalues[sorted_pvalues <= seq_along(sorted_pvalues)/length(sorted_pvalues)*alpha_value]
            )
          ]
        )
      )
    ] <- 1
    
    perf_binomial_ben_hoch <- bind_cols(
      method = paste("Binomial", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Benjamini-Hochberg",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_ben_hoch, true.ess = true_ess)
    )
    
    
    est_ess_bon_holm <- rep(0, length(sorted_pvalues))
    est_ess_bon_holm[
      as.numeric(
        names(
          sorted_pvalues[
            1:min(seq_along(sorted_pvalues)[
              sorted_pvalues > alpha_value/(length(sorted_pvalues)-seq_along(sorted_pvalues)+1)
            ]
            )-1
          ]
        )
      )
    ] <- 1
    
    
    perf_binomial_bon_holm <- bind_cols(
      method = paste("Binomial", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Bonferroni-Holm",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_bon_holm, true.ess = true_ess)
    )
    
    
    
    est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))
    
    perf_binomial_bon <- bind_cols(
      method = paste("Binomial", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Bonferroni",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_bon, true.ess = true_ess)
    )
    
    out <- rbind(perf_binomial_ben_hoch, perf_binomial_bon_holm, perf_binomial_bon)
    
    out
    
  })
  performance_binomial <- do.call(rbind, performance_binomial)
  
  
  performance_connis <- lapply(seq_along(ws), function(i){
    
    pvalues <- results_ConNIS[[i]]$p_value
    
    names(pvalues) <- seq_along(pvalues)
    sorted_pvalues <- sort(pvalues)
    
    est_ess_ben_hoch <- rep(0, length(sorted_pvalues))
    est_ess_ben_hoch[
      as.numeric(
        names(
          pvalues[
            pvalues <= max(
              sorted_pvalues[sorted_pvalues <= seq_along(sorted_pvalues)/length(sorted_pvalues)*alpha_value]
            )
          ]
        )
      )
    ] <- 1
    
    perf_connis_ben_hoch <- bind_cols(
      method = paste("ConNIS", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Benjamini-Hochberg",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_ben_hoch, true.ess = true_ess)
    )
    
    
    est_ess_bon_holm <- rep(0, length(sorted_pvalues))
    est_ess_bon_holm[
      as.numeric(
        names(
          sorted_pvalues[
            1:min(seq_along(sorted_pvalues)[
              sorted_pvalues > alpha_value/(length(sorted_pvalues)-seq_along(sorted_pvalues)+1)
            ]
            )-1
          ]
        )
      )
    ] <- 1
    
    
    perf_connis_bon_holm <- bind_cols(
      method = paste("ConNIS", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Bonferroni-Holm",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_bon_holm, true.ess = true_ess)
    )
    
    
    
    est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))
    
    perf_connis_bon <- bind_cols(
      method = paste("ConNIS", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Bonferroni",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_bon, true.ess = true_ess)
    )
    
    out <- rbind(perf_connis_ben_hoch, perf_connis_bon_holm, perf_connis_bon)
    
    out
    
  })
  performance_connis <- 
    do.call(rbind, performance_connis)
  
  
  
  performance_expvsgamma <- lapply(seq_along(ts), function(i){
    threshold <- ts[i]
    
    est_ess <- results_ExpVsGamma[[i]]$essential
    
    perf_expvsgamma <- bind_cols(
      method = paste("Exp. vs. Gamma", sep=""), 
      alpha_value = NA,
      adjustment = NA,
      weighting = NA,
      post_prob = NA,
      threshold_log2_expgamma = threshold,
      classification_performance(est.ess = est_ess, true.ess = true_ess)
    )
    
    out <- perf_expvsgamma
    
    out
    
  })
  performance_expvsgamma <- do.call(rbind, performance_expvsgamma)
  
  
  performance_geometric <- lapply(seq_along(ws), function(i){
    
    pvalues <- results_Geometric[[i]]$p_value
    
    names(pvalues) <- seq_along(pvalues)
    sorted_pvalues <- sort(pvalues)
    
    est_ess_ben_hoch <- rep(0, length(sorted_pvalues))
    est_ess_ben_hoch[
      as.numeric(
        names(
          pvalues[
            pvalues <= max(
              sorted_pvalues[sorted_pvalues <= seq_along(sorted_pvalues)/length(sorted_pvalues)*alpha_value]
            )
          ]
        )
      )
    ] <- 1
    
    perf_geometric_ben_hoch <- bind_cols(
      method = paste("Geometric", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Benjamini-Hochberg",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_ben_hoch, true.ess = true_ess)
    )
    
    
    est_ess_bon_holm <- rep(0, length(sorted_pvalues))
    est_ess_bon_holm[
      as.numeric(
        names(
          sorted_pvalues[
            1:min(seq_along(sorted_pvalues)[
              sorted_pvalues > alpha_value/(length(sorted_pvalues)-seq_along(sorted_pvalues)+1)
            ]
            )-1
          ]
        )
      )
    ] <- 1
    
    
    perf_geometric_bon_holm <- bind_cols(
      method = paste("Geometric", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Bonferroni-Holm",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_bon_holm, true.ess = true_ess)
    )
    
    
    
    est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))
    
    perf_geometric_bon <- bind_cols(
      method = paste("Geometric", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Bonferroni",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_bon, true.ess = true_ess)
    )
    
    out <- rbind(perf_geometric_ben_hoch, perf_geometric_bon_holm, perf_geometric_bon)
    
    out
    
  })
  performance_geometric <- 
    do.call(rbind, performance_geometric)
  
  
  performance_insdens <- lapply(seq_along(rs), function(r){
    
    threshold <- rs[r]
    
    est_ess <- results_InsDens[[r]]$essential
    
    perf_insdens <- bind_cols(
      method = paste("InsDens", sep=""), 
      alpha_value = NA,
      adjustment = NA,
      weighting = NA,
      post_prob = threshold,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess, true.ess = true_ess)
    )
    
    out <- perf_insdens
    
    out
    
  })
  performance_insdens <- do.call(rbind, performance_insdens)
  
  
  performance_tn5gaps <- lapply(seq_along(ws), function(i){
    
    pvalues <- results_Tn5Gaps[[i]]$p_value
    
    names(pvalues) <- seq_along(pvalues)
    sorted_pvalues <- sort(pvalues)
    
    est_ess_ben_hoch <- rep(0, length(sorted_pvalues))
    est_ess_ben_hoch[
      as.numeric(
        names(
          pvalues[
            pvalues <= max(
              sorted_pvalues[sorted_pvalues <= seq_along(sorted_pvalues)/length(sorted_pvalues)*alpha_value]
            )
          ]
        )
      )
    ] <- 1
    
    perf_tn5gaps_ben_hoch <- bind_cols(
      method = paste("Tn5Gaps", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Benjamini-Hochberg",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_ben_hoch, true.ess = true_ess)
    )
    
    
    est_ess_bon_holm <- rep(0, length(sorted_pvalues))
    est_ess_bon_holm[
      as.numeric(
        names(
          sorted_pvalues[
            1:min(seq_along(sorted_pvalues)[
              sorted_pvalues > alpha_value/(length(sorted_pvalues)-seq_along(sorted_pvalues)+1)
            ]
            )-1
          ]
        )
      )
    ] <- 1
    
    
    perf_tn5gaps_bon_holm <- bind_cols(
      method = paste("Tn5Gaps", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Bonferroni-Holm",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_bon_holm, true.ess = true_ess)
    )
    
    
    
    est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))
    
    perf_tn5gaps_bon <- bind_cols(
      method = paste("Tn5Gaps", sep=""), 
      alpha_value = alpha_value,
      adjustment = "Bonferroni",
      weighting = ws[i],
      post_prob = NA,
      threshold_log2_expgamma = NA,
      classification_performance(est.ess = est_ess_bon, true.ess = true_ess)
    )
    
    out <- rbind(perf_tn5gaps_ben_hoch, perf_tn5gaps_bon_holm, perf_tn5gaps_bon)
    
    out
    
  })
  performance_tn5gaps <- do.call(rbind, performance_tn5gaps)
  
  
  
  performance_all <- 
    bind_rows(
      performance_binomial,
      performance_connis,
      performance_expvsgamma,
      performance_geometric,
      performance_insdens,
      performance_tn5gaps
    )
  
  performance_all$trimming_start <- trimming_start
  performance_all$trimming_end <- trimming_end
  
  performance_all
  
})


performances <- do.call(rbind, performances)


saveRDS(performances, file="./performance/Performance_example_simu1.RDS")



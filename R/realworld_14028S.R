library(tidyverse)
library(insdens)
library(MASS)
library(gmp)

setwd("./R")

# https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.01723/full#supplementary-material
# https://genomics.lbl.gov/supplemental/rbarseq/html/Keio/
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099820#pone-0099820-t001

source("functions.R")


### Mandal et all data
all_genes_14028S <- readRDS("./14028s_data/Typhimurium14028S_genes.RDS")
# some locusIDs miss a leading 0:
all_genes_14028S$locusId[928:1000] <- paste("STM14_0", 927:999, sep="")
# some loucsIDs are duplicated because of two different gene_syn; remove duplicate
all_genes_14028S <- all_genes_14028S[!duplicated(all_genes_14028S$locusId),]

# Original IS data from the authors
all_IS_raw <- readxl::read_xlsx("./14028s_data/Data_Sheet_2.xlsx", skip = 4)


# true essentials of Kaio library
true_ess_kaio <- readRDS("./bw25113_data/essential_genes_kaio.RDS")

# knockouts of Porwollik
# no and MGD?
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099820#pone-0099820-t001
all_genes_porwollik <- readxl::read_xlsx("./14028s_data/pone.0099820.s001.xlsx", skip = 1)

all_genes_porwollik <- all_genes_porwollik[1:5585,]

# Keep only geens withlocusIDs in both
genes_in_both <- all_genes_14028S$locusId[all_genes_14028S$locusId %in% all_genes_porwollik$`Target gene`]

all_genes_porwollik <- all_genes_porwollik[ all_genes_porwollik$`Target gene` %in% genes_in_both,]
all_genes_porwollik <- all_genes_porwollik[-which(duplicated(all_genes_porwollik$`Target gene`)),]

all_genes_14028S <- all_genes_14028S[all_genes_14028S$locusId %in% genes_in_both,]



true_ess_pro <- 
  na.omit(all_genes_porwollik$`gene symbol`[
    which(
      na.omit(all_genes_porwollik$`Mutant status KanR (* indicates present in more than one location)` == "deleted in an MGD" & 
                all_genes_porwollik$`Mutant status CmR (* indicates present in more than one location)` == "deleted in an MGD") |
        na.omit(all_genes_porwollik$`Mutant status KanR (* indicates present in more than one location)` == "no" & 
                  all_genes_porwollik$`Mutant status CmR (* indicates present in more than one location)` == "no"))])

true_ess <- union(true_ess_kaio, true_ess_pro)

true_ess <- true_ess[true_ess %in% na.omit(all_genes_14028S$gene) | 
                            true_ess %in% na.omit(all_genes_14028S$gene_syn)]

### noise thresholds
read_count_thresholds <- c(0,1,3,5,10)

#### trimming
list_trimming_start <- c(0, 0.05, 0.1)
list_trimming_end <- c(0, 0.05, 0.1)

alpha_value <- 0.05
ws <- ws <- c(0.01, seq(0.05,1,0.05))
ts <- c(seq(0.1,2, 0.1),3:12)
rs <- c(seq(0.01,0.1, 0.01), seq(0.2,0.9,0.1), seq(0.91,0.99,0.01))


performances_read_count_thresholds <- 
  lapply(read_count_thresholds, function(read_count_threshold){
    
    # set read count threshold
    all_IS <- tibble(pos=unlist(all_IS_raw[which(all_IS_raw[,2] > read_count_threshold),1]))
    
    
    performances_trimming <- lapply(seq_along(list_trimming_start), function(n){
      
      trimming_start <- list_trimming_start[n]
      trimming_end <- list_trimming_end[n]
      
      gene_data <- lapply(1:nrow(all_genes_14028S), function(i){
        
        gene_i_start <- all_genes_14028S$gene_start[i] + 
          floor(trimming_start * (all_genes_14028S$gene_end[i] - all_genes_14028S$gene_start[i]+1))
        
        gene_i_end <- all_genes_14028S$gene_end[i] - 
          floor(trimming_end * (all_genes_14028S$gene_end[i] - all_genes_14028S$gene_start[i]+1))
        
        IS_gene_i <- 
          all_IS$pos[which(gene_i_start <= all_IS$pos &
                             gene_i_end >= all_IS$pos)]
        
        num_ins <- length(IS_gene_i)
        
        gene <- all_genes_14028S$gene[i]
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
        
        tibble(entry = all_genes_14028S$locusId[i],
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
      
      
      true_ess <- as.numeric(gene_data$gene %in% true_ess)
      
      
      # number of insertions
      ins_sites <- 
        sort(unique(all_IS$pos))
      
      ## insdens
      ins.densities.per.gene <- gene_data %>% dplyr::select(gene, insdens)
      ins.densities.per.gene$insdens[ins.densities.per.gene$insdens == 0] <- 0.0000001
      
      
      filenname <- "./tmpData/MG1655_insdens.csv"
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
      
      
      print("exp")
      ## expvsgamma
      
      results_ExpVsGamma <- lapply(ts, function(t){
        
        ExpVsGamma(gene.names = gene_data$gene, 
                   gene.starts = gene_data$gene_start, 
                   gene.stops = gene_data$gene_end, 
                   num.ins.per.gene = gene_data$num_ins, 
                   log2threshold = t
        )
        
      })
      
      
      
      print("bin")
      ### binomial 
      
      results_Binomial <- lapply(ws, function(w){
        
        Binomial(ins.positions = ins_sites, 
                 gene.names = gene_data$gene, 
                 gene.starts = gene_data$gene_start, 
                 gene.stops = gene_data$gene_end, 
                 num.ins.per.gene = gene_data$num_ins, 
                 genome.length = 4870265, 
                 weighting = w)
        
      })
      
      
      ### ConNIS
      print("con")
      results_ConNIS <- lapply(ws, function(w){
        
        ConNIS(ins.positions = ins_sites, 
               gene.names = gene_data$gene, 
               gene.starts = gene_data$gene_start, 
               gene.stops = gene_data$gene_end, 
               num.ins.per.gene = gene_data$num_ins, 
               genome.length = 4870265, 
               weighting = w)
        
      })
      
      
      ### Geometric
      print("geo")
      results_Geometric <- lapply(ws, function(w){
        
        Geometric(ins.positions = ins_sites, 
                  gene.names = gene_data$gene, 
                  gene.starts = gene_data$gene_start, 
                  gene.stops = gene_data$gene_end, 
                  num.ins.per.gene = gene_data$num_ins, 
                  genome.length = 4870265, 
                  weighting = w)
        
      })
      
      # Tn5Gaps
      
      ins_sites <- 
        sort(unique(
          c(1, 
            unique(all_IS$pos), 
            4870265)))
      
      print("tn5")
      results_Tn5Gaps <- lapply(ws, function(w){
        
        Tn5Gaps(ins.positions = ins_sites, 
                gene.names = gene_data$gene, 
                gene.starts = gene_data$gene_start, 
                gene.stops = gene_data$gene_end, 
                genome.length = 4870265, 
                weighting = w)
        
      })
      
      
      
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
      performance_all$read_count_threshold <- read_count_threshold
      
      performance_all
      
    })
    
    
    performances_trimming <- do.call(rbind, performances_trimming)
    
    
    
  })

performances_read_count_thresholds <- do.call(rbind, performances_read_count_thresholds)

saveRDS(performances_read_count_thresholds, file="./performance/Performance_14028S_realWorld.RDS")



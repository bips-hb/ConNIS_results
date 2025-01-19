library(tidyverse)
library(insdens)
library(MASS)
library(gmp)
library(parallel)

setwd("./R")

# Goodall et al. (2018)
# https://journals.asm.org/doi/10.1128/mbio.02096-17
# https://genomics.lbl.gov/supplemental/rbarseq/html/Keio/

source("functions.R")
all_genes_BW25113_raw <- tibble(
  read.csv("./bw25113_data/bw25113_genbank.txt", sep="\n", header = F))

all_genes_BW25113_raw <- all_genes_BW25113_raw[grep(">", all_genes_BW25113_raw$V1),]

all_genes_BW25113 <- lapply(1:nrow(all_genes_BW25113_raw), function(i){
  tmp_str <- as.character(unlist(all_genes_BW25113_raw[i,]))
  gene <- stringr::str_match(tmp_str, "\\[gene=\\s*(.*?)\\s*\\]")[,2]
  locus_tag <- stringr::str_match(tmp_str, "\\[locus_tag=\\s*(.*?)\\s*\\]")[,2]
  protein <- stringr::str_match(tmp_str, "\\[protein=\\s*(.*?)\\s*\\]")[,2]
  position <- stringr::str_match(tmp_str, "\\[location=\\s*(.*?)\\s*\\]")[,2]
  
  if(length(grep("join", position))>0){
    
    insert_comp <- FALSE
    if(length(grep("complement", position)) > 0){
      left_comp <- "complement("
      right_comp <- ")"
      insert_comp <- TRUE
    }
    
    position <- 
      stringr::str_match(tmp_str, "join\\(\\s*(.*?)\\s*\\)")[,2]
    
    position <- 
      str_replace_all(position, ",", "..")
    
    position <- paste(
      min(str_split(position, "\\.\\.")[[1]]), 
      "..",
      max(str_split(position, "\\.\\.")[[1]]),
      sep=""
    )  
    
    if(insert_comp){
      position <- 
        paste(left_comp,
              position,
              right_comp,
              sep="")
    }
    
  }
  
  
  if(length(grep("complement", position)) > 0){
    Orientation <- "minus"
    position <- str_remove(position, "complement\\(")
    position <- str_remove(position, "\\)")
  }
  
  if(length(grep("..", position)) > 0){
    gene_start <- as.numeric(str_split(position, "\\.\\.")[[1]][1])
    gene_end <- as.numeric(str_split(position, "\\.\\.")[[1]][2])
  }
  
  if(length(grep(" ", position)) > 0){
    gene_start <- as.numeric(str_split(position, " ")[[1]][1])
    gene_end <- as.numeric(str_split(position, " ")[[1]][2])
  }
  
  gene_length <- gene_end-gene_start+1
  
  tibble(
    gene,
    locus_tag,
    gene_start,
    gene_end,
    gene_length,
    protein
  )
})

all_genes_BW25113 <- do.call(rbind,all_genes_BW25113)

# need to remove "YdfJ" because of false start stop annotation (its length is over 20000)
all_genes_BW25113 <- all_genes_BW25113[all_genes_BW25113$gene_length < 10000,] 

true_ess <- readRDS("./bw25113_data/essential_genes_kaio.RDS")
true_ess <- true_ess[true_ess %in% all_genes_BW25113$gene]

#Is from Goodall
all_IS_raw <- read.csv("./bw25113_data/Ecoli_BW25113_chlor-tn_position_depth_count.txt")
all_IS_full <- tibble(pos=which(all_IS_raw[,1] >0),
                      n=all_IS_raw[which(all_IS_raw[,1] >0),1],
                      strand = NA)
# Clsuter size
nc <- 50

# Number of Simulation runs:
n_sim <- 100


### noise thresholds
noise_thresholds <- c(0,1,3)

#### trimming
list_trimming_start <- c(0, 0.05)
list_trimming_end <- c(0, 0.05)

# significance level
alpha_value <- 0.05

ws <- seq(0.1, 1, 0.1)
ts <- c(2:12)
rs <- c(0.01, seq(0.1, 0.9, 0.1), 0.99)



for(numIS in c(100000, 200000, 400000)){
  
  for(noise_threshold in noise_thresholds){
    
    cl <- makeCluster(nc, type="PSOCK", outfile = "")
    
    clusterExport(cl,
                  list(
                    "all_genes_BW25113",
                    "all_IS_raw",
                    "all_IS_full",
                    "noise_threshold",
                    "list_trimming_start",
                    "list_trimming_end",
                    "alpha_value",
                    "ws",
                    "ts",
                    "rs",
                    "numIS",
                    "true_ess",
                    "bc_prob_bigZ",
                    "classification_performance",
                    "bc_freq_bigZ",
                    "pgumbel",
                    "Binomial",
                    "ConNIS",
                    "ExpVsGamma",
                    "Geometric",
                    "Tn5Gaps"
                  ))
    
    clusterEvalQ(cl, {
      library(tidyverse)
      library(insdens)
      library(MASS)
      library(gmp)}
    )
    
    
    BW25113_simulation <- parLapplyLB(cl, 1:n_sim, function(sim_run){
      
      set.seed(123*numIS+sim_run)
      # subsample of IS
      all_IS <- all_IS_full[sort(sample(1:nrow(all_IS_full), numIS, replace = F)), ]
      
      
      # Only IS with more than noise_threshold insertion count are used
      all_IS <- all_IS[all_IS$n > noise_threshold,]
      
      performances_trimming <- lapply(seq_along(list_trimming_start), function(n){
        
        trimming_start <- list_trimming_start[n]
        trimming_end <- list_trimming_end[n]
        
        gene_data <- lapply(1:nrow(all_genes_BW25113), function(i){
          
          gene_i_start <- all_genes_BW25113$gene_start[i] + 
            floor(trimming_start * (all_genes_BW25113$gene_end[i] - all_genes_BW25113$gene_start[i]+1))
          
          gene_i_end <- all_genes_BW25113$gene_end[i] - 
            floor(trimming_end * (all_genes_BW25113$gene_end[i] - all_genes_BW25113$gene_start[i]+1))
          
          IS_gene_i <- 
            all_IS$pos[which(gene_i_start <= all_IS$pos &
                               gene_i_end >= all_IS$pos)]
          
          num_ins <- length(IS_gene_i)
          
          gene <- all_genes_BW25113$gene[i]
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
          
          tibble(entry = all_genes_BW25113$locusId[i],
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
        
        
        filenname <- paste("./tmpData/MG1655_insdens_", sim_run, "_", numIS, ".csv", sep="")
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
                   genome.length =  4631469, 
                   weighting = w)
          
        })
        
        
        ### ConNIS
        
        results_ConNIS <- lapply(ws, function(w){
          
          ConNIS(ins.positions = ins_sites,
                 gene.names = gene_data$gene,
                 gene.starts = gene_data$gene_start,
                 gene.stops = gene_data$gene_end,
                 num.ins.per.gene = gene_data$num_ins,
                 genome.length =  4631469,
                 weighting = w)
          
        })
        
        
        ### Geometric
        
        results_Geometric <- lapply(ws, function(w){
          
          Geometric(ins.positions = ins_sites,
                    gene.names = gene_data$gene,
                    gene.starts = gene_data$gene_start,
                    gene.stops = gene_data$gene_end,
                    num.ins.per.gene = gene_data$num_ins,
                    genome.length =  4631469,
                    weighting = w)
          
        })
        
        # Tn5Gaps
        
        ins_sites <- 
          sort(unique(
            c(1, 
              unique(all_IS$pos), 
              4631469)))
        
        
        results_Tn5Gaps <- lapply(ws, function(w){
          
          Tn5Gaps(ins.positions = ins_sites,
                  gene.names = gene_data$gene,
                  gene.starts = gene_data$gene_start,
                  gene.stops = gene_data$gene_end,
                  genome.length =  4631469,
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
        performance_all$noise_threshold <- noise_threshold
        
        performance_all
        
      })
      
      
      performances_trimming <- do.call(rbind, performances_trimming)
      
      performances_trimming$sim_run <- sim_run
      performances_trimming$noise_threshold <- noise_threshold
      performances_trimming
      
      
      
    })
    
    BW25113_simulation <- do.call(rbind, BW25113_simulation)
    stopCluster(cl)
    saveRDS(BW25113_simulation, 
            file=paste("./performance/Performance_BW25113_simulation",
                       "_numIS_", numIS, 
                       "_noise_threshold_", noise_threshold,
                       ".RDS", sep=""))
  
    }
  
}





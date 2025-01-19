library(tidyverse)
library(insdens)
library(MASS)
library(gmp)
library(parallel)

setwd("./R")

source("functions.R")

### Simulation data: 400,000 IS, 2% noise and lambda=0.7  
simulated_data <- readRDS("./simulatedData/SimuData_Tnseq_GeneClusters_ess_ORF_uniform_numLoci_400000_lambda_0.7_bp_per_wave_2e+06_sine_scaling_factor_1.3_NBnc_100_NBdis_1_NBp_0.3_tecnoise_0.02.RDS")

set.seed(99243)
random_sim_run <- sample(1:100, 1)
tnseq_data <- simulated_data[[random_sim_run]]$tnseqData[, c(1,4,5)]
names(tnseq_data)<- c("gene", "gene_start", "gene_end")

all_IS <- tibble(pos=simulated_data[[random_sim_run]]$observed_all_unique_loci)
  

original_data <- tnseq_data

nc <- 50
cl <- makeCluster(nc, type="PSOCK", outfile = "")

trimming_end <- trimming_start <- 0.05 
alpha_value <- 0.05
n_draws <- 500
ws <- ws <- c(0.01, seq(0.05,1,0.05))
ts <- c(seq(0.1,2, 0.1),3:12)
rs <- c(seq(0.01,0.1, 0.01), seq(0.2,0.9,0.1), seq(0.91,0.99,0.01))

clusterExport(cl,
              list(
                "original_data",
                "all_IS",
                "trimming_start",
                "trimming_end",
                "alpha_value",
                "ws",
                "ts", 
                "rs",
                "n_draws",
                "classification_performance",
                "prob_seq_misses",
                "freq_seq_misses",
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


# Stability Binomial
out_perm_binomial <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)/2)))))
  
  gene_data <- lapply(1:nrow(original_data), function(i){
    
    gene_i_start <- original_data$gene_start[i] + 
      floor(trimming_start * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    gene_i_end <- original_data$gene_end[i] - 
      floor(trimming_end * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    IS_gene_i <- 
      all_IS_perm$pos[which(gene_i_start <= all_IS_perm$pos &
                              gene_i_end >= all_IS_perm$pos)]
    
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
    
    tibble(gene = gene,
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
  
  ins_sites <- 
    sort(unique(all_IS_perm$pos))
  
  
  results_Binomial <- lapply(ws, function(w){
    
    Binomial(ins.positions = ins_sites,
             gene.names = gene_data$gene,
             gene.starts = gene_data$gene_start,
             gene.stops = gene_data$gene_end,
             num.ins.per.gene = gene_data$num_ins,
             genome.length =  4641652,
             weighting = w)
    
  })
  
  est_binomial <- lapply(seq_along(ws), function(i){
    
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
    
    est_ess_ben_hoch
    
    est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))
    est_ess_bon
  })
  
  est_binomial
  
})

stability_binomial <- sapply(seq_along(ws), function(j){
  a <- lapply(out_perm_binomial, function(i){i[[j]]})
  
  sum(na.omit((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws)))/sum(na.omit(Reduce("+", a)>0))
  
})

saveRDS(out_perm_binomial, "./results/subsamples_binomial_example_simu5.RDS")


# stability ConNIS
out_perm_connis <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)/2)))))
  
  gene_data <- lapply(1:nrow(original_data), function(i){
    
    gene_i_start <- original_data$gene_start[i] + 
      floor(trimming_start * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    gene_i_end <- original_data$gene_end[i] - 
      floor(trimming_end * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    IS_gene_i <- 
      all_IS_perm$pos[which(gene_i_start <= all_IS_perm$pos &
                              gene_i_end >= all_IS_perm$pos)]
    
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
    
    tibble(gene = gene,
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
  
  ins_sites <- 
    sort(unique(all_IS_perm$pos))
  
  
  results_ConNIS <- lapply(ws, function(w){
    
    ConNIS(ins.positions = ins_sites, 
           gene.names = gene_data$gene, 
           gene.starts = gene_data$gene_start, 
           gene.stops = gene_data$gene_end, 
           num.ins.per.gene = gene_data$num_ins, 
           genome.length =  4641652, 
           weighting = w)
    
  })
  
  est_connis <- lapply(seq_along(ws), function(i){
    
    pvalues <- results_ConNIS[[i]]$p_value
    
    names(pvalues) <- seq_along(pvalues)
    sorted_pvalues <- sort(pvalues)
    # 
    #   est_ess_ben_hoch <- rep(0, length(sorted_pvalues))
    #   est_ess_ben_hoch[
    #     as.numeric(
    #       names(
    #         pvalues[
    #           pvalues <= max(
    #             sorted_pvalues[sorted_pvalues <= seq_along(sorted_pvalues)/length(sorted_pvalues)*alpha_value]
    #           )
    #         ]
    #       )
    #     )
    #   ] <- 1
    # 
    #   est_ess_ben_hoch
    
    est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))
    est_ess_bon
  })
  
  est_connis
  
})

stability_connis <- sapply(seq_along(ws), function(j){
  a <- lapply(out_perm_connis, function(i){i[[j]]})
  
  sum(na.omit((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws)))/sum(na.omit(Reduce("+", a)>0))
  
})

saveRDS(out_perm_connis, "./results/subsamples_connis_example_simu5.RDS")


# stability Exp. vs Gamma
out_perm_expvsgamma <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)/2)))))
  
  gene_data <- lapply(1:nrow(original_data), function(i){
    
    gene_i_start <- original_data$gene_start[i] + 
      floor(trimming_start * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    gene_i_end <- original_data$gene_end[i] - 
      floor(trimming_end * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    IS_gene_i <- 
      all_IS_perm$pos[which(gene_i_start <= all_IS_perm$pos &
                              gene_i_end >= all_IS_perm$pos)]
    
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
    
    tibble(gene = gene,
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
  
  ins_sites <- 
    sort(unique(all_IS_perm$pos))
  
  results_ExpVsGamma <- lapply(ts, function(t){
    
    ExpVsGamma(gene.names = gene_data$gene,
               gene.starts = gene_data$gene_start,
               gene.stops = gene_data$gene_end,
               num.ins.per.gene = gene_data$num_ins,
               log2threshold = t
    )
    
  })
  
  est_expvsgamma <- lapply(results_ExpVsGamma, function(i){
    as.numeric(unlist(i[,2]))
  })
  
  est_expvsgamma
  
})

stability_expvsgamma <- sapply(seq_along(ts), function(j){
  a <- lapply(out_perm_expvsgamma, function(i){i[[j]]})
  
  sum(na.omit((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws)))/sum(na.omit(Reduce("+", a)>0))
  
})

saveRDS(out_perm_expvsgamma, "./results/subsamples_expvsgamma_example_simu5.RDS")


# stability Geometric
out_perm_geometric <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)/2)))))
  
  gene_data <- lapply(1:nrow(original_data), function(i){
    
    gene_i_start <- original_data$gene_start[i] + 
      floor(trimming_start * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    gene_i_end <- original_data$gene_end[i] - 
      floor(trimming_end * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    IS_gene_i <- 
      all_IS_perm$pos[which(gene_i_start <= all_IS_perm$pos &
                              gene_i_end >= all_IS_perm$pos)]
    
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
    
    tibble(gene = gene,
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
  
  ins_sites <- 
    sort(unique(all_IS_perm$pos))
  results_Geometric <- lapply(ws, function(w){
    
    Geometric(ins.positions = ins_sites, 
              gene.names = gene_data$gene, 
              gene.starts = gene_data$gene_start, 
              gene.stops = gene_data$gene_end, 
              num.ins.per.gene = gene_data$num_ins, 
              genome.length =  4641652, 
              weighting = w)
    
  })
  
  est_geometric <- lapply(seq_along(ws), function(i){
    
    pvalues <- results_Geometric[[i]]$p_value
    
    names(pvalues) <- seq_along(pvalues)
    sorted_pvalues <- sort(pvalues)
    
    # est_ess_ben_hoch <- rep(0, length(sorted_pvalues))
    # est_ess_ben_hoch[
    #   as.numeric(
    #     names(
    #       pvalues[
    #         pvalues <= max(
    #           sorted_pvalues[sorted_pvalues <= seq_along(sorted_pvalues)/length(sorted_pvalues)*alpha_value]
    #         )
    #       ]
    #     )
    #   )
    # ] <- 1
    # 
    # est_ess_ben_hoch
    
    est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))
    est_ess_bon
    
  })
  
  est_geometric
  
})

stability_geometric <- sapply(seq_along(ws), function(j){
  a <- lapply(out_perm_geometric, function(i){i[[j]]})
  
  sum(na.omit((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws)))/sum(na.omit(Reduce("+", a)>0))
  
})

saveRDS(out_perm_geometric, "./results/subsamples_geometric_example_simu5.RDS")


# stability InsDens
out_perm_insdens <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)/2)))))
  
  gene_data <- lapply(1:nrow(original_data), function(i){
    
    gene_i_start <- original_data$gene_start[i] + 
      floor(trimming_start * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    gene_i_end <- original_data$gene_end[i] - 
      floor(trimming_end * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    IS_gene_i <- 
      all_IS_perm$pos[which(gene_i_start <= all_IS_perm$pos &
                              gene_i_end >= all_IS_perm$pos)]
    
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
    
    tibble(gene = gene,
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
  
  ins_sites <- 
    sort(unique(all_IS_perm$pos))
  
  ins.densities.per.gene <- gene_data %>% dplyr::select(gene, insdens)
  ins.densities.per.gene$insdens[ins.densities.per.gene$insdens == 0] <- 0.0000001
  
  
  filenname <- paste("./tmpData/insdens_", draw_i, ".csv", sep="")
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
  
  est_insdens <- lapply(results_InsDens, function(i){
    as.numeric(unlist(i[,2]))
  })
  
  est_insdens
  
})

stability_insdens <- sapply(seq_along(rs), function(j){
  a <- lapply(out_perm_insdens, function(i){i[[j]]})
  
  sum(na.omit((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws)))/sum(na.omit(Reduce("+", a)>0))
  
})

saveRDS(out_perm_insdens, "./results/subsamples_insdens_example_simu5.RDS")


# stability Tn5Gaps
out_perm_tn5gaps <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)/2)))))
  
  gene_data <- lapply(1:nrow(original_data), function(i){
    
    gene_i_start <- original_data$gene_start[i] + 
      floor(trimming_start * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    gene_i_end <- original_data$gene_end[i] - 
      floor(trimming_end * (original_data$gene_end[i] - original_data$gene_start[i]+1))
    
    IS_gene_i <- 
      all_IS_perm$pos[which(gene_i_start <= all_IS_perm$pos &
                              gene_i_end >= all_IS_perm$pos)]
    
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
    
    tibble(gene = gene,
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
  
  ins_sites <- 
    sort(unique(all_IS_perm$pos))
  
  results_Tn5Gaps <- lapply(ws, function(w){
    
    Tn5Gaps(ins.positions = ins_sites, 
            gene.names = gene_data$gene, 
            gene.starts = gene_data$gene_start, 
            gene.stops = gene_data$gene_end, 
            genome.length =  4641652, 
            weighting = w)
    
  })
  
  est_tn5gaps <- lapply(seq_along(ws), function(i){
    
    pvalues <- results_Tn5Gaps[[i]]$p_value
    
    names(pvalues) <- seq_along(pvalues)
    sorted_pvalues <- sort(pvalues)
    
    # est_ess_ben_hoch <- rep(0, length(sorted_pvalues))
    # est_ess_ben_hoch[
    #   as.numeric(
    #     names(
    #       pvalues[
    #         pvalues <= max(
    #           sorted_pvalues[sorted_pvalues <= seq_along(sorted_pvalues)/length(sorted_pvalues)*alpha_value]
    #         )
    #       ]
    #     )
    #   )
    # ] <- 1
    # 
    # est_ess_ben_hoch
    
    est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))
    est_ess_bon
  })
  
  est_tn5gaps
  
})

stability_tn5gaps <- sapply(seq_along(ws), function(j){
  a <- lapply(out_perm_tn5gaps, function(i){i[[j]]})
  
  sum(na.omit((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws)))/sum(na.omit(Reduce("+", a)>0))
  
})

saveRDS(out_perm_tn5gaps, "./results/subsamples_tn5gaps_example_simu5.RDS")


stopCluster(cl)


stabilities <- list(
  stability_binomial,
  stability_connis,
  stability_expvsgamma,
  stability_geometric,
  stability_insdens,
  stability_tn5gaps
)

names(stabilities) <- 
  c("stability_binomial",
    "stability_connis",
    "stability_expvsgamma",
    "stability_geometric",
    "stability_insdens",
    "stability_tn5gaps")

names(stabilities[[1]]) <- as.character(ws)
names(stabilities[[2]]) <- as.character(ws)
names(stabilities[[3]]) <- as.character(ts)
names(stabilities[[4]]) <- as.character(ws)
names(stabilities[[5]]) <- as.character(rs)
names(stabilities[[6]]) <- as.character(ws)


saveRDS(stabilities, "./performance/stabilities_example_simu5.RDS")

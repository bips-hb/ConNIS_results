library(tidyverse)
library(insdens)
library(MASS)
library(gmp)
library(parallel)

setwd("./R")

source("functions.R")

### Mandal et all data
all_genes_14028S <- readRDS("./14028s_data/Typhimurium14028S_genes.RDS")
# some locusIDs miss a leading 0:
all_genes_14028S$locusId[928:1000] <- paste("STM14_0", 927:999, sep="")
# some loucsIDs are duplicated because of two different gene_syn; remove duplicate
all_genes_14028S <- all_genes_14028S[!duplicated(all_genes_14028S$locusId),]


# Original data from the authors

all_IS_raw <- readxl::read_xlsx("./14028s_data/Data_Sheet_2.xlsx", skip = 4)


# manual inspection showed possibility for spurious IS; omit IS with 1 read count
read_count_threshold <- 1
all_IS <- tibble(pos=unlist(all_IS_raw[which(all_IS_raw[,2] > read_count_threshold),1]))



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


original_data <- all_genes_14028S

genome_length <- 4870265

nc <- 50
cl <- makeCluster(nc, type="PSOCK", outfile = "")

frac_subsamples <- 0.5
trimming_end <- trimming_start <- 0.05 
alpha_value <- 0.05
n_draws <- 500
ws <- ws <- c(0.01, seq(0.05,1,0.05))
ts <- c(1:12)
rs <- c(seq(0.1,0.9,0.1), seq(0.91,0.99,0.01))

clusterExport(cl,
              list(
                "original_data",
                "genome_length",
                "frac_subsamples",
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
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)*frac_subsamples)))))
  
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
             genome.length = genome_length,
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
  
  sum((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws))/sum(Reduce("+", a)>0)
  
})

saveRDS(out_perm_binomial, "./results/subsamples_binomial_14028s.RDS")


# stability ConNIS
out_perm_connis <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)*frac_subsamples)))))
  
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
           genome.length = genome_length, 
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
  
  sum((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws))/sum(Reduce("+", a)>0)
  
})

saveRDS(out_perm_connis, "./results/subsamples_connis_14028s.RDS")


# stability Exp. vs Gamma
out_perm_expvsgamma <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)*frac_subsamples)))))
  
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
  
  sum((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws))/sum(Reduce("+", a)>0)
  
})

saveRDS(out_perm_expvsgamma, "./results/subsamples_expvsgamma_14028s.RDS")


# stability Geometric
out_perm_geometric <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)*frac_subsamples)))))
  
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
              genome.length = genome_length, 
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
  
  sum((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws))/sum(Reduce("+", a)>0)
  
})

saveRDS(out_perm_geometric, "./results/subsamples_geometric_14028s.RDS")


# stability InsDens
out_perm_insdens <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)*frac_subsamples)))))
  
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
  
  sum((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws))/sum(Reduce("+", a)>0)
  
})

saveRDS(out_perm_insdens, "./results/subsamples_insdens_14028s.RDS")


# stability Tn5Gaps
out_perm_tn5gaps <- parLapply(cl, 1:n_draws, function(draw_i){
  
  print(draw_i)
  
  set.seed(draw_i)
  all_IS_perm <- tibble(pos=sort(as.numeric(sample(unlist(all_IS), floor(nrow(all_IS)*frac_subsamples)))))
  
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
            genome.length = genome_length, 
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
  
  sum((Reduce("+", a)/n_draws)*(1-Reduce("+", a)/n_draws), na.rm=T)/sum(Reduce("+", a)>0, na.rm=T)
  
})

saveRDS(out_perm_tn5gaps, "./results/subsamples_tn5gaps_14028s.RDS")


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

saveRDS(stabilities, "./performance/stabilities_14028s.RDS")

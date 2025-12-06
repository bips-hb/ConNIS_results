library(tidyverse)
library(insdens)
library(MASS)
library(gmp)
library(ggvenn)

setwd("./R")

# https://journals.asm.org/doi/10.1128/mbio.00306-15
# https://genomics.lbl.gov/supplemental/rbarseq/html/Keio/

source("functions.R")



######## E. Coli #######

### Wetmore es all data
all_genes_BW25113 <- read_csv("./bw25113_data/genes_wetmore_et_al.csv")

all_genes_BW25113 <- tibble(all_genes_BW25113)

colnames(all_genes_BW25113)[c(5:6,8)] <-
  c("gene_start", "gene_end", "gene")

all_genes_BW25113 <- all_genes_BW25113[all_genes_BW25113$gene_start <= 4631469, ]


# read insertion data from Wetmore et al. (2015) 
# Original source: https://genomics.lbl.gov/supplemental/rbarseq/html/Keio/all.poolcount

all_IS_raw <- read_tsv("./bw25113_data/all.poolcount")


# true essentials of Kaio library
true_ess <- readRDS("./bw25113_data/essential_genes_kaio.RDS")
true_ess <- true_ess[true_ess %in% all_genes_BW25113$gene]
true_ess_names <- true_ess

### noise thresholds
read_count_threshold <- 0

#### trimming
trimming_start <- 0.05
trimming_end <- 0.05

alpha_value <- 0.05
ws <- ws <- c(0.01, seq(0.05,1,0.05))
ts <- c(seq(0.1,2, 0.1),3:12)
rs <- c(seq(0.01,0.1, 0.01), seq(0.2,0.9,0.1), seq(0.91,0.99,0.01))


# Only the T0 data unconstrained in LB medium is used, i.e. sets with IDs IT001, IT002,
# IT049 and IT050 (see https://genomics.lbl.gov/supplemental/rbarseq/html/Keio/exps)
# use different thresholds for denoising
all_IS <- 
  tibble(pos = sort(
    unique(
      unlist(
        all_IS_raw[which(
          (all_IS_raw$Keio_ML9_set1.IT001 + 
             all_IS_raw$Keio_ML9_set1.IT002 + 
             all_IS_raw$Keio_ML9_set1.IT049 + 
             all_IS_raw$Keio_ML9_set1.IT050) > read_count_threshold),"pos"]
      )
    )
  )
  )



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


filenname <- "./tmpData/BW25113_insdens.csv"
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
           genome.length = 4631469, 
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
         genome.length = 4631469, 
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
            genome.length = 4631469, 
            weighting = w)
  
})

# Tn5Gaps

ins_sites <- 
  sort(unique(
    c(1, 
      unique(all_IS$pos), 
      4631469)))

print("tn5")
results_Tn5Gaps <- lapply(ws, function(w){
  
  Tn5Gaps(ins.positions = ins_sites, 
          gene.names = gene_data$gene, 
          gene.starts = gene_data$gene_start, 
          gene.stops = gene_data$gene_end, 
          genome.length = 4631469, 
          weighting = w)
  
})


# Selected Genes by Binomial
# from generate_stability_tables.R
index_stability <- which(round(ws,3) == round(0.5,3))

pvalues <- results_Binomial[[index_stability]]$p_value

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

genes_ben_hoch <- gene_data$gene[est_ess_ben_hoch==1]


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

genes_bon_holm <- gene_data$gene[est_ess_bon_holm==1]


est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))

genes_bon <- gene_data$gene[est_ess_bon==1]

selected_genes_binomial <- 
  list(
    genes_bon=genes_bon,
    genes_bon_holm=genes_bon_holm,
    genes_ben_hoch=genes_ben_hoch
  )


# Selected Genes by ConNIS
# from generate_stability_tables.R
index_stability <- which(round(ws,3) == round(0.15,3))

pvalues <- results_ConNIS[[index_stability]]$p_value

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

genes_ben_hoch <- gene_data$gene[est_ess_ben_hoch==1]


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

genes_bon_holm <- gene_data$gene[est_ess_bon_holm==1]


est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))

genes_bon <- gene_data$gene[est_ess_bon==1]

selected_genes_connis <- 
  list(
    genes_bon=genes_bon,
    genes_bon_holm=genes_bon_holm,
    genes_ben_hoch=genes_ben_hoch
  )




# Selected Genes by Exp vs Gamma
# from generate_stability_tables.R
index_stability <- which(round(ts,3) == round(9,3))

est_ess <- results_ExpVsGamma[[index_stability]]$essential
selected_genes_expvsgamma <- gene_data$gene[est_ess==1]



# Selected Genes by Geometric
# from generate_stability_tables.R
index_stability <- which(round(ws,3) == round(0.4,3))

pvalues <- results_Geometric[[index_stability]]$p_value

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

genes_ben_hoch <- gene_data$gene[est_ess_ben_hoch==1]


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

genes_bon_holm <- gene_data$gene[est_ess_bon_holm==1]


est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))

genes_bon <- gene_data$gene[est_ess_bon==1]

selected_genes_geometric <- 
  list(
    genes_bon=genes_bon,
    genes_bon_holm=genes_bon_holm,
    genes_ben_hoch=genes_ben_hoch
  )



# Selected Genes by InsDens
# from generate_stability_tables.R
index_stability <- which(round(rs,3) == round(0.3,3))

est_ess <- results_InsDens[[index_stability]]$essential
selected_genes_insdens <- gene_data$gene[est_ess==1]


# Selected Genes by Tn5Gaps
# from generate_stability_tables.R
index_stability <- which(round(ws,3) == round(0.85,3))

pvalues <- results_Tn5Gaps[[index_stability]]$p_value

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

genes_ben_hoch <- gene_data$gene[est_ess_ben_hoch==1]


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

genes_bon_holm <- gene_data$gene[est_ess_bon_holm==1]


est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))

genes_bon <- gene_data$gene[est_ess_bon==1]

selected_genes_tn5gaps <- 
  list(
    genes_bon=genes_bon,
    genes_bon_holm=genes_bon_holm,
    genes_ben_hoch=genes_ben_hoch
  )

saveRDS(
  list(Binomial = selected_genes_binomial,
       ConNIS = selected_genes_connis,
       ExpVsGamma = selected_genes_expvsgamma,
       Geometric = selected_genes_geometric,
       InsDens = selected_genes_insdens,   
       Tn5Gaps = selected_genes_tn5gaps),
  file="results/selected_genes_tuned_EcoliBW25113.RDS")

genes_found_union_methods <- unique(c(
  selected_genes_binomial[[3]],
  selected_genes_connis[[3]],
  selected_genes_expvsgamma,
  selected_genes_geometric[[3]],
  selected_genes_insdens,   
  selected_genes_tn5gaps[[3]]
))

genes_sorted <- gene_data$gene[gene_data$gene %in% genes_found_union_methods]

comparison_genes_found <- 
  tibble(gene = genes_sorted,
         Binomial = genes_sorted %in% selected_genes_binomial[[3]],
         ConNIS = genes_sorted %in% selected_genes_connis[[1]],
         ExpVsGamma = genes_sorted %in% selected_genes_expvsgamma,
         Geometric = genes_sorted %in% selected_genes_geometric[[1]],
         Tn5Gaps = genes_sorted %in% selected_genes_tn5gaps[[3]])




write.csv2(comparison_genes_found, file = "table_selected_genes_tuned_EcoliBW25113.csv")

which(comparison_genes_found$ConNIS == TRUE & 
        comparison_genes_found$Binomial == FALSE &
        comparison_genes_found$ExpVsGamma == FALSE &
        comparison_genes_found$Geometric == FALSE &
        comparison_genes_found$Tn5Gaps == FALSE)

# nadD only in ConNIS, an essential gene  https://journals.asm.org/doi/10.1128/mbio.00747-13#:~:text=The%20respective%20enzymes%2C%20NaMN%20adenylyltransferase,tuberculosis%20(23–25).


venn_connis_vs_binomial <- ggvenn(list(
  ConNIS = selected_genes_connis[[3]],
  Binomial = selected_genes_binomial[[3]],
  Reference = true_ess_names), 
  show_percentage = FALSE, 
  text_size = 7, 
  stroke_size = 1.5, 
  set_name_size = 7,
  fill_color = c("#117733", "#6699CC", "white"))

venn_connis_vs_expvsgamma <- ggvenn(list(
  ConNIS = selected_genes_connis[[3]],
  "Exp. vs Gamma" = selected_genes_expvsgamma,
  Reference = true_ess_names), 
  show_percentage = FALSE, 
  text_size = 7, 
  stroke_size = 1.5, 
  set_name_size = 7,
  fill_color = c("#117733", "#CC6677", "white"))

venn_connis_vs_geometric <- ggvenn(list(
  ConNIS = selected_genes_connis[[3]],
  Geometric = selected_genes_geometric[[3]],
  Reference = true_ess_names), 
  show_percentage = FALSE, 
  text_size = 7, 
  stroke_size = 1.5, 
  set_name_size = 7,
  fill_color = c("#117733", "#7f7f7f", "white"))

venn_connis_vs_insdens <- ggvenn(list(
  ConNIS = selected_genes_connis[[3]],
  InsDens = selected_genes_insdens,
  Reference = true_ess_names), 
  show_percentage = FALSE, 
  text_size = 7, 
  stroke_size = 1.5, 
  set_name_size = 7,
  fill_color = c("#117733", "#DDCC77", "white"))

venn_connis_vs_tn5gaps <- ggvenn(list(
  ConNIS = selected_genes_connis[[3]],
  Tn5Gaps = selected_genes_tn5gaps[[1]],
  Reference = true_ess_names), 
  show_percentage = FALSE, 
  text_size = 7, 
  stroke_size = 1.5, 
  set_name_size = 7,
  fill_color = c("#117733", "#9651A0", "white"))

output_plot <- ggarrange(venn_connis_vs_binomial, 
                         venn_connis_vs_expvsgamma, 
                         venn_connis_vs_geometric,
                         venn_connis_vs_insdens,
                         venn_connis_vs_tn5gaps,
                         nrow = 1, ncol=5,
                         align='v', 
                         labels=c("A", "B", "C", "D", "E"),
                         font.label = list(size = 3,color= "#525252")
) 

save_plot(filename = "./plots/venn_connis_vs_all.pdf",
          plot =output_plot, dpi =600, base_height = 5, base_asp = 5)


gene_sets <- list(
  Binomial = selected_genes_binomial[[3]],
  ConNIS = selected_genes_connis[[3]],
  ExpVsGamma = selected_genes_expvsgamma,
  Geometric = selected_genes_geometric[[3]],
  InsDens = selected_genes_insdens,
  Tn5Gaps = selected_genes_tn5gaps[[3]]
)

method_colors <- 
  c(Binomial = "#6699CC", 
    ConNIS = "#117733", 
    ExpVsGamma = "#CC6677", 
    Geometric = "#7f7f7f", 
    InsDens = "#DDCC77",  
    Tn5Gaps = "#9651A0")

p_venns <- plot_pairwise_venn_matrix(gene_sets, method_colors, show_elements = FALSE,
                                     global_set   = true_ess_names,
                                     global_name  = "True Gene Set",
                                     global_color = "white",
                                     show_percentage = FALSE,
                                     venn_set_name_size = 1.8,
                                     venn_text_size     = 1.8,
                                     label_size         = 2.5,
                                     label_row_rel      = 0.08,  # relative Höhe der Kopfzeile
                                     label_col_rel      = 0.08,   # relative Breite der linken Spalte
                                     inner_margin = 0.1
)

save_plot(filename = "./plots/venns_EcoliBW25113.pdf",
          plot =p_venns, dpi =600, base_height = 10, base_asp = 1)






######## ECOLI MG1655 ############


# https://journals.asm.org/doi/10.1128/aac.00452-24
# https://genomics.lbl.gov/supplemental/rbarseq/html/Keio/


source("functions.R")



######## E. Coli #######

### Ma et al. data
all_genes_MG1655 <- read_csv("./MG1655_data/MG1655_Input.tradis_gene_insert_sites.csv.all.csv")

colnames(all_genes_MG1655)[c(2,4,5)] <-
  c("gene", "gene_start", "gene_end" )



# read insertion data from Ma et al. (2024)
# Original data from the authors
all_IS_raw <- read.csv("./MG1655_data/MG1655input.insert_site_plot", 
                       header = F, sep = " ")

# true essentials of PEC (Profiling of E. Coli Chromosome) data 
# https://shigen.nig.ac.jp/ecoli/pec/
true_ess <- read_table("./MG1655_data/M1655_PEC.txt")
true_ess <- as.character(unlist(true_ess[,2]))
true_ess <- true_ess[true_ess %in% all_genes_MG1655$gene]
true_ess_names <- true_ess

### noise thresholds
read_count_threshold <- 0

#### trimming
trimming_start <- 0.05
trimming_end <- 0.05

alpha_value <- 0.05
ws <- ws <- c(0.01, seq(0.05,1,0.05))
ts <- c(seq(0.1,2, 0.1),3:12)
rs <- c(seq(0.01,0.1, 0.01), seq(0.2,0.9,0.1), seq(0.91,0.99,0.01))

# combine + and - strang and set noise threshold
all_IS <- tibble(pos=which(c(all_IS_raw[,1]) + c(all_IS_raw[,2]) > read_count_threshold))


gene_data <- lapply(1:nrow(all_genes_MG1655), function(i){
  
  gene_i_start <- all_genes_MG1655$gene_start[i] + 
    floor(trimming_start * (all_genes_MG1655$gene_end[i] - all_genes_MG1655$gene_start[i]+1))
  
  gene_i_end <- all_genes_MG1655$gene_end[i] - 
    floor(trimming_end * (all_genes_MG1655$gene_end[i] - all_genes_MG1655$gene_start[i]+1))
  
  IS_gene_i <- 
    all_IS$pos[which(gene_i_start <= all_IS$pos &
                       gene_i_end >= all_IS$pos)]
  
  num_ins <- length(IS_gene_i)
  
  gene <- all_genes_MG1655$gene[i]
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
  
  tibble(entry = all_genes_MG1655$locusId[i],
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
           genome.length = 4641652, 
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
         genome.length = 4641652, 
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
            genome.length = 4641652, 
            weighting = w)
  
})

# Tn5Gaps

ins_sites <- 
  sort(unique(
    c(1, 
      unique(all_IS$pos), 
      4641652)))

print("tn5")
results_Tn5Gaps <- lapply(ws, function(w){
  
  Tn5Gaps(ins.positions = ins_sites, 
          gene.names = gene_data$gene, 
          gene.starts = gene_data$gene_start, 
          gene.stops = gene_data$gene_end, 
          genome.length = 4641652, 
          weighting = w)
  
})


# Selected Genes by Binomial
# from generate_stability_tables.R
index_stability <- which(round(ws,3) == round(0.45,3))

pvalues <- results_Binomial[[index_stability]]$p_value

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

genes_ben_hoch <- gene_data$gene[est_ess_ben_hoch==1]


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

genes_bon_holm <- gene_data$gene[est_ess_bon_holm==1]


est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))

genes_bon <- gene_data$gene[est_ess_bon==1]

selected_genes_binomial <- 
  list(
    genes_bon=genes_bon,
    genes_bon_holm=genes_bon_holm,
    genes_ben_hoch=genes_ben_hoch
  )


# Selected Genes by ConNIS
# from generate_stability_tables.R
index_stability <- which(round(ws,3) == round(0.25,3))

pvalues <- results_ConNIS[[index_stability]]$p_value

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

genes_ben_hoch <- gene_data$gene[est_ess_ben_hoch==1]


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

genes_bon_holm <- gene_data$gene[est_ess_bon_holm==1]


est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))

genes_bon <- gene_data$gene[est_ess_bon==1]

selected_genes_connis <- 
  list(
    genes_bon=genes_bon,
    genes_bon_holm=genes_bon_holm,
    genes_ben_hoch=genes_ben_hoch
  )




# Selected Genes by Exp vs Gamma
# from generate_stability_tables.R
index_stability <- which(round(ts,3) == round(1,3))

est_ess <- results_ExpVsGamma[[index_stability]]$essential
selected_genes_expvsgamma <- gene_data$gene[est_ess==1]



# Selected Genes by Geometric
# from generate_stability_tables.R
index_stability <- which(round(ws,3) == round(0.35,3))

pvalues <- results_Geometric[[index_stability]]$p_value

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

genes_ben_hoch <- gene_data$gene[est_ess_ben_hoch==1]


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

genes_bon_holm <- gene_data$gene[est_ess_bon_holm==1]


est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))

genes_bon <- gene_data$gene[est_ess_bon==1]

selected_genes_geometric <- 
  list(
    genes_bon=genes_bon,
    genes_bon_holm=genes_bon_holm,
    genes_ben_hoch=genes_ben_hoch
  )



# Selected Genes by InsDens
# from generate_stability_tables.R
index_stability <- which(round(rs,3) == round(0.1,3))

est_ess <- results_InsDens[[index_stability]]$essential
selected_genes_insdens <- gene_data$gene[est_ess==1]


# Selected Genes by Tn5Gaps
# from generate_stability_tables.R
index_stability <- which(round(ws,3) == round(0.85,3))

pvalues <- results_Tn5Gaps[[index_stability]]$p_value

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

genes_ben_hoch <- gene_data$gene[est_ess_ben_hoch==1]


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

genes_bon_holm <- gene_data$gene[est_ess_bon_holm==1]


est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))

genes_bon <- gene_data$gene[est_ess_bon==1]

selected_genes_tn5gaps <- 
  list(
    genes_bon=genes_bon,
    genes_bon_holm=genes_bon_holm,
    genes_ben_hoch=genes_ben_hoch
  )

saveRDS(
  list(Binomial = selected_genes_binomial,
       ConNIS = selected_genes_connis,
       ExpVsGamma = selected_genes_expvsgamma,
       Geometric = selected_genes_geometric,
       InsDens = selected_genes_insdens,   
       Tn5Gaps = selected_genes_tn5gaps),
  file="results/selected_genes_tuned_EcoliMG1655.RDS")

genes_found_union_methods <- unique(c(
  selected_genes_binomial[[3]],
  selected_genes_connis[[3]],
  selected_genes_expvsgamma,
  selected_genes_geometric[[3]],
  selected_genes_insdens,   
  selected_genes_tn5gaps[[3]]
))

genes_sorted <- gene_data$gene[gene_data$gene %in% genes_found_union_methods]

comparison_genes_found <- 
  tibble(gene = genes_sorted,
         Binomial = genes_sorted %in% selected_genes_binomial[[3]],
         ConNIS = genes_sorted %in% selected_genes_connis[[1]],
         ExpVsGamma = genes_sorted %in% selected_genes_expvsgamma,
         Geometric = genes_sorted %in% selected_genes_geometric[[1]],
         Tn5Gaps = genes_sorted %in% selected_genes_tn5gaps[[3]])




write.csv2(comparison_genes_found, file = "table_selected_genes_tuned_EcoliMG1655.csv")

which(comparison_genes_found$ConNIS == TRUE & 
        comparison_genes_found$Binomial == FALSE &
        comparison_genes_found$ExpVsGamma == FALSE &
        comparison_genes_found$Geometric == FALSE &
        comparison_genes_found$Tn5Gaps == FALSE)


gene_sets <- list(
  Binomial = selected_genes_binomial[[3]],
  ConNIS = selected_genes_connis[[3]],
  ExpVsGamma = selected_genes_expvsgamma,
  Geometric = selected_genes_geometric[[3]],
  InsDens = selected_genes_insdens,
  Tn5Gaps = selected_genes_tn5gaps[[3]]
)

method_colors <- 
  c(Binomial = "#6699CC", 
    ConNIS = "#117733", 
    ExpVsGamma = "#CC6677", 
    Geometric = "#7f7f7f", 
    InsDens = "#DDCC77",  
    Tn5Gaps = "#9651A0")

p_venns <- plot_pairwise_venn_matrix(gene_sets, method_colors, show_elements = FALSE,
                                     global_set   = true_ess_names,
                                     global_name  = "True Gene Set",
                                     global_color = "white",
                                     show_percentage = FALSE,
                                     venn_set_name_size = 1.8,
                                     venn_text_size     = 1.8,
                                     label_size         = 2.5,
                                     label_row_rel      = 0.08,  # relative Höhe der Kopfzeile
                                     label_col_rel      = 0.08,   # relative Breite der linken Spalte
                                     inner_margin = 0.1
)

save_plot(filename = "./plots/venns_EcoliMG1655.pdf",
          plot =p_venns, dpi =600, base_height = 10, base_asp = 1)









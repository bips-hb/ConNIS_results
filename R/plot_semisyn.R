library(tidyverse)
library(ggpubr)
library(cowplot)

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

# use subsample sizes of semi-synthetic simulation
subsample_sizes <- c(50000,100000,200000,400000)

is_subsamples <- lapply(subsample_sizes, function(subsample_size){
  set.seed(subsample_size)
  sort(sample(is_pos, subsample_size, replace = FALSE))
})

gene_wise_insertion_densities_per_subsetsize <-
  lapply(seq_along(subsample_sizes), function(j){
    ins_dense_per_gene_subsample <- sapply(1:nrow(gene_data), function(i) {
      sum(is_subsamples[[j]] >= gene_data$start[i] &
            is_subsamples[[j]] <= gene_data$end[i])/
        (gene_data$end[i]-gene_data$start[i]+1)
    })
  })


# load references essetial genes based on ExpVsGamma and the IS of Goodall, 2018
ref_ess_gene <- readRDS("./bw25113_data/ref_ess_semisyn.RDS")

# load results of semi synthetic
semi_syn_results <- readRDS(file="./performance/Results_semisyn_BW25113.RDS")

# load performances of semi syntetic
semi_syn_performances <- readRDS(file="./performance/Performance_semisyn_BW25113.RDS")

# best MCC of each method per subsample size
best_MCCs <- semi_syn_performances %>%
  group_by(method, subsample_size) %>%
  slice_max(MCC)

# create an index number for the subsample sizes
best_MCCs$subsample_index_number <- rep(1:4, length(unique(best_MCCs$method)))

# semi_syn_results is a list of 4 lists (each list for one subsample sizes),
# each list containing 4 lists (for four methods), each containing 10
# lists, each list representing the results of one weighting value

plots <- lapply(1:nrow(best_MCCs), function(i){

  all_results_per_subsample_and_method <-
    semi_syn_results[[
      best_MCCs[i,]$subsample_index_number]][[
        paste("results_", best_MCCs[i,]$method, sep="")]]

  all_results_per_subsample_and_method <-
    do.call(rbind, all_results_per_subsample_and_method)

  results_optimal_tuning <-
    all_results_per_subsample_and_method[
      all_results_per_subsample_and_method[,3] == best_MCCs[i,]$tuning,]

  if(best_MCCs[i,]$method %in% c("ExpVsGamma","InsDens")){
    est_ess_genes <- results_optimal_tuning$gene[
      results_optimal_tuning$essential ==1 ]
    tuning <- NULL
    if(best_MCCs[i,]$method == c("ExpVsGamma")){
      tuning_type <- "log2 threshold"
    }else{
      tuning_type <- "posterior probability"
    }

  }else{
    est_ess_genes <- results_optimal_tuning$gene[
      p.adjust(results_optimal_tuning$p_value, method = "hochberg") <= 0.05]
    tuning <- best_MCCs[i,]$tuning
    tuning_type <- "weight value"
  }

  TP_genes <- est_ess_genes[est_ess_genes %in% ref_ess_gene]

  FP_genes <- est_ess_genes[!(est_ess_genes %in% ref_ess_gene)]

  FN_genes <- ref_ess_gene[!(ref_ess_gene %in% est_ess_genes)]

  color_tps <- rep(NA, nrow(gene_data))
  color_tps[which(gene_data$gene %in% TP_genes)]  <- "#2CA68F"

  color_fps <- rep(NA, nrow(gene_data))
  color_fps[which(gene_data$gene %in% FP_genes)]  <- "#FF348F"

  color_fns <- rep(NA, nrow(gene_data))
  color_fns[which(gene_data$gene %in% FN_genes)]  <- "#7A6CB0"

  plot_data <-
    bind_cols(gene_data,
              is_density = gene_wise_insertion_densities_per_subsetsize[[best_MCCs[i,]$subsample_index_number]]
    )

  p <- ggplot(plot_data,
         aes(x=start, y=is_density))+
    geom_point(color="#B8B8B8", size=0.8) +
    geom_hline(yintercept = mean(plot_data$is_density),
               color="#8C8C8C", size=1)+
    geom_hline(yintercept = mean(plot_data$is_density)*tuning,
               color="#414141", size=1, linetype="dashed")+
    geom_point(data=plot_data,
               aes(x=start,y=is_density), color=color_fps, size=0.8)+
    geom_point(data=plot_data,
               aes(x=start,y=is_density), color=color_fns, size=0.8)+
    geom_point(data=plot_data,
               aes(x=start,y=is_density), color=color_tps, size=0.8)+
    annotate("rect", xmin = 3390000, xmax = 4550000,
             ymin = max(plot_data$is_density)*0.725,
             ymax = max(plot_data$is_density)*0.925,
             fill="white")+
    annotate("text", x = 4000000, y=max(plot_data$is_density)*0.9, size=8/.pt,
             label = paste0("MCC: ", round(best_MCCs[i,]$MCC, 2)), color = "#4A4A4A")+
    annotate("text", x = 4000000, y=max(plot_data$is_density)*0.85, size=8/.pt,
             label = paste0("TP: ", length(TP_genes)), color = "#2CA68F")+
    annotate("text", x = 4000000, y=max(plot_data$is_density)*0.8, size=8/.pt,
             label = paste0("FP: ", length(FP_genes)), color = "#FF348F")+
    annotate("text", x = 4000000, y=max(plot_data$is_density)*0.75, size=8/.pt,
             label = paste0("FN: ", length(FN_genes)), color = "#7A6CB0")+
    xlab("Genomic postion")+
    ylab("Gene-wise insertion density")+
    ggtitle(paste("Subsample size: ", as.integer(best_MCCs[i,]$subsample_size),
                    ", ", tuning_type,": ", best_MCCs[i,]$tuning, sep="")) +
    theme_minimal() +
    theme(legend.position="bottom",
          plot.title=element_text( face='italic', size=6, hjust = 0.5),
          legend.title=element_blank(),
          legend.text=element_text(size=rel(0.6)),
          strip.text = element_text(size=rel(0.6)),
          axis.text.x=element_text(size = rel(0.8)),
          axis.text.y=element_text(size = rel(0.8)),
          axis.title.x = element_text(size = rel(0.6)),
          axis.title.y = element_text(size = rel(0.6)))

  p

})

#plots for Binomial
plots_binomial <- ggarrange(
  plots[[1]],
  plots[[2]],
  plots[[3]],
  plots[[4]],
  nrow = 2, ncol=2,
  labels=LETTERS[1:4],
  font.label = list(size = 12,color= "#525252"),
  align='v'
)

plots_binomial <- annotate_figure(plots_binomial,
                top = text_grob("Binomial", size = 14))

save_plot(filename = "./plots/semisyn_performance_Binomial.pdf",
          plot =plots_binomial, dpi =600, base_height = 6, base_asp = 4/4 )


#plots for ConNIS
plots_connis <- ggarrange(
  plots[[5]],
  plots[[6]],
  plots[[7]],
  plots[[8]],
  nrow = 2, ncol=2,
  labels=LETTERS[1:4],
  font.label = list(size = 12,color= "#525252"),
  align='v'
)

plots_connis <- annotate_figure(plots_connis,
                                  top = text_grob("ConNIS", size = 14))

save_plot(filename = "./plots/semisyn_performance_ConNIS.pdf",
          plot =plots_connis, dpi =600, base_height = 6, base_asp = 4/4 )

#plots for Exp vs. Gamma
plots_expvsgamma <- ggarrange(
  plots[[9]],
  plots[[10]],
  plots[[11]],
  plots[[12]],
  nrow = 2, ncol=2,
  labels=LETTERS[1:4],
  font.label = list(size = 12,color= "#525252"),
  align='v'
)

plots_expvsgamma <- annotate_figure(plots_expvsgamma,
                                    top = text_grob("Exp. vs. Gamma", size = 14))

save_plot(filename = "./plots/semisyn_performance_ExpVsGamma.pdf",
          plot =plots_expvsgamma, dpi =600, base_height = 6, base_asp = 4/4 )


#plots for Geometric
plots_geometric <- ggarrange(
  plots[[13]],
  plots[[14]],
  plots[[15]],
  plots[[16]],
  nrow = 2, ncol=2,
  labels=LETTERS[1:4],
  font.label = list(size = 12,color= "#525252"),
  align='v'
)

plots_geometric <- annotate_figure(plots_geometric,
                                  top = text_grob("Geometric", size = 14))

save_plot(filename = "./plots/semisyn_performance_Geometric.pdf",
          plot =plots_geometric, dpi =600, base_height = 6, base_asp = 4/4 )


#plots for InsDens
plots_insdens <- ggarrange(
  plots[[17]],
  plots[[18]],
  plots[[19]],
  plots[[20]],
  nrow = 2, ncol=2,
  labels=LETTERS[1:4],
  font.label = list(size = 12,color= "#525252"),
  align='v'
)

plots_insdens <- annotate_figure(plots_insdens,
                                  top = text_grob("InsDens", size = 14))
save_plot(filename = "./plots/semisyn_performance_InsDens.pdf",
          plot =plots_insdens, dpi =600, base_height = 6, base_asp = 4/4 )

#plots for Tn5Gaps
plots_tn5gaps <- ggarrange(
  plots[[21]],
  plots[[22]],
  plots[[23]],
  plots[[24]],
  nrow = 2, ncol=2,
  labels=LETTERS[1:4],
  font.label = list(size = 12,color= "#525252"),
  align='v'
)

plots_tn5gaps <- annotate_figure(plots_tn5gaps,
                                 top = text_grob("Tn5Gaps", size = 14))
save_plot(filename = "./plots/semisyn_performance_Tn5Gaps.pdf",
          plot =plots_tn5gaps, dpi =600, base_height = 6, base_asp = 4/4 )



plots_comparison <- lapply(subsample_sizes, function(i){
  p <- ggplot(semi_syn_performances %>% filter(subsample_size==i),
         aes(x=TP+FP, y=MCC, color=method)) +
    geom_vline(xintercept = length(ref_ess_gene), linetype="dashed", alpha = 0.7, color="orange") +
    scale_color_manual(values =
                         c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
    geom_line(size=1) +
    xlab("Number of genes labeled as 'essential'") +
    theme_minimal() +
    ggtitle(paste("Subsample size: ", as.integer(i), sep="")) +
    theme(legend.position="bottom",
          plot.title=element_text( face='italic', size=6, hjust = 0.5),
          legend.title=element_blank(),
          legend.text=element_text(size=rel(0.6)),
          strip.text = element_text(size=rel(0.6)),
          axis.text.x=element_text(size = rel(0.8)),
          axis.text.y=element_text(size = rel(0.8)),
          axis.title.x = element_text(size = rel(0.6)),
          axis.title.y = element_text(size = rel(0.6))) +
    coord_cartesian(ylim = c(0, 1))

    p
})

p_comp_subset_sizes <- ggarrange(
  plots_comparison[[1]],
  plots_comparison[[2]],
  plots_comparison[[3]],
  plots_comparison[[4]],
  nrow = 2, ncol=2,
  labels=LETTERS[1:4],
  font.label = list(size = 12,color= "#525252"),
  align='v',
  common.legend = T,
  legend = "bottom"
)

save_plot(filename = "./plots/semisyn_comparison.pdf",
          plot =p_comp_subset_sizes, dpi =600, base_height = 6, base_asp = 4/4 )

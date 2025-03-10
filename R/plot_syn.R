# Plots for synthetic settings
## Major Synthetic Setting 1 A)
input <- list()
input$unique_loci_syn <- 200000
input$lambda_syn <- 0.75
input$wave_size_syn <- 2000000
input$tecnoise_syn <- 0
input$trim_value <- 0.1
input$metric <- "MCC"

results_sim <- lapply(c("Binomial", "ConNIS", "ExpGamma", "Geometric", "InsDens", "Tn5Gaps"), function(m){
  readRDS(
    paste("./performance/Performance_",
          m, 
          "_ess_ORF_", "uniform",#input$ess_ORF,
          "_numLoci_", format(input$unique_loci_syn, scientific = F),
          "_lambda_",input$lambda_syn,
          "_bp_per_wave_", "2e+06",
          "_sine_scaling_factor_", "1.3",
          "_NBnc_", 100, #input$NegBinomial_num_cluster,
          "_NBdis_", 1, #input$NegBinomial_dispersion,
          "_NBp_", 0.3, #input$NegBinomial_p,
          "_tecnoise_", input$tecnoise_syn,
          "_trimmed_", input$trim_value,
          ".RDS", sep="")
  )
})
results_sim <- do.call(bind_rows, results_sim)


results_to_analyze <-
  results_sim %>%
  group_by(method, 
           sim.number,
           weighting,
           threshold_log2_expgamma,
           post_prob,
           adjustment,
           alpha_value
  ) %>%
  mutate(tuning = sum(
    c(
      threshold_log2_expgamma,
      post_prob,
      weighting
    ), na.rm=T) 
  )

get_optimal_tuning <- results_to_analyze %>% 
  group_by(method, tuning, adjustment) %>% 
  summarise(mean_Metric = mean(!!sym(input$metric), na.rm = T)) %>% 
  group_by(method) %>% 
  filter(mean_Metric==max(mean_Metric, na.rm = T)) %>% 
  distinct(method, mean_Metric, .keep_all = TRUE)

x_min <- 0
x_max <- 1

best_results <- lapply(1:nrow(get_optimal_tuning), function(i){
  results_to_analyze %>% 
    filter(method == get_optimal_tuning$method[i] & 
             tuning == get_optimal_tuning$tuning[i] &
             (adjustment == get_optimal_tuning$adjustment[i] | 
                is.na(adjustment)))
})

best_results <- do.call(rbind, best_results)

best_results$adjustment_short <- NULL

best_results <- best_results %>% mutate(adjustment_short = case_when(
  adjustment == "Benjamini-Hochberg" ~ "Benj.-Hoch.",
  adjustment == "Bonferroni-Holm" ~ "Bonf.-Holm",
  adjustment == "Bonferroni" ~ "Bonf."
))

best_results <- best_results %>%
  arrange(desc(method))

best_results <- 
  melt(data = best_results, 
       id.vars = "method",
       measure.vars =  c("MCC", "Precision", "TPR"), 
       variable.name = "Metric")

unweigthed_results_no_melt <- results_to_analyze %>% 
  group_by(method) %>% 
  filter((tuning == 1 & method == "Binomial" & adjustment == "Benjamini-Hochberg") |
           (tuning == 1 & method == "ConNIS" & adjustment == "Bonferroni-Holm") |
           (threshold_log2_expgamma == 2 & method == "Exp. vs. Gamma" ) |
           (tuning == 1 & method == "Geometric" & adjustment == "Bonferroni-Holm") |
           (post_prob == 0.1 & method == "InsDens") |
           (tuning == 1 & method == "Tn5Gaps" & adjustment == "Benjamini-Hochberg")
  )

unweigthed_results <- unweigthed_results_no_melt %>%
  arrange(desc(method))

unweigthed_results <- 
  melt(data = unweigthed_results, 
       id.vars = "method",
       measure.vars =  c("MCC", "Precision", "TPR"), 
       variable.name = "Metric")

df_scales <- data.frame(
  Panel = c("MCC", "Precision", "TPR"),
  xmin = c(0, 0, 0),
  xmax = c(1, 1, 1),
  n = c(4, 4, 4)
)
df_scales <- split(df_scales, df_scales$Panel)

scales <- lapply(df_scales, function(x) {
  scale_x_continuous(limits = c(x$xmin, x$xmax), n.breaks = x$n)
})

p1_rl_max_mcc_precision_tpr <- ggplot(
  best_results,
  aes(y=fct_inorder(method), x= value, fill=method)) + 
  # geom_vline(xintercept = 0.5, color="black",linetype="dotted", alpha=0.5)+
  geom_vline(xintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
  geom_density_ridges2(alpha=0.5, quantile_lines=TRUE, 
                       quantile_fun=function(MCC,...){mean(MCC)}, 
                       rel_min_height = 0.005 ) +
  geom_density_ridges2(alpha=0, quantile_lines=TRUE,linewidth = 0.5, linetype = "dashed",
                       quantile_fun=function(MCC,...){mean(MCC)},
                       rel_min_height = 0.005, data=unweigthed_results, aes(y=fct_inorder(method), x= value, color=method)) +
  facet_grid(~Metric , scales = "free", labeller = as_labeller(c("MCC" = "MCC", "Precision" = "Corresponding Precision", "TPR" = "Corresponding TPR"))) +
  ggh4x::facetted_pos_scales(
    x = scales
  )+
  theme(legend.position="none") +
  theme(axis.text.x=element_text(size=rel(1)),
        axis.text.y=element_text(size=rel(1)),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  scale_color_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))


best_results_mcc <- best_results %>% filter(Metric == "MCC")
unweigthed_results_mcc <- unweigthed_results %>% filter(Metric == "MCC")
p1_rl_max_mcc <- ggplot(
  best_results_mcc,
  aes(y=fct_inorder(method), x= value, fill=method)) + 
  # geom_vline(xintercept = 0.5, color="black",linetype="dotted", alpha=0.5)+
  geom_vline(xintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  # geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
  geom_density_ridges2(alpha=0.5, quantile_lines=TRUE, 
                       quantile_fun=function(MCC,...){mean(MCC)}, 
                       rel_min_height = 0.005 ) +
  geom_density_ridges2(alpha=0, quantile_lines=TRUE,linewidth = 0.5, linetype = "dashed",
                       quantile_fun=function(MCC,...){mean(MCC)},
                       rel_min_height = 0.005, data=unweigthed_results_mcc, aes(y=fct_inorder(method), x= value, color=method)) +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6), angle = 0, vjust = 0, hjust=1),
        axis.title.y=element_blank())+
  scale_fill_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  scale_color_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  xlab("MCC")+
  xlim(0,1)+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS1: 200,000 IS and min. 75% insertion-free")



plot_data <- results_to_analyze %>%
  filter(
    (method == "Binomial" & adjustment == "Benjamini-Hochberg") |
      (method == "ConNIS" & adjustment == "Bonferroni-Holm") |
      (method == "Exp. vs. Gamma" & is.na(adjustment)) |
      (method == "Geometric" & adjustment == "Bonferroni-Holm") |
      (method == "InsDens" & is.na(adjustment)) |
      (method == "Tn5Gaps" & adjustment == "Benjamini-Hochberg")
  )

y_min <- 0
y_max <- 1


# data_unweigted_range <- 
#   bind_cols(unweigthed_results_no_melt %>% 
#               group_by(method) %>% 
#               mutate(N_min=TP+FP) %>% 
#               slice(which.min(N_min)) %>%  
#               select(MCC_min = MCC, N_min), 
#             unweigthed_results %>% 
#               group_by(method) %>% 
#               mutate(N_max=TP+FP) %>% 
#               slice(which.max(N_max)) %>%  
#               select(MCC_max=MCC, N_max) ) %>% 
#   select(method=method...1, MCC_min, MCC_max, N_min, N_max)



p1_setsize <- ggplot(plot_data %>%mutate(N=TP+FP), 
                     aes(x=N, y=!!sym(input$metric), color=method)) +
  # geom_rect(data = unweigthed_results_no_melt,inherit.aes = FALSE,
  #           aes(fill=method, xmin = N_min, ymin = 0, xmax = N_max, ymax = 1), alpha=0.3)+
  # scale_fill_manual(values =
  #                     c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_hline(yintercept = 0, color="black",linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = mean(plot_data$TP+plot_data$FN), color="orange",linetype="dashed", alpha=0.7)+
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE, size = 0.75) +
  xlab("Number of genes labeled 'essential'")+
  coord_cartesian(ylim=c(y_min,y_max))+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS1: 200,000 IS and min. 75% insertion-free") +
  geom_point(data = unweigthed_results_no_melt %>%
               group_by(method) %>% reframe(MCC=mean(MCC), N=mean(TP+FP)),
             aes(x=N, y=MCC), colour = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"), size=2)
  # geom_segment(data = unweigthed_results_no_melt, 
  #              aes(x = N_min, y = -0.05, xend = N_min, yend = 0.05), size=0.8
  #              )+
  # geom_segment(data = unweigthed_results_no_melt, 
  #              aes(x = N_max, y = -0.05, xend = N_max, yend = 0.05), size=0.8
  # )
  
  

p1_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE, size = 0.75) +
  xlab("Recall")+
  coord_cartesian(ylim=c(0,1), xlim=c(0,1))+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS1: 200,000 IS and min. 75% insertion-free")+
  geom_point(data = unweigthed_results_no_melt %>%
               group_by(method) %>% reframe(Precision=mean(Precision), TPR=mean(TPR)),
             aes(x=TPR, y=Precision), 
             colour = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"), size=2)


## Major Synthetic Setting 1 B)
input <- list()
input$unique_loci_syn <- 400000
input$lambda_syn <- 0.85
input$wave_size_syn <- 2000000
input$tecnoise_syn <- 0.02
input$trim_value <- 0.1
input$metric <- "MCC"

results_sim <- lapply(c("Binomial", "ConNIS", "ExpGamma", "Geometric", "InsDens", "Tn5Gaps"), function(m){
  readRDS(
    paste("./performance/Performance_",
          m, 
          "_ess_ORF_", "uniform",#input$ess_ORF,
          "_numLoci_", format(input$unique_loci_syn, scientific = F),
          "_lambda_",input$lambda_syn,
          "_bp_per_wave_", "2e+06",
          "_sine_scaling_factor_", "1.3",
          "_NBnc_", 100, #input$NegBinomial_num_cluster,
          "_NBdis_", 1, #input$NegBinomial_dispersion,
          "_NBp_", 0.3, #input$NegBinomial_p,
          "_tecnoise_", input$tecnoise_syn,
          "_trimmed_", input$trim_value,
          ".RDS", sep="")
  )
})
results_sim <- do.call(bind_rows, results_sim)


results_to_analyze <-
  results_sim %>%
  group_by(method, 
           sim.number,
           weighting,
           threshold_log2_expgamma,
           post_prob,
           adjustment,
           alpha_value
  ) %>%
  mutate(tuning = sum(
    c(
      threshold_log2_expgamma,
      post_prob,
      weighting
    ), na.rm=T) 
  )

get_optimal_tuning <- results_to_analyze %>% 
  group_by(method, tuning, adjustment) %>% 
  summarise(mean_Metric = mean(!!sym(input$metric), na.rm = T)) %>% 
  group_by(method) %>% 
  filter(mean_Metric==max(mean_Metric, na.rm = T)) %>% 
  distinct(method, mean_Metric, .keep_all = TRUE)

x_min <- 0
x_max <- 1

best_results <- lapply(1:nrow(get_optimal_tuning), function(i){
  results_to_analyze %>% 
    filter(method == get_optimal_tuning$method[i] & 
             tuning == get_optimal_tuning$tuning[i] &
             (adjustment == get_optimal_tuning$adjustment[i] | 
                is.na(adjustment)))
})

best_results <- do.call(rbind, best_results)

best_results$adjustment_short <- NULL

best_results <- best_results %>% mutate(adjustment_short = case_when(
  adjustment == "Benjamini-Hochberg" ~ "Benj.-Hoch.",
  adjustment == "Bonferroni-Holm" ~ "Bonf.-Holm",
  adjustment == "Bonferroni" ~ "Bonf."
))

best_results <- best_results %>%
  arrange(desc(method))

best_results <- 
  melt(data = best_results, 
       id.vars = "method",
       measure.vars =  c("MCC", "Precision", "TPR"), 
       variable.name = "Metric")

unweigthed_results_no_melt <- results_to_analyze %>% 
  group_by(method) %>% 
  filter((tuning == 1 & method == "Binomial" & adjustment == "Benjamini-Hochberg") |
           (tuning == 1 & method == "ConNIS" & adjustment == "Bonferroni-Holm") |
           (threshold_log2_expgamma == 2 & method == "Exp. vs. Gamma" ) |
           (tuning == 1 & method == "Geometric" & adjustment == "Bonferroni-Holm") |
           (post_prob == 0.1 & method == "InsDens") |
           (tuning == 1 & method == "Tn5Gaps" & adjustment == "Benjamini-Hochberg")
  )

unweigthed_results <- unweigthed_results_no_melt %>%
  arrange(desc(method))

unweigthed_results <- 
  melt(data = unweigthed_results, 
       id.vars = "method",
       measure.vars =  c("MCC", "Precision", "TPR"), 
       variable.name = "Metric")

df_scales <- data.frame(
  Panel = c("MCC", "Precision", "TPR"),
  xmin = c(0, 0, 0),
  xmax = c(1, 1, 1),
  n = c(4, 4, 4)
)
df_scales <- split(df_scales, df_scales$Panel)

scales <- lapply(df_scales, function(x) {
  scale_x_continuous(limits = c(x$xmin, x$xmax), n.breaks = x$n)
})

p2_rl_max_mcc_precision_tpr <- ggplot(
  best_results,
  aes(y=fct_inorder(method), x= value, fill=method)) + 
  # geom_vline(xintercept = 0.5, color="black",linetype="dotted", alpha=0.5)+
  geom_vline(xintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
  geom_density_ridges2(alpha=0.5, quantile_lines=TRUE, 
                       quantile_fun=function(MCC,...){mean(MCC)}, 
                       rel_min_height = 0.005 ) +
  geom_density_ridges2(alpha=0, quantile_lines=TRUE,linewidth = 0.5, linetype = "dashed",
                       quantile_fun=function(MCC,...){mean(MCC)},
                       rel_min_height = 0.005, data=unweigthed_results, aes(y=fct_inorder(method), x= value, color=method)) +
  facet_grid(~Metric , scales = "free", labeller = as_labeller(c("MCC" = "MCC", "Precision" = "Corresponding Precision", "TPR" = "Corresponding TPR"))) +
  ggh4x::facetted_pos_scales(
    x = scales
  )+
  theme(legend.position="none") +
  theme(axis.text.x=element_text(size=rel(1)),
        axis.text.y=element_text(size=rel(1)),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  scale_color_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))


best_results_mcc <- best_results %>% filter(Metric == "MCC")
unweigthed_results_mcc <- unweigthed_results %>% filter(Metric == "MCC")
p2_rl_max_mcc <- ggplot(
  best_results_mcc,
  aes(y=fct_inorder(method), x= value, fill=method)) + 
  # geom_vline(xintercept = 0.5, color="black",linetype="dotted", alpha=0.5)+
  geom_vline(xintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  # geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
  geom_density_ridges2(alpha=0.5, quantile_lines=TRUE, 
                       quantile_fun=function(MCC,...){mean(MCC)}, 
                       rel_min_height = 0.005 ) +
  geom_density_ridges2(alpha=0, quantile_lines=TRUE,linewidth = 0.5, linetype = "dashed",
                       quantile_fun=function(MCC,...){mean(MCC)},
                       rel_min_height = 0.005, data=unweigthed_results_mcc, aes(y=fct_inorder(method), x= value, color=method)) +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6), angle = 0, vjust = 0, hjust=1),
        axis.title.y=element_blank())+
  scale_fill_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  scale_color_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  xlab("MCC")+
  xlim(0,1)+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS2: 400,000 IS, min. 85% insertion-free and 2% noise IS")



plot_data <- results_to_analyze %>%
  filter(
    (method == "Binomial" & adjustment == "Benjamini-Hochberg") |
      (method == "ConNIS" & adjustment == "Bonferroni-Holm") |
      (method == "Exp. vs. Gamma" & is.na(adjustment)) |
      (method == "Geometric" & adjustment == "Bonferroni-Holm") |
      (method == "InsDens" & is.na(adjustment)) |
      (method == "Tn5Gaps" & adjustment == "Benjamini-Hochberg")
  )

y_min <- 0
y_max <- 1


# data_unweigted_range <- 
#   bind_cols(unweigthed_results_no_melt %>% 
#               group_by(method) %>% 
#               mutate(N_min=TP+FP) %>% 
#               slice(which.min(N_min)) %>%  
#               select(MCC_min = MCC, N_min), 
#             unweigthed_results %>% 
#               group_by(method) %>% 
#               mutate(N_max=TP+FP) %>% 
#               slice(which.max(N_max)) %>%  
#               select(MCC_max=MCC, N_max) ) %>% 
#   select(method=method...1, MCC_min, MCC_max, N_min, N_max)



p2_setsize <- ggplot(plot_data %>%mutate(N=TP+FP), 
                     aes(x=N, y=!!sym(input$metric), color=method)) +
  # geom_rect(data = unweigthed_results_no_melt,inherit.aes = FALSE,
  #           aes(fill=method, xmin = N_min, ymin = 0, xmax = N_max, ymax = 1), alpha=0.3)+
  # scale_fill_manual(values =
  #                     c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_hline(yintercept = 0, color="black",linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = mean(plot_data$TP+plot_data$FN), color="orange",linetype="dashed", alpha=0.7)+
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE, size = 0.75) +
  xlab("Number of genes labeled 'essential'")+
  coord_cartesian(ylim=c(y_min,y_max))+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS2: 400,000 IS, min. 85% insertion-free and 2% noise IS") +
  geom_point(data = unweigthed_results_no_melt %>%
               group_by(method) %>% reframe(MCC=mean(MCC), N=mean(TP+FP)),
             aes(x=N, y=MCC), colour = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"), size=2)
# geom_segment(data = unweigthed_results_no_melt, 
#              aes(x = N_min, y = -0.05, xend = N_min, yend = 0.05), size=0.8
#              )+
# geom_segment(data = unweigthed_results_no_melt, 
#              aes(x = N_max, y = -0.05, xend = N_max, yend = 0.05), size=0.8
# )



p2_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE, size = 0.75) +
  xlab("Recall")+
  coord_cartesian(ylim=c(0,1), xlim=c(0,1))+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS2: 400,000 IS, min. 85% insertion-free and 2% noise IS")+
  geom_point(data = unweigthed_results_no_melt %>%
               group_by(method) %>% reframe(Precision=mean(Precision), TPR=mean(TPR)),
             aes(x=TPR, y=Precision), 
             colour = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"), size=2)

## Major Synthetic Setting 1 C)
input <- list()
input$unique_loci_syn <- 50000
input$lambda_syn <- 0.85
input$wave_size_syn <- 2000000
input$tecnoise_syn <- 0
input$trim_value <- 0.1
input$metric <- "MCC"

results_sim <- lapply(c("Binomial", "ConNIS", "ExpGamma", "Geometric", "InsDens", "Tn5Gaps"), function(m){
  readRDS(
    paste("./performance/Performance_",
          m, 
          "_ess_ORF_", "uniform",#input$ess_ORF,
          "_numLoci_", format(input$unique_loci_syn, scientific = F),
          "_lambda_",input$lambda_syn,
          "_bp_per_wave_", "2e+06",
          "_sine_scaling_factor_", "1.3",
          "_NBnc_", 100, #input$NegBinomial_num_cluster,
          "_NBdis_", 1, #input$NegBinomial_dispersion,
          "_NBp_", 0.3, #input$NegBinomial_p,
          "_tecnoise_", input$tecnoise_syn,
          "_trimmed_", input$trim_value,
          ".RDS", sep="")
  )
})
results_sim <- do.call(bind_rows, results_sim)


results_to_analyze <-
  results_sim %>%
  group_by(method, 
           sim.number,
           weighting,
           threshold_log2_expgamma,
           post_prob,
           adjustment,
           alpha_value
  ) %>%
  mutate(tuning = sum(
    c(
      threshold_log2_expgamma,
      post_prob,
      weighting
    ), na.rm=T) 
  )

get_optimal_tuning <- results_to_analyze %>% 
  group_by(method, tuning, adjustment) %>% 
  summarise(mean_Metric = mean(!!sym(input$metric), na.rm = T)) %>% 
  group_by(method) %>% 
  filter(mean_Metric==max(mean_Metric, na.rm = T)) %>% 
  distinct(method, mean_Metric, .keep_all = TRUE)

x_min <- 0
x_max <- 1

best_results <- lapply(1:nrow(get_optimal_tuning), function(i){
  results_to_analyze %>% 
    filter(method == get_optimal_tuning$method[i] & 
             tuning == get_optimal_tuning$tuning[i] &
             (adjustment == get_optimal_tuning$adjustment[i] | 
                is.na(adjustment)))
})

best_results <- do.call(rbind, best_results)

best_results$adjustment_short <- NULL

best_results <- best_results %>% mutate(adjustment_short = case_when(
  adjustment == "Benjamini-Hochberg" ~ "Benj.-Hoch.",
  adjustment == "Bonferroni-Holm" ~ "Bonf.-Holm",
  adjustment == "Bonferroni" ~ "Bonf."
))

best_results <- best_results %>%
  arrange(desc(method))

best_results <- 
  melt(data = best_results, 
       id.vars = "method",
       measure.vars =  c("MCC", "Precision", "TPR"), 
       variable.name = "Metric")

unweigthed_results_no_melt <- results_to_analyze %>% 
  group_by(method) %>% 
  filter((tuning == 1 & method == "Binomial" & adjustment == "Benjamini-Hochberg") |
           (tuning == 1 & method == "ConNIS" & adjustment == "Bonferroni-Holm") |
           (threshold_log2_expgamma == 2 & method == "Exp. vs. Gamma" ) |
           (tuning == 1 & method == "Geometric" & adjustment == "Bonferroni-Holm") |
           (post_prob == 0.1 & method == "InsDens") |
           (tuning == 1 & method == "Tn5Gaps" & adjustment == "Benjamini-Hochberg")
  )

unweigthed_results <- unweigthed_results_no_melt %>%
  arrange(desc(method))

unweigthed_results <- 
  melt(data = unweigthed_results, 
       id.vars = "method",
       measure.vars =  c("MCC", "Precision", "TPR"), 
       variable.name = "Metric")

df_scales <- data.frame(
  Panel = c("MCC", "Precision", "TPR"),
  xmin = c(0, 0, 0),
  xmax = c(1, 1, 1),
  n = c(4, 4, 4)
)
df_scales <- split(df_scales, df_scales$Panel)

scales <- lapply(df_scales, function(x) {
  scale_x_continuous(limits = c(x$xmin, x$xmax), n.breaks = x$n)
})

p3_rl_max_mcc_precision_tpr <- ggplot(
  best_results,
  aes(y=fct_inorder(method), x= value, fill=method)) + 
  # geom_vline(xintercept = 0.5, color="black",linetype="dotted", alpha=0.5)+
  geom_vline(xintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
  geom_density_ridges2(alpha=0.5, quantile_lines=TRUE, 
                       quantile_fun=function(MCC,...){mean(MCC)}, 
                       rel_min_height = 0.005 ) +
  geom_density_ridges2(alpha=0, quantile_lines=TRUE,linewidth = 0.5, linetype = "dashed",
                       quantile_fun=function(MCC,...){mean(MCC)},
                       rel_min_height = 0.005, data=unweigthed_results, aes(y=fct_inorder(method), x= value, color=method)) +
  facet_grid(~Metric , scales = "free", labeller = as_labeller(c("MCC" = "MCC", "Precision" = "Corresponding Precision", "TPR" = "Corresponding TPR"))) +
  ggh4x::facetted_pos_scales(
    x = scales
  )+
  theme(legend.position="none") +
  theme(axis.text.x=element_text(size=rel(1)),
        axis.text.y=element_text(size=rel(1)),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  scale_color_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))


best_results_mcc <- best_results %>% filter(Metric == "MCC")
unweigthed_results_mcc <- unweigthed_results %>% filter(Metric == "MCC")
p3_rl_max_mcc <- ggplot(
  best_results_mcc,
  aes(y=fct_inorder(method), x= value, fill=method)) + 
  # geom_vline(xintercept = 0.5, color="black",linetype="dotted", alpha=0.5)+
  geom_vline(xintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  # geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
  geom_density_ridges2(alpha=0.5, quantile_lines=TRUE, 
                       quantile_fun=function(MCC,...){mean(MCC)}, 
                       rel_min_height = 0.005 ) +
  geom_density_ridges2(alpha=0, quantile_lines=TRUE,linewidth = 0.5, linetype = "dashed",
                       quantile_fun=function(MCC,...){mean(MCC)},
                       rel_min_height = 0.005, data=unweigthed_results_mcc, aes(y=fct_inorder(method), x= value, color=method)) +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6), angle = 0, vjust = 0, hjust=1),
        axis.title.y=element_blank())+
  scale_fill_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  scale_color_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  xlab("MCC")+
  xlim(0,1)+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS3: 50,000 IS, min. 85% insertion-free and 2% noise IS")



plot_data <- results_to_analyze %>%
  filter(
    (method == "Binomial" & adjustment == "Benjamini-Hochberg") |
      (method == "ConNIS" & adjustment == "Bonferroni-Holm") |
      (method == "Exp. vs. Gamma" & is.na(adjustment)) |
      (method == "Geometric" & adjustment == "Bonferroni-Holm") |
      (method == "InsDens" & is.na(adjustment)) |
      (method == "Tn5Gaps" & adjustment == "Benjamini-Hochberg")
  )

y_min <- 0
y_max <- 1


# data_unweigted_range <- 
#   bind_cols(unweigthed_results_no_melt %>% 
#               group_by(method) %>% 
#               mutate(N_min=TP+FP) %>% 
#               slice(which.min(N_min)) %>%  
#               select(MCC_min = MCC, N_min), 
#             unweigthed_results %>% 
#               group_by(method) %>% 
#               mutate(N_max=TP+FP) %>% 
#               slice(which.max(N_max)) %>%  
#               select(MCC_max=MCC, N_max) ) %>% 
#   select(method=method...1, MCC_min, MCC_max, N_min, N_max)



p3_setsize <- ggplot(plot_data %>%mutate(N=TP+FP), 
                     aes(x=N, y=!!sym(input$metric), color=method)) +
  # geom_rect(data = unweigthed_results_no_melt,inherit.aes = FALSE,
  #           aes(fill=method, xmin = N_min, ymin = 0, xmax = N_max, ymax = 1), alpha=0.3)+
  # scale_fill_manual(values =
  #                     c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_hline(yintercept = 0, color="black",linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = mean(plot_data$TP+plot_data$FN), color="orange",linetype="dashed", alpha=0.7)+
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE, size = 0.75) +
  xlab("Number of genes labeled 'essential'")+
  coord_cartesian(ylim=c(y_min,y_max))+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS3: 50,000 IS and min. 85% insertion-free") +
  geom_point(data = unweigthed_results_no_melt %>%
               group_by(method) %>% reframe(MCC=mean(MCC), N=mean(TP+FP)),
             aes(x=N, y=MCC), colour = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"), size=2)
# geom_segment(data = unweigthed_results_no_melt, 
#              aes(x = N_min, y = -0.05, xend = N_min, yend = 0.05), size=0.8
#              )+
# geom_segment(data = unweigthed_results_no_melt, 
#              aes(x = N_max, y = -0.05, xend = N_max, yend = 0.05), size=0.8
# )



p3_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE, size = 0.75) +
  xlab("Recall")+
  coord_cartesian(ylim=c(0,1), xlim=c(0,1))+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS3: 50,000 IS and min. 85% insertion-free")+
  geom_point(data = unweigthed_results_no_melt %>%
               group_by(method) %>% reframe(Precision=mean(Precision), TPR=mean(TPR)),
             aes(x=TPR, y=Precision), 
             colour = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"), size=2)


## Major Synthetic Setting 1 D)
input <- list()
input$unique_loci_syn <- 200000
input$lambda_syn <- 0.8
input$num_hot_spots <- 0
input$num_cold_spots <- 25
input$tecnoise_syn <- 0
input$trim_value <- 0.1
input$metric <- "MCC"

results_sim <- lapply(c("Binomial", "ConNIS", "ExpGamma", "Geometric", "InsDens", "Tn5Gaps"), function(m){
  readRDS(
    paste("./performance/Performance_",
          m,
          "_ess_ORF_", "uniform",#input$ess_ORF,
          "_numLoci_", format(input$unique_loci, scientific = F),
          "_lambda_",input$lambda_syn,
          "_num_hot_spots_", input$num_hot_spots,
          "_hot_spot_size_", 10000, #input$hot_spot_size,
          "_num_cold_spots_", input$num_cold_spots,
          "_cold_spot_size_", 10000, #input$cold_spot_size,
          "_NBnc_", 100, #input$NegBinomial_num_cluster,
          "_NBdis_", 1, #input$NegBinomial_dispersion,
          "_NBp_", 0.3, #input$NegBinomial_p,
          "_tecnoise_", input$tecnoise_syn,
          "_trimmed_", input$trim_value,
          ".RDS", sep="")
  )
})
results_sim <- do.call(bind_rows, results_sim)


results_to_analyze <-
  results_sim %>%
  group_by(method, 
           sim.number,
           weighting,
           threshold_log2_expgamma,
           post_prob,
           adjustment,
           alpha_value
  ) %>%
  mutate(tuning = sum(
    c(
      threshold_log2_expgamma,
      post_prob,
      weighting
    ), na.rm=T) 
  )

get_optimal_tuning <- results_to_analyze %>% 
  group_by(method, tuning, adjustment) %>% 
  summarise(mean_Metric = mean(!!sym(input$metric), na.rm = T)) %>% 
  group_by(method) %>% 
  filter(mean_Metric==max(mean_Metric, na.rm = T)) %>% 
  distinct(method, mean_Metric, .keep_all = TRUE)

x_min <- 0
x_max <- 1

best_results <- lapply(1:nrow(get_optimal_tuning), function(i){
  results_to_analyze %>% 
    filter(method == get_optimal_tuning$method[i] & 
             tuning == get_optimal_tuning$tuning[i] &
             (adjustment == get_optimal_tuning$adjustment[i] | 
                is.na(adjustment)))
})

best_results <- do.call(rbind, best_results)

best_results$adjustment_short <- NULL

best_results <- best_results %>% mutate(adjustment_short = case_when(
  adjustment == "Benjamini-Hochberg" ~ "Benj.-Hoch.",
  adjustment == "Bonferroni-Holm" ~ "Bonf.-Holm",
  adjustment == "Bonferroni" ~ "Bonf."
))

best_results <- best_results %>%
  arrange(desc(method))

best_results <- 
  melt(data = best_results, 
       id.vars = "method",
       measure.vars =  c("MCC", "Precision", "TPR"), 
       variable.name = "Metric")

unweigthed_results_no_melt <- results_to_analyze %>% 
  group_by(method) %>% 
  filter((tuning == 1 & method == "Binomial" & adjustment == "Benjamini-Hochberg") |
           (tuning == 1 & method == "ConNIS" & adjustment == "Bonferroni-Holm") |
           (threshold_log2_expgamma == 2 & method == "Exp. vs. Gamma" ) |
           (tuning == 1 & method == "Geometric" & adjustment == "Bonferroni-Holm") |
           (post_prob == 0.1 & method == "InsDens") |
           (tuning == 1 & method == "Tn5Gaps" & adjustment == "Benjamini-Hochberg")
  )

unweigthed_results <- unweigthed_results_no_melt %>%
  arrange(desc(method))

unweigthed_results <- 
  melt(data = unweigthed_results, 
       id.vars = "method",
       measure.vars =  c("MCC", "Precision", "TPR"), 
       variable.name = "Metric")

df_scales <- data.frame(
  Panel = c("MCC", "Precision", "TPR"),
  xmin = c(0, 0, 0),
  xmax = c(1, 1, 1),
  n = c(4, 4, 4)
)
df_scales <- split(df_scales, df_scales$Panel)

scales <- lapply(df_scales, function(x) {
  scale_x_continuous(limits = c(x$xmin, x$xmax), n.breaks = x$n)
})

p4_rl_max_mcc_precision_tpr <- ggplot(
  best_results,
  aes(y=fct_inorder(method), x= value, fill=method)) + 
  # geom_vline(xintercept = 0.5, color="black",linetype="dotted", alpha=0.5)+
  geom_vline(xintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
  geom_density_ridges2(alpha=0.5, quantile_lines=TRUE, 
                       quantile_fun=function(MCC,...){mean(MCC)}, 
                       rel_min_height = 0.005 ) +
  geom_density_ridges2(alpha=0, quantile_lines=TRUE,linewidth = 0.5, linetype = "dashed",
                       quantile_fun=function(MCC,...){mean(MCC)},
                       rel_min_height = 0.005, data=unweigthed_results, aes(y=fct_inorder(method), x= value, color=method)) +
  facet_grid(~Metric , scales = "free", labeller = as_labeller(c("MCC" = "MCC", "Precision" = "Corresponding Precision", "TPR" = "Corresponding TPR"))) +
  ggh4x::facetted_pos_scales(
    x = scales
  )+
  theme(legend.position="none") +
  theme(axis.text.x=element_text(size=rel(1)),
        axis.text.y=element_text(size=rel(1)),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_fill_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  scale_color_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))


best_results_mcc <- best_results %>% filter(Metric == "MCC")
unweigthed_results_mcc <- unweigthed_results %>% filter(Metric == "MCC")
p4_rl_max_mcc <- ggplot(
  best_results_mcc,
  aes(y=fct_inorder(method), x= value, fill=method)) + 
  # geom_vline(xintercept = 0.5, color="black",linetype="dotted", alpha=0.5)+
  geom_vline(xintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  # geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
  geom_density_ridges2(alpha=0.5, quantile_lines=TRUE, 
                       quantile_fun=function(MCC,...){mean(MCC)}, 
                       rel_min_height = 0.005 ) +
  geom_density_ridges2(alpha=0, quantile_lines=TRUE,linewidth = 0.5, linetype = "dashed",
                       quantile_fun=function(MCC,...){mean(MCC)},
                       rel_min_height = 0.005, data=unweigthed_results_mcc, aes(y=fct_inorder(method), x= value, color=method)) +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6), angle = 0, vjust = 0, hjust=1),
        axis.title.y=element_blank())+
  scale_fill_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  scale_color_manual(
    values = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"))+
  xlab("MCC")+
  xlim(0,1)+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS4: 200,000 IS, min. 80% insertion-free and 25 cold spots")



plot_data <- results_to_analyze %>%
  filter(
    (method == "Binomial" & adjustment == "Benjamini-Hochberg") |
      (method == "ConNIS" & adjustment == "Bonferroni-Holm") |
      (method == "Exp. vs. Gamma" & is.na(adjustment)) |
      (method == "Geometric" & adjustment == "Bonferroni-Holm") |
      (method == "InsDens" & is.na(adjustment)) |
      (method == "Tn5Gaps" & adjustment == "Benjamini-Hochberg")
  )

y_min <- 0
y_max <- 1


# data_unweigted_range <- 
#   bind_cols(unweigthed_results_no_melt %>% 
#               group_by(method) %>% 
#               mutate(N_min=TP+FP) %>% 
#               slice(which.min(N_min)) %>%  
#               select(MCC_min = MCC, N_min), 
#             unweigthed_results %>% 
#               group_by(method) %>% 
#               mutate(N_max=TP+FP) %>% 
#               slice(which.max(N_max)) %>%  
#               select(MCC_max=MCC, N_max) ) %>% 
#   select(method=method...1, MCC_min, MCC_max, N_min, N_max)



p4_setsize <- ggplot(plot_data %>%mutate(N=TP+FP), 
                     aes(x=N, y=!!sym(input$metric), color=method)) +
  # geom_rect(data = unweigthed_results_no_melt,inherit.aes = FALSE,
  #           aes(fill=method, xmin = N_min, ymin = 0, xmax = N_max, ymax = 1), alpha=0.3)+
  # scale_fill_manual(values =
  #                     c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_hline(yintercept = 0, color="black",linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = mean(plot_data$TP+plot_data$FN), color="orange",linetype="dashed", alpha=0.7)+
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE, size = 0.75) +
  xlab("Number of genes labeled 'essential'")+
  coord_cartesian(ylim=c(y_min,y_max))+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS4: 200,000 IS, min. 80% insertion-free and 25 cold spots") +
  geom_point(data = unweigthed_results_no_melt %>%
               group_by(method) %>% reframe(MCC=mean(MCC), N=mean(TP+FP)),
             aes(x=N, y=MCC), colour = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"), size=2)
# geom_segment(data = unweigthed_results_no_melt, 
#              aes(x = N_min, y = -0.05, xend = N_min, yend = 0.05), size=0.8
#              )+
# geom_segment(data = unweigthed_results_no_melt, 
#              aes(x = N_max, y = -0.05, xend = N_max, yend = 0.05), size=0.8
# )



p4_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE, size = 0.75) +
  xlab("Recall")+
  coord_cartesian(ylim=c(0,1), xlim=c(0,1))+
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title=element_text( face='italic', size=6, hjust = 0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  ggtitle("\nSS4: 200,000 IS, min. 80% insertion-free and 25 cold spots")+
  geom_point(data = unweigthed_results_no_melt %>%
               group_by(method) %>% reframe(Precision=mean(Precision), TPR=mean(TPR)),
             aes(x=TPR, y=Precision), 
             colour = c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0"), size=2)






  plots_syn_results_mcc_subsetsize <- 
    ggarrange(
      p1_setsize,
      p2_setsize,
      p3_setsize,
      p4_setsize,
      nrow = 2, ncol=2,
      labels=LETTERS[1:4],
      font.label = list(size = 12,color= "#525252"),
      align='v', 
      common.legend = T, legend = "bottom", legend.grob = get_legend(p4_tpr_vs_precision)
    )

  
  plots_syn_results_tpr_vs_precision <- 
    ggarrange(
      p1_tpr_vs_precision,
      p2_tpr_vs_precision,
      p3_tpr_vs_precision,
      p4_tpr_vs_precision,
      nrow = 2, ncol=2,
      labels=LETTERS[1:4],
      font.label = list(size = 12,color= "#525252"),
      align='v', 
      common.legend = T, legend = "bottom", legend.grob = get_legend(p4_tpr_vs_precision)
    )
  
  plots_syn_best_mcc <- 
    ggarrange(
      p1_rl_max_mcc,
      p2_rl_max_mcc,
      p3_rl_max_mcc,
      p4_rl_max_mcc,
      nrow = 2, ncol=2,
      labels=LETTERS[1:4],
      font.label = list(size = 12,color= "#525252"),
      align='v', 
      common.legend = T, legend = "bottom", legend.grob = get_legend(p4_tpr_vs_precision)
    )
  
  
  save_plot(filename = "./plots/syn_results_mcc_subsetsize.pdf",
            plot =plots_syn_results_mcc_subsetsize, dpi =600, base_height = 6, base_asp = 4/4 )
  
  save_plot(filename = "./plots/syn_results_tpr_vs_precision.pdf",
            plot =plots_syn_results_tpr_vs_precision, dpi =600, base_height = 6, base_asp = 4/4 )
  
  save_plot(filename = "./plots/syn_results_best_mcc.pdf",
            plot =plots_syn_best_mcc, dpi =600, base_height = 6, base_asp = 4/4 )
    
  
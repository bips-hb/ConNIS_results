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

unweigthed_results <- results_to_analyze %>% 
  group_by(method) %>% 
  filter((tuning == 1 & method == "Binomial" & adjustment == "Benjamini-Hochberg") |
           (tuning == 1 & method == "ConNIS" & adjustment == "Bonferroni-Holm") |
           (threshold_log2_expgamma == 2 & method == "Exp. vs. Gamma" ) |
           (tuning == 1 & method == "Geometric" & adjustment == "Bonferroni-Holm") |
           (post_prob == 0.1 & method == "InsDens") |
           (tuning == 1 & method == "Tn5Gaps" & adjustment == "Bonferroni-Holm")
  )

unweigthed_results <- unweigthed_results %>%
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
  xlim(0,1)


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

p1_setsize <- ggplot(plot_data, 
                     aes(x=TP+FP, y=!!sym(input$metric), color=method)) +
  geom_hline(yintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  geom_vline(xintercept = mean(plot_data$TP+plot_data$FN), color="orange",linetype="dashed", alpha=0.7)+
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE) +
  xlab("Number of genes labeled 'essential'")+
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim=c(y_min,y_max))


p1_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE) +
  xlab("True positve rate")+
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1))


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

unweigthed_results <- results_to_analyze %>% 
  group_by(method) %>% 
  filter((tuning == 1 & method == "Binomial" & adjustment == "Benjamini-Hochberg") |
           (tuning == 1 & method == "ConNIS" & adjustment == "Bonferroni-Holm") |
           (threshold_log2_expgamma == 2 & method == "Exp. vs. Gamma" ) |
           (tuning == 1 & method == "Geometric" & adjustment == "Bonferroni-Holm") |
           (post_prob == 0.1 & method == "InsDens") |
           (tuning == 1 & method == "Tn5Gaps" & adjustment == "Bonferroni-Holm")
  )

unweigthed_results <- unweigthed_results %>%
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
  # geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
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
  xlim(0,1)

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

p2_setsize <- ggplot(plot_data, 
                     aes(x=TP+FP, y=!!sym(input$metric), color=method)) +
  geom_hline(yintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  geom_vline(xintercept = mean(plot_data$TP+plot_data$FN), color="orange",linetype="dashed", alpha=0.7)+
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE) +
  xlab("Number of genes labeled 'essential'")+
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim=c(y_min,y_max))

p2_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE) +
  xlab("True positve rate")+
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1))


## Major Synthetic Setting 1 C)
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

unweigthed_results <- results_to_analyze %>% 
  group_by(method) %>% 
  filter((tuning == 1 & method == "Binomial" & adjustment == "Benjamini-Hochberg") |
           (tuning == 1 & method == "ConNIS" & adjustment == "Bonferroni-Holm") |
           (threshold_log2_expgamma == 2 & method == "Exp. vs. Gamma" ) |
           (tuning == 1 & method == "Geometric" & adjustment == "Bonferroni-Holm") |
           (post_prob == 0.1 & method == "InsDens") |
           (tuning == 1 & method == "Tn5Gaps" & adjustment == "Bonferroni-Holm")
  )

unweigthed_results <- unweigthed_results %>%
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
  # geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
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
  xlim(0,1)

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

p3_setsize <- ggplot(plot_data, 
                     aes(x=TP+FP, y=!!sym(input$metric), color=method)) +
  geom_hline(yintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  geom_vline(xintercept = mean(plot_data$TP+plot_data$FN), color="orange",linetype="dashed", alpha=0.7)+
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE) +
  xlab("Number of genes labeled 'essential'")+
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim=c(y_min,y_max))

p3_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE) +
  xlab("True positve rate")+
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1))


## Major Synthetic Setting 1 D)
input <- list()
input$unique_loci_syn <- 200000
input$lambda_syn <- 0.8
input$num_hot_spots <- 25
input$num_cold_spots <- 0
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

unweigthed_results <- results_to_analyze %>% 
  group_by(method) %>% 
  filter((tuning == 1 & method == "Binomial" & adjustment == "Benjamini-Hochberg") |
           (tuning == 1 & method == "ConNIS" & adjustment == "Bonferroni-Holm") |
           (threshold_log2_expgamma == 2 & method == "Exp. vs. Gamma" ) |
           (tuning == 1 & method == "Geometric" & adjustment == "Bonferroni-Holm") |
           (post_prob == 0.1 & method == "InsDens") |
           (tuning == 1 & method == "Tn5Gaps" & adjustment == "Bonferroni-Holm")
  )

unweigthed_results <- unweigthed_results %>%
  arrange(desc(method))

unweigthed_results <- 
  melt(data = unweigthed_results, 
       id.vars = "method",
       measure.vars =  c("MCC", "Precision", "TPR"), 
       variable.name = "Metric")

df_scales <- data.frame(
  Panel = c("MCC", "Precision", "TPR"),
  xmin = c(-1, -1, -1),
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
  # geom_vline(xintercept = -0.5, color="black",linetype="dotted", alpha=0.5) +
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
  xlim(0,1)


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

p4_setsize <- ggplot(plot_data, 
                     aes(x=TP+FP, y=!!sym(input$metric), color=method)) +
  geom_hline(yintercept = 0, color="black",linetype="dotted", alpha=0.5) +
  geom_vline(xintercept = mean(plot_data$TP+plot_data$FN), color="orange",linetype="dashed", alpha=0.7)+
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE) +
  xlab("Number of genes labeled 'essential'")+
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim=c(y_min,y_max))

p4_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_smooth(se=TRUE) +
  xlab("True positve rate")+
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1))


# plots_rl_max_mcc_precision_tpr <- ggarrange(p1_rl_max_mcc_precision_tpr + rremove("legend.title"),  
#                                             p2_rl_max_mcc_precision_tpr + rremove("legend.title"), 
#                                             p3_rl_max_mcc_precision_tpr + rremove("legend.title"), 
#                                             p4_rl_max_mcc_precision_tpr + rremove("legend.title"),
#                                             nrow = 4, ncol=1,
#                                             align='v', labels=c(
#                                               paste(LETTERS[1], c(""), sep=""), 
#                                               paste(LETTERS[2], c(""), sep=""),
#                                               paste(LETTERS[3], c(""), sep=""),
#                                               paste(LETTERS[4], c(""), sep="")
#                                             ),
#                                             common.legend = T, legend = "none"
# ) 
# 
# 
# 
# plots_subsetsize_tpr_precision  <- 
#   ggarrange( p1_setsize + rremove("legend.title") , 
#              p1_tpr_vs_precision + rremove("legend.title") , 
#              p2_setsize + rremove("legend.title") , 
#              p2_tpr_vs_precision + rremove("legend.title") ,
#              p3_setsize + rremove("legend.title") ,
#              p3_tpr_vs_precision + rremove("legend.title") ,
#              p4_setsize + rremove("legend.title") ,
#              p4_tpr_vs_precision + rremove("legend.title") ,
#              nrow = 4, ncol=2,
#              labels=c(LETTERS[1], "", LETTERS[2], "", LETTERS[3], "", LETTERS[4], ""),
#              align='v', 
#              common.legend = T, legend = "bottom"
#   ) 
# 
# 
# plots_syn_results  <- 
#   ggarrange(
#     p1_rl_max_mcc + rremove("legend.title") , 
#     p1_setsize + rremove("legend.title") , 
#     p1_tpr_vs_precision + rremove("legend.title") , 
#     p2_rl_max_mcc + rremove("legend.title") ,
#     p2_setsize + rremove("legend.title") ,
#     p2_tpr_vs_precision + rremove("legend.title") ,
#     p3_rl_max_mcc + rremove("legend.title") , 
#     p3_setsize + rremove("legend.title") ,
#     p3_tpr_vs_precision + rremove("legend.title") ,
#     p4_rl_max_mcc + rremove("legend.title") , 
#     p4_setsize + rremove("legend.title") ,
#     p4_tpr_vs_precision + rremove("legend.title") ,
#     nrow = 4, ncol=3,
#     labels=c(LETTERS[1], "", "", 
#              LETTERS[2], "", "", 
#              LETTERS[3], "", "",
#              LETTERS[4], "", ""
#     ),
#     #align='v', 
#     common.legend = T, legend = "bottom", legend.grob = get_legend(p4_tpr_vs_precision)
#   ) 
# 
# save_plot(filename = "./plots/syn_best_mcc.pdf",
#           plot =plots_rl_max_mcc_precision_tpr, dpi =600, base_height = 9.5, base_asp = 1)
# 
# 
# save_plot(filename = "./plots/syn_subsetsize.pdf",
#           plot =plots_subsetsize_tpr_precision, dpi =600, base_height = 9.5*2/3, base_asp = 2/3)

  
  plot_first_line <- ggarrange(p1_rl_max_mcc + rremove("legend.title") , 
                               p1_setsize + rremove("legend.title") , 
                               p1_tpr_vs_precision + rremove("legend.title") ,
                        nrow = 1, ncol=3,
                        #align='v', 
                        legend = "none"
  ) 
  
  plot_second_line <- ggarrange(p2_rl_max_mcc + rremove("legend.title") , 
                                   p2_setsize + rremove("legend.title") , 
                                   p2_tpr_vs_precision + rremove("legend.title") ,
                                   nrow = 1, ncol=3,
                                   #align='v', 
                                legend = "none"
  ) 
  
  plot_third_line <- ggarrange(p3_rl_max_mcc + rremove("legend.title") , 
                                p3_setsize + rremove("legend.title") , 
                                p3_tpr_vs_precision + rremove("legend.title") ,
                                nrow = 1, ncol=3,
                                #align='v', 
                               legend = "none"
  ) 
  
  plot_fourth_line <- ggarrange(p4_rl_max_mcc + rremove("legend.title") , 
                                p4_setsize + rremove("legend.title") , 
                                p4_tpr_vs_precision + rremove("legend.title") ,
                                nrow = 1, ncol=3,
                                #align='v', 
                                legend = "none"
  ) 
  
  plot_first_line <- annotate_figure(plot_first_line, top=text_grob("         SS1: 200,000 IS and min. 75% insertion", face = "italic", size = 10))
  plot_second_line <- annotate_figure(plot_second_line, top=text_grob("         SS1: 400,000 IS, min. 85% insertion free and 2% noise", face = "italic", size = 10))
  plot_third_line <- annotate_figure(plot_third_line, top=text_grob("         SS2: 200,000 IS, min. 80% insertion free and 25 cold spots", face = "italic", size = 10))
  plot_fourth_line <- annotate_figure(plot_fourth_line, top=text_grob("         SS2: 200,000 IS, min. 80% insertion free 25 and hot spots", face = "italic", size = 10))
  
  plots_syn_results  <- 
    ggarrange(
      plot_first_line,
      plot_second_line,
      plot_third_line,
      plot_fourth_line,
      nrow = 4, ncol=1,
      labels=LETTERS[1:4],
      align='v', 
      common.legend = T, legend = "bottom", legend.grob = get_legend(p4_tpr_vs_precision)
    )
  

  # save_plot(filename = "./plots/syn_results.pdf",
  #           plot =plots_syn_results, dpi =600, base_height = 9.5, base_asp = 3/4 )
  
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
  
  
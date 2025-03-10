# Real data application of Wetmore et al. (2015) with Keio library as gold standard
performances <- 
  readRDS("./performance/Performance_BW25113_realWorld.RDS")

performances_all <- 
  performances %>%
  group_by(method,  
           post_prob, 
           threshold_log2_expgamma, 
           weighting,
           adjustment,
           trimming_start,
           trimming_end) %>%
  mutate(tuning = 
           sum(c(
             threshold_log2_expgamma, 
             weighting, 
             post_prob), na.rm=T),
         N=TP+FP)

y_min <- 0
y_max <- 1

plot_data <- 
  performances_all %>% 
  filter(
    (
      (trimming_start == 0.05 & read_count_threshold == 0 & adjustment == "Benjamini-Hochberg" & method == "Binomial")| 
        (trimming_start == 0.05 & read_count_threshold == 0 & adjustment == "Bonferroni" & method == "ConNIS") | 
        (trimming_start == 0.05 & read_count_threshold == 0 & is.na(adjustment) & method == "Exp. vs. Gamma") | 
        (trimming_start == 0.05 & read_count_threshold == 0 & adjustment == "Bonferroni" & method == "Geometric") | 
        (trimming_start == 0.05 & read_count_threshold == 0 & adjustment == "Benjamini-Hochberg" & method == "Tn5Gaps") | 
        (trimming_start == 0.05 & read_count_threshold == 0 & is.na(adjustment) & method == "InsDens") 
    )
  )

plot_data_no_weight <- 
  plot_data %>% 
  filter(
    (method == "Binomial" & weighting == 1) |
      (method == "ConNIS" & weighting == 1) |
      (method == "Exp. vs. Gamma" & threshold_log2_expgamma == 2) |
      (method == "Geometric" & weighting == 1) |
      (method == "InsDens" & round(post_prob,2) == round(0.1,2)) |
      (method == "Tn5Gaps" & weighting == 1) 
  )

p1_mcc_subsetsize <- ggplot(plot_data, 
                            aes(x=N, y=MCC, color=method)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha=0.5) +
  geom_vline(xintercept = 254,  linetype="dashed", alpha = 0.7, color="orange") +
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_line(size=0.75) +
  geom_point(data=plot_data_no_weight, aes(x=N, y=MCC, color=method), size=2) +
  xlab("Number of genes labeled as 'essential'") +
  theme_minimal() +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim = c(y_min, y_max))

p1_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha=0.5) +
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_line(size=0.75) +
  geom_point(data=plot_data_no_weight, aes(x=TPR, y=Precision, color=method), size=2) +
  xlab("Recall") +
  theme_minimal() +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim = c(y_min, y_max), xlim = c(0, 1))


# Ecoli Realdata Ma et al. data
performances <- 
  readRDS("./performance/Performance_MG1655_realWorld.RDS")

performances_all <- 
  performances %>%
  group_by(method,  
           post_prob, 
           threshold_log2_expgamma, 
           weighting,
           adjustment,
           trimming_start,
           trimming_end) %>%
  mutate(tuning = 
           sum(c(
             threshold_log2_expgamma, 
             weighting, 
             post_prob), na.rm=T),
         N=TP+FP)

y_min <- 0
y_max <- 1


plot_data <- 
  performances_all %>% 
  filter(
    (
      (trimming_start == 0.05 & read_count_threshold == 0 & adjustment == "Benjamini-Hochberg" & method == "Binomial")| 
        (trimming_start == 0.05 & read_count_threshold == 0 & adjustment == "Bonferroni" & method == "ConNIS") | 
        (trimming_start == 0.05 & read_count_threshold == 0 & is.na(adjustment) & method == "Exp. vs. Gamma") | 
        (trimming_start == 0.05 & read_count_threshold == 0 & adjustment == "Bonferroni" & method == "Geometric") | 
        (trimming_start == 0.05 & read_count_threshold == 0 & adjustment == "Benjamini-Hochberg" & method == "Tn5Gaps") | 
        (trimming_start == 0.05 & read_count_threshold == 0 & is.na(adjustment) & method == "InsDens") 
    )
  )

plot_data_no_weight <- 
  plot_data %>% 
  filter(
    (method == "Binomial" & weighting == 1) |
      (method == "ConNIS" & weighting == 1) |
      (method == "Exp. vs. Gamma" & threshold_log2_expgamma == 2) |
      (method == "Geometric" & weighting == 1) |
      (method == "InsDens" & round(post_prob,2) == round(0.1,2)) |
      (method == "Tn5Gaps" & weighting == 1) 
  )

p2_mcc_subsetsize <- ggplot(plot_data, 
                            aes(x=N, y=MCC, color=method)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha=0.5) +
  geom_vline(xintercept = 300,  linetype="dashed", alpha = 0.7, color="orange") +
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_line(size=0.75) +
  geom_point(data=plot_data_no_weight, aes(x=N, y=MCC, color=method), size=2) +
  xlab("Number of genes labeled as 'essential'") +
  theme_minimal() +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim = c(y_min, y_max))

p2_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha=0.5) +
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_line(size=0.75) +
  geom_point(data=plot_data_no_weight, aes(x=TPR, y=Precision, color=method), size=2) +
  xlab("Recall") +
  theme_minimal() +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim = c(y_min, y_max), xlim = c(0, 1))

# Salmonella enterica Serovar Typhimurium real world by Mandal
performances <- 
  readRDS("./performance/Performance_14028S_realWorld.RDS")

performances_all <- 
  performances %>%
  group_by(method,  
           post_prob, 
           threshold_log2_expgamma, 
           weighting,
           adjustment,
           trimming_start,
           trimming_end) %>%
  mutate(tuning = 
           sum(c(
             threshold_log2_expgamma, 
             weighting, 
             post_prob), na.rm=T),
         N=TP+FP)

y_min <- 0
y_max <- 1


plot_data <- 
  performances_all %>% 
  filter(
    (
      (trimming_start == 0.05 & read_count_threshold == 1 & adjustment == "Benjamini-Hochberg" & method == "Binomial")| 
        (trimming_start == 0.05 & read_count_threshold == 1 & adjustment == "Bonferroni" & method == "ConNIS") | 
        (trimming_start == 0.05 & read_count_threshold == 1 & is.na(adjustment) & method == "Exp. vs. Gamma") | 
        (trimming_start == 0.05 & read_count_threshold == 1 & adjustment == "Bonferroni" & method == "Geometric") | 
        (trimming_start == 0.05 & read_count_threshold == 1 & adjustment == "Benjamini-Hochberg" & method == "Tn5Gaps") | 
        (trimming_start == 0.05 & read_count_threshold == 1 & is.na(adjustment) & method == "InsDens") 
    )
  )

plot_data_no_weight <- 
  plot_data %>% 
  filter(
    (method == "Binomial" & weighting == 1) |
      (method == "ConNIS" & weighting == 1) |
      (method == "Exp. vs. Gamma" & threshold_log2_expgamma == 2) |
      (method == "Geometric" & weighting == 1) |
      (method == "InsDens" & round(post_prob,2) == round(0.1,2)) |
      (method == "Tn5Gaps" & weighting == 1) 
  )

p3_mcc_subsetsize <- ggplot(plot_data, 
                            aes(x=N, y=MCC, color=method)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha=0.5) +
  geom_vline(xintercept = 461, linetype="dashed", alpha = 0.7, color="orange") +
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_line(size=0.75) +
  geom_point(data=plot_data_no_weight, aes(x=N, y=MCC, color=method), size=2) +
  xlab("Number of genes labeled as 'essential'") +
  theme_minimal() +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim = c(-0.5, y_max))

p3_mcc_subsetsize_window <- 
  ggplot(plot_data, 
       aes(x=N, y=MCC, color=method)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha=0.5) +
  geom_vline(xintercept = 461, linetype="dashed", alpha = 0.7, color="orange", size=1) +
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_line(size=2) +
  geom_point(data=plot_data_no_weight, aes(x=N, y=MCC, color=method), size=4) +
  xlab("Number of genes labeled as 'essential'") +
    theme_minimal() +
    theme(legend.position="none", 
          plot.title=element_blank(),
          legend.title=element_blank(), 
          legend.text=element_blank(),
          strip.text = element_text(size=rel(2)),
          axis.text.x=element_text(size = rel(2)),
          axis.text.y=element_text(size = rel(2)),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    coord_cartesian(xlim = c(200, 750), ylim=c(0.35,0.6)) +
  scale_y_continuous(breaks=c(0.4, 0.5, 0.6))

p3_tpr_vs_precision <- ggplot(plot_data, 
                              aes(x=TPR, y=Precision, color=method)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha=0.5) +
  scale_color_manual(values =
                       c("#6699CC", "#117733", "#CC6677", "#7f7f7f", "#DDCC77",  "#9651A0")) +
  geom_line(size=0.75) +
  geom_point(data=plot_data_no_weight, aes(x=TPR, y=Precision, color=method), size=2) +
  xlab("Recall") +
  theme_minimal() +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text=element_text(size=rel(0.6)),
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.6)),         
        axis.title.y = element_text(size = rel(0.6))) +
  coord_cartesian(ylim = c(y_min, y_max), xlim = c(0, 1))


# output_plot <- ggarrange(p1_mcc_subsetsize + rremove("legend.title") , 
#                          p1_tpr_vs_precision + rremove("legend.title") ,
#                          p2_mcc_subsetsize + rremove("legend.title") , 
#                          p2_tpr_vs_precision + rremove("legend.title") ,
#                          p3_mcc_subsetsize + rremove("legend.title") , 
#                          p3_tpr_vs_precision + rremove("legend.title") ,
#                          nrow = 3, ncol=2,
#                          align='', labels=c("A", "", "B", "", "C", ""),
#                          common.legend = T, legend = "bottom"
# ) 



top_plot <- ggarrange(p1_mcc_subsetsize + rremove("legend.title") , 
                      p1_tpr_vs_precision + rremove("legend.title") ,
                      nrow = 1, ncol=2,
                      align='v', legend = "none"
) 

middle_plot <- ggarrange( p2_mcc_subsetsize + rremove("legend.title") , 
                          p2_tpr_vs_precision + rremove("legend.title"), 
                      nrow = 1, ncol=2,
                      align='v', legend = "none"
) 

bottom_plot <- ggarrange(p3_mcc_subsetsize + rremove("legend.title") , 
                         p3_tpr_vs_precision + rremove("legend.title"),
                      nrow = 1, ncol=2,
                      align='v', legend = "none"
) 

top_plot <- annotate_figure(top_plot, top=text_grob("\nE. coli BW25113", face = "italic", size = 6))
middle_plot <- annotate_figure(middle_plot, top=text_grob("\nE. coli K-12 MG1655", face = "italic", size = 6))
bottom_plot <- annotate_figure(bottom_plot, top=text_grob("\nS. Typhimurium 14028S", face = "italic", size = 6))


output_plot <- ggarrange(top_plot, 
          middle_plot, 
          bottom_plot,
          nrow = 3, ncol=1,
          align='v', 
          labels=c("A", "B", "C"),
          font.label = list(size = 10,color= "#525252"),
          common.legend = T, 
          legend.grob = get_legend(p3_mcc_subsetsize), 
          legend = "bottom"
) 

save_plot(filename = "./plots/real_data.pdf",
          plot =output_plot, dpi =600, base_height = 9.5*(3/4), base_asp = (2/3))

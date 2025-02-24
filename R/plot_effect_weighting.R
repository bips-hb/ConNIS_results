# Plots for the effect of weighting

for(j in c("Binomial", "ConNIS", "Geometric", "Tn5Gaps")){
  
  selected_method <- j
  corrections <- c("Bonferroni-Holm", "Benjamini-Hochberg")
  
  plots <- lapply(seq_along(corrections), function(i){
    
    subplot_letters <- list(c("A", "B"), c("C", "D"))[[i]]
    
    selected_correction <- corrections[i]
    
    performances_14028S <- 
      readRDS("./performance/Performance_14028S_realWorld.RDS")
    
    performances_all_14028S <- 
      performances_14028S %>%
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
             N=TP+FP) %>% 
      filter(trimming_start == 0.05 & 
               read_count_threshold == 1 & 
               adjustment == selected_correction & 
               method == selected_method)
    
    performances_all_14028S$FNR <- 
      performances_all_14028S$FN/(performances_all_14028S$FN + 
                                    performances_all_14028S$TP)
    
    performances_all_14028S$dataset <- 
      "S. Typhimurium 14028S"
    
    performances_BW25113 <- 
      readRDS("./performance/Performance_BW25113_realWorld.RDS")
    
    performances_all_BW25113 <- 
      performances_BW25113 %>%
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
             N=TP+FP) %>% 
      filter(trimming_start == 0.05 & 
               read_count_threshold == 0 & 
               adjustment == selected_correction & 
               method == selected_method)
    
    performances_all_BW25113$FNR <- 
      performances_all_BW25113$FN/(performances_all_BW25113$FN + 
                                     performances_all_BW25113$TP)
    
    performances_all_BW25113$dataset <- 
      "E. coli BW25113"
    
    
    performances_MG1655 <- 
      readRDS("./performance/Performance_MG1655_realWorld.RDS")
    
    performances_all_MG1655 <- 
      performances_MG1655 %>%
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
             N=TP+FP) %>% 
      filter(trimming_start == 0.05 & 
               read_count_threshold == 0 & 
               adjustment == selected_correction & 
               method == selected_method)
    
    performances_all_MG1655$FNR <- 
      performances_all_MG1655$FN/(performances_all_MG1655$FN + 
                                    performances_all_MG1655$TP)
    
    performances_all_MG1655$dataset <- 
      "E. coli MG1655"
    
    
    FNR_plot_data <- 
      bind_rows(
        performances_all_BW25113,
        performances_all_MG1655,
        performances_all_14028S)
    
    FNR_plot_data$FDR <- 
      FNR_plot_data$FP/(FNR_plot_data$FP + FNR_plot_data$TP)
    
    FNR_vs_FPR <- ggplot(FNR_plot_data,
                         aes(x=FNR, y=FDR, color=dataset, label=weighting)) +
      geom_point(alpha=1) +
      geom_line()+
      scale_color_manual(values = c("#38CCFF", "#FFD850", "#E914B6"))+
      geom_label(size=2, nudge_y = 0.001, nudge_x = 0.03, show.legend = FALSE) +
      theme(legend.position="bottom", 
            legend.title=element_blank(),
            plot.title = element_text(hjust = 0.5, size=8)) +
      guides(fill=guide_legend(nrow=1,byrow=FALSE)) + 
      coord_cartesian(xlim=c(0,1), ylim=c(0,1))+
      ggtitle(paste(selected_correction, "correction"))
    
    FNR_vs_MCC <- ggplot(FNR_plot_data,
                         aes(x=FNR, y=MCC, color=dataset, label=weighting)) +
      geom_point(alpha=1) +
      geom_line()+
      scale_color_manual(values =
                           c("#38CCFF", "#FFD850", "#E914B6"))+
      geom_label(size=2, nudge_y = 0.025, nudge_x = -0.02, show.legend = FALSE) +
      theme(legend.position="bottom", 
            legend.title=element_blank(),
            plot.title = element_text(hjust = 0.5, size=8)) +
      guides(fill=guide_legend(nrow=1,byrow=FALSE)) + 
      coord_cartesian(xlim=c(0,1), ylim=c(0,1))+
      ggtitle(paste(selected_correction, "correction"))
    
    
    
    list(FNR_vs_FPR, FNR_vs_MCC)
    
    
  })
  
  
  output_plot <- ggarrange(plots[[1]][[1]], 
                           plots[[1]][[2]],
                           plots[[2]][[1]],
                           plots[[2]][[2]],
                           nrow = 2, ncol=2,
                           align='v', 
                           labels=LETTERS[1:4],
                           font.label = list(size = 12,color= "#525252"),
                           common.legend = T, 
                           legend = "bottom"
  ) 
  
  output_plot <- annotate_figure(output_plot, top = text_grob(selected_method,  size = 14))
  save_plot(filename = paste("./plots/effect_weighting_", selected_method, ".pdf", sep=""),
            plot =output_plot, dpi =600, base_height = 6, base_asp = (1))
  
  
}



# Analysis of perfomance of Binomial (in combination with BC approach)
source("./functions.R")


if(distortion == "sin"){
  simulated_data <- readRDS(paste("./simulatedData/SimuData_Tnseq_GeneClusters",
                                  "_ess_ORF_", ess_ORF,
                                  "_numLoci_", format(unique_loci, scientific = F),
                                  "_lambda_",lambda,
                                  "_bp_per_wave_", bp_per_wave , 
                                  "_sine_scaling_factor_", sine_scaling_factor,
                                  "_NBnc_", NegBinomial_num_cluster,
                                  "_NBdis_", NegBinomial_dispersion,
                                  "_NBp_", NegBinomial_p,
                                  "_tecnoise_", technical_noise,
                                  ".RDS", sep=""))
  
  results_binomial <- readRDS( 
    paste("./results/Results_binomial",
          "_ess_ORF_", ess_ORF,
          "_numLoci_", format(unique_loci, scientific = F),
          "_lambda_",lambda,
          "_bp_per_wave_", bp_per_wave , 
          "_sine_scaling_factor_", sine_scaling_factor,
          "_NBnc_", NegBinomial_num_cluster,
          "_NBdis_", NegBinomial_dispersion,
          "_NBp_", NegBinomial_p,
          "_tecnoise_", technical_noise,
          "_trimmed_", data_trimming,
          ".RDS", sep="")
  )
  
}else{
  simulated_data <- readRDS(paste("./simulatedData/SimuData_Tnseq_GeneClusters",
                                  "_ess_ORF_", ess_ORF,
                                  "_numLoci_", format(unique_loci, scientific = F),
                                  "_lambda_",lambda,
                                  "_num_hot_spots_", num_hot_spots,
                                  "_hot_spot_size_", hot_spot_size,
                                  "_num_cold_spots_", num_cold_spots,
                                  "_cold_spot_size_", cold_spot_size,
                                  "_NBnc_", NegBinomial_num_cluster,
                                  "_NBdis_", NegBinomial_dispersion,
                                  "_NBp_", NegBinomial_p,
                                  "_tecnoise_", technical_noise,
                                  ".RDS", sep=""))
  
  results_binomial <- readRDS( 
    paste("./results/Results_binomial",
          "_ess_ORF_", ess_ORF,
          "_numLoci_", format(unique_loci, scientific = F),
          "_lambda_",lambda,
          "_num_hot_spots_", num_hot_spots,
          "_hot_spot_size_", hot_spot_size,
          "_num_cold_spots_", num_cold_spots,
          "_cold_spot_size_", cold_spot_size,
          "_NBnc_", NegBinomial_num_cluster,
          "_NBdis_", NegBinomial_dispersion,
          "_NBp_", NegBinomial_p,
          "_tecnoise_", technical_noise,
          "_trimmed_", data_trimming,
          ".RDS", sep="")
  )
}


performance_sim <- lapply(seq_along(simulated_data), function(sim_run){
  
  alpha_value <- 0.05
  
  sim_run_data <- simulated_data[[sim_run]]$tnseqData
  
  sim_run_results <- results_binomial[[sim_run]]
  
  weights <- unique(sim_run_results$weight_insdens)
  
  performance_per_weight <- lapply(weights, function(w){
    
    truth_ess <- 
        as.numeric(sim_run_data[, "true_essentiality"] == "essential")
    
    
      sim_run_prob_binomial <- 
        sim_run_results[round(sim_run_results$weight_insdens, 2) == round(w, 2) , ]
      sim_run_prob_binomial$p_value[sim_run_prob_binomial$p_value==Inf] <- 1
      sim_run_prob_binomial$p_value[sim_run_prob_binomial$p_value==-Inf] <- 0
      
      
      pvalues <- unlist(sim_run_prob_binomial$p_value)
      
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
        weighting = w,
        post_prob = NA,
        threshold_log2_expgamma = NA,
        data_trinning = data_trimming,
        classification_performance(est.ess = est_ess_ben_hoch, true.ess = truth_ess),
        sim.number = sim_run)
      
      
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
        weighting = w,
        post_prob = NA,
        threshold_log2_expgamma = NA,
        data_trinning = data_trimming,
        classification_performance(est.ess = est_ess_bon_holm, true.ess = truth_ess),
        sim.number = sim_run)
      
      
      
      est_ess_bon <- as.numeric(pvalues <= alpha_value/length(pvalues))
      
      perf_binomial_bon <- bind_cols(
        method = paste("Binomial", sep=""), 
        alpha_value = alpha_value,
        adjustment = "Bonferroni",
        weighting = w,
        post_prob = NA,
        threshold_log2_expgamma = NA,
        data_trinning = data_trimming,
        classification_performance(est.ess = est_ess_bon, true.ess = truth_ess),
        sim.number = sim_run)
      
      out <- rbind(perf_binomial_ben_hoch, perf_binomial_bon_holm, perf_binomial_bon)

      out

    })
    performance_per_weight <- do.call(rbind, performance_per_weight)
})

performance_sim <- do.call(rbind, performance_sim)


if(distortion == "sin"){
  saveRDS(performance_sim, paste("./performance/Performance_Binomial",
                                 "_ess_ORF_", ess_ORF,
                                 "_numLoci_", format(unique_loci, scientific = F),
                                 "_lambda_",lambda,
                                 "_bp_per_wave_", bp_per_wave , 
                                 "_sine_scaling_factor_", sine_scaling_factor,
                                 "_NBnc_", NegBinomial_num_cluster,
                                 "_NBdis_", NegBinomial_dispersion,
                                 "_NBp_", NegBinomial_p,
                                 "_tecnoise_", technical_noise,
                                 "_trimmed_", data_trimming,
                                 ".RDS", sep="")
  )
  
}else{
  saveRDS(performance_sim, paste("./performance/Performance_Binomial",
                                 "_ess_ORF_", ess_ORF,
                                 "_numLoci_", format(unique_loci, scientific = F),
                                 "_lambda_",lambda,
                                 "_num_hot_spots_", num_hot_spots,
                                 "_hot_spot_size_", hot_spot_size,
                                 "_num_cold_spots_", num_cold_spots,
                                 "_cold_spot_size_", cold_spot_size,
                                 "_NBnc_", NegBinomial_num_cluster,
                                 "_NBdis_", NegBinomial_dispersion,
                                 "_NBp_", NegBinomial_p,
                                 "_tecnoise_", technical_noise,
                                 "_trimmed_", data_trimming,
                                 ".RDS", sep="")
  )
}

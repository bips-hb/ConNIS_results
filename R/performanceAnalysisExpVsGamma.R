
# Analysis of perfomance of Insdens and non-disrupted gene regions (in combination mit BC apporach)
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

  results_exp_vs_gamma <- readRDS( 
    paste("./results/Results_exp_vs_gamma",
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

  results_exp_vs_gamma <- readRDS( 
    paste("./results/Results_exp_vs_gamma",
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
  
  
  sim_run_data <- simulated_data[[sim_run]]$tnseqData
  
  sim_run_results <- results_exp_vs_gamma[[sim_run]]
  
  thresholds <- unique(sim_run_results$log2threshold)
  
  
  perfomance_log_fold_threshold <- lapply(thresholds, function(log_fold_threshold){
    
    truth_ess <- 
      as.numeric(sim_run_data[, "true_essentiality"] == "essential")
    
    sim_run_results_log2_threshold <- 
      sim_run_results %>% filter(log2threshold == log_fold_threshold)
    
    est_ess <- sim_run_results_log2_threshold$essential

    perf_expvsgamma <- bind_cols(
      method = paste("Exp. vs. Gamma", sep=""), 
      alpha_value = NA,
      adjustment = NA,
      weighting = NA,
      post_prob = NA,
      threshold_log2_expgamma = log_fold_threshold,
      data_trinning = data_trimming,
      classification_performance(est.ess = est_ess, true.ess = truth_ess),
      sim.number = sim_run)

    
    
  })
  
  perfomance_log_fold_threshold <- 
    do.call(rbind, perfomance_log_fold_threshold)
  
  perfomance_log_fold_threshold
  
})

performance_sim <- do.call(rbind, performance_sim)



if(distortion == "sin"){
  saveRDS(performance_sim, paste("./performance/Performance_ExpGamma",
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
  saveRDS(performance_sim, paste("./performance/Performance_ExpGamma",
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

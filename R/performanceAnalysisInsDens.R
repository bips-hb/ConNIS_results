# Analysis of perfomance of Insdens (in combination mit BC apporach)
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

  results_insdens <- readRDS( 
    paste("./results/Results_insdens",
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

  results_insdens <- readRDS( 
    paste("./results/Results_insdens",
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


performance_sim <- lapply(seq_along(results_insdens), function(sim_run){

  sim_run_data <- simulated_data[[sim_run]]$tnseqData
  
  sim_run_results <- results_insdens[[sim_run]][[1]]
  
  post_probs <- c(0.01, seq(0.1, 0.9, 0.1), 0.99)
  
  performance_post_probs <- lapply(post_probs, function(p){

    truth_ess <- 
      as.numeric(sim_run_data[, "true_essentiality"] == "essential")
    
    est_ess <- as.numeric(sim_run_results$post_prob >= p)
    
    perf_expvsgamma <- bind_cols(
      method = paste("InsDens", sep=""), 
      alpha_value = NA,
      adjustment = NA,
      weighting = NA,
      post_prob = p,
      threshold_log2_expgamma = NA,
      data_trinning = data_trimming,
      classification_performance(est.ess = est_ess, true.ess = truth_ess),
      sim.number = sim_run)
    

  })
  
  performance_post_probs <- do.call(rbind, performance_post_probs)

  performance_post_probs

})


performance_sim <- do.call(rbind, performance_sim)


if(distortion == "sin"){
  saveRDS(performance_sim, paste("./performance/Performance_InsDens",
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
  saveRDS(performance_sim, paste("./performance/Performance_InsDens",
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

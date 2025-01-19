# script for apllying the Method Insdens

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
}



cl <- makeCluster(nc, type="PSOCK") 

clusterExport(cl, 
  list("simulated_data",
       "data_trimming"
       ))

clusterEvalQ(cl, {
  library(tidyverse)
  library(gmp)
  library(insdens)}
  )


results_insdens <- parLapply(cl, seq_along(simulated_data), function(sim_run){
  
  # get genome length
  max_bp <- simulated_data[[sim_run]]$max_bp
  
  # get unque insertion sites
  ins_sites <- simulated_data[[sim_run]]$observed_all_unique_loci
  
  # get data
  sim_run_tnseqData <- simulated_data[[sim_run]]$tnseqData
  
  # create trimmed data:
  sim_run_tnseqData_trimmed <- sim_run_tnseqData
  # trimm by first data_trimming/2%
  sim_run_tnseqData_trimmed$right_bp <- 
    round(sim_run_tnseqData_trimmed$right_bp - sim_run_tnseqData_trimmed$gene.length * data_trimming/2)
  # trimm by last data_trimming/2%
  sim_run_tnseqData_trimmed$left_bp <- 
    round(sim_run_tnseqData_trimmed$left_bp + sim_run_tnseqData_trimmed$gene.length * data_trimming/2)
  # new length
  sim_run_tnseqData_trimmed$gene.length <- 
    sim_run_tnseqData_trimmed$right_bp - sim_run_tnseqData_trimmed$left_bp +1
  # new is count
  ins.per.gene.trimmed <- sapply(seq_along(sim_run_tnseqData_trimmed$gene), function(i){
    
    length(
      ins_sites[ins_sites>=sim_run_tnseqData_trimmed$left_bp[i] &
                  ins_sites<=sim_run_tnseqData_trimmed$right_bp[i]
      ]
    )
    
  })
  sim_run_tnseqData_trimmed$num_ins <- ins.per.gene.trimmed
  #new insertion density per gene
  sim_run_tnseqData_trimmed$ins_index <- sim_run_tnseqData_trimmed$num_ins/sim_run_tnseqData_trimmed$gene.length
  
  # start analysis
  
  #######
  # insdens 
  ins.densities.per.gene <- sim_run_tnseqData_trimmed$ins_index
  
  ins.densities.per.gene[ins.densities.per.gene == 0] <- 0.0000001
  ins.densities.per.gene[ins.densities.per.gene >= 1] <- 1

  InsDensTnseqData <- 
    tibble(gene = sim_run_tnseqData_trimmed$gene, 
           insdens = ins.densities.per.gene)
  filenname <- 
    paste("./tmpData/InsDensTnseqData_SimuRun_", 
          sim_run,
          ".csv",
          sep="")
  write.csv(x = InsDensTnseqData, 
            file= filenname, 
            row.names = F)
  
  mcmc_out <- NULL
  try(
    mcmc_out <- 
    post_prob(
      file_path = filenname,
      print_prop_sd = F,
      print_it = F)
  )
  # if insdens can not be apllied
  if(is.null(mcmc_out)){
    cat(
      unique.loci, 
      NegBinomial.num.cluster,
      NegBinomial.dispersion,
      NegBinomial.mu,
      jsim_run,
      "\n", 
      file = 'settings_error_insdens.txt', sep = ',', append = TRUE)
  }
  
  mcmc_out
  
})

stopCluster(cl)

if(distortion == "sin"){
  saveRDS(results_insdens, 
          file= paste("./results/Results_insdens",
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
  saveRDS(results_insdens, 
          file= paste("./results/Results_insdens",
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

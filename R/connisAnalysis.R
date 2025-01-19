# script for aplying our Consecutive Insertions method 

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
                   "csc",
                   "bc_prob_bigZ",
                   "cluster_determiner",
                   "classification_performance",
                   "bc_freq_bigZ",
                   "ConNIS",
                   "data_trimming"))

clusterEvalQ(cl, {
  library(tidyverse)
  library(gmp)}
)


results_connis <- parLapply(cl, seq_along(simulated_data), function(sim_run){

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
  out <- lapply(seq(0.1, 1, 0.1), function(w){
    
    # Use ConNIS Function
    ConNIS(ins.positions = ins_sites, 
             gene.names = sim_run_tnseqData_trimmed$gene, 
             gene.starts = sim_run_tnseqData_trimmed$left_bp, 
             gene.stops = sim_run_tnseqData_trimmed$right_bp, 
             num.ins.per.gene = sim_run_tnseqData_trimmed$num_ins, 
             genome.length = max_bp, 
             weighting = w)
    
  })
  
  do.call(rbind, out)
  
  
})


stopCluster(cl)

if(distortion == "sin"){
  saveRDS(results_connis, 
          file= paste("./results/Results_connis",
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
  saveRDS(results_connis, 
          file= paste("./results/Results_connis",
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

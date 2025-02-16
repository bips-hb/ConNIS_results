# Main simulation script
setwd("./R")

library(parallel)
library(tidyverse)
library(gmp)
library(insdens)

# number of unique insertion sites (IS) 
list_unique_loci <- 400000#c(50000, 100000, 200000, 400000)

# type of insertion free sections in essential genes (only uniform is available)
ess_ORF <- "uniform"  
# uniform parameter values U(lambda/1)
list_lambda <- 0.85#c(0.7, 0.75, 0.8, 0.85)

if(ess_ORF == "uniform" & any(list_lambda)> 1){
  stop("STOP: for uniform lambda needs to be smaller 1")
}

# determination which simulation setting to use
distortion <- "sin" #"spots" # "sin"  

# sin scaling values and wave length for simualtion setting 1 (SS1)
list_sine_scaling_factor <- 1.3
list_bp_per_wave <- 2000000

# number of hot&cold spots an sizes for simulation setting 2 (SS2)
if(distortion == "spots"){
  list_num_hot_spots <- c(25)
  list_hot_spot_size <- c(10000)
  list_num_cold_spots <- c(0, 25)
  list_cold_spot_size <- c(10000)
}else{
  list_num_hot_spots <- c(0)
  list_hot_spot_size <- c(10000)
  list_num_cold_spots <- c(0)
  list_cold_spot_size <- c(10000)
  
}


# paramter values for the negative Binomial distribution to determine the essential genes
# (NegBinomial_num_cluster is the number of observed clusters of essential genes)
list_NegBinomial_num_cluster <- 100
list_NegBinomial_dispersion <- 1 
list_NegBinomial_p <- 0.3 

# unique IS as technical noise (number of unique IS times a technical noise value) 
list_technical_noise <- 0.02# c(0, 0.02)

# total proportion of trimming of genes (half of the values for the start of genes, half for the end)
list_data_trimming <- 0.1#c(0, 0.1)

# number of simualtion runs per setting
num_simu <- 100

# number of workers
nc <- 50

for(unique_loci in list_unique_loci){
  
  for(lambda in list_lambda){
    
    for(sine_scaling_factor in list_sine_scaling_factor){
      
      for(bp_per_wave in list_bp_per_wave){
        
        for(num_hot_spots in list_num_hot_spots){
          
          for(hot_spot_size in list_hot_spot_size){
            
            for(num_cold_spots in list_num_cold_spots){
              
              for(cold_spot_size in list_cold_spot_size){
                
                for(NegBinomial_num_cluster in list_NegBinomial_num_cluster){
                  
                  for(NegBinomial_dispersion in list_NegBinomial_dispersion){
                    
                    for(NegBinomial_p in list_NegBinomial_p){
                      
                      for(technical_noise in list_technical_noise){
                        
                        # Data simulation
                        
                        cat("################################\n")
                        cat(paste("Starting data simulation with\n",
                                  unique_loci, " unique loci\n",
                                  lambda, " lambda\n",
                                  distortion, " distortion\n",
                                  num_hot_spots, " num hot spots\n",
                                  hot_spot_size, " hot spot size\n",
                                  num_cold_spots, " num cold spots\n",
                                  cold_spot_size, " cold spot size\n",
                                  NegBinomial_num_cluster, " cluster (draws NB-Dis)\n",
                                  NegBinomial_dispersion, " NB-dispersion\n",
                                  NegBinomial_p, " NB-prob\n",
                                  technical_noise, " Technical noise\n",
                                  sep=""))
                        
                        source("./dataSimulation.R")
                        
                        for(data_trimming in list_data_trimming){
                          
                          # Method Analysis
                          
                          cat("################################\n")
                          cat(paste("Starting Method Analysis for\n",
                                    unique_loci, " unique loci\n",
                                    lambda, " lambda\n",
                                    distortion, " distortion\n",
                                    num_hot_spots, " num hot spots\n",
                                    hot_spot_size, " hot spot size\n",
                                    num_cold_spots, " num cold spots\n",
                                    cold_spot_size, " cold spot size\n",
                                    NegBinomial_num_cluster, " cluster (draws NB-Dis)\n",
                                    NegBinomial_dispersion, " NB-dispersion\n",
                                    NegBinomial_p, " NB-prob\n",
                                    technical_noise, " Technical noise\n",
                                    data_trimming, " trimmed\n",
                                    sep=""))
                          
                          cat(paste("Binomial\n"))
                          source("./binomialAnalysis.R")
                          cat(paste("ConNIS\n"))
                          source("./connisAnalysis.R")
                          cat(paste("ExpVsGamma\n"))
                          source("./expVsGammaAnalysis.R")
                          cat(paste("Geometric\n"))
                          source("./geometricAnalysis.R")
                          cat(paste("InsDens\n"))
                          source("./insdensAnalysis.R")
                          cat(paste("Tn5gap\n"))
                          source("./tn5gapAnalysis.R")
                          
                          # Performance Analysis
                          
                          cat("################################\n")
                          cat(paste("Starting performance analysis for\n",
                                    unique_loci, " unique loci\n",
                                    lambda, " lambda\n",
                                    distortion, " distortion\n",
                                    num_hot_spots, " num hot spots\n",
                                    hot_spot_size, " hot spot size\n",
                                    num_cold_spots, " num cold spots\n",
                                    cold_spot_size, " cold spot size\n",
                                    NegBinomial_num_cluster, " cluster (draws NB-Dis)\n",
                                    NegBinomial_dispersion, " NB-dispersion\n",
                                    NegBinomial_p, " NB-prob\n",
                                    technical_noise, " Technical noise\n",
                                    data_trimming, " trimmed\n",
                                    sep=""))
                          
                          cat(paste("Binomial\n"))
                          source("./performanceAnalysisBinomial.R")
                          cat(paste("ConNIS\n"))
                          source("./performanceAnalysisConNIS.R")
                          cat(paste("ExpVsGamma\n"))
                          source("./performanceAnalysisExpVsGamma.R")
                          cat(paste("Geometric\n"))
                          source("./performanceAnalysisGeometric.R")
                          cat(paste("InsDens\n"))
                          source("./performanceAnalysisInsDens.R")
                          cat(paste("Tn5gap\n"))
                          source("./performanceAnalysisTn5Gaps.R")
                          
                        }
                        
                      }
                      
                    }    
                    
                  }
                  
                }
                
              }
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
}




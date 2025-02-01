# This script simualtes data based on a specific user selected data-structure like a real world data
# use Keio E.coli for simualtion of data
# https://www.embopress.org/doi/full/10.1038/msb4100050

# Packages
library(tidyverse)
library(parallel)

source("./functions.R")

#------------------------------------------------------------------------------
# Data Preparation
#------------------------------------------------------------------------------

# load data 
ecoli <- readxl::read_xls("./data_for_synthetic_data_generation/msb4100050-sup-0004.xls", skip = 2)

# we use the W3110 wildtype mutant (we only need the first 2 and 7:10 columns) 
ecoli <- ecoli[,c(1:3, 5:6)]

# change names
names(ecoli) <- c(
  "ECK_number",
  "gene",
  "JW_id",
  "left_bp",
  "right_bp"
)

# there are unknown genes in the data, coded as "none"; we have to make them 
# unique:
ecoli$ECK_number[ecoli$ECK_number == "none"] <- 
  paste("unknown_", 1:sum(ecoli$gene == "none"), sep="")

ecoli$gene[ecoli$gene == "none"] <- 
  paste("unknown_", 1:sum(ecoli$gene == "none"), sep="")

# get gene names appearing multiple times
mult.genes <- names(which(table(ecoli$gene)>1))

# remove multiple entry 
for(i in mult.genes){
  tmp_gene <- ecoli %>% filter(gene %in% i)
  one_gene <- tmp_gene[1,]
  one_gene$right_bp <- max(tmp_gene$right_bp)
  index_multiple_gene <- which(ecoli$gene %in% i)
  ecoli <- ecoli[-c(index_multiple_gene[-1]),]
  ecoli[ecoli$gene == i, ] <- one_gene
  rm(list=c("tmp_gene", "one_gene", "index_multiple_gene"))
}
rm(i)
rm(mult.genes)

# get min and max basepair of ECOLI
min_bp <- 1
max_bp <- 4646332 # see https://www.ncbi.nlm.nih.gov/nuccore/AP009048.1

# Create a new data frame based on the (in this case) Ecoli data.
# This data set will be the standard form for our simulation, that is you can
# use any other kind of data set in the beginning or construct your own (different
# genes length etc) as long as you have the following columns that follow our
# name scheme

tnseqData <- ecoli
rm(ecoli)

### start generation of no-name dataset
# new gene names:
tnseqData$gene <- paste("gene_", 1:nrow(tnseqData), sep="")

names(tnseqData) <- c("UGI1", "gene", "UGI2", "left_bp", "right_bp") 
tnseqData$UGI1 <- paste("UGI1_", 1:nrow(tnseqData), sep="") 
tnseqData$UGI2 <- paste("UGI2_", 1:nrow(tnseqData), sep="") 

tnseqData <- tnseqData[,c("gene", "UGI1", "UGI2", "left_bp", "right_bp")]



# get the gene length
tnseqData$gene.length <- tnseqData$right_bp - tnseqData$left_bp + 1 

num_genes <- nrow(tnseqData)

#------------------------------------------------------------------------------
# Parallel Processing Setup
#------------------------------------------------------------------------------

cl <- makeCluster(nc, type="PSOCK")

clusterExport(cl,
  list("tnseqData",
        "ess_ORF",
       "num_genes",
       "unique_loci",
       "lambda", 
       "distortion",
       "num_hot_spots", 
       "hot_spot_size", 
       "num_cold_spots", 
       "cold_spot_size", 
       "sine_scaling_factor",
       "bp_per_wave",
       "sine_prob_scaling",
       "NegBinomial_num_cluster", 
       "NegBinomial_dispersion", 
       "NegBinomial_p", 
       "max_bp",
       "technical_noise",
  "num_simu"))

clusterEvalQ(cl, {
library(tidyverse)
library(gmp)
library(insdens)}
)

simulated_data <- parLapply(cl, 1:num_simu, function(sim_run){
  # set seed for drawing cluster sizes (NegBinomial_num_cluster samples) from
  # a negative binomial distribution
  set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
             num_hot_spots * 1000 + num_cold_spots * 1000000)
  cluster.sizes.sim <- 
    rnbinom(NegBinomial_num_cluster, 
            NegBinomial_dispersion, 
            NegBinomial_p)+1 # +1 to avoid clusters of sizes 0
  
  # number of (non)-essential genes
  num_ess_genes <- sum(cluster.sizes.sim)
  num.noness.genes <- num_genes - num_ess_genes
  
  # We mix sections with non-essentials and essentials
  # To make the generation more straight forward, we determine first where the
  # essentials should be inserted by generating a sequence of non-essentials
  # and essentials where the latter will be filled with multiple essentials
  # set.seed(sim_run)
  set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
             num_hot_spots * 1000 + num_cold_spots * 1000000)
  tmp.gene.true_essentiality <- 
    sample(
      c(
        rep("non-essential", num_genes -sum(cluster.sizes.sim)),
        rep("essentials", length(cluster.sizes.sim))
      )
    )
  
  # Next, we generate a vector "gene.true_essentiality" telling if 
  # genes are essential or non-essential
  gene.true_essentiality <- NULL
  cluster.counter <- 1 
  for(i in tmp.gene.true_essentiality){
    
    if(i == "essentials"){
      
      gene.true_essentiality <- 
        c(gene.true_essentiality, 
          rep("essential", cluster.sizes.sim[cluster.counter])
        )
      
      cluster.counter <- cluster.counter+1
      
    }else{
      
      gene.true_essentiality <- c(gene.true_essentiality, "non-essential")
      
    }
  }
  rm(i)
  
  # we label each gene with its "true_essentiality"
  tnseqData$true_essentiality <- gene.true_essentiality
  
  
  # each potential insertion site of the genome gets the same "starting" probability, 
  # i.e. 1 divided by the length of the genome
  all_my_probs <- rep(1/max_bp, max_bp)
  
  # next, we determine if the distortion along the genome is sin or with hot/coldspots
  if(distortion == "sin"){
    scaling_factors <- 
      sine_prob_scaling( bp.per.wave = bp_per_wave , sine.scaling.factor = sine_scaling_factor )
    
    num_cold_spots <- 0
    cold_spot_size <- 0
    num_hot_spots <- 0
    hot_spot_size <- 0
    
    all_my_probs <- all_my_probs * scaling_factors
    
  }else{
    # only if at least some cold/hotspots are on the genome. Otherwise no distortion
    if(num_cold_spots + num_hot_spots > 0){
      
      genome <- 1:(max_bp-max(cold_spot_size, hot_spot_size))
      
      set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
                 num_hot_spots * 1000 + num_cold_spots * 1000000)
      
      # sample a sequence of hot- and coldspots
      seq_spots <- 
        sample(
          c(
            rep("H", num_hot_spots),
            rep("C", num_cold_spots)
          ), 
          size = num_cold_spots + num_hot_spots, 
          replace = F)
      
      spots_data <- tibble()  
      
      delete_genome <- NULL
      work_genome <- genome  
      counter <- 0
      
      # for each spot we sample a starting point
      for(spot in seq_spots){
        
        counter <- counter + 1
        set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
                   num_hot_spots * 1000 + num_cold_spots * 1000000 + 11111 * counter)
        
        spot_start <- sample(work_genome, 1)
        
        spot_end <- spot_start + (max(cold_spot_size, hot_spot_size)-1)
        
        spots_data <- rbind(
          spots_data,
          tibble(spot_kind = spot,
                 spot_start = spot_start,
                 spot_end = spot_end))
        
        remove_bp <- 
          max(1,(spot_start - (max(cold_spot_size, hot_spot_size)))) : 
          (spot_end+(max(cold_spot_size, hot_spot_size)))
        
        delete_genome <- c(delete_genome, remove_bp)
        
        work_genome <- genome[-delete_genome]
        
      }
      
      rm(spot)
      
      for(spot in 1:nrow(spots_data)){
        
        if(spots_data$spot_kind[spot] == "C"){
          all_my_probs[spots_data$spot_start[spot] : spots_data$spot_end[spot]] <- 
            all_my_probs[spots_data$spot_start[spot] : spots_data$spot_end[spot]] * 0.1
        }else{
          all_my_probs[spots_data$spot_start[spot] : spots_data$spot_end[spot]] <- 
            all_my_probs[spots_data$spot_start[spot] : spots_data$spot_end[spot]] * 5
        }
        
      }
      rm(spot)
      
      
    }
    
  }
  
  
  # We add variance to each potential insertion site by adding values from a
  # normal distribution
  for(i in 1:nrow(tnseqData)){
    set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
               num_hot_spots * 1000 + num_cold_spots * 1000000 + 1000*i)
    if(i == 1){
      all_my_probs[1: (tnseqData$left_bp[i])-1 ] <-  
        all_my_probs[1: tnseqData$left_bp[i]-1 ] + 
        rnorm(1, 0, 1/(max_bp*3))
      all_my_probs[tnseqData$left_bp[i] : (tnseqData$right_bp[i])] <-  
        all_my_probs[tnseqData$left_bp[i] : (tnseqData$right_bp[i]) ] + 
        rnorm(1, 0, 1/(max_bp*3))
    }else if(i == nrow(tnseqData)){
      all_my_probs[tnseqData$left_bp[i] : (tnseqData$right_bp[i])] <-  
        all_my_probs[tnseqData$left_bp[i] : (tnseqData$right_bp[i]) ] + 
        rnorm(1, 0, 1/(max_bp*3))
      all_my_probs[(tnseqData$right_bp[i-1] +1) : (tnseqData$left_bp[i]-1)] <-  
        all_my_probs[(tnseqData$right_bp[i-1] +1 ) : (tnseqData$left_bp[i]-1)] + 
        rnorm(1, 0, 1/(max_bp*3))
      all_my_probs[(tnseqData$right_bp[i] +1) : max_bp] <-  
        all_my_probs[(tnseqData$right_bp[i] +1) : max_bp] + 
        rnorm(1, 0, 1/(max_bp*3))
    }else{
      all_my_probs[tnseqData$left_bp[i] : (tnseqData$right_bp[i])] <-  
        all_my_probs[tnseqData$left_bp[i] : (tnseqData$right_bp[i]) ] + 
        rnorm(1, 0, 1/(max_bp*3))
      all_my_probs[(tnseqData$right_bp[i-1] +1) : (tnseqData$left_bp[i]-1)] <-  
        all_my_probs[(tnseqData$right_bp[i-1] +1 ) : (tnseqData$left_bp[i]-1)] + 
        rnorm(1, 0, 1/(max_bp*3))
    }
    
  }
  #
  
  # Let only a fraction of size lambda of essential genes have insertion sites
  # with probability 0 (i.e. being insertion free)
  for(i in 1:nrow(tnseqData)){
    
    if(tnseqData[i,"true_essentiality"] == "essential"){
      
      set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
                 num_hot_spots * 1000 + num_cold_spots * 1000000 + i*100)
      lambda_per_gene <- runif(1, lambda, 1)
      
      if(lambda_per_gene == 1){
        all_my_probs[as.numeric(tnseqData[i,"left_bp"]) : as.numeric(tnseqData[i,"right_bp"])] <- 0
      }
      if(lambda_per_gene < 1){
        set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
                   num_hot_spots * 1000 + num_cold_spots * 1000000 + i*100)
        # set random start of essential part within an essential gene
        start_bp_for_essential_gene_part <- sample(0:floor(as.numeric(tnseqData[i,"gene.length"])*(1-lambda_per_gene)),1)
        
        all_my_probs[
          as.numeric(tnseqData[i,"left_bp"]+ start_bp_for_essential_gene_part) : 
            (as.numeric(tnseqData[i,"left_bp"]) + 
               start_bp_for_essential_gene_part +
               floor(as.numeric(tnseqData[i,"gene.length"])*lambda_per_gene))] <- 0
      }
      
      
    }
    
  }
  rm(i)
  
  # all negative probabilies due to the normal variance term are set to 0
  all_my_probs[all_my_probs < 0] <- 0
  
  # normalize to get final probabilites for each potential insertion site
  all_my_probs <- all_my_probs/sum(all_my_probs)
  
  # sample observed insertion sites 
  # we have to sample 'unique_loci * (1+0.05/(400000/unique_loci))' i.e. a bit 
  # more insertion sites than given by 'unique_loci'. This is neccessary because
  # drawing without replacement is not feasible for over 4,000,000 different 
  # individual probabilities in the sample function. Thus, we use drawing with
  # repalcment but need to drasw a bit more due to replicates
  set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
             num_hot_spots * 1000 + num_cold_spots * 1000000)
  observed_all_unique_loci <- 
    sort(
      unique(sample(1:max_bp, unique_loci * (1+0.05/(400000/unique_loci)), prob = all_my_probs, replace = T))
    ) 
  
  # technical noise are insertion sites that can be occur everywhere on the genome:
  if(technical_noise > 0){
    set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
               num_hot_spots * 1000 + num_cold_spots * 1000000)
    ins_tec_noise <- sample(c(1:max_bp)[-observed_all_unique_loci],
                            technical_noise*unique_loci, replace = F)
    observed_all_unique_loci <- sort(unique(c(observed_all_unique_loci,ins_tec_noise)))
  }
  
  # count insertions per gene
  ins.per.gene <- sapply(seq_along(tnseqData$gene), function(i){
    
    length(
      observed_all_unique_loci[observed_all_unique_loci>=tnseqData$left_bp[i] &
                                 observed_all_unique_loci<=tnseqData$right_bp[i]
      ]
    )
    
  })
  
  # add to tnseqdata
  tnseqData$num_ins <- ins.per.gene
  tnseqData$ins_index <- tnseqData$num_ins/tnseqData$gene.length
  
  # Return results
  out <- list(
    tnseqData = tnseqData,
    sim_run = sim_run,
    cluster.sizes.sim = cluster.sizes.sim,
    NegBinomial_num_cluster = NegBinomial_num_cluster,
    NegBinomial_dispersion = NegBinomial_dispersion,
    NegBinomial_p = NegBinomial_p, 
    unique_loci = unique_loci,
    observed_all_unique_loci =observed_all_unique_loci,
    num_cold_spots = num_cold_spots,
    cold_spot_size = cold_spot_size,
    num_hot_spots = num_hot_spots,
    hot_spot_size = hot_spot_size,
    bp.per.wave = bp_per_wave , 
    sine.scaling.factor = sine_scaling_factor,
    max_bp = max_bp
  )
  
  out 
})

stopCluster(cl)

#------------------------------------------------------------------------------
# Save Results
#------------------------------------------------------------------------------

if(distortion == "sin"){
  saveRDS(simulated_data, 
      paste("./simulatedData/SimuData_Tnseq_GeneClusters",
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
  saveRDS(simulated_data, 
        paste("./simulatedData/SimuData_Tnseq_GeneClusters",
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


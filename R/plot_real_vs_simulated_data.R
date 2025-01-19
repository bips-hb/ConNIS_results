# Examples of IS Distributions along the genome

## Real data

### Wetmore es all data
# https://journals.asm.org/doi/10.1128/mbio.00306-15
# Escherichia coli K-12 strain BW25113

length_genome <- 4631469

all_genes <- read_csv("./bw25113_data/genes_wetmore_et_al.csv")

all_genes <- tibble(all_genes)

colnames(all_genes)[c(5:6,8)] <-
  c("gene_start", "gene_end", "gene")

all_genes <- all_genes[all_genes$gene_start <= 4631469, ]

# read insertion data from Wetmore et al. (2015) 
# Original source: https://genomics.lbl.gov/supplemental/rbarseq/html/Keio/all.poolcount

all_IS_raw <- read_tsv("./bw25113_data/all.poolcount")

# Only the first four samples of set1 are T0 data unconstrained in LB medium
all_IS <- 
  tibble(pos = sort(
    unique(
      (all_IS_raw %>% 
         filter(Keio_ML9_set1.IT001 >0 | 
                  Keio_ML9_set1.IT002  >0 | 
                  Keio_ML9_set1.IT003  > 0 | 
                  Keio_ML9_set1.IT004) 
      )$pos
    )
  )
  )


gene_data <- lapply(1:nrow(all_genes), function(i){
  
  IS_gene_i <- all_IS$pos[which(all_genes$gene_start[i] <= all_IS$pos &
                                  all_genes$gene_end[i] >= all_IS$pos)]
  
  num_ins <- length(IS_gene_i)
  
  gene <- all_genes$gene[i]
  if(is.null(gene)){
    gene <- NA
  }
  
  
  if(all_genes$gene_start[i] < min(IS_gene_i) |
     length(IS_gene_i)==0){
    start_is <- all_genes$gene_start[i]-1
  }else{
    start_is <- min(IS_gene_i)
  }
  
  if(all_genes$gene_end[i] > max(IS_gene_i) |
     length(IS_gene_i)==0){
    end_is <- all_genes$gene_end[i]+1
  }else{
    end_is <- max(IS_gene_i)
  }
  
  tibble(entry = all_genes$locusId[i],
         gene = gene,
         gene_start = all_genes$gene_start[i],
         gene_end = all_genes$gene_end[i],
         length = all_genes$gene_end[i] - all_genes$gene_start[i] +1,
         num_ins = num_ins,
         nondis = max(diff(unique(c(start_is, IS_gene_i, end_is)))-1)/
           (all_genes$gene_end[i] - all_genes$gene_start[i] + 1),
         insdens = num_ins/(all_genes$gene_end[i] - all_genes$gene_start[i] + 1)
  )
  
})


gene_data <- do.call(rbind, gene_data)

gene_data$nondis[is.na(gene_data$nondis)] <- 1
gene_data$insdens[is.na(gene_data$insdens)] <- 0

gene_data$num_ins[is.na(gene_data$num_ins)] <- 0

gene_data <- 
  gene_data[!duplicated(gene_data),]

gene_data$gene[is.na(gene_data$gene)] <- 
  paste("unnamed_gene_", seq_along(gene_data$gene[is.na(gene_data$gene)]),
        sep="")


p1 <- ggplot(gene_data,
             aes(x=round((gene_start+gene_end)/2), y=insdens))+
  geom_point(alpha=0.2, color="cyan4", size=0.5)+
  theme(legend.position = "none") +
  ylab("Gene-wise insertion density") +
  xlab("Genome postion") +
  geom_smooth(color="firebrick")+
  geom_hline(yintercept = nrow(all_IS)/length_genome, color="darkorange",linetype="dashed")+
  theme(legend.position="none",
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.title.y = element_text(size = rel(0.6))) +
  ylim(0,0.15)

rm(gene_data)

### Goodall et al. 2018
# https://journals.asm.org/doi/10.1128/mbio.02096-17

length_genome <- 4631469

all_genes_raw <- tibble(
  read.csv("./bw25113_data/bw25113_genbank.txt", sep="\n", header = F))

all_genes_raw <- all_genes_raw[grep(">", all_genes_raw$V1),]

all_genes <- lapply(1:nrow(all_genes_raw), function(i){
  tmp_str <- as.character(unlist(all_genes_raw[i,]))
  gene <- stringr::str_match(tmp_str, "\\[gene=\\s*(.*?)\\s*\\]")[,2]
  locus_tag <- stringr::str_match(tmp_str, "\\[locus_tag=\\s*(.*?)\\s*\\]")[,2]
  protein <- stringr::str_match(tmp_str, "\\[protein=\\s*(.*?)\\s*\\]")[,2]
  position <- stringr::str_match(tmp_str, "\\[location=\\s*(.*?)\\s*\\]")[,2]
  
  if(length(grep("join", position))>0){
    
    insert_comp <- FALSE
    if(length(grep("complement", position)) > 0){
      left_comp <- "complement("
      right_comp <- ")"
      insert_comp <- TRUE
    }
    
    position <- 
      stringr::str_match(tmp_str, "join\\(\\s*(.*?)\\s*\\)")[,2]
    
    position <- 
      str_replace_all(position, ",", "..")
    
    position <- paste(
      min(str_split(position, "\\.\\.")[[1]]), 
      "..",
      max(str_split(position, "\\.\\.")[[1]]),
      sep=""
    )  
    
    if(insert_comp){
      position <- 
        paste(left_comp,
              position,
              right_comp,
              sep="")
    }
    
  }
  
  
  if(length(grep("complement", position)) > 0){
    Orientation <- "minus"
    position <- str_remove(position, "complement\\(")
    position <- str_remove(position, "\\)")
  }
  
  if(length(grep("..", position)) > 0){
    gene_start <- as.numeric(str_split(position, "\\.\\.")[[1]][1])
    gene_end <- as.numeric(str_split(position, "\\.\\.")[[1]][2])
  }
  
  if(length(grep(" ", position)) > 0){
    gene_start <- as.numeric(str_split(position, " ")[[1]][1])
    gene_end <- as.numeric(str_split(position, " ")[[1]][2])
  }
  
  gene_length <- gene_end-gene_start+1
  
  tibble(
    gene,
    locus_tag,
    gene_start,
    gene_end,
    gene_length,
    protein
  )
})

all_genes <- do.call(rbind,all_genes)

# need to remove "YdfJ" because of false start stop annotation (its length is over 20000)
all_genes <- all_genes[all_genes$gene_length < 10000,] 


all_genes <- all_genes[all_genes$gene_length >= 100,] 

#IS from Goodall
all_IS_raw <- read.csv("./bw25113_data/Ecoli_BW25113_chlor-tn_position_depth_count.txt")
all_IS <- tibble(pos=which(all_IS_raw[,1] >0),
                 n=all_IS_raw[which(all_IS_raw[,1] >0),1],
                 strand = NA)


gene_data <- lapply(1:nrow(all_genes), function(i){
  
  IS_gene_i <- all_IS$pos[which(all_genes$gene_start[i] <= all_IS$pos &
                                  all_genes$gene_end[i] >= all_IS$pos)]
  
  num_ins <- length(IS_gene_i)
  
  gene <- all_genes$gene[i]
  if(is.null(gene)){
    gene <- NA
  }
  
  
  if(all_genes$gene_start[i] < min(IS_gene_i) |
     length(IS_gene_i)==0){
    start_is <- all_genes$gene_start[i]-1
  }else{
    start_is <- min(IS_gene_i)
  }
  
  if(all_genes$gene_end[i] > max(IS_gene_i) |
     length(IS_gene_i)==0){
    end_is <- all_genes$gene_end[i]+1
  }else{
    end_is <- max(IS_gene_i)
  }
  
  tibble(entry = all_genes$locusId[i],
         gene = gene,
         gene_start = all_genes$gene_start[i],
         gene_end = all_genes$gene_end[i],
         length = all_genes$gene_end[i] - all_genes$gene_start[i] +1,
         num_ins = num_ins,
         nondis = max(diff(unique(c(start_is, IS_gene_i, end_is)))-1)/
           (all_genes$gene_end[i] - all_genes$gene_start[i] + 1),
         insdens = num_ins/(all_genes$gene_end[i] - all_genes$gene_start[i] + 1)
  )
  
})


gene_data <- do.call(rbind, gene_data)

gene_data$nondis[is.na(gene_data$nondis)] <- 1
gene_data$insdens[is.na(gene_data$insdens)] <- 0

gene_data$num_ins[is.na(gene_data$num_ins)] <- 0

gene_data <- 
  gene_data[!duplicated(gene_data),]

gene_data$gene[is.na(gene_data$gene)] <- 
  paste("unnamed_gene_", seq_along(gene_data$gene[is.na(gene_data$gene)]),
        sep="")

p2 <- ggplot(gene_data,
             aes(x=round((gene_start+gene_end)/2), y=insdens))+
  geom_point(alpha=0.2, color="cyan4", size=0.5)+
  theme(legend.position = "none") +
  ylab("Gene-wise insertion density") +
  xlab("Genome postion") +
  geom_smooth(color="firebrick")+
  geom_hline(yintercept = nrow(all_IS)/length_genome, color="darkorange",linetype="dashed")+
  theme(legend.position="none",
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.title.y = element_text(size = rel(0.6))) 

rm(gene_data)

### Pechter et al. (2016)
# https://journals.asm.org/doi/10.1128/jb.00771-15
# Rhodopseudomonas palustris strain CGA009
length_genome <- 5467641

rpa_TA_sites <- readxl::read_xlsx("./rpa_CGA009/zjb999093952sd3.xlsx", skip = 2)
rpa_genes <- readxl::read_xlsx("./rpa_CGA009/zjb999093952sd4.xlsx", skip = 3)
rpa_shared_ess_genes <- readxl::read_xlsx("./rpa_CGA009/zjb999093952sd4.xlsx", sheet = 2, skip=1)

all_IS <- tibble(pos=sort(unique(rpa_TA_sites$Position)))

gene_data <- lapply(1:nrow(rpa_genes), function(i){
  
  IS_gene_i <- all_IS$pos[which(rpa_genes$Start[i] <= all_IS$pos &
                                  rpa_genes$Stop[i] >= all_IS$pos)]
  
  num_ins <- length(IS_gene_i)
  
  gene <- rpa_genes$Gene[i]
  if(is.null(gene)){
    gene <- NA
  }
  
  tibble(entry = rpa_genes$`Locus Tag`[i],
         gene = gene,
         gene_start = rpa_genes$Start[i],
         gene_end = rpa_genes$Stop[i],
         length = rpa_genes$Stop[i] - rpa_genes$Start[i] +1,
         num_ins = num_ins,
         nondis = max(diff(c(rpa_genes$Start[i], IS_gene_i, rpa_genes$Stop[i]))+1)/
           (rpa_genes$Stop[i] - rpa_genes$Start[i] + 1),
         insdens = num_ins/(rpa_genes$Stop[i] - rpa_genes$Start[i] + 1)
  )
  
})


gene_data <- do.call(rbind, gene_data)

gene_data$nondis[is.na(gene_data$nondis)] <- 1
gene_data$insdens[is.na(gene_data$insdens)] <- 0

gene_data$num_ins[is.na(gene_data$num_ins)] <- 0

gene_data$gene[is.na(gene_data$gene)] <- 
  paste("unnamed_gene_", seq_along(gene_data$gene[is.na(gene_data$gene)]),
        sep="")

p3 <- ggplot(gene_data,
             aes(x=round((gene_start+gene_end)/2), y=insdens))+
  geom_point(alpha=0.2, color="cyan4", size=0.5)+
  theme(legend.position = "none") +
  ylab("Gene-wise insertion density") +
  xlab("Genome postion") +
  geom_smooth(color="firebrick")+
  theme(legend.position="none",
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.title.y = element_text(size = rel(0.6))) +
  geom_hline(yintercept = nrow(all_IS)/length_genome, color="darkorange",linetype="dashed")

rm(gene_data)

### Bleem et al. (2023)
# https://www.sciencedirect.com/science/article/pii/S2211124723008586

length_genome <- 4199332

all_IS_raw <- readxl::read_xlsx("./syk6_data/mmc2.xlsx")

all_IS_raw <- all_IS_raw %>% dplyr::select(pos, strand)

all_IS <- all_IS_raw[!duplicated(all_IS_raw$pos),]

all_IS <- all_IS %>% arrange(pos)

all_genes_SYK6 <- read.csv("./syk6_data/All-genes-of-S.-sp.-SYK-6.txt", sep="\t")
all_genes_SYK6 <- tibble(all_genes_SYK6)
colnames(all_genes_SYK6)[c(3:4)] <- 
  c("gene_start", "gene_end")

plasmid_genes <- readRDS("./syk6_data/syk6_plasmid_gene_names.RDS")

genome_SYK6 <- all_genes_SYK6[!c(all_genes_SYK6$Gene.Name %in% plasmid_genes),]
genome_SYK6 <- genome_SYK6 %>% arrange(gene_start)

plasmid_SYK6 <- all_genes_SYK6[c(all_genes_SYK6$Gene.Name %in% plasmid_genes),]
plasmid_SYK6 <- plasmid_SYK6 %>% arrange(gene_start)

gene_data <- lapply(1:nrow(genome_SYK6), function(i){
  
  IS_gene_i <- all_IS$pos[which(genome_SYK6$gene_start[i] <= all_IS$pos &
                                  genome_SYK6$gene_end[i] >= all_IS$pos)]
  
  num_ins <- length(IS_gene_i)
  
  gene <- genome_SYK6$gene[i]
  if(is.null(gene)){
    gene <- NA
  }
  
  tibble(entry = genome_SYK6$Gene.Name[i],
         gene = gene,
         gene_start = genome_SYK6$gene_start[i],
         gene_end = genome_SYK6$gene_end[i],
         length = genome_SYK6$gene_end[i] - genome_SYK6$gene_start[i] +1,
         num_ins = num_ins,
         nondis = max(diff(c(genome_SYK6$gene_start[i], IS_gene_i, genome_SYK6$gene_end[i]))+1)/
           (genome_SYK6$gene_end[i] - genome_SYK6$gene_start[i] + 1),
         insdens = num_ins/(genome_SYK6$gene_end[i] - genome_SYK6$gene_start[i] + 1)
  )
  
})


gene_data <- do.call(rbind, gene_data)

gene_data$nondis[is.na(gene_data$nondis)] <- 1
gene_data$insdens[is.na(gene_data$insdens)] <- 0

gene_data$num_ins[is.na(gene_data$num_ins)] <- 0

gene_data$gene[is.na(gene_data$gene)] <- 
  paste("unnamed_gene_", seq_along(gene_data$gene[is.na(gene_data$gene)]),
        sep="")


p4 <- ggplot(gene_data,
             aes(x=round((gene_start+gene_end)/2), y=insdens))+
  geom_point(alpha=0.2, color="cyan4", size=0.5)+
  theme(legend.position = "none") +
  ylab("Gene-wise insertion density") +
  xlab("Genome postion") +
  geom_smooth(color="firebrick")+
  theme(legend.position="none",
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.title.y = element_text(size = rel(0.6))) +
  geom_hline(yintercept = nrow(all_IS)/length_genome, color="darkorange",linetype="dashed")

rm(gene_data)

### Simualted data

### Ecoli data preperation
# load data 
ecoli <- readxl::read_xls("./data_for_synthetic_data_generation/msb4100050-sup-0004.xls", skip = 2)

# we use the w3110 wildtype mutant (we only need the first 2 and 7:10 columns) 
# but could also use the mg1655 wildtype (columns 1:2, 7:10) 
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

# get min and max basepair of ECOLI (youe can set different ones for your data):
min_bp <- 1
max_bp <- 4646332 # see https://www.ncbi.nlm.nih.gov/nuccore/AP009048.1

# We know create a new dataframe based on the (in this case) the ecoli data
# this data set will be the standard form for our simulation, that is you can
# use any other kind of data set in the beginning or construct your own (different
# genes length etc) as long as you have the following columns

tnseqData <- ecoli
rm(ecoli)

### start generation of no-name dataset
# new gene names:
tnseqData$gene <- paste("gene_", 1:nrow(tnseqData), sep="")

names(tnseqData) <- c("UGI1", "gene", "UGI2", "left_bp", "right_bp") 
tnseqData$UGI1 <- paste("UGI1_", 1:nrow(tnseqData), sep="") 
tnseqData$UGI2 <- paste("UGI2_", 1:nrow(tnseqData), sep="") 

tnseqData <- tnseqData[,c("gene", "UGI1", "UGI2", "left_bp", "right_bp")]

tnseqData$gene.length <- tnseqData$right_bp - tnseqData$left_bp + 1 

num_genes <- nrow(tnseqData)

length_genome <- 4646332

source("functions.R")


##### p5

sim_run <- 1
unique_loci <- 100000
ess_ORF <- "uniform"
lambda <- 0.85
distortion <- "sin"
sine_scaling_factor <- 1.3
bp_per_wave <- 2000000
NegBinomial_num_cluster <- 100
NegBinomial_dispersion <- 1
NegBinomial_p <- 0.3
num_hot_spots <- 0
num_cold_spots <- 0
technical_noise <- 0.02

# set seed for drawing cluster sizes (NegBinomial_num_cluster samples) from
# a negative binomial distribution
set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
           num_hot_spots * 1000 + num_cold_spots * 1000000)
cluster.sizes.sim <- 
  rnbinom(NegBinomial_num_cluster, 
          NegBinomial_dispersion, 
          NegBinomial_p)+1. # +1 to avoid clusters of sizes 0

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
          all_my_probs[spots_data$spot_start[spot] : spots_data$spot_end[spot]] * 3
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

# Let only a fraction of size lambda of essential genes only have insertion sites
# with probabilty 0 (i.e. being insertion free)
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
      # set random start of essential part within a essential gene
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

# all negative probabilies due to the normal variance. term are set to 0
all_my_probs[all_my_probs < 0] <- 0

# normalize to get final probabilites for each potential insertion site
all_my_probs <- all_my_probs/sum(all_my_probs)

# sample observed insertion sites 
# we have to (1+0.05/(400000/unique_loci) to sample a bit more insertion sites to
# get ~unique insertion sites due to drawing with replacement (wothout replacement is 
# not possible for over 4,000,000 different probabilities in the sample function)
set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
           num_hot_spots * 1000 + num_cold_spots * 1000000)
observed_all_unique_loci <- 
  sort(
    unique(sample(1:max_bp, unique_loci * (1+0.05/(400000/unique_loci)), prob = all_my_probs, replace = T))
  ) 

# technical noise are insertion sites that can be everywhere on the genome:
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


tnseqData$num_ins <- ins.per.gene
tnseqData$ins_index <- tnseqData$num_ins/tnseqData$gene.length

all_IS <- tibble(pos=observed_all_unique_loci)

borders_genome <- seq(1, length_genome, 1000)

densities <- sapply(1:(length(borders_genome)-1), function(i){
  sum(
    all_IS$pos >= borders_genome[i] &
      all_IS$pos <= borders_genome[i+1])/
    (borders_genome[i+1]-borders_genome[i]+1)
})

gene_densities <- 
  tibble(pos=borders_genome[-length(borders_genome)],
         gene_insdens = densities)

p5 <- ggplot(tnseqData,
             aes(x=round((left_bp + right_bp)/2), y=ins_index))+
  geom_point(alpha=0.2, color="mediumpurple3", size=0.5)+
  theme(legend.position = "none") +
  ylab("Gene-wise insertion density") +
  xlab("Genome postion") +
  geom_smooth(color="firebrick")+
  theme(legend.position="none",
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.title.y = element_text(size = rel(0.6))) +
  geom_hline(yintercept = nrow(all_IS)/length_genome, color="darkorange",linetype="dashed")



#### p6

sim_run <- 2
unique_loci <- 400000
ess_ORF <- "uniform"
lambda <- 0.7
distortion <- "sin"
sine_scaling_factor <- 1.3
bp_per_wave <- 2000000
NegBinomial_num_cluster <- 100
NegBinomial_dispersion <- 1
NegBinomial_p <- 0.3
num_hot_spots <- 0
num_cold_spots <- 0
technical_noise <- 0

# set seed for drawing cluster sizes (NegBinomial_num_cluster samples) from
# a negative binomial distribution
set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
           num_hot_spots * 1000 + num_cold_spots * 1000000)
cluster.sizes.sim <- 
  rnbinom(NegBinomial_num_cluster, 
          NegBinomial_dispersion, 
          NegBinomial_p)+1. # +1 to avoid clusters of sizes 0

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
          all_my_probs[spots_data$spot_start[spot] : spots_data$spot_end[spot]] * 3
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

# Let only a fraction of size lambda of essential genes only have insertion sites
# with probabilty 0 (i.e. being insertion free)
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
      # set random start of essential part within a essential gene
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

# all negative probabilies due to the normal variance. term are set to 0
all_my_probs[all_my_probs < 0] <- 0

# normalize to get final probabilites for each potential insertion site
all_my_probs <- all_my_probs/sum(all_my_probs)

# sample observed insertion sites 
# we have to (1+0.05/(400000/unique_loci) to sample a bit more insertion sites to
# get ~unique insertion sites due to drawing with replacement (wothout replacement is 
# not possible for over 4,000,000 different probabilities in the sample function)
set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
           num_hot_spots * 1000 + num_cold_spots * 1000000)
observed_all_unique_loci <- 
  sort(
    unique(sample(1:max_bp, unique_loci * (1+0.05/(400000/unique_loci)), prob = all_my_probs, replace = T))
  ) 

# technical noise are insertion sites that can be everywhere on the genome:
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

tnseqData$num_ins <- ins.per.gene
tnseqData$ins_index <- tnseqData$num_ins/tnseqData$gene.length

all_IS <- tibble(pos=observed_all_unique_loci)


p6 <- ggplot(tnseqData,
             aes(x=round((left_bp + right_bp)/2), y=ins_index))+
  geom_point(alpha=0.2, color="mediumpurple3", size=0.5)+
  theme(legend.position = "none") +
  ylab("Gene-wise insertion density") +
  xlab("Genome postion") +
  geom_smooth(color="firebrick")+
  theme(legend.position="none",
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.title.y = element_text(size = rel(0.6))) +
  geom_hline(yintercept = nrow(all_IS)/length_genome, color="darkorange",linetype="dashed")




##### p7

sim_run <- 3
unique_loci <- 200000
ess_ORF <- "uniform"
lambda <- 0.8
distortion <- "uniform"
sine_scaling_factor <- 1.3
bp_per_wave <- 2000000
NegBinomial_num_cluster <- 100
NegBinomial_dispersion <- 1
NegBinomial_p <- 0.3
num_hot_spots <- 0
num_cold_spots <- 25
cold_spot_size <- 10000
hot_spot_size <- 10000
technical_noise <- 0

# set seed for drawing cluster sizes (NegBinomial_num_cluster samples) from
# a negative binomial distribution
set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
           num_hot_spots * 1000 + num_cold_spots * 1000000)
cluster.sizes.sim <- 
  rnbinom(NegBinomial_num_cluster, 
          NegBinomial_dispersion, 
          NegBinomial_p)+1. # +1 to avoid clusters of sizes 0

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
          all_my_probs[spots_data$spot_start[spot] : spots_data$spot_end[spot]] * 3
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

# Let only a fraction of size lambda of essential genes only have insertion sites
# with probabilty 0 (i.e. being insertion free)
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
      # set random start of essential part within a essential gene
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

# all negative probabilies due to the normal variance. term are set to 0
all_my_probs[all_my_probs < 0] <- 0

# normalize to get final probabilites for each potential insertion site
all_my_probs <- all_my_probs/sum(all_my_probs)

# sample observed insertion sites 
# we have to (1+0.05/(400000/unique_loci) to sample a bit more insertion sites to
# get ~unique insertion sites due to drawing with replacement (wothout replacement is 
# not possible for over 4,000,000 different probabilities in the sample function)
set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
           num_hot_spots * 1000 + num_cold_spots * 1000000)
observed_all_unique_loci <- 
  sort(
    unique(sample(1:max_bp, unique_loci * (1+0.05/(400000/unique_loci)), prob = all_my_probs, replace = T))
  ) 

# technical noise are insertion sites that can be everywhere on the genome:
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


tnseqData$num_ins <- ins.per.gene
tnseqData$ins_index <- tnseqData$num_ins/tnseqData$gene.length

all_IS <- tibble(pos=observed_all_unique_loci)


p7 <- ggplot(tnseqData,
             aes(x=round((left_bp + right_bp)/2), y=ins_index))+
  geom_point(alpha=0.2, color="mediumpurple3", size=0.5)+
  theme(legend.position = "none") +
  ylab("Gene-wise insertion density") +
  xlab("Genome postion") +
  geom_smooth(color="firebrick")+
  theme(legend.position="none",
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.title.y = element_text(size = rel(0.6))) +
  geom_hline(yintercept = nrow(all_IS)/length_genome, color="darkorange",linetype="dashed")


##### p8

sim_run <- 4
unique_loci <- 400000
ess_ORF <- "uniform"
lambda <- 0.7
distortion <- "uniform"
sine_scaling_factor <- 1.3
bp_per_wave <- 2000000
NegBinomial_num_cluster <- 100
NegBinomial_dispersion <- 1
NegBinomial_p <- 0.3
num_hot_spots <- 25
num_cold_spots <- 25
cold_spot_size <- 10000
hot_spot_size <- 10000
technical_noise <- 0.02


# set seed for drawing cluster sizes (NegBinomial_num_cluster samples) from
# a negative binomial distribution
set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
           num_hot_spots * 1000 + num_cold_spots * 1000000)
cluster.sizes.sim <- 
  rnbinom(NegBinomial_num_cluster, 
          NegBinomial_dispersion, 
          NegBinomial_p)+1. # +1 to avoid clusters of sizes 0

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
          all_my_probs[spots_data$spot_start[spot] : spots_data$spot_end[spot]] * 3
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

# Let only a fraction of size lambda of essential genes only have insertion sites
# with probabilty 0 (i.e. being insertion free)
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
      # set random start of essential part within a essential gene
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

# all negative probabilies due to the normal variance. term are set to 0
all_my_probs[all_my_probs < 0] <- 0

# normalize to get final probabilites for each potential insertion site
all_my_probs <- all_my_probs/sum(all_my_probs)

# sample observed insertion sites 
# we have to (1+0.05/(400000/unique_loci) to sample a bit more insertion sites to
# get ~unique insertion sites due to drawing with replacement (wothout replacement is 
# not possible for over 4,000,000 different probabilities in the sample function)
set.seed(sim_run + unique_loci + lambda *1000000 + NegBinomial_dispersion * 1000000 +
           num_hot_spots * 1000 + num_cold_spots * 1000000)
observed_all_unique_loci <- 
  sort(
    unique(sample(1:max_bp, unique_loci * (1+0.05/(400000/unique_loci)), prob = all_my_probs, replace = T))
  ) 

# technical noise are insertion sites that can be everywhere on the genome:
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


tnseqData$num_ins <- ins.per.gene
tnseqData$ins_index <- tnseqData$num_ins/tnseqData$gene.length

all_IS <- tibble(pos=observed_all_unique_loci)

p8 <- ggplot(tnseqData,
             aes(x=round((left_bp + right_bp)/2), y=ins_index))+
  geom_point(alpha=0.2, color="mediumpurple3", size=0.5)+
  theme(legend.position = "none") +
  ylab("Gene-wise insertion density") +
  xlab("Genome postion") +
  geom_smooth(color="firebrick")+
  theme(legend.position="none",
        strip.text = element_text(size=rel(0.6)),
        axis.text.x=element_text(size = rel(0.6)),
        axis.text.y=element_text(size = rel(0.6)),
        axis.title.x = element_text(size = rel(0.6)),
        axis.title.y = element_text(size = rel(0.6))) + 
  geom_hline(yintercept = nrow(all_IS)/length_genome, color="darkorange",linetype="dashed")


real_vs_simulated_data <- cowplot::plot_grid(
  p1, p2, p3, p4, p5, p6, p7, p8, 
  nrow = 4, byrow = T,
  labels = as.list(LETTERS[1:8]),
  label_size = 11,
  align = "v"
)


save_plot(filename = "./plots/real_vs_simulated_data.pdf",
          plot =  real_vs_simulated_data, dpi =600, base_height = 3.71*2, base_asp = 1.618/2)

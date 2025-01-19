# Functions used for the Simualtion study and real world application of the draft
# 'A novel analysis method and a sub-sampling tuning procedure to improve the 
# detection of essential genes in Tn5 libraries'

Binomial <- function(ins.positions, 
                     gene.names, 
                     gene.starts, 
                     gene.stops, 
                     num.ins.per.gene, 
                     genome.length, 
                     weighting=1){
  
  if(!length(unique(
    c(length(gene.names), 
      length(gene.starts), 
      length(gene.stops), 
      length(num.ins.per.gene))))==1){
    stop("Different lengths of gene.names, gene.starts, gene.stops and num.ins.per.gene")
  }
  
  ins_sites <- sort(ins.positions)
  
  observed_genome_insertion_densitiy <- 
    length(ins_sites)/genome.length
  
  results_per_gene <- lapply(seq_along(gene.names), function(i){
    
    gene_i_start <- 
      gene.starts[i]
    
    gene_i_stop <- 
      gene.stops[i]
    
    gene_i_length <- gene_i_stop - gene_i_start + 1
    gene_i_num_ins <- num.ins.per.gene[i]
    
    p_value <- sum(dbinom(0:gene_i_num_ins,
                          gene_i_length, 
                          observed_genome_insertion_densitiy * weighting))
    
    if(is.na(p_value)){
      p_value <- 1
    }
    
    tibble(gene = gene.names[i],
           p_value = p_value,
           weight_insdens = weighting)
    
    
  })
  
  results_per_gene <- do.call(rbind, results_per_gene)
  results_per_gene
  
  
}


ConNIS <- function(ins.positions, 
                   gene.names, 
                   gene.starts, 
                   gene.stops, 
                   num.ins.per.gene, 
                   genome.length, 
                   weighting=1){
  
  if(!length(unique(
    c(length(gene.names), 
      length(gene.starts), 
      length(gene.stops), 
      length(num.ins.per.gene))))==1){
    stop("Different lengths of gene.names, gene.starts, gene.stops and num.ins.per.gene")
  }
  
  ins_sites <- sort(ins.positions)
  
  observed_genome_insertion_densitiy <- 
    length(ins_sites)/genome.length
  
  results_per_gene <- lapply(seq_along(gene.names), function(i){
    
    gene_i_start <- 
      gene.starts[i]
    
    gene_i_stop <- 
      gene.stops[i]
    
    gene_i_length <- gene_i_stop - gene_i_start + 1
    
    expected_num_IS_gene_i <- 
      floor(gene_i_length * observed_genome_insertion_densitiy)
    
    gene_i_num_ins <- num.ins.per.gene[i]
    
    max_gap <- max(
      diff(
        unique(
          c(
            gene_i_start-1,
            ins_sites[ins_sites >= gene_i_start & ins_sites <= gene_i_stop],
            gene_i_stop+1
          ) 
        )
      )
    )
    max_gap <- min(gene_i_length, max_gap)
    
    if(max_gap > 1000){
      probs_bc <- 
        prob_seq_misses(
          gene_i_length, 
          gene_i_length-
            expected_num_IS_gene_i*weighting)
      p_value <- 
        sum(probs_bc[max_gap:(gene_i_length-
                                expected_num_IS_gene_i*weighting)])
      
      if(gene_i_length == gene_i_num_ins){
        p_value <- 1
      }
      if(max_gap >= gene_i_length-expected_num_IS_gene_i*weighting){
        p_value <- probs_bc[length(probs_bc)]
      }
      
    }else{
      
      if(gene_i_length == gene_i_num_ins){
        p_value <- 1
      }
      if(max_gap >= gene_i_length-expected_num_IS_gene_i*weighting){
        p_value <- 
          as.numeric(chooseZ(gene_i_length - (gene_i_length-expected_num_IS_gene_i*weighting) - 1,                            
                             gene_i_length- expected_num_IS_gene_i*weighting - 
                               (gene_i_length-expected_num_IS_gene_i*weighting)) /
          chooseZ(gene_i_length - 1, 
                             gene_i_length- expected_num_IS_gene_i*weighting - 1))
      }else{ # this is a spped-up for gap sizes of <=1000
        
        probs_bc <- 
          sapply(1: (max_gap-1), function(s){
            as.numeric(chooseZ(gene_i_length - s - 1, 
                               gene_i_length- expected_num_IS_gene_i*weighting - s) /
                         chooseZ(gene_i_length - 1, 
                                 gene_i_length- expected_num_IS_gene_i*weighting - 1))
          })
        
        p_value <- 1-sum(probs_bc)
      }
      
    }
    
    tibble(gene = gene.names[i],
           p_value = p_value,
           weighting = weighting)
    
  })
  
  results_per_gene <- do.call(rbind, results_per_gene)
  results_per_gene
}


ExpVsGamma <- function(gene.names, 
                       gene.starts, 
                       gene.stops, 
                       num.ins.per.gene,
                       log2threshold=2){
  
  if(!length(unique(
    c(length(gene.names), 
      length(gene.starts), 
      length(gene.stops))))==1){
    stop("Different lengths of gene.names, gene.starts and gene.stops")
  }
  
  gene_insdensity <- num.ins.per.gene/(gene.stops-gene.starts+1)
  
  
  # https://github.com/sanger-pathogens/Bio-Tradis/blob/master/bin/tradis_essentiality.R
  
  # original: ii <- STM_baseline$ins_index
  ii <- gene_insdensity
  
  
  #identify second maxima
  h <- hist(ii, breaks=200,plot=FALSE)
  maxindex <- which.max(h$density[10:length(h$density)])
  maxval <- h$mids[maxindex+3]
  
  # print pdf of loess curve and later on, histogram
  # pdf(paste(input, "QC_and_changepoint_plots", "pdf", sep = ".")) # we don't need a plot
  
  #find inter-mode minima with loess
  nG <- length(ii) # original: length(STM_baseline$read_count)
  r <- floor(maxval *2000)
  I = ii < r / 2000
  h1 = hist(ii[I],breaks=(0:r/2000), plot=FALSE)
  lo <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density
  #plot(h1$density, main="Density") # we don't need a plot
  #lines(predict(lo),col='red',lwd=2) # we don't need a plot
  m = h1$mids[which.min(predict(lo))]
  I1 = ((ii < m)&(ii >= 0))
  
  h = hist(ii, breaks="FD",plot=FALSE) 
  I2 = ((ii >= m)&(ii < h$mids[max(which(h$counts>5))]))
  f1 = (sum(I1) + sum(ii == 0))/nG
  f2 = (sum(I2))/nG
  
  if(all(I1==FALSE)){
    d1 <- list(
      estimate =NA,  
      sd = NA, 
      vcov = NA, 
      n = NA, 
      loglik = NA 
    )
  }else{
    d1 = fitdistr(ii[I1], "exponential")
  }
  
  if(all(I2==FALSE)){
    d2 <- list(
      estimate =NA,  
      sd = NA, 
      vcov = NA, 
      n = NA, 
      loglik = NA 
    )
  }else{
    d2 = fitdistr(ii[I2], "gamma") #fit curves
  }
  
  # print pdf of histogram
  #pdf("Loess_and_changepoint_estimation.pdf")
  
  #plots
  # hist(ii,breaks="FD", xlim=c(0,max(ii)), freq=FALSE,xlab="Insertion index", main="Gamma fits")
  # lines(0:200/500, f1*dgamma(0:200/500, 1, d1$estimate[1])) # was [2]
  # lines(0:200/500, f2*dgamma(0:200/500, d2$estimate[1], d2$estimate[2]))
  # print changepoint
  
  lower <- max(which(log((pgamma(1:1000/10000, d2$e[1],d2$e[2])*
                            (1-pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)))/
                           (pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)*
                              (1-pgamma(1:1000/10000, d2$e[1],d2$e[2]))) , base=2) < -log2threshold))
  
  upper <- min(which(log((pgamma(1:1000/10000, d2$e[1],d2$e[2])*
                            (1-pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)))/
                           (pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)*
                              (1-pgamma(1:1000/10000, d2$e[1],d2$e[2]))) , base=2) > log2threshold))
  
  essen <- lower/10000
  ambig <- upper/10000
  
  est_ess <- which(ii < essen)
  est_unclear <- which(ii >= essen & ii < ambig)
  est_non_ess <- seq_along(ii)[-c(est_ess, est_unclear)]
  
  est_essentiality <- rep("non-essential", length(ii))
  est_essentiality[est_ess] <- "essential"
  est_essentiality[est_unclear] <- "unclear"
  
  est_ess_exp_gamma <- as.numeric(est_essentiality == "essential")
  
  tibble(gene = gene.names,
         essential = est_ess_exp_gamma,
         log2threshold = log2threshold)
  
}


Geometric <- function(ins.positions, 
                      gene.names, 
                      gene.starts, 
                      gene.stops, 
                      num.ins.per.gene, 
                      genome.length, 
                      weighting=1){
  
  if(!length(unique(
    c(length(gene.names), 
      length(gene.starts), 
      length(gene.stops), 
      length(num.ins.per.gene))))==1){
    stop("Different lengths of gene.names, gene.starts, gene.stops and num.ins.per.gene")
  }
  
  ins_sites <- sort(ins.positions)
  
  observed_genome_insertion_densitiy <- 
    length(ins_sites)/genome.length
  
  results_per_gene <- lapply(seq_along(gene.names), function(i){
    
    gene_i_start <- 
      gene.starts[i]
    
    gene_i_stop <- 
      gene.stops[i]
    
    gene_i_length <- gene_i_stop - gene_i_start + 1
    
    expected_num_IS_gene_i <- 
      floor(gene_i_length * observed_genome_insertion_densitiy)
    
    gene_i_num_ins <- num.ins.per.gene[i]
    
    max_gap <- max(
      diff(
        unique(
          c(
            gene_i_start-1,
            ins_sites[ins_sites >= gene_i_start & ins_sites <= gene_i_stop],
            gene_i_stop+1
          ) 
        )
      )
    )
    max_gap <- min(gene_i_length, max_gap)
    
    p_value <- 1 - pgeom(max_gap, observed_genome_insertion_densitiy * weighting)
    
    tibble(gene = gene.names[i],
           p_value = p_value,
           weighting = weighting) 
  })
  
  results_per_gene <- do.call(rbind, results_per_gene)
  results_per_gene
}


Tn5Gaps <- function(ins.positions, 
                    gene.names, 
                    gene.starts, 
                    gene.stops, 
                    genome.length, 
                    weighting=1){
  
  if(!length(unique(
    c(length(gene.names), 
      length(gene.starts), 
      length(gene.stops))))==1){
    stop("Different lengths of gene.names, gene.starts and gene.stops")
  }
  
  ins_sites <- 
    sort(unique(
      c(1, 
        unique(ins.positions), 
        genome.length)))
  
  pins <- length(unique(ins.positions))/genome.length * weighting
  
  pnon <- 1 - pins
  
  results_per_gene <- lapply(seq_along(gene.names), function(i){
    
    gene_i_start <- 
      gene.starts[i]
    
    gene_i_stop <- 
      gene.stops[i]
    
    gene_i_length <- gene_i_stop - gene_i_start + 1
    
    lower_ins_index <- 
      max(
        which(ins_sites <= gene_i_start)
      )
    
    upper_ins_index <- 
      min(
        which(ins_sites >= gene_i_stop)
      )
    if(upper_ins_index == Inf){
      upper_ins_index <- lower_ins_index
    }
    
    gaps_gene_i <- diff(ins_sites[lower_ins_index:upper_ins_index])
    
    max_gap_gene_i <- unique(max(gaps_gene_i))          
    
    ins_sites_overlapping_gene_i <- 
      ins_sites[lower_ins_index:upper_ins_index]
    
    overlap_sizes_gene_i <- sapply(1:(length(ins_sites_overlapping_gene_i)-1), function(j){
      min(ins_sites_overlapping_gene_i[j+1], 
          gene_i_stop 
      ) -
        max(ins_sites_overlapping_gene_i[j], gene_i_start) -1
    })
    
    max_overlap_gene_i <- unique(max(overlap_sizes_gene_i))
    
    gap_size_for_p_value <- 
      gaps_gene_i[overlap_sizes_gene_i %in% max_overlap_gene_i][1]
    
    # pvalue from tn5gaps.pv of transit software
    # https://github.com/mad-lab/transit/blob/master/src/pytransit/analysis/tn5gaps.py
    p_value <- 
      1 - pgumbel(
        k = gap_size_for_p_value, 
        mu = log(length(ins_sites)* pins, 1/pnon), 
        s = 1/log(1/pnon))
    
    
    tibble(gene = gene.names[i],
           p_value = p_value,
           weighting = weighting)
    
  })
  
  results_per_gene <- do.call(rbind, results_per_gene)
  results_per_gene
}


#' uses GMP
# Function to get freqences of sequences of misses of lenght 1,2,...,s,....,k  
#' @param n : number of positions
#' @param k : number of misses
#' @return : number of possibilities to see sizes s, s = 1...t

freq_seq_misses <- function(n, k){
  
  lapply(1:k, function(s){ 
    (n-k+1) * chooseZ(n-s-1, k-s)
  })
  
}


# Function to get probabilities of sequences of misses of lenght 1,2,...,s,....,k 
#' @param n : number of positions
#' @param k : number of misses
#' @return : probability of cluster sizes s, s = 1...k

prob_seq_misses <- function(n, k){
  
  counts <- freq_seq_misses(n, k)
  
  sum_approx_freq <- Reduce("+", counts)
  prob_clusters_sizes <- 
    sapply(seq_along(counts), function(i){
      as.numeric(counts[[i]]/sum_approx_freq)
    }
    )
  prob_clusters_sizes
  
}


#' Function to get the Expectation of the cluster size in the bucket cluster problem 
#' @param n : number of positions
#' @param k : number of successes
#' @return : expectation value
seq_misses_expectation <- function(n, k){
  
   as.numeric(
     chooseZ(n, k-1) /
       chooseZ(n-1, k-1)
     )
    
}



#' Function to get the Variance of the of sequence of misses  
#' @param n : number of positions
#' @param k : number of successes
#' @return : expectation value
seq_misses_variance <- function(n, k){
  
  as.numeric(
    (chooseZ(n, k-1) + 2*  chooseZ(n, k-2)) /
      chooseZ(n-1, k-1)  - 
      (chooseZ(n, k-1)   /
         chooseZ(n-1, k-1) ) ^2
  )
  
}


#' Probability to randomly place t balls into n bins with no bin 
#' containing more than one ball, so that at least m consecutive bins are empty is
#' https://math.stackexchange.com/questions/3034654/arranging-n-balls-in-k-bins-so-that-m-consecutive-bins-are-empty
emptyBins <- function(n, k, m){
  
  choose(n, k)^-1 * 
    sum(
      sapply(1:floor(n/m), function(j){
        choose(n-m*j, k) * 
          choose(k+1, j) *
          (-1)^(j+1)
      })
    )
}

#' We have to fill in n-t to get the probabilty of having at least m consecutive
#' bins filled
# emptyBins(n, n-t, 3)




#### function for plotting possibilities
#' @param n
#' @param t

generate_plot_matrix <- function(n, k){
  
  if( n > 20){
    stop("Please you smaller n for plotting and exakt search")
  }
  
  as.matrix(
    do.call(rbind,
            lapply(combn(n, k, simplify = FALSE), function(i){
              as.numeric(c(1:n) %in% i)
            })
    )
  )
}
  


plot_combs <- function(n, k, s = NULL, m = NULL, sort.for.s = FALSE, truncate.non.s=FALSE){
  
  # m is a user specified matrix
  
  if(is.null(m)){
    plot_matrix <- generate_plot_matrix(n, k)
    # reverse order of rows
    plot_matrix <- plot_matrix[nrow(plot_matrix):1, ]
  }else{
    plot_matrix <- m
  }
  
  
  
  
  if(!is.null(s)){
    
    if(s > 1){
      
      count_list <- lapply(1:nrow(plot_matrix), function(i){
        
        actual_sequence <- plot_matrix[i,]
        
        out_sequence <- actual_sequence
        
        counter <- 0
        previous <- FALSE
        for(j in seq_along(actual_sequence)){
          
          if(actual_sequence[j] == 1 &
             j < length(actual_sequence)){
            counter <- counter + 1
            previous <- TRUE
            
          }else if(actual_sequence[j] == 1 &
                   j==length(actual_sequence) 
          ){
            if(previous == TRUE){
              out_sequence[(j):(j-counter)] <- counter+1
            }else{
              out_sequence[j] <- 1
            }
            
          }else if(actual_sequence[j] == 0 &
                   previous == TRUE){
            out_sequence[(j-1):(j-counter)] <- counter
            counter <- 0
          }else if(actual_sequence[j] == 0 &
                   previous == FALSE){
            counter <- 0
          }
          
        }
        out_sequence
        
        
      })
      
      count_matrix <- do.call(rbind, count_list)
      
      
      count_matrix[count_matrix != s & count_matrix > 0] <- 1
      
      count_matrix[count_matrix == s ] <- 2
      
    }else{
      
      count_list <- lapply(1:nrow(plot_matrix), function(i){
        
        actual_sequence <- plot_matrix[i,]
        
        out_sequence <- actual_sequence
        
        for(j in seq_along(actual_sequence)){

          if(j == 1 & 
             actual_sequence[j] == 1 &
             actual_sequence[j+1] == 0){
            out_sequence[j] <- 2
            
          }else if(j > 1 & 
                   j < length(actual_sequence)){
                     if(actual_sequence[j-1] == 0 &
                        actual_sequence[j] == 1 &
                        actual_sequence[j+1] == 0){
                       
                       out_sequence[j] <- 2
                       
                     }
                   
          
            
            
          }else if(j == length(actual_sequence)){
            if(actual_sequence[j] == 1 &
               actual_sequence[j-1] == 0){
              out_sequence[j] <- 2
            }
          } 
          
        }
        out_sequence
        
        
      })
      
      count_matrix <- do.call(rbind, count_list)
      
      
      
    }
    
    if(sort.for.s == TRUE){
      position_first_s_cluster <- sapply(1:nrow(count_matrix), function(i){
        first_found <- which(cumsum(count_matrix[i,] == 2) == 1)
        if(length(first_found) == 0){
          first_found <- 0
        }
        first_found
      })
      
      sort_matrix <- 
        cbind(count_matrix,
              position_first_s_cluster)
      
      
      
      count_matrix <- 
        matrix(sort_matrix[order(sort_matrix[,ncol(sort_matrix)], decreasing = T),
                           -ncol(sort_matrix)], ncol=ncol(count_matrix))
      
    }
    
    if(truncate.non.s == TRUE){
      kepp_combs <- 
        sapply(1:nrow(count_matrix),function(i){
          any(count_matrix[i,] == 2)})
      count_matrix<- count_matrix[kepp_combs,]
    }
    
    
    pic <- ggplot(melt(count_matrix), 
                 aes(x = Var2, y = Var1, fill = as.factor(value))) + 
      geom_tile(color="black") +
      scale_fill_manual(values = c("white", "red", "blue")) +
      theme_bw() +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()) + 
      coord_fixed()
    
    
    
  }else{
    
    pic <- ggplot(melt(plot_matrix == 1), 
           aes(x = Var2, y = Var1, fill = value)) + 
      geom_tile(color="black") +
      scale_fill_manual(values = c("white", "red")) +
      theme_bw() +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()) + 
      coord_fixed()
    
  }
  
  pic
  
  
}



true_freq_s <- function(n, k, s){
  
  plot_matrix <- generate_plot_matrix(n, k)
  
  out <- lapply(1:nrow(plot_matrix), function(i){
    
    csc(which(rep(1,n) == plot_matrix[i,]))
    
  })
  
  out <- do.call(c, out)
  
  freq_s <- sum(s == out)
  prob_s <- freq_s/sum(sapply(1:k, function(i){
    sum(i == out)
  }))
  tibble::tibble(
    n=n,
    k=k,
    s=s, 
    freq_s=freq_s,
    prob_s=prob_s)
}



# function for performanve evaluation

classification_performance <- function(est.ess, true.ess){
  
  TP <- sum(which(est.ess == 1) %in% 
              which(true.ess == 1) )
  FP <- sum(est.ess)-TP
  TN <- sum(which(est.ess == 0) %in% 
              which(true.ess == 0) )
  FN <- sum(true.ess == 1) - TP
  
  FPR <- FP/(FP+TN)
  TPR <- TP/(TP+FN)
  TNR <- TN/(TN+FP)
  Precision <- TP/(TP+FP)
  Accuracy <- (TP +TN) / (TP + TN + FP + FN)
  bACC <- (TPR +TNR)/2
  F1 <- 2*TP/(2*TP+FP+FN)
  MCC <- (TP*TN - FP*FN)/
    sqrt( (TP+FP) * (TP + FN) * (TN + FP ) * (TN + FN) )
  
  tibble(TP, FP, TN, FN, FPR, TPR, TNR, Precision, Accuracy, bACC, F1, MCC)
}


cluster_determiner <- function(x){
  out <- NULL
  counter <- 0
  for(i in 1:length(x)){
    if(x[i] == 1 ){
      counter <- counter + 1
      out <- c(out, counter)
      out[i:(i-counter+1)] <- counter
    }else{
      counter <- 0
      out <- c(out, 0)
    }
  }
  out
}


## function to generate scaling factors based on sine function
# bp.per.wave determines the number of base pairs for one complete wave
# scaling.factor determines the factor we wish to se for the highes upscaling 
# and downscaling; has to be >1 

sine_prob_scaling <- function(bp.per.wave, sine.scaling.factor){
  
  tmp_multiplier <- bp.per.wave/100000
  
  # from max to min length
  track_length <- 25000 * tmp_multiplier
  
  rate1 <- 0.00002 / tmp_multiplier
  rate2_range <- list(0, sine.scaling.factor - 1)
  
  x1 <- 1:(track_length)
  amp <- (rate2_range[[2]] - rate2_range[[1]])/2
  y1 <- amp*cos(2*pi*rate1*x1) + amp
  # plot(x1, y1, type='l')
  
  y_first <- y1+1
  y_second <- (1-(rev(y_first)-1) / sine.scaling.factor)
  y <- c(y_first, y_second)
  y <- c(y, rev(y))
  
  # plot(seq_along(y), y, type='l')
  y
}

#gumbel density
pgumbel <- function(k, mu, s){
  exp(
    - exp( 
      (mu - k) / s
    )
  ) 
}

# function to remove the all values before the first peak
remove_before_1st_peak <- function(stab_values, tuning_values, reverse=FALSE){
  
  # Reverse the order of stab_values and tuning_values if reverse is TRUE
  if(reverse == TRUE){
    stab_values <- rev(stab_values)
    tuning_values <- rev(tuning_values)
  }
  
  # Identify indices of NA in stab_values
  indices_NA <- which(is.na(stab_values))
  
  # Remove NA from stab_values and corresponding tuning_values
  if(length(indices_NA) > 0){
    stab_values_no_NA <- stab_values[-c(indices_NA)]
    tuning_values_no_NA <- tuning_values[-c(indices_NA)]
  }else{
    stab_values_no_NA <- stab_values
    tuning_values_no_NA <- tuning_values
  }
  
  # Calculate the differences between consecutive stability values
  diff_stab_values <- diff(stab_values_no_NA)
  
  if(all(diff_stab_values >= 0) | all(diff_stab_values <= 0)){
    
    out <- tibble(stab_values=as.numeric(stab_values_no_NA),
                  tuning_values=as.numeric(tuning_values_no_NA))
    
  }else{
    
    # Initialize dismiss_indices to store indices of non-monotonic values
    dismiss_indices <- 0
    for(i in seq_along(diff_stab_values)){
      if(diff_stab_values[i] >= 0 & max(dismiss_indices) == (i-1)){
        dismiss_indices <- c(dismiss_indices, i)
      }
    }
    
    # Remove non-monotonic values if any
    if(length(dismiss_indices) > 1){
      stab_values_leftpeak <- stab_values_no_NA[-dismiss_indices]
      tuning_values_leftpeak <- tuning_values_no_NA[-dismiss_indices]
    } else {
      stab_values_leftpeak <- stab_values_no_NA
      tuning_values_leftpeak <- tuning_values_no_NA
    }
    
    
    out <- tibble(stab_values=as.numeric(stab_values_leftpeak),
                  tuning_values=as.numeric(tuning_values_leftpeak))
  }
  
  out
  
}

# get the tunint value with the most stable results
get_most_stable <- function(stab_values, tuning_values, reverse = F){
  
  tmp <- remove_before_1st_peak(stab_values=stab_values, tuning_values=tuning_values, reverse = reverse)
  
  tmp_stab_values <- tmp$stab_values
  tmp_tuning_values <- tmp$tuning_values
  
  index_most_stable <- min(as.numeric(which(tmp_stab_values == min(tmp_stab_values))))
  
  out <- tmp[index_most_stable,]
  names(out) <- c("stability", "tuning_value")
  out 
}

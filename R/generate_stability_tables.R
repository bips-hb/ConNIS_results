# Best possible value vs selected
## BW25113
source("functions.R")

performance_metric <- "MCC"
beta <- 0.05

performances <- 
  readRDS("./performance/Performance_BW25113_realWorld.RDS")

adjust_bin <- "Benjamini-Hochberg"
adjust_con <- "Benjamini-Hochberg"
adjust_geo <- "Benjamini-Hochberg"
adjust_tn5 <- "Benjamini-Hochberg"

performances_all <- 
  performances %>%
  filter(read_count_threshold == 0 &
           trimming_start == 0.05 ) %>%
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
         N=TP+FP)

max_performance_metric <- 
  performances_all %>% 
  filter(
    
      (
        (adjustment == adjust_con & method == "ConNIS") | 
          (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
          (adjustment == adjust_geo & method == "Geometric") | 
          (adjustment == adjust_bin & method == "Binomial")| 
          is.na(adjustment)
      )
  ) %>%
  group_by(method) %>%
  slice(which.max(!!sym(performance_metric)))

min_performance_metric <- 
  performances_all %>% 
  filter(
      (
        (adjustment == adjust_con & method == "ConNIS") | 
          (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
          (adjustment == adjust_geo & method == "Geometric") | 
          (adjustment == adjust_bin & method == "Binomial")| 
          is.na(adjustment)
      )
  ) %>%
  group_by(method) %>%
  slice(which.min(!!sym(performance_metric)))

unweighted_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(1,3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(1,3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(2,3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(1,3)) |
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(0.1,3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(1,3)) 
      )
  )

stabilities <- readRDS("./performance/stabilities_bw25113.RDS")

stability_threshold <- 
  do.call(rbind,
          lapply(seq_along(stabilities), function(i){
            
            if(i %in% c(3,5)){
              stab <- get_most_stable(stabilities[[i]], as.numeric(names(stabilities[[i]])), 
                                       reverse = T)
            }else{
              
                if(i == 1){
                  if(adjust_bin == "Benjamini-Hochberg"){
                    stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                  }
                  if(adjust_bin == "Bonferroni"){
                    stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                  }
                  else{
                    stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                  }
                }
              
              if(i == 2){
                if(adjust_con == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_con == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 4){
                if(adjust_geo == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_geo == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 6){
                if(adjust_tn5 == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_tn5 == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
            }
            stab
          })
  )

stability_threshold$method  <- names(stabilities)

selected_threshold_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(stability_threshold$tuning_value[1],3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(stability_threshold$tuning_value[2],3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(stability_threshold$tuning_value[3],3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(stability_threshold$tuning_value[4],3)) | 
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(stability_threshold$tuning_value[5],3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(stability_threshold$tuning_value[6],3)) 
      )
  )

stabiltiy_BW25113 <-
  as_tibble(
    rbind(
      paste("$", 
            as.numeric(unlist(round(selected_threshold_performance_metric[,performance_metric],2))),
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(max_performance_metric[,performance_metric],2))), 
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(min_performance_metric[,performance_metric],2))),
            "$", sep=""),
      paste("$", 
            c(
              as.character(unlist(round(unweighted_performance_metric[,performance_metric],2)))), 
            "$",
            sep="")
    )
  )
  

names(stabiltiy_BW25113) <- 
  c("Binomial",
    "ConNIS", 
    "Exp. vs. Gamma", 
    "Geometric", 
    "InsDens", 
    "Tn5Gaps")

rm(list=c("performances", 
          "stability_threshold", 
          "max_performance_metric",
          "min_performance_metric",
          "performances_all",
          "selected_threshold_performance_metric",
          "unweighted_performance_metric",
          "stabilities", "adjust_bin", "adjust_con", "adjust_geo", "adjust_tn5"))

## MG1655

performances <- 
  readRDS("./performance/Performance_MG1655_realWorld.RDS")

adjust_bin <- "Benjamini-Hochberg"
adjust_con <- "Bonferroni"
adjust_geo <- "Bonferroni"
adjust_tn5 <- "Benjamini-Hochberg"

performances_all <- 
  performances %>%
  filter(read_count_threshold == 0 &
           trimming_start == 0.05 ) %>%
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
         N=TP+FP)

max_performance_metric <- 
  performances_all %>% 
  filter(
      (
        (adjustment == adjust_con & method == "ConNIS") | 
          (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
          (adjustment == adjust_geo & method == "Geometric") | 
          (adjustment == adjust_bin & method == "Binomial")| 
          is.na(adjustment)
      )
  ) %>%
  group_by(method) %>%
  slice(which.max(!!sym(performance_metric)))

min_performance_metric <- 
  performances_all %>% 
  filter(
      (
        (adjustment == adjust_con & method == "ConNIS") | 
          (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
          (adjustment == adjust_geo & method == "Geometric") | 
          (adjustment == adjust_bin & method == "Binomial")| 
          is.na(adjustment)
      )
  ) %>%
  group_by(method) %>%
  slice(which.min(!!sym(performance_metric)))

unweighted_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(1,3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(1,3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(2,3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(1,3)) |
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(0.1,3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(1,3)) 
      )
  )

stabilities <- readRDS("./performance/stabilities_MG1655.RDS")

stability_threshold <- 
  do.call(rbind,
          lapply(seq_along(stabilities), function(i){
            
            if(i %in% c(3,5)){
              stab <- get_most_stable(stabilities[[i]], as.numeric(names(stabilities[[i]])), 
                                      reverse = T)
            }else{
              
              if(i == 1){
                if(adjust_bin == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_bin == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 2){
                if(adjust_con == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_con == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 4){
                if(adjust_geo == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_geo == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 6){
                if(adjust_tn5 == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_tn5 == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
            }
            stab
          })
  )

stability_threshold$method  <- names(stabilities)

selected_threshold_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(stability_threshold$tuning_value[1],3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(stability_threshold$tuning_value[2],3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(stability_threshold$tuning_value[3],3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(stability_threshold$tuning_value[4],3)) | 
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(stability_threshold$tuning_value[5],3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(stability_threshold$tuning_value[6],3)) 
      )
  )

stabiltiy_MG1655 <-
  as_tibble(
    rbind(
      paste("$", 
            as.numeric(unlist(round(selected_threshold_performance_metric[,performance_metric],2))),
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(max_performance_metric[,performance_metric],2))), 
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(min_performance_metric[,performance_metric],2))),
            "$", sep=""),
      paste("$", 
            c(
              as.character(unlist(round(unweighted_performance_metric[,performance_metric],2)))), 
            "$",
            sep="")
    )
  )
  
names(stabiltiy_MG1655) <- 
  c("Binomial",
    "ConNIS", 
    "Exp. vs. Gamma", 
    "Geometric", 
    "InsDens", 
    "Tn5Gaps")

rm(list=c("performances", 
          "stability_threshold", 
          "max_performance_metric",
          "min_performance_metric",
          "performances_all",
          "selected_threshold_performance_metric",
          "unweighted_performance_metric",
          "stabilities", "adjust_bin", "adjust_con", "adjust_geo", "adjust_tn5"))

## 14028s

performances <- 
  readRDS("./performance/Performance_14028S_realWorld.RDS")

adjust_bin <- "Benjamini-Hochberg"
adjust_con <- "Bonferroni"
adjust_geo <- "Bonferroni"
adjust_tn5 <- "Benjamini-Hochberg"

performances_all <- 
  performances %>%
  filter(read_count_threshold == 1 &
           trimming_start == 0.05 ) %>%
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
         N=TP+FP)

max_performance_metric <- 
  performances_all %>% 
  filter(
      (
        (adjustment == adjust_con & method == "ConNIS") | 
          (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
          (adjustment == adjust_geo & method == "Geometric") | 
          (adjustment == adjust_bin & method == "Binomial")| 
          is.na(adjustment)
      )
  ) %>%
  group_by(method) %>%
  slice(which.max(!!sym(performance_metric)))

min_performance_metric <- 
  performances_all %>% 
  filter(
      (
        (adjustment == adjust_con & method == "ConNIS") | 
          (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
          (adjustment == adjust_geo & method == "Geometric") | 
          (adjustment == adjust_bin & method == "Binomial")| 
          is.na(adjustment)
      )
  ) %>%
  group_by(method) %>%
  slice(which.min(!!sym(performance_metric)))

unweighted_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(1,3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(1,3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(2,3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(1,3)) |
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(0.1,3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(1,3)) 
      )
  )

stabilities <- readRDS("./performance/stabilities_14028s.RDS")

stability_threshold <- 
  do.call(rbind,
          lapply(seq_along(stabilities), function(i){
            
            if(i %in% c(3,5)){
              stab <- get_most_stable(stabilities[[i]], as.numeric(names(stabilities[[i]])), 
                                      reverse = T)
            }else{
              
              if(i == 1){
                if(adjust_bin == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_bin == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 2){
                if(adjust_con == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_con == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 4){
                if(adjust_geo == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_geo == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 6){
                if(adjust_tn5 == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_tn5 == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
            }
            stab
          })
  )

stability_threshold$method  <- names(stabilities)

selected_threshold_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(stability_threshold$tuning_value[1],3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(stability_threshold$tuning_value[2],3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(stability_threshold$tuning_value[3],3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(stability_threshold$tuning_value[4],3)) | 
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(stability_threshold$tuning_value[5],3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(stability_threshold$tuning_value[6],3)) 
      )
  )


stabiltiy_14028S <-
  as_tibble(
    rbind(
      paste("$", 
            as.numeric(unlist(round(selected_threshold_performance_metric[,performance_metric],2))),
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(max_performance_metric[,performance_metric],2))), 
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(min_performance_metric[,performance_metric],2))),
            "$", sep=""),
      paste("$", 
            c(
              as.character(unlist(round(unweighted_performance_metric[,performance_metric],2)))), 
            "$",
            sep="")
    )
  )

names(stabiltiy_14028S) <- 
  c("Binomial",
    "ConNIS", 
    "Exp. vs. Gamma", 
    "Geometric", 
    "InsDens", 
    "Tn5Gaps")

rm(list=c("performances", 
          "stability_threshold", 
          "max_performance_metric",
          "min_performance_metric",
          "performances_all",
          "selected_threshold_performance_metric",
          "unweighted_performance_metric",
          "stabilities", "adjust_bin", "adjust_con", "adjust_geo", "adjust_tn5"))

## example simu1
performances <- 
  readRDS("./performance/Performance_example_simu1.RDS")

adjust_bin <- "Benjamini-Hochberg"
adjust_con <- "Bonferroni-Holm"
adjust_geo <- "Bonferroni-Holm"
adjust_tn5 <- "Benjamini-Hochberg"

performances_all <- 
  performances %>%
  filter(trimming_start == 0.05) %>%
  group_by(method,  
           post_prob, 
           threshold_log2_expgamma, 
           weighting,
           adjustment) %>%
  mutate(tuning = 
           sum(c(
             threshold_log2_expgamma, 
             weighting, 
             post_prob), na.rm=T),
         N=TP+FP)

max_performance_metric <- 
  performances_all %>% 
  filter(
      (
        (adjustment == adjust_con & method == "ConNIS") | 
          (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
          (adjustment == adjust_geo & method == "Geometric") | 
          (adjustment == adjust_bin & method == "Binomial")| 
          is.na(adjustment)
      )
  ) %>%
  group_by(method) %>%
  slice(which.max(!!sym(performance_metric)))

min_performance_metric <- 
  performances_all %>% 
  filter(
      (
        (adjustment == adjust_con & method == "ConNIS") | 
          (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
          (adjustment == adjust_geo & method == "Geometric") | 
          (adjustment == adjust_bin & method == "Binomial")| 
          is.na(adjustment)
      )
  ) %>%
  group_by(method) %>%
  slice(which.min(!!sym(performance_metric)))

unweighted_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(1,3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(1,3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(2,3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(1,3)) |
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(0.1,3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(1,3)) 
      )
  )

stabilities <- readRDS("./performance/stabilities_example_simu1.RDS")


stability_threshold <- 
  do.call(rbind,
          lapply(seq_along(stabilities), function(i){
            
            if(i %in% c(3,5)){
              stab <- get_most_stable(stabilities[[i]], as.numeric(names(stabilities[[i]])), 
                                      reverse = T)
            }else{
              
              if(i == 1){
                if(adjust_bin == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_bin == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 2){
                if(adjust_con == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_con == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 4){
                if(adjust_geo == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_geo == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 6){
                if(adjust_tn5 == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_tn5 == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
            }
            stab
          })
  )

stability_threshold$method  <- names(stabilities)

selected_threshold_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(stability_threshold$tuning_value[1],3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(stability_threshold$tuning_value[2],3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(stability_threshold$tuning_value[3],3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(stability_threshold$tuning_value[4],3)) | 
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(stability_threshold$tuning_value[5],3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(stability_threshold$tuning_value[6],3)) 
      )
  )

stabiltiy_example_simu1 <-
  as_tibble(
    rbind(
      paste("$", 
            as.numeric(unlist(round(selected_threshold_performance_metric[,performance_metric],2))),
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(max_performance_metric[,performance_metric],2))), 
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(min_performance_metric[,performance_metric],2))),
            "$", sep=""),
      paste("$", 
            c(
              as.character(unlist(round(unweighted_performance_metric[,performance_metric],2)))), 
            "$",
            sep="")
    )
  )

names(stabiltiy_example_simu1) <- 
  c("Binomial",
    "ConNIS", 
    "Exp. vs. Gamma", 
    "Geometric", 
    "InsDens", 
    "Tn5Gaps")

rm(list=c("performances", 
          "stability_threshold", 
          "max_performance_metric",
          "min_performance_metric",
          "performances_all",
          "selected_threshold_performance_metric",
          "unweighted_performance_metric",
          "stabilities", "adjust_bin", "adjust_con", "adjust_geo", "adjust_tn5"))

## example simu2
performances <- 
  readRDS("./performance/Performance_example_simu2.RDS")

adjust_bin <- "Benjamini-Hochberg"
adjust_con <- "Bonferroni-Holm"
adjust_geo <- "Bonferroni-Holm"
adjust_tn5 <- "Benjamini-Hochberg"

performances_all <- 
  performances %>%
  filter(trimming_start == 0.05) %>%
  group_by(method,  
           post_prob, 
           threshold_log2_expgamma, 
           weighting,
           adjustment) %>%
  mutate(tuning = 
           sum(c(
             threshold_log2_expgamma, 
             weighting, 
             post_prob), na.rm=T),
         N=TP+FP)

max_performance_metric <- 
  performances_all %>% 
  filter(
    (
      (adjustment == adjust_con & method == "ConNIS") | 
        (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
        (adjustment == adjust_geo & method == "Geometric") | 
        (adjustment == adjust_bin & method == "Binomial")| 
        is.na(adjustment)
    )
  ) %>%
  group_by(method) %>%
  slice(which.max(!!sym(performance_metric)))

min_performance_metric <- 
  performances_all %>% 
  filter(
    (
      (adjustment == adjust_con & method == "ConNIS") | 
        (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
        (adjustment == adjust_geo & method == "Geometric") | 
        (adjustment == adjust_bin & method == "Binomial")| 
        is.na(adjustment)
    )
  ) %>%
  group_by(method) %>%
  slice(which.min(!!sym(performance_metric)))

unweighted_performance_metric <- 
  performances_all %>% 
  filter(
    (
      
      (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(1,3))|
        (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(1,3)) | 
        (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(2,3)) |
        (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(1,3)) |
        (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(0.1,3)) |
        (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(1,3)) 
    )
  )

stabilities <- readRDS("./performance/stabilities_example_simu2.RDS")


stability_threshold <- 
  do.call(rbind,
          lapply(seq_along(stabilities), function(i){
            
            if(i %in% c(3,5)){
              stab <- get_most_stable(stabilities[[i]], as.numeric(names(stabilities[[i]])), 
                                      reverse = T)
            }else{
              
              if(i == 1){
                if(adjust_bin == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_bin == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 2){
                if(adjust_con == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_con == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 4){
                if(adjust_geo == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_geo == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 6){
                if(adjust_tn5 == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_tn5 == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
            }
            stab
          })
  )

stability_threshold$method  <- names(stabilities)

selected_threshold_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(stability_threshold$tuning_value[1],3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(stability_threshold$tuning_value[2],3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(stability_threshold$tuning_value[3],3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(stability_threshold$tuning_value[4],3)) | 
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(stability_threshold$tuning_value[5],3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(stability_threshold$tuning_value[6],3)) 
      )
  )

stabiltiy_example_simu2 <-
  as_tibble(
    rbind(
      paste("$", 
            as.numeric(unlist(round(selected_threshold_performance_metric[,performance_metric],2))),
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(max_performance_metric[,performance_metric],2))), 
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(min_performance_metric[,performance_metric],2))),
            "$", sep=""),
      paste("$", 
            c(
              as.character(unlist(round(unweighted_performance_metric[,performance_metric],2)))), 
            "$",
            sep="")
    )
  )



names(stabiltiy_example_simu2) <- 
  c("Binomial",
    "ConNIS", 
    "Exp. vs. Gamma", 
    "Geometric", 
    "InsDens", 
    "Tn5Gaps")

rm(list=c("performances", 
          "stability_threshold", 
          "max_performance_metric",
          "min_performance_metric",
          "performances_all",
          "selected_threshold_performance_metric",
          "unweighted_performance_metric",
          "stabilities", "adjust_bin", "adjust_con", "adjust_geo", "adjust_tn5"))

## example simu3
performances <- 
  readRDS("./performance/Performance_example_simu3.RDS")

adjust_bin <- "Benjamini-Hochberg"
adjust_con <- "Bonferroni-Holm"
adjust_geo <- "Bonferroni-Holm"
adjust_tn5 <- "Benjamini-Hochberg"

performances_all <- 
  performances %>%
  filter(trimming_start == 0.05) %>%
  group_by(method,  
           post_prob, 
           threshold_log2_expgamma, 
           weighting,
           adjustment) %>%
  mutate(tuning = 
           sum(c(
             threshold_log2_expgamma, 
             weighting, 
             post_prob), na.rm=T),
         N=TP+FP)

max_performance_metric <- 
  performances_all %>% 
  filter(
    (
      (adjustment == adjust_con & method == "ConNIS") | 
        (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
        (adjustment == adjust_geo & method == "Geometric") | 
        (adjustment == adjust_bin & method == "Binomial")| 
        is.na(adjustment)
    )
  ) %>%
  group_by(method) %>%
  slice(which.max(!!sym(performance_metric)))

min_performance_metric <- 
  performances_all %>% 
  filter(
    (
      (adjustment == adjust_con & method == "ConNIS") | 
        (adjustment == adjust_tn5 & method == "Tn5Gaps") | 
        (adjustment == adjust_geo & method == "Geometric") | 
        (adjustment == adjust_bin & method == "Binomial")| 
        is.na(adjustment)
    )
  ) %>%
  group_by(method) %>%
  slice(which.min(!!sym(performance_metric)))

unweighted_performance_metric <- 
  performances_all %>% 
  filter(
    (
      
      (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(1,3))|
        (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(1,3)) | 
        (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(2,3)) |
        (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(1,3)) |
        (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(0.1,3)) |
        (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(1,3)) 
    )
  )

stabilities <- readRDS("./performance/stabilities_example_simu3.RDS")


stability_threshold <- 
  do.call(rbind,
          lapply(seq_along(stabilities), function(i){
            
            if(i %in% c(3,5)){
              stab <- get_most_stable(stabilities[[i]], as.numeric(names(stabilities[[i]])), 
                                      reverse = T)
            }else{
              
              if(i == 1){
                if(adjust_bin == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_bin == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 2){
                if(adjust_con == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_con == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 4){
                if(adjust_geo == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_geo == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
              if(i == 6){
                if(adjust_tn5 == "Benjamini-Hochberg"){
                  stab <- get_most_stable(stabilities[[i]][[1]], as.numeric(names(stabilities[[i]][[1]])))
                }
                if(adjust_tn5 == "Bonferroni"){
                  stab <- get_most_stable(stabilities[[i]][[2]], as.numeric(names(stabilities[[i]][[2]])))
                }
                else{
                  stab <- get_most_stable(stabilities[[i]][[3]], as.numeric(names(stabilities[[i]][[3]])))
                }
              }
              
            }
            stab
          })
  )

stability_threshold$method  <- names(stabilities)

selected_threshold_performance_metric <- 
  performances_all %>% 
  filter(
      (
        
        (adjustment == adjust_bin & method == "Binomial" & round(tuning,3) == round(stability_threshold$tuning_value[1],3))|
          (adjustment == adjust_con & method == "ConNIS" & round(tuning,3) == round(stability_threshold$tuning_value[2],3)) | 
          (is.na(adjustment) & method == "Exp. vs. Gamma" & round(tuning,3) == round(stability_threshold$tuning_value[3],3)) |
          (adjustment == adjust_geo & method == "Geometric" & round(tuning,3) == round(stability_threshold$tuning_value[4],3)) | 
          (is.na(adjustment) & method == "InsDens" & round(tuning,3) == round(stability_threshold$tuning_value[5],3)) |
          (adjustment == adjust_tn5 & method == "Tn5Gaps" & round(tuning,3) == round(stability_threshold$tuning_value[6],3)) 
      )
  )


stabiltiy_example_simu3 <-
  as_tibble(
    rbind(
      paste("$", 
            as.numeric(unlist(round(selected_threshold_performance_metric[,performance_metric],2))),
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(max_performance_metric[,performance_metric],2))), 
            "$",sep=""),
      paste("$", 
            as.numeric(unlist(round(min_performance_metric[,performance_metric],2))),
            "$", sep=""),
      paste("$", 
            c(
              as.character(unlist(round(unweighted_performance_metric[,performance_metric],2)))), 
            "$",
            sep="")
    )
  )



names(stabiltiy_example_simu3) <- 
  c("Binomial",
    "ConNIS", 
    "Exp. vs. Gamma", 
    "Geometric", 
    "InsDens", 
    "Tn5Gaps")

rm(list=c("performances", 
          "stability_threshold", 
          "max_performance_metric",
          "min_performance_metric",
          "performances_all",
          "selected_threshold_performance_metric",
          "unweighted_performance_metric",
          "stabilities", "adjust_bin", "adjust_con", "adjust_geo", "adjust_tn5"))

write_csv(
  bind_rows(
    stabiltiy_BW25113,
    stabiltiy_MG1655,
    stabiltiy_14028S,
    stabiltiy_example_simu1,
    stabiltiy_example_simu2,
    stabiltiy_example_simu3
  ), file = paste("table_stabilities", performance_metric, ".csv", sep=""), 
  col_names = T
)

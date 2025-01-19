# get gene data for Salmonella Typhimurium (S. Typhimurium) 14028s 
# https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.01723/full

library(tidyverse)
library(stringr)
genes_and_positions <- read.csv("./14028s_data/genes_and_positions.txt", header = F)


indices <- which(str_detect(genes_and_positions[,1], "  gene  ") )
next_lines <- diff(c(indices, nrow(genes_and_positions)))-1

genes_ <- lapply(seq_along(indices), function(i){
  
  act_index <- indices[i]
  tmp_value <- genes_and_positions[act_index,]
  
 
    
    tmp_value_no_gene <- str_remove_all(tmp_value, "gene")
    tmp_value_no_ws <- str_remove_all(tmp_value_no_gene, " ")
    tmp_value_no_comp <- str_remove_all(tmp_value_no_ws, "complement\\(")
    tmp_value_no_par <- str_remove_all(tmp_value_no_comp, "\\)")
    
    gene_start <- as.numeric(str_split_1(tmp_value_no_par, pattern = "\\.\\.")[1])
    gene_end <- as.numeric(str_split_1(tmp_value_no_par, pattern = "\\.\\.")[2])
    
    
    out_gene <- sapply((act_index+1) : (act_index+next_lines[i]), function(j){
      
      tmp_value2 <- genes_and_positions[j,]
      
      gene <- locus_tag <- gene_syn <- NA
      
      if(str_detect(tmp_value2, "/gene=")){
        tmp_value2_no_gene <- str_remove_all(tmp_value2, "/gene=")
        tmp_value2_no_ws <- str_remove_all(tmp_value2_no_gene, " ")
        gene <- tmp_value2_no_ws
      }
      gene
      
    } )
    
    gene <- unique(out_gene)
    if(length(gene)>1 & any(is.na(gene))){
      gene <- gene[!is.na(gene)]
    }
    
    out_locus_tag <- sapply((act_index+1) : (act_index+next_lines[i]), function(j){
      
      tmp_value2 <- genes_and_positions[j,]
      
      gene <- locus_tag <- gene_syn <- NA
      
      if(str_detect(tmp_value2, "/locus_tag=")){
        tmp_value2_no_gene <- str_remove_all(tmp_value2, "/locus_tag=")
        tmp_value2_no_ws <- str_remove_all(tmp_value2_no_gene, " ")
        locus_tag <- tmp_value2_no_ws
      }
      locus_tag
      
    } )
    locus_tag <- unique(out_locus_tag)
    if(length(locus_tag)>1 & any(is.na(locus_tag))){
      locus_tag <- locus_tag[!is.na(locus_tag)]
    }
    
    
    out_gene_syn <- sapply((act_index+1) : (act_index+next_lines[i]), function(j){
      
      tmp_value2 <- genes_and_positions[j,]
      
      gene <- locus_tag <- gene_syn <- NA
      
      if(str_detect(tmp_value2, "/gene_synonym=")){
        tmp_value2_no_gene <- str_remove_all(tmp_value2, "/gene_synonym=")
        tmp_value2_no_ws <- str_remove_all(tmp_value2_no_gene, " ")
        gene_syn <- tmp_value2_no_ws
      }
      
      gene_syn
      
    } )
    gene_syn <- unique(out_gene_syn)
    if(length(gene_syn)>1 & any(is.na(gene_syn))){
      gene_syn <- gene_syn[!is.na(gene_syn)]
    }
    
    tibble(gene, gene_syn, locusId = locus_tag, gene_start, gene_end)
  
  
})

genes <- do.call(rbind, genes)

saveRDS(genes, "./14028s_data/Typhimurium14028S_genes.RDS")

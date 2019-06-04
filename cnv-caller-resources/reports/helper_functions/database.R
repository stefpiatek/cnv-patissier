library(DBI)
library(RSQLite)
library(tidyverse)
library(stringr)



get_callers <- function(con){
  callers <- con %>%
    tbl("callers") %>%
    collect() %>%
    mutate(name = str_remove(name, "_case")) %>%
    mutate(name =str_replace(name, "cnvkit", "CNV kit")) %>%
    mutate(name =str_replace(name, "codex2", "CODEX2")) %>%
    mutate(name =str_replace(name, "copywriter", "CopywriteR")) %>%
    mutate(name =str_replace(name, "decon", "DECoN")) %>%
    mutate(name =str_replace(name, "exome-depth", "ExomeDepth")) %>%
    mutate(name =str_replace(name, "excavator2", "EXCAVATOR2")) %>%
    mutate(name =str_replace(name, "gatk", "GATK gCNV")) %>%
    mutate(name =str_replace(name, "panelcn_mops", "panelcn.MOPS")) %>%
    mutate(name =str_replace(name, "savvycnv", "SavvyCNV")) %>%
    mutate(name =str_replace(name, "xhmm", "XHMM")) 
    
  
  
  return(callers)
}

get_callers_list <- function(db_file){
  con <- dbConnect(SQLite(), dbname= db_file)
  callers <- get_callers(con) 
  dbDisconnect(con)
  
  return(callers)
}

get_true_positives <- function(db_file){
  con <- dbConnect(SQLite(), dbname= db_file)
  callers <- get_callers(con)
  
  genes <- con %>%
    tbl("genes") 
  
  positives <- genes %>%
    left_join(tbl(con, "samples"), by = c("id" = "gene_id"), suffix = c("_gene", "")) %>%
    filter(result_type == "positive")
  
  true_positives <- positives %>%
    left_join(tbl(con, "called_cnvs"), by=c("id" = "sample_id"),  suffix = c("", "_called_cnv")) %>%
    collect() %>%
    mutate(sample_id = id) %>%
    filter(!(capture == "MiniPanel" & name_gene == "BRCA1" & start == 41197636 & end == 41276193 & cnv_id == 66)) %>%
    filter(!(capture == "MiniPanel" & name_gene == "BRCA2" & start == 32890510 & end == 32972962 & cnv_id == 115)) %>%
    select(-id) %>%
    group_by(caller_id, sample_id, name_gene) %>%
    slice(which.min(id_called_cnv)) %>%
    ungroup() %>%
    left_join(callers, by = c("caller_id" = "id"), suffix = c("", "_caller"))
  
  dbDisconnect(con)
  
  return(true_positives)
}

get_positives <- function(db_file){
  con <- dbConnect(SQLite(), dbname= db_file)
  genes <- con %>%
    tbl("genes") 
  
  positives <- genes %>%
    left_join(tbl(con, "samples"), by = c("id" = "gene_id"), suffix = c("_gene", "")) %>%
    filter(result_type == "positive") %>%
    collect() 
  
  dbDisconnect(con)
  

  return(positives)
}

get_false_positives <- function(db_file){
  con <- dbConnect(SQLite(), dbname= db_file)
  
  callers <- get_callers(con) 
  
  genes <- con %>%
    tbl("genes") 
  
  negatives <-  genes %>%
    left_join(tbl(con, "samples"), by = c("id" = "gene_id"), suffix = c("_gene", "")) %>%
    filter(result_type == "normal") 
  
  false_positives <- negatives %>%
    left_join(tbl(con, "called_cnvs"), by=c("id" = "sample_id"),  suffix = c("", "_called_cnv")) %>%
    collect() %>%
    mutate(sample_id = id) %>%
    filter(!(capture == "MiniPanel" & name_gene == "BRCA1" & start == 41197636 & end == 41276193 & cnv_id == 66)) %>%
    filter(!(capture == "MiniPanel" & name_gene == "BRCA2" & start == 32890510 & end == 32972962 & cnv_id == 115)) %>%
    filter(!(capture == "ICR" & name_gene != "BRCA1" & name_gene != "BRCA2")) %>%
    select(-id) %>%
    group_by(caller_id, sample_id, name_gene) %>%
    slice(which.min(id_called_cnv)) %>%
    ungroup() %>%
    left_join(callers, by = c("caller_id" = "id"), suffix = c("", "_caller")) 

  dbDisconnect(con)
  

  return(false_positives)
}

get_negatives <- function(db_file){
  con <- dbConnect(SQLite(), dbname= db_file)
  genes <- con %>%
    tbl("genes") 
  
  negatives <-  genes %>%
    left_join(tbl(con, "samples"), by = c("id" = "gene_id"), suffix = c("_gene", "")) %>%
    filter(result_type == "normal") %>% 
    filter(!(capture == "ICR" & name_gene != "BRCA1" & name_gene != "BRCA2")) %>%
    collect() 
  
  dbDisconnect(con)
  
  return(negatives)
}




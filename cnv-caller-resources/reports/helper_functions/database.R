library(DBI)
library(RSQLite)
library(tidyverse)
library(stringr)



get_callers <- function(con){
  callers <- con %>%
    tbl("callers") %>%
    collect() %>%
    mutate(name = str_remove(name, "_case")) 
  
  
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
    filter(!(capture == "MiniPanel" & name_gene == "BRCA1" & cnv_id == 37)) %>%
    filter(!(capture == "MiniPanel" & name_gene == "BRCA2" & cnv_id == 115)) %>%
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

get_false_negatives <- function(db_file){
  con <- dbConnect(SQLite(), dbname= db_file)
  
  callers <- get_callers(con) 
  
  genes <- con %>%
    tbl("genes") 
  
  negatives <-  genes %>%
    left_join(tbl(con, "samples"), by = c("id" = "gene_id"), suffix = c("_gene", "")) %>%
    filter(result_type == "normal") 
  
  false_negatives <- negatives %>%
    left_join(tbl(con, "called_cnvs"), by=c("id" = "sample_id"),  suffix = c("", "_called_cnv")) %>%
    collect() %>%
    mutate(sample_id = id) %>%
    filter(!(capture == "MiniPanel" & name_gene == "BRCA1" & cnv_id == 37)) %>%
    filter(!(capture == "MiniPanel" & name_gene == "BRCA2" & cnv_id == 115)) %>%
    select(-id) %>%
    group_by(caller_id, sample_id, name_gene) %>%
    slice(which.min(id_called_cnv)) %>%
    ungroup() %>%
    left_join(callers, by = c("caller_id" = "id"), suffix = c("", "_caller")) 

  dbDisconnect(con)
  

  return(false_negatives)
}

get_negatives <- function(db_file){
  con <- dbConnect(SQLite(), dbname= db_file)
  genes <- con %>%
    tbl("genes") 
  
  negatives <-  genes %>%
    left_join(tbl(con, "samples"), by = c("id" = "gene_id"), suffix = c("_gene", "")) %>%
    filter(result_type == "normal") %>% 
    collect() 
  
  dbDisconnect(con)
  
  return(negatives)
}




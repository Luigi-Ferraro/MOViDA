library(readr)
library(tidyverse)
library(readxl)
library(httr)
library(jsonlite)
library(furrr)
library(magrittr)
library(yaGST)
library(org.Hs.eg.db)
library(annotate)
plan(multisession, workers = 100)


get_enrich <- function(experiment,expression_data, GO){
  rankedList <- expression_data[[experiment]]
  names(rankedList) <- expression_data$GENE_SYMBOLS
  rankedList <- sort(rankedList,decreasing = T)
  ans <- GO %>% map_dbl(possibly(function(x) mwwGST(rankedList, x, minLenGeneSet = 0, verbose = FALSE)$nes, 0))
  #ans <- lapply(GO, function(x) mwwGST(rankedList, x, minLenGeneSet = 0, verbose = FALSE))
  ans
}


get_enrich_p <- function(experiment,expression_data,GO){
  rankedList <- expression_data[[experiment]]
  names(rankedList) <- expression_data$GENE_SYMBOLS
  rankedList <- sort(rankedList,decreasing = T)
  ans <- GO %>% map_dbl(possibly(function(x) mwwGST(rankedList, x, minLenGeneSet = 0, verbose = FALSE)$p.value, 0))
  ans
}



#gene_ontologies <- read_delim("gene_ontologies.txt", "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
#expression_data <- read_delim("/storage/qnap_vol1/data/PUBLIC/CCLE/expression/GDSC_Cell_line_RMA_proc_basalExp.csv", 
#                                                                  "\t", escape_double = FALSE, trim_ws = TRUE)
gene_ontologies <- read_csv("code_prep_data/gene_ontologies.txt", 
                            col_names = FALSE)
expression_data <- as.data.frame(dataNorm)
expression_data$GENE_SYMBOLS <- rownames(expression_data)

fun <- possibly(.f = function(x) getSYMBOL(get(x, org.Hs.egGO2ALLEGS),"org.Hs.eg"),otherwise = NA)
GO2 <- gene_ontologies$X1 %>% map(fun)
names(GO2) <- gene_ontologies$X1

mapped <- GO2 %>% map_dbl(length) %>% enframe()
all(mapped$value > 0)

enrich_df <- colnames(expression_data)[-ncol(expression_data)] %>% future_map_dfr(get_enrich, expression_data, GO2)
enrich_df$DepMap_ID <- colnames(expression_data)[-ncol(expression_data)]
enrich_df <- enrich_df %>% dplyr::select(DepMap_ID,everything())

enrich_df <- as.data.frame(enrich_df)
names(enrich_df)[names(enrich_df) == colnames(enrich_df)[1]] <- "CCLE_Name"
#enrich_df <- enrich_df %>% rename(CCLE_Name = colnames(enrich_df)[1])
write.table(enrich_df,"./data_ge_ont/GO_enrich_gbm.tsv", sep = "\t", col.names = T, row.names = F, quote = F)


# enrich_p_df <- colnames(expression_data)[-c(1:2)] %>% future_map_dfr(get_enrich_p,expression_data,GO2)
# enrich_p_df$DepMap_ID <- colnames(expression_data)[-c(1:2)]
# enrich_p_df <- enrich_p_df %>% dplyr::select(DepMap_ID,everything())
# 
# enrich_p_df_ann <- enrich_p_df %>% left_join(dplyr::select(all_pairs,X1,DepMap_ID)) %>% dplyr::select(X1,DepMap_ID,everything())
# write.table(enrich_p_df_ann,"GO_enrich_p_ann_GDSC.tab", sep = "\t", col.names = T, row.names = F, quote = F)


save.image(file = "GDSC_GO.Rdata")

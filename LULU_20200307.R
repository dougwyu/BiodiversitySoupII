## ------------------------------------------------------------------------------
# set some general options
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) # if you want to pass arguments to R from the shell (bash) command line
print(args)

otutablefilename <- args[1]


## ------------------------------------------------------------------------------
library(tidyverse)
# devtools::install_github("tobiasgf/lulu", force = TRUE)
library(lulu)
library(knitr)


## ---- eval=FALSE, include=FALSE------------------------------------------------
## setwd("~/Xiaoyangmiseqdata/BiodiversitySoupII/data/seqs/folderA/Filter_min1PCRs_min1copies_A")
## otutablefilename <- "../data/seqs/folderA/Filter_min1PCRs_min1copies_A/table_BioSoupII_A_97.txt"
## matchlistfilename <- "../data/seqs/folderA/Filter_min1PCRs_min1copies_A/matchlist.txt"


## ------------------------------------------------------------------------------
otutable <- read_tsv(otutablefilename) %>% 
    rename(OTU_ID = "#OTU ID") %>% 
    column_to_rownames(var = "OTU_ID")

matchlist <- read_tsv("matchlist.txt", col_names = FALSE)


## ------------------------------------------------------------------------------
curated_result <- lulu(otutable, matchlist)
curated_table <- rownames_to_column(curated_result$curated_table, var = "OTU_ID")
write_tsv(curated_table, path = "table_BioSoupII_97_lulu.txt")


## ----eval=FALSE, include=FALSE-------------------------------------------------
## setwd("~/Xiaoyangmiseqdata/BiodiversitySoupII/scripts")
## purl("LULU_20200307.Rmd")


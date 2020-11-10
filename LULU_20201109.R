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


## ---- eval=FALSE, include=FALSE------------------------------------------------
## test code
# pcr <- "A"
# setwd(glue::glue("~/Xiaoyangmiseqdata/BiodiversitySoupII_repo/data/seqs/folder{pcr}/Filter_min1PCRs_min2copies_{pcr}"))
# otutablefilename <- glue::glue("table_BioSoupII_{pcr}_97.txt")
# matchlistfilename <- "matchlist.txt"

## ------------------------------------------------------------------------------
otutable <- read_tsv(otutablefilename) %>% 
    rename(OTU_ID = "#OTU ID") %>% 
    column_to_rownames(var = "OTU_ID")

matchlist <- read.table("matchlist.txt", header=FALSE,
                        as.is=TRUE, stringsAsFactors=FALSE)

## ------------------------------------------------------------------------------
curated_result <- lulu(otutable, matchlist, minimum_match = 97.1) # default minimum_match = 84% similarity, but such a low value only works if there are many samples, so that otu co-occurrence gives a strong signal that two otus are parent and child.

curated_table <- rownames_to_column(curated_result$curated_table, 
                                    var = "OTU_ID"
                                    )

write_tsv(curated_table, file = "table_BioSoupII_97_lulu.txt")


## ----eval=FALSE, include=FALSE-------------------------------------------------
## setwd("~/Xiaoyangmiseqdata/BiodiversitySoupII/scripts")
## knitr::purl("LULU_20200307.Rmd")


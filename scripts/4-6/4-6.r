 
library(readr)

# importing results of Panther

setwd("/mnt/hdd0/univ/thesis-support-material/data/4-6/all_vs_50kb/panther/output/allDEgenes/raw/")

BP_all <- read_delim("BP.txt", 
                     "\t", escape_double = FALSE, col_types = cols(`allDEgenes.txt (FDR)` = col_double()), 
                     trim_ws = TRUE, skip = 6)

MF_all <- read_delim("MF.txt", 
                     "\t", escape_double = FALSE, col_types = cols(`allDEgenes.txt (FDR)` = col_double()), 
                     trim_ws = TRUE, skip = 6)

CC_all <- read_delim("CC.txt", 
                     "\t", escape_double = FALSE, col_types = cols(`allDEgenes.txt (FDR)` = col_double()), 
                     trim_ws = TRUE, skip = 6)

BP_50KB <- read_delim("../../DE50KB/raw/BP.txt", 
                      "\t", escape_double = FALSE, trim_ws = TRUE, skip = 6)

library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(annotables)

library(readr)
DE_tp53_sh1 <- read_csv("sh1.txt", 
                        col_names = FALSE)

DE_LINC <- read_delim("lincko.txt", 
                      "\t", escape_double = FALSE, col_names = FALSE,
                      trim_ws = TRUE)

universe <- grch37$ensgene

ego <- enrichGO(gene = as.character(DE_LINC$X1), 
                universe = as.character(universe),
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "CC", # BP, MF, CC
                pvalueCutoff = 1E-3,
                pAdjustMethod = "none", 
                readable = TRUE)

cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "~/univ/lab2/05-10/clusterProfiler/linc_ko/pval0001/CC.csv") # BP, MF, CC

dotplot(ego, showCategory=500)

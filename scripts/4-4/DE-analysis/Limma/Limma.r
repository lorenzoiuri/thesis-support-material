
basedir <- "/mnt/hdd0/univ/thesis-support-material/data/"

library(edgeR)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)

# GUIDE:
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# importing data
setwd(paste0(basedir, "4-4/featureCounts-output/"))

filename = "counts.txt"
#filename = "filtered-only-lncrnas/counts_only_lncRNA.txt"

counts_full <- read_delim(filename, "\t", 
                          escape_double = FALSE, trim_ws = TRUE, skip = 1)

counts_simple_full <- as.matrix(counts_full[,c(7,8,9,10,11,12)])

# renaming the rows and columns

# note: bc*-1 are the shNT samples
#       bc*-2 are the sh1 samples
#       bc*-3 are the sh2 samples

colnames(counts_simple_full) <- c("bc1-1", "bc1-2", "bc1-3", "bc2-1", "bc2-2", "bc2-3")
rownames(counts_simple_full) <- counts_full$Geneid

# set adjusted pvalue and LFC cutoffs
padj.cutoff <- 0.01
lfc.cutoff <- 0.58

# ----- shNT vs sh-1 ----

# keeping only shNT and sh1 samples
counts <- counts_simple_full[,c(1,4,2,5)]

condition <- factor(c(rep("control", 2), rep("sh1", 2)))

sampleInfo <- data.frame(row.names = colnames(counts), condition)

group <- sampleInfo$condition

# creating object
d0 <- DGEList(counts)

# compute normalization factors
d0 <- calcNormFactors(d0)

# filtering low expressed genes
# we remove the genes whose cpm normalized count is < 1 
# in every sample
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = F)

# model fitting
fit <- lmFit(y, mm)

# setting the comparison between control and sh1
contr <- makeContrasts(groupsh1 - groupcontrol, levels = colnames(coef(fit)))

# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

result.table <- topTable(tmp, sort.by = "P", n = Inf)

DE.genes.sh1.limma <- result.table %>%
                        filter(adj.P.Val < padj.cutoff &
                               abs(logFC) > lfc.cutoff ) %>%
                        arrange(-logFC) %>%
                        mutate(gene = rownames(.)) %>%
                        select(gene, adj.P.Val, logFC)

DE.genes.sh1.limma <- as_tibble(DE.genes.sh1.limma)

# -----------------------
rm(counts)
rm(condition)
rm(sampleInfo)
rm(group)
rm(d0)
rm(d)
rm(mm)
rm(drop)
rm(y)
rm(fit)
rm(contr)
rm(tmp)
rm(result.table)

# ----- shNT vs sh2 ----

# keeping only shNT and sh2 samples
counts <- counts_simple_full[,c(1,4,3,6)]

condition <- factor(c(rep("control", 2), rep("sh2", 2)))

sampleInfo <- data.frame(row.names = colnames(counts), condition)

group <- sampleInfo$condition

# creating object
d0 <- DGEList(counts)

# compute normalization factors
d0 <- calcNormFactors(d0)

# filtering low expressed genes
# we remove the genes whose cpm normalized count is < 1 
# in every sample
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = F)

# model fitting
fit <- lmFit(y, mm)

# setting the comparison between control and sh1
contr <- makeContrasts(groupsh2 - groupcontrol, levels = colnames(coef(fit)))

# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

result.table <- topTable(tmp, sort.by = "P", n = Inf)

DE.genes.sh2.limma <- result.table %>%
    filter(adj.P.Val < padj.cutoff &
               abs(logFC) > lfc.cutoff ) %>%
    arrange(-logFC) %>%
    mutate(gene = rownames(.)) %>%
    select(gene, adj.P.Val, logFC)

DE.genes.sh2.limma <- as_tibble(DE.genes.sh2.limma)

# ----------------------------------

# a value of LFC > 0 means that the gene is MORE EXPRESSED in the sh1 (or sh2)
# samples (the gene is upregulated upon tp53-mut silencing)

rm(counts)
rm(condition)
rm(sampleInfo)
rm(group)
rm(d0)
rm(d)
rm(mm)
rm(drop)
rm(y)
rm(fit)
rm(contr)
rm(tmp)
rm(result.table)

rm(counts_full)
rm(counts_simple_full)
rm(cutoff)
rm(lfc.cutoff)
rm(padj.cutoff)
rm(filename)

# -----------------------------------

#DE.genes.sh1.limma.upreg <- DE.genes.sh1.limma %>% filter(logFC > 0)
#DE.genes.sh1.limma.downreg <- DE.genes.sh1.limma %>% filter(logFC < 0)

#DE.genes.sh2.limma.upreg <- DE.genes.sh2.limma %>% filter(logFC > 0)
#DE.genes.sh2.limma.downreg <- DE.genes.sh2.limma %>% filter(logFC < 0)


basedir <- "/mnt/hdd0/univ/thesis-support-material/data/"

library(edgeR)
library(dplyr)
library(tibble)
library(readr)
library(HTSFilter)

# GUIDE used:
# http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html

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

# ----- shNT vs sh1 ----

# keeping only shNT and sh1 samples
counts <- counts_simple_full[,c(1,4,2,5)]

condition <- factor(c(rep("control", 2), rep("sh1", 2)))

sampleInfo <- data.frame(row.names = colnames(counts), condition)

# creating edger data object
dgeFull <- DGEList(counts, group=sampleInfo$condition)

# remove genes with zero counts for all samples
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)

# estimate the normalization factors
dgeFull <- calcNormFactors(dgeFull, method="TMM")

# estimate common and tagwise dispersion
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)

# perform an exact test for the difference in expression between the two conditions
# and get top DE genes
dgeTest <- exactTest(dgeFull)
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))

# perform independent filtering, removing poorly expressed genes to
# reduce problem of multiple testing
# note: DESEQ2 does this filtering automatically
filtData <- HTSFilter(dgeFull)$filteredData
dgeTestFilt <- exactTest(filtData)
resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))

# keep only the DE genes that satisfy the cutoffs
DE.genes.sh1.edger <- resFilt$table[resFilt$table$FDR < padj.cutoff &
                            abs(resFilt$table$logFC) > lfc.cutoff,]

DE.genes.sh1.edger <- DE.genes.sh1.edger[order(DE.genes.sh1.edger$logFC),]

DE.genes.sh1.edger <- DE.genes.sh1.edger %>%
                        mutate(gene = rownames(.)) %>%
                        select(gene, FDR, logFC, PValue)

DE.genes.sh1.edger <- as_tibble(DE.genes.sh1.edger)

# -----------------------

rm(counts)
rm(dgeFull)
rm(dgeTest)
rm(dgeTestFilt)
rm(filtData)
rm(resNoFilt)
rm(resFilt)
rm(sampleInfo)
rm(condition)

# ----- shNT vs sh2 ----

# keeping only shNT and sh2 samples
counts <- counts_simple_full[,c(1,4,3,6)]

condition <- factor(c(rep("control", 2), rep("sh2", 2)))

sampleInfo <- data.frame(row.names = colnames(counts), condition)

# creating edger data object
dgeFull <- DGEList(counts, group=sampleInfo$condition)

# remove genes with zero counts for all samples
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)

# estimate the normalization factors
dgeFull <- calcNormFactors(dgeFull, method="TMM")

# estimate common and tagwise dispersion
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)

# perform an exact test for the difference in expression between the two conditions
# and get top DE genes
dgeTest <- exactTest(dgeFull)
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))

# perform independent filtering, removing poorly expressed genes to
# reduce problem of multiple testing
# note: DESEQ2 does this filtering automatically
filtData <- HTSFilter(dgeFull)$filteredData
dgeTestFilt <- exactTest(filtData)
resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))

# keep only the DE genes that satisfy the cutoffs
DE.genes.sh2.edger <- resFilt$table[resFilt$table$FDR < padj.cutoff &
                                 abs(resFilt$table$logFC) > lfc.cutoff,]

DE.genes.sh2.edger <- DE.genes.sh2.edger[order(DE.genes.sh2.edger$logFC),]

DE.genes.sh2.edger <- DE.genes.sh2.edger %>%
                        mutate(gene = rownames(.)) %>%
                        select(gene, FDR, logFC, PValue)

DE.genes.sh2.edger <- as_tibble(DE.genes.sh2.edger)

# ------------------------------

# a value of LFC > 0 means that the gene is MORE EXPRESSED in the sh1 (or sh2)
# samples (the gene is upregulated upon tp53-mut silencing)

rm(counts)
rm(dgeFull)
rm(dgeTest)
rm(dgeTestFilt)
rm(filtData)
rm(resNoFilt)
rm(resFilt)
rm(sampleInfo)
rm(condition)

rm(counts_full)
rm(counts_simple_full)
rm(lfc.cutoff)
rm(padj.cutoff)
rm(filename)

# ------------------------------

#DE.genes.sh1.edger.upreg <- DE.genes.sh1.edger %>% filter(logFC > 0)
#DE.genes.sh1.edger.downreg <- DE.genes.sh1.edger %>% filter(logFC < 0)

#DE.genes.sh2.edger.upreg <- DE.genes.sh2.edger %>% filter(logFC > 0)
#DE.genes.sh2.edger.downreg <- DE.genes.sh2.edger %>% filter(logFC < 0)

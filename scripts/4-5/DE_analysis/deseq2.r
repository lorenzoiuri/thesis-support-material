
library(dplyr)
library(tibble)
library(DESeq2)
library(readr)
library(ggrepel)

# importing data

setwd("/mnt/hdd0/univ/thesis-support-material/data/4-5/question2/1-DE_histones/")

counts_full <- read_delim("counts.txt", "\t", 
                          escape_double = FALSE, trim_ws = TRUE, skip = 1)

counts_simple_full <- as.matrix(counts_full[,c(7,8,9,10,11,12)])

colnames(counts_simple_full) <- c("bc1-1", "bc1-2", "bc1-3", "bc2-1", "bc2-2", "bc2-3")
rownames(counts_simple_full) <- counts_full$Geneid

padj.cutoff <- 0.05
lfc.cutoff <- 0.26

# ----- shNT vs sh-1 ----

counts <- counts_simple_full[,c(1,4,2,5)]

condition <- factor(c(rep("control", 2), rep("sh1", 2)))

meta <- data.frame(row.names = colnames(counts), condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~condition)

dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)

contrast.sh1 <- c("condition", "sh1", "control")

res_table_sh1_unshrunken <- results(dds, contrast = contrast.sh1, alpha = padj.cutoff)
res_table_sh1_shrunken <- lfcShrink(dds, contrast = contrast.sh1, res=res_table_sh1_unshrunken, type = "normal")

res_table_tb_sh1 <- res_table_sh1_shrunken %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()

DE.genes.sh1.deseq2 <- res_table_tb_sh1 %>%
                        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>%
                        arrange(-abs(log2FoldChange))

DE.genes.sh1.deseq2 <- DE.genes.sh1.deseq2 %>%
                        select(gene, padj, log2FoldChange, pvalue)

#res_table_tb_sh1 %>%
#    filter(!is.na(padj)) %>%
#    mutate(DE = (padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)) %>%
#    mutate(DEALOT = (-log10(padj) > 50 & abs(log2FoldChange) > 2)) %>%
#    mutate(text = ifelse(DEALOT, gene, NA)) %>%
#    ggplot(mapping = aes(x=log2FoldChange, 
#                         y=-log10(padj), 
#                         col=DE,
#                         label=text)) +
#        geom_point() + 
#        theme_minimal() +
#        geom_text_repel() +
#        scale_color_manual(values=c("blue", "black", "red")) +
#        geom_vline(xintercept=c(-0.58, 0.58), col="red") +
#        geom_hline(yintercept=-log10(0.01), col="red") + 
#        ggtitle("Differentially expressed genes in shNT vs sh-1, by Deseq2")

# ----------------------

rm(dds)
rm(meta)
rm(counts)
rm(normalized_counts)
rm(condition)
rm(res_table_sh1_shrunken)
rm(res_table_sh1_unshrunken)
rm(res_table_tb_sh1)
rm(contrast.sh1)

# ----- shNT vs sh-2 ----

counts <- counts_simple_full[,c(1,4,3,6)]

condition <- factor(c(rep("control", 2), rep("sh2", 2)))

meta <- data.frame(row.names = colnames(counts), condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~condition)

dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)

contrast.sh2 <- c("condition", "sh2", "control")

res_table_sh2_unshrunken <- results(dds, contrast = contrast.sh2, alpha = padj.cutoff)
res_table_sh2_shrunken <- lfcShrink(dds, contrast = contrast.sh2, res=res_table_sh2_unshrunken, type = "normal")

res_table_tb_sh2 <- res_table_sh2_shrunken %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()

DE.genes.sh2.deseq2 <- res_table_tb_sh2 %>%
    filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>%
    arrange(-abs(log2FoldChange))

DE.genes.sh2.deseq2 <- DE.genes.sh2.deseq2 %>%
    select(gene, padj, log2FoldChange, pvalue)

# --------------------------------

rm(dds)
rm(meta)
rm(counts)
rm(normalized_counts)
rm(condition)
rm(res_table_sh2_shrunken)
rm(res_table_sh2_unshrunken)
rm(res_table_tb_sh2)
rm(contrast.sh2)

rm(counts_full)
rm(counts_simple_full)
rm(lfc.cutoff)
rm(padj.cutoff)

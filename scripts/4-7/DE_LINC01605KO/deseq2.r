
library(dplyr)
library(tibble)
library(DESeq2)
library(readr)


# importing data

setwd("/mnt/hdd0/univ/thesis-support-material/data/4-7/DE_genes")
counts_full <- read_delim("featurecounts/counts_LINC01605KO.txt", "\t", 
                          escape_double = FALSE, trim_ws = TRUE, 
                          skip = 1)

counts_simple_full <- as.matrix(counts_full[,c(7,8,9,10,11,12,13,14)])

colnames(counts_simple_full) <- c("IR1", "IR2", "IR4", "IR5", "IR6", "IR7", "IR8", "IR10")
rownames(counts_simple_full) <- counts_full$Geneid

padj.cutoff <- 0.01
lfc.cutoff <- 0.58

#---------------------------------

counts <- counts_simple_full

condition <- factor(c(rep("WT", 4), rep("KO", 4)))

meta <- data.frame(row.names = colnames(counts), condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~condition)

dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)

write.table(normalized_counts, file = "deseq2_normalized_counts.txt",
            row.names = T, col.names = T, quote = F, sep = "\t")

contrast <- c("condition", "KO", "WT")

res_table_unshrunken <- results(dds, contrast = contrast, alpha = padj.cutoff)
res_table_shrunken <- lfcShrink(dds, contrast = contrast, res=res_table_unshrunken)

res_table_tb <- res_table_shrunken %>%
                    data.frame() %>%
                    rownames_to_column(var="gene") %>%
                    as_tibble()

DE.genes.deseq2 <- res_table_tb %>%
                        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>%
                        arrange(-abs(log2FoldChange))

DE.genes.deseq2 <- DE.genes.deseq2 %>%
                        select(gene, padj, log2FoldChange, pvalue)

# ----------------------

rm(dds)
rm(meta)
rm(counts)
rm(normalized_counts)
rm(condition)
rm(res_table_shrunken)
rm(res_table_unshrunken)
rm(res_table_tb)
rm(contrast)

rm(counts_full)
rm(counts_simple_full)
rm(lfc.cutoff)
rm(padj.cutoff)

# ---------------------

write.table(DE.genes.deseq2$gene, file = "DE_DESeq2_LINC01605ko.txt", 
            row.names = F, col.names = F, sep = "\t", quote = F)

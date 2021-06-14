
basedir <- "/mnt/hdd0/univ/thesis-support-material/data/"

library(dplyr)
library(tibble)
library(DESeq2)
library(readr)
library(ggrepel)

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

# setting the value for the filtering parameters

padj.cutoff <- 0.01
lfc.cutoff <- 0.58

# ----- shNT vs sh1 ----

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

# --- volcano plot locating the genes of the GO term "regulation of nuclear division"

# temp <- res_table_tb_sh1 %>% mutate(gene = substr(gene, 0, 15))
# 
# regulation_nuclear_division <- read_csv("../aux/regulation_nuclear_division.txt",
#                                         col_names = FALSE)
# 
# gene_symbol_Enseml_mapping <- read_delim("../../4-5/gene-symbol-Enseml-mapping.txt",
#                                          "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
# 
# a <- temp %>%
#             filter(!is.na(padj)) %>%
#             mutate(DE = (padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)) %>%
#             mutate(target = ifelse(gene %in% regulation_nuclear_division$X1 , T,F)) %>%
#             inner_join(gene_symbol_Enseml_mapping, by = c("gene" = "X1")) %>%
#             mutate(text = ifelse(target, X2, NA))
# a %>%
# ggplot(mapping = aes(x=log2FoldChange,
#                          y=-log10(padj),
#                          shape=DE,
#                          col=target,
#                          label = text)) +
#         geom_point(alpha = 0.7) +
#         geom_point(data = (a %>% filter(gene %in% regulation_nuclear_division$X1))) +
#         theme_minimal() +
#         geom_text_repel(max.overlaps = 35) +
#         scale_color_manual(values=c("lightblue", "plum4")) +
#         xlim(c(-2,4.2)) +
#         ylim(c(0, 75)) +
#         theme(legend.position = "none") +
#         geom_vline(xintercept=c(-1*lfc.cutoff, lfc.cutoff), col="grey") +
#         geom_hline(yintercept=-log10(padj.cutoff), col="grey") +
#         labs(title = "Differentially expressed genes involved in the regulation of nuclear division process",
#              x = "log2 fold change",
#              y = "-log10(pvalue)")
# 
# rm(temp, a, regulation_nuclear_division, gene_symbol_Enseml_mapping)

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

# ----- shNT vs sh2 ----

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

# a value of LFC > 0 means that the gene is MORE EXPRESSED in the sh-1 (or sh-2)
# samples (the gene is upregulated upon tp53-mut silencing)
# example:
# gene ENSG00000039068.19_8 has LFC > 5, 
# the raw counts for that gene are:
#
# bc1-1 bc1-2 bc1-3 bc2-1 bc2-2 bc2-3 
#    37  2612  1547    36  3093  1759 

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
rm(filename)

# ----------------------------------

#DE.genes.sh1.deseq2.upreg <- DE.genes.sh1.deseq2 %>% filter(log2FoldChange > 0)
#DE.genes.sh1.deseq2.downreg <- DE.genes.sh1.deseq2 %>% filter(log2FoldChange < 0)

#DE.genes.sh2.deseq2.upreg <- DE.genes.sh2.deseq2 %>% filter(log2FoldChange > 0)
#DE.genes.sh2.deseq2.downreg <- DE.genes.sh2.deseq2 %>% filter(log2FoldChange < 0)

#setwd("/mnt/hdd0/univ/lab2/rnaseq/DE_analysis/Deseq2")
#write.table(DE.genes.sh1.deseq2, file = "DE_sh1_DESEQ2.txt",
#            col.names = T, row.names = F, sep = "\t", quote = F)
#write.table(DE.genes.sh2.deseq2, file = "DE_sh2_DESEQ2.txt",
#            col.names = T, row.names = F, sep = "\t", quote = F)

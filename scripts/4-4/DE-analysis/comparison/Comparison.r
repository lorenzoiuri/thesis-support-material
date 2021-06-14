
basedir <- "/mnt/hdd0/univ/thesis-support-material/data/"

# it's necessary to have execute the script DESeq2.r, EdgeR.r and Limma.r and to have kept the produced objects

# creating a Venn diagram visualization (sh1)
library("ggplot2")
library("ggVennDiagram")
library("egg")
library("grid")

x <- list(
    DESeq2 = DE.genes.sh1.deseq2$gene,
    EdgeR = DE.genes.sh1.edger$gene,
    Limma = DE.genes.sh1.limma$gene
)

a <- ggVennDiagram(x[c(1,2)]) +
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange")

b <- ggVennDiagram(x[c(1,3)]) +
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange")

c <- ggVennDiagram(x[c(2,3)])+
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange")

d <- ggVennDiagram(x) +
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
define_region <- function(row, col){
    viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(a, vp = define_region(row = 1, col = 1))
print(b, vp = define_region(row = 2, col = 1))
print(c, vp = define_region(row = 3, col = 1))
print(d, vp = define_region(row = 1:3, col = 2))

rm(a,b,c,d,x, define_region)


# creating a Venn diagram visualization (sh2)

x <- list(
    Deseq2 = DE.genes.sh2.deseq2$gene,
    EdgeR = DE.genes.sh2.edger$gene,
    Limma = DE.genes.sh2.limma$gene
)

a <- ggVennDiagram(x[c(1,2)]) +
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange")

b <- ggVennDiagram(x[c(1,3)]) +
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange")

c <- ggVennDiagram(x[c(2,3)])+
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange")

d <- ggVennDiagram(x) +
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
define_region <- function(row, col){
    viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(a, vp = define_region(row = 1, col = 1))
print(b, vp = define_region(row = 2, col = 1))
print(c, vp = define_region(row = 3, col = 1))
print(d, vp = define_region(row = 1:3, col = 2))

rm(a,b,c,d,x, define_region)


# number of genes DE in both experiments

x <- list(
    "shNT vs sh1" = DE.genes.sh1.deseq2$gene,
    "shNT vs sh2" = DE.genes.sh2.deseq2$gene
)

a <- ggVennDiagram(x[c(1,2)]) +
        ggplot2::scale_fill_gradient(low="lightyellow",high = "orange") + 
        theme(legend.position = "none")

x <- list(
    "shNT vs sh1" = DE.genes.sh1.edger$gene,
    "shNT vs sh2" = DE.genes.sh2.edger$gene
)

b <- ggVennDiagram(x[c(1,2)]) +
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange") + 
    theme(legend.position = "none")

x <- list(
    "shNT vs sh1" = DE.genes.sh1.limma$gene,
    "shNT vs sh2" = DE.genes.sh2.limma$gene
)

c <- ggVennDiagram(x[c(1,2)]) +
    ggplot2::scale_fill_gradient(low="lightyellow",high = "orange") + 
    theme(legend.position = "none")

ggarrange(a, b, c, ncol = 3, nrow = 1, labels = c("DESeq2","EdgeR","Limma"))

rm(a,b,c,x)

# -------------------------------

# biotype of the commonly reported DE genes


# common genes in deseq2 and edger (sh1)
common.deseq2_edger.sh1 <- DE.genes.sh1.deseq2 %>%
                            inner_join(DE.genes.sh1.edger, by = c("gene" = "gene")) %>%
                            select(gene)

# common genes in deseq2 and limma (sh1)
common.deseq2_limma.sh1 <- DE.genes.sh1.deseq2 %>%
                            inner_join(DE.genes.sh1.limma, by = c("gene" = "gene")) %>%
                            select(gene)

# common genes in edger and limma (sh1)
common.edger_limma.sh1 <- DE.genes.sh1.edger %>%
                            inner_join(DE.genes.sh1.limma, by = c("gene" = "gene")) %>%
                            select(gene)

# common genes in all 3 (sh1)
common.sh1 <- common.deseq2_edger.sh1 %>% 
                inner_join(common.deseq2_limma.sh1) %>% 
                inner_join(common.edger_limma.sh1)

rm(common.deseq2_edger.sh1)
rm(common.deseq2_limma.sh1)
rm(common.edger_limma.sh1)

# common genes in deseq2 and edger (sh2)
common.deseq2_edger.sh2 <- DE.genes.sh2.deseq2 %>%
                            inner_join(DE.genes.sh2.edger, by = c("gene" = "gene")) %>%
                            select(gene)

# common genes in deseq2 and limma (sh2)
common.deseq2_limma.sh2 <- DE.genes.sh2.deseq2 %>%
                            inner_join(DE.genes.sh2.limma, by = c("gene" = "gene")) %>%
                            select(gene)

# common genes in edger and limma (sh2)
common.edger_limma.sh2 <- DE.genes.sh2.edger %>%
                            inner_join(DE.genes.sh2.limma, by = c("gene" = "gene")) %>%
                            select(gene)

# common genes in all 3 (sh2)
common.sh2 <- common.deseq2_edger.sh2 %>% 
                inner_join(common.deseq2_limma.sh2) %>% 
                inner_join(common.edger_limma.sh2)

rm(common.deseq2_edger.sh2)
rm(common.deseq2_limma.sh2)
rm(common.edger_limma.sh2)

# assigning a biotype to the common DE genes and counting how
# many are protein_coding and how many are lncRNA

library(readr)
library(stringr)
library(dplyr)

setwd(paste0(basedir, "4-4"))

gene_name_type <- read_delim("gtf_only_genes-name_type.txt", 
                             ";", 
                             escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

gene_name_type <- gene_name_type %>%
                    mutate(gene = str_remove_all(X1, "gene_id")) %>%
                    mutate(gene = str_remove_all(gene, "\"")) %>%
                    mutate(gene = str_remove_all(gene, " ")) %>%
                    mutate(biotype = str_remove_all(X2, "gene_type")) %>%
                    mutate(biotype = str_remove_all(biotype, "\"")) %>%
                    mutate(biotype = str_remove_all(biotype, " ")) %>%
                    select(gene, biotype)

# temp.sh1(sh2) contains the common DE genes in Deseq2, EdgeR and Limma
# and their biotype

temp.sh1 <- common.sh1 %>%
    inner_join(gene_name_type, by = c("gene" = "gene") )

temp.sh2 <- common.sh2 %>%
    inner_join(gene_name_type, by = c("gene" = "gene") )

# number of genes of every biotype in the Gencode V35 annotation
annotation_biotype <- gene_name_type %>% 
                        group_by(biotype) %>% 
                        summarise(total = n()) %>% 
                        arrange(-total)

# common.shX.biotype contains, for every biotype, the number of commonly detected DE genes of that biotype
common.sh1.biotype <- temp.sh1 %>% 
                        group_by(biotype) %>% 
                        summarise(n=n()) %>%
                        inner_join(annotation_biotype, by = c("biotype" = "biotype")) %>%
                        arrange(-n) %>%
                        mutate('% of DE' = n/nrow(temp.sh1) * 100) %>%
                        mutate('% in annotation' = n/total * 100)

common.sh2.biotype <- temp.sh2 %>% 
                        group_by(biotype) %>% 
                        summarise(n=n()) %>%
                        inner_join(annotation_biotype, by = c("biotype" = "biotype")) %>%
                        arrange(-n) %>%
                        mutate('% of DE' = n/nrow(temp.sh2) * 100) %>%
                        mutate('% in annotation' = n/total * 100)

# export common DE genes
write.table(common.sh1, file = "common-DE-genes/DE_common_sh1.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

write.table(common.sh2, file = "common-DE-genes/DE_common_sh2.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

rm(common.sh1, common.sh2)
rm(gene_name_type, temp.sh1, temp.sh2, common.sh1.biotype, common.sh2.biotype)

# ---------------------------------------------

# genes lncRNA DE

# DE genes by Deseq2 in sh1 (sh2) vs shNT that are lncRNAs
DE.genes.sh1.deseq2 %>% 
    inner_join(gene_name_type, by = c("gene" = "gene")) %>% 
    filter(biotype == "lncRNA") %>%
    nrow()

DE.genes.sh2.deseq2 %>% 
    inner_join(gene_name_type, by = c("gene" = "gene")) %>% 
    filter(biotype == "lncRNA") %>%
    nrow()

# DE genes by EdgeR in sh1 (sh2) vs shNT that are lncRNAs
DE.genes.sh1.edger %>% 
    inner_join(gene_name_type, by = c("gene" = "gene")) %>% 
    filter(biotype == "lncRNA") %>%
    nrow()

DE.genes.sh2.edger %>% 
    inner_join(gene_name_type, by = c("gene" = "gene")) %>% 
    filter(biotype == "lncRNA") %>%
    nrow()

# DE genes by Limma in sh1 (sh2) vs shNT that are lncRNAs
DE.genes.sh1.limma %>% 
    inner_join(gene_name_type, by = c("gene" = "gene")) %>% 
    filter(biotype == "lncRNA") %>%
    nrow()

DE.genes.sh2.limma %>% 
    inner_join(gene_name_type, by = c("gene" = "gene")) %>% 
    filter(biotype == "lncRNA") %>%
    nrow()

# ------------------------------------------


basedir <- "/mnt/hdd0/univ/thesis-support-material/data/"

setwd(paste0(basedir, "4-5"))

library(readr)
library(dplyr)
library(ggplot2)

n_h3k27ac <- 58201

data.frame(sites = rep(c("p53", "random 1", "random 2", "random 3") , 2),
           value = c(7256, 353, 376, 344)) %>%
    ggplot(mapping = aes(x = sites, y = value, fill = sites)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    geom_text(aes(label = value), vjust=-0.5, size=5, color = "black", position = position_dodge(0.9)) +
    geom_text(aes(label = paste0(format(round(value/n_h3k27ac*100, 1), nsmall = 1), "%")), 
              vjust=1.5, size=4, color = "white", position = position_dodge(0.9)) +
    theme_minimal() + 
    theme(legend.position = c(0.9,0.85)) +
    labs(title = "Number of H3K27ac sites containing a p53 binding site", x = "", y = "Sites")


# for each h3k27ac region on a p53 binding site: the coordinates of region and site and the name, logfc, distance of
# the closest DE gene

gene_symbol_Enseml_mapping <- read_delim("gene-symbol-Enseml-mapping.txt", 
                                         "\t", escape_double = FALSE, col_names = FALSE, 
                                         trim_ws = TRUE)

closest_DEgene_sh1 <- read_delim("question1/2-closest_DE_gene/closest-gene/closest_DEgene_sh1.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_h3k27ac_p53 <- read_delim("question1/2-closest_DE_gene/intersections-h3k27ac_p53/intersections.txt", 
                                        "\t", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE)

DESeq2_DE_sh1 <- read_delim("question1/aux/DESeq2_DE_sh1.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(intersections_h3k27ac_p53) <- c("chr_h3k27ac", "start_h3k27ac", "end_h3k27ac", "name_h3k27ac", 
                                         "chr_p53", "start_p53", "end_p53", "name_p53", "overlap")

colnames(closest_DEgene_sh1) <- c("chr_h3k27ac", "start_h3k27ac", "end_h3k27ac", "name_h3k27ac", 
                                  "chr_gene", "start_gene", "end_gene", "name_gene", "nd", "strand", "distance")

temp <- closest_DEgene_sh1 %>%
            select("name_h3k27ac", "chr_gene", "start_gene", "end_gene", "name_gene", "distance")

table1 <- intersections_h3k27ac_p53 %>%
    inner_join(temp, by = c("name_h3k27ac" = "name_h3k27ac")) %>%
    inner_join(DESeq2_DE_sh1, by = c("name_gene" = "gene")) %>%
    select("chr_h3k27ac", "start_h3k27ac", "end_h3k27ac", "start_p53", "end_p53", "name_gene", "log2FoldChange", "distance") %>%
    mutate(name_gene = substr(name_gene, 0, 15)) %>%
    inner_join(gene_symbol_Enseml_mapping, by = c("name_gene" = "X1")) %>%
    rename("gene name" = "X2") %>%
    select("chr_h3k27ac", "start_h3k27ac", "end_h3k27ac", "start_p53", "end_p53", "gene name", "log2FoldChange", "distance")

rm(temp)

# LINC01605: 
#chr8:37278859-37411701
b1 <- 37278859
b2 <- 37411701
table1 %>% filter(chr_h3k27ac == "chr8") %>% filter(start_h3k27ac >= b1 & end_h3k27ac <= b2)

# --------------------

# marking which of the genes within 1MB from the h3k27ac regions on a p53 site are DE

genes_within_1MB <- read_table2("question1/3-genes_within_1MB_from_1/a-genes_within_1MB/genes_within_1MB.txt", 
                                col_names = FALSE)

# importing the DE genes in sh1 vs shNT and in sh2 vs shNT, then keeping only those DE in both experiments
DE_common_sh1 <- read_delim("../4-4/sites-on-DE-genes/genes/DE_common_sh1_namepos.bed", 
                            "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

DE_common_sh2 <- read_delim("../4-4/sites-on-DE-genes/genes/DE_common_sh2_namepos.bed", 
                            "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

DE_sh1sh2 <- (DE_common_sh1 %>% inner_join(DE_common_sh2))


genes_within_1MB <- genes_within_1MB %>% 
                    mutate(isDE = X8 %in% DE_sh1sh2$X4)

colnames(genes_within_1MB) <- c("chr_h", "start_h", "end_h", "region_name", 
                                "chr_g", "start_g", "end_g", "gene_name", 
                                "nd", "gene_strand", "distance", "isGeneDE" )

a <- genes_within_1MB %>% 
        dplyr::group_by(chr_h, start_h, end_h, region_name) %>% 
        dplyr::summarise(nGenes = n(), nGenesDE = sum(isGeneDE[isGeneDE])) %>%
        dplyr::mutate(fractionOfDEPercent = 100 * nGenesDE / nGenes) %>%
        dplyr::filter(nGenes >= 10) %>%
        dplyr::arrange(-fractionOfDEPercent)

# regions on LINC01605's locus
a %>% filter(chr_h == "chr8", start_h %in% c(37371529,37373258,37397879,37401725))

# -- graph

nodes <-  a %>%
            select(region_name, chr_h, start_h, end_h, fractionOfDEPercent) %>%
            mutate(ycoord = as.integer(substr(chr_h, 4, 10)))

DE_sh1sh2 <- DE_sh1sh2 %>%
                dplyr::mutate(ycoord = as.integer(substr(X1, 4, 99)))

nodes %>%
    filter(fractionOfDEPercent > 3.0) %>%
    ggplot() +
        geom_hline(mapping = aes(yintercept = ycoord), color = "gray80") +
        geom_point(mapping = aes(x = start_h, y = ycoord, size = fractionOfDEPercent, color = as.factor(ycoord))) + 
        #geom_point(mapping = aes(x = start_h, y = ycoord,size = fractionOfDEPercent), colour = "black", shape = 1) +
        scale_color_discrete(guide = F) + 
        scale_size_continuous(range = c(2,12), name = "% of DE genes within 1MB") +
        geom_point(data = DE_sh1sh2, mapping = aes(x = X2, y = ycoord), color = "black", shape = 4) +
        theme_classic() +
        theme(legend.position = c(0.9, 0.9) ) +
        labs(title = "H3K27ac regions having DE genes within 1M base pairs") + 
        xlab("Position within the chromosome") + 
        ylab("Chromosome number")

# --------------------

# marking which of the genes within 50KB from the h3k27ac regions on a p53 site are DE

genes_within_50KB <- read_table2("question1/4-genes_within_50KB_from_1/genes_within_50KB/genes_within_50KB.txt", 
                                  col_names = FALSE)

genes_within_50KB <- genes_within_50KB %>% 
                     mutate(isDE = X8 %in% DE_sh1sh2$X4)

colnames(genes_within_50KB) <- c("chr_h", "start_h", "end_h", "region_name", 
                                "chr_g", "start_g", "end_g", "gene_name", 
                                "nd", "gene_strand", "distance", "isGeneDE")


#write.table(DE_sh1sh2$X4, file = "../4-6/all_vs_50kb/input/raw/allDEgenes.txt",
#            quote = F, row.names = F, col.names = F, sep = "\t")

#DE50KB <- unique((genes_within_50KB %>% dplyr::filter(isGeneDE == T) %>% dplyr::select(gene_name))$gene_name)

#write.table(DE50KB, file = "../4-6/all_vs_50kb/input/raw/DE50KB.txt",
#            quote = F, row.names = F, col.names = F, sep = "\t")

# -----------------------------------------------
# considering now the DE genes close to a DE h3k27ac region with p53 binding

closest_DEgene_sh1 <- read_delim("question2/3-closest_DE_gene/closest_DEgene_sh1.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

colnames(closest_DEgene_sh1) <- c("chr_h3k27ac", "start_h3k27ac", "end_h3k27ac", "name_h3k27ac", 
                                  "chr_gene", "start_gene", "end_gene", "name_gene", "nd", "strand", "distance")

temp <- closest_DEgene_sh1 %>%
    select("name_h3k27ac", "chr_gene", "start_gene", "end_gene", "name_gene", "distance")

table1 <- intersections_h3k27ac_p53 %>%
    inner_join(temp, by = c("name_h3k27ac" = "name_h3k27ac")) %>%
    inner_join(DESeq2_DE_sh1, by = c("name_gene" = "gene")) %>%
    select("chr_h3k27ac", "start_h3k27ac", "end_h3k27ac", "start_p53", "end_p53", "name_gene", "log2FoldChange", "distance") %>%
    mutate(name_gene = substr(name_gene, 0, 15)) %>%
    inner_join(gene_symbol_Enseml_mapping, by = c("name_gene" = "X1")) %>%
    rename("gene name" = "X2") %>%
    select("chr_h3k27ac", "start_h3k27ac", "end_h3k27ac", "start_p53", "end_p53", "gene name", "log2FoldChange", "distance")

table1 %>% filter(chr_h3k27ac == "chr8") %>% filter(start_h3k27ac >= b1 & end_h3k27ac <= b2)

# marking which of the genes within 1MB from the h3k27ac regions on a p53 site are DE

genes_within_1MB_sh1 <- read_table2("question2/4-genes_within_1MB_from_2/genes_within_1MB_sh1.txt", 
                                col_names = FALSE)

genes_within_1MB_sh2 <- read_table2("question2/4-genes_within_1MB_from_2/genes_within_1MB_sh2.txt", 
                                    col_names = FALSE)

genes_within_1MB_sh1 <- genes_within_1MB_sh1 %>% 
    mutate(isDE = X8 %in% DE_common_sh1$X4)

genes_within_1MB_sh2 <- genes_within_1MB_sh2 %>% 
    mutate(isDE = X8 %in% DE_common_sh2$X4)

colnames(genes_within_1MB_sh1) <- c("chr_h", "start_h", "end_h", "region_name", 
                                    "chr_g", "start_g", "end_g", "gene_name", 
                                    "nd", "gene_strand", "distance", "isGeneDE" )

colnames(genes_within_1MB_sh2) <- c("chr_h", "start_h", "end_h", "region_name", 
                                    "chr_g", "start_g", "end_g", "gene_name", 
                                    "nd", "gene_strand", "distance", "isGeneDE" )

a <- genes_within_1MB_sh1 %>% 
    dplyr::group_by(chr_h, start_h, end_h, region_name) %>% 
    dplyr::summarise(nGenes = n(), nGenesDE = sum(isGeneDE[isGeneDE])) %>%
    dplyr::mutate(fractionOfDEPercent = 100 * nGenesDE / nGenes) %>%
    dplyr::filter(nGenes >= 10) %>%
    dplyr::arrange(-fractionOfDEPercent)


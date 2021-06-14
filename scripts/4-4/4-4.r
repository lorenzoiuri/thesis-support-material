
basedir <- "/mnt/hdd0/univ/thesis-support-material/data/"

# number of reads mapped to the TP53 and LINC01605 genes

library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd(paste0(basedir, "4-4/"))

counts_full <- read_delim("featureCounts-output/counts.txt", "\t", 
                          escape_double = FALSE, trim_ws = TRUE, skip = 1)

counts_simple_full <- as.matrix(counts_full[,c(7,8,9,10,11,12)])

# renaming the rows and columns

# note: bc*-1 are the shNT samples
#       bc*-2 are the sh1 samples
#       bc*-3 are the sh2 samples

colnames(counts_simple_full) <- c("bc1-1", "bc1-2", "bc1-3", "bc2-1", "bc2-2", "bc2-3")
rownames(counts_simple_full) <- counts_full$Geneid

counts_simple_full['ENSG00000141510.18_9',] # TP53
counts_simple_full['ENSG00000253161.5_9',] # LINC01605

# -----------------------

a <- data.frame(exp = c(rep("sh1 vs shNT" , 2) , rep("sh2 vs shNT" , 2)),
           annotation = rep(c("only lncRNA", "full") , 2),
           value = c(21, 34, 25, 30)) %>%
        ggplot(mapping = aes(x = exp, y = value, fill = annotation)) +
            geom_bar(stat = "identity", position = "dodge") +
            scale_fill_manual(values = c("gold3", "darkslateblue")) +
            geom_text(aes(label = value), vjust=-0.3, size=4, color = "black", position = position_dodge(0.9))+
            theme_minimal() + 
            theme(legend.position = "none") +
            labs(title = "DESeq2", x = "", y = "lncRNA DE genes")

b <- data.frame(exp = c(rep("sh1 vs shNT" , 2) , rep("sh2 vs shNT" , 2)),
                annotation = rep(c("only lncRNA", "full") , 2),
                value = c(24, 37, 28, 62)) %>%
    ggplot(mapping = aes(x = exp, y = value, fill = annotation)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("gold3", "darkslateblue")) +
    geom_text(aes(label = value), vjust=-0.3, size=4, color = "black", position = position_dodge(0.9))+
    theme_minimal() + 
    theme(legend.position = "none") +
    labs(title = "EdgeR", x = "", y = "lncRNA DE genes")

c <- data.frame(exp = c(rep("sh1 vs shNT" , 2) , rep("sh2 vs shNT" , 2)),
                annotation = rep(c("only lncRNA", "full") , 2),
                value = c(7, 73, 3, 76)) %>%
    ggplot(mapping = aes(x = exp, y = value, fill = annotation)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("gold3", "darkslateblue")) +
    geom_text(aes(label = value), vjust=-0.3, size=4, color = "black", position = position_dodge(0.9))+
    theme_minimal() + 
    theme(legend.position = "right") +
    labs(title = "Limma", x = "", y = "lncRNA DE genes")

ggarrange(a, b, c, ncol = 3, nrow = 1)
rm(a,b,c)

# -------------------------------------------------------

# Interactions between differentially expressed genes and mutant p53}

# - p53 sites

intersections_p53sites_DEsh1 <- read_delim("sites-on-DE-genes/p53-sites/intersections_p53sites_DEsh1.txt", 
                                           "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
intersections_p53sites_DEsh2 <- read_delim("sites-on-DE-genes/p53-sites/intersections_p53sites_DEsh2.txt", 
                                           "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# p53 sites intersecting at least one DE gene (sh1)
intersections_p53sites_DEsh1 %>% dplyr::group_by(X4) %>% dplyr::count() %>% nrow()

# p53 sites intersecting at least one DE gene (sh2)
intersections_p53sites_DEsh2 %>% dplyr::group_by(X4) %>% dplyr::count() %>% nrow()

# - random sites sample 1

intersections_rs1_DEsh1 <- read_delim("sites-on-DE-genes/random-sites/sample1/intersections_rs1_DEsh1.txt", 
                                      "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_rs1_DEsh2 <- read_delim("sites-on-DE-genes/random-sites/sample1/intersections_rs1_DEsh2.txt", 
                                      "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_rs1_DEsh1 %>% dplyr::group_by(X4) %>% dplyr::count() %>% nrow()
intersections_rs1_DEsh2 %>% dplyr::group_by(X4) %>% dplyr::count() %>% nrow()

# - random sites sample 2

intersections_rs2_DEsh1 <- read_delim("sites-on-DE-genes/random-sites/sample2/intersections_rs2_DEsh1.txt", 
                                      "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_rs2_DEsh2 <- read_delim("sites-on-DE-genes/random-sites/sample2/intersections_rs2_DEsh2.txt", 
                                      "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_rs2_DEsh1 %>% dplyr::group_by(X4) %>% dplyr::count() %>% nrow()
intersections_rs2_DEsh2 %>% dplyr::group_by(X4) %>% dplyr::count() %>% nrow()

# - random sites sample 3

intersections_rs3_DEsh1 <- read_delim("sites-on-DE-genes/random-sites/sample3/intersections_rs3_DEsh1.txt", 
                                      "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_rs3_DEsh2 <- read_delim("sites-on-DE-genes/random-sites/sample3/intersections_rs3_DEsh2.txt", 
                                      "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_rs3_DEsh1 %>% dplyr::group_by(X4) %>% dplyr::count() %>% nrow()
intersections_rs3_DEsh2 %>% dplyr::group_by(X4) %>% dplyr::count() %>% nrow()

# - plots

# p53 sites on DE genes (sh1): 607
# random sites on DE genes (sh1): 265 282 268 - 2.25x

# p53 sites on DE genes (sh2): 685
# random sites on DE genes (sh2): 263 268 245 - 4.4x
# 1200x550
data.frame(exp = c(rep("sh1 vs shNT" , 4) , rep("sh2 vs shNT" , 4)),
           sites = rep(c("p53", "random 1", "random 2", "random 3") , 2),
           value = c(607, 265, 282, 268, 685, 263, 268, 245)) %>%
    ggplot(mapping = aes(x = exp, y = value, fill = sites)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    geom_text(aes(label = value), vjust=1.6, size=5, color = "white", position = position_dodge(0.9)) +
    geom_text(aes(label = paste0(format(round(value/11028*100, 1), nsmall = 1), "%")), 
                  vjust=4, size=4, color = "white", position = position_dodge(0.9)) +
    theme_minimal() + 
    theme(legend.position = c(0.9, 0.85)) + 
    #theme(legend.position = "right") +
    labs(title = "Number of p53 and random sites (out of 11028) located on DE genes", x = "", y = "sites")


# -- sites not on DE genes but close to one --

closest_p53_sh1 <- read_delim("sites-on-DE-genes/p53-sites/sites-not-on-DE-genes/closest_DEgene_sh1.txt", 
                              "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
closest_p53_sh2 <- read_delim("sites-on-DE-genes/p53-sites/sites-not-on-DE-genes/closest_DEgene_sh2.txt", 
                              "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

closest_p53_sh1 <- closest_p53_sh1 %>% filter(X8 != ".")
closest_p53_sh2 <- closest_p53_sh2 %>% filter(X8 != ".")

closest_rs1_sh1 <- read_delim("sites-on-DE-genes/random-sites/sample1/sites-not-on-DE-genes/closest_DEgene_sh1.txt", 
                              "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
closest_rs1_sh2 <- read_delim("sites-on-DE-genes/random-sites/sample1/sites-not-on-DE-genes/closest_DEgene_sh2.txt", 
                              "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

closest_rs1_sh1 <- closest_rs1_sh1 %>% filter(X8 != ".")
closest_rs1_sh2 <- closest_rs1_sh2 %>% filter(X8 != ".")

closest_rs2_sh1 <- read_delim("sites-on-DE-genes/random-sites/sample2/sites-not-on-DE-genes/closest_DEgene_sh1.txt", 
                              "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
closest_rs2_sh2 <- read_delim("sites-on-DE-genes/random-sites/sample2/sites-not-on-DE-genes/closest_DEgene_sh2.txt", 
                              "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

closest_rs2_sh1 <- closest_rs2_sh1 %>% filter(X8 != ".")
closest_rs2_sh2 <- closest_rs2_sh2 %>% filter(X8 != ".")

closest_rs3_sh1 <- read_delim("sites-on-DE-genes/random-sites/sample3/sites-not-on-DE-genes/closest_DEgene_sh1.txt", 
                              "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
closest_rs3_sh2 <- read_delim("sites-on-DE-genes/random-sites/sample3/sites-not-on-DE-genes/closest_DEgene_sh2.txt", 
                              "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

closest_rs3_sh1 <- closest_rs3_sh1 %>% filter(X8 != ".")
closest_rs3_sh2 <- closest_rs3_sh2 %>% filter(X8 != ".")

# number of not sites not on DE genes but within threshold to one (in percent)

threshold <- 2000

target <- closest_p53_sh1
a1 <- target %>% filter((X11 >= 0 & X11 <= threshold) | (X11 < 0 & X11 >= -1*threshold)) %>%
    nrow()/nrow(target) * 100

target <- closest_p53_sh2
a2 <- target %>% filter((X11 >= 0 & X11 <= threshold) | (X11 < 0 & X11 >= -1*threshold)) %>%
    nrow()/nrow(target) * 100

target <- closest_rs1_sh1
b1 <- target %>% filter((X11 >= 0 & X11 <= threshold) | (X11 < 0 & X11 >= -1*threshold)) %>%
    nrow()/nrow(target) * 100

target <- closest_rs1_sh2
b2 <- target %>% filter((X11 >= 0 & X11 <= threshold) | (X11 < 0 & X11 >= -1*threshold)) %>%
    nrow()/nrow(target) * 100

target <- closest_rs2_sh1
c1 <- target %>% filter((X11 >= 0 & X11 <= threshold) | (X11 < 0 & X11 >= -1*threshold)) %>%
    nrow()/nrow(target) * 100

target <- closest_rs2_sh2
c2 <- target %>% filter((X11 >= 0 & X11 <= threshold) | (X11 < 0 & X11 >= -1*threshold)) %>%
    nrow()/nrow(target) * 100

target <- closest_rs3_sh1
d1 <- target %>% filter((X11 >= 0 & X11 <= threshold) | (X11 < 0 & X11 >= -1*threshold)) %>%
    nrow()/nrow(target) * 100

target <- closest_rs3_sh2
d2 <- target %>% filter((X11 >= 0 & X11 <= threshold) | (X11 < 0 & X11 >= -1*threshold)) %>%
    nrow()/nrow(target) * 100

# 1200x550
data.frame(exp = c(rep("sh1 vs shNT" , 4) , rep("sh2 vs shNT" , 4)),
           sites = rep(c("p53", "random 1", "random 2", "random 3") , 2),
           value = c(a1, b1, c1, d1, a2, b2, c2, d2)) %>%
    ggplot(mapping = aes(x = exp, y = value, fill = sites)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    geom_text(aes(label = paste0(format(round(value, 2), nsmall = 1), "%")), 
              vjust=-0.5, size=4, color = "black", position = position_dodge(0.9)) +
    theme_minimal() + 
    theme(legend.position = c(0.9, 0.85)) + 
    #theme(legend.position = "right") +
    labs(title = "Percentage of sites not on DE genes and within 2000 bp from one", 
         x = "", y = "percentage")

# distance distribution

distances <- data.frame(
    value = c(closest_p53_sh1$X11, closest_rs1_sh1$X11, closest_rs2_sh1$X11, closest_rs3_sh1$X11),
    sites = c(rep("p53", length(closest_p53_sh1$X11)), 
              rep("random 1", length(closest_rs1_sh1$X11)), 
              rep("random 2", length(closest_rs2_sh1$X11)), 
              rep("random 3", length(closest_rs3_sh1$X11))) ) %>%
    mutate(value = if_else(value < 0, -1*value, value))

distances %>% filter(sites == "p53") %>% summary() # median 653k
distances %>% filter(sites == "random 1") %>% summary() # median 1685k

a <- distances %>%
    ggplot(mapping = aes(x = sites, y = value, fill = sites)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    theme_minimal() + 
    theme(legend.position = "none") +
    labs(title = "Distance distribution to the closest DE gene", 
         x = "", y = "base pairs")

b <- distances %>% filter(value < 5000000) %>%
    ggplot(mapping = aes(x = sites, y = value, fill = sites)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    theme_minimal() + 
    theme(legend.position = "right") +
    labs(title = "Distance distribution to the closest DE gene, \n values limited to 5M base pairs", 
         x = "", y = "base pairs")

# 1200x550
ggarrange(a, b, ncol = 2, nrow = 1)
rm(a,b)

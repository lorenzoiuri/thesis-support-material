
basedir <- "/mnt/hdd0/univ/thesis-support-material/data/"

# producing plot for first figure

linc <- data.frame(experiment = factor(c(rep("WT",4), rep("KO",4)), levels = c("WT", "KO")),
                   sample = factor(rep(c("sample 1", "sample 2", "sample 3", "sample 4"), 2)),
                   value = c(33,167,202,138,6,9,0,7))

tp53 <- data.frame(experiment = factor(c(rep("shNT",2), rep("sh1",2), rep("sh2",2)), levels = c("shNT", "sh1", "sh2")),
                   sample = factor(rep(c("sample 1", "sample 2"), 3)),
                   value = c(588,647,137,182,245,215))


library(ggplot2)
library(dplyr)

library(RColorBrewer)

a <- linc %>%
    ggplot(mapping = aes(fill = sample, y = value, x = experiment)) + 
    theme_classic() +
    geom_bar(position = "dodge", stat = "identity") + 
    scale_fill_brewer(palette = "Paired") +
    labs(title = "LINC01605 expression levels upon LINC01605 knockout") +
    theme(legend.position = "top",
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line( size=.05, color="black"),
          panel.grid.minor.y = element_line( size=.05, color="black"),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.title.y = element_blank()
    )

b <- tp53 %>%
    ggplot(mapping = aes(fill = sample, y = value, x = experiment)) + 
    theme_classic() +
    geom_bar(position = "dodge", stat = "identity") + 
    scale_fill_brewer(palette = "Paired") +
    labs(title = "TP53 expression levels upon mutTP53 silencing (shRNA)") +
    theme(legend.position = "top",
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line( size=.05, color="black"),
          panel.grid.minor.y = element_line( size=.05, color="black"),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.title.y = element_blank()
    )

library(ggpubr)

# 1000x500
ggarrange(b, a, ncol = 2, nrow = 1)
rm(a,b)
rm(tp53, linc)

# -----------------------------------

# counting DE genes common in the TP53 silencing and in the LINC01605 KO experiments

library(readr)
library(ggVennDiagram)
library(ggplot2)

setwd(paste0(basedir, "4-7/DE_genes"))

DE_DESeq2_TP53sil_common_sh12 <- read_csv("DE_DESeq2_TP53sil_common_sh12.txt", col_names = FALSE)$X1

DE_DESeq2_LINC01605ko <- read_csv("DE_DESeq2_LINC01605ko.txt", col_names = FALSE)$X1

x <- list("DE genes in TP53 silencing" = DE_DESeq2_TP53sil_common_sh12,
          "DE genes in LINC01605 knockout" = DE_DESeq2_LINC01605ko)

ggVennDiagram(x) + ggplot2::scale_fill_gradient(low="lightyellow", high = "orange")

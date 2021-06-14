
basedir <- "/mnt/hdd0/univ/thesis-support-material/data/4-7/"

plotDots <- function(data,
                     plot_title = "",
                     subtitle = "") {
    data %>%
        ggplot(mapping=aes(x=ES, y=description, color=logpval, size=SIZE)) +
        geom_point() +
        scale_color_continuous(low = "green", high = "red") +
        theme_bw() +
        labs(x = "Enrichment score (ES)", y = "", 
             color = "-log10(pvalue)",
             size = "Gene set size",
             title = plot_title,
             subtitle = subtitle) +
        theme(text = element_text(size=14, family = "Helvetica")) +
        guides(size = guide_legend(order = 1),
               color_continuos = guide_legend(order = 2))
}

library(readr)

setwd(basedir)

# --- LINC KO ---

gsea_wt_lincKO <- read_delim("GSEA_KEGG/LINC01605KO/gsea_report_for_WT_1622032517847.tsv", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

gsea_ko_lincKO <- read_delim("GSEA_KEGG/LINC01605KO/gsea_report_for_KO_1622032517847.tsv", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)


library(dplyr)
library(stringr)

temp1 <- gsea_ko_lincKO %>% filter(`NOM p-val` <= 0.01)
temp2 <- gsea_wt_lincKO %>% filter(`NOM p-val` <= 0.01)

linc <- temp1 %>% dplyr::union(temp2) %>%
    mutate(description = stringr::str_to_sentence(gsub('_', ' ', gsub('kegg_', ' ', tolower(NAME))))) %>%
    mutate(logpval = -1 * log10(`NOM p-val`)) %>%
    select(description, SIZE, ES, logpval) %>%
    mutate(test = max(logpval))

library(ggplot2)

plot_title <- "Enriched KEGG pathways in LINC01605 KO experiment"
subtitle <- "ES > 0: pathway enriched in WT samples \nES < 0: pathway enriched in KO samples"

# ES > 0: pathway enriched in WT samples
# ES < 0: pathway enriched in KO samples
# 880x780
plotDots(linc, plot_title, subtitle)

# --- TP53 sh1 ---

gsea_sh1_tp53sh1 <- read_delim("GSEA_KEGG/TP53/sh1/gsea_report_for_sh1_1622033140136.tsv", 
                                "\t", escape_double = FALSE)

gsea_shNT_tp53sh1 <- read_delim("GSEA_KEGG/TP53/sh1/gsea_report_for_shNT_1622033140136.tsv", 
                               "\t", escape_double = FALSE)

temp1 <- gsea_sh1_tp53sh1 %>% filter(`NOM p-val` <= 0.01)
temp2 <- gsea_shNT_tp53sh1 %>% filter(`NOM p-val` <= 0.01)

sh1 <- temp1 %>% dplyr::union(temp2) %>%
    mutate(description = stringr::str_to_sentence(gsub('_', ' ', gsub('kegg_', ' ', tolower(NAME))))) %>%
    mutate(logpval = -1 * log10(`NOM p-val`)) %>%
    select(description, SIZE, ES, logpval)

plot_title <- "Enriched KEGG pathways in TP53 sh1 experiment"
subtitle <- "ES > 0: pathway enriched in non-silenced samples \nES < 0: pathway enriched in silenced samples"

# ES > 0: pathway enriched in shNT samples
# ES < 0: pathway enriched in sh1 samples
# 880x780
plotDots(sh1, plot_title, subtitle)

# --- TP53 sh2 ---

gsea_sh2_tp53sh2 <- read_delim("GSEA_KEGG/TP53/sh2/gsea_report_for_sh2_1622033289948.tsv", 
                               "\t", escape_double = FALSE)

gsea_shNT_tp53sh2 <- read_delim("GSEA_KEGG/TP53/sh2/gsea_report_for_shNT_1622033289948.tsv", 
                                "\t", escape_double = FALSE)

temp1 <- gsea_sh2_tp53sh2 %>% filter(`NOM p-val` <= 0.01)
temp2 <- gsea_shNT_tp53sh2 %>% filter(`NOM p-val` <= 0.01)

sh2 <- temp1 %>% dplyr::union(temp2) %>%
    mutate(description = stringr::str_to_sentence(gsub('_', ' ', gsub('kegg_', ' ', tolower(NAME))))) %>%
    mutate(logpval = -1 * log10(`NOM p-val`)) %>%
    select(description, SIZE, ES, logpval)

plot_title <- "Enriched KEGG pathways in TP53 sh2 experiment"
subtitle <- "ES > 0: pathway enriched in non-silenced samples \nES < 0: pathway enriched in silenced samples"

# ES > 0: pathway enriched in shNT samples
# ES < 0: pathway enriched in sh1 samples
# 880x780
plotDots(sh2, plot_title, subtitle)  

# --- common pathways ---

# 11 common pathways
sh1sh2 <- sh1 %>% inner_join(sh2, by = c("description" = "description")) %>%
    filter(description != " Ribosome") %>%
    mutate(SIZE = SIZE.x,
           ES = (ES.x + ES.y) / 2,
           logpval = (logpval.x + logpval.y) / 2)

plot_title <- "Enriched KEGG pathways in both TP53 \nsh1 and sh2 experiments"
subtitle <- "ES > 0: pathway enriched in non-silenced samples \nES < 0: pathway enriched in silenced samples"
# 880x780
plotDots(sh1sh2, plot_title, subtitle)

# "empty"
sh1sh2 %>% inner_join(linc, by = c("description" = "description"))

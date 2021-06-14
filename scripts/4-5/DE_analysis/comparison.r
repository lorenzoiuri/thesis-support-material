

# histonic regions commonly considered DE (sh1)
common.sh1 <- DE.genes.sh1.deseq2$gene %>% 
                dplyr::intersect(DE.genes.sh1.edger$gene) %>%
                dplyr::intersect(DE.genes.sh1.limma$gene)

# histonic regions commonly considered DE (sh2)
common.sh2 <- DE.genes.sh2.deseq2$gene %>% 
                dplyr::intersect(DE.genes.sh2.edger$gene) %>%
                dplyr::intersect(DE.genes.sh2.limma$gene)

# importing the file that associated the name of a histonic region to its location
library(readr)
histone_bed <- read_delim("/mnt/hdd0/univ/thesis-support-material/data/4-5/h3k27ac_mda231.bed", 
                          "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


# getting the location for the DE histonic regions

histone_bed_sh1 <- histone_bed %>% 
                    dplyr::inner_join( data.frame(name = common.sh1), 
                                       by = c("X4" = "name"))

histone_bed_sh2 <- histone_bed %>% 
                    dplyr::inner_join( data.frame(name = common.sh2), 
                                       by = c("X4" = "name"))

# ---- venn diagrams ----

library("ggplot2")
library("ggVennDiagram")
library("egg")
library("grid")

printVennDiagrams <- function(x) {
    
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

}

# --sh1

x <- list(
    Deseq2 = DE.genes.sh1.deseq2$gene,
    EdgeR = DE.genes.sh1.edger$gene,
    Limma = DE.genes.sh1.limma$gene
)

printVennDiagrams(x)

# --sh2

x <- list(
    Deseq2 = DE.genes.sh2.deseq2$gene,
    EdgeR = DE.genes.sh2.edger$gene,
    Limma = DE.genes.sh2.limma$gene
)

printVennDiagrams(x)

# ---- exporting to bed ----

setwd("/mnt/hdd0/univ/thesis-support-material/data/4-5/question2/1-DE_histones/")
# 
write.table(histone_bed_sh1, file = "DE_histones_sh1.bed", 
            sep = "\t", quote = F, row.names = F, col.names = F)
# 
write.table(histone_bed_sh2, file = "DE_histones_sh2.bed", 
            sep = "\t", quote = F, row.names = F, col.names = F)




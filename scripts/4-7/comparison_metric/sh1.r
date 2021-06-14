
basedir <- "/mnt/hdd0/univ/thesis-support-material/"

setwd(basedir)

source("scripts/4-7/comparison_metric/function.r")
#source("function_mean.r")

# ----------

library(readr)
library(dplyr)

exp <- "sh1" # sh1 | linc_ko
pval <- "00001"

# --- BP ---

setwd(paste0(basedir, "data/4-7/functional_analysis_results/"))
GO <- "BP"

cp <- read_csv(paste0("clusterProfiler/", exp, "/pval", pval, "/revigo/", GO, ".csv"))
panther <- read_csv(paste0("Panther/", exp, "/pval", pval, "/revigo/", GO, ".csv"))
gsea <- read_csv(paste0("GSEA/", exp, "/pval", pval, "/revigo/", GO, ".csv"))

rs350 <- read_csv(paste0("random_samples_goterms/sample350t/revigo/", GO, ".csv"))
rs200 <- read_csv(paste0("random_samples_goterms/sample200t/revigo/", GO, ".csv"))

# the values in bracket are obtained using the "function_mean"

# the path must end with trailing slash!
files_path = paste0(basedir, "scripts/4-7/comparison_metric/")

# -- cp
computeSimilarity(cp, panther, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.94 (46+15/65) [0.91]
computeSimilarity(cp, gsea, GOType = GO, preComputedMatrix = "cp-gsea-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.79 (24+27/65) [0.77]
computeSimilarity(cp, rs350, GOType = GO, preComputedMatrix = "cp-rs350-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.54 (0+35/65) [0.50]
computeSimilarity(cp, rs200, GOType = GO, preComputedMatrix = "cp-rs200-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.52 (1+32/65) [0.49]

# -- panther
computeSimilarity(panther, cp, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.78 (48+52/129) [0.71]
computeSimilarity(panther, gsea, GOType = GO, preComputedMatrix = "panther-gsea-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F,
                  returnCommonClusters = T) # 0.85 (59+51/129) [0.83]
computeSimilarity(panther, rs350, GOType = GO, preComputedMatrix = "panther-rs350-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.55 (1+70/129) [0.51]
computeSimilarity(panther, rs200, GOType = GO, preComputedMatrix = "panther-rs200-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.53 (1+68/129) [0.50]

# -- gsea
computeSimilarity(gsea, cp, GOType = GO, preComputedMatrix = "gsea-cp-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.58 (24+74/169) [0.52]
computeSimilarity(gsea, panther, GOType = GO, preComputedMatrix = "gsea-panther-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.75 (59+68/169) [0.71]
computeSimilarity(gsea, rs350, GOType = GO, preComputedMatrix = "gsea-rs350-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.56 (4+90/169) [0.50]
computeSimilarity(gsea, rs200, GOType = GO, preComputedMatrix = "gsea-rs200-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.51 (3+84/169) [0.46]

# -- rs350
computeSimilarity(rs350, cp, GOType = GO, preComputedMatrix = "rs350-cp-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.23 (0+37/159) [0.21]
computeSimilarity(rs350, panther, GOType = GO, preComputedMatrix = "rs350-panther-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.28 (1+43/159) [0.25]
computeSimilarity(rs350, gsea, GOType = GO, preComputedMatrix = "rs350-gsea-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.30 (4+44/159) [0.28]
computeSimilarity(rs350, rs200, GOType = GO, preComputedMatrix = "rs350-rs200-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.24 (2+36/159) [0.23]

# -- rs200
computeSimilarity(rs200, cp, GOType = GO, preComputedMatrix = "rs200-cp-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.24 (1+24/104) [0.22]
computeSimilarity(rs200, panther, GOType = GO, preComputedMatrix = "rs200-panther-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.28 (1+28/104) [0.26]
computeSimilarity(rs200, gsea, GOType = GO, preComputedMatrix = "rs200-gsea-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.31 (3+29/104) [0.30]
computeSimilarity(rs200, rs350, GOType = GO, preComputedMatrix = "rs200-rs350-BP-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.28 (2+27/104) [0.25]


# ------- corrplot 

library(ggcorrplot)

r.cp = c(1, 0.94, 0.79, 0.54, 0.52)
r.panther = c(0.78, 1, 0.79, 0.54, 0.52)
r.gsea = c(0.58, 0.75, 1, 0.56, 0.51)
r.rs350 = c(0.23, 0.28, 0.30, 1, 0.24)
r.rs200 = c(0.24, 0.28, 0.31, 0.28, 1)

names <- c("cp", "panther", "gsea", "rs350", "rs200")

m <- matrix(data = c(r.cp, r.panther, r.gsea, r.rs350, r.rs200),
            nrow = length(names), ncol = length(names))

colnames(m) <- names
rownames(m) <- names

# 700x700
ggcorrplot(m, outline.col = "white",
           ggtheme = ggplot2::theme_minimal,
           colors = c("white", "burlywood3", "brown3"),
           lab = T,
           show.legend = T,
           title = "Values of the comparison metric \nTP53 sh1 silencing - \"Biological process\" GO terms")






# --- MF ---

setwd(paste0(basedir, "data/4-7/functional_analysis_results/"))
GO <- "MF"
 
cp <- read_csv(paste0("clusterProfiler/", exp, "/pval", pval, "/revigo/", GO, ".csv"))
panther <- read_csv(paste0("Panther/", exp, "/pval", pval, "/revigo/", GO, ".csv"))
gsea <- read_csv(paste0("GSEA/", exp, "/pval", pval, "/revigo/", GO, ".csv"))

rs350 <- read_csv(paste0("random_samples_goterms/sample350t/revigo/", GO, ".csv"))
rs200 <- read_csv(paste0("random_samples_goterms/sample200t/revigo/", GO, ".csv"))

# -- cp
computeSimilarity(cp, panther, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = T, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.81 (4+1.7/7)
computeSimilarity(cp, gsea, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.75 (3+2.3/7)
computeSimilarity(cp, rs350, GOType = GO, preComputedMatrix = "cp-rs350-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = T) # 0.10 (0+0.7/7)
computeSimilarity(cp, rs200, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.15 (0+1/7)

# -- panther
computeSimilarity(panther, cp, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = T, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.60 (4+4.3/14)
computeSimilarity(panther, gsea, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.82 (9+2.6/14)
computeSimilarity(panther, rs350, GOType = GO, preComputedMatrix = "panther-rs350-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.23 (0+3.2/14)
computeSimilarity(panther, rs200, GOType = GO, preComputedMatrix = "panther-rs200-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.15 (0+2/14)

# -- gsea
computeSimilarity(gsea, cp, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.37 (3+7.2/28)
computeSimilarity(gsea, panther, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.53 (9+5.8/28)
computeSimilarity(gsea, rs350, GOType = GO, preComputedMatrix = "gsea-rs350-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.24 (1+5.7/28)
computeSimilarity(gsea, rs200, GOType = GO, preComputedMatrix = "gsea-rs200-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.22 (0+6.1/28)

# -- rs350
computeSimilarity(rs350, cp, GOType = GO, preComputedMatrix = "rs350-cp-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.02 (0+1.2/57)
computeSimilarity(rs350, panther, GOType = GO, preComputedMatrix = "rs350-panther-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.07 (0+4/57) 
computeSimilarity(rs350, gsea, GOType = GO, preComputedMatrix = "rs350-gsea-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.11 (1+5.5/57)
computeSimilarity(rs350, rs200, GOType = GO, preComputedMatrix = "rs350-rs200-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.11 (1+5.2/57)


# -- rs200
computeSimilarity(rs200, cp, GOType = GO, preComputedMatrix = "rs200-cp-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.05 (0+1.7/35)
computeSimilarity(rs200, panther, GOType = GO, preComputedMatrix = "rs200-panther-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.09 (0+3/35)
computeSimilarity(rs200, gsea, GOType = GO, preComputedMatrix = "rs200-gsea-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.15 (0+5/35)
computeSimilarity(rs200, rs350, GOType = GO, preComputedMatrix = "rs200-rs350-MF-sh1", naiveMatrixCompute = T,
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.16 (1+4.5/35)

# ------- corrplot 

library(ggcorrplot)

r.cp = c(1, 0.81, 0.75, 0.10, 0.15)
r.panther = c(0.6, 1, 0.82, 0.23, 0.15)
r.gsea = c(0.37, 0.53, 1, 0.24, 0.22)
r.rs350 = c(0.02, 0.07, 0.11, 1, 0.11)
r.rs200 = c(0.05, 0.09, 0.15, 0.16, 1)

names <- c("cp", "panther", "gsea", "rs350", "rs200")

m <- matrix(data = c(r.cp, r.panther, r.gsea, r.rs350, r.rs200),
            nrow = length(names), ncol = length(names))

colnames(m) <- names
rownames(m) <- names

# 700x700
ggcorrplot(m, outline.col = "white",
           ggtheme = ggplot2::theme_minimal,
           colors = c("white", "burlywood3", "brown3"),
           lab = T,
           show.legend = T,
           title = "Values of the comparison metric \nTP53 sh1 silencing - \"Molecular function\" GO terms")

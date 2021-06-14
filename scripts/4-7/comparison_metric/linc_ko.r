
basedir <- "/mnt/hdd0/univ/thesis-support-material/"

setwd(basedir)
source("scripts/4-7/comparison_metric/function.r")

# ----------

library(readr)
library(dplyr)

exp <- "linc_ko" # sh1 | linc_ko
pval <- "0001"

# --- BP ---

setwd(paste0(basedir, "data/4-7/functional_analysis_results/"))
GO <- "BP"

cp <- read_csv(paste0("clusterProfiler/", exp, "/pval", pval, "/revigo/", GO, ".csv"))
panther <- read_csv(paste0("Panther/", exp, "/pval", pval, "/revigo/", GO, ".csv"))
gsea <- read_csv(paste0("GSEA/", exp, "/pval", pval, "/revigo/", GO, ".csv"))
rs50 <- read_csv(paste0("random_samples_goterms/sample50t/revigo/", GO, ".csv"))


# the path must end with trailing slash!
files_path = paste0(basedir, "scripts/4-7/comparison_metric/")

# -- cp
computeSimilarity(cp, panther, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.92 (6+2.3/9)
computeSimilarity(cp, gsea, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.93 (6+2.3/9)
computeSimilarity(cp, rs50, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.4 (0+3.6/9)

# -- panther
computeSimilarity(panther, cp, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.66 (6+3.2/14)
computeSimilarity(panther, gsea, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.81 (7+4.4/14)
computeSimilarity(panther, rs50, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.38 (0+5.3/14)

# -- gsea
computeSimilarity(gsea, cp, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.32 (6+8.55/45)
computeSimilarity(gsea, panther, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.37 (5+11.8/45)
computeSimilarity(gsea, rs50, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.4 (1+17/45)

# -- rs50
computeSimilarity(rs50, cp, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.17 (0+5.2/30)
computeSimilarity(rs50, panther, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.19 (0+5.7/30)
computeSimilarity(rs50, gsea, GOType = GO, preComputedMatrix = NULL, naiveMatrixCompute = F, 
                  path = files_path,
                  experiment = exp, exportMatrix = F) # 0.30 (1+7.9/30)

# ------ corrplot

library(ggcorrplot)

r.cp = c(1, 0.92, 0.93, 0.4)
r.panther = c(0.66, 1, 0.81, 0.38)
r.gsea = c(0.32, 0.37, 1, 0.4)
r.rs50 = c(0.17, 0.19, 0.30, 1)

names <- c("cp", "panther", "gsea", "rs50")

m <- matrix(data = c(r.cp, r.panther, r.gsea, r.rs50),
            nrow = length(names), ncol = length(names))

colnames(m) <- names
rownames(m) <- names

ggcorrplot(m, outline.col = "white",
           ggtheme = ggplot2::theme_minimal,
           colors = c("white", "burlywood3", "brown3"),
           lab = T,
           show.legend = T,
           title = "Values of the comparison metric \nLINC01605 knockout - \"Biological process\" GO terms")



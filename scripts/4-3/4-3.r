
basedir <- "/mnt/hdd0/univ/thesis-support-material/data"

setwd(paste0(basedir, "/4-3/"))

library(readr)
library(dplyr)
library(ggplot2)

# importing gencode_v35 and gencode_v24 annotation

gencode_v35 <- read_delim("gencode/gencode-v35.bed", 
                          "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

gencode_v24 <- read_delim("gencode/gencode-v24.bed", 
                          "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# --- 4.3.1 ---

# number of transcripts in annotation
gencode_v35 %>% nrow(.)
gencode_v24 %>% nrow(.)

# number of genes in annotation
gencode_v35 %>% group_by(X5) %>% count() %>% nrow(.)
gencode_v24 %>% group_by(X5) %>% count() %>% nrow(.)

# importing the files with the intersections between the p53 sites and the annotations

intersections_p53_v35 <- read_delim("bedtools-intersect/with-v35/intersections-p53-v35.bed", 
                                    "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_p53_v24 <- read_delim("bedtools-intersect/with-v24/intersections-p53-v24.bed", 
                                    "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# importing the files with the p53 sites not on genes

p53_no_intersections_v35 <- read_delim("bedtools-intersect/with-v35/p53-no-intersections.bed", 
                                       "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

p53_no_intersections_v24 <- read_delim("bedtools-intersect/with-v24/p53-no-intersections.bed", 
                                       "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# first approach

intersections_p53_v35 %>% group_by(X4,X9,X10) %>% count()
intersections_p53_v24 %>% group_by(X4,X9,X10) %>% count()

# second approach (see merge_transcripts.r)

merged_v35 <- read_delim("gencode/merged-annotation/merged-v35.bed", 
                         "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

merged_v24 <- read_delim("gencode/merged-annotation/merged-v24.bed", 
                         "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# intersections between the p53 sites and the merged transcripts

intersections_p53_merged_v35 <- read_delim("bedtools-intersect/with-merged-v35/intersections-p53-v35.bed", 
                                           "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_p53_merged_v24 <- read_delim("bedtools-intersect/with-merged-v24/intersections-p53-v24.bed", 
                                           "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# intersections, by biotype, v35
intersections_p53_merged_v35 %>% group_by(X9) %>% count() %>% arrange(-n)

# genes by biotype, in the merged v35 annotation
(merged_v35 %>% group_by(X4,X5) %>% count()) %>% group_by(X5) %>% count() %>% arrange(-n)

# transcripts of the gene snoU13
gencode_v35 %>% filter(X5 == "snoU13") %>% group_by(X1) %>% summarise(count = n(), avglen = mean(X3-X2))

# average length in the merged v35 annotation, of the genes by biotype
merged_v35 %>% group_by(X5) %>% summarise(count = n(), avglen = mean(X3-X2)) %>% arrange(-avglen)

# number of intersections between sites and genes, after snoRNA and misc_RNA removal
intersections_p53_merged_v35 %>% filter(X9 != "snoRNA", X9 != "misc_RNA")

intersections_p53_merged_v24 %>% filter(X9 != "snoRNA", X9 != "misc_RNA")


# --- 4.3.2 ---

# importing the estimated mutant p53 sites
p53_sites <- read_delim("sites/p53-sites.bed",
                        "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

p53_sites <- p53_sites %>% mutate(len = X3-X2)

# information about the lengths distribution
summary(p53_sites$len)

# 98th percentiles is <= 1000
percents <- seq(from = 0, to = 1, by = 0.01)
quantile(p53_sites$len, probs = percents)['98%'] 
rm(percents)

p53_sites_filtered <- p53_sites %>% filter(len <= 1000)

# boxplot length distribution
# export in 330x550
data.frame(sample = "p53", length = p53_sites$len) %>%
ggplot(aes(x = sample, y = length, fill = sample)) + 
    geom_boxplot() + 
    scale_fill_manual(values = c("coral3")) +
    theme_minimal() +
    theme(legend.position="none") + 
    labs(title = "Length distribution of the \n estimated p53 binding sites", x = "sample", y = "base pairs")

# using GAMLSS
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
ex_fit <- fitDist(p53_sites_filtered$len, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)
summary(ex_fit)

# random samples generation
random_sample_1 <- rGG(n = length(p53_sites_filtered$len), 
                         mu = exp(ex_fit$mu.coefficients), 
                         sigma = exp(ex_fit$sigma.coefficients), 
                         nu = ex_fit$nu.coefficients )

random_sample_2 <- rGG(n = length(p53_sites_filtered$len), 
                       mu = exp(ex_fit$mu.coefficients), 
                       sigma = exp(ex_fit$sigma.coefficients), 
                       nu = ex_fit$nu.coefficients )

random_sample_3 <- rGG(n = length(p53_sites_filtered$len), 
                       mu = exp(ex_fit$mu.coefficients), 
                       sigma = exp(ex_fit$sigma.coefficients), 
                       nu = ex_fit$nu.coefficients )

# removal of Inf values
random_sample_1 <- random_sample_1[random_sample_1 != Inf]
random_sample_2 <- random_sample_2[random_sample_2 != Inf]
random_sample_3 <- random_sample_3[random_sample_3 != Inf]

# comparison of the distribution in the samples, boxplots
df <- data.frame(sample = "p53", length = p53_sites_filtered$len)
df <- rbind(df, data.frame(sample = "random 1", length = random_sample_1))
df <- rbind(df, data.frame(sample = "random 2", length = random_sample_2))
df <- rbind(df, data.frame(sample = "random 3", length = random_sample_3))

# export at 530x550
ggplot(df, aes(x = sample, y = length, fill = sample)) + 
    geom_boxplot() + 
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    theme_minimal() +
    theme(legend.position="none") + 
    labs(title = "Length distribution of the estimated p53 binding sites \n compared to the randomly generated sites", 
         x = "sample", y = "base pairs")

# comparison of the distribution in the samples, density plot
# 630x550
ggplot(df, aes(x = length, fill = sample)) + 
    geom_density(alpha = 0.4, position = "identity") + 
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    theme_minimal() +
    theme(legend.position = c(0.9, 0.85)) + 
    labs(title = "Length distribution of the estimated p53 binding sites \n compared to the randomly generated sites",
         x = "base pairs")

rm(df)

# quantile - quantile plot
# 650x650
qqplot(p53_sites_filtered$len, random_sample_1, 
       main = "Comparison of the length distributions' quantiles: \n binding sites vs random sites from the fitted model",
       xlab = "Lengths of the estimated binding sites",
       ylab = "Lengths of the sites from the random samples")
abline(0, 1, col = 'red')


# --- 4.3.3 ---

# fraction of genes with at least one binding site
p1 <- (intersections_p53_v35 %>% group_by(X9,X10) %>% count()) %>%
    group_by(X10) %>% count() %>% arrange(-n)

p2 <- (gencode_v35 %>% group_by(X5,X6) %>% count()) %>%
    group_by(X6) %>% count() %>% arrange(-n)

p1 %>% inner_join(p2, by = c("X10" = "X6")) %>% mutate(ratio = n.x/n.y)

rm(p1)
rm(p2)

# comparison with the random samples

random_sample1 <- read_delim("sites/random-sites-samples/sample1.bed", 
                             "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
random_sample2 <- read_delim("sites/random-sites-samples/sample2.bed", 
                             "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
random_sample3 <- read_delim("sites/random-sites-samples/sample3.bed", 
                             "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# importing the results of bedtools intersect on the samples and the annotation

intersections_sample1 <- read_delim("bedtools-intersect/random-samples/sample1/bedtools_intersect.bed", 
                                    "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sites_not_on_genes_sample1 <- read_delim("bedtools-intersect/random-samples/sample1/random_no_intersections.bed", 
                                         "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_sample2 <- read_delim("bedtools-intersect/random-samples/sample2/bedtools_intersect.bed", 
                                    "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sites_not_on_genes_sample2 <- read_delim("bedtools-intersect/random-samples/sample2/random_no_intersections.bed", 
                                         "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

intersections_sample3 <- read_delim("bedtools-intersect/random-samples/sample3/bedtools_intersect.bed", 
                                    "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sites_not_on_genes_sample3 <- read_delim("bedtools-intersect/random-samples/sample3/random_no_intersections.bed", 
                                         "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# number of intersections
intersections_p53_v35 %>% group_by(X4,X9,X10) %>% count() %>% nrow() # 10368
intersections_sample1 %>% group_by(X4,X9,X10) %>% count() %>% nrow() # 6481
intersections_sample2 %>% group_by(X4,X9,X10) %>% count() %>% nrow() # 6554
intersections_sample3 %>% group_by(X4,X9,X10) %>% count() %>% nrow() # 6657

# number of genes intersected at least once
intersections_p53_v35 %>% group_by(X9,X10) %>% count() %>% nrow() # 8021
intersections_sample1 %>% group_by(X9,X10) %>% count() %>% nrow() # 4573
intersections_sample2 %>% group_by(X9,X10) %>% count() %>% nrow() # 4652
intersections_sample3 %>% group_by(X9,X10) %>% count() %>% nrow() # 4742

# number of sites not on genes
p53_no_intersections_v35 %>% nrow()   # 2899
sites_not_on_genes_sample1 %>% nrow() # 5214
sites_not_on_genes_sample2 %>% nrow() # 5124
sites_not_on_genes_sample3 %>% nrow() # 5084

# graphs
# 500x400
data.frame(sample = c("p53", "random 1", "random 2", "random 3"),
           value = c(10368,6481,6554,6657)  ) %>%
ggplot(mapping = aes(x = sample, y = value, fill = sample)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    geom_text(aes(label = value), vjust=1.6, size=4, color = "white")+
    theme_minimal() +
    labs(title = "Number of intersections between sites and genes", x = "sites sample", y = "intersections") + 
    theme(legend.position = "none")

# number of intersected genes by at least one region
# 500x400
data.frame(sample = c("p53", "random 1", "random 2", "random 3"),
                 value = c(8021,4573,4652,4742) ) %>%
ggplot(mapping = aes(x = sample, y = value, fill = sample)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    geom_text(aes(label = value), vjust=1.6, size=4, color = "white")+
    theme_minimal() +
    labs(title = "Number of genes intersected by at least one site", x = "sites sample", y = "genes") + 
    theme(legend.position = "none")

# number of sites not intersecting genes
# 500x400
data.frame(sample = c("p53", "random 1", "random 2", "random 3"),
                 value = c(2899,5214,5124,5084) ) %>%
ggplot(mapping = aes(x = sample, y = value, fill = sample)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    geom_text(aes(label = value), vjust=1.6, size=4, color = "white")+
    theme_minimal() +
    labs(title = "Number of sites not on genes (out of 11028)", x = "sites sample", y = "sites") + 
    theme(legend.position = "none")

# number of genes by biotype containing at least one p53 binding site or a site from the random samples

tot <- intersections_p53_v35 %>% group_by(X9,X10) %>% count() %>% nrow()
p1 <- (intersections_p53_v35 %>% group_by(X9,X10) %>% count()) %>% 
        group_by(X10) %>% count() %>% arrange(-n) %>%
        mutate(ratio = n/tot)

tot <- intersections_sample1 %>% group_by(X9,X10) %>% count() %>% nrow()
p2 <- (intersections_sample1 %>% group_by(X9,X10) %>% count()) %>% 
        group_by(X10) %>% count() %>% arrange(-n) %>%
        mutate(ratio = n/tot)

tot <- intersections_sample2 %>% group_by(X9,X10) %>% count() %>% nrow()
p3 <- (intersections_sample2 %>% group_by(X9,X10) %>% count()) %>% 
        group_by(X10) %>% count() %>% arrange(-n) %>%
        mutate(ratio = n/tot)

tot <- intersections_sample3 %>% group_by(X9,X10) %>% count() %>% nrow()
p4 <- (intersections_sample3 %>% group_by(X9,X10) %>% count()) %>% 
        group_by(X10) %>% count() %>% arrange(-n) %>%
        mutate(ratio = n/tot)

p1 %>% 
    inner_join(p2, by = c("X10" = "X10")) %>%
    inner_join(p3, by = c("X10" = "X10")) %>%
    inner_join(p4, by = c("X10" = "X10"))

rm(tot,p1,p2,p3,p4)

# ----------

# distribution of distances to the closest transcript

# importing distributions

p53_sites_closest_dw <- read_delim("bedtools-closest/p53-sites/closest_genes_downstream.txt", 
                                   "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
p53_sites_closest_up <- read_delim("bedtools-closest/p53-sites/closest_genes_upstream.txt", 
                                   "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

sample1_closest_dw <- read_delim("bedtools-closest/sample1/closest_genes_downstream.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sample1_closest_up <- read_delim("bedtools-closest/sample1/closest_genes_upstream.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

sample2_closest_dw <- read_delim("bedtools-closest/sample2/closest_genes_downstream.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sample2_closest_up <- read_delim("bedtools-closest/sample2/closest_genes_upstream.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

sample3_closest_dw <- read_delim("bedtools-closest/sample3/closest_genes_downstream.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
sample3_closest_up <- read_delim("bedtools-closest/sample3/closest_genes_upstream.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

p53_sites_closest_dw <- p53_sites_closest_dw %>% filter(X8 != ".")
p53_sites_closest_up <- p53_sites_closest_up %>% filter(X8 != ".")

sample1_closest_dw <- sample1_closest_dw %>% filter(X8 != ".")
sample1_closest_up <- sample1_closest_up %>% filter(X8 != ".")

sample2_closest_dw <- sample2_closest_dw %>% filter(X8 != ".")
sample2_closest_up <- sample2_closest_up %>% filter(X8 != ".")

sample3_closest_dw <- sample3_closest_dw %>% filter(X8 != ".")
sample3_closest_up <- sample3_closest_up %>% filter(X8 != ".")

# distance distributions

p53_dw <- p53_sites_closest_dw$X13
p53_up <- -1*(p53_sites_closest_up$X13)

r1_dw <- sample1_closest_dw$X13
r1_up <- -1*(sample1_closest_up$X13)

r2_dw <- sample2_closest_dw$X13
r2_up <- -1*(sample2_closest_up$X13)

r3_dw <- sample3_closest_dw$X13
r3_up <- -1*(sample3_closest_up$X13)

threshold <- 400000

dist_down <- data.frame(sample = "p53",
                        distance = p53_dw[p53_dw <= threshold])

dist_down <- rbind(dist_down, data.frame(sample = "random 1",
                                         distance = r1_dw[r1_dw <= threshold]))

dist_down <- rbind(dist_down, data.frame(sample = "random 2",
                                         distance = r2_dw[r2_dw <= threshold]))

dist_down <- rbind(dist_down, data.frame(sample = "random 3",
                                         distance = r3_dw[r3_dw <= threshold]))

a <- ggplot(dist_down, aes(x = sample, y = distance, fill = sample)) + 
    geom_boxplot() + 
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    theme_minimal() +
    theme(legend.position="none") + 
    labs(title = "Distance distribution to the closest downstream gene \n (distances limited to 400 000)", 
         x = "sample", y = "base pairs")

ggplot(dist_down, mapping = aes(x = distance, fill = sample)) + 
    geom_density(alpha = 0.5, position = "identity") +
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    theme_minimal() +
    theme(legend.position = c(0.9, 0.85))

dist_up <- data.frame(sample = "p53",
                        distance = p53_up[p53_up <= threshold])

dist_up <- rbind(dist_up, data.frame(sample = "random 1",
                                         distance = r1_up[r1_up <= threshold]))

dist_up <- rbind(dist_up, data.frame(sample = "random 2",
                                         distance = r2_up[r2_up <= threshold]))

dist_up <- rbind(dist_up, data.frame(sample = "random 3",
                                         distance = r3_up[r3_up <= threshold]))

b <- ggplot(dist_up, aes(x = sample, y = distance, fill = sample)) + 
    geom_boxplot() + 
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    theme_minimal() +
    theme(legend.position="none") +
    labs(title = "Distance distribution to the closest upstream gene \n (distances limited to 400 000)", 
         x = "sample", y = "base pairs")

ggplot(dist_up, mapping = aes(x = distance, fill = sample)) + 
    geom_density(alpha = 0.5, position = "identity") +
    scale_fill_manual(values = c("coral3", "dodgerblue1", "dodgerblue2", "dodgerblue3")) +
    theme_minimal() +
    theme(legend.position = c(0.9, 0.85))

library(ggpubr)
# 550x1100
ggarrange(a, b, ncol = 2, nrow = 1)

rm(dist_down)
rm(dist_up)
rm(threshold)
rm(a,b)

# biotype of the closest trascript

p53_sites_closest_dw <- p53_sites_closest_dw %>% 
    group_by(X10) %>% count() %>% arrange(-n) %>% 
    mutate(ratio = n/nrow(p53_sites_closest_dw))

p53_sites_closest_up <- p53_sites_closest_up %>% 
    group_by(X10) %>% count() %>% arrange(-n) %>% 
    mutate(ratio = n/nrow(p53_sites_closest_up))

sample1_closest_dw <- sample1_closest_dw %>% 
    group_by(X10) %>% count() %>% arrange(-n) %>% 
    mutate(ratio = n/nrow(sample1_closest_dw))

sample1_closest_up <- sample1_closest_up %>% 
    group_by(X10) %>% count() %>% arrange(-n) %>% 
    mutate(ratio = n/nrow(sample1_closest_up))

sample2_closest_dw <- sample2_closest_dw %>% 
    group_by(X10) %>% count() %>% arrange(-n) %>% 
    mutate(ratio = n/nrow(sample2_closest_dw))

sample2_closest_up <- sample2_closest_up %>% 
    group_by(X10) %>% count() %>% arrange(-n) %>% 
    mutate(ratio = n/nrow(sample2_closest_up))

sample3_closest_dw <- sample3_closest_dw %>% 
    group_by(X10) %>% count() %>% arrange(-n) %>% 
    mutate(ratio = n/nrow(sample3_closest_dw))

sample3_closest_up <- sample3_closest_up %>% 
    group_by(X10) %>% count() %>% arrange(-n) %>% 
    mutate(ratio = n/nrow(sample3_closest_up))

# Considering the sites not located on genes, this table shows, for each biotype, the fraction of sites having
# gene with such biotype as the closest one. The estimated p53 binding sites have more frequently a protein coding
# gene as the closest one. 

p53_sites_closest_dw %>% 
    inner_join(sample1_closest_dw, by = c("X10" = "X10")) %>% 
    inner_join(sample2_closest_dw, by = c("X10" = "X10")) %>% 
    inner_join(sample3_closest_dw, by = c("X10" = "X10")) %>%
    dplyr::select(X10, ratio.x, ratio.y, ratio.x.x, ratio.y.y)

p53_sites_closest_up %>% 
    inner_join(sample1_closest_up, by = c("X10" = "X10")) %>% 
    inner_join(sample2_closest_up, by = c("X10" = "X10")) %>% 
    inner_join(sample3_closest_up, by = c("X10" = "X10")) %>%
    dplyr::select(X10, ratio.x, ratio.y, ratio.x.x, ratio.y.y)

rm(p53_sites_closest_dw, p53_sites_closest_up)
rm(sample1_closest_dw, sample1_closest_up)
rm(sample2_closest_dw, sample2_closest_up)
rm(sample3_closest_dw, sample3_closest_up)

# ------------------------

# Considerations about normalization with respect to gene length

#Genes, by biotype, containing at least one site

tot_v35 <- (gencode_v35 %>% group_by(X5,X6) %>% summarise(m = n())) %>% group_by(X6) %>% count() %>% arrange(-n)
p53 <- (intersections_p53_v35 %>% group_by(X9,X10) %>% summarise(m = n())) %>% group_by(X10) %>% count() %>% arrange(-n)
r1 <- (intersections_sample1 %>% group_by(X9,X10) %>% summarise(m = n())) %>% group_by(X10) %>% count() %>% arrange(-n)
r2 <- (intersections_sample2 %>% group_by(X9,X10) %>% summarise(m = n())) %>% group_by(X10) %>% count() %>% arrange(-n)
r3 <- (intersections_sample3 %>% group_by(X9,X10) %>% summarise(m = n())) %>% group_by(X10) %>% count() %>% arrange(-n)

tot_v35 %>% 
    inner_join(p53, by = c("X6" = "X10")) %>%
    inner_join(r1, by = c("X6" = "X10")) %>%
    inner_join(r2, by = c("X6" = "X10")) %>%
    inner_join(r3, by = c("X6" = "X10")) %>%
    mutate(ratio_p53 = n.y/n.x) %>%
    mutate(ratio_r1 = n.x.x/n.x) %>%
    mutate(ratio_r2 = n.y.y/n.x) %>% 
    mutate(ratio_r3 = n/n.x) %>%
    dplyr::select(X6, n.x, n.y, ratio_p53, n.x.x, ratio_r1, n.y.y, ratio_r2, n, ratio_r3)

rm(p53,r1,r2,r3)

# gene lengths
merged_v35 %>% 
    mutate(len = X3-X2) %>%
    group_by(X5) %>% 
    summarise(count = n(), mean_len = mean(len), median_len = median(len)) %>%
    arrange(-count)

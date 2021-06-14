
basedir <- "/mnt/hdd0/univ/thesis-support-material/data/4-7/"

library("ViSEAGO")
library(dplyr)
library(readr)

# Genes DE in LINC KO
DE_linc01605ko <- read_delim(paste0(basedir, "DE_genes/no_suffix/DE_DESeq2_LINC01605ko.txt"), 
                             "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# ALL the genes in the annotation
all_genes_id_no_suffix <- read_delim(paste0(basedir, "DE_genes/no_suffix/all_genes_id_no_suffix.txt"), 
                                     "\t", escape_double = FALSE, col_names = FALSE, 
                                     trim_ws = TRUE)

# format the gene codes how viseago wants
selection <- (DE_linc01605ko %>% mutate(code = as.integer(sub("ENSG", "", X1))))$code
background <- (all_genes_id_no_suffix %>% mutate(code = as.integer(sub("ENSG", "", X1))))$code

# download GO annotation 
Bioconductor<-ViSEAGO::Bioconductor2GO()
#ViSEAGO::available_organisms(Bioconductor)
myGENE2GO<-ViSEAGO::annotate(
    "org.Hs.eg.db",
    Bioconductor )

# -- BP --

# perform GO terms enrichment analysis on BP terms
BP <- ViSEAGO::create_topGOdata(
        geneSel=selection,
        allGenes=background,
        gene2GO=myGENE2GO, 
        ont="BP",
        nodeSize=5 )

classic <- topGO::runTest(
    BP,
    algorithm ="classic",
    statistic = "fisher" )

# show the resulting GO terms
BP_sResults <- ViSEAGO::merge_enrich_terms(
    Input=list(
        condition=c("BP","classic")
     ) )
ViSEAGO::show_table(BP_sResults)

#ViSEAGO::show_table(
#    BP_sResults,
#    "BP_sResults.xls"
#)

# cluster the GO terms by similarity
myGOs <- ViSEAGO::build_GO_SS(
    gene2GO=myGENE2GO,
    enrich_GO_terms=BP_sResults
)

myGOs<-ViSEAGO::compute_SS_distances(
    myGOs,
    distance="Wang"
)

# visualize
ViSEAGO::MDSplot(myGOs)

# print MDSplot
#ViSEAGO::MDSplot(
#    myGOs,
#    file="mdsplot1.png" )

# 750x800
Wang_clusters_wardD2 <- ViSEAGO::GOterms_heatmap(
    myGOs,
    showIC=F,
    showGOlabels=T,
    GO.tree=list(
        tree=list(
            distance="Wang",
            aggreg.method="ward.D2"
        ),
        cut=list(
            dynamic=list(
                pamStage=TRUE,
                pamRespectsDendro=TRUE,
                deepSplit=2,
                minClusterSize =2
            )
        )
    ),
    samples.tree=NULL
)


# Display the clusters-heatmap
ViSEAGO::show_heatmap(
    Wang_clusters_wardD2,
    "GOterms"
)

# print the clusters-heatmap
#ViSEAGO::show_heatmap(
#    Wang_clusters_wardD2,
#    "GOterms",
#    "cluster_heatmap_Wang_wardD2.png" )


# Display the clusters-heatmap table
ViSEAGO::show_table(Wang_clusters_wardD2)

# display colored MDSplot
ViSEAGO::MDSplot(
    Wang_clusters_wardD2,
    "GOterms" )

# ------------------------

# calculate semantic similarities between clusters of GO terms
Wang_clusters_wardD2 <- ViSEAGO::compute_SS_distances(
    Wang_clusters_wardD2,
    distance=c("max", "avg","rcmax", "BMA")
)

ViSEAGO::MDSplot(
    Wang_clusters_wardD2,
    "GOclusters"
)

# GOclusters heatmap
Wang_clusters_wardD2 <- ViSEAGO::GOclusters_heatmap(
    Wang_clusters_wardD2,
    tree=list(
        distance="BMA",
        aggreg.method="ward.D2"
    )
)

# display the GOClusters heatmap
ViSEAGO::show_heatmap(
    Wang_clusters_wardD2,
    "GOclusters"
)

# -- MF --

# perform GO terms enrichment analysis on MF terms
MF <- ViSEAGO::create_topGOdata(
    geneSel=selection,
    allGenes=background,
    gene2GO=myGENE2GO, 
    ont="MF",
    nodeSize=5 )

classic<-topGO::runTest(
    MF,
    algorithm ="classic",
    statistic = "ks" )

# show the resulting GO terms
MF_sResults <- ViSEAGO::merge_enrich_terms(
    Input=list(
        condition=c("MF","classic")
    ) )
ViSEAGO::show_table(MF_sResults)

#ViSEAGO::show_table(
#    BP_sResults,
#    "BP_sResults.xls"
#)

# cluster the GO terms by similarity
myGOs <- ViSEAGO::build_GO_SS(
    gene2GO=myGENE2GO,
    enrich_GO_terms=MF_sResults
)

myGOs<-ViSEAGO::compute_SS_distances(
    myGOs,
    distance="Wang"
)

# visualize
ViSEAGO::MDSplot(myGOs)

# print MDSplot
#ViSEAGO::MDSplot(
#    myGOs,
#    file="mdsplot1.png" )

Wang_clusters_wardD2 <- ViSEAGO::GOterms_heatmap(
    myGOs,
    showIC=F,
    showGOlabels=TRUE,
    GO.tree=list(
        tree=list(
            distance="Wang",
            aggreg.method="ward.D2"
        ),
        cut=list(
            dynamic=list(
                pamStage=TRUE,
                pamRespectsDendro=TRUE,
                deepSplit=2,
                minClusterSize =2
            )
        )
    ),
    samples.tree=NULL
)


# Display the clusters-heatmap
ViSEAGO::show_heatmap(
    Wang_clusters_wardD2,
    "GOterms"
)

# print the clusters-heatmap
#ViSEAGO::show_heatmap(
#    Wang_clusters_wardD2,
#    "GOterms",
#    "cluster_heatmap_Wang_wardD2.png" )


# Display the clusters-heatmap table
ViSEAGO::show_table(Wang_clusters_wardD2)

# display colored MDSplot
ViSEAGO::MDSplot(
    Wang_clusters_wardD2,
    "GOterms" )

# ------------------------

# calculate semantic similarities between clusters of GO terms
Wang_clusters_wardD2 <- ViSEAGO::compute_SS_distances(
    Wang_clusters_wardD2,
    distance=c("max", "avg","rcmax", "BMA")
)

ViSEAGO::MDSplot(
    Wang_clusters_wardD2,
    "GOclusters"
)

# GOclusters heatmap
Wang_clusters_wardD2 <- ViSEAGO::GOclusters_heatmap(
    Wang_clusters_wardD2,
    tree=list(
        distance="BMA",
        aggreg.method="ward.D2"
    )
)

# display the GOClusters heatmap
ViSEAGO::show_heatmap(
    Wang_clusters_wardD2,
    "GOclusters"
)

# -- CC --

# perform GO terms enrichment analysis on MF terms
CC <- ViSEAGO::create_topGOdata(
    geneSel=selection,
    allGenes=background,
    gene2GO=myGENE2GO, 
    ont="CC",
    nodeSize=5 )

classic<-topGO::runTest(
    CC,
    algorithm ="classic",
    statistic = "fisher" )

# show the resulting GO terms
CC_sResults <- ViSEAGO::merge_enrich_terms(
    Input=list(
        condition=c("CC","classic")
    ) )
ViSEAGO::show_table(CC_sResults)

#ViSEAGO::show_table(
#    BP_sResults,
#    "BP_sResults.xls"
#)

# cluster the GO terms by similarity
myGOs <- ViSEAGO::build_GO_SS(
    gene2GO=myGENE2GO,
    enrich_GO_terms=CC_sResults
)

myGOs<-ViSEAGO::compute_SS_distances(
    myGOs,
    distance="Wang"
)

# visualize
ViSEAGO::MDSplot(myGOs)

# print MDSplot
#ViSEAGO::MDSplot(
#    myGOs,
#    file="mdsplot1.png" )

Wang_clusters_wardD2 <- ViSEAGO::GOterms_heatmap(
    myGOs,
    showIC=TRUE,
    showGOlabels=TRUE,
    GO.tree=list(
        tree=list(
            distance="Wang",
            aggreg.method="ward.D2"
        ),
        cut=list(
            dynamic=list(
                pamStage=TRUE,
                pamRespectsDendro=TRUE,
                deepSplit=2,
                minClusterSize =2
            )
        )
    ),
    samples.tree=NULL
)


# Display the clusters-heatmap
ViSEAGO::show_heatmap(
    Wang_clusters_wardD2,
    "GOterms"
)

# print the clusters-heatmap
#ViSEAGO::show_heatmap(
#    Wang_clusters_wardD2,
#    "GOterms",
#    "cluster_heatmap_Wang_wardD2.png" )


# Display the clusters-heatmap table
ViSEAGO::show_table(Wang_clusters_wardD2)

# display colored MDSplot
ViSEAGO::MDSplot(
    Wang_clusters_wardD2,
    "GOterms" )

# ------------------------

# calculate semantic similarities between clusters of GO terms
Wang_clusters_wardD2 <- ViSEAGO::compute_SS_distances(
    Wang_clusters_wardD2,
    distance=c("max", "avg","rcmax", "BMA")
)

ViSEAGO::MDSplot(
    Wang_clusters_wardD2,
    "GOclusters"
)

# GOclusters heatmap
Wang_clusters_wardD2 <- ViSEAGO::GOclusters_heatmap(
    Wang_clusters_wardD2,
    tree=list(
        distance="BMA",
        aggreg.method="ward.D2"
    )
)

# display the GOClusters heatmap
ViSEAGO::show_heatmap(
    Wang_clusters_wardD2,
    "GOclusters"
)


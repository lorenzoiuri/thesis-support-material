
# using mean as linkage criterion

# GOType must be "BP", "MF", "CC"
computeSimilarity <- function(first, 
                              second, 
                              GOType, 
                              preComputedMatrix = NULL, 
                              naiveMatrixCompute = F,
                              precomputed_hs_path = "/mnt/hdd0/univ/thesis-support-material/data/4-7/comparison_metric/precomputed_hs",
                              precomputed_matrix_path = "/mnt/hdd0/univ/thesis-support-material/data/4-7/comparison_metric/precomputed_matrix",
                              exportMatrix = F,
                              experiment){

    # setting frequency threshold to discard uninformative GO terms
    freq.threshold <- 16
    
    p <- first %>% 
        dplyr::filter(Frequency < freq.threshold) %>%
        dplyr::mutate(code = as.numeric(substr(TermID, 4, 999))) %>%
        dplyr::mutate(Representative = ifelse(Representative == -1, code, Representative))
    
    g <- second %>% 
        dplyr::filter(Frequency < freq.threshold) %>%
        dplyr::mutate(code = as.numeric(substr(TermID, 4, 999))) %>%
        dplyr::mutate(Representative = ifelse(Representative == -1, code, Representative))
    
    # computing first part of the score
    print("computing the first part of the score...")
    
    count <- 0
    not_in_g <- p %>% dplyr::filter(TermID == 1) # empty
    for(i in 1:nrow(p)){
        # only consider the leader
        if(p[i,]$Eliminated == F){
            # considering now all the terms in its cluster
            code <- p[i,]$code
            friends <- p %>% dplyr::filter(Representative == !!code)
            #friends <- p[p$Representative == code, ]
            # checking if in g exists one of the terms in the cluster of i
            res <- friends %>% dplyr::inner_join(g, by = c("TermID" = "TermID")) %>% nrow()
            if (res > 0) {
                count <- count + 1
            } else {
                not_in_g <- rbind(not_in_g, p[i,])
            }
            rm(code, friends, res)
        }
    }
    
    part_a <- count
    rm(count)
    print("done")
    
    # if preComputedMatrix is NULL, then the similarities between the GO terms need to be computed
    # otherwise they are imported from the preComputedMatrix file path
    
    simMatrix <- NULL
    # GO IDs of the terms belonging to a p cluster not in g
    go1 <- (p %>% dplyr::filter(Representative %in% not_in_g$Representative))$TermID
    # GO IDs of the terms in g
    go2 <- g$TermID
    
    if(length(go1) == 0){
        part_b <- 0
    } else {
        if (is.null(preComputedMatrix)){
            print("computing similarity matrix...")
            suppressPackageStartupMessages(suppressWarnings(suppressMessages(library(GOSemSim))))
            
            hsGO <- NULL
            setwd(precomputed_hs_path)
            if(GOType == "BP"){
                hsGO <- readRDS(file = "bp.rds")
            } else if(GOType == "MF"){
                hsGO <- readRDS(file = "mf.rds")
            } else if(GOType == "CC"){
                hsGO <- readRDS(file = "cc.rds")
            } else {
                print("Error: Unknown GOType")
                exit()
            }
                
            #hsGO <- godata('org.Hs.eg.db', ont=GOType, computeIC=T)
        
            # computing the preComputedMatrix  between GO terms
            simMatrix <- NULL
            if (naiveMatrixCompute){
                print("using inefficient way to compute similarity Matrix")
                simMatrix <- matrix(data = 0, nrow = length(go1), ncol = length(go2))
                
                for (i in 1:nrow(simMatrix)) {
                    for (j in 1:ncol(simMatrix)){
                        simMatrix[i,j] <- goSim(go1[i], go2[j], semData=hsGO, measure="Rel")
                    }
                    print(paste0("computing row ", i, " of ", nrow(simMatrix)))
                }
                
                if(exportMatrix == T){
                    setwd(precomputed_matrix_path)
                    print("exporting computed matrix")
                    a <- deparse(substitute(first))
                    b <- deparse(substitute(second))
                    saveRDS(simMatrix, file = paste0(a,"-", b, "-", GOType, "-", experiment))
                    rm(a,b)
                }
                
            } else {
                print("using efficient way to compute similarity Matrix")
                simMatrix <- mgoSim(go1, go2, semData=hsGO, measure="Rel", combine=NULL)
            }
        } else {
            print("importing similarities matrix...")
            setwd(precomputed_matrix_path)
            simMatrix <- readRDS(file = preComputedMatrix)
        }
        simMatrix[is.na(simMatrix)] <- 0
        rownames(simMatrix) <- go1
        colnames(simMatrix) <- go2
        print("done")
        rm(go1, go2)
        #print(simMatrix)
        
        # computing the second part of the score
        print("computing the second part of the score...")
        
        sim.threshold <- 0.1
        score <- 0
        g_clusters <- g %>% filter(Eliminated == F)
        for(i in 1:nrow(not_in_g)){
            code <- not_in_g[i,]$code
            friends <- p %>% filter(Representative == !!code)
            max_mean_similarity <- 0
            max_g_cluster_code <- 0
            for(j in 1:nrow(g_clusters)){
                g_cluster_code <- g_clusters[j,]$code
                g_cluster_friends <- g %>% filter(Representative == !!g_cluster_code)
                similarity_values <- c()
                for(k in 1:nrow(friends)){
                    for(h in 1:nrow(g_cluster_friends)){
                        value <- simMatrix[friends[[k,1]], g_cluster_friends[[h,1]] ]
                        if (is.na(value)){
                            #print("warning: NA generated, replaced with 0")
                            value <- 0
                        }
                        similarity_values <- c(similarity_values, value)
                    }
                }
                mean_similarity <- mean(similarity_values)
                if (mean_similarity > max_mean_similarity) {
                    max_mean_similarity <- mean_similarity
                    max_g_cluster_code <- g_cluster_code
                }
            }
            if(max_mean_similarity >= sim.threshold){
                score <- score + max_mean_similarity
                #print(paste0(code, " ~ ", max_g_cluster_code, " : ", max_mean_similarity))
            }
            #print(paste0((i)/nrow(not_in_g)*100, " %..."))
        }
        part_b <- score
        rm(score)
        print("done")
        rm(similarity_values, mean_similarity, i, j, k, h)
        rm(sim.threshold, value)
        rm(friends, code, max_mean_similarity)
        rm(g_cluster_code, g_cluster_friends, g_clusters, max_g_cluster_code)
    }
    
    scale.factor <- 1
    representatives <- (p %>% dplyr::filter(Eliminated == F) %>% count())$n
    print(paste0("first part: ", part_a, " / ", representatives))
    print(paste0("second part: ", part_b, " / ", representatives))
    score <- (part_a + part_b)/representatives
    print(paste0("final score: ", score))
    return(score)
}

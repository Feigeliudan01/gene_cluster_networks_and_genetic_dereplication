#!/usr/bin/env Rscript
library(optparse)

library(igraph)
library(parallel)
library(ggtree)

getwd()


log_con <- file("test_run.log")

args<-commandArgs(TRUE)


print("Working on file")
print(args[1])
filename = args[1]
print("Changing directory to")
print(args[2])
setwd(args[2])


aspDf <- read.table(filename, header = TRUE,  stringsAsFactors = FALSE, sep = '\t', quote = "\n")


clusterBlastAll <- read.csv("clusterBlastAll.csv", header = TRUE, stringsAsFactors = FALSE)

clusterBlast <- clusterBlastAll
clusterBlast <- clusterBlast[clusterBlast$pident_score != 0, c('q_clust_id', 'h_clust_id', 'pident_score')]

names(clusterBlast)[3] <- 'weight'

print("Testing aspDf")
head(aspDf)

cat("Missing gene clusters in clusterblast query", file = log_con, append = TRUE)
cat(paste(setdiff(clusterBlast$q_clust_id, aspDf$cluster_id), collapse = ",") , file = log_con, append = TRUE)
setdiff(clusterBlast$q_clust_id, aspDf$cluster_id)
cat("Missing gene clusters in downloaded data frame", file = log_con, append = TRUE)
cat(paste(setdiff(aspDf$cluster_id, clusterBlast$q_clust_id), collapse = ",") , file = log_con, append = TRUE)
missingOnes <- setdiff(aspDf$cluster_id, clusterBlast$q_clust_id)
print(missingOnes)
cat(paste(setdiff(clusterBlast$h_clust_id, aspDf$cluster_id), collapse = ","), file = log_con, append = TRUE)

clusterBlast <- clusterBlast[(clusterBlast$q_clust_id %in% aspDf$cluster_id & clusterBlast$h_clust_id %in% aspDf$cluster_id),]
print("Checking clusterBlast")


clusterDups <- clusterBlast[duplicated(clusterBlast[,c('q_clust_id', 'h_clust_id')]),c('q_clust_id', 'h_clust_id')]
cat(paste(clusterDups, collapse = ","), file = log_con, append = TRUE)

sum(clusterBlast$q_clust_id == clusterBlast$h_clust_id)
clusterDups

write.csv( clusterBlast, file = "clusterBlast.csv", row.names = FALSE )


set.seed(123)

cat("Info about clusterBLast", file = log_con, append = TRUE)

# clusterBlastSub <- clusterBlast[clusterBlast$weight >20,]
cat("Minimum clusterBlast weight", file = log_con, append = TRUE)
cat(min(clusterBlast$weight), file = log_con, append = TRUE)

cat("ClusterBlast weight over 100", file = log_con, append = TRUE)
# cat(clusterBlastAll[clusterBlastAll$pident_score > 100,], file = log_con, append = TRUE)


print("Creating Network")
present_clusters <- unique(aspDf$cluster_id[aspDf$cluster_id %in% clusterBlast$q_clust_id | aspDf$cluster_id %in% clusterBlast$h_clust_id])
g.raw <- graph_from_data_frame(d=clusterBlast, vertices= present_clusters, directed=FALSE) # Decided on undirected graph because I it's a bidirectional relationship for all anyway
g <- igraph::simplify(g.raw, edge.attr.comb = list(weight = "mean")) #remove.multiple = TRUE)
walks <- cluster_walktrap(as.undirected(g), steps = 1)

if(max(as.numeric(sizes(walks))) > length(unique(aspDf$org_id)) ){
    print("Network families are too large, needs more clustering")
    print("There are too many gene clusters per family, since they are exceeding the number of orgs (we assume that duplications will form another family, hence another round of clustering)")

   superNetworks <- as.numeric(which(sizes(walks) > length(unique(aspDf$org_id))))

     firstNetworks <- lapply(names(groups(walks))[!(names(groups(walks)) %in% superNetworks)],function(x){
         sub <- induced.subgraph(g, vids = walks[[x]])
         })

     secondNetwork <- lapply(names(groups(walks))[names(groups(walks)) %in% superNetworks],function(x){
         sub <- induced.subgraph(g, vids = walks[[x]])
         walkMore <- cluster_walktrap(sub, steps = 1)
         lapply(names(groups(walkMore)), function(x){
             induced.subgraph(g, vids = walkMore[[x]])
             } )
         })

     testg <- induced_subgraph(g, vids = walks[[superNetworks[[1]]]])

     secondNetwork <- unlist(secondNetwork, recursive = FALSE)
     cNetwork <- c(firstNetworks, secondNetwork) # complete network
     sizes <- sapply(cNetwork, function(x){length(V(x))})
          counter = 0
 clusterFamDf <- lapply(cNetwork, function(x){
     counter <<- counter + 1
     clusters <- V(x)$name
     data.frame(cluster_id = clusters, clusterFam = rep(counter, length(clusters)),
     stringsAsFactors = FALSE)
     })


 clusterFamDf <- do.call(rbind, clusterFamDf)
 clusterFamDf$cluster_id <- as.character(clusterFamDf$cluster_id)

famLimit <-  max(clusterFamDf$clusterFam)
names(clusterFamDf)
missingClusters <- setdiff(unique(aspDf$cluster_id), clusterFamDf$cluster_id)
missingDf <- data.frame(cluster_id = missingClusters, clusterFam = seq(from = famLimit+1, to = famLimit+length(missingClusters), by = 1))

clusterFamDf <- rbind(clusterFamDf, missingDf)

viFams <- merge(aspDf, clusterFamDf, by = c("cluster_id"), all.x = TRUE)


}else{

    print("No further clustering")
    coms <- communities(walks)
    tmpFam <- lapply(names(coms), function(x){
        data.frame(cluster_id = coms[[x]], clusterFam = (x), stringsAsFactors = FALSE)
        })

    famDf <- do.call(rbind, tmpFam)

    str(famDf)
    head(aspDf)
    viFams <- merge(aspDf, famDf, by = c("cluster_id"), all.x = TRUE)
    str(viFams)

}

write.table(viFams, paste0(gsub(".tsv","", filename),"_c.tsv"), row.names = FALSE, sep = "\t")



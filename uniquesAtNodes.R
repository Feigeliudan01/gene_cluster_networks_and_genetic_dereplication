#!/usr/bin/env Rscript
library(ape)
library(ggtree)

args<-commandArgs(TRUE)
print(args[1])
setwd(args[2])
filename = args[1]
treeFile = args[3] 

aspDf <- read.table(filename, stringsAsFactors = FALSE, header = TRUE, sep = "\t")

dummyTree <- read.tree(treeFile) # Tree needs new organism ids

p <- ggtree(dummyTree, branch.length = "none")

getUniquesAtNodes <- function(p, tree, clusters){
  "function changed to return df, changed function to take name"

  nodeContent <- lapply(p$data$node, function(x){
    get.offspring.tip(tree, p$data$node[x])
  })

  names(nodeContent) <- p$data$node

  orgsByFam <- split(clusters$name, clusters$clusterFam)

  # Comparing orgs at nodes and in families

  ndcUnique <- lapply(nodeContent, function(x){
    pattern <- sapply(orgsByFam, function(y){ all(x  %in% y) & all(y %in% x)})
    as.integer(names(orgsByFam)[pattern])
  })


  tmpUniques <- lapply(names(ndcUnique), function(x){
    tmpFams <- as.character(ndcUnique[[x]])
    if(length(tmpFams) == 0){ tmpFams <- "none"}
    data.frame(clusterFam = tmpFams, nodePosition = as.character(x) )
  })

  unNodeDf <- do.call(rbind, tmpUniques)

  return(unNodeDf)
}


uFamsAtNodes <- getUniquesAtNodes(p, dummyTree, aspDf)

viFamsNodes <- merge(aspDf, uFamsAtNodes, by = "clusterFam", all.x = TRUE)

write.table(viFamsNodes, paste0(gsub(".tsv","", filename),"_n.tsv"), row.names = TRUE)

uNodes <- lapply(split(uFamsAtNodes$clusterFam, uFamsAtNodes$nodePosition),function(x){
    if(x == 'none'){res = 0
        }else{
            res <- length(x)
        }
        res
    })

uNodes <- data.frame(counts = do.call(rbind, uNodes), node = as.integer(names(uNodes)))

write.table(uNodes, sprintf("%s_nodePosition.tsv", gsub(".tsv","", filename)), row.names = FALSE)

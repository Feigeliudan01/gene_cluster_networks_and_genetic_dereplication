setwd("../nigri_test") # Needs to be changed to output directory of main pipeline
source("../smModule/gvizPreprocessor.R")
sm <- read.table("sm_data_nigri_test_c_n.tsv", header = TRUE, stringsAsFactors = FALSE)

sm$gff_start <- as.numeric(sm$gff_start)
sm$gff_end <- as.numeric(sm$gff_end)

sm$pub_name <- paste(substr(sm$genus,1,1), sm$real_name)

malfOrgs <-  c("Asplac1", "Aspni_NRRL3_1","Asptu1", "Aspni_DSM_1", "Aspbr1", "Aspni7")


smList <- split(sm, sm$clusterFam)

resultList <- lapply(smList, function(x){
  tmp <- subset(x, sm_short == "NRPS")
  
  NrpsScore <- sum((tmp$gff_end - tmp$gff_start) > 14400)
  
  if(is.na(NrpsScore)){    NrpsScore <- 0}
  
  orgsPresent = length(intersect(x$name, malfOrgs))/ length(malfOrgs)
    
  data.frame(interproScore = sum(grepl("Glutathione|disulphide", x$ipr_desc)),
             rightNrps = NrpsScore,
             orgScore = orgsPresent
             )
  
})

result <- do.call(rbind, resultList)

result$clusterFam <- rownames(result)

candidates <- subset(result, subset = (orgScore > 0.8 & rightNrps > 6 & interproScore > 0))

candidates <- candidates[base::order(candidates$interproScore, decreasing = TRUE),]

bestCandidate <- candidates[1,]$clusterFam

clust <- subset(sm, clusterFam == bestCandidate) 

for(i in unique(clust$cluster_id)){
  print(i)
}

exchangeList <- list(c("DNA-binding|Transcription factor","Transcription factor"),
     c("Methyltransferase", "Methyltransferase"),
       c("Glutathione S-transferase", "Glutathione S-transferase"),
     c("pyridine nucleotide-disulphide oxidoreductase", "pyridine nucleotide-disulphide oxidoreductase"),
     c("Aminotransferase", "Aminotransferase")
     )


clust$ipr_desc[grepl("AMP-dependent synthetase/ligase",
                     clust$ipr_desc ) &
                 !grepl("AMP-dependent synthetase/ligase,Phosphopantetheine-binding,Condensation domain",
                        clust$ipr_desc)] <- "AMP-dependent synthetase/ligase-like"

for(i in exchangeList){
  clust$ipr_desc[grepl(i[1], clust$ipr_desc)] <- i[2]
}

clustList <- aspGviz(clust)


annFrame <- unique(
  subset(
    do.call(rbind, clustList),
    select = c(feature, annotation))
)
annFrame <- annFrame[
  order(
    match(annFrame$feature, LETTERS702)
  )
  ,]





annFrame



mytracks <- createTrackList(clustList)


pdf("malfClusters.pdf", width = 7, height = 5)
plotTracks(mytracks, groupAnnotation = "group" ,
           extend.left = 10000,
           shape = "arrow",
           fontfamily.group = "Helvetica",
           fontface.group = 3,
           A = "#FFCC00", B = "#660066", E = "#CC3399",
           H = "#66CC00", P = "#3366CC" , Q = "#33CCFF",
           U = "#FF0099", Y ="#990000", AJ = "#660099",
           AK = "#FF0000", N = "#FF9900"
           )
dev.off()

pdf("malfClustersWONames.pdf", width = 12, height = 8)
plotTracks(mytracks,
           shape = "arrow",
           fontfamily.group = "Helvetica",
           fontface.group = 3,
           A = "#FFCC00", B = "#660066", E = "#CC3399",
           H = "#66CC00", P = "#3366CC" , Q = "#33CCFF",
           U = "#FF0099", Y ="#990000", AJ = "#660099",
           AK = "#FF0000", N = "#FF9900"
)
dev.off()

selFeatures <- c("A", "B", "E" ,"H", "P", "Q", "U" , "Y" , "AJ", "AK", "N")



x11()

gvizCols <- c("#FFCC00", "#660066", "#CC3399",
              "#66CC00","#3366CC" , "#33CCFF",
              "#FF0099", "#990000", "#660099",
              "#FF0000","#FF9900")




pdf("legend.pdf")
plot.new()
legend("topleft", legend= annFrame$annotation[annFrame$feature %in% selFeatures], pt.cex = 2, pch = 21,   pt.bg=gvizCols)
dev.off()
# Load clusterBlast and show how similar they are

clusterBlast <- read.csv("clusterBlast.csv", stringsAsFactors = FALSE)
library(igraph)
g.raw <- graph_from_data_frame(d=clusterBlast, vertices= unique(sm$cluster_id), directed=FALSE) # Decided on undirected graph because I it's a bidirectional relationship for all anyway
g <- igraph::simplify(g.raw, edge.attr.comb = list(weight = "mean"))


malfNetwork <- induced.subgraph(g, vids = subset(sm, subset = clusterFam == bestCandidate)$cluster_id)

ceb <- cluster_edge_betweenness(malfNetwork)

names(membership(ceb))
as.numeric(membership(ceb))
plot(ceb, malfNetwork)

library(ggraph)
library(RColorBrewer)

V(malfNetwork)
V(malfNetwork)$clustMember <- c("core", "outer", "outer")[as.numeric(membership(ceb))]

cols <- colorRampPalette(brewer.pal(9,"Blues"))(100)
bcols <- brewer.pal(9,"Blues")

display.brewer.all()

gp <- ggraph(malfNetwork) + 
  geom_edge_link(aes(colour = weight), edge_width = .3) + 
  geom_node_point(size = 3, aes(color = clustMember))+
  scale_edge_color_gradient(low = bcols[2], high = bcols[8])+
  theme_graph(base_family = "Helvetica")



ggsave("malfGraph2col.pdf", plot = gp, width = 7, height = 7)

ggsave("malfGraph15col.pdf", plot = gp, width = 4.5, height = 2.4)

ggsave("malfGraph1col.pdf", plot = gp, width = 3.42, height = 2.4)



library(Gviz)

LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
LETTERS702 <- LETTERS702[c(-6, -7)]
MYLETTERS <- LETTERS[c(-6,-7)]




modCluster <- function(df){
  
  # Normalizes cluster coordinates and turns clusters if the cluster
  # backbone is located on the '-' strand.
  # Input: Dataframe with gff coordinates and strand.
  # Changes direction of cluster if necessary (see 'change')
  
  if(sum(c("protein_id", "gff_start", "gff_end", "gff_strand") %in%  colnames(df)) != 4 ){
    stop("protein_id and gff information must be provided")}
  
  df$start <- df$gff_start - min(df$gff_start)
  df$end <- df$start + df$gff_end-df$gff_start
  df$strand <- df$gff_strand
  
  bb = strsplit(unique(df$group), '_')[[1]][2] # reads bb protein from clsuter_id
  
  change = sum(df$protein_id == bb & df$strand == '-')
  
  if(change){
    
    old_orientation = as.numeric(df$strand == '-')
    df$strand <- c('-','+')[old_orientation+1] # Changing strand
    
    diff <- df$end-df$start
    df$start <- abs(df$end - max(df$end))
    df$end <- df$start + diff
    
    df <- df[rev(rownames(df)),]
    
  }
  return(df)
}


addFeatureColumn <- function(cdf){
  # Adding an additional column to gene cluster df
  # for better diplay in plots
  
  cdf$annotation <- cdf$ipr_desc
  cdf$annotation[cdf$sm_short != 'none'] <- cdf$sm_short[cdf$sm_short != 'none']
  # Not used right now
  
  annTemp <- as.data.frame(table(cdf$annotation[cdf$annotation != 'none']))
  annTemp$Var1 <-  as.character(annTemp$Var1)
  annOrder <- annTemp$Var1[order(annTemp$Freq, decreasing = TRUE)]
  cdf$annotation <- factor(cdf$annotation, levels = annOrder)
  
  if(length(unique(cdf$annotation)) < length(LETTERS702)){
    cdf$feature <- LETTERS702[cdf$annotation]
    cdf$feature[cdf$sm_short != 'none'] <- as.character(cdf$annotation[cdf$sm_short != 'none'])
    # cdf$feature[cdf$annotation == 'none'] <- 'none'
    return(cdf)
    
  }else{
    print("WARNING! Not enough letters")
  }
}


aspGviz <- function(sm) {
  # Needs data frame with gene cluster data (ipr and gff and sm_short)
  # Returns list of clusters for AnnotationTrack
  
  if(sum(is.na(sm$gff_start)) > 0 ){
    "Warning: There are gene clusters without gff coordinates"
  }
  
  sm <- subset(sm, !is.na(gff_start) & !is.na(gff_end))
  
  
  sm$group <- sm$cluster_id
  sm <- addFeatureColumn(sm)
  
  perCluster <- split(sm, sm$group)
  
  smList <- lapply(perCluster, function(x){
    x <- modCluster(x)
  })
  # sm <- do.call(rbind, newPosition)
  
  return(smList)
  
}


getNiceLegend <- function(clustList){
  
  if(class(clustList)!="list"){stop("Input must be list")}
  if(sum(c("feature", "annotation") %in%  colnames(clustList[[1]])) !=2 ){
    stop("Input must contain feature and annotation columns")}
  
  # Selects featues and their description to create
  # an overview of features to print to a legend
  
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
  
  return(annFrame)
}


createTrackList <- function(clustList){
  # Automatically creates AnnotationTracks
  # of gene clusters
  if(sum(c("feature", "complete_name") %in%  colnames(clustList)) ){
    stop("Input must contain feature and complete_name columns")}
  
  mytracks <- GenomeAxisTrack() # Initialize track
  
  aTrack <- lapply(seq(1,length(clustList)), function(x){ # create one track per cluster
    
    AnnotationTrack(clustList[[x]],
                    name = names(clustList)[x],
                    feature = clustList[[x]]$feature,
                    group = clustList[[x]]$complete_name) # change to pub_name or complete Name
  })
  
  mytracks <- append(mytracks, aTrack) # append to initialized track
  
  return(mytracks)
}

#  setwd("/home/seth/aspmineSM/NigriPaperSet")
# #
#  sm <- read.table("smurfOrgIprMibigFamsNodes.tsv", header = TRUE, stringsAsFactors = FALSE)
#  sm <- subset(sm, clusterFam == 2, select = c("clusterFam", "cluster_id", "protein_id","name","complete_name", "sm_short","clust_size","gff_start","gff_end","gff_strand","ipr_desc"))
# colnames(sm)
#
# write.table(sm, "syntenyTest.tsv", row.names = FALSE)
#
#
# sm <- read.table("/home/seth/asptoolbox/misc/syntenyTest.tsv",
#                  stringsAsFactors = FALSE,
#                  header = TRUE)
#
#
#  clustList <- aspGviz(sm)
#
#
# # library(Gviz)
# mytracks <- GenomeAxisTrack()
#
# aTrack <- lapply(seq(1,length(clustList)), function(x){
#   AnnotationTrack(clustList[[x]], name = names(clustList)[x], feature = clustList[[x]]$feature, group = clustList[[x]]$complete_name)
# })
#
#  mytracks <- append(mytracks, aTrack)
#  plotTracks(mytracks, groupAnnotation = "group" , featureAnnotation = "feature",A = "#FFCC00", B = "#660066", C = "#CC3399", D = "#66CC00", E = "#3366CC" , H = "#33CCFF", I = "#FF0099", J ="#990000", K = "#660099", L = "#FF0000", PKS = "#FF9900")
#
#  print(getNiceLegend(clustList), row.names = FALSE)
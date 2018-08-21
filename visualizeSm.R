#!/usr/bin/env Rscript
library(gplots)
library(RColorBrewer)
library(ape)
library(ggtree)
library(ggstance)

##############################
# USE AS RSCRIPT visualizeSm.R fileName setName treeFile
# From the secMet folder if needed



#############
# ARGUMENTS
####
args<-commandArgs(TRUE)
fileName = args[1]
setName = args[2]
treeFile = args[3]


setwd(setName)

sm <- read.table(fileName, header = TRUE, stringsAsFactors = FALSE )

sm$complete_name <- paste0(substr(sm$genus,1,1),'. ', sm$real_name)
sm$compound <- gsub( " biosynthetic gene cluster", "", sm$compound)
sm$compound[sm$compound == "none"] <- NA

subset(sm, compound == "Patulin")
patulinFam <- subset(sm, compound == "Patulin")$clusterFam
aculinic_acid_protein_id <- subset(sm, (clusterFam == patulinFam & name == "Aspac1" & sm_short != "none"))$protein_id
sm[sm$protein_id == aculinic_acid_protein_id & sm$name == "Aspac1",]$compound <- "Aculinic acid"

# "4,4'-piperazine-2,5-diyldimethyl-bis-phenol" replaced by piperazine*
sm[grepl("piperazine",sm$compound),]$compound <- "Piperazine*"
sm <- subset(sm, sm$name!="Aspcal1")



sm[(sm$protein_id ==5447 & sm$name == "Aspnid1"),]$compound <- "Nidulanin A"
# http://www.aspergillusgenome.org/cgi-bin/locus.pl?locus=AN0150 
# mpdG AN0150
sm[(sm$protein_id ==6652 & sm$name == "Aspnid1"),]$compound <- "Emodin"


tree <- read.tree(paste0("../",treeFile))

p <- ggtree(tree, branch.length = 'none')
dfNames <- unique(sm$name)
treeNames <-  unique(subset(p$data, isTip, select = label))$label

tree_select <- ape::drop.tip(tree, setdiff(treeNames, dfNames))

p_select <- ggtree(tree_select)


tipsY <- subset(p_select$data, isTip, select = c(label, y))
tipOrder <- tipsY$label[order(tipsY$y, decreasing = TRUE)]


nameTrans <- unique(subset(sm, select = c(name, complete_name)))
tipsComplete <- nameTrans$complete_name[match(tipOrder, nameTrans$name)]

unique(sm$name)
tree$tip.label

orgsMissingInTree <- setdiff(unique(sm$name), tree$tip.label)
orgsMissingInDf <- setdiff(tree$tip.label, unique(sm$name))



#####################
# Family distribution
######
familyFrequencies <- as.data.frame(table(unique(subset(sm, select = c(clusterFam, name)))$clusterFam))
library(ggplot2)
ggplot(familyFrequencies) +
  geom_histogram(aes(Freq), binwidth = 1, colour = "black")+
  stat_bin(binwidth=1, geom="text", colour="black", size=3.5,
           aes(x = Freq, label=..count..), vjust = -1)+
  ggtitle("Number of SMGC families with clusters in x number of orgs")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Number of orgs in family")+
  ylab("Number of families")

ggsave("Family size distribution.pdf", device = "pdf")

length(unique(sm$cluster_id))
length(unique(sm$name))
table(familyFrequencies$Freq)
subset(sm, is.na(sm$clusterFam))
max(sm$clusterFam[!is.na(sm$clusterFam)])




############################
# HEATMAP SHARED GENE CLUSTERS BETWEEN ORGS


orgsVsOrgsRaw <- split(unique(sm[,c("clusterFam","complete_name")])$clusterFam,
                       unique(sm[,c("clusterFam","complete_name")])$complete_name )

orgsVsOrgs <- sapply(orgsVsOrgsRaw, function(x){
  sapply(orgsVsOrgsRaw, function(y){
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    round(length(intersect(unique(x),
                           unique(y)))*2/sum(length(unique(x)),
                                             length(unique(y)) )  *100, digits = 2)
  })
})


write.table(orgsVsOrgs, sprintf("SMGC_Heatmap_orgVsOrg.tsv", fileName))


pdf(sprintf("orgsVsOrgs_%s.pdf", fileName), width = 16, height = 16, pointsize = 18)
heatmap.2(orgsVsOrgs,
          dendrogram = 'both',
          margins=c(12,9),
          #breaks = c(0,0.5,1),
          #na.color = 'white',
          sepcolor = "grey",
          colsep = 1:ncol(orgsVsOrgs),
          rowsep = 1:nrow(orgsVsOrgs),
          col = c("grey",rev(brewer.pal(9,"Spectral"))), #brewer.pal(10,"Oranges"),
          trace = 'none',
          xlab = NULL,
          # breaks = seq(0, 100, 20),
          breaks = seq(0, 100, 10),
          #  cexRow = 1,
          main = "Shared gene clusters in percent",
          # labRow = FALSE,
          # color key and info
          key = TRUE
)
dev.off()


#################
# Compound heatmap
###########


    compFamDf <- aggregate(compound ~ clusterFam ,
                           data = subset(sm, subset = !(is.na(compound))),
                           FUN = function(x){
                             base::paste(unique(x), collapse = ', ')
                           })

    compFamDf <- merge(
      unique(subset(sm, select = c("complete_name", "clusterFam"))),
      compFamDf,
      by = "clusterFam")



    predHm <- unclass(table(compFamDf[,c("complete_name", "compound")]))
    predHm[predHm > 0 ] = 1

     predHmSorted <- predHm[tipsComplete,]
     predHmSorted


     # remove metabolites which are only found in one organism
    compOnlyFoundOnce <- colnames(predHmSorted)[colSums(predHmSorted) == 1]


     pdf("compoundHeatmapAll.pdf", width = 10, height = 16, pointsize = 18)
     heatmap.2(predHmSorted[,!(colnames(predHmSorted) %in% compOnlyFoundOnce )],
               dendrogram = 'col',
               Rowv = FALSE,
               main = 'Dereplication of\nSMGC families by mibig',
               margins=c(12,9),
               na.color = 'white',
               col = c("white", brewer.pal(9,"Blues")[9]),
               trace = 'none',
               key = TRUE,
               #   tracecol = 'black',
               # linecol = NULL,
               sepcolor = "grey",
               colsep = 1:ncol(predHm),
               rowsep = 1:nrow(predHm)#,
     )
     dev.off()

write.table(compFamDf, "compoundHeatmapDf.tsv")





###############
# TREE COUNT of unique SMGC per branch
######


if(length(orgsMissingInDf) > 0){
  print("Warning, the following orgs are missing in the data frame:")
  print(orgsMissingInDf)
}

# if(length(orgsMissingInTree) >0){
#   stop("There are organisms from the dataframe missing in the tree:",
#        orgsMissingInTree)
# }

smTree <- unique(subset(sm,
       subset = sm_short != "none",
       select = c("name", "real_name", "complete_name")))


smTree <- rbind(smTree, data.frame(name = c("Neucr2", "Sacce1"), real_name = c("crassa","cerevisiae"), complete_name = c("N. crassa", "S. cerevisiae")))



# unique(sm$nodePosition)



famPerNode <- unique(subset(sm, subset = !is.na(nodePosition), select = c(clusterFam, nodePosition)))
famPerNode <- as.data.frame(table(famPerNode$nodePosition))
colnames(famPerNode) = c("node", "famCount")


p_select$data$famCount <- merge(data.frame(node =p_select$data$node), famPerNode, by = "node", all.x = TRUE)$famCount

psm <- p_select

# psm + geom_tiplab(aes(label = complete_name)) + geom_text2(aes(label = famCount), vjust = -.3, hjust = 1.5) + xlim(0,20)

p_select + geom_tiplab() + geom_text2(aes(label = famCount), vjust = -.3, hjust = 1.5) + xlim(0,20)



#######
# TREE WITH SM DATA
####

smOv <- unique(subset(sm,
                      subset = sm_short != "none",
                      select = c("name", "protein_id", "sm_short")))


smOv <- aggregate(cluster_id ~ name + sm_short,
                data = subset(sm, sm_short != 'none'),
                FUN = length)
names(smOv)[3] <- "sm_count"



orgsVsOrgsSorted <- orgsVsOrgs[tipsComplete,]


smOvRn <- aggregate(cluster_id ~ complete_name + sm_short,
                  data = subset(sm, sm_short != 'none'),
                  FUN = length)
names(smOvRn)[3] <- "sm_count"
smOvRn$complete_name <- factor(smOvRn$complete_name, levels=rev(tipsComplete))
smOvRn

ggplot(smOvRn)+ geom_barh(aes(y = complete_name, x = sm_count, fill = sm_short), color = 'black', stat = "identity")
length(tipsComplete)
tipsComplete[1:7]

ggsave("sm_bars.pdf", height = 11)




pdf("orgsVsOrgsMod.pdf", width = 12, height = 12, pointsize = 14)
heatmap.2(orgsVsOrgsSorted,
          dendrogram = 'col',
          Rowv = FALSE,
          Colv = TRUE,
          margins=c(12,10),
          #breaks = c(0,0.5,1),
          #na.color = 'white',
          sepcolor = "grey",
          colsep = 1:ncol(orgsVsOrgs),
          rowsep = 1:nrow(orgsVsOrgs),
          col = c("grey",rev(brewer.pal(9,"Spectral"))), #brewer.pal(10,"Oranges"),
          trace = 'none',
          # xlab = NULL,
          # breaks = seq(0, 100, 20),
          breaks = seq(0, 100, 10),
          #  cexRow = 1,
          #  main = "Shared gene clusters in percent"
          # labRow = FALSE,
          # color key and info
          key = TRUE#,
          # key.xlab = "Shared SMGC \n families (%)",
          # key.ylab = "Count",
          # keysize = 1.5
          # srtCol = 45
          # lwid=c(0.1,4),
          # lhei=c(0.1,1)
)
dev.off()
sm$complete_name
sm_selected_tree <- unique(subset(sm, select = c(name, complete_name) ))
p_select <- p_select %<+% sm_selected_tree


svg("treesMod.svg", width = 4, height = 10, pointsize = 14)
p_select +  geom_tiplab(aes(label = complete_name), align = TRUE) + xlim(0, .6)
dev.off()

p_select  + xlim(0, .2)+   theme(plot.background = element_blank(),
                                 panel.border = element_blank(),
                                 panel.background = element_blank()) + geom_treescale()
ggsave("treesMod.pdf", width = 4, height = 10, units = "in")


d <- p_select$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)


svg("treeBS.svg", width = 7, height = 7, pointsize = 14)
p_select + geom_text2(data = d, aes(label = label))  +  geom_tiplab(aes(label = complete_name), align = TRUE) + xlim(0, .25)
dev.off()


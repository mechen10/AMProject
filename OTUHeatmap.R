#!/bin/Rscript

# FOR AM PROJECT

library("optparse")
library("gplots")
########################### OPT PARSE #################################
option_list = list(
  make_option(c("-i", "--inputFolder"), type="character",
              help="Folder with reduced OTU tables that are only core OTUs for each treatment type"),
  make_option(c("-t", "--otuTable"), type="character",
              help="OTU Table input"),
  make_option(c("-m", "--metadata"), type="character",
              help="Metadata fp")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

inputFolder = opt$inputFolder
otuTableFP = opt$otuTableFP
metadataFP = opt$metadata

########################### FOR TESTING #################################
setwd('/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis')
inputFolder <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/ALLMORPH_cores'
otuTableFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/OTU_Table_text.txt'
metadataFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/TEMP_frombotclust/MF_nochlpmito_m1000.txt'

# setwd('/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis')
# inputFolder <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/CORE/ALLMORPH_cores'
# otuTableFP <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/OTU_Table_text.txt'
# metadataFP <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/OTU_MP_filt/MF_nochlpmito_m1000.txt'


########################### LOAD DATA #################################
system("mkdir OTUHEATMAP")

OTUTable <- read.delim(paste0(otuTableFP)
                       , header = TRUE
                       , skip = 1
                       , row.names = 1)
taxonomyNames <- as.data.frame(OTUTable[,ncol(OTUTable)])
OTUTable <- OTUTable[,-ncol(OTUTable)]

# Make OTU Table relative abundance
OTUTable.RelAbund <- OTUTable
colSumsOTUTable <- colSums(OTUTable)
for (i in 1:ncol(OTUTable)) {
  for (j in 1:nrow(OTUTable)) {
    OTUTable.RelAbund[j,i] <- OTUTable[j,i]/colSumsOTUTable[i]
  }
}


MF <- read.delim(paste0(metadataFP)
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE)
rownames(MF) <- gsub("-",".", rownames(MF))

ListOfFiles <- system(paste0('ls ',inputFolder), intern = TRUE)

P.Files <- ListOfFiles[grep("P_(BL|CR|FB|H2O1440|H2O20)", ListOfFiles)]
H.Files <- ListOfFiles[grep("H_(BL|CR|FB|W)", ListOfFiles)]

for (i in H.Files) {
  assign(i, read.delim(paste0(inputFolder,"/",i)
                       , stringsAsFactors = FALSE
                       , header = TRUE
                       , row.names = 1
                       , strip.white = TRUE))
}

for (i in P.Files) {
  assign(i, read.delim(paste0(inputFolder,"/",i)
                       , stringsAsFactors = FALSE
                       , header = TRUE
                       , row.names = 1
                       , strip.white = TRUE))
}
########################### START #################################

################ HAKAI ######################
# Get all OTUs for Hakai
allOTUs.H <- c()
for (i in H.Files) {
  allOTUs.H <- c(allOTUs.H, rownames(get(i)))
}
allOTUs.H <- unique(allOTUs.H)

# Get all samples from OTU Table
allSamples <- c()
for (i in H.Files) {
  allSamples <- c(allSamples, colnames(get(i)))
}

OTUTable.H.filt <- OTUTable.RelAbund[allOTUs.H,allSamples]
MF.H.filt <- MF[allSamples,]


OTUTemp <- aggregate(t(OTUTable.H.filt), by = list(MF.H.filt$ColRep), FUN = mean)
OTUTable.filt.col <- t(data.frame(OTUTemp, row.names = 1))
rownames(OTUTable.filt.col) <- rownames(OTUTable.H.filt)

# Now, order them

H.Files.time <- colnames(OTUTable.filt.col)[c(grep('(B|L|R|W)20', colnames(OTUTable.filt.col))
                                              , grep('(B|L|R|W)60', colnames(OTUTable.filt.col))
                                              , grep('(B|L|R|W)360', colnames(OTUTable.filt.col))
                                              , grep('(B|L|R|W)720', colnames(OTUTable.filt.col))
                                              , grep('(B|L|R|W)5760', colnames(OTUTable.filt.col))
                                            )]
H.Files.morphtime <- H.Files.time[c(grep('FB', H.Files.time)
                                    ,grep('BL', H.Files.time)
                                    ,grep('CR', H.Files.time)
                                    , grep('W', H.Files.time))]

OTUTable.filt.col <- OTUTable.filt.col[,H.Files.morphtime]


########## Make core/not core matrix############
OTUTable.core <- OTUTable.filt.col
for (i in 1:nrow(OTUTable.core)) {
  for (j in 1:ncol(OTUTable.core)) {
    FileName <- H.Files[grep(colnames(OTUTable.core)[j], H.Files)]
    if (rownames(OTUTable.core)[i] %in% rownames(get(FileName))) {
      OTUTable.core[i,j] <- 1
    } else {
      OTUTable.core[i,j] <- 0
      
    }
  }
}

########## Presence/Absence matrix ############
OTUTable.presabs <- OTUTable.filt.col
for (i in 1:nrow(OTUTable.presabs)) {
  for (j in 1:ncol(OTUTable.presabs)) {
    if (OTUTable.presabs[i,j] > 0) {
      OTUTable.presabs[i,j] <- 1
    } else {
      OTUTable.presabs[i,j] <- 0
      
    }
  }
}

################################# PLOT ####
colortemp <- colorRampPalette(c("white","darkgreen"))
xlabels <- c("20 m"
  , "1 h"
  , "6 h"
  , "12 h"
  , "4 d"
  , "20 m"
  , "1 h"
  , "6 h"
  , "12 h"
  , "4 d"
  , "20 m"
  , "1 h"
  , "6 h"
  , "12 h"
  , "4 d"
  , "6 h"
  , "4 d")


hc <- hclust(as.dist(1-cor(t(OTUTable.filt.col))))
hc.v <- hclust(as.dist(1-cor(OTUTable.filt.col)))

pdf("./OTUHEATMAP/OTUHeatmap_Hakai_scaled.pdf", pointsize = 14)
par(fig = c(0,1,0,1))
heatmap.2(as.matrix(OTUTable.filt.col)
          , Colv = NA
          , Rowv = as.dendrogram(hc)
          , cexCol = 1
          , labRow = NA
          , labCol = xlabels
          , col = colortemp(10)
          , margins = c(10,0.5)
          , density.info= "none"
          , trace = "none"
          # , breaks = c(0,0.00001,0.0001,0.001,0.01,0.1)
          , dendrogram = 'none'
          , ColSideColors = c(rep("salmon",5),rep("dodgerblue",5),rep("darkorchid",5), rep("darkblue",2))
          , scale = "row"
)
par(fig = c(0,1,0,1))
legend("topright"
       , legend = c("Finely branched","Bladed","Crustose","Water")
       , col = c("salmon","dodgerblue","darkorchid","darkblue")
       , pch = 19
       , cex = 0.7
       )
dev.off()
# 
# pdf("./OTUHEATMAP/OTUHeatmap_Hakai_scaledZERO.pdf", pointsize = 14)
# par(fig = c(0,1,0,1))
# heatmap.2(as.matrix(OTUTable.filt.col)
#           , Colv = NA
#           , Rowv = as.dendrogram(hc)
#           , cexCol = 1
#           , labRow = NA
#           , labCol = xlabels
#           , col = colortemp(10)
#           , margins = c(10,0.5)
#           , density.info= "none"
#           , trace = "none"
#           # , breaks = c(0,0.00001,0.0001,0.001,0.01,0.1)
#           , dendrogram = 'none'
#           , ColSideColors = c(rep("salmon",5),rep("dodgerblue",5),rep("darkorchid",5), rep("darkblue",2))
#           , scale = "row"
#           , breaks = seq(0,3,length.out = 11)
# )
# par(fig = c(0,1,0,1))
# legend("topright"
#        , legend = c("Finely branched","Bladed","Crustose","Water")
#        , col = c("salmon","dodgerblue","darkorchid","darkblue")
#        , pch = 19
#        , cex = 0.7
# )
# dev.off()

pdf("./OTUHEATMAP/OTUHeatmap_Hakai_coreorot.pdf", pointsize = 14)
heatmap.2(OTUTable.core
          , Colv = NA
          , Rowv = as.dendrogram(hc)
          , dendrogram = 'none'
          , cexCol = 1
          , labRow = NA
          , labCol = xlabels
          , col = colortemp(2)
          , margins = c(10,0.5)
          , ColSideColors = c(rep("salmon",5),rep("dodgerblue",5),rep("darkorchid",5), rep("darkblue",2))
          , scale = "none"
          , density.info= "none"
          , trace = "none"
)
par(fig = c(0,1,0,1))
legend("topright"
       , legend = c("Finely branched","Bladed","Crustose","Water")
       , col = c("salmon","dodgerblue","darkorchid","darkblue")
       , pch = 19
       , cex = 0.7
)
dev.off()

pdf("./OTUHEATMAP/OTUHeatmap_Hakai_presabs.pdf", pointsize = 14)
par(fig = c(0,1,0,1))
heatmap.2(OTUTable.presabs
          , Colv = NA
          , Rowv = as.dendrogram(hc)
          , cexCol = 1
          , labRow = NA
          , labCol = xlabels
          , col = colortemp(2)
          , margins = c(10,0.5)
          , density.info= "none"
          , trace = "none"
          # , breaks = c(0,0.00001,0.0001,0.001,0.01,0.1)
          , dendrogram = 'none'
          , ColSideColors = c(rep("salmon",5),rep("dodgerblue",5),rep("darkorchid",5), rep("darkblue",2))
          , scale = "none"
)
par(fig = c(0,1,0,1))
legend("topright"
       , legend = c("Finely branched","Bladed","Crustose","Water")
       , col = c("salmon","dodgerblue","darkorchid","darkblue")
       , pch = 19
       , cex = 0.7
)
dev.off()

# For looking at turnover
OTUTable.H.presabs <- OTUTable.presabs
################ PORTMOODY ######################
# Get all OTUs for PM
allOTUs.P <- c()
for (i in P.Files) {
  allOTUs.P <- c(allOTUs.P, rownames(get(i)))
}
allOTUs.P <- unique(allOTUs.P)

# Get all samples from OTU Table
allSamples <- c()
for (i in P.Files) {
  allSamples <- c(allSamples, colnames(get(i)))
}
# For P only-- get rid of X
allSamples <- gsub("X","", allSamples)
colnames(OTUTable.RelAbund) <- gsub("X","",colnames(OTUTable.RelAbund))

OTUTable.P.filt <- OTUTable.RelAbund[allOTUs.P,allSamples]
MF.P.filt <- MF[allSamples,]


OTUTemp <- aggregate(t(OTUTable.P.filt), by = list(MF.P.filt$ColRep), FUN = mean)
OTUTable.filt.col <- t(data.frame(OTUTemp, row.names = 1))
rownames(OTUTable.filt.col) <- rownames(OTUTable.P.filt)

# Now, order them

P.Files.time <- colnames(OTUTable.filt.col)[c(grep('(B|L|R|O)20', colnames(OTUTable.filt.col))
                                              , grep('(B|L|R|O)60', colnames(OTUTable.filt.col))
                                              , grep('(B|L|R|O)180', colnames(OTUTable.filt.col))
                                              , grep('(B|L|R|O)360', colnames(OTUTable.filt.col))
                                              , grep('(B|L|R|O)720', colnames(OTUTable.filt.col))
                                              , grep('(B|L|R|O)1440', colnames(OTUTable.filt.col))
)]
P.Files.morphtime <- P.Files.time[c(grep('FB', P.Files.time)
                                    ,grep('BL', P.Files.time)
                                    ,grep('CR', P.Files.time)
                                    , grep('H2O', P.Files.time))]

OTUTable.filt.col <- OTUTable.filt.col[,P.Files.morphtime]

########## Make core/not core matrix############
OTUTable.core <- OTUTable.filt.col
for (i in 1:nrow(OTUTable.core)) {
  for (j in 1:ncol(OTUTable.core)) {
    FileName <- P.Files[grep(colnames(OTUTable.core)[j], P.Files)]
    if (rownames(OTUTable.core)[i] %in% rownames(get(FileName))) {
      OTUTable.core[i,j] <- 1
    } else {
      OTUTable.core[i,j] <- 0
      
    }
  }
}

########## Presence/Absence matrix ############
OTUTable.presabs <- OTUTable.filt.col
for (i in 1:nrow(OTUTable.presabs)) {
  for (j in 1:ncol(OTUTable.presabs)) {
    if (OTUTable.presabs[i,j] > 0) {
      OTUTable.presabs[i,j] <- 1
    } else {
      OTUTable.presabs[i,j] <- 0
      
    }
  }
}


################################# PLOT ####
colortemp <- colorRampPalette(c("white","darkgreen"))
xlabels <- c("20 m"
             , "1 h"
             , "3 h"
             , "6 h"
             , "12 h"
             , "1 d"
             , "20 m"
             , "1 h"
             , "3 h"
             , "6 h"
             , "12 h"
             , "1 d"
             , "20 m"
             , "1 h"
             , "3 h"
             , "6 h"
             , "12 h"
             , "1 d"
             , "20 m"
             , "1 d")

hc <- hclust(as.dist(1-cor(t(OTUTable.filt.col))))
hc.v <- hclust(as.dist(1-cor(OTUTable.filt.col)))

pdf("./OTUHEATMAP/OTUHeatmap_PM_scaled.pdf", pointsize = 14)
par(fig = c(0,1,0,1))
heatmap.2(as.matrix(OTUTable.filt.col)
          , Colv = NA
          , Rowv = as.dendrogram(hc)
          , cexCol = 1
          , labRow = NA
          , labCol = xlabels
          , col = colortemp(10)
          , margins = c(10,0.5)
          , density.info= "none"
          , trace = "none"
          # , breaks = c(0,0.00001,0.0001,0.001,0.01,0.1)
          , dendrogram = 'none'
          , ColSideColors = c(rep("salmon",6),rep("dodgerblue",6),rep("darkorchid",6), rep("darkblue",2))
          , scale = "row"
)
par(fig = c(0,1,0,1))
legend("topright"
       , legend = c("Finely branched","Bladed","Crustose","Water")
       , col = c("salmon","dodgerblue","darkorchid","darkblue")
       , pch = 19
       , cex = 0.7
)
dev.off()

# 
# pdf("./OTUHEATMAP/OTUHeatmap_PM_scaledZERO.pdf", pointsize = 14)
# par(fig = c(0,1,0,1))
# heatmap.2(as.matrix(OTUTable.filt.col)
#           , Colv = NA
#           , Rowv = as.dendrogram(hc)
#           , cexCol = 1
#           , labRow = NA
#           , labCol = xlabels
#           , col = colortemp(10)
#           , margins = c(10,0.5)
#           , density.info= "none"
#           , trace = "none"
#           # , breaks = c(0,0.00001,0.0001,0.001,0.01,0.1)
#           , dendrogram = 'none'
#           , ColSideColors = c(rep("salmon",6),rep("dodgerblue",6),rep("darkorchid",6), rep("darkblue",2))
#           , scale = "row"
#           , breaks = seq(0,4,length.out = 11)
# )
# par(fig = c(0,1,0,1))
# legend("topright"
#        , legend = c("Finely branched","Bladed","Crustose","Water")
#        , col = c("salmon","dodgerblue","darkorchid","darkblue")
#        , pch = 19
#        , cex = 0.7
# )
# dev.off()

pdf("./OTUHEATMAP/OTUHeatmap_PM_coreornot.pdf", pointsize = 14)
par(fig = c(0,1,0,1))
heatmap.2(as.matrix(OTUTable.core)
          , Colv = NA
          , Rowv = as.dendrogram(hc)
          , cexCol = 1
          , labRow = NA
          , labCol = xlabels
          , col = colortemp(2)
          , margins = c(10,0.5)
          , density.info= "none"
          , trace = "none"
          # , breaks = c(0,0.00001,0.0001,0.001,0.01,0.1)
          , dendrogram = 'none'
          , ColSideColors = c(rep("salmon",6),rep("dodgerblue",6),rep("darkorchid",6), rep("darkblue",2))
          , scale = "none"
)
par(fig = c(0,1,0,1))
legend("topright"
       , legend = c("Finely branched","Bladed","Crustose","Water")
       , col = c("salmon","dodgerblue","darkorchid","darkblue")
       , pch = 19
       , cex = 0.7
)
dev.off()

pdf("./OTUHEATMAP/OTUHeatmap_PM_presabs.pdf", pointsize = 14)
par(fig = c(0,1,0,1))
heatmap.2(as.matrix(OTUTable.presabs)
          , Colv = NA
          , Rowv = as.dendrogram(hc)
          , cexCol = 1
          , labRow = NA
          , labCol = xlabels
          , col = colortemp(2)
          , margins = c(10,0.5)
          , density.info= "none"
          , trace = "none"
          # , breaks = c(0,0.00001,0.0001,0.001,0.01,0.1)
          , dendrogram = 'none'
          , ColSideColors = c(rep("salmon",6),rep("dodgerblue",6),rep("darkorchid",6), rep("darkblue",2))
          , scale = "none"
)
par(fig = c(0,1,0,1))
legend("topright"
       , legend = c("Finely branched","Bladed","Crustose","Water")
       , col = c("salmon","dodgerblue","darkorchid","darkblue")
       , pch = 19
       , cex = 0.7
)
dev.off()

OTUTable.P.presabs <- OTUTable.presabs
########### Finding presence/absence for turnover #############

# # HAKAI
# H.samples <- rownames(MF.H.filt)
# OTUTable.H <- OTUTable.RelAbund[allOTUs.H,H.samples]
# 
# OTUTable.H.presabs <- OTUTable.H
# for (r in 1:nrow(OTUTable.H.presabs)) {
#   for (c in 1:ncol(OTUTable.H.presabs)) {
#     if (OTUTable.H.presabs[r,c] > 0.001 ) {
#       OTUTable.H.presabs[r,c] <- 1
#     } else {
#       OTUTable.H.presabs[r,c] <- 0
#     }
#   }
# }
# #PM
# 
# P.samples <- rownames(MF.P.filt)
# OTUTable.P <- OTUTable.RelAbund[allOTUs.P,P.samples]
# 
# OTUTable.P.presabs <- OTUTable.P
# for (r in 1:nrow(OTUTable.P.presabs)) {
#   for (c in 1:ncol(OTUTable.P.presabs)) {
#     if (OTUTable.P.presabs[r,c] > 0.001 ) {
#       OTUTable.P.presabs[r,c] <- 1
#     } else {
#       OTUTable.P.presabs[r,c] <- 0
#     }
#   }
# }

#HAKAI

timesList.H <- c("20","60","360","720","5760")
turnOver.H <- matrix(nrow = 4, ncol = length(timesList.H))
rownames(turnOver.H) <- c("total","lost","gained", "retained")
colnames(turnOver.H) <- timesList.H

turnOver.H.FB <- matrix(nrow = 4, ncol = length(timesList.H))
rownames(turnOver.H.FB) <- c("total","lost","gained","retained")
colnames(turnOver.H.FB) <- timesList.H

turnOver.H.BL <- matrix(nrow = 4, ncol = length(timesList.H))
rownames(turnOver.H.BL) <- c("total","lost","gained","retained")
colnames(turnOver.H.BL) <- timesList.H

turnOver.H.CR <- matrix(nrow = 4, ncol = length(timesList.H))
rownames(turnOver.H.CR) <- c("total","lost","gained","retained")
colnames(turnOver.H.CR) <- timesList.H

allold <- c()
FBold <- c()
BLold <- c()
CRold <- c()
i <- 2
for (i in 1:length(timesList.H)) {
  ttemp <- timesList.H[i]
  FBonly <- OTUTable.H.presabs[,grep(paste0("FB", ttemp), colnames(OTUTable.H.presabs))]
  BLonly <- OTUTable.H.presabs[,grep(paste0("BL", ttemp), colnames(OTUTable.H.presabs))]
  CRonly <- OTUTable.H.presabs[,grep(paste0("CR", ttemp), colnames(OTUTable.H.presabs))]
  allonly <- OTUTable.H.presabs[,grep(paste0("(FB|BL|CR)[.]", ttemp,"[.]"), colnames(OTUTable.H.presabs))]
  
  # Filter by which ones are actually there
  FBonly <- FBonly[which(FBonly != 0)]
  BLonly <- BLonly[which(BLonly != 0)]
  CRonly <- CRonly[which(CRonly != 0)]
  allonly <- allonly[which(rowSums(allonly) != 0),]
  
  # Get total OTU count for this time
  FBcount <- length(FBonly)
  BLcount <- length(BLonly)
  CRcount <- length(CRonly)
  allcount <- nrow(allonly)
  
  # Get "old" list of OTUs and compare which ones are the same
  FBlost <- length(FBold)-sum(FBold %in% names(FBonly))
  BLlost <- length(BLold)-sum(BLold %in% names(BLonly))
  CRlost <- length(CRold)-sum(CRold %in% names(CRonly))
  alllost <- length(allold)-sum(allold %in% rownames(allonly))
  
  FBgain <- FBcount-sum(names(FBonly) %in% FBold)
  BLgain <- BLcount-sum(names(BLonly) %in% BLold)
  CRgain <- CRcount-sum(names(CRonly) %in% CRold)
  allgain <- allcount-sum(rownames(allonly) %in% allold)
  
  # Get "retained" OTUs between each
  
  FBretain <- sum(names(FBonly) %in% FBold)
  BLretain <- sum(names(BLonly) %in% BLold)
  CRretain <- sum(names(CRonly) %in% CRold)
  allretain <- sum(rownames(allonly) %in% allold)
  
  
  # Load into matrices
  turnOver.H.FB["total", paste0(ttemp)] <- FBcount
  turnOver.H.FB["lost", paste0(ttemp)] <- round(FBlost/FBcount,2)
  turnOver.H.FB["gained", paste0(ttemp)] <- round(FBgain/FBcount,2)
  turnOver.H.FB["retained", paste0(ttemp)] <- round(FBretain/FBcount,2)
  
  turnOver.H.BL["total", paste0(ttemp)] <- BLcount
  turnOver.H.BL["lost", paste0(ttemp)] <- round(BLlost/BLcount,2)
  turnOver.H.BL["gained", paste0(ttemp)] <- round(BLgain/BLcount,2)
  turnOver.H.BL["retained", paste0(ttemp)] <- round(BLretain/BLcount,2)
  
  turnOver.H.CR["total", paste0(ttemp)] <- CRcount
  turnOver.H.CR["lost", paste0(ttemp)] <- round(CRlost/CRcount,2)
  turnOver.H.CR["gained", paste0(ttemp)] <- round(CRgain/CRcount,2)
  turnOver.H.CR["retained", paste0(ttemp)] <- round(CRretain/CRcount,2)
  
  turnOver.H["total", paste0(ttemp)] <- allcount
  turnOver.H["lost", paste0(ttemp)] <- round(alllost/allcount,2)
  turnOver.H["gained", paste0(ttemp)] <- round(allgain/allcount,2)
  turnOver.H["retained", paste0(ttemp)] <- round(allretain/allcount,2)
  
  # Finally, set "old" as the new ones

  FBold <- names(which(FBonly != 0))
  BLold <- names(which(BLonly != 0))
  CRold <- names(which(CRonly != 0))
  allold <- names(which(rowSums(allonly) != 0))
}

# PM


timesList.P <- c("20","60","180","360","720","1440")
turnOver.P <- matrix(nrow = 4, ncol = length(timesList.P))
rownames(turnOver.P) <- c("total","lost","gained", "retained")
colnames(turnOver.P) <- timesList.P

turnOver.P.FB <- matrix(nrow = 4, ncol = length(timesList.P))
rownames(turnOver.P.FB) <- c("total","lost","gained","retained")
colnames(turnOver.P.FB) <- timesList.P

turnOver.P.BL <- matrix(nrow = 4, ncol = length(timesList.P))
rownames(turnOver.P.BL) <- c("total","lost","gained","retained")
colnames(turnOver.P.BL) <- timesList.P

turnOver.P.CR <- matrix(nrow = 4, ncol = length(timesList.P))
rownames(turnOver.P.CR) <- c("total","lost","gained","retained")
colnames(turnOver.P.CR) <- timesList.P

allold <- c()
FBold <- c()
BLold <- c()
CRold <- c()
for (i in 1:length(timesList.P)) {
  ttemp <- timesList.P[i]
  FBonly <- OTUTable.P.presabs[,grep(paste0("FB",ttemp), colnames(OTUTable.P.presabs))]
  BLonly <- OTUTable.P.presabs[,grep(paste0("BL",ttemp), colnames(OTUTable.P.presabs))]
  CRonly <- OTUTable.P.presabs[,grep(paste0("CR",ttemp), colnames(OTUTable.P.presabs))]
  allonly <- OTUTable.P.presabs[,grep(paste0("(FB|BL|CR)",ttemp), colnames(OTUTable.P.presabs))]
  
  # Filter by which ones are actually there
  FBonly <- FBonly[which(FBonly != 0)]
  BLonly <- BLonly[which(BLonly != 0)]
  CRonly <- CRonly[which(CRonly != 0)]
  allonly <- allonly[which(rowSums(allonly) != 0),]
  
  # Get total OTU count for this time
  FBcount <- length(FBonly)
  BLcount <- length(BLonly)
  CRcount <- length(CRonly)
  allcount <- nrow(allonly)
  
  # Get "old" list of OTUs and compare which ones are the same
  FBlost <- length(FBold)-sum(FBold %in% names(FBonly))
  BLlost <- length(BLold)-sum(BLold %in% names(BLonly))
  CRlost <- length(CRold)-sum(CRold %in% names(CRonly))
  alllost <- length(allold)-sum(allold %in% rownames(allonly))
  
  FBgain <- FBcount-sum(names(FBonly) %in% FBold)
  BLgain <- BLcount-sum(names(BLonly) %in% BLold)
  CRgain <- CRcount-sum(names(CRonly) %in% CRold)
  allgain <- allcount-sum(rownames(allonly) %in% allold)
  
  # Get "retained" OTUs between each
  
  FBretain <- sum(names(FBonly) %in% FBold)
  BLretain <- sum(names(BLonly) %in% BLold)
  CRretain <- sum(names(CRonly) %in% CRold)
  allretain <- sum(rownames(allonly) %in% allold)
  
  
  # Load into matrices
  turnOver.P.FB["total", paste0(ttemp)] <- FBcount
  turnOver.P.FB["lost", paste0(ttemp)] <- round(FBlost/FBcount,2)
  turnOver.P.FB["gained", paste0(ttemp)] <- round(FBgain/FBcount,2)
  turnOver.P.FB["retained", paste0(ttemp)] <- round(FBretain/FBcount,2)
  
  turnOver.P.BL["total", paste0(ttemp)] <- BLcount
  turnOver.P.BL["lost", paste0(ttemp)] <- round(BLlost/BLcount,2)
  turnOver.P.BL["gained", paste0(ttemp)] <- round(BLgain/BLcount,2)
  turnOver.P.BL["retained", paste0(ttemp)] <- round(BLretain/BLcount,2)
  
  turnOver.P.CR["total", paste0(ttemp)] <- CRcount
  turnOver.P.CR["lost", paste0(ttemp)] <- round(CRlost/CRcount,2)
  turnOver.P.CR["gained", paste0(ttemp)] <- round(CRgain/CRcount,2)
  turnOver.P.CR["retained", paste0(ttemp)] <- round(CRretain/CRcount,2)
  
  turnOver.P["total", paste0(ttemp)] <- allcount
  turnOver.P["lost", paste0(ttemp)] <- round(alllost/allcount,2)
  turnOver.P["gained", paste0(ttemp)] <- round(allgain/allcount,2)
  turnOver.P["retained", paste0(ttemp)] <- round(allretain/allcount,2)
  
  # Finally, set "old" as the new ones
  
  FBold <- names(which(FBonly != 0))
  BLold <- names(which(BLonly != 0))
  CRold <- names(which(CRonly != 0))
  allold <- names(which(rowSums(allonly) != 0))
}

retained <- matrix(nrow = 6, ncol = 6)
rownames(retained) <- c("H.FB","H.BL","H.CR","P.FB","P.BL","P.CR")
colnames(retained) <- c("60","180","360","720","1440","5760")
for (m in c("FB","CR","BL")) {
  for (t in c("60","180","360","720","1440")) {
    TempP <- get(paste0("turnOver.P.",m))["retained",t]
    retained[paste0("P.",m), t] <- TempP
  }
  for (t in c("60","360","720","5760")) {
    TempH <- get(paste0("turnOver.H.",m))["retained",t]
    retained[paste0("H.",m), t] <- TempH
  }
}




pdf("OTUHEATMAP/retainedOTUs.pdf", pointsize = 14, width = 7, height = 5)
par(fig = c(0,0.8,0,1), mar = c(6,4,2,4))
plot(c(1,2,3,4,5,6),NULL
     , xlim = c(1,6)
     , ylim = c(0,0.35)
     , pch= ""
     , xaxt = "n"
     , xlab = ""
     , ylab = "Percent community turnover"
)
title(xlab = "Time"
      , line = 4)
axis(side = 1
     , at = c(1,2,3,4,5,6)
     , labels = c("1 hour","3 hours","6 hours","12 hours","1 day","4 days")
     , las = 2)
points(1-retained[1,] ~ c(1,2,3,4,5,6)
       , col = "salmon"
       , pch = 19)
points(1-retained[2,] ~ c(1,2,3,4,5,6)
        , col = "dodgerblue"
       , pch = 19)
points(1-retained[3,] ~ c(1,2,3,4,5,6)
        , col = "darkorchid4"
       , pch = 19)
points(1-retained[4,] ~ c(1,2,3,4,5,6)
        , col = "salmon"
       , pch = 21)
points(1-retained[5,] ~ c(1,2,3,4,5,6)
        , col = "dodgerblue"
       , pch = 21)
points(1-retained[6,] ~ c(1,2,3,4,5,6)
        , col = "darkorchid4"
       , pch = 21)
lines(1-na.exclude(retained[1,]) ~ c(1,2,3,4,5,6)[!is.na(retained[1,])] 
      , col = "salmon"
      , lty = "solid"
      , lwd = 2
      )
lines(1-na.exclude(retained[2,]) ~ c(1,2,3,4,5,6)[!is.na(retained[2,])] 
      , col = "dodgerblue"
      , lty = "solid"
      , lwd = 2
      )
lines(1-na.exclude(retained[3,]) ~ c(1,2,3,4,5,6)[!is.na(retained[3,])] 
      , col = "darkorchid4"
      , lty = "solid"
      , lwd = 2
      )
lines(1-na.exclude(retained[4,]) ~ c(1,2,3,4,5,6)[!is.na(retained[4,])] 
      , col = "salmon"
      , lty = "dotted"
      , lwd = 2
      )
lines(1-na.exclude(retained[5,]) ~ c(1,2,3,4,5,6)[!is.na(retained[5,])] 
      , col = "dodgerblue"
      , lty = "dotted"
      , lwd = 2
      )
lines(1-na.exclude(retained[6,]) ~ c(1,2,3,4,5,6)[!is.na(retained[6,])] 
      , col = "darkorchid4"
      , lty = "dotted"
      , lwd = 2
      )
par(fig = c(0.7,1,0,1), mar = c(4,0,4,0), new = TRUE)
plot(0, NULL
     , pch = ""
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = ""
     , ylab = ""
     , bty = 'n'
     )
legend( "topleft"
        , legend = c("FB","BL","CR")
        , pch = 22
        , pt.bg = c("salmon","dodgerblue","darkorchid4")
        , col = c("salmon","dodgerblue","darkorchid4"))
legend( "left"
        , legend = c("Reed Point","Hakai")
        , lty = c("dotted","solid")
        , lwd = c(2,2)
        , pch = c(21,19)
        )
dev.off()

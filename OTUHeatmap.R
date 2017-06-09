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
# setwd('/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis')
# inputFolder <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/ALLMORPH_cores'
# otuTableFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/OTU_Table_text.txt'
# metadataFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/TEMP_frombotclust/MF_nochlpmito_m1000.txt'

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

#!/bin/Rscript

library("optparse")
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

P.Files <- ListOfFiles[grep("P_(BL|CR|FB|H201440|H2O20)", ListOfFiles)]
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

# Reorder OTUs


# Make water core first
water.Files <- H.Files[grep("W", H.Files)]
coreOTUs.W <- c()
for (i in water.Files) {
  coreOTUs.W <- c(coreOTUs.W, rownames(get(i)))
}
coreOTUs.W <- unique(coreOTUs.W)

# Core across all shapes?
AM.Files <- H.Files[grep("CR|BL|FB", H.Files)][-grep("5760", H.Files)]
maxcount <- length(AM.Files)
coreOTUs.morphs <- c()
for (i in AM.Files) {
  coreOTUs.morphs <- c(coreOTUs.morphs, rownames(get(i)))
}
coreOTUs.morphs.all <- c()
for (i in 1:length(table(coreOTUs.morphs)/maxcount)) {
  if ((table(coreOTUs.morphs)/maxcount)[i] > 0.9) {
    coreOTUs.morphs.all <- c(coreOTUs.morphs.all, rownames(table(coreOTUs.morphs))[i])
  }
}

# Core across each shape
# FB
FB.Files <- H.Files[grep("FB", H.Files)]
maxcount <- length(FB.Files)
coreOTUs.FB <- c()
for (i in FB.Files) {
  coreOTUs.FB <- c(coreOTUs.FB, rownames(get(i)))
}
coreOTUs.FB.all <- c()
coreOTUs.FB.leftover <- c()
for (i in 1:length(table(coreOTUs.FB)/maxcount)) {
  if ((table(coreOTUs.FB)/maxcount)[i] > 0.9) {
    coreOTUs.FB.all <- c(coreOTUs.FB.all, rownames(table(coreOTUs.FB))[i])
  }
}
# FB - 5760 without
FB.Files.nolast <- FB.Files[-grep("5760", FB.Files)]
maxcount <- length(FB.Files.nolast)
coreOTUs.FB.nolast <- c()
for (i in FB.Files.nolast) {
  coreOTUs.FB.nolast <- c(coreOTUs.FB.nolast, rownames(get(i)))
}
coreOTUs.FB.all.nolast <- c()
coreOTUs.FB.nolast.leftover <- c()
for (i in 1:length(table(coreOTUs.FB.nolast)/maxcount)) {
  if ((table(coreOTUs.FB.nolast)/maxcount)[i] > 0.9) {
    coreOTUs.FB.all.nolast <- c(coreOTUs.FB.all.nolast, rownames(table(coreOTUs.FB.nolast))[i])
  } 
}
FB.unique.final <- unique(c(coreOTUs.FB.all, coreOTUs.FB.all.nolast, coreOTUs.FB))

# BL
BL.Files <- H.Files[grep("BL", H.Files)]
maxcount <- length(BL.Files)
coreOTUs.BL <- c()
for (i in BL.Files) {
  coreOTUs.BL <- c(coreOTUs.BL, rownames(get(i)))
}
coreOTUs.BL.all <- c()
coreOTUs.BL.leftover <- c()
for (i in 1:length(table(coreOTUs.BL)/maxcount)) {
  if ((table(coreOTUs.BL)/maxcount)[i] > 0.9) {
    coreOTUs.BL.all <- c(coreOTUs.BL.all, rownames(table(coreOTUs.BL))[i])
  }
}
# BL - 5760 without
BL.Files.nolast <- BL.Files[-grep("5760", BL.Files)]
maxcount <- length(BL.Files.nolast)
coreOTUs.BL.nolast <- c()
for (i in BL.Files.nolast) {
  coreOTUs.BL.nolast <- c(coreOTUs.BL.nolast, rownames(get(i)))
}
coreOTUs.BL.all.nolast <- c()
coreOTUs.BL.nolast.leftover <- c()
for (i in 1:length(table(coreOTUs.BL.nolast)/maxcount)) {
  if ((table(coreOTUs.BL.nolast)/maxcount)[i] > 0.9) {
    coreOTUs.BL.all.nolast <- c(coreOTUs.BL.all.nolast, rownames(table(coreOTUs.BL.nolast))[i])
  } 
}
BL.unique.final <- unique(c(coreOTUs.BL.all, coreOTUs.BL.all.nolast, coreOTUs.BL))



# CR
CR.Files <- H.Files[grep("CR", H.Files)]
maxcount <- length(CR.Files)
coreOTUs.CR <- c()
for (i in CR.Files) {
  coreOTUs.CR <- c(coreOTUs.CR, rownames(get(i)))
}
coreOTUs.CR.all <- c()
coreOTUs.CR.leftover <- c()
for (i in 1:length(table(coreOTUs.CR)/maxcount)) {
  if ((table(coreOTUs.CR)/maxcount)[i] > 0.9) {
    coreOTUs.CR.all <- c(coreOTUs.CR.all, rownames(table(coreOTUs.CR))[i])
  }
}
# CR - 5760 without
CR.Files.nolast <- CR.Files[-grep("5760", CR.Files)]
maxcount <- length(CR.Files.nolast)
coreOTUs.CR.nolast <- c()
for (i in CR.Files.nolast) {
  coreOTUs.CR.nolast <- c(coreOTUs.CR.nolast, rownames(get(i)))
}
coreOTUs.CR.all.nolast <- c()
coreOTUs.CR.nolast.leftover <- c()
for (i in 1:length(table(coreOTUs.CR.nolast)/maxcount)) {
    if ((table(coreOTUs.CR.nolast)/maxcount)[i] > 0.9) {
      coreOTUs.CR.all.nolast <- c(coreOTUs.CR.all.nolast, rownames(table(coreOTUs.CR.nolast))[i])
    } 
}
CR.unique.final <- unique(c(coreOTUs.CR.all, coreOTUs.CR.all.nolast, coreOTUs.CR))

# Only 5760
H.Files.lastonly <- H.Files[grep("(BL|FB|CR)5760", H.Files)]
maxcount <- length(H.Files.lastonly)
coreOTUs.lastonly <- c()
for (i in H.Files.lastonly) {
  coreOTUs.lastonly <- c(coreOTUs.lastonly, rownames(get(i)))
}
coreOTUs.lastonly.all <- c()
for (i in 1:length(table(coreOTUs.lastonly)/maxcount)) {
  if ((table(coreOTUs.lastonly)/maxcount)[i] > 0.9) {
    coreOTUs.lastonly.all <- c(coreOTUs.lastonly.all, rownames(table(coreOTUs.lastonly))[i])
  }
}


# Now, sort OTUTable.filt.col OTUs
ORDERED.otus <- unique(c(
      coreOTUs.lastonly
     , coreOTUs.morphs.all
     , CR.unique.final
     , BL.unique.final
     , FB.unique.final
     , coreOTUs.W 
))

OTUTable.filt.col <- OTUTable.filt.col[ORDERED.otus,]

################################# PLOT ####
colortemp <- colorRampPalette(c("white","darkgreen"))
xlabels <- c("Finely Branched: 20 m"
  , "Finely Branched: 1 h"
  , "Finely Branched: 6 h"
  , "Finely Branched: 12 h"
  , "Finely Branched: 4 d"
  , "Bladed: 20 m"
  , "Bladed: 1 h"
  , "Bladed: 6 h"
  , "Bladed: 12 h"
  , "Bladed: 4 d"
  , "Crustose: 20 m"
  , "Crustose: 1 h"
  , "Crustose: 6 h"
  , "Crustose: 12 h"
  , "Crustose: 4 d"
  , "Water: 6 h"
  , "Water: 4 d")

pdf(paste0("./OTUHEATMAP/OTUHeatmap_","Hakai",".pdf"), pointsize = 14)
heatmap(as.matrix(OTUTable.filt.col)
        , Colv = NA
        , Rowv = NA
        , cexCol = 0.7
        , labRow = NA
        , labCol = xlabels
        , col = colortemp(10)
        , margins = c(10,0.5)
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

# Reorder OTUs


# Make water core first
water.Files <- P.Files[grep("H2O", P.Files)]
coreOTUs.W <- c()
for (i in water.Files) {
  coreOTUs.W <- c(coreOTUs.W, rownames(get(i)))
}
coreOTUs.W <- unique(coreOTUs.W)

# Core across all shapes?
AM.Files <- P.Files[grep("CR|BL|FB", P.Files)]
maxcount <- length(AM.Files)
coreOTUs.morphs <- c()
for (i in AM.Files) {
  coreOTUs.morphs <- c(coreOTUs.morphs, rownames(get(i)))
}
coreOTUs.morphs.all <- c()
for (i in 1:length(table(coreOTUs.morphs)/maxcount)) {
  if ((table(coreOTUs.morphs)/maxcount)[i] > 0.9) {
    coreOTUs.morphs.all <- c(coreOTUs.morphs.all, rownames(table(coreOTUs.morphs))[i])
  }
}

# Core across each shape
# FB
FB.Files <- P.Files[grep("FB", P.Files)]
maxcount <- length(FB.Files)
coreOTUs.FB <- c()
for (i in FB.Files) {
  coreOTUs.FB <- c(coreOTUs.FB, rownames(get(i)))
}
coreOTUs.FB.all <- c()
coreOTUs.FB.leftover <- c()
for (i in 1:length(table(coreOTUs.FB)/maxcount)) {
  if ((table(coreOTUs.FB)/maxcount)[i] > 0.9) {
    coreOTUs.FB.all <- c(coreOTUs.FB.all, rownames(table(coreOTUs.FB))[i])
  }
}
# # FB - 5760 without
# FB.Files.nolast <- FB.Files[-grep("5760", FB.Files)]
# maxcount <- length(FB.Files.nolast)
# coreOTUs.FB.nolast <- c()
# for (i in FB.Files.nolast) {
#   coreOTUs.FB.nolast <- c(coreOTUs.FB.nolast, rownames(get(i)))
# }
# coreOTUs.FB.all.nolast <- c()
# coreOTUs.FB.nolast.leftover <- c()
# for (i in 1:length(table(coreOTUs.FB.nolast)/maxcount)) {
#   if ((table(coreOTUs.FB.nolast)/maxcount)[i] > 0.9) {
#     coreOTUs.FB.all.nolast <- c(coreOTUs.FB.all.nolast, rownames(table(coreOTUs.FB.nolast))[i])
#   } 
# }
FB.unique.final <- unique(c(coreOTUs.FB.all, coreOTUs.FB))

# BL
BL.Files <- P.Files[grep("BL", P.Files)]
maxcount <- length(BL.Files)
coreOTUs.BL <- c()
for (i in BL.Files) {
  coreOTUs.BL <- c(coreOTUs.BL, rownames(get(i)))
}
coreOTUs.BL.all <- c()
coreOTUs.BL.leftover <- c()
for (i in 1:length(table(coreOTUs.BL)/maxcount)) {
  if ((table(coreOTUs.BL)/maxcount)[i] > 0.9) {
    coreOTUs.BL.all <- c(coreOTUs.BL.all, rownames(table(coreOTUs.BL))[i])
  }
}
# # BL - 5760 without
# BL.Files.nolast <- BL.Files[-grep("5760", BL.Files)]
# maxcount <- length(BL.Files.nolast)
# coreOTUs.BL.nolast <- c()
# for (i in BL.Files.nolast) {
#   coreOTUs.BL.nolast <- c(coreOTUs.BL.nolast, rownames(get(i)))
# }
# coreOTUs.BL.all.nolast <- c()
# coreOTUs.BL.nolast.leftover <- c()
# for (i in 1:length(table(coreOTUs.BL.nolast)/maxcount)) {
#   if ((table(coreOTUs.BL.nolast)/maxcount)[i] > 0.9) {
#     coreOTUs.BL.all.nolast <- c(coreOTUs.BL.all.nolast, rownames(table(coreOTUs.BL.nolast))[i])
#   } 
# }
BL.unique.final <- unique(c(coreOTUs.BL.all, coreOTUs.BL))



# CR
CR.Files <- P.Files[grep("CR", P.Files)]
maxcount <- length(CR.Files)
coreOTUs.CR <- c()
for (i in CR.Files) {
  coreOTUs.CR <- c(coreOTUs.CR, rownames(get(i)))
}
coreOTUs.CR.all <- c()
coreOTUs.CR.leftover <- c()
for (i in 1:length(table(coreOTUs.CR)/maxcount)) {
  if ((table(coreOTUs.CR)/maxcount)[i] > 0.9) {
    coreOTUs.CR.all <- c(coreOTUs.CR.all, rownames(table(coreOTUs.CR))[i])
  }
}
# # CR - 5760 without
# CR.Files.nolast <- CR.Files[-grep("5760", CR.Files)]
# maxcount <- length(CR.Files.nolast)
# coreOTUs.CR.nolast <- c()
# for (i in CR.Files.nolast) {
#   coreOTUs.CR.nolast <- c(coreOTUs.CR.nolast, rownames(get(i)))
# }
# coreOTUs.CR.all.nolast <- c()
# coreOTUs.CR.nolast.leftover <- c()
# for (i in 1:length(table(coreOTUs.CR.nolast)/maxcount)) {
#   if ((table(coreOTUs.CR.nolast)/maxcount)[i] > 0.9) {
#     coreOTUs.CR.all.nolast <- c(coreOTUs.CR.all.nolast, rownames(table(coreOTUs.CR.nolast))[i])
#   } 
# }
CR.unique.final <- unique(c(coreOTUs.CR.all, coreOTUs.CR))

# # Only 5760
# P.Files.lastonly <- P.Files[grep("(BL|FB|CR)5760", P.Files)]
# maxcount <- length(P.Files.lastonly)
# coreOTUs.lastonly <- c()
# for (i in P.Files.lastonly) {
#   coreOTUs.lastonly <- c(coreOTUs.lastonly, rownames(get(i)))
# }
# coreOTUs.lastonly.all <- c()
# for (i in 1:length(table(coreOTUs.lastonly)/maxcount)) {
#   if ((table(coreOTUs.lastonly)/maxcount)[i] > 0.9) {
#     coreOTUs.lastonly.all <- c(coreOTUs.lastonly.all, rownames(table(coreOTUs.lastonly))[i])
#   }
# }


# Now, sort OTUTable.filt.col OTUs
ORDERED.otus <- unique(c(
  # coreOTUs.lastonly
   coreOTUs.morphs.all
  , CR.unique.final
  , BL.unique.final
  , FB.unique.final
  , coreOTUs.W 
))

OTUTable.filt.col <- OTUTable.filt.col[ORDERED.otus,]

################################# PLOT ####
colortemp <- colorRampPalette(c("white","darkgreen"))
xlabels <- c("Finely Branched: 20 m"
             , "Finely Branched: 1 h"
             , "Finely Branched: 3 h"
             , "Finely Branched: 6 h"
             , "Finely Branched: 12 h"
             , "Finely Branched: 1 d"
             , "Bladed: 20 m"
             , "Bladed: 1 h"
             , "Bladed: 3 h"
             , "Bladed: 6 h"
             , "Bladed: 12 h"
             , "Bladed: 1 d"
             , "Crustose: 20 m"
             , "Crustose: 1 h"
             , "Crustose: 3 h"
             , "Crustose: 6 h"
             , "Crustose: 12 h"
             , "Crustose: 1 d"
             , "Water: 20 m"
             , "Water: 1 d")

pdf(paste0("./OTUHEATMAP/OTUHeatmap_","PM",".pdf"), pointsize = 14)
heatmap(as.matrix(OTUTable.filt.col)
        , Colv = NA
        , Rowv = NA
        , cexCol = 0.7
        , labRow = NA
        , labCol = xlabels
        , col = colortemp(10)
        , margins = c(10,0.5)
)
dev.off()

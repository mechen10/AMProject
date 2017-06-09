#!/bin/Rscript

# FOR AM PROJECT

library("optparse")
########################### OPT PARSE #################################
option_list = list(
  make_option(c("-o", "--otutable"), type="character",
              help="Full OTU Table"),
  make_option(c("-c", "--condensedOTUTable"),
              help="OTU Table at desired taxa level to collapse", type="character"),
  make_option(c("-m", "--mappingfile"),
              help="Mappingfile", type="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

OTUFP = opt$otutable
condensedFP = opt$condensedOTUTable
MFPWD = opt$mappingfile
########################### FOR TESTING #################################
# 
# setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/")
# OTUFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/OTU_Table_text.txt'
# condensedFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/OTU_L4.txt'
# keepUnColFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/TOKEEP.txt'
# MFPWD <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/MF_nochlpmito_m1000.txt'


# setwd("/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis")
# OTUFP <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/OTU_Table_text.txt'
# condensedFP <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/summarize_taxa/OTU_Table_nochlpmito_m1000_L4.txt'
# MFPWD <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/OTU_MP_filt/MF_nochlpmito_m1000.txt'
# 

########################### LOADING #################################
system("mkdir TAXASUMMARIES")
# Make Taxa summaries

OTUTable <- read.delim(paste0(OTUFP)
                       , header = TRUE
                       , skip = 1
                       , row.names = 1)
taxonomyNames <- as.data.frame(OTUTable[,ncol(OTUTable)])
rownames(taxonomyNames) <- rownames(OTUTable)

OTUTable.3 <- read.delim(paste0(condensedFP)
                       , header = TRUE
                       , skip = 1
                       , row.names = 1)

OTUTable <- OTUTable[,-ncol(OTUTable)]
# OTUTable.3 <- OTUTable.3[,-ncol(OTUTable.3)]

colnames(OTUTable) <- gsub('X', '',colnames(OTUTable))
colnames(OTUTable.3) <- gsub('X','',colnames(OTUTable.3))

# Get rid of FB.20.6 and FB.20.9
OTUTable <- OTUTable[,-grep("FB.20.6", colnames(OTUTable))]
OTUTable <- OTUTable[,-grep("FB.20.9", colnames(OTUTable))]

OTUTable.3 <- OTUTable.3[,-grep("FB.20.6", colnames(OTUTable.3))]
OTUTable.3 <- OTUTable.3[,-grep("FB.20.9", colnames(OTUTable.3))]


MF <- read.delim(paste0(MFPWD)
                       , header = TRUE
                       , row.names = 1
                       , stringsAsFactors = FALSE)
rownames(MF) <- gsub("-",".", rownames(MF))

# Make sure they are mutually in each other

MF <- MF[rownames(MF) %in% colnames(OTUTable),]
OTUTable <- OTUTable[,colnames(OTUTable) %in% rownames(MF)]
OTUTable.3 <- OTUTable.3[,colnames(OTUTable.3) %in% rownames(MF)]

# Split hakai and non hakai samples
# Get all samples I want; remove the rest
MF.inclWater <- MF
MF.P.inclWater <- MF
OTUTable.inclWater <- OTUTable
OTUTable.P.inclWater <- OTUTable
OTUTable.3.inclWater <- OTUTable.3
OTUTable.3.P.inclWater <- OTUTable.3
for (ROW in 1:nrow(MF)) {
  typeMorph <- MF[ROW, 'Morph']
  typeType <- MF[ROW, 'Type']
  rownameTemp <- rownames(MF)[ROW]
  if ((!(typeType == "H")) | (!(typeMorph %in% c('CR',"BL","FB","W","H2O")))) {
      MF.inclWater <- MF.inclWater[-grep(paste0("^",rownameTemp,"$"), rownames(MF.inclWater)),]
      OTUTable.inclWater <- OTUTable.inclWater[,-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(OTUTable.inclWater))]
      OTUTable.3.inclWater <- OTUTable.3.inclWater[,-grep(paste0("^", rownames(MF)[ROW], "$"), colnames(OTUTable.3.inclWater))]
  } 
  if ((!(typeType == "P")) | (!(typeMorph %in% c('CR',"BL","FB","W","H2O")))) {
    MF.P.inclWater <- MF.P.inclWater[-grep(paste0("^",rownameTemp,"$"), rownames(MF.P.inclWater)),]
    OTUTable.P.inclWater <- OTUTable.P.inclWater[,-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(OTUTable.P.inclWater))]
    OTUTable.3.P.inclWater <- OTUTable.3.P.inclWater[,-grep(paste0("^", rownames(MF)[ROW], "$"), colnames(OTUTable.3.P.inclWater))]
  }  
}



# Now make relative abundance

OTUTable.RelAbund <- OTUTable.inclWater
colSumsOTUTable <- colSums(OTUTable.inclWater)
for (i in 1:ncol(OTUTable.inclWater)) {
  for (j in 1:nrow(OTUTable.inclWater)) {
    OTUTable.RelAbund[j,i] <- OTUTable.inclWater[j,i]/colSumsOTUTable[i]
  }
}

OTUTable.P.RelAbund <- OTUTable.P.inclWater
colSumsOTUTable <- colSums(OTUTable.P.inclWater)
for (i in 1:ncol(OTUTable.P.inclWater)) {
  for (j in 1:nrow(OTUTable.P.inclWater)) {
    OTUTable.P.RelAbund[j,i] <- OTUTable.P.inclWater[j,i]/colSumsOTUTable[i]
  }
}

OTUTable.3.RelAbund <- OTUTable.3.inclWater
OTUTable.3.P.RelAbund <- OTUTable.3.P.inclWater


# Get list of random colors
set.seed(10)
nonGreyColors <- colors()[-grep("gr(a|e)y|darkolivegreen|white|light", colors())]
randomColors <- sample(nonGreyColors, nrow(OTUTable.3.RelAbund), replace = FALSE)
randomColorsEdited <- randomColors

# filter colors so that small groups are 'grey'
for (OTU in 1:nrow(OTUTable.3.RelAbund)) {
  if (max(c(as.numeric(OTUTable.3.RelAbund[OTU,]),as.numeric(OTUTable.3.P.RelAbund[OTU,]))) <= 0.05){
    randomColorsEdited[OTU] <- "white"
  }
}
# Then insert specific color for Oleispira
randomColorsEdited <- c("darkolivegreen",randomColorsEdited)
randomColors <- c("darkolivegreen", randomColors)

######### HAKAI #########

# Now, find OTU that is the abundant Oleispira
allOleispira <- rownames(taxonomyNames)[grep("Oleispira", taxonomyNames[,1])]
olsp <- names(which.max(rowSums(OTUTable.RelAbund[allOleispira,])))

# Make sure the .3 and not .3 have same order and same number of things

OTUTable.3.RelAbund <- OTUTable.3.RelAbund[,unlist(lapply(colnames(OTUTable.RelAbund), function(x) {
  grep(paste0("^",x,"$"), colnames(OTUTable.3.RelAbund))
}))]


# subtract relative abundance of this from the OTUTable.3 group
# Gammaproteobacteria = group3
# Oceanospirillales = group 4

row3Olsp <-grep("Oceanospirillales", rownames(OTUTable.3.RelAbund))
rowOlsp <- grep(paste0(olsp), rownames(OTUTable.RelAbund))
newGamma <- OTUTable.3.RelAbund[row3Olsp,] - OTUTable.RelAbund[rowOlsp,]
OTUTable.3.RelAbundFinal <- OTUTable.3.RelAbund
OTUTable.3.RelAbundFinal[row3Olsp,] <- newGamma
OTUTable.3.RelAbundFinal <- rbind(OTUTable.RelAbund[rowOlsp,],OTUTable.3.RelAbundFinal)

# Grep out FB and split by time
FB.only.OTU <- OTUTable.3.RelAbundFinal[,grep("FB", colnames(OTUTable.3.RelAbundFinal))]
# Get metadata
MF.inclWater$Time <- factor(MF.inclWater$Time)
# Get order
MF.inclWater <- MF.inclWater[with(MF.inclWater, order(Time,Rep)),]
FB.order <- rownames(MF.inclWater)[grep("FB", rownames(MF.inclWater))]
FB.only.OTU <- FB.only.OTU[,FB.order]

# Grep out BL and split by time
BL.only.OTU <- OTUTable.3.RelAbundFinal[,grep("BL", colnames(OTUTable.3.RelAbundFinal))]
# Get metadata
MF.inclWater$Time <- factor(MF.inclWater$Time)
# Get order
MF.inclWater <- MF.inclWater[with(MF.inclWater, order(Time,Rep)),]
BL.order <- rownames(MF.inclWater)[grep("BL", rownames(MF.inclWater))]
BL.only.OTU <- BL.only.OTU[,BL.order]

# Grep out CR and split by time
CR.only.OTU <- OTUTable.3.RelAbundFinal[,grep("CR", colnames(OTUTable.3.RelAbundFinal))]
# Get metadata
MF.inclWater$Time <- factor(MF.inclWater$Time)
# Get order
MF.inclWater <- MF.inclWater[with(MF.inclWater, order(Time,Rep)),]
CR.order <- rownames(MF.inclWater)[grep("CR", rownames(MF.inclWater))]
CR.only.OTU <- CR.only.OTU[,CR.order]

# Grep out water for Hakai
W.only.OTU <- OTUTable.3.RelAbundFinal[,grep("W", colnames(OTUTable.3.RelAbundFinal))]
W.order <- rownames(MF.inclWater)[grep("W", rownames(MF.inclWater))]
W.only.OTU <- W.only.OTU[,W.order]

# Count number of 20,60,360,720, 5670's- FB
count.FB.20 <- length(grep(".20.", colnames(FB.only.OTU), fixed = TRUE))
count.FB.60 <- length(grep(".60.", colnames(FB.only.OTU), fixed = TRUE))
count.FB.360 <- length(grep(".360.", colnames(FB.only.OTU), fixed = TRUE))
count.FB.720 <- length(grep(".720.", colnames(FB.only.OTU), fixed = TRUE))
count.FB.5760 <- length(grep(".5760.", colnames(FB.only.OTU), fixed = TRUE))

# Count number of 20,60,360,720, 5670's- BL
count.BL.20 <- length(grep(".20.", colnames(BL.only.OTU), fixed = TRUE))
count.BL.60 <- length(grep(".60.", colnames(BL.only.OTU), fixed = TRUE))
count.BL.360 <- length(grep(".360.", colnames(BL.only.OTU), fixed = TRUE))
count.BL.720 <- length(grep(".720.", colnames(BL.only.OTU), fixed = TRUE))
count.BL.5760 <- length(grep(".5760.", colnames(BL.only.OTU), fixed = TRUE))

# Count number of 20,60,360,720, 5670's- CR
count.CR.20 <- length(grep(".20.", colnames(CR.only.OTU), fixed = TRUE))
count.CR.60 <- length(grep(".60.", colnames(CR.only.OTU), fixed = TRUE))
count.CR.360 <- length(grep(".360.", colnames(CR.only.OTU), fixed = TRUE))
count.CR.720 <- length(grep(".720.", colnames(CR.only.OTU), fixed = TRUE))
count.CR.5760 <- length(grep(".5760.", colnames(CR.only.OTU), fixed = TRUE))

# Count number of 20,60,360,720, 5670's- CR
count.W.1 <- length(grep(".360", colnames(W.only.OTU), fixed = TRUE))
count.W.2 <- length(grep(".5760", colnames(W.only.OTU), fixed = TRUE))

pdf(paste0("./TAXASUMMARIES/Taxasummaries","Hakai",".pdf"), pointsize = 14)
par(fig = c(0,0.6,0.65,0.95), oma = c(0,3,0,0), mar = c(1,1,0,0))
barplot(as.matrix(FB.only.OTU)
        , col = randomColors
        , space = c(rep(0,count.FB.20)
                    ,0.5
                    ,rep(0,count.FB.60-1)
                    ,0.5
                    ,rep(0,count.FB.360-1)
                    ,0.5
                    ,rep(0,count.FB.720-1)
                    ,0.5
                    ,rep(0,count.FB.5760-1))
        , xaxt = 'n'
        , yaxt = 'n'
        , cex.axis = 0.5
        ,border = FALSE
        )
axis(1
     , labels = c("20m","1h","6h","12h","4d")
     , at = c(4,14,24,34,44)
     , tick = FALSE
     , cex = 0.1
     , las = 1
     , line = -1)
title(ylab = "FinelyBranched"
      , line = 0)
par(fig = c(0,0.6,0.35,0.65), oma = c(0,3,0,0), mar = c(1,1,0,0), new = TRUE)
barplot(as.matrix(BL.only.OTU)
        , col = randomColors
        , space = c(rep(0,count.BL.20)
                    ,0.5
                    ,rep(0,count.BL.60-1)
                    ,0.5
                    ,rep(0,count.BL.360-1)
                    ,0.5
                    ,rep(0,count.BL.720-1)
                    ,0.5
                    ,rep(0,count.BL.5760-1))
        , xaxt = 'n'
        , yaxt = 'n'
        , cex.axis = 0.5
        , border = FALSE
)
axis(1
     , labels = c("20m","1h","6h","12h","4d")
     , at = c(4,14,24,34,44)
     , tick = FALSE
     , cex = 0.1
     , las = 1
     , line = -1)
title(ylab = "Bladed"
      , line = 0)
par(fig = c(0,0.6,0.05,0.35), oma = c(0,3,0,0), mar = c(1,1,0,0),new = TRUE)
barplot(as.matrix(CR.only.OTU)
        , col = randomColors
        , space = c(rep(0,count.CR.20)
                    ,0.5
                    ,rep(0,count.CR.60-1)
                    ,0.5
                    ,rep(0,count.CR.360-1)
                    ,0.5
                    ,rep(0,count.CR.720-1)
                    ,0.5
                    ,rep(0,count.CR.5760-1))
        , xaxt = 'n'
        , yaxt = 'n'
        , cex.axis = 0.5
        , border = FALSE
)
axis(1
     , labels = c("20m","1h","6h","12h","4d")
     , at = c(4,14,24,32,40)
     , tick = FALSE
     , cex = 0.1
     , las = 1
     , line = -1)
title(ylab = "Crustose"
      , line = 0)
par(fig = c(0.6,0.9,0.65,0.95), oma = c(0,3,0,0), mar = c(1,1,0,0),new = TRUE)
barplot(as.matrix(W.only.OTU)
        , col = randomColors
        , space = c(rep(0,count.W.1)
                    ,0.5
                    ,rep(0,count.W.2-1)
                    )
        , xaxt = 'n'
        , yaxt = 'n'
        , cex.axis = 0.5
        , border = FALSE
)
axis(1
     , labels = c("6h","4d")
     , at = c(1,5)
     , tick = FALSE
     , cex = 0.1
     , las = 1
     , line = -1)
title(ylab = "Water"
      , line = 0)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0),new = TRUE)
plot(0,0
     , pch = ''
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = 'n'
     , bty = 'n')
title(ylab = "RELATIVE ABUNDANCE"
      , line = -2)



# legend("right"
#        , pch = 19
#        , legend = rownames(CR.only.OTU)
#        , col = randomColors
#        , cex = 0.4)
dev.off()
########## PM ##########

# Now, find OTU that is the abundant Oleispira
allOleispira <- rownames(taxonomyNames)[grep("Oleispira", taxonomyNames[,1])]
olsp <- names(which.max(rowSums(OTUTable.P.RelAbund[allOleispira,])))

# Make sure the .3 and not .3 have same order and same number of things

OTUTable.3.P.RelAbund <- OTUTable.3.P.RelAbund[,unlist(lapply(colnames(OTUTable.P.RelAbund), function(x) {
  grep(paste0("^",x,"$"), colnames(OTUTable.3.P.RelAbund))
}))]

# subtract relative abundance of this from the OTUTable.3 group
# Gammaproteobacteria = group3
# Oceanospirillales = group 4

row3Olsp <-grep("Oceanospirillales", rownames(OTUTable.3.P.RelAbund))
rowOlsp <- grep(paste0(olsp), rownames(OTUTable.P.RelAbund))
newGamma <- OTUTable.3.P.RelAbund[row3Olsp,] - OTUTable.P.RelAbund[rowOlsp,]
OTUTable.3.P.RelAbundFinal <- OTUTable.3.P.RelAbund
OTUTable.3.P.RelAbundFinal[row3Olsp,] <- newGamma
OTUTable.3.P.RelAbundFinal <- rbind(OTUTable.P.RelAbund[rowOlsp,],OTUTable.3.P.RelAbundFinal)

# Grep out FB and split by time
FB.only.OTU <- OTUTable.3.P.RelAbundFinal[,grep("FB", colnames(OTUTable.3.P.RelAbundFinal))]
# Get metadata
MF.P.inclWater$Time <- factor(MF.P.inclWater$Time)
# Get order
MF.P.inclWater <- MF.P.inclWater[with(MF.P.inclWater, order(Time,Rep)),]
FB.order <- rownames(MF.P.inclWater)[grep("FB", rownames(MF.P.inclWater))]
FB.only.OTU <- FB.only.OTU[,FB.order]

# Grep out BL and split by time
BL.only.OTU <- OTUTable.3.P.RelAbundFinal[,grep("BL", colnames(OTUTable.3.P.RelAbundFinal))]
# Get metadata
MF.P.inclWater$Time <- factor(MF.P.inclWater$Time)
# Get order
MF.P.inclWater <- MF.P.inclWater[with(MF.P.inclWater, order(Time,Rep)),]
BL.order <- rownames(MF.P.inclWater)[grep("BL", rownames(MF.P.inclWater))]
BL.only.OTU <- BL.only.OTU[,BL.order]

# Grep out CR and split by time
CR.only.OTU <- OTUTable.3.P.RelAbundFinal[,grep("CR", colnames(OTUTable.3.P.RelAbundFinal))]
# Get metadata
MF.P.inclWater$Time <- factor(MF.P.inclWater$Time)
# Get order
MF.P.inclWater <- MF.P.inclWater[with(MF.P.inclWater, order(Time,Rep)),]
CR.order <- rownames(MF.P.inclWater)[grep("CR", rownames(MF.P.inclWater))]
CR.only.OTU <- CR.only.OTU[,CR.order]

# Grep out water for PM
W.only.OTU <- OTUTable.3.P.RelAbundFinal[,grep("H2O", colnames(OTUTable.3.P.RelAbundFinal))]
W.order <- rownames(MF.P.inclWater)[grep("H2O", rownames(MF.P.inclWater))]
W.only.OTU <- W.only.OTU[,W.order]

# Count number of 20,60,360,720, 5670's - FB
count.FB.20 <- length(grep("^20[.]", colnames(FB.only.OTU)))
count.FB.60 <- length(grep("^60[.]", colnames(FB.only.OTU)))
count.FB.180 <- length(grep("^180[.]", colnames(FB.only.OTU)))
count.FB.360 <- length(grep("^360[.]", colnames(FB.only.OTU)))
count.FB.720 <- length(grep("^720[.]", colnames(FB.only.OTU)))
count.FB.1440 <- length(grep("^1440[.]", colnames(FB.only.OTU)))

# Count number of 20,60,360,720, 5670's - BL
count.BL.20 <- length(grep("^20[.]", colnames(BL.only.OTU)))
count.BL.60 <- length(grep("^60[.]", colnames(BL.only.OTU)))
count.BL.180 <- length(grep("^180[.]", colnames(BL.only.OTU)))
count.BL.360 <- length(grep("^360[.]", colnames(BL.only.OTU)))
count.BL.720 <- length(grep("^720[.]", colnames(BL.only.OTU)))
count.BL.1440 <- length(grep("^1440[.]", colnames(BL.only.OTU)))

# Count number of 20,60,360,720, 5670's - CR
count.CR.20 <- length(grep("^20[.]", colnames(CR.only.OTU)))
count.CR.60 <- length(grep("^60[.]", colnames(CR.only.OTU)))
count.CR.180 <- length(grep("^180[.]", colnames(CR.only.OTU)))
count.CR.360 <- length(grep("^360[.]", colnames(CR.only.OTU)))
count.CR.720 <- length(grep("^720[.]", colnames(CR.only.OTU)))
count.CR.1440 <- length(grep("^1440[.]", colnames(CR.only.OTU)))

# Count number of 20,60,360,720, 5670's- CR
count.W.1 <- length(grep("20.", colnames(W.only.OTU), fixed = TRUE))
count.W.2 <- length(grep("1440.", colnames(W.only.OTU), fixed = TRUE))


pdf(paste0("./TAXASUMMARIES/Taxasummaries","PM",".pdf"), pointsize = 14)
par(fig = c(0,0.6,0.65,0.95), oma = c(0,3,0,0), mar = c(1,1,0,0))
barplot(as.matrix(FB.only.OTU)
        , col = randomColors
        , space = c(rep(0,count.FB.20)
                    ,0.5
                    ,rep(0,count.FB.60-1)
                    ,0.5
                    ,rep(0,count.FB.180-1)
                    ,0.5
                    ,rep(0,count.FB.360-1)
                    ,0.5
                    ,rep(0,count.FB.720-1)
                    ,0.5
                    ,rep(0,count.FB.1440-1))
        , xaxt = 'n'
        , yaxt = 'n'
        , cex.axis = 0.5
        , border = FALSE
)
axis(1
     , labels = c("20m","1h","3h","6h","12h","24h")
     , at = c(1.5,5,8,11.5,14.5,18)
     , tick = FALSE
     , cex = 0.1
     , las = 1
     , line = -1)
title(ylab = "FinelyBranched"
      , line = 0)
par(fig = c(0,0.6,0.35,0.65), oma = c(0,3,0,0), mar = c(1,1,0,0), new = TRUE)
barplot(as.matrix(BL.only.OTU)
        , col = randomColors
        , space = c(rep(0,count.BL.20)
                    ,0.5
                    ,rep(0,count.BL.60-1)
                    ,0.5
                    ,rep(0,count.BL.180-1)
                    ,0.5
                    ,rep(0,count.BL.360-1)
                    ,0.5
                    ,rep(0,count.BL.720-1)
                    ,0.5
                    ,rep(0,count.BL.1440-1))
        , xaxt = 'n'
        , yaxt = 'n'
        , cex.axis = 0.5
        , border = FALSE
)
axis(1
     , labels = c("20m","1h","3h","6h","12h","24h")
     , at = c(1.5,5,8,11.5,15,18.5)
     , tick = FALSE
     , cex = 0.1
     , las = 1
     , line = -1)
title(ylab = "Bladed"
      , line = 0)
par(fig = c(0,0.6,0.05,0.35), oma = c(0,3,0,0), mar = c(1,1,0,0),new = TRUE)
barplot(as.matrix(CR.only.OTU)
        , col = randomColors
        , space = c(rep(0,count.CR.20)
                    ,0.5
                    ,rep(0,count.CR.60-1)
                    ,0.5
                    ,rep(0,count.CR.180-1)
                    ,0.5
                    ,rep(0,count.CR.360-1)
                    ,0.5
                    ,rep(0,count.CR.720-1)
                    ,0.5
                    ,rep(0,count.CR.1440-1))
        , xaxt = 'n'
        , yaxt = 'n'
        , cex.axis = 0.5
        , border = FALSE
)
axis(1
     , labels = c("20m","1h","3h","6h","12h","24h")
     , at = c(1.5,5,8,11.5,15,18.5)
     , tick = FALSE
     , cex = 0.1
     , las = 1
     , line = -1)
title(ylab = "Crustose"
      , line = 0)
par(fig = c(0.6,0.9,0.65,0.95), oma = c(0,3,0,0), mar = c(1,1,0,0),new = TRUE)
barplot(as.matrix(W.only.OTU)
        , col = randomColors
        , space = c(rep(0,count.W.1)
                    ,0.5
                    ,rep(0,count.W.2-1)
        )
        , xaxt = 'n'
        , yaxt = 'n'
        , cex.axis = 0.5
        , border = FALSE
)
axis(1
     , labels = c("20m","1d")
     , at = c(1,4)
     , tick = FALSE
     , cex = 0.1
     , las = 1
     , line = -1)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0),new = TRUE)
plot(0,0
     , pch = ''
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = 'n'
     , bty = 'n')
title(ylab = "RELATIVE ABUNDANCE"
      , line = -2)

dev.off()

########## Make Legend ############

# Get rownames of taxa that are above 5%, using 'white'.
over5 <- rownames(OTUTable.3.RelAbundFinal)[-grep('white', randomColorsEdited)]
oleispiraName <- taxonomyNames[grep(over5[1], rownames(taxonomyNames)),]
over5[1] <- paste0(oleispiraName)

over5colors <- randomColorsEdited[-grep('white', randomColorsEdited)]



pdf("./TAXASUMMARIES/LEGEND_taxasummaries.pdf", pointsize = 14, height = 6, width = 7.5)
par(mar = c(0,0,0,0))
plot(0,0
     , pch = ''
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = ''
     , ylab = ''
     , bty = 'n')
legend('center'
       , legend = rev(over5)
       , pch = 22
       , col = rev(over5colors)
       , pt.bg = rev(over5colors)
       , cex = 0.7
       , pt.cex = 2
       , bty = 'n'
)
dev.off()
# 
# ######### SHORT VERSION ##############
######### HAKAI ##########
# Collapse by replicate
OTUTable.RelAbund.trans <- t(OTUTable.RelAbund)

colRepNames <- c()
for (i in 1:ncol(OTUTable.RelAbund)) {
  colRepNames <- c(colRepNames, MF.inclWater[colnames(OTUTable.RelAbund)[i],"ColRep"])
}

rownames(OTUTable.RelAbund.trans) <- colRepNames
OTUTable.RelAbund.trans.agg <- aggregate(OTUTable.RelAbund.trans, by = list(row.names(OTUTable.RelAbund.trans)), FUN = mean)
OTUTable.RelAbund.trans.agg.sd <- aggregate(OTUTable.RelAbund.trans, by = list(row.names(OTUTable.RelAbund.trans)), FUN = sd)

headersTemp <- OTUTable.RelAbund.trans.agg[,1]
OTUTable.RelAbund.trans.agg <- OTUTable.RelAbund.trans.agg[,-1]
OTUTable.RelAbund.trans.agg.sd <- OTUTable.RelAbund.trans.agg.sd[,-1]
rownames(OTUTable.RelAbund.trans.agg) <- headersTemp
rownames(OTUTable.RelAbund.trans.agg.sd) <- headersTemp

OTUTable.RelAbund.trans.agg <- t(OTUTable.RelAbund.trans.agg)
OTUTable.RelAbund.trans.agg.sd <- t(OTUTable.RelAbund.trans.agg.sd)

# colnames(OTUTable.RelAbund.trans.agg) <- OTUTable.RelAbund.trans.agg[1,]
# OTUTable.RelAbund.trans.agg <- OTUTable.RelAbund.trans.agg[-1,]

# Order

orderPosition <- c()
for (i in c("FB","BL","CR","W")) {
  typeList <- colnames(OTUTable.RelAbund.trans.agg)[grep(i, colnames(OTUTable.RelAbund.trans.agg))]
  for (j in c("20","60","360","720","5760")) {
    timeList <- typeList[grep(paste0(i,j), typeList, fixed = TRUE)]
    orderPosition <- c(orderPosition, which(colnames(OTUTable.RelAbund.trans.agg) %in% timeList))
        }
}


OTUTable.RelAbund.trans.agg <- OTUTable.RelAbund.trans.agg[,orderPosition]
OTUTable.RelAbund.trans.agg.sd <- OTUTable.RelAbund.trans.agg.sd[,orderPosition]


# Get legend names

grThanFive <- c()
for (i in 1:nrow(OTUTable.RelAbund.trans.agg)) {
  grThanFive <- c(grThanFive, any(OTUTable.RelAbund.trans.agg[i,] >= 0.01))
}

namesKeep <- rownames(OTUTable.RelAbund.trans.agg)[grThanFive]

positionKeep <- unlist(sapply(namesKeep, function(x) {grep(paste0("^",x,"$"), rownames(taxonomyNames))}))

taxonomyNamesForPlot <- taxonomyNames[positionKeep,]
taxonomyNamesForPlotSplit <-sapply(taxonomyNamesForPlot, function(x) {
  tempString <- gsub("; __", "SPLIT", x)
  strsplit(tempString, split = "SPLIT")
  })

taxonomyNamesSplit <- c()
for (i in taxonomyNamesForPlotSplit) {
  if (is.na(i[3])) {
    Class <- "Unassigned"
  } else {
    Class <- i[3]
  }

  if (is.na(i[6])) {
    Genus <- "Unassigned"
  } else {
    Genus <- i[6]
  }

  if (is.na(i[7])) {
    Species <- "Unassigned"
  } else {
    Species <- i[7]
  }

  taxonomyNamesSplit <- c(taxonomyNamesSplit, paste0(Class, "--",Genus,"_",Species))
}


coloursMatched <- cbind(names(positionKeep), randomColors[1:length(positionKeep)])
colorsPlot <- cbind(rownames(OTUTable.RelAbund.trans.agg), rep("grey", nrow(OTUTable.RelAbund.trans.agg)))
count <- 1
for (i in 1:nrow(colorsPlot)) {
  if (colorsPlot[i,1] %in% coloursMatched[,1]) {
    colorsPlot[i,2] <- coloursMatched[count,2]
    count <- count + 1
  }
}

legendColors <- coloursMatched[,2]


pdf("./TAXASUMMARIES/TaxaSummaries_Hakai_OTU.pdf", pointsize = 14, width = 7, height = 5)
par(fig = c(0,1,0,1), mar = c(2,2,2,2))
barplot(as.matrix(OTUTable.RelAbund.trans.agg)
        , col = colorsPlot[,2]
        , border = NA
        , las = 2
        , space = c(0,0,0,0,0,2,0,0,0,0,2,0,0,0,0,2,0)
        , xaxt = "n"
        , yaxt = "n"

        )
title(ylab="Relative Abundance", line=0, cex.lab=1.2)
axis(1
     , las = 2
     , labels = c("20 m","1 h","6 h","12 h","4 d","",""
                  , "20 m","1 h","6 h","12 h","4 d","",""
                  , "20 m","1 h","6 h","12 h","4 d","",""
                  , "6 h","4 d", "")
     , at = seq(0.5,23.5, by = 1)
     , tick = FALSE
     , line = -1
    )
axis(1
     , labels = c("Finely Br.","Bladed", "Crustose", "Water")
     , at = c(2,9,16.5,21.5)
     , tick = FALSE
     , line = 1
     , cex = 0.5
)

dev.off()


########### PM #############
# Collapse by replicate
OTUTable.P.RelAbund.trans <- t(OTUTable.P.RelAbund)

colRepNames <- c()
for (i in 1:ncol(OTUTable.P.RelAbund)) {
  colRepNames <- c(colRepNames, MF.P.inclWater[colnames(OTUTable.P.RelAbund)[i],"ColRep"])
}

rownames(OTUTable.P.RelAbund.trans) <- colRepNames
OTUTable.P.RelAbund.trans.agg <- aggregate(OTUTable.P.RelAbund.trans, by = list(row.names(OTUTable.P.RelAbund.trans)), FUN = mean)
OTUTable.P.RelAbund.trans.agg.sd <- aggregate(OTUTable.P.RelAbund.trans, by = list(row.names(OTUTable.P.RelAbund.trans)), FUN = sd)

headersTemp <- OTUTable.P.RelAbund.trans.agg[,1]
OTUTable.P.RelAbund.trans.agg <- OTUTable.P.RelAbund.trans.agg[,-1]
OTUTable.P.RelAbund.trans.agg.sd <- OTUTable.P.RelAbund.trans.agg.sd[,-1]
rownames(OTUTable.P.RelAbund.trans.agg) <- headersTemp
rownames(OTUTable.P.RelAbund.trans.agg.sd) <- headersTemp

OTUTable.P.RelAbund.trans.agg <- t(OTUTable.P.RelAbund.trans.agg)
OTUTable.P.RelAbund.trans.agg.sd <- t(OTUTable.P.RelAbund.trans.agg.sd)

# colnames(OTUTable.P.RelAbund.trans.agg) <- OTUTable.P.RelAbund.trans.agg[1,]
# OTUTable.P.RelAbund.trans.agg <- OTUTable.P.RelAbund.trans.agg[-1,]

# Order

orderPosition <- c()
for (i in c("FB","BL","CR","H2O")) {
  typeList <- colnames(OTUTable.P.RelAbund.trans.agg)[grep(i, colnames(OTUTable.P.RelAbund.trans.agg))]
  for (j in c("20","60","180","360","720","1440")) {
    timeList <- typeList[grep(paste0(i,j), typeList, fixed = TRUE)]
    orderPosition <- c(orderPosition, which(colnames(OTUTable.P.RelAbund.trans.agg) %in% timeList))
  }
}


OTUTable.P.RelAbund.trans.agg <- OTUTable.P.RelAbund.trans.agg[,orderPosition]
OTUTable.P.RelAbund.trans.agg.sd <- OTUTable.P.RelAbund.trans.agg.sd[,orderPosition]



pdf("./TAXASUMMARIES/TaxaSummaries_PM_OTU.pdf", pointsize = 14, width = 7, height = 5)
par(fig = c(0,1,0,1), mar = c(2,2,2,2))
barplot(as.matrix(OTUTable.P.RelAbund.trans.agg)
        , col = colorsPlot[,2]
        , border = NA
        , las = 2
        , space = c(0,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,0,2,0)
        , xaxt = "n"
        , yaxt = "n"
        
)
title(ylab="Relative Abundance", line=0, cex.lab=1.2)
axis(1
     , las = 2
     , labels = c("20 m","1 h","3 h","6 h","12 h","1 d","",""
                  , "20 m","1 h","3 h","6 h","12 h","1 d","",""
                  , "20 m","1 h","3 h","6 h","12 h","1 d","",""
                  , "20 m","1 d", "")
     , at = seq(0.5,26.5, by = 1)
     , tick = FALSE
     , line = -1
)
axis(1
     , labels = c("Finely Br.","Bladed", "Crustose", "Water")
     , at = c(3,11,19.5,25.5)
     , tick = FALSE
     , line = 1
     , cex = 0.5
)

dev.off()

########## Make Legend for OTUs ############

pdf("./TAXASUMMARIES/LEGEND_taxasummaries_OTUs.pdf", pointsize = 14, height = 10, width = 7)
par(mar = c(0,0,0,0))
plot(0,0
     , pch = ''
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = ''
     , ylab = ''
     , bty = 'n')
legend("center"
       , legend = rev(taxonomyNamesSplit)
       , col = rev(legendColors)
       , pch = 15
       , cex = 0.65
       , pt.cex = 1.5
       , bty = 'n')

dev.off()


######## EXTRA-- for figuring out taxa
#
# taxonomyNamesSplitwOTUID <- taxonomyNamesSplit
# for (i in 1:length(taxonomyNamesSplitwOTUID)) {
#   taxonomyNamesSplitwOTUID[i] <- paste0(taxonomyNamesSplitwOTUID[i], "_", as.character(names(positionKeep)[i]))
# }
#
#
#
# jpeg("TaxaSummaries_morelegend.jpeg", pointsize = 32, width = 1500, height = 1500)
# par(fig = c(0,0.6,0,1), mar = c(2,2,2,2))
# barplot(as.matrix(OTUTable.RelAbund.trans.agg)
#         , col = randomColors
#         , border = NA
#         , las = 2
#         , space = c(0,0,0,0,0,2,0,0,0,0,2,0,0,0,0,2,0)
#         , xaxt = "n"
#         , yaxt = "n"
#
# )
# title(ylab="Relative Abundance", line=0, cex.lab=1.2)
# axis(1
#      , las = 2
#      , labels = c("20","60","360","720","5760","",""
#                   , "20", "60","360","720","5760","",""
#                   , "20","60","360","720","5760","",""
#                   , "360","5760", "")
#      , at = seq(0.5,23.5, by = 1)
#      , tick = FALSE
#      , line = -1
# )
# axis(1
#      , labels = c("Crust","Bladed", "Finely Br.", "Water")
#      , at = c(2,9,16.5,21.5)
#      , tick = FALSE
#      , line = 1
#      , cex = 0.5
# )
# par(fig = c(0.5,1,0,1), mar = c(0,0,0,0), new = TRUE)
# plot(0,0
#      , pch = ""
#      , xaxt = "n"
#      , yaxt = "n"
#      , xlab = ""
#      , ylab = ""
#      , bty = "n"
# )
# legend("center"
#        , legend = rev(taxonomyNamesSplitwOTUID)
#        , col = rev(legendColors)
#        , pch = 15
#        , cex = 0.5
#        , bty = 'n')
# dev.off()


# 
# # Isolating individual taxa to see abundances
# # xlab
# xlabelsTime <- factor(c("20","60","360","720","5760"), levels = c("20","60","360","720","5760"))
# 
# # Find error bars for each mean
# OTUTable.RelAbund.trans.agg.sd.filtered <- OTUTable.RelAbund.trans.agg.sd[positionKeep,]
# 
# OTUTable.RelAbund.trans.agg.sd.filtered["82837",]
# 
# colorOleispira <-  randomColors[positionKeep[grep(paste0("^","82837","$"), names(positionKeep))]]
# 
# # OLEISPIRA
# jpeg("RelAbund_Oleispira.jpeg")
# par(mar = c(5,5,4,4))
# plot(xlabelsTime,NULL
#      , pch = ""
#      , ylim = c(0,max(OTUTable.RelAbund.trans.agg["82837",]))
#      , xlab = "Time (Minutes)"
#      , ylab = "Relative Abundance"
#      , main = "Relative abundance of Oleispira across time"
#      , cex.lab = 1.5
#      , cex.main = 1.5
#      , xaxt = "n")
# axis(1, at=1:5, labels=xlabelsTime)
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82837",][1:5]
#        , lty = 1
#        , col = colorOleispira
#        , lwd = 3)
# arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
#        , x1 = c(0.99,1.99,2.99,3.99,4.99)
#        , y0 = c(OTUTable.RelAbund.trans.agg["82837",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["82837",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorOleispira
#        )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82837",][6:10]
#        , lty = 2
#        , col = colorOleispira
#        , lwd = 3)
# arrows(x0 = c(1,2,3,4,5)
#        , x1 = c(1,2,3,4,5)
#        , y0 = c(OTUTable.RelAbund.trans.agg["82837",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["82837",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorOleispira
# )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82837",][11:15]
#        , lty = 3
#        , col = colorOleispira
#        , lwd = 3)
# arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
#        , x1 = c(1.01,2.01,3.01,4.01,5.01)
#        , y0 = c(OTUTable.RelAbund.trans.agg["82837",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["82837",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorOleispira
# )
# legend("topleft"
#        , legend = c("Finely Br.","Bladed","Crustose")
#        , lty = c(1,2,3)
#        , col = c(colorOleispira,colorOleispira,colorOleispira)
#        , lwd = 3
#        )
# dev.off()
# 
# 
# 
# # PSUEDOMONAS FLUORESCENS
# colorfluorescens <-  randomColors[positionKeep[grep(paste0("^","82692","$"), names(positionKeep))]]
# 
# jpeg("RelAbund_Pseudomonas_fluorenscens.jpeg")
# par(mar = c(5,5,4,4))
# plot(xlabelsTime,NULL
#      , pch = ""
#      , ylim = c(0,max(OTUTable.RelAbund.trans.agg["82692",]))
#      , xlab = "Time (Minutes)"
#      , ylab = "Relative Abundance"
#      , main = "Relative abundance of Pseudomonas_fluorenscens across time"
#      , cex.lab = 1.5
#      , cex.main = 1.5
#      , xaxt = "n")
# axis(1, at=1:5, labels=xlabelsTime)
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82692",][1:5]
#       , lty = 1
#       , col = colorfluorescens
#       , lwd = 3)
# arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
#        , x1 = c(0.99,1.99,2.99,3.99,4.99)
#        , y0 = c(OTUTable.RelAbund.trans.agg["82692",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["82692",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorfluorescens
# )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82692",][6:10]
#       , lty = 2
#       , col = colorfluorescens
#       , lwd = 3)
# arrows(x0 = c(1,2,3,4,5)
#        , x1 = c(1,2,3,4,5)
#        , y0 = c(OTUTable.RelAbund.trans.agg["82692",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["82692",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorfluorescens
# )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82692",][11:15]
#       , lty = 3
#       , col = colorfluorescens
#       , lwd = 3)
# arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
#        , x1 = c(1.01,2.01,3.01,4.01,5.01)
#        , y0 = c(OTUTable.RelAbund.trans.agg["82692",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["82692",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorfluorescens
# )
# legend("topright"
#        , legend = c("Finely Br.","Bladed","Crustose")
#        , lty = c(1,2,3)
#        , col = c(colorfluorescens,colorfluorescens,colorfluorescens)
#        , lwd = 3
# )
# dev.off()
# 
# 
# # Pseudoalteromonas_unassigned
# colorPseudoalteromonas <-  randomColors[positionKeep[grep(paste0("^","91303","$"), names(positionKeep))]]
# 
# jpeg("RelAbund_Pseudoalteromonas_unassigned.jpeg")
# par(mar = c(5,5,4,4))
# plot(xlabelsTime,NULL
#      , pch = ""
#      , ylim = c(0,max(OTUTable.RelAbund.trans.agg["91303",]))
#      , xlab = "Time (Minutes)"
#      , ylab = "Relative Abundance"
#      , main = "Relative abundance of Pseudoalteromonas across time"
#      , cex.lab = 1.5
#      , cex.main = 1.5
#      , xaxt = "n")
# axis(1, at=1:5, labels=xlabelsTime)
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["91303",][1:5]
#       , lty = 1
#       , col = colorPseudoalteromonas
#       , lwd = 3)
# arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
#        , x1 = c(0.99,1.99,2.99,3.99,4.99)
#        , y0 = c(OTUTable.RelAbund.trans.agg["91303",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["91303",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorPseudoalteromonas
# )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["91303",][6:10]
#       , lty = 2
#       , col = colorPseudoalteromonas
#       , lwd = 3)
# arrows(x0 = c(1,2,3,4,5)
#        , x1 = c(1,2,3,4,5)
#        , y0 = c(OTUTable.RelAbund.trans.agg["91303",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["91303",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorPseudoalteromonas
# )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["91303",][11:15]
#       , lty = 3
#       , col = colorPseudoalteromonas
#       , lwd = 3)
# arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
#        , x1 = c(1.01,2.01,3.01,4.01,5.01)
#        , y0 = c(OTUTable.RelAbund.trans.agg["91303",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["91303",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorPseudoalteromonas
# )
# legend("topright"
#        , legend = c("Finely Br.","Bladed","Crustose")
#        , lty = c(1,2,3)
#        , col = c(colorPseudoalteromonas,colorPseudoalteromonas,colorPseudoalteromonas)
#        , lwd = 3
# )
# dev.off()
# 
# 
# 
# # Achromobacter
# colorAchromobacter <-  randomColors[positionKeep[grep(paste0("^","77068","$"), names(positionKeep))]]
# 
# jpeg("RelAbund_Achromobacter.jpeg")
# par(mar = c(5,5,4,4))
# plot(xlabelsTime,NULL
#      , pch = ""
#      , ylim = c(0,max(OTUTable.RelAbund.trans.agg["77068",]))
#      , xlab = "Time (Minutes)"
#      , ylab = "Relative Abundance"
#      , main = "Relative abundance of Achromobacter across time"
#      , cex.lab = 1.5
#      , cex.main = 1.5
#      , xaxt = "n")
# axis(1, at=1:5, labels=xlabelsTime)
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["77068",][1:5]
#       , lty = 1
#       , col = colorAchromobacter
#       , lwd = 3)
# arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
#        , x1 = c(0.99,1.99,2.99,3.99,4.99)
#        , y0 = c(OTUTable.RelAbund.trans.agg["77068",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["77068",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorAchromobacter
# )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["77068",][6:10]
#       , lty = 2
#       , col = colorAchromobacter
#       , lwd = 3)
# arrows(x0 = c(1,2,3,4,5)
#        , x1 = c(1,2,3,4,5)
#        , y0 = c(OTUTable.RelAbund.trans.agg["77068",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["77068",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorAchromobacter
# )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["77068",][11:15]
#       , lty = 3
#       , col = colorAchromobacter
#       , lwd = 3)
# arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
#        , x1 = c(1.01,2.01,3.01,4.01,5.01)
#        , y0 = c(OTUTable.RelAbund.trans.agg["77068",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["77068",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorAchromobacter
# )
# legend("bottomleft"
#        , legend = c("Finely Br.","Bladed","Crustose")
#        , lty = c(1,2,3)
#        , col = c(colorAchromobacter,colorAchromobacter,colorAchromobacter)
#        , lwd = 3
# )
# dev.off()
# 
# 
# 
# # Gammaproteobacteria_unassigned
# colorGammaproteobacteria_unassigned <-  randomColors[positionKeep[grep(paste0("^","92169","$"), names(positionKeep))]]
# 
# jpeg("RelAbund_Gammaproteobacteria_unassigned.jpeg")
# par(mar = c(5,5,4,4))
# plot(xlabelsTime,NULL
#      , pch = ""
#      , ylim = c(0,max(OTUTable.RelAbund.trans.agg["92169",]))
#      , xlab = "Time (Minutes)"
#      , ylab = "Relative Abundance"
#      , main = "Relative abundance of Gammaproteobacteria_unassigned across time"
#      , cex.lab = 1.5
#      , cex.main = 1.5
#      , xaxt = "n")
# axis(1, at=1:5, labels=xlabelsTime)
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["92169",][1:5]
#       , lty = 1
#       , col = colorGammaproteobacteria_unassigned
#       , lwd = 3)
# arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
#        , x1 = c(0.99,1.99,2.99,3.99,4.99)
#        , y0 = c(OTUTable.RelAbund.trans.agg["92169",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["92169",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorGammaproteobacteria_unassigned
# )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["92169",][6:10]
#       , lty = 2
#       , col = colorGammaproteobacteria_unassigned
#       , lwd = 3)
# arrows(x0 = c(1,2,3,4,5)
#        , x1 = c(1,2,3,4,5)
#        , y0 = c(OTUTable.RelAbund.trans.agg["92169",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["92169",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorGammaproteobacteria_unassigned
# )
# lines(xlabelsTime,OTUTable.RelAbund.trans.agg["92169",][11:15]
#       , lty = 3
#       , col = colorGammaproteobacteria_unassigned
#       , lwd = 3)
# arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
#        , x1 = c(1.01,2.01,3.01,4.01,5.01)
#        , y0 = c(OTUTable.RelAbund.trans.agg["92169",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , y1 = c(OTUTable.RelAbund.trans.agg["92169",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
#        , angle = 90
#        , code = 3
#        , length = 0
#        , col = colorGammaproteobacteria_unassigned
# )
# legend("topleft"
#        , legend = c("Finely Br.","Bladed","Crustose")
#        , lty = c(1,2,3)
#        , col = c(colorGammaproteobacteria_unassigned,colorGammaproteobacteria_unassigned,colorGammaproteobacteria_unassigned)
#        , lwd = 3
# )
# dev.off()
# 
# 
# 
# ########## COLLAPSE BY LEVEL ################
# 
# OTUTable.RelAbund

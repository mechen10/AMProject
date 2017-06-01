#!/bin/Rscript

setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/TAXASUMMARIES")
# Make Taxa summaries

OTUTable <- read.delim("OTU_table_Text_real.txt"
                       , header = TRUE
                       , skip = 1
                       , row.names = 1)
taxonomyNames <- as.data.frame(OTUTable[,ncol(OTUTable)])
rownames(taxonomyNames) <- rownames(OTUTable)

OTUTable <- OTUTable[,-ncol(OTUTable)]
colnames(OTUTable) <- gsub('X', '',colnames(OTUTable))

# Get rid of FB.20.6 and FB.20.9
OTUTable <- OTUTable[,-grep("FB.20.6", colnames(OTUTable))]
OTUTable <- OTUTable[,-grep("FB.20.9", colnames(OTUTable))]

MF <- read.delim("MF_All_m1000.txt"
                       , header = TRUE
                       , row.names = 1
                       , stringsAsFactors = FALSE)
rownames(MF) <- gsub("-",".", rownames(MF))

# Make sure they are mutually in each other

MF <- MF[rownames(MF) %in% colnames(OTUTable),]
OTUTable <- OTUTable[,colnames(OTUTable) %in% rownames(MF)]

# Get rid of all non-hakai samples
# Get all samples I want; remove the rest
MF.inclWater <- MF
OTUTable.inclWater <- OTUTable
for (ROW in 1:nrow(MF)) {
  typeMorph <- MF[ROW, 'Morph']
  typeType <- MF[ROW, 'Type']
  rownameTemp <- rownames(MF)[ROW]
  if (!(typeMorph %in% c('CR',"BL","FB","W"))) {
    MF.inclWater <- MF.inclWater[-grep(paste0("^",rownameTemp,"$"), rownames(MF.inclWater)),]
    OTUTable.inclWater <- OTUTable.inclWater[,-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(OTUTable.inclWater))]
  } else if (!(typeType == "H")) {
    MF.inclWater <- MF.inclWater[-grep(paste0("^",rownameTemp,"$"), rownames(MF.inclWater)),]
    OTUTable.inclWater <- OTUTable.inclWater[,-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(OTUTable.inclWater))]
  }
}

# Get rid of zeros
toDelete <- c()
for (i in 1:nrow(OTUTable.inclWater)) {
  if (sum(OTUTable.inclWater[i,]) <= 5) {
    toDelete <- c(toDelete, i)
  }
}

OTUTable.inclWater <- OTUTable.inclWater[-toDelete,]

# Now make relative abundance

OTUTable.RelAbund <- OTUTable.inclWater
colSumsOTUTable <- colSums(OTUTable.inclWater)
for (i in 1:ncol(OTUTable.inclWater)) {
  for (j in 1:nrow(OTUTable.inclWater)) {
    OTUTable.RelAbund[j,i] <- OTUTable.inclWater[j,i]/colSumsOTUTable[i]
  }
}

# Get list of random colors
set.seed(5)
nonGreyColors <- colors()[-grep("gr(a|e)y", colors())]
randomColors <- sample(nonGreyColors, nrow(OTUTable.RelAbund), replace = TRUE)

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



legendColors <- randomColors[positionKeep]

jpeg("TaxaSummaries.jpeg", pointsize = 32, width = 1000, height = 1500)
par(fig = c(0,1,0.4,1), mar = c(2,2,2,2))
barplot(as.matrix(OTUTable.RelAbund.trans.agg)
        , col = randomColors
        , border = NA
        , las = 2
        , space = c(0,0,0,0,0,2,0,0,0,0,2,0,0,0,0,2,0)
        , xaxt = "n"
        , yaxt = "n"

        )
title(ylab="Relative Abundance", line=0, cex.lab=1.2)
axis(1
     , las = 2
     , labels = c("20","60","360","720","5760","",""
                  , "20", "60","360","720","5760","",""
                  , "20","60","360","720","5760","",""
                  , "360","5760", "")
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
par(fig = c(0,1,0,0.4), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n"
     )
legend("center"
       , legend = rev(taxonomyNamesSplit)
       , col = rev(legendColors)
       , pch = 15
       , cex = 0.8
       , bty = 'n')
dev.off()

OTUTable.RelAbund.trans.agg["82837",]




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



# Isolating individual taxa to see abundances
# xlab
xlabelsTime <- factor(c("20","60","360","720","5760"), levels = c("20","60","360","720","5760"))

# Find error bars for each mean
OTUTable.RelAbund.trans.agg.sd.filtered <- OTUTable.RelAbund.trans.agg.sd[positionKeep,]

OTUTable.RelAbund.trans.agg.sd.filtered["82837",]

colorOleispira <-  randomColors[positionKeep[grep(paste0("^","82837","$"), names(positionKeep))]]

# OLEISPIRA
jpeg("RelAbund_Oleispira.jpeg")
par(mar = c(5,5,4,4))
plot(xlabelsTime,NULL
     , pch = ""
     , ylim = c(0,max(OTUTable.RelAbund.trans.agg["82837",]))
     , xlab = "Time (Minutes)"
     , ylab = "Relative Abundance"
     , main = "Relative abundance of Oleispira across time"
     , cex.lab = 1.5
     , cex.main = 1.5
     , xaxt = "n")
axis(1, at=1:5, labels=xlabelsTime)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82837",][1:5]
       , lty = 1
       , col = colorOleispira
       , lwd = 3)
arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
       , x1 = c(0.99,1.99,2.99,3.99,4.99)
       , y0 = c(OTUTable.RelAbund.trans.agg["82837",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["82837",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorOleispira
       )
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82837",][6:10]
       , lty = 2
       , col = colorOleispira
       , lwd = 3)
arrows(x0 = c(1,2,3,4,5)
       , x1 = c(1,2,3,4,5)
       , y0 = c(OTUTable.RelAbund.trans.agg["82837",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["82837",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorOleispira
)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82837",][11:15]
       , lty = 3
       , col = colorOleispira
       , lwd = 3)
arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
       , x1 = c(1.01,2.01,3.01,4.01,5.01)
       , y0 = c(OTUTable.RelAbund.trans.agg["82837",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["82837",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorOleispira
)
legend("topleft"
       , legend = c("Finely Br.","Bladed","Crustose")
       , lty = c(1,2,3)
       , col = c(colorOleispira,colorOleispira,colorOleispira)
       , lwd = 3
       )
dev.off()



# PSUEDOMONAS FLUORESCENS
colorfluorescens <-  randomColors[positionKeep[grep(paste0("^","82692","$"), names(positionKeep))]]

jpeg("RelAbund_Pseudomonas_fluorenscens.jpeg")
par(mar = c(5,5,4,4))
plot(xlabelsTime,NULL
     , pch = ""
     , ylim = c(0,max(OTUTable.RelAbund.trans.agg["82692",]))
     , xlab = "Time (Minutes)"
     , ylab = "Relative Abundance"
     , main = "Relative abundance of Pseudomonas_fluorenscens across time"
     , cex.lab = 1.5
     , cex.main = 1.5
     , xaxt = "n")
axis(1, at=1:5, labels=xlabelsTime)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82692",][1:5]
      , lty = 1
      , col = colorfluorescens
      , lwd = 3)
arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
       , x1 = c(0.99,1.99,2.99,3.99,4.99)
       , y0 = c(OTUTable.RelAbund.trans.agg["82692",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["82692",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorfluorescens
)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82692",][6:10]
      , lty = 2
      , col = colorfluorescens
      , lwd = 3)
arrows(x0 = c(1,2,3,4,5)
       , x1 = c(1,2,3,4,5)
       , y0 = c(OTUTable.RelAbund.trans.agg["82692",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["82692",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorfluorescens
)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["82692",][11:15]
      , lty = 3
      , col = colorfluorescens
      , lwd = 3)
arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
       , x1 = c(1.01,2.01,3.01,4.01,5.01)
       , y0 = c(OTUTable.RelAbund.trans.agg["82692",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["82692",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorfluorescens
)
legend("topright"
       , legend = c("Finely Br.","Bladed","Crustose")
       , lty = c(1,2,3)
       , col = c(colorfluorescens,colorfluorescens,colorfluorescens)
       , lwd = 3
)
dev.off()


# Pseudoalteromonas_unassigned
colorPseudoalteromonas <-  randomColors[positionKeep[grep(paste0("^","91303","$"), names(positionKeep))]]

jpeg("RelAbund_Pseudoalteromonas_unassigned.jpeg")
par(mar = c(5,5,4,4))
plot(xlabelsTime,NULL
     , pch = ""
     , ylim = c(0,max(OTUTable.RelAbund.trans.agg["91303",]))
     , xlab = "Time (Minutes)"
     , ylab = "Relative Abundance"
     , main = "Relative abundance of Pseudoalteromonas across time"
     , cex.lab = 1.5
     , cex.main = 1.5
     , xaxt = "n")
axis(1, at=1:5, labels=xlabelsTime)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["91303",][1:5]
      , lty = 1
      , col = colorPseudoalteromonas
      , lwd = 3)
arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
       , x1 = c(0.99,1.99,2.99,3.99,4.99)
       , y0 = c(OTUTable.RelAbund.trans.agg["91303",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["91303",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorPseudoalteromonas
)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["91303",][6:10]
      , lty = 2
      , col = colorPseudoalteromonas
      , lwd = 3)
arrows(x0 = c(1,2,3,4,5)
       , x1 = c(1,2,3,4,5)
       , y0 = c(OTUTable.RelAbund.trans.agg["91303",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["91303",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorPseudoalteromonas
)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["91303",][11:15]
      , lty = 3
      , col = colorPseudoalteromonas
      , lwd = 3)
arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
       , x1 = c(1.01,2.01,3.01,4.01,5.01)
       , y0 = c(OTUTable.RelAbund.trans.agg["91303",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["91303",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorPseudoalteromonas
)
legend("topright"
       , legend = c("Finely Br.","Bladed","Crustose")
       , lty = c(1,2,3)
       , col = c(colorPseudoalteromonas,colorPseudoalteromonas,colorPseudoalteromonas)
       , lwd = 3
)
dev.off()



# Achromobacter
colorAchromobacter <-  randomColors[positionKeep[grep(paste0("^","77068","$"), names(positionKeep))]]

jpeg("RelAbund_Achromobacter.jpeg")
par(mar = c(5,5,4,4))
plot(xlabelsTime,NULL
     , pch = ""
     , ylim = c(0,max(OTUTable.RelAbund.trans.agg["77068",]))
     , xlab = "Time (Minutes)"
     , ylab = "Relative Abundance"
     , main = "Relative abundance of Achromobacter across time"
     , cex.lab = 1.5
     , cex.main = 1.5
     , xaxt = "n")
axis(1, at=1:5, labels=xlabelsTime)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["77068",][1:5]
      , lty = 1
      , col = colorAchromobacter
      , lwd = 3)
arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
       , x1 = c(0.99,1.99,2.99,3.99,4.99)
       , y0 = c(OTUTable.RelAbund.trans.agg["77068",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["77068",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorAchromobacter
)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["77068",][6:10]
      , lty = 2
      , col = colorAchromobacter
      , lwd = 3)
arrows(x0 = c(1,2,3,4,5)
       , x1 = c(1,2,3,4,5)
       , y0 = c(OTUTable.RelAbund.trans.agg["77068",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["77068",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorAchromobacter
)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["77068",][11:15]
      , lty = 3
      , col = colorAchromobacter
      , lwd = 3)
arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
       , x1 = c(1.01,2.01,3.01,4.01,5.01)
       , y0 = c(OTUTable.RelAbund.trans.agg["77068",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["77068",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorAchromobacter
)
legend("bottomleft"
       , legend = c("Finely Br.","Bladed","Crustose")
       , lty = c(1,2,3)
       , col = c(colorAchromobacter,colorAchromobacter,colorAchromobacter)
       , lwd = 3
)
dev.off()



# Gammaproteobacteria_unassigned
colorGammaproteobacteria_unassigned <-  randomColors[positionKeep[grep(paste0("^","92169","$"), names(positionKeep))]]

jpeg("RelAbund_Gammaproteobacteria_unassigned.jpeg")
par(mar = c(5,5,4,4))
plot(xlabelsTime,NULL
     , pch = ""
     , ylim = c(0,max(OTUTable.RelAbund.trans.agg["92169",]))
     , xlab = "Time (Minutes)"
     , ylab = "Relative Abundance"
     , main = "Relative abundance of Gammaproteobacteria_unassigned across time"
     , cex.lab = 1.5
     , cex.main = 1.5
     , xaxt = "n")
axis(1, at=1:5, labels=xlabelsTime)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["92169",][1:5]
      , lty = 1
      , col = colorGammaproteobacteria_unassigned
      , lwd = 3)
arrows(x0 = c(0.99,1.99,2.99,3.99,4.99)
       , x1 = c(0.99,1.99,2.99,3.99,4.99)
       , y0 = c(OTUTable.RelAbund.trans.agg["92169",][1:5] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["92169",][1:5] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][1:5]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorGammaproteobacteria_unassigned
)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["92169",][6:10]
      , lty = 2
      , col = colorGammaproteobacteria_unassigned
      , lwd = 3)
arrows(x0 = c(1,2,3,4,5)
       , x1 = c(1,2,3,4,5)
       , y0 = c(OTUTable.RelAbund.trans.agg["92169",][6:10] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["92169",][6:10] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][6:10]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorGammaproteobacteria_unassigned
)
lines(xlabelsTime,OTUTable.RelAbund.trans.agg["92169",][11:15]
      , lty = 3
      , col = colorGammaproteobacteria_unassigned
      , lwd = 3)
arrows(x0 = c(1.01,2.01,3.01,4.01,5.01)
       , x1 = c(1.01,2.01,3.01,4.01,5.01)
       , y0 = c(OTUTable.RelAbund.trans.agg["92169",][11:15] - OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , y1 = c(OTUTable.RelAbund.trans.agg["92169",][11:15] + OTUTable.RelAbund.trans.agg.sd.filtered["82837",][11:15]/2)
       , angle = 90
       , code = 3
       , length = 0
       , col = colorGammaproteobacteria_unassigned
)
legend("topleft"
       , legend = c("Finely Br.","Bladed","Crustose")
       , lty = c(1,2,3)
       , col = c(colorGammaproteobacteria_unassigned,colorGammaproteobacteria_unassigned,colorGammaproteobacteria_unassigned)
       , lwd = 3
)
dev.off()


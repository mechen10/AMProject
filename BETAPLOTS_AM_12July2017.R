#!/bin/RScript
set.seed(10)
# MAKE BETA PLOTS
library("MASS")
library("vegan")
library("nlme")
library("optparse")
library("xtable")
########################### OPT PARSE #################################


option_list = list(
  make_option(c("-b", "--braycurtisdm"), type="character",
              help="Distance Matrix- Braycurtis"),
  make_option(c("-u", "--unweightedunifracdm"),
              help="Distance Matrix- Unweighted Unifrac", type="character"),
  make_option(c("-w", "--weightedunifracdm"), 
              help="Distance Matrix- Weighted Unifrac", type="character"),
  make_option(c("-m", "--mappingfile"),
              help="Mappingfile", type="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

BCPWD = opt$braycurtisdm
UWUFPWD = opt$unweightedunifracdm
WUFPWD = opt$weightedunifracdm
MFPWD = opt$mappingfile
########################### FOR TESTING #################################

# BCPWD <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/ANALYSIS_ALPHABETATAXA/beta_div/bray_curtis_dm.txt'
# UWUFPWD <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/ANALYSIS_ALPHABETATAXA/beta_div/unweighted_unifrac_dm.txt'
# WUFPWD <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/ANALYSIS_ALPHABETATAXA/beta_div/weighted_unifrac_dm.txt'
# MFPWD <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/OTU_MP_filt/MF_nochlpmito_m1000.txt'
# setwd("/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis")
# 
BCPWD <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/TEMP_frombotclust/dm/bray_curtis_dm.txt'
UWUFPWD <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/TEMP_frombotclust/dm/unweighted_unifrac_dm.txt'
WUFPWD <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/TEMP_frombotclust/dm/weighted_unifrac_dm.txt'
MFPWD <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/TEMP_frombotclust/MF_nochlpmito_m1000.txt'
setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis")

########################### LOADING #################################
 

dm.UWUF <- read.delim(paste0(UWUFPWD)
                      , header = TRUE
                      , row.names = 1
                      , stringsAsFactors = FALSE)
colnames(dm.UWUF) <- rownames(dm.UWUF)

dm.WUF <- read.delim(paste0(WUFPWD)
                     , header = TRUE
                     , row.names = 1
                     , stringsAsFactors = FALSE)
colnames(dm.WUF) <- rownames(dm.WUF)

dm.BC <- read.delim(paste0(BCPWD)
                    , header = TRUE
                    , row.names = 1
                    , stringsAsFactors = FALSE)
colnames(dm.BC) <- rownames(dm.BC)


MF <- read.delim(paste0(MFPWD)
                 , stringsAsFactors = FALSE
                 , header = TRUE
                 , row.names= 1
                 , strip.white = TRUE)



# REMOVE OUTLIERS
# CR-20-5
# 60.CR.3 
# FB-20-6
# FB-20-9

dm.UWUF <- dm.UWUF[-grep("^CR-20-5$", rownames(dm.UWUF)),-grep("^CR-20-5$", colnames(dm.UWUF))]
dm.WUF <- dm.WUF[-grep("^CR-20-5$", rownames(dm.WUF)),-grep("^CR-20-5$", colnames(dm.WUF))]
dm.BC <- dm.BC[-grep("^CR-20-5$", rownames(dm.BC)),-grep("^CR-20-5$", colnames(dm.BC))]


dm.UWUF <- dm.UWUF[-grep("^60.CR.3$", rownames(dm.UWUF)),-grep("^60.CR.3$", colnames(dm.UWUF))]
dm.WUF <- dm.WUF[-grep("^60.CR.3$", rownames(dm.WUF)),-grep("^60.CR.3$", colnames(dm.WUF))]
dm.BC <- dm.BC[-grep("^60.CR.3$", rownames(dm.BC)),-grep("^60.CR.3$", colnames(dm.BC))]


dm.UWUF <- dm.UWUF[-grep("^FB-20-6$", rownames(dm.UWUF)),-grep("^FB-20-6$", colnames(dm.UWUF))]
dm.WUF <- dm.WUF[-grep("^FB-20-6$", rownames(dm.WUF)),-grep("^FB-20-6$", colnames(dm.WUF))]
dm.BC <- dm.BC[-grep("^FB-20-6$", rownames(dm.BC)),-grep("^FB-20-6$", colnames(dm.BC))]

dm.UWUF <- dm.UWUF[-grep("^FB-20-9$", rownames(dm.UWUF)),-grep("^FB-20-9$", colnames(dm.UWUF))]
dm.WUF <- dm.WUF[-grep("^FB-20-9$", rownames(dm.WUF)),-grep("^FB-20-9$", colnames(dm.WUF))]
dm.BC <- dm.BC[-grep("^FB-20-9$", rownames(dm.BC)),-grep("^FB-20-9$", colnames(dm.BC))]


# Make sure everything in dm is in MF

dm.UWUF <- dm.UWUF[unlist(lapply(rownames(MF), function(x) {
  grep(paste0("^",x,"$"), rownames(dm.UWUF))
})), unlist(lapply(rownames(MF), function(x) {
  grep(paste0("^",x,"$"), colnames(dm.UWUF))
}))]

dm.WUF <- dm.WUF[unlist(lapply(rownames(MF), function(x) {
  grep(paste0("^",x,"$"), rownames(dm.WUF))
})), unlist(lapply(rownames(MF), function(x) {
  grep(paste0("^",x,"$"), colnames(dm.WUF))
}))]


dm.BC <- dm.BC[unlist(lapply(rownames(MF), function(x) {
  grep(paste0("^",x,"$"), rownames(dm.BC))
})), unlist(lapply(rownames(MF), function(x) {
  grep(paste0("^",x,"$"), colnames(dm.BC))
}))]


# MAKE ONLY HAKAI (but save pm for later)
dm.UWUF.P <- dm.UWUF[unlist(sapply(rownames(dm.UWUF), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "P"
})), unlist(sapply(colnames(dm.UWUF), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "P"
}))]

dm.UWUF <- dm.UWUF[unlist(sapply(rownames(dm.UWUF), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "H"
})), unlist(sapply(colnames(dm.UWUF), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "H"
}))]


dm.WUF.P <- dm.WUF[unlist(sapply(rownames(dm.WUF), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "P"
})), unlist(sapply(colnames(dm.WUF), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "P"
}))]

dm.WUF <- dm.WUF[unlist(sapply(rownames(dm.WUF), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "H"
})), unlist(sapply(colnames(dm.WUF), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "H"
}))]

dm.BC.P <- dm.BC[unlist(sapply(rownames(dm.BC), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "P"
})), unlist(sapply(colnames(dm.BC), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "P"
}))]

dm.BC <- dm.BC[unlist(sapply(rownames(dm.BC), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "H"
})), unlist(sapply(colnames(dm.BC), function(x) {
  MF[grep(paste0("^",x,"$"), rownames(MF)), 'Type'] == "H"
}))]


# Make sure everything in MF is in dm (only need to use one dm)
MF.P <- MF[unlist(lapply(rownames(dm.UWUF.P), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
})),]

MF <- MF[unlist(lapply(rownames(dm.UWUF), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
})),]

# Get morph out
MF.morphkeep <- MF
dm.UWUF.morphonly <- dm.UWUF
dm.WUF.morphonly <- dm.WUF
dm.BC.morphonly <- dm.BC
# ROW <- 67
for (ROW in 1:nrow(MF)) {
  typeMorph <- MF[ROW, 'Morph']
  typeTime <- MF[ROW, 'Time']
  rownameTemp <- rownames(MF)[ROW]
  if ((!(typeMorph %in% c('CR',"BL","FB"))) | ((typeTime %in% c("5760")))) {
    MF.morphkeep <- MF.morphkeep[-grep(paste0("^",rownameTemp,"$"), rownames(MF.morphkeep)),]
    
    dm.UWUF.morphonly <- dm.UWUF.morphonly[-grep(paste0("^",rownames(MF)[ROW],"$"), rownames(dm.UWUF.morphonly)),-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(dm.UWUF.morphonly))]
    dm.WUF.morphonly <- dm.WUF.morphonly[-grep(paste0("^",rownames(MF)[ROW],"$"), rownames(dm.WUF.morphonly)),-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(dm.WUF.morphonly))]
    dm.BC.morphonly <- dm.BC.morphonly[-grep(paste0("^",rownames(MF)[ROW],"$"), rownames(dm.BC.morphonly)),-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(dm.BC.morphonly))]
    }
}

# Get at 5760
MF.morphkeep.5760 <- MF.morphkeep[-grep("5760", rownames(MF.morphkeep)),]

MF.P.morphkeep <- MF.P
dm.UWUF.P.morphonly <- dm.UWUF.P
dm.WUF.P.morphonly <- dm.WUF.P
dm.BC.P.morphonly <- dm.BC.P
# ROW <- 67
for (ROW in 1:nrow(MF.P)) {
  typeMorph <- MF.P[ROW, 'Morph']
  typeTime <- MF.P[ROW, 'Time']
  rownameTemp <- rownames(MF.P)[ROW]
  if ((!(typeMorph %in% c('CR',"BL","FB"))) | ((typeTime %in% c("5760")))) {
    MF.P.morphkeep <- MF.P.morphkeep[-grep(paste0("^",rownameTemp,"$"), rownames(MF.P.morphkeep)),]
    
    dm.UWUF.P.morphonly <- dm.UWUF.P.morphonly[-grep(paste0("^",rownames(MF.P)[ROW],"$"), rownames(dm.UWUF.P.morphonly)),-grep(paste0("^",rownames(MF.P)[ROW],"$"), colnames(dm.UWUF.P.morphonly))]
    dm.WUF.P.morphonly <- dm.WUF.P.morphonly[-grep(paste0("^",rownames(MF.P)[ROW],"$"), rownames(dm.WUF.P.morphonly)),-grep(paste0("^",rownames(MF.P)[ROW],"$"), colnames(dm.WUF.P.morphonly))]
    dm.BC.P.morphonly <- dm.BC.P.morphonly[-grep(paste0("^",rownames(MF.P)[ROW],"$"), rownames(dm.BC.P.morphonly)),-grep(paste0("^",rownames(MF.P)[ROW],"$"), colnames(dm.BC.P.morphonly))]
  }
}

# Get non-water out
MF.inclWater <- MF
dm.UWUF.inclWater <- dm.UWUF
dm.WUF.inclWater <- dm.WUF
dm.BC.inclWater <- dm.BC
for (ROW in 1:nrow(MF)) {
  typeMorph <- MF[ROW, 'Morph']
  rownameTemp <- rownames(MF)[ROW]
  if (!(typeMorph %in% c('CR',"BL","FB","W"))) {
    MF.inclWater <- MF.inclWater[-grep(paste0("^",rownameTemp,"$"), rownames(MF.inclWater)),]
    
    dm.UWUF.inclWater <- dm.UWUF.inclWater[-grep(paste0("^",rownames(MF)[ROW],"$"), rownames(dm.UWUF.inclWater)),-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(dm.UWUF.inclWater))]
    dm.WUF.inclWater <- dm.WUF.inclWater[-grep(paste0("^",rownames(MF)[ROW],"$"), rownames(dm.WUF.inclWater)),-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(dm.WUF.inclWater))]
    dm.BC.inclWater <- dm.BC.inclWater[-grep(paste0("^",rownames(MF)[ROW],"$"), rownames(dm.BC.inclWater)),-grep(paste0("^",rownames(MF)[ROW],"$"), colnames(dm.BC.inclWater))]
  }
}

MF.P.inclWater <- MF.P
dm.UWUF.P.inclWater <- dm.UWUF.P
dm.WUF.P.inclWater <- dm.WUF.P
dm.BC.P.inclWater <- dm.BC.P
for (ROW in 1:nrow(MF.P)) {
  typeMorph <- MF.P[ROW, 'Morph']
  rownameTemp <- rownames(MF.P)[ROW]
  if (!(typeMorph %in% c('CR',"BL","FB","W"))) {
    MF.P.inclWater <- MF.P.inclWater[-grep(paste0("^",rownameTemp,"$"), rownames(MF.P.inclWater)),]
    
    dm.UWUF.P.inclWater <- dm.UWUF.P.inclWater[-grep(paste0("^",rownames(MF.P)[ROW],"$"), rownames(dm.UWUF.P.inclWater)),-grep(paste0("^",rownames(MF.P)[ROW],"$"), colnames(dm.UWUF.P.inclWater))]
    dm.WUF.P.inclWater <- dm.WUF.P.inclWater[-grep(paste0("^",rownames(MF.P)[ROW],"$"), rownames(dm.WUF.P.inclWater)),-grep(paste0("^",rownames(MF.P)[ROW],"$"), colnames(dm.WUF.P.inclWater))]
    dm.BC.P.inclWater <- dm.BC.P.inclWater[-grep(paste0("^",rownames(MF.P)[ROW],"$"), rownames(dm.BC.P.inclWater)),-grep(paste0("^",rownames(MF.P)[ROW],"$"), colnames(dm.BC.P.inclWater))]
  }
}




system("mkdir BETAPLOTS_H")
system('mkdir ./BETAPLOTS_H/individualtests/')
system("mkdir BETAPLOTS_P")
system('mkdir ./BETAPLOTS_P/individualtests/')
system("mkdir ./BETAPLOTS_LATEX")

############ ******** ALL METRICS ********** #############
metrics <- c("UWUF","WUF","BC")
fullMetrics <- c("Unweighted Unifrac","Weighted Unifrac","Bray-Curtis")

for (metric in metrics) {
  ######## --HAKAI-- ###########
  ### NMDS #####
  assign(paste0("NMDS.",metric,".morphonly"), isoMDS(as.matrix(get(paste0("dm.",metric,".morphonly")))
                                                     # , y = cmdscale(as.matrix(get(paste0("dm.",metric,".morphonly"))),2)
                                                     ))
  assign(paste0("NMDS.",metric,".all"), isoMDS(as.matrix(get(paste0("dm.",metric,".inclWater")))
                                               # , y = cmdscale(as.matrix(get(paste0("dm.",metric,".inclWater"))), 2)
                                               ))
  
  ###### STATS ##########
  MF.morphkeep <- MF.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
  MF.morphkeep$Morph <- factor(MF.morphkeep$Morph, levels = c('CR','BL','FB'))
  MF.morphkeep$Time <- factor(MF.morphkeep$Time, levels = c('20','60','360','720','5760'))
  MF.morphkeep$Type <- factor(MF.morphkeep$Type, levels = c('H'))
  
  assign(paste0("ANOVA.",metric,".morphtime"), adonis(get(paste0("dm.",metric,".morphonly")) ~ Time*Morph, data = MF.morphkeep, by = "margin"))
  capture.output(get(paste0("ANOVA.",metric,".morphtime")), file = paste0("BETAPLOTS_H/adonis_", metric,"_Hakai.txt"))
  
  # Dispersion across time and between morphs
  dist.morphonly <- as.dist(get(paste0("dm.",metric,".inclWater"))[-grep("W", rownames(get(paste0("dm.",metric,".inclWater")))), -grep("W", colnames(get(paste0("dm.",metric,".inclWater"))))])
  MF.incl5760 <- MF.inclWater[-grep("W", rownames(MF.inclWater)),]
  assign(paste0("betadisp.",metric,".time"), betadisper(d = dist.morphonly, group = MF.incl5760$Time))
  assign(paste0("betadisp.",metric,".morph"), betadisper(d = dist.morphonly, group = MF.incl5760$Morph))
  
  capture.output(get(paste0("betadisp.",metric,".time")), file = paste0("BETAPLOTS_H/betadispTime_", metric, "_Hakai.txt"))
  capture.output(get(paste0("betadisp.",metric,".morph")), file = paste0("BETAPLOTS_H/betadispMorph_", metric, "_Hakai.txt"))
  
  
  # Grep values for each timepoint and morph, and then plot dispersion within each
  TP <- levels(factor(MF.inclWater$Time))
  MorphTypes <- levels(factor(MF.morphkeep$Morph))
  
  assign(paste0("FB.BL.",metric,".ALL"), matrix(ncol = 2, nrow = 0))
  assign(paste0("FB.CR.",metric,".ALL"), matrix(ncol = 2, nrow = 0))
  assign(paste0("CR.BL.",metric,".ALL"), matrix(ncol = 2, nrow = 0))
  for (t in 1:length(TP)) {
    # Timepoint
    TPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[t],"$"), MF.inclWater$Time)]
    toremove <- grep("W", TPNames)
    if (length(toremove) > 0) {
      TPNames <- TPNames[-toremove]
    }
    # Between morphs
    assign(paste0("dm.",metric,".",TP[t]), get(paste0("dm.",metric,".inclWater"))[TPNames,TPNames])
    FB.BL.values <- as.vector(as.dist(get(paste0("dm.",metric,".",TP[t]))[grep("FB", rownames(get(paste0("dm.",metric,".",TP[t])))), grep("BL", colnames(get(paste0("dm.",metric,".",TP[t]))))]))
    FB.CR.values <- as.vector(as.dist(get(paste0("dm.",metric,".",TP[t]))[grep("FB", rownames(get(paste0("dm.",metric,".",TP[t])))), grep("CR", colnames(get(paste0("dm.",metric,".",TP[t]))))]))
    CR.BL.values <- as.vector(as.dist(get(paste0("dm.",metric,".",TP[t]))[grep("CR", rownames(get(paste0("dm.",metric,".",TP[t])))), grep("BL", colnames(get(paste0("dm.",metric,".",TP[t]))))]))
    
    # Add to matrix
    assign(paste0("FB.BL.",metric,".ALL"), rbind(get(paste0("FB.BL.",metric,".ALL"))
                     , cbind(as.numeric(FB.BL.values), rep(TP[t], length(FB.BL.values))))
           )
    assign(paste0("FB.CR.",metric,".ALL"), rbind(get(paste0("FB.CR.",metric,".ALL"))
                                                 , cbind(as.numeric(FB.CR.values), rep(TP[t], length(FB.CR.values))))
    )
    assign(paste0("CR.BL.",metric,".ALL"), rbind(get(paste0("CR.BL.",metric,".ALL"))
                                                 , cbind(as.numeric(CR.BL.values), rep(TP[t], length(CR.BL.values))))
    )
    
    # ANOVA
    assign(paste0("MF.",TP[t],".only"), MF.inclWater[match(rownames(get(paste0("dm.",metric,".",TP[t]))),rownames(MF.inclWater)),])
    assign(paste0("ANOVA.",metric,".",TP[t],".only"), adonis(get(paste0("dm.",metric,".",TP[t])) ~ Morph, data = get(paste0("MF.",TP[t],".only")), by = "margin"))
    get(paste0("ANOVA.",metric,".",TP[t],".only"))
  }
  
  # Plot the distances
  assign(paste0("FB.BL.",metric,".ALL.mean"), aggregate(as.numeric(get(paste0("FB.BL.",metric,".ALL"))[,1]), by = list(factor(get(paste0("FB.BL.",metric,".ALL"))[,2], levels = c("20","60","360","720","5760"))), mean))
  assign(paste0("FB.BL.",metric,".ALL.sd"), aggregate(as.numeric(get(paste0("FB.BL.",metric,".ALL"))[,1]), by = list(factor(get(paste0("FB.BL.",metric,".ALL"))[,2], levels = c("20","60","360","720","5760"))), sd))

  assign(paste0("CR.BL.",metric,".ALL.mean"), aggregate(as.numeric(get(paste0("CR.BL.",metric,".ALL"))[,1]), by = list(factor(get(paste0("CR.BL.",metric,".ALL"))[,2], levels = c("20","60","360","720","5760"))), mean))
  assign(paste0("CR.BL.",metric,".ALL.sd"), aggregate(as.numeric(get(paste0("CR.BL.",metric,".ALL"))[,1]), by = list(factor(get(paste0("CR.BL.",metric,".ALL"))[,2], levels = c("20","60","360","720","5760"))), sd))
  
  assign(paste0("FB.CR.",metric,".ALL.mean"), aggregate(as.numeric(get(paste0("FB.CR.",metric,".ALL"))[,1]), by = list(factor(get(paste0("FB.CR.",metric,".ALL"))[,2], levels = c("20","60","360","720","5760"))), mean))
  assign(paste0("FB.CR.",metric,".ALL.sd"), aggregate(as.numeric(get(paste0("FB.CR.",metric,".ALL"))[,1]), by = list(factor(get(paste0("FB.CR.",metric,".ALL"))[,2], levels = c("20","60","360","720","5760"))), sd))

  
  ######## PLOT DISP ###############

  maxy <- ceiling(max(unlist(lapply(c("FB.BL.","CR.BL.","FB.CR."), function(x) {
    max(get(paste0(x,metric,".ALL.mean"))[,2] + get(paste0(x,metric,".ALL.sd"))[,2])
  })))*100)/100
  miny <- floor(min(unlist(lapply(c("FB.BL.","CR.BL.","FB.CR."), function(x) {
    min(get(paste0(x,metric,".ALL.mean"))[,2] - get(paste0(x,metric,".ALL.sd"))[,2])
  })))*100)/100
  
  ylimits <- c(miny, maxy)
  
  # xvalues <- log(FB.BL.UWUF.ALL.mean[,1])
  xvalues <- as.character(get(paste0("FB.BL.",metric,".ALL.mean"))[,1])
  pdf(paste0("BETAPLOTS_H/DispOverTime_H_",metric,".pdf"), pointsize = 14)
  par(fig = c(0,0.8,0,1))
  plot(xvalues, NULL
       , main = "Dispersion of morphologies across time"
       , xlab = "Time"
       , ylab = paste0("Distance (",fullMetrics[which(metrics %in% metric)],")")
       , ylim = ylimits
       , xaxt = 'n')
  axis(side = 1, at = c(1,2,3,4,5), labels = c('20 min','1 h', '6 h', '12 h','4 d'))
  points(get(paste0("FB.BL.",metric,".ALL.mean"))[,2] ~ c(1,2,3,4,5)
         , type = 'l'
         , col = "purple"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5)*0.99
         , x1 = c(1,2,3,4,5)*0.99
         , y0 = c(get(paste0("FB.BL.",metric,".ALL.mean"))[,2] - get(paste0("FB.BL.",metric,".ALL.sd"))[,2]/2)
         , y1 = c(get(paste0("FB.BL.",metric,".ALL.mean"))[,2] + get(paste0("FB.BL.",metric,".ALL.sd"))[,2]/2)
         , angle = 90
         , code = 3
         , length = 0.03)
  points(get(paste0("FB.CR.",metric,".ALL.mean"))[,2] ~ c(1,2,3,4,5)
         , type = 'l'
         , col = "darkgreen"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5)*1
         , x1 = c(1,2,3,4,5)*1
         , y0 = c(get(paste0("FB.CR.",metric,".ALL.mean"))[,2] - get(paste0("FB.CR.",metric,".ALL.sd"))[,2]/2)
         , y1 = c(get(paste0("FB.CR.",metric,".ALL.mean"))[,2] + get(paste0("FB.CR.",metric,".ALL.sd"))[,2]/2)
         , angle = 90
         , code = 3
         , length = 0.03)
  points(get(paste0("CR.BL.",metric,".ALL.mean"))[,2] ~ c(1,2,3,4,5)
         , type = 'l'
         , col = "grey"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5)*1.01
         , x1 = c(1,2,3,4,5)*1.01
         , y0 = c(get(paste0("CR.BL.",metric,".ALL.mean"))[,2] - get(paste0("CR.BL.",metric,".ALL.sd"))[,2]/2)
         , y1 = c(get(paste0("CR.BL.",metric,".ALL.mean"))[,2] + get(paste0("CR.BL.",metric,".ALL.sd"))[,2]/2)
         , angle = 90
         , code = 3
         , length = 0.03)
  par(fig = c(0.7,1,0,1), mar = c(5,0,5,0), new = TRUE)
  plot(0,0
       , pch = ""
       , xlab = ""
       , ylab = ""
       , xaxt = "n"
       , yaxt = "n"
       , bty = "n")
  legend("center"
         , legend = c("FB:CR", "FB:BL", "CR:BL")
         , lty = 1
         , col = c("darkgreen","purple","grey")
         , lwd = 2)
  dev.off()
  
  ####### PLOT MORPH ############
  # NO 5760 UWUF
  
  # MAKE POLYGONS for plotting
  for (morph in c("FB","BL","CR")) {
    assign(paste0("NMDS.",metric,".",morph), get(paste0("NMDS.",metric,".morphonly"))$points[grep(paste0(morph), rownames(get(paste0("NMDS.",metric,".morphonly"))$points)),])
    assign(paste0("NMDS.",metric,".",morph,".chull"), chull(get(paste0("NMDS.",metric,".",morph))))
    assign(paste0("NMDS.",metric,".",morph,".chull"), c(get(paste0("NMDS.",metric,".",morph,".chull")), get(paste0("NMDS.",metric,".",morph,".chull"))[1]))
    
  }
  
  MorphColours <- c("darkorchid4","dodgerblue","salmon") 
  
  pdf(paste0("BETAPLOTS_H/NMDS_H_",metric,"_Morph.pdf"), pointsize = 14)
  par(fig = c(0,0.8,0,1))
  plot(get(paste0("NMDS.",metric,".morphonly"))$points
       , main = "NMDS of Artificial Seaweed Shapes"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.morphkeep$Morph)]
       , sub = round(get(paste0("NMDS.",metric,".morphonly"))$stress/100,2)
       , xlab = "NMDS 1"
       , ylab = "NMDS 2"
       , cex = 2
  )
  lines(get(paste0("NMDS.",metric,".CR"))[get(paste0("NMDS.",metric,".CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".BL"))[get(paste0("NMDS.",metric,".BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".FB"))[get(paste0("NMDS.",metric,".FB.chull")),]
        , col = MorphColours[3])
  par(fig = c(0.75,1,0,1), mar = c(5,0,5,0), new = TRUE)
  plot(0,0
       , pch = ''
       , xlab = ''
       , ylab = ''
       , xaxt = 'n'
       , yaxt = 'n'
       , bty = 'n')
  legend("topleft"
         , pch = 21
         , col = "black"
         , pt.bg = MorphColours
         , legend = c("Crustose","Bladed","Finely Branched")
         , cex = 0.8
  )
  dev.off()
  
  ####### PLOT TIME ############
  TimeColours <- c("grey","lightblue","blue","darkblue")
  
  pdf(paste0("BETAPLOTS_H/NMDS_H_",metric,"_Time.pdf"), pointsize = 14)
  par(fig = c(0,0.8,0,1))
  plot(get(paste0("NMDS.",metric,".morphonly"))$points
       , main = "NMDS of Artificial Seaweed Shapes"
       , pch = 21
       , col = "black"
       , bg = TimeColours[factor(MF.morphkeep$Time)]
       , sub = round(get(paste0("NMDS.",metric,".morphonly"))$stress/100,2)
       , xlab = "NMDS 1"
       , ylab = "NMDS 2"
       , cex = 2
  )
  # lines(NMDS.UWUF.CR[NMDS.UWUF.CR.chull,]
  #       , col = "red")
  # lines(NMDS.UWUF.BL[NMDS.UWUF.BL.chull,]
  #       , col = "magenta")
  # lines(NMDS.UWUF.FB[NMDS.UWUF.FB.chull,]
  #       , col = "blue")
  par(fig = c(0.75,1,0,1), mar = c(5,0,5,0), new = TRUE)
  plot(0,0
       , pch = ''
       , xlab = ''
       , ylab = ''
       , xaxt = 'n'
       , yaxt = 'n'
       , bty = 'n')
  legend("topleft"
         , pch = 21
         , col = "black"
         , pt.bg = TimeColours
         , legend = c("20 minutes","1 hour","6 hours","12 hours")
  )
  dev.off()
  
  
  
  # # TIME-- ONE at a time
  # 
  # for (t in TP) {
  #   # Per time point
  #   assign(paste0("NMDS.",metric,".",t,".only"), isoMDS(as.matrix(get(paste0("dm.",metric,".",t)))
  #                                                       # , y = cmdscale(as.matrix(get(paste0("dm.",metric,".",t))), 2)
  #                                                       ))
  #   for (morph in c("FB","BL","CR")) {
  #     assign(paste0("NMDS.",metric,".",t,".only.",morph), get(paste0("NMDS.",metric,".",t,".only"))$points[grep(morph, rownames(get(paste0("NMDS.",metric,".",t,".only"))$points)),])
  #     assign(paste0("NMDS.",metric,".",t,".only.",morph,".chull"), chull(get(paste0("NMDS.",metric,".",t,".only.",morph))))
  #     assign(paste0("NMDS.",metric,".",t,".only.",morph,".chull"), c(get(paste0("NMDS.",metric,".",t,".only.",morph,".chull")), get(paste0("NMDS.",metric,".",t,".only.",morph,".chull"))[1]))
  #   }
  # }
  
  # TIME-- ONE at a time
  
  for (t in TP) {
    # Per time point
    if (t == "5760") {
      assign(paste0("NMDS.",metric,".",t,".only.FULL"), isoMDS(as.matrix(get(paste0("dm.",metric,".",t))), k = 2
                                                     ))
      assign(paste0("NMDS.",metric,".",t,".only"), get(paste0("NMDS.",metric,".",t,".only.FULL"))$points)
      assign(paste0("NMDS.",metric,".",t,".only.stress"), get(paste0("NMDS.",metric,".",t,".only.FULL"))$stress)
    } else {
      assign(paste0("NMDS.",metric,".",t,".only"), get(paste0("NMDS.",metric,".morphonly"))$points[grep(paste0("-",t,"-"), rownames(get(paste0("NMDS.",metric,".morphonly"))$points)),]
      )
    }

    for (morph in c("FB","BL","CR")) {
      assign(paste0("NMDS.",metric,".",t,".only.",morph), get(paste0("NMDS.",metric,".",t,".only"))[grep(morph, rownames(get(paste0("NMDS.",metric,".",t,".only")))),])
      assign(paste0("NMDS.",metric,".",t,".only.",morph,".chull"), chull(get(paste0("NMDS.",metric,".",t,".only.",morph))))
      assign(paste0("NMDS.",metric,".",t,".only.",morph,".chull"), c(get(paste0("NMDS.",metric,".",t,".only.",morph,".chull")), get(paste0("NMDS.",metric,".",t,".only.",morph,".chull"))[1]))
    }
  }
  
  ############ COMBO DISP BETA ################
  # Disp and beta through time combined
  xvalues <- as.character(get(paste0("FB.BL.",metric,".ALL.mean"))[,1])
  
  # Change factor of individual plots
  MF.20.only$Morph <- factor(MF.20.only$Morph, levels = c("CR","BL","FB"))
  MF.60.only$Morph <- factor(MF.60.only$Morph, levels = c("CR","BL","FB"))
  MF.360.only$Morph <- factor(MF.360.only$Morph, levels = c("CR","BL","FB"))
  MF.720.only$Morph <- factor(MF.720.only$Morph, levels = c("CR","BL","FB"))
  MF.5760.only$Morph <- factor(MF.5760.only$Morph, levels = c("CR","BL","FB"))
  
  # Get significance for each timepoint
  listSig <- c()
  for (t in TP) {
    tempSig <- get(paste0("ANOVA.",metric,".",t,".only"))$aov.tab$`Pr(>F)`[1]
    if (as.numeric(tempSig) <= 0.001) {
      listSig <- c(listSig, "***")
    } else if (as.numeric(tempSig) <= 0.01) {
      listSig <- c(listSig, "**")
    } else if (as.numeric(tempSig) <= 0.05) {
      listSig <- c(listSig,"*")
    } else {
      listSig <- c(listSig, " ")
    }
  }
  names(listSig) <- TP
  
  pdf(paste0("BETAPLOTS_H/COMBO_H_dispbeta_",metric,".pdf"), pointsize = 14, width = 10, height = 6)
  par(fig = c(0,0.8,0.3,1))
  plot(xvalues, NULL
       , main = "Dispersion of morphologies across time"
       , xlab = ''
       , ylab = paste0("Distance (",fullMetrics[which(metrics %in% metric)],")")
       , ylim = ylimits
       , xaxt = 'n')
  axis(side = 1
       , at = c(1,2,3,4,5)
       , labels = c("20 min","1 h","6 h","12 h","4 d")
       , las = 2)
  points(get(paste0("FB.BL.",metric,".ALL.mean"))[,2] ~ c(1,2,3,4,5)
         , type = 'l'
         , col = "purple"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5)*0.99
         , x1 = c(1,2,3,4,5)*0.99
         , y0 = c(get(paste0("FB.BL.",metric,".ALL.mean"))[,2] - get(paste0("FB.BL.",metric,".ALL.sd"))[,2])
         , y1 = c(get(paste0("FB.BL.",metric,".ALL.mean"))[,2] + get(paste0("FB.BL.",metric,".ALL.sd"))[,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  points(get(paste0("FB.CR.",metric,".ALL.mean"))[,2] ~ c(1,2,3,4,5)
         , type = 'l'
         , col = "darkgreen"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5)*1
         , x1 = c(1,2,3,4,5)*1
         , y0 = c(get(paste0("FB.CR.",metric,".ALL.mean"))[,2] - get(paste0("FB.CR.",metric,".ALL.sd"))[,2])
         , y1 = c(get(paste0("FB.CR.",metric,".ALL.mean"))[,2] + get(paste0("FB.CR.",metric,".ALL.sd"))[,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  points(get(paste0("CR.BL.",metric,".ALL.mean"))[,2] ~ c(1,2,3,4,5)
         , type = 'l'
         , col = "grey"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5)*1.01
         , x1 = c(1,2,3,4,5)*1.01
         , y0 = c(get(paste0("CR.BL.",metric,".ALL.mean"))[,2] - get(paste0("CR.BL.",metric,".ALL.sd"))[,2])
         , y1 = c(get(paste0("CR.BL.",metric,".ALL.mean"))[,2] + get(paste0("CR.BL.",metric,".ALL.sd"))[,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  title(xlab = paste0("Stress: ", round(get(paste0("NMDS.",metric,".morphonly"))$stress/100,2))
        , line = 3)
  par(fig = c(0.7,1,0.5,1), mar = c(2,0,0,0),oma = c(0,0,0,0), new = TRUE)
  plot(0,0
       , pch = ""
       , xlab = ""
       , ylab = ""
       , xaxt = "n"
       , yaxt = "n"
       , bty = "n")
  legend("center"
         , legend = c("FB:CR", "FB:BL", "CR:BL")
         , lty = 1
         , col = c("darkgreen","purple","grey")
         , lwd = 2 
  )
  # EACH GETS 0.15 SPACE TOTAL;
  par(oma = c(1,0.1,1,0.1), mar = c(2,0,3,0), fig = c(0.049,0.1932,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".20.only"))#$points
       # , main = "20 Minutes"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.20.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
  )
  title(xlab = paste0(listSig["20"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ",round(get(paste0("NMDS.",metric,".20.only"))$stress)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)
  # title(xlab = paste0("Df = ", get(paste0("ANOVA.",metric,".20.only"))$aov.tab[1]$Df[1],",",get(paste0("ANOVA.",metric,".20.only"))$aov.tab[1]$Df[3]) 
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.20.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.20.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  
  lines(get(paste0("NMDS.",metric,".20.only.CR"))[get(paste0("NMDS.",metric,".20.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".20.only.BL"))[get(paste0("NMDS.",metric,".20.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".20.only.FB"))[get(paste0("NMDS.",metric,".20.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(2,0,3,0), fig = c(0.2032,0.3474,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".60.only"))#$points
       # , main = "1 Hour"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.60.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
  )
  title(xlab = paste0(listSig["60"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ",round(get(paste0("NMDS.",metric,".60.only"))$stress)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)  
  # title(xlab = paste0("Df = ", ANOVA.UWUF.60.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.60.only$aov.tab[1]$Df[3]) 
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.60.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.60.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".60.only.CR"))[get(paste0("NMDS.",metric,".60.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".60.only.BL"))[get(paste0("NMDS.",metric,".60.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".60.only.FB"))[get(paste0("NMDS.",metric,".60.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(2,0,3,0), fig = c(0.3574,0.5016,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".360.only"))#$points
       # , main = "6 Hours"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.360.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
  )
  title(xlab = paste0(listSig["360"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ",round(get(paste0("NMDS.",metric,".360.only"))$stress)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)  
  # title(xlab = paste0("Df = ", ANOVA.UWUF.360.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.360.only$aov.tab[1]$Df[3]) 
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.360.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.360.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".360.only.CR"))[get(paste0("NMDS.",metric,".360.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".360.only.BL"))[get(paste0("NMDS.",metric,".360.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".360.only.FB"))[get(paste0("NMDS.",metric,".360.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(2,0,3,0), fig = c(0.5116,0.6558,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".720.only"))#$points
       # , main = "12 Hours"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.720.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
  )
  title(xlab = paste0(listSig["720"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ",round(get(paste0("NMDS.",metric,".720.only"))$stress)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)
  # title(xlab = paste0("Df = ", ANOVA.UWUF.720.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.720.only$aov.tab[1]$Df[3]) 
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.720.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.720.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".720.only.CR"))[get(paste0("NMDS.",metric,".720.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".720.only.BL"))[get(paste0("NMDS.",metric,".720.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".720.only.FB"))[get(paste0("NMDS.",metric,".720.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(2,0,3,0), fig = c(0.6658,0.81,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".5760.only"))#$points
       # , main = "12 Hours"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.5760.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
   )
  title(xlab = paste0(listSig["5760"])
        , line = 1.25
        , cex.lab = 3)
  title(sub = paste0("Stress: ",round(get(paste0("NMDS.",metric,".5760.only.stress"))/100,2))
        , line = 1
        , cex.sub = 0.75)
  # title(xlab = paste0("Df = ", ANOVA.UWUF.5760.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.5760.only$aov.tab[1]$Df[3]) 
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.5760.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.5760.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".5760.only.CR"))[get(paste0("NMDS.",metric,".5760.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".5760.only.BL"))[get(paste0("NMDS.",metric,".5760.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".5760.only.FB"))[get(paste0("NMDS.",metric,".5760.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(2,0,2,0), fig = c(0.775,1,0,0.4), new = TRUE)
  plot(0,0
       , pch = ''
       , bty = 'n'
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = ''
       , ylab = ''
  )
  legend("center"
         , legend = levels(MF.morphkeep$Morph)
         , pch = 21
         , col = "black"
         , pt.bg = MorphColours
  )
  #STOP
  dev.off()
  
  ########### BETADISP#############
  assign(paste0("betadisp.H.",metric,".FB"), get(paste0("betadisp.",metric,".time"))$distances[grep("FB", names(get(paste0("betadisp.",metric,".time"))$distances))])
  assign(paste0("betadisp.H.",metric,".BL"), get(paste0("betadisp.",metric,".time"))$distances[grep("BL", names(get(paste0("betadisp.",metric,".time"))$distances))])
  assign(paste0("betadisp.H.",metric,".CR"), get(paste0("betadisp.",metric,".time"))$distances[grep("CR", names(get(paste0("betadisp.",metric,".time"))$distances))])
  
  MF.H.FB <- MF.incl5760[grep("FB", rownames(MF.incl5760)),]
  MF.H.BL <- MF.incl5760[grep("BL", rownames(MF.incl5760)),]
  MF.H.CR <- MF.incl5760[grep("CR", rownames(MF.incl5760)),]
  
  betadisp.FB.H.agg <- aggregate(get(paste0("betadisp.H.",metric,".FB")), by = list(MF.H.FB$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
  betadisp.BL.H.agg <- aggregate(get(paste0("betadisp.H.",metric,".BL")), by = list(MF.H.BL$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
  betadisp.CR.H.agg <- aggregate(get(paste0("betadisp.H.",metric,".CR")), by = list(MF.H.CR$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
  
  assign(paste0("betadisp.",metric,".time.forstat"), cbind(Distance = get(paste0("betadisp.",metric,".time"))$distances, MF.incl5760[unlist(lapply(names(get(paste0("betadisp.",metric,".time"))$distances), function(x) {grep(paste0("^",x,"$"), rownames(MF.incl5760))})),]))
  assign(paste0("ANOVA.betadisp.",metric), anova(lm(Distance ~ Time*Morph, data = get(paste0("betadisp.",metric,".time.forstat")))))
  capture.output(get(paste0("ANOVA.betadisp.",metric)), file = paste0("./BETAPLOTS_H/ANOVA.betadisp.",metric,".txt"))
  
  xvalues <- c("20","60","360","720","5760")
  
  maxy <- ceiling(max(unlist(lapply(c(".FB.",".BL.",".CR."), function(x) {
    max(get(paste0("betadisp",x,"H.agg"))[,2][,1] + get(paste0("betadisp",x,"H.agg"))[,2][,2])
  })))*100)/100
  miny <- round(min(unlist(lapply(c(".FB.",".BL.",".CR."), function(x) {
    min(get(paste0("betadisp",x,"H.agg"))[,2][,1] - get(paste0("betadisp",x,"H.agg"))[,2][,2])
  })))*100)/100
  
  ylimits <- c(miny, maxy)
  pdf(paste0("./BETAPLOTS_H/BetaDisp_H_",metric,"_eachmorph.pdf"),pointsize = 14)
  plot(xvalues, NULL
       , main = "Dispersion of morphologies across time"
       , xlab = 'Time'
       , ylab = paste0("Distance (",fullMetrics[which(metrics %in% metric)],")")
       , ylim = ylimits
       , xaxt = 'n')
  axis(side = 1
       , at = c(1,2,3,4,5)
       , labels = c("20 min","1 h","6 h","12 h","4 d")
       , las = 2)
  points(betadisp.FB.H.agg[,2][,1] ~ c(1,2,3,4,5)
         , type = 'l'
         , col = "salmon"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5)*0.99
         , x1 = c(1,2,3,4,5)*0.99
         , y0 = c(betadisp.FB.H.agg[,2][,1] - betadisp.FB.H.agg[,2][,2])
         , y1 = c(betadisp.FB.H.agg[,2][,1] + betadisp.FB.H.agg[,2][,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  points(betadisp.BL.H.agg[,2][,1] ~ c(1,2,3,4,5)
         , type = 'l'
         , col = "dodgerblue"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5)*1
         , x1 = c(1,2,3,4,5)*1
         , y0 = c(betadisp.BL.H.agg[,2][,1] - betadisp.BL.H.agg[,2][,2])
         , y1 = c(betadisp.BL.H.agg[,2][,1] + betadisp.BL.H.agg[,2][,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  points(betadisp.CR.H.agg[,2][,1] ~ c(1,2,3,4,5)
         , type = 'l'
         , col = "darkorchid4"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5)*1.01
         , x1 = c(1,2,3,4,5)*1.01
         , y0 = c(betadisp.CR.H.agg[,2][,1] - betadisp.CR.H.agg[,2][,2])
         , y1 = c(betadisp.CR.H.agg[,2][,1] + betadisp.CR.H.agg[,2][,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  par(fig = c(0.7,1,0.5,1), mar = c(2,0,0,0),oma = c(0,0,0,0), new = TRUE)
  plot(0,0
       , pch = ""
       , xlab = ""
       , ylab = ""
       , xaxt = "n"
       , yaxt = "n"
       , bty = "n")
  legend("center"
         , legend = c("FB", "BL", "CR")
         , lty = 1
         , col = c("salmon","dodgerblue","darkorchid4")
         , lwd = 2 
  )
  dev.off()
  
  
  ############ PLOT 5760 ################

  
  pdf(paste0("BETAPLOTS_H/NMDS_",metric,"_5760Only.pdf"), pointsize = 14)
  par(fig = c(0,0.75,0,1))
  plot(get(paste0("NMDS.",metric,".5760.only"))#$points
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.5760.only$Morph)]
       , xlab = "NMDS 1"
       , ylab = "NMDS 2"
       , sub = paste0("Stress: ", round(get(paste0("NMDS.",metric,".5760.only.stress"))/100,2))
       , cex = 1.5
       , main = "NMDS of community composition (4 days)")
  lines(get(paste0("NMDS.",metric,".5760.only.CR"))[get(paste0("NMDS.",metric,".5760.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".5760.only.BL"))[get(paste0("NMDS.",metric,".5760.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".5760.only.FB"))[get(paste0("NMDS.",metric,".5760.only.FB.chull")),]
        , col = MorphColours[3])
  par(fig = c(0.7,1,0,1), mar = c(0,0,0,0), new = TRUE)
  plot(0,0
       , xlab = ""
       , ylab = ""
       , xaxt = "n"
       , yaxt = "n"
       , bty = "n"
       , pch = ""
  )
  legend("center"
         , legend = c("Crust","Blade","Finely Br.")
         , pch = 21
         , col= "black"
         , pt.bg = c(MorphColours[1],MorphColours[2],MorphColours[3]))
  dev.off()
  
  
  ######## --PM-- ###########
  ### NMDS #####

  assign(paste0("NMDS.",metric,".P.morphonly"), isoMDS(as.matrix(get(paste0("dm.",metric,".P.morphonly")))
                                                       # , y = cmdscale(as.matrix(get(paste0("dm.",metric,".P.morphonly"))), 2)
                                                       ))
  assign(paste0("NMDS.",metric,".P.all"), isoMDS(as.matrix(get(paste0("dm.",metric,".P.inclWater")))
                                                 # , y = cmdscale(as.matrix(get(paste0("dm.",metric,".P.inclWater"))), 2)
                                                 ))
  
  ###### STATS ##########
  MF.P.morphkeep <- MF.P.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
  MF.P.morphkeep$Morph <- factor(MF.P.morphkeep$Morph, levels = c('CR','BL','FB'))
  MF.P.morphkeep$Time <- factor(MF.P.morphkeep$Time, levels = c('20','60','180','360','720','1440'))
  MF.P.morphkeep$Type <- factor(MF.P.morphkeep$Type, levels = c('P'))
  
  assign(paste0("ANOVA.",metric,".P.morphtime"), adonis(get(paste0("dm.",metric,".P.morphonly")) ~ Time*Morph, data = MF.P.morphkeep, by = "margin"))
  capture.output(get(paste0("ANOVA.",metric,".P.morphtime")), file = paste0("BETAPLOTS_P/adonis_", metric,"_PM.txt"))

  # Dispersion across time and between morphs
  dist.P.morphonly <- as.dist(get(paste0("dm.",metric,".P.morphonly")))
  assign(paste0("betadisp.",metric,".P.time"), betadisper(d = dist.P.morphonly, group = MF.P.morphkeep$Time))
  assign(paste0("betadisp.",metric,".P.morph"), betadisper(d = dist.P.morphonly, group = MF.P.morphkeep$Morph))
  
  capture.output(get(paste0("betadisp.",metric,".P.time")), file = paste0("BETAPLOTS_P/betadispTime_", metric, "_PM.txt"))
  capture.output(get(paste0("betadisp.",metric,".P.morph")), file = paste0("BETAPLOTS_P/betadispMorph_", metric, "_PM.txt"))
  
  
  # Grep values for each timepoint and morph, and then plot dispersion within each
  TP <- levels(factor(MF.P.morphkeep$Time))
  MorphTypes <- levels(factor(MF.P.morphkeep$Morph))
  
  
  assign(paste0("FB.BL.",metric,".P.ALL"), matrix(ncol = 2, nrow = 0))
  assign(paste0("FB.CR.",metric,".P.ALL"), matrix(ncol = 2, nrow = 0))
  assign(paste0("CR.BL.",metric,".P.ALL"), matrix(ncol = 2, nrow = 0))
  for (t in 1:length(TP)) {
    # Timepoint
    TPNames <- rownames(MF.P.inclWater)[grep(paste0("^",TP[t],"$"), MF.P.inclWater$Time)]
    toremove <- grep("W", TPNames)
    if (length(toremove) > 0) {
      TPNames <- TPNames[-toremove]
    }
    # Between morphs
    assign(paste0("dm.",metric,".P.",TP[t]), get(paste0("dm.",metric,".P.inclWater"))[TPNames,TPNames])
    FB.BL.values <- as.vector(as.dist(get(paste0("dm.",metric,".P.",TP[t]))[grep("FB", rownames(get(paste0("dm.",metric,".P.",TP[t])))), grep("BL", colnames(get(paste0("dm.",metric,".P.",TP[t]))))]))
    FB.CR.values <- as.vector(as.dist(get(paste0("dm.",metric,".P.",TP[t]))[grep("FB", rownames(get(paste0("dm.",metric,".P.",TP[t])))), grep("CR", colnames(get(paste0("dm.",metric,".P.",TP[t]))))]))
    CR.BL.values <- as.vector(as.dist(get(paste0("dm.",metric,".P.",TP[t]))[grep("CR", rownames(get(paste0("dm.",metric,".P.",TP[t])))), grep("BL", colnames(get(paste0("dm.",metric,".P.",TP[t]))))]))
    
    # Add to matrix
    assign(paste0("FB.BL.",metric,".P.ALL"), rbind(get(paste0("FB.BL.",metric,".P.ALL"))
                                                 , cbind(as.numeric(FB.BL.values), rep(TP[t], length(FB.BL.values))))
    )
    assign(paste0("FB.CR.",metric,".P.ALL"), rbind(get(paste0("FB.CR.",metric,".P.ALL"))
                                                 , cbind(as.numeric(FB.CR.values), rep(TP[t], length(FB.CR.values))))
    )
    assign(paste0("CR.BL.",metric,".P.ALL"), rbind(get(paste0("CR.BL.",metric,".P.ALL"))
                                                 , cbind(as.numeric(CR.BL.values), rep(TP[t], length(CR.BL.values))))
    )
    
    # ANOVA
    assign(paste0("MF.P.",TP[t],".only"), MF.P.inclWater[match(rownames(get(paste0("dm.",metric,".P.",TP[t]))),rownames(MF.P.inclWater)),])
    assign(paste0("ANOVA.",metric,".P.",TP[t],".only"), adonis(get(paste0("dm.",metric,".P.",TP[t])) ~ Morph, data = get(paste0("MF.P.",TP[t],".only")), by = "margin"))
  }
  
  
  # Plot the distances
  assign(paste0("FB.BL.",metric,".P.ALL.mean"), aggregate(as.numeric(get(paste0("FB.BL.",metric,".P.ALL"))[,1]), by = list(factor(get(paste0("FB.BL.",metric,".P.ALL"))[,2], levels = c("20","60","180","360","720","1440"))), mean))
  assign(paste0("FB.BL.",metric,".P.ALL.sd"), aggregate(as.numeric(get(paste0("FB.BL.",metric,".P.ALL"))[,1]), by = list(factor(get(paste0("FB.BL.",metric,".P.ALL"))[,2], levels = c("20","60","180","360","720","1440"))), sd))
  
  assign(paste0("CR.BL.",metric,".P.ALL.mean"), aggregate(as.numeric(get(paste0("CR.BL.",metric,".P.ALL"))[,1]), by = list(factor(get(paste0("CR.BL.",metric,".P.ALL"))[,2], levels = c("20","60","180","360","720","1440"))), mean))
  assign(paste0("CR.BL.",metric,".P.ALL.sd"), aggregate(as.numeric(get(paste0("CR.BL.",metric,".P.ALL"))[,1]), by = list(factor(get(paste0("CR.BL.",metric,".P.ALL"))[,2], levels = c("20","60","180","360","720","1440"))), sd))
  
  assign(paste0("FB.CR.",metric,".P.ALL.mean"), aggregate(as.numeric(get(paste0("FB.CR.",metric,".P.ALL"))[,1]), by = list(factor(get(paste0("FB.CR.",metric,".P.ALL"))[,2], levels = c("20","60","180","360","720","1440"))), mean))
  assign(paste0("FB.CR.",metric,".P.ALL.sd"), aggregate(as.numeric(get(paste0("FB.CR.",metric,".P.ALL"))[,1]), by = list(factor(get(paste0("FB.CR.",metric,".P.ALL"))[,2], levels = c("20","60","180","360","720","1440"))), sd))
  
  get(paste0("FB.CR.",metric,".P.ALL.sd"))
  ######## PLOT DISP ###############

  maxy <- ceiling(max(unlist(lapply(c("FB.BL.","CR.BL.","FB.CR."), function(x) {
    naList <- which(is.na(get(paste0(x,metric,".P.ALL.sd"))[,2]))
    sdTemp <- get(paste0(x,metric,".P.ALL.sd"))[,2]
    if (length(naList) >=1) {
      sdTemp[naList] <- 0
    }
    max(get(paste0(x,metric,".P.ALL.mean"))[,2] + sdTemp, na.rm = TRUE)
  })))*100)/100
  miny <- floor(min(unlist(lapply(c("FB.BL.","CR.BL.","FB.CR."), function(x) {
    naList <- which(is.na(get(paste0(x,metric,".P.ALL.sd"))[,2]))
    sdTemp <- get(paste0(x,metric,".P.ALL.sd"))[,2]
    if (length(naList) >=1) {
      sdTemp[naList] <- 0
    }
    min(get(paste0(x,metric,".P.ALL.mean"))[,2] - sdTemp, na.rm = TRUE)
  })))*100)/100
  
  ylimits <- c(miny, maxy)
  
  # xvalues <- log(FB.BL.UWUF.P.ALL.mean[,1])
  xvalues <- as.character(get(paste0("FB.BL.",metric,".P.ALL.mean"))[,1])
  pdf(paste0("BETAPLOTS_P/DispOverTime_P_",metric,".pdf"), pointsize = 14)
  par(fig = c(0,0.8,0,1))
  plot(xvalues, NULL
       , main = "Dispersion of morphologies across time"
       , xlab = "Time"
       , ylab = paste0("Distance (",fullMetrics[which(metrics %in% metric)],")")
       , ylim = ylimits
       , xaxt = 'n')
  axis(side = 1, at = c(1,2,3,4,5,6), labels = c('20 min','1 h','3 h','6 h', '12 h', '24 h'))
  points(get(paste0("FB.BL.",metric,".P.ALL.mean"))[,2] ~ c(1,2,3,4,5,6)
         , type = 'l'
         , col = "purple"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5,6)*0.99
         , x1 = c(1,2,3,4,5,6)*0.99
         , y0 = c(get(paste0("FB.BL.",metric,".P.ALL.mean"))[,2] - get(paste0("FB.BL.",metric,".P.ALL.sd"))[,2]/2)
         , y1 = c(get(paste0("FB.BL.",metric,".P.ALL.mean"))[,2] + get(paste0("FB.BL.",metric,".P.ALL.sd"))[,2]/2)
         , angle = 90
         , code = 3
         , length = 0.03)
  points(get(paste0("FB.CR.",metric,".P.ALL.mean"))[,2] ~ c(1,2,3,4,5,6)
         , type = 'l'
         , col = "darkgreen"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5,6)*1
         , x1 = c(1,2,3,4,5,6)*1
         , y0 = c(get(paste0("FB.CR.",metric,".P.ALL.mean"))[,2] - get(paste0("FB.CR.",metric,".P.ALL.sd"))[,2]/2)
         , y1 = c(get(paste0("FB.CR.",metric,".P.ALL.mean"))[,2] + get(paste0("FB.CR.",metric,".P.ALL.sd"))[,2]/2)
         , angle = 90
         , code = 3
         , length = 0.03)
  points(get(paste0("CR.BL.",metric,".P.ALL.mean"))[,2] ~ c(1,2,3,4,5,6)
         , type = 'l'
         , col = "grey"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5,6)*1.01
         , x1 = c(1,2,3,4,5,6)*1.01
         , y0 = c(get(paste0("CR.BL.",metric,".P.ALL.mean"))[,2] - get(paste0("CR.BL.",metric,".P.ALL.sd"))[,2]/2)
         , y1 = c(get(paste0("CR.BL.",metric,".P.ALL.mean"))[,2] + get(paste0("CR.BL.",metric,".P.ALL.sd"))[,2]/2)
         , angle = 90
         , code = 3
         , length = 0.03)
  par(fig = c(0.7,1,0,1), mar = c(5,0,5,0), new = TRUE)
  plot(0,0
       , pch = ""
       , xlab = ""
       , ylab = ""
       , xaxt = "n"
       , yaxt = "n"
       , bty = "n")
  legend("center"
         , legend = c("FB:CR", "FB:BL", "CR:BL")
         , lty = 1
         , col = c("darkgreen","purple","grey")
         , lwd = 2)
  dev.off()
  
  ####### PLOT MORPH ############

  
  # MAKE POLYGONS for plotting
  for (morph in c("CR","BL","FB")) {
    assign(paste0("NMDS.",metric,".P.",morph), get(paste0("NMDS.",metric,".P.morphonly"))$points[grep(paste0(morph), rownames(get(paste0("NMDS.",metric,".P.morphonly"))$points)),])
    
    assign(paste0("NMDS.",metric,".P.",morph,".chull"), chull(get(paste0("NMDS.",metric,".P.",morph))))
    assign(paste0("NMDS.",metric,".P.",morph,".chull"), c(get(paste0("NMDS.",metric,".P.",morph,".chull")), get(paste0("NMDS.",metric,".P.",morph,".chull"))[1]))
    
  }
  
  
  MorphColours <- c("darkorchid4","dodgerblue","salmon") 
  
  pdf(paste0("BETAPLOTS_P/NMDS_P_",metric,"_Morph.pdf"), pointsize = 14)
  par(fig = c(0,0.8,0,1))
  plot(get(paste0("NMDS.",metric,".P.morphonly"))$points
       , main = "NMDS of Artificial Seaweed Shapes"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.P.morphkeep$Morph)]
       , sub = round(get(paste0("NMDS.",metric,".P.morphonly"))$stress/100,2)
       , xlab = "NMDS 1"
       , ylab = "NMDS 2"
       , cex = 2
  )
  
  lines(get(paste0("NMDS.",metric,".P.CR"))[get(paste0("NMDS.",metric,".P.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".P.BL"))[get(paste0("NMDS.",metric,".P.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".P.FB"))[get(paste0("NMDS.",metric,".P.FB.chull")),]
        , col = MorphColours[3])
  par(fig = c(0.75,1,0,1), mar = c(5,0,5,0), new = TRUE)
  plot(0,0
       , pch = ''
       , xlab = ''
       , ylab = ''
       , xaxt = 'n'
       , yaxt = 'n'
       , bty = 'n')
  legend("topleft"
         , pch = 21
         , col = "black"
         , pt.bg = MorphColours
         , legend = c("Crustose","Bladed","Finely Branched")
         , cex = 0.8
  )
  dev.off()
  
  ####### PLOT TIME ############
  TimeColours <- c("grey","lightblue","deepskyblue","blue","darkblue", "purple")
  
  pdf(paste0("BETAPLOTS_P/NMDS_P_",metric,"_Time.pdf"), pointsize = 14)
  par(fig = c(0,0.8,0,1))
  plot(get(paste0("NMDS.",metric,".P.morphonly"))$points[,1]
       , main = "NMDS of Artificial Seaweed Shapes"
       , pch = 21
       , col = "black"
       , bg = TimeColours[factor(MF.P.morphkeep$Time)]
       , sub = round(get(paste0("NMDS.",metric,".P.morphonly"))$stress/100,2)
       , xlab = "NMDS 1"
       , ylab = "NMDS 2"
       , cex = 2
  )
  # lines(NMDS.UWUF.P.CR[NMDS.UWUF.P.CR.chull,]
  #       , col = "red")
  # lines(NMDS.UWUF.P.BL[NMDS.UWUF.P.BL.chull,]
  #       , col = "magenta")
  # lines(NMDS.UWUF.P.FB[NMDS.UWUF.P.FB.chull,]
  #       , col = "blue")
  par(fig = c(0.75,1,0,1), mar = c(5,0,5,0), new = TRUE)
  plot(0,0
       , pch = ''
       , xlab = ''
       , ylab = ''
       , xaxt = 'n'
       , yaxt = 'n'
       , bty = 'n')
  legend("topleft"
         , pch = 21
         , col = "black"
         , pt.bg = TimeColours
         , legend = c("20 minutes","1 hour","3 hours","6 hours","12 hours", "24 hours")
  )
  dev.off()

  # # TIME-- ONE at a time
  # for (t in TP) {
  #   # Per time point
  #   dm.temp <- as.dist(get(paste0("dm.",metric,".P.",t)))
  #   assign(paste0("NMDS.",metric,".P.",t,".only"), cmdscale(dm.temp, k = 2))
  #   for (morph in c("FB","BL","CR")) {
  #     assign(paste0("NMDS.",metric,".P.",t,".only.",morph), get(paste0("NMDS.",metric,".P.",t,".only"))[grep(morph, rownames(get(paste0("NMDS.",metric,".P.",t,".only")))),])
  #     assign(paste0("NMDS.",metric,".P.",t,".only.",morph,".chull"), chull(get(paste0("NMDS.",metric,".P.",t,".only.",morph))))
  #     assign(paste0("NMDS.",metric,".P.",t,".only.",morph,".chull"), c(get(paste0("NMDS.",metric,".P.",t,".only.",morph,".chull")), get(paste0("NMDS.",metric,".P.",t,".only.",morph,".chull"))[1]))
  #   }
  # }
  
  # TIME-- ONE at a time
  for (t in TP) {
    # Per time point
    
    assign(paste0("NMDS.",metric,".P.",t,".only"), get(paste0("NMDS.",metric,".P.morphonly"))$points[grep(paste0("^",t,"."), rownames(get(paste0("NMDS.",metric,".P.morphonly"))$points)),])
    for (morph in c("FB","BL","CR")) {
      assign(paste0("NMDS.",metric,".P.",t,".only.",morph), get(paste0("NMDS.",metric,".P.",t,".only"))[grep(morph, rownames(get(paste0("NMDS.",metric,".P.",t,".only")))),])
      assign(paste0("NMDS.",metric,".P.",t,".only.",morph,".chull"), chull(get(paste0("NMDS.",metric,".P.",t,".only.",morph))))
      assign(paste0("NMDS.",metric,".P.",t,".only.",morph,".chull"), c(get(paste0("NMDS.",metric,".P.",t,".only.",morph,".chull")), get(paste0("NMDS.",metric,".P.",t,".only.",morph,".chull"))[1]))
    }
  }

 
  ############ COMBO DISP BETA ################
  xvalues <- as.character(get(paste0("FB.BL.",metric,".P.ALL.mean"))[,1])
  
  
  # Change factor of individual plots
  MF.P.20.only$Morph <- factor(MF.P.20.only$Morph, levels = c("CR","BL","FB"))
  MF.P.60.only$Morph <- factor(MF.P.60.only$Morph, levels = c("CR","BL","FB"))
  MF.P.180.only$Morph <- factor(MF.P.180.only$Morph, levels = c("CR","BL","FB"))
  MF.P.360.only$Morph <- factor(MF.P.360.only$Morph, levels = c("CR","BL","FB"))
  MF.P.720.only$Morph <- factor(MF.P.720.only$Morph, levels = c("CR","BL","FB"))
  MF.P.1440.only$Morph <- factor(MF.P.1440.only$Morph, levels = c("CR","BL","FB"))
  
  # Get significance for each timepoint
  listSig <- c()
  for (t in TP) {
    tempSig <- get(paste0("ANOVA.",metric,".P.",t,".only"))$aov.tab$`Pr(>F)`[1]
    if (as.numeric(tempSig) <= 0.001) {
      listSig <- c(listSig, "***")
    } else if (as.numeric(tempSig) <= 0.01) {
      listSig <- c(listSig, "**")
    } else if (as.numeric(tempSig) <= 0.05) {
      listSig <- c(listSig,"*")
    } else {
      listSig <- c(listSig, " ")
    }
  }
  names(listSig) <- TP
  
  
  pdf(paste0("BETAPLOTS_P/COMBO_P_dispbeta_",metric,".pdf"), pointsize = 14, width = 10, height = 6)
  par(fig = c(0,0.8,0.3,1))
  plot(xvalues, NULL
       , main = "Dispersion of morphologies across time"
       , xlab = ''
       , ylab = paste0("Distance (",fullMetrics[which(metrics %in% metric)],")")
       , ylim = ylimits
       , xaxt = 'n')
  axis(side = 1
       , at = c(1,2,3,4,5,6)
       , labels = c("20 min","1 h","3 h","6 h","12 h","24 h")
       , las = 2)
  points(get(paste0("FB.BL.",metric,".P.ALL.mean"))[,2] ~ c(1,2,3,4,5,6)
         , type = 'l'
         , col = "purple"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5,6)*0.99
         , x1 = c(1,2,3,4,5,6)*0.99
         , y0 = c(get(paste0("FB.BL.",metric,".P.ALL.mean"))[,2] - get(paste0("FB.BL.",metric,".P.ALL.sd"))[,2])
         , y1 = c(get(paste0("FB.BL.",metric,".P.ALL.mean"))[,2] + get(paste0("FB.BL.",metric,".P.ALL.sd"))[,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  points(get(paste0("FB.CR.",metric,".P.ALL.mean"))[,2] ~ c(1,2,3,4,5,6)
         , type = 'l'
         , col = "darkgreen"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5,6)*1
         , x1 = c(1,2,3,4,5,6)*1
         , y0 = c(get(paste0("FB.CR.",metric,".P.ALL.mean"))[,2] - get(paste0("FB.CR.",metric,".P.ALL.sd"))[,2])
         , y1 = c(get(paste0("FB.CR.",metric,".P.ALL.mean"))[,2] + get(paste0("FB.CR.",metric,".P.ALL.sd"))[,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  points(get(paste0("CR.BL.",metric,".P.ALL.mean"))[,2] ~ c(1,2,3,4,5,6)
         , type = 'l'
         , col = "grey"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4)*1.01
         , x1 = c(1,2,3,4)*1.01
         , y0 = c(get(paste0("CR.BL.",metric,".P.ALL.mean"))[,2] - get(paste0("CR.BL.",metric,".P.ALL.sd"))[,2])
         , y1 = c(get(paste0("CR.BL.",metric,".P.ALL.mean"))[,2] + get(paste0("CR.BL.",metric,".P.ALL.sd"))[,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  title(xlab = paste0("Stress: ", round(get(paste0("NMDS.",metric,".P.morphonly"))$stress/100,2))
        , line = 3)
  par(fig = c(0.7,1,0.2,1), mar = c(2,0,0,0),oma = c(0,0,0,0), new = TRUE)
  plot(0,0
       , pch = ""
       , xlab = ""
       , ylab = ""
       , xaxt = "n"
       , yaxt = "n"
       , bty = "n")
  legend("center"
         , legend = c("FB:CR", "FB:BL", "CR:BL")
         , lty = 1
         , col = c("darkgreen","purple","grey")
         , lwd = 2 
  )
  par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.05,0.158333,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".P.20.only"))#$points
       # , main = "20 Minutes"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.P.20.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
 )
  title(xlab = paste0(listSig["20"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ", round(get(paste0("NMDS.",metric,".P.20.only"))$stress,2)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)
  # title(xlab = paste0("Df = ", ANOVA.UWUF.P.20.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.20.only$aov.tab[1]$Df[3])
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.20.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.P.20.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".P.20.only.CR"))[get(paste0("NMDS.",metric,".P.20.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".P.20.only.BL"))[get(paste0("NMDS.",metric,".P.20.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".P.20.only.FB"))[get(paste0("NMDS.",metric,".P.20.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.178333,0.28666,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".P.60.only"))#$points
       # , main = "1 Hour"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.P.60.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
  )
  title(xlab = paste0(listSig["60"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ", round(get(paste0("NMDS.",metric,".P.60.only"))$stress,2)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)
  # title(xlab = paste0("Df = ", ANOVA.UWUF.P.60.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.60.only$aov.tab[1]$Df[3])
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.60.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.P.60.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".P.60.only.CR"))[get(paste0("NMDS.",metric,".P.60.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".P.60.only.BL"))[get(paste0("NMDS.",metric,".P.60.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".P.60.only.FB"))[get(paste0("NMDS.",metric,".P.60.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.30666,0.414999,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".P.180.only"))#$points
       # , main = "6 Hours"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.P.180.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
  )
  title(xlab = paste0(listSig["180"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ", round(get(paste0("NMDS.",metric,".P.180.only"))$stress,2)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)
  # title(xlab = paste0("Df = ", ANOVA.UWUF.P.180.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.180.only$aov.tab[1]$Df[3])
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.180.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.P.180.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".P.180.only.CR"))[get(paste0("NMDS.",metric,".P.180.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".P.180.only.BL"))[get(paste0("NMDS.",metric,".P.180.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".P.180.only.FB"))[get(paste0("NMDS.",metric,".P.180.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.434999,0.543333,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".P.360.only"))#$points
       # , main = "6 Hours"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.P.360.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
  )
  title(xlab = paste0(listSig["360"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ", round(get(paste0("NMDS.",metric,".P.360.only"))$stress,2)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)
  # title(xlab = paste0("Df = ", ANOVA.UWUF.P.360.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.360.only$aov.tab[1]$Df[3])
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.360.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.P.360.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".P.360.only.CR"))[get(paste0("NMDS.",metric,".P.360.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".P.360.only.BL"))[get(paste0("NMDS.",metric,".P.360.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".P.360.only.FB"))[get(paste0("NMDS.",metric,".P.360.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.565555,0.671666,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".P.720.only"))#$points
       # , main = "12 Hours"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.P.720.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
  )
  title(xlab = paste0(listSig["720"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ", round(get(paste0("NMDS.",metric,".P.720.only"))$stress,2)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)
  # title(xlab = paste0("Df = ", ANOVA.UWUF.P.720.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.720.only$aov.tab[1]$Df[3])
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.720.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.P.720.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".P.720.only.CR"))[get(paste0("NMDS.",metric,".P.720.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".P.720.only.BL"))[get(paste0("NMDS.",metric,".P.720.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".P.720.only.FB"))[get(paste0("NMDS.",metric,".P.720.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.691666,0.79999,0,0.4), new = TRUE)
  plot(get(paste0("NMDS.",metric,".P.1440.only"))#$points
       # , main = "12 Hours"
       , pch = 21
       , col = "black"
       , bg = MorphColours[factor(MF.P.1440.only$Morph)]
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = 'n'
       , ylab = ''
  )
  title(xlab = paste0(listSig["1440"])
        , line = 1.25
        , cex.lab = 3)
  # title(xlab = paste0("Stress: ", round(get(paste0("NMDS.",metric,".P.1440.only"))$stress,2)/100)
  #       , line = 1.5
  #       , cex.lab = 0.75)
  # title(xlab = paste0("Df = ", ANOVA.UWUF.P.1440.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.1440.only$aov.tab[1]$Df[3])
  #       , line = 0.25
  #       , cex.axis = 0.5)
  # title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.1440.only$aov.tab[5]$R2[1], digits = 2))) 
  #       , line= 1
  #       , cex.axis = 0.5)
  # title(xlab = paste0("p = ", ANOVA.UWUF.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
  #       , line = 1.75
  #       , cex.axis = 0.5)
  lines(get(paste0("NMDS.",metric,".P.1440.only.CR"))[get(paste0("NMDS.",metric,".P.1440.only.CR.chull")),]
        , col = MorphColours[1])
  lines(get(paste0("NMDS.",metric,".P.1440.only.BL"))[get(paste0("NMDS.",metric,".P.1440.only.BL.chull")),]
        , col = MorphColours[2])
  lines(get(paste0("NMDS.",metric,".P.1440.only.FB"))[get(paste0("NMDS.",metric,".P.1440.only.FB.chull")),]
        , col = MorphColours[3])
  par(oma = c(1,0.1,1,0.1), mar = c(2,0,2,0), fig = c(0.775,1,0,0.4), new = TRUE)
  plot(0,0
       , pch = ''
       , bty = 'n'
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = ''
       , ylab = ''
  )
  legend("center"
         , legend = levels(MF.P.morphkeep$Morph)
         , pch = 21
         , col = "black"
         , pt.bg = MorphColours
  )
  #STOP
  dev.off()
  
  # Combine individual stats and print out
  
  allindividualTests <- c(paste0("ANOVA.",metric,".P.20.only")
                          ,paste0("ANOVA.",metric,".P.60.only")
                          ,paste0("ANOVA.",metric,".P.180.only")
                          ,paste0("ANOVA.",metric,".P.360.only")
                          ,paste0("ANOVA.",metric,".P.720.only")
                          ,paste0("ANOVA.",metric,".P.1440.only")
                          
  )
  for (i in allindividualTests) {
    # print(get(i))
    capture.output(get(i), file = paste0("./BETAPLOTS_P/individualtests/",i,".txt"))
  }
  
  
  ########### BETADISP#############

  assign(paste0("betadisp.P.",metric,".FB"), get(paste0("betadisp.",metric,".P.time"))$distances[grep("FB", names(get(paste0("betadisp.",metric,".P.time"))$distances))])
  assign(paste0("betadisp.P.",metric,".BL"), get(paste0("betadisp.",metric,".P.time"))$distances[grep("BL", names(get(paste0("betadisp.",metric,".P.time"))$distances))])
  assign(paste0("betadisp.P.",metric,".CR"), get(paste0("betadisp.",metric,".P.time"))$distances[grep("CR", names(get(paste0("betadisp.",metric,".P.time"))$distances))])
  
  MF.P.FB <- MF.P.morphkeep[grep("FB", rownames(MF.P.morphkeep)),]
  MF.P.BL <- MF.P.morphkeep[grep("BL", rownames(MF.P.morphkeep)),]
  MF.P.CR <- MF.P.morphkeep[grep("CR", rownames(MF.P.morphkeep)),]
  
  betadisp.FB.P.agg <- aggregate(get(paste0("betadisp.P.",metric,".FB")), by = list(MF.P.FB$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
  betadisp.BL.P.agg <- aggregate(get(paste0("betadisp.P.",metric,".BL")), by = list(MF.P.BL$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
  betadisp.CR.P.agg <- aggregate(get(paste0("betadisp.P.",metric,".CR")), by = list(MF.P.CR$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
  
  assign(paste0("betadisp.",metric,".P.time.forstat"), cbind(Distance = get(paste0("betadisp.",metric,".P.time"))$distances
                                                             , MF.P.morphkeep[unlist(lapply(names(get(paste0("betadisp.",metric,".P.time"))$distances)
                                                                                          , function(x) {grep(paste0("^",x,"$"), rownames(MF.P.morphkeep))})),]))
  

  assign(paste0("ANOVA.betadisp.",metric,".P"), anova(lm(Distance ~ Time*Morph, data = get(paste0("betadisp.",metric,".P.time.forstat")))))
  capture.output(get(paste0("ANOVA.betadisp.",metric,".P")), file = paste0("./BETAPLOTS_P/ANOVA.betadisp.",metric,".txt"))
  
  
  maxy <- ceiling(max(unlist(lapply(c("FB","BL","CR"), function(x) {
    naList <- which(is.na(get(paste0("betadisp.",x,".P.agg"))[,2][,2]))
    sdTemp <- get(paste0("betadisp.",x,".P.agg"))[,2][,2]
    if (length(naList) >=1) {
      sdTemp[naList] <- 0
    }
    max(get(paste0("betadisp.",x,".P.agg"))[,2][,1] + sdTemp, na.rm = TRUE)
  })))*100)/100
  miny <- floor(min(unlist(lapply(c("FB","BL","CR"), function(x) {
    naList <- which(is.na(get(paste0("betadisp.",x,".P.agg"))[,2][,2]))
    sdTemp <- get(paste0("betadisp.",x,".P.agg"))[,2][,2]
    if (length(naList) >=1) {
      sdTemp[naList] <- 0
    }
    min(get(paste0("betadisp.",x,".P.agg"))[,2][,1] - sdTemp, na.rm = TRUE)
  })))*100)/100
  ylimits <- c(miny, maxy)
  
  pdf(paste0("./BETAPLOTS_P/BetaDisp_P_",metric,"_eachmorph.pdf"),pointsize = 14)
  plot(xvalues, NULL
       , main = "Dispersion of morphologies across time"
       , xlab = 'Time'
       , ylab = paste0("Distance (",fullMetrics[which(metrics %in% metric)],")")
       , ylim = ylimits
       , xaxt = 'n')
  axis(side = 1
       , at = c(1,2,3,4,5,6)
       , labels = c("20 min","1 h","3 h","6 h","12 h", "1 d")
       , las = 2)
  points(betadisp.FB.P.agg[,2][,1] ~ c(1,2,3,4,5,6)
         , type = 'l'
         , col = "salmon"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5,6)*0.99
         , x1 = c(1,2,3,4,5,6)*0.99
         , y0 = c(betadisp.FB.P.agg[,2][,1] - betadisp.FB.P.agg[,2][,2])
         , y1 = c(betadisp.FB.P.agg[,2][,1] + betadisp.FB.P.agg[,2][,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  points(betadisp.BL.P.agg[,2][,1] ~ c(1,2,3,4,5,6)
         , type = 'l'
         , col = "dodgerblue"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5,6)*1
         , x1 = c(1,2,3,4,5,6)*1
         , y0 = c(betadisp.BL.P.agg[,2][,1] - betadisp.BL.P.agg[,2][,2])
         , y1 = c(betadisp.BL.P.agg[,2][,1] + betadisp.BL.P.agg[,2][,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  points(betadisp.CR.P.agg[,2][,1] ~ c(1,2,3,4,5,6)
         , type = 'l'
         , col = "darkorchid4"
         , lwd = 2
         , lty = 1)
  arrows(x0 = c(1,2,3,4,5,6)*1.01
         , x1 = c(1,2,3,4,5,6)*1.01
         , y0 = c(betadisp.CR.P.agg[,2][,1] - betadisp.CR.P.agg[,2][,2])
         , y1 = c(betadisp.CR.P.agg[,2][,1] + betadisp.CR.P.agg[,2][,2])
         , angle = 90
         , code = 3
         , length = 0.03)
  par(fig = c(0.7,1,0.5,1), mar = c(2,0,0,0),oma = c(0,0,0,0), new = TRUE)
  plot(0,0
       , pch = ""
       , xlab = ""
       , ylab = ""
       , xaxt = "n"
       , yaxt = "n"
       , bty = "n")
  legend("center"
         , legend = c("FB", "BL", "CR")
         , lty = 1
         , col = c("salmon","dodgerblue","darkorchid4")
         , lwd = 2 
  )
  dev.off()
  
  ################## Combine individual stats and print out ####################
  ### DO BETADISP FOR OVERALL ##
  assign(paste0("anova.betadisp.",metric,".morph"), anova(get(paste0("betadisp.",metric,".morph"))))
  assign(paste0("anova.betadisp.",metric,".time"), anova(get(paste0("betadisp.",metric,".time"))))
  assign(paste0("anova.betadisp.P.",metric,".morph"), anova(get(paste0("betadisp.",metric,".P.morph"))))
  assign(paste0("anova.betadisp.P.",metric,".time"), anova(get(paste0("betadisp.",metric,".P.time"))))
  
  # This is PERMADISP-- don't need a table, will quote in text
  capture.output(get(paste0("anova.betadisp.",metric,".morph")), file = paste0("BETAPLOTS_H/anova.betadisp.",metric,".morph.txt"))
  capture.output(get(paste0("anova.betadisp.",metric,".time")), file = paste0("BETAPLOTS_H/anova.betadisp.",metric,".time.txt"))
  capture.output(get(paste0("anova.betadisp.P.",metric,".morph")), file = paste0("BETAPLOTS_P/anova.betadisp.",metric,".P.morph.txt"))
  capture.output(get(paste0("anova.betadisp.P.",metric,".time")), file = paste0("BETAPLOTS_P/anova.betadisp.",metric,".P.time.txt"))
  
  ### NOW DO EACH TIME POINT ##
  ## This is a table with both Hakai and PM for the supplementary figures
  permdisp.morphology.across.time <- matrix(ncol = 8, nrow = 2)
  rownames(permdisp.morphology.across.time) <- c("P","H")
  colnames(permdisp.morphology.across.time) <- c("20","60","180","360","720","1440","5760", "Overall")
  for (t in c("20","60","180","360","720","1440")) {
    assign(paste0("betadisp.P.",t), betadisper(as.dist(get(paste0("dm.",metric,".P.",t))), group = get(paste0("MF.P.",t,".only"))$Morph))
    assign(paste0("anova.betadisp.P.",t), anova(get(paste0("betadisp.P.",t))))
    
    if (get(paste0("anova.betadisp.P.",t))$`Pr(>F)`[1] < 0.001) {
      ptemp <- paste0("\\bf{<0.001}")
    } else if (get(paste0("anova.betadisp.P.",t))$`Pr(>F)`[1] <= 0.05) {
      ptemp <- paste0("\\bf{", signif(get(paste0("anova.betadisp.P.",t))$`Pr(>F)`[1],2),"}")
    } else {
      ptemp <- signif(get(paste0("anova.betadisp.P.",t))$`Pr(>F)`[1],2)
    }
    ftemp <- get(paste0("anova.betadisp.P.",t))$`F value`[1]
    dftemp <- paste0(get(paste0("anova.betadisp.P.",t))$Df[1],",",get(paste0("anova.betadisp.P.",t))$Df[2])
    
    toPaste <- paste0(ptemp, " ($F_{",dftemp,"}$=",round(ftemp,2),")")
    permdisp.morphology.across.time["P", paste0(t)] <- toPaste
  }
  for (t in c("20","60","360","720","5760")) {
    assign(paste0("betadisp.",t), betadisper(as.dist(get(paste0("dm.",metric,".",t))), group = get(paste0("MF.",t,".only"))$Morph))
    assign(paste0("anova.betadisp.",t), anova(get(paste0("betadisp.",t))))
    
    if (get(paste0("anova.betadisp.",t))$`Pr(>F)`[1] < 0.001) {
      ptemp <- paste0("\\bf{<0.001}")
    } else if (get(paste0("anova.betadisp.",t))$`Pr(>F)`[1] <= 0.05) {
      ptemp <- paste0("\\bf{", signif(get(paste0("anova.betadisp.",t))$`Pr(>F)`[1],2),"}")
    } else {
      ptemp <- signif(get(paste0("anova.betadisp.",t))$`Pr(>F)`[1],2)
    }
    
    ftemp <- get(paste0("anova.betadisp.",t))$`F value`[1]
    dftemp <- paste0(get(paste0("anova.betadisp.",t))$Df[1],",",get(paste0("anova.betadisp.",t))$Df[2])
    
    toPaste <- paste0(ptemp, " ($F_{",dftemp,"}$=",round(ftemp,2),")")
    permdisp.morphology.across.time["H", paste0(t)] <- toPaste
    
  }
  # Get overall
  if (get(paste0("anova.betadisp.P.",metric,".morph"))$`Pr(>F)`[1] < 0.001) {
    ptemp <- paste0("\\bf{<0.001}")
  } else if (get(paste0("anova.betadisp.P.",metric,".morph"))$`Pr(>F)`[1] <= 0.05) {
    ptemp <- paste0("\\bf{", signif(get(paste0("anova.betadisp.P.",metric,".morph"))$`Pr(>F)`[1],2),"}")
  } else {
    ptemp <- signif(get(paste0("anova.betadisp.P.",metric,".morph"))$`Pr(>F)`[1],2)
  }
  ftemp <- get(paste0("anova.betadisp.P.",metric,".morph"))$`F value`[1]
  dftemp <- paste0(get(paste0("anova.betadisp.P.",metric,".morph"))$Df[1],",",get(paste0("anova.betadisp.P.",metric,".morph"))$Df[2])
  toPaste <- paste0(ptemp," ($F_{",dftemp,"}$=", round(ftemp,2),")")
  permdisp.morphology.across.time["P","Overall"] <- toPaste
  
  if (get(paste0("anova.betadisp.",metric,".morph"))$`Pr(>F)`[1] < 0.001) {
    ptemp <- paste0("\\bf{<0.001}")
    
  } else if (get(paste0("anova.betadisp.",metric,".morph"))$`Pr(>F)`[1] <= 0.05) {
    ptemp <- paste0("\\bf{", signif(get(paste0("anova.betadisp.",metric,".morph"))$`Pr(>F)`[1],2),"}")
  } else {
    ptemp <- signif(get(paste0("anova.betadisp.",metric,".morph"))$`Pr(>F)`[1],2)
  }
  
  ftemp <- get(paste0("anova.betadisp.",metric,".morph"))$`F value`[1]
  dftemp <- paste0(get(paste0("anova.betadisp.",metric,".morph"))$Df[1],",",get(paste0("anova.betadisp.",metric,".morph"))$Df[2])
  toPaste <- paste0(ptemp, " ($F_{",dftemp,"}$=",round(ftemp,2),")")
  permdisp.morphology.across.time["H","Overall"] <- toPaste
  
  colnames(permdisp.morphology.across.time) <- c("20 minutes","1 hour","3 hours","6 hours","12 hours","1 day","4 days", "Overall")
  
  tempfile1 <- permdisp.morphology.across.time
  for (r in 1:nrow(tempfile1)) {
    for (c in 1:ncol(tempfile1)) {
      if (is.na(tempfile1[r,c])) {
        permdisp.morphology.across.time[r,c] <- "-"
      }}}
  
  
  # Make double header
  permdisp.morphology.across.time <- cbind(Site = c("Reed Point","Hakai"), permdisp.morphology.across.time)
  assign(paste0("permdisp.morphology.across.time.",metric), permdisp.morphology.across.time)

  capture.output(print(xtable(get(paste0("permdisp.morphology.across.time.",metric)))
                       , include.rownames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x}
                       )
                 , file = paste0("BETAPLOTS_LATEX/permdisp.morph.across.time.",metric,".txt"))
  
  ### MAKE BETA DIV TABLES-- metrics separately but H and P together; extras will go in supp
  
  anova.morphology.across.time <- matrix(ncol = 8, nrow = 2)
  colnames(anova.morphology.across.time) <- c("20","60","180","360","720","1440","5760", "Overall")
  rownames(anova.morphology.across.time) <- c("P","H")
  for (t in c("20","60","180","360","720","1440")) {
    
    if (get(paste0("ANOVA.",metric,".P.",t,".only"))$aov.tab$`Pr(>F)`[1] < 0.001) {
      ptemp <- paste0("\\bf{<0.001}")
    } else if (get(paste0("ANOVA.",metric,".P.",t,".only"))$aov.tab$`Pr(>F)`[1] <= 0.05) {
      ptemp <- paste0("\\bf{", signif(get(paste0("ANOVA.",metric,".P.",t,".only"))$aov.tab$`Pr(>F)`[1],2),"}")
    } else {
      ptemp <- signif(get(paste0("ANOVA.",metric,".P.",t,".only"))$aov.tab$`Pr(>F)`[1],2)
    }
    rtemp <- get(paste0("ANOVA.",metric,".P.",t,".only"))$aov.tab$`R2`[1]
    dftemp <- paste0(get(paste0("ANOVA.",metric,".P.",t,".only"))$aov.tab$`Df`[1],",",get(paste0("ANOVA.",metric,".P.",t,".only"))$aov.tab$`Df`[3])
    toPaste <- paste0(ptemp, " ($R^2$=",round(rtemp,2),",$df$=",dftemp,")")
    
    anova.morphology.across.time["P",paste0(t)] <- toPaste
  }
  for (t in c("20","60","360","720","5760")) {
    if (get(paste0("ANOVA.",metric,".",t,".only"))$aov.tab$`Pr(>F)`[1] < 0.001) {
      ptemp <- paste0("\\bf{<0.001}")
    } else if (get(paste0("ANOVA.",metric,".",t,".only"))$aov.tab$`Pr(>F)`[1] <= 0.05) {
      ptemp <- paste0("\\bf{", signif(get(paste0("ANOVA.",metric,".",t,".only"))$aov.tab$`Pr(>F)`[1],2),"}")
    } else {
      ptemp <- signif(get(paste0("ANOVA.",metric,".",t,".only"))$aov.tab$`Pr(>F)`[1],2)
    }
    rtemp <- get(paste0("ANOVA.",metric,".",t,".only"))$aov.tab$`R2`[1]
    dftemp <- paste0(get(paste0("ANOVA.",metric,".",t,".only"))$aov.tab$`Df`[1],",",get(paste0("ANOVA.",metric,".",t,".only"))$aov.tab$`Df`[3])
    toPaste <- paste0(ptemp, " ($R^2$=",round(rtemp,2),",$df$=",dftemp,")")
    
    anova.morphology.across.time["H",paste0(t)] <- toPaste
  }
  # Do overall P
  if (get(paste0("ANOVA.",metric,".P.morphtime"))$aov.tab$`Pr(>F)`[2] < 0.001) {
    ptemp <- paste0("\\bf{<0.001}")
  } else if (get(paste0("ANOVA.",metric,".P.morphtime"))$aov.tab$`Pr(>F)`[2] <= 0.05) {
    ptemp <- paste0("\\bf{", signif(get(paste0("ANOVA.",metric,".P.morphtime"))$aov.tab$`Pr(>F)`[2],2),"}")
  } else {
    ptemp <- signif(get(paste0("ANOVA.",metric,".P.morphtime"))$aov.tab$`Pr(>F)`[2],2)
  }
  rtemp <- get(paste0("ANOVA.",metric,".P.morphtime"))$aov.tab$`R2`[2]
  dftemp <- paste0(get(paste0("ANOVA.",metric,".P.morphtime"))$aov.tab$`Df`[2],",",get(paste0("ANOVA.",metric,".P.morphtime"))$aov.tab$`Df`[4])
  toPaste <- paste0(ptemp, " ($R^2$=",round(rtemp,2),",$df$=",dftemp,")")
  anova.morphology.across.time["P","Overall"] <- toPaste
  
  # Do overall H
  if (get(paste0("ANOVA.",metric,".morphtime"))$aov.tab$`Pr(>F)`[2] < 0.001) {
    ptemp <- paste0("\\bf{<0.001}")
  } else if (get(paste0("ANOVA.",metric,".morphtime"))$aov.tab$`Pr(>F)`[2] <= 0.05) {
    ptemp <- paste0("\\bf{", signif(get(paste0("ANOVA.",metric,".morphtime"))$aov.tab$`Pr(>F)`[2],2),"}")
  } else {
    ptemp <- signif(get(paste0("ANOVA.",metric,".morphtime"))$aov.tab$`Pr(>F)`[2],2)
  }
  rtemp <- get(paste0("ANOVA.",metric,".morphtime"))$aov.tab$`R2`[2]
  dftemp <- paste0(get(paste0("ANOVA.",metric,".morphtime"))$aov.tab$`Df`[2],",",get(paste0("ANOVA.",metric,".morphtime"))$aov.tab$`Df`[4])
  toPaste <- paste0(ptemp, " ($R^2$=",round(rtemp,2),",$df$=",dftemp,")")
  anova.morphology.across.time["H","Overall"] <- toPaste
  
  colnames(anova.morphology.across.time) <- c("20 minutes","1 hour","3 hours","6 hours","12 hours","1 day","4 days", "Overall")
  
  tempfile1 <- anova.morphology.across.time
  for (r in 1:nrow(tempfile1)) {
    for (c in 1:ncol(tempfile1)) {
      if (is.na(tempfile1[r,c])) {
        anova.morphology.across.time[r,c] <- "-"
      } 
    }
  }
  
  # Make double header
  anova.morphology.across.time <- cbind(Site = c("Reed Point","Hakai"),Factors = c("Morph","Morph"), anova.morphology.across.time)
  assign(paste0("anova.morphology.across.time.",metric), anova.morphology.across.time)

  capture.output(print(xtable(get(paste0("anova.morphology.across.time.",metric)))
                       , include.rownames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x}
                       )
                 , file = paste0("BETAPLOTS_LATEX/anova.morph.across.time.",metric,".txt"))
  
  ### DO FB/CR/BL TEST FOR EACH
  # Make table
  pairwiseAdonis.all <- matrix(ncol = 8,nrow = 6)
  colnames(pairwiseAdonis.all) <- c("20","60","180","360","720","1440", "5760","Overall")
  rownames(pairwiseAdonis.all) <- c("pFB:BL","pFB:CR","pBL:CR","hFB:BL","hFB:CR","hBL:CR" )
  listMorphs <- c("FB","BL","CR")
  for (m in 1:(length(listMorphs)-1)) {
    for (n in (m+1):length(listMorphs)) {
      
      #  Reed point
      tempMFALL.P <- MF.P.inclWater[grep(paste0(listMorphs[m],"|",listMorphs[n]), MF.P.inclWater$Morph),]
      tempDMALL.P <- get(paste0("dm.",metric,".P.inclWater"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.",metric,".P.inclWater")))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.",metric,".P.inclWater"))))]
      tempAdonisALL.P <- adonis(tempDMALL.P ~ Morph, data = tempMFALL.P, by = "marginal")
      ptemp <- tempAdonisALL.P$aov.tab$`Pr(>F)`[1]
      ptemp <- p.adjust(ptemp, method = "fdr",n = 3)
      
      if (ptemp < 0.001) {
        ptemp <- paste0("\\bf{<0.001}")
      } else if (ptemp <= 0.05) {
        ptemp <- paste0("\\bf{",signif(ptemp,2),"}")
      } else {
        ptemp <- signif(ptemp,2)
      }
      toPaste <- paste0( ptemp
                         ," ($R^2$=", round(tempAdonisALL.P$aov.tab$R2[1],digits = 2)
                         ,", $df$=", paste0(tempAdonisALL.P$aov.tab$Df[1],",",tempAdonisALL.P$aov.tab$Df[3])
                         , ")")
      pairwiseAdonis.all[paste0("p",listMorphs[m],":",listMorphs[n]),"Overall"] <- toPaste
      
      # Hakai
      tempMFALL <- MF.inclWater[grep(paste0(listMorphs[m],"|",listMorphs[n]), MF.inclWater$Morph),]
      tempDMALL <- get(paste0("dm.",metric,".inclWater"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.",metric,".inclWater")))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.",metric,".inclWater"))))]
      tempAdonisALL <- adonis(tempDMALL ~ Morph, data = tempMFALL, by = "marginal")
      ptemp <- tempAdonisALL$aov.tab$`Pr(>F)`[1]
      ptemp <- p.adjust(ptemp, method = "fdr", n = 3)
      
      if (ptemp < 0.001) {
        ptemp <- paste0("\\bf{<0.001}")
      } else if (ptemp <= 0.05) {
        ptemp <- paste0("\\bf{",signif(ptemp,2),"}")
      } else {
        ptemp <- signif(ptemp,2)
      }
      toPaste <- paste0( ptemp
                         ," ($R^2$=", round(tempAdonisALL$aov.tab$R2[1],digits = 2)
                         ,", $df$=", paste0(tempAdonisALL$aov.tab$Df[1],",",tempAdonisALL$aov.tab$Df[3])
                         , ")")
      pairwiseAdonis.all[paste0("h",listMorphs[m],":",listMorphs[n]),"Overall"] <- toPaste
      
      for (t in c("20","60","180","360","720","1440")) {
        tempMF <- get(paste0("MF.P.",t,".only"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), get(paste0("MF.P.",t,".only"))$Morph),]
        tempDM <- get(paste0("dm.",metric,".P.",t))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.",metric,".P.",t)))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.",metric,".P.",t))))]
        tempAdonis <- adonis(tempDM ~ Morph, data = tempMF, by = "marginal")
        ptemp <- tempAdonis$aov.tab$`Pr(>F)`[1]
        ptemp <- p.adjust(ptemp, method = "fdr", n = 3)
        
        if (ptemp < 0.001) {
          ptemp <- paste0("\\bf{<0.001}")
        } else if (ptemp <= 0.05) {
          ptemp <- paste0("\\bf{",signif(ptemp,2),"}")
        } else {
          ptemp <- signif(ptemp,2)
        }
        toPaste <- paste0( ptemp
                           ," ($R^2$=", round(tempAdonis$aov.tab$R2[1],digits = 2)
                           ,", $df$=", paste0(tempAdonis$aov.tab$Df[1],",",tempAdonis$aov.tab$Df[3])
                           , ")")
        pairwiseAdonis.all[paste0("p",listMorphs[m],":",listMorphs[n]),paste0(t)] <- toPaste
      }
      for (t in c("20","60","360","720","5760")) {
        tempMF <- get(paste0("MF.",t,".only"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), get(paste0("MF.",t,".only"))$Morph),]
        tempDM <- get(paste0("dm.",metric,".",t))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.",metric,".",t)))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.",metric,".",t))))]
        tempAdonis <- adonis(tempDM ~ Morph, data = tempMF, by = "marginal")
        ptemp <- tempAdonis$aov.tab$`Pr(>F)`[1]
        ptemp <- p.adjust(ptemp, method = "fdr", n = 3)
        
        if (ptemp < 0.001) {
          ptemp <- paste0("\\bf{<0.001}")
        } else if (ptemp <= 0.05) {
          ptemp <- paste0("\\bf{",signif(ptemp,2),"}")
        } else {
          ptemp <- signif(ptemp,2)
        }
        toPaste <- paste0( ptemp
                           ," ($R^2$=", round(tempAdonis$aov.tab$R2[1],digits = 2)
                           ,", $df$=", paste0(tempAdonis$aov.tab$Df[1],",",tempAdonis$aov.tab$Df[3])
                           , ")")
        pairwiseAdonis.all[paste0("h",listMorphs[m],":",listMorphs[n]),paste0(t)] <- toPaste
      }
    }
  }
  
  # Overall
  
  
  
  # Get rid of NAs
  for (r in 1:nrow(pairwiseAdonis.all)) {
    for (c in 1:ncol(pairwiseAdonis.all)) {
      if (is.na(pairwiseAdonis.all[r,c])) {
        pairwiseAdonis.all[r,c] <- "-"
      }
    }
  }
  
  # Change Rownames
  assign(paste0("pairwiseAdonis.all.",metric), cbind(Sites = c("Reed Point", " ","  ","Hakai","   ","    "),Factors = c("FB:BL","FB:CR","BL:CR","FB:BL","FB:CR","BL:CR"), pairwiseAdonis.all))

  capture.output(print(xtable(get(paste0("pairwiseAdonis.all.",metric)))
                       , include.rownames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x}
                       )
                 , file = paste0("BETAPLOTS_LATEX/pairwiseAdonis.",metric,".txt"))
  
  
}

############ PRINTING LATEX-DEP STATS ############
# COMBINING BETA DISP
permdisp.morphology.across.time.MASTER <- matrix(ncol = 10, nrow = 6)
colnames(permdisp.morphology.across.time.MASTER) <- colnames(permdisp.morphology.across.time.UWUF)
rownames(permdisp.morphology.across.time.MASTER) <- c("pBC","pWUF","pUWUF", "hBC","hWUF","hUWUF")
for (met in c("BC","WUF","UWUF")) {
  PMtemp <- get(paste0("permdisp.morphology.across.time.",met))[1,]
  HKtemp <- get(paste0("permdisp.morphology.across.time.",met))[2,]
  
  permdisp.morphology.across.time.MASTER[paste0("p",met),] <- PMtemp
  permdisp.morphology.across.time.MASTER[paste0("h",met),] <- HKtemp
  
}
# Change rownames

permdisp.morphology.across.time.MASTER.edit <- cbind(rep(c("Bray-Curtis","Weighted Unifrac","Un-weighted Unifrac"),2), permdisp.morphology.across.time.MASTER)
rownames(permdisp.morphology.across.time.MASTER.edit) <- c("Reed Point"," ","  ","Hakai","   ","    ")

capture.output(print(xtable(permdisp.morphology.across.time.MASTER.edit)
                     , math.style.negative = TRUE
                     , math.style.exponents = TRUE
                     , sanitize.text.function = function(x) {x}
                     )
               , file = paste0("BETAPLOTS_LATEX/permdisp.morphology.across.time.MASTER.txt"))





# COMBINING BETA ANOVAs
anova.morphology.across.time.MASTER <- matrix(ncol = 9, nrow = 6)
colnames(anova.morphology.across.time.MASTER) <- colnames(anova.morphology.across.time.UWUF)
rownames(anova.morphology.across.time.MASTER) <- c("pBC","pWUF","pUWUF", "hBC","hWUF","hUWUF")
for (met in c("BC","WUF","UWUF")) {
  PMtemp <- get(paste0("anova.morphology.across.time.",met))[1,]
  HKtemp <- get(paste0("anova.morphology.across.time.",met))[2,]
  
  anova.morphology.across.time.MASTER[paste0("p",met),] <- PMtemp
  anova.morphology.across.time.MASTER[paste0("h",met),] <- HKtemp
  
}
# Change rownames

anova.morphology.across.time.MASTER.edit <- cbind(rep(c("Bray-Curtis","Weighted Unifrac","Un-weighted Unifrac"),2), anova.morphology.across.time.MASTER)
rownames(anova.morphology.across.time.MASTER.edit) <- c("Reed Point"," ","  ","Hakai","   ","    ")

capture.output(print(xtable(anova.morphology.across.time.MASTER.edit)
                            , math.style.negative = TRUE
                            , math.style.exponents = TRUE
                            , sanitize.text.function = function(x) {x}
                            )
                     , file = paste0("BETAPLOTS_LATEX/anova.morphology.across.time.MASTER.txt"))


# Combine for figure table

tempCombo <- rbind(anova.morphology.across.time.BC[1,]
      , pairwiseAdonis.all.BC[1:3,]
      , anova.morphology.across.time.BC[2,]
      , pairwiseAdonis.all.BC[4:6,])

LATEX_adonis_all_and_pairwise <- cbind(Sites = c("Reed Point","","","","Hakai","","","")
                                       , "Test type" = c("PERMANOVA", "" ,"Pairwise PERMANOVA","","PERMANOVA","","Pairwise PERMANOVA","")
                                       , tempCombo[,2:ncol(tempCombo)])
LATEX_adonis_all_and_pairwise <- gsub("^","\\\\makecell{", LATEX_adonis_all_and_pairwise)
LATEX_adonis_all_and_pairwise <- gsub("$","}", LATEX_adonis_all_and_pairwise)
LATEX_adonis_all_and_pairwise <- gsub("(","\\\\(", LATEX_adonis_all_and_pairwise, fixed = TRUE)


capture.output(print(xtable(LATEX_adonis_all_and_pairwise)
                     , include.rownames = FALSE
                     , math.style.negative = TRUE
                     , math.style.exponents = TRUE
                     , sanitize.text.function = function(x) {x}
                     )
               , file = paste0("BETAPLOTS_LATEX/LATEX_adonis_all_and_pairwise_BC.txt"))




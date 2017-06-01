#!/bin/RScript
set.seed(3)
# MAKE BETA PLOTS
library("MASS")
library("vegan")
library("nlme")
library("optparse")
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

BCPWD <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/ANALYSIS_ALPHABETATAXA/beta_div/bray_curtis_dm.txt'
UWUFPWD <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/ANALYSIS_ALPHABETATAXA/beta_div/unweighted_unifrac_dm.txt'
WUFPWD <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/ANALYSIS_ALPHABETATAXA/beta_div/weighted_unifrac_dm.txt'
MFPWD <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/OTU_MP_filt/MF_nochlpmito_m1000.txt'
setwd("/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis")
  
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
                 , row.names= 1)



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

MF <- MF[unlist(lapply(rownames(dm.UWUF), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
})),]

MF.P <- MF[unlist(lapply(rownames(dm.UWUF.P), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
})),]


# Get morph out-- also, 5760
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
system("mkdir BETAPLOTS_P")

#*** NOTE: ONLY WORKING ON UWUF RIGHT NOW; DO OTHERS AFTER
####### UWUF ############# 
metric <- "UWUF"

### NMDS #####

NMDS.UWUF.morphonly <- isoMDS(as.matrix(dm.UWUF.morphonly), y = cmdscale(as.matrix(dm.UWUF.morphonly), 2))
NMDS.UWUF.all <- isoMDS(as.matrix(dm.UWUF.inclWater), y = cmdscale(as.matrix(dm.UWUF.inclWater), 2))

###### STATS ##########
MF.morphkeep <- MF.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
dm.UWUF.morphonly.H <- dm.UWUF.morphonly[grep(".", rownames(dm.UWUF.morphonly), fixed = TRUE),grep(".", colnames(dm.UWUF.morphonly), fixed = TRUE)]
dm.UWUF.morphonly.P <- dm.UWUF.morphonly[grep("-", rownames(dm.UWUF.morphonly), fixed = TRUE),grep("-", colnames(dm.UWUF.morphonly), fixed = TRUE)]

ANOVA.UWUF.morphtime <- adonis(dm.UWUF.morphonly ~ Time*Morph, data = MF.morphkeep, by = "margin")
capture.output(ANOVA.UWUF.morphtime, file = paste0("BETAPLOTS/adonis_", metric,"_Hakai.txt"))
ANOSIM.UWUF.morphtime <- anosim(dm.UWUF.morphonly, grouping = MF.morphkeep$Morph)

# Dispersion across time and between morphs
dist.UWUF.morphonly <- as.dist(dm.UWUF.morphonly)
betadisp.UWUF.time <- betadisper(d = dist.UWUF.morphonly, group = MF.morphkeep$Time)
betadisp.UWUF.morph <- betadisper(d = dist.UWUF.morphonly, group = MF.morphkeep$Morph)

capture.output(betadisp.UWUF.time, file = paste0("BETAPLOTS/betadispTime_", metric, "_Hakai.txt"))
capture.output(betadisp.UWUF.morph, file = paste0("BETAPLOTS/betadispMorph_", metric, "_Hakai.txt"))


# Grep values for each timepoint and morph, and then plot dispersion within each
TP <- levels(factor(MF.morphkeep$Time))
MorphTypes <- levels(factor(MF.morphkeep$Morph))
# First Timepoint
FirstTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[1],"$"), MF.morphkeep$Time)]
# Between morphs
dm.UWUF.FirstTP <- dm.UWUF.morphonly[FirstTPNames,FirstTPNames]
FB.BL.Firstvalues <- as.vector(as.dist(dm.UWUF.FirstTP[grep("FB", rownames(dm.UWUF.FirstTP)), grep("BL", colnames(dm.UWUF.FirstTP))]))
FB.CR.Firstvalues <- as.vector(as.dist(dm.UWUF.FirstTP[grep("FB", rownames(dm.UWUF.FirstTP)), grep("CR", colnames(dm.UWUF.FirstTP))]))
CR.BL.Firstvalues <- as.vector(as.dist(dm.UWUF.FirstTP[grep("CR", rownames(dm.UWUF.FirstTP)), grep("BL", colnames(dm.UWUF.FirstTP))]))
# ANOVA
dm.UWUF.20 <- dm.UWUF.morphonly[grep("-20-", rownames(dm.UWUF.morphonly)), grep("-20-", colnames(dm.UWUF.morphonly))]
MF.UWUF.20.only <- MF.morphkeep[grep("-20-", rownames(MF.morphkeep)),]
ANOVA.UWUF.20.only <- adonis(dm.UWUF.20 ~ Morph, data = MF.UWUF.20.only, by = "margin")

# Second Timepoine
SecondTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[2],"$"), MF.morphkeep$Time)]
# Between morphs
dm.UWUF.SecondTP <- dm.UWUF.morphonly[SecondTPNames,SecondTPNames]
FB.BL.Secondvalues <- as.vector(as.dist(dm.UWUF.SecondTP[grep("FB", rownames(dm.UWUF.SecondTP)), grep("BL", colnames(dm.UWUF.SecondTP))]))
FB.CR.Secondvalues <- as.vector(as.dist(dm.UWUF.SecondTP[grep("FB", rownames(dm.UWUF.SecondTP)), grep("CR", colnames(dm.UWUF.SecondTP))]))
CR.BL.Secondvalues <- as.vector(as.dist(dm.UWUF.SecondTP[grep("CR", rownames(dm.UWUF.SecondTP)), grep("BL", colnames(dm.UWUF.SecondTP))]))
# ANOVA
dm.UWUF.60 <- dm.UWUF.morphonly[grep("-60-", rownames(dm.UWUF.morphonly)), grep("-60-", colnames(dm.UWUF.morphonly))]
MF.UWUF.60.only <- MF.morphkeep[grep("-60-", rownames(MF.morphkeep)),]
ANOVA.UWUF.60.only <- adonis(dm.UWUF.60 ~ Morph, data = MF.UWUF.60.only, by = "margin")

# Third Timepoint
ThirdTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[3],"$"), MF.morphkeep$Time)]
# Between morphs
dm.UWUF.ThirdTP <- dm.UWUF.morphonly[ThirdTPNames,ThirdTPNames]
FB.BL.Thirdvalues <- as.vector(as.dist(dm.UWUF.ThirdTP[grep("FB", rownames(dm.UWUF.ThirdTP)), grep("BL", colnames(dm.UWUF.ThirdTP))]))
FB.CR.Thirdvalues <- as.vector(as.dist(dm.UWUF.ThirdTP[grep("FB", rownames(dm.UWUF.ThirdTP)), grep("CR", colnames(dm.UWUF.ThirdTP))]))
CR.BL.Thirdvalues <- as.vector(as.dist(dm.UWUF.ThirdTP[grep("CR", rownames(dm.UWUF.ThirdTP)), grep("BL", colnames(dm.UWUF.ThirdTP))]))
# ANOVA
dm.UWUF.360 <- dm.UWUF.morphonly[grep("-360-", rownames(dm.UWUF.morphonly)), grep("-360-", colnames(dm.UWUF.morphonly))]
MF.UWUF.360.only <- MF.morphkeep[grep("-360-", rownames(MF.morphkeep)),]
ANOVA.UWUF.360.only <- adonis(dm.UWUF.360 ~ Morph, data = MF.UWUF.360.only, by = "margin")


# Fourth Timepoint
FourthTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[4],"$"), MF.morphkeep$Time)]
# Between morphs
dm.UWUF.FourthTP <- dm.UWUF.morphonly[FourthTPNames,FourthTPNames]
FB.BL.Fourthvalues <- as.vector(as.dist(dm.UWUF.FourthTP[grep("FB", rownames(dm.UWUF.FourthTP)), grep("BL", colnames(dm.UWUF.FourthTP))]))
FB.CR.Fourthvalues <- as.vector(as.dist(dm.UWUF.FourthTP[grep("FB", rownames(dm.UWUF.FourthTP)), grep("CR", colnames(dm.UWUF.FourthTP))]))
CR.BL.Fourthvalues <- as.vector(as.dist(dm.UWUF.FourthTP[grep("CR", rownames(dm.UWUF.FourthTP)), grep("BL", colnames(dm.UWUF.FourthTP))]))
# ANOVA
dm.UWUF.720 <- dm.UWUF.morphonly[grep("-720-", rownames(dm.UWUF.morphonly)), grep("-720-", colnames(dm.UWUF.morphonly))]
MF.UWUF.720.only <- MF.morphkeep[grep("-720-", rownames(MF.morphkeep)),]
ANOVA.UWUF.720.only <- adonis(dm.UWUF.720 ~ Morph, data = MF.UWUF.720.only, by = "margin")


# Combine into single tables

FB.BL.UWUF.ALL <- cbind(as.numeric(c(FB.BL.Firstvalues
                                     , FB.BL.Secondvalues
                                     , FB.BL.Thirdvalues
                                     , FB.BL.Fourthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.BL.Firstvalues))
                                       , rep(TP[2], length(FB.BL.Secondvalues))
                                       , rep(TP[3], length(FB.BL.Thirdvalues))
                                       , rep(TP[4], length(FB.BL.Fourthvalues)))))
FB.CR.UWUF.ALL <- cbind(as.numeric(c(FB.CR.Firstvalues
                                     , FB.CR.Secondvalues
                                     , FB.CR.Thirdvalues
                                     , FB.CR.Fourthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.CR.Firstvalues))
                                       , rep(TP[2], length(FB.CR.Secondvalues))
                                       , rep(TP[3], length(FB.CR.Thirdvalues))
                                       , rep(TP[4], length(FB.CR.Fourthvalues)))))
CR.BL.UWUF.ALL <- cbind(as.numeric(c(CR.BL.Firstvalues
                                     , CR.BL.Secondvalues
                                     , CR.BL.Thirdvalues
                                     , CR.BL.Fourthvalues))
                        , as.numeric(c(rep(TP[1], length(CR.BL.Firstvalues))
                                       , rep(TP[2], length(CR.BL.Secondvalues))
                                       , rep(TP[3], length(CR.BL.Thirdvalues))
                                       , rep(TP[4], length(CR.BL.Fourthvalues)))))

FB.BL.UWUF.lm <- lm(FB.BL.UWUF.ALL[,1] ~ log(FB.BL.UWUF.ALL[,2]))
summary(FB.BL.UWUF.lm)

CR.BL.UWUF.lm <- lm(CR.BL.UWUF.ALL[,1] ~ log(CR.BL.UWUF.ALL[,2]))
summary(CR.BL.UWUF.lm)

FB.CR.UWUF.lm <- lm(FB.CR.UWUF.ALL[,1] ~ log(FB.CR.UWUF.ALL[,2]))
summary(FB.CR.UWUF.lm)

# Plot the distances

FB.BL.UWUF.ALL.mean <- aggregate(FB.BL.UWUF.ALL, by = list(FB.BL.UWUF.ALL[,2]), mean)
FB.BL.UWUF.ALL.sd <- aggregate(FB.BL.UWUF.ALL, by = list(FB.BL.UWUF.ALL[,2]), sd)

CR.BL.UWUF.ALL.mean <- aggregate(CR.BL.UWUF.ALL, by = list(CR.BL.UWUF.ALL[,2]), mean)
CR.BL.UWUF.ALL.sd <- aggregate(CR.BL.UWUF.ALL, by = list(CR.BL.UWUF.ALL[,2]), sd)

FB.CR.UWUF.ALL.mean <- aggregate(FB.CR.UWUF.ALL, by = list(FB.CR.UWUF.ALL[,2]), mean)
FB.CR.UWUF.ALL.sd <- aggregate(FB.CR.UWUF.ALL, by = list(FB.CR.UWUF.ALL[,2]), sd)

######## PLOT DISP ###############
# ylimits <- c(min(FB.BL.UWUF.ALL[,1],FB.CR.UWUF.ALL[,1],CR.BL.UWUF.ALL[,1]), max(FB.BL.UWUF.ALL[,1],FB.CR.UWUF.ALL[,1],CR.BL.UWUF.ALL[,1]))
ylimits <- c(0.45,0.7)

# xvalues <- log(FB.BL.UWUF.ALL.mean[,1])
xvalues <- as.character(FB.BL.UWUF.ALL.mean[,1])
jpeg(paste0("BETAPLOTS/DispOverTime_",metric,".jpeg"))
par(fig = c(0,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1, at = c(1,2,3,4), labels = xvalues)
points(FB.BL.UWUF.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*0.99
       , x1 = c(1,2,3,4)*0.99
       , y0 = c(FB.BL.UWUF.ALL.mean[,2] - FB.BL.UWUF.ALL.sd[,2]/2)
       , y1 = c(FB.BL.UWUF.ALL.mean[,2] + FB.BL.UWUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.UWUF.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1
       , x1 = c(1,2,3,4)*1
       , y0 = c(FB.CR.UWUF.ALL.mean[,2] - FB.CR.UWUF.ALL.sd[,2]/2)
       , y1 = c(FB.CR.UWUF.ALL.mean[,2] + FB.CR.UWUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.UWUF.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1.01
       , x1 = c(1,2,3,4)*1.01
       , y0 = c(CR.BL.UWUF.ALL.mean[,2] - CR.BL.UWUF.ALL.sd[,2]/2)
       , y1 = c(CR.BL.UWUF.ALL.mean[,2] + CR.BL.UWUF.ALL.sd[,2]/2)
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
       , col = c("darkgreen","purple","grey"))
dev.off()

####### PLOT MORPH ############
# NOCHLP UWUF
MF.morphkeep$Morph <- factor(MF.morphkeep$Morph, levels = c('CR','BL','FB'))
MF.morphkeep$Time <- factor(MF.morphkeep$Time, levels = c('20','60','180','360','720','1440'))
MF.morphkeep$Type <- factor(MF.morphkeep$Type, levels = c('P','H'))
MorphColours <- c("red","magenta","blue")


# MAKE POLYGONS for plotting
NMDS.UWUF.CR <- NMDS.UWUF.morphonly$points[grep("CR", rownames(NMDS.UWUF.morphonly$points)),]
NMDS.UWUF.BL <- NMDS.UWUF.morphonly$points[grep("BL", rownames(NMDS.UWUF.morphonly$points)),]
NMDS.UWUF.FB <- NMDS.UWUF.morphonly$points[grep("FB", rownames(NMDS.UWUF.morphonly$points)),]

NMDS.UWUF.CR.chull <- chull(NMDS.UWUF.CR)
NMDS.UWUF.CR.chull <- c(NMDS.UWUF.CR.chull, NMDS.UWUF.CR.chull[1])

NMDS.UWUF.BL.chull <- chull(NMDS.UWUF.BL)
NMDS.UWUF.BL.chull <- c(NMDS.UWUF.BL.chull, NMDS.UWUF.BL.chull[1])

NMDS.UWUF.FB.chull <- chull(NMDS.UWUF.FB)
NMDS.UWUF.FB.chull <- c(NMDS.UWUF.FB.chull, NMDS.UWUF.FB.chull[1])


jpeg(paste0("BETAPLOTS/NMDS_",metric,"_Morph.jpeg"))
par(fig = c(0,0.8,0,1))
plot(NMDS.UWUF.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.morphkeep$Morph)]
     , sub = round(NMDS.UWUF.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
     )
lines(NMDS.UWUF.CR[NMDS.UWUF.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.BL[NMDS.UWUF.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.FB[NMDS.UWUF.FB.chull,]
      , col = "blue")
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
       )
dev.off()

####### PLOT TIME ############
TimeColours <- c("grey","lightblue","blue","darkblue")

jpeg(paste0("BETAPLOTS/NMDS_",metric,"_Time.jpeg"))
par(fig = c(0,0.8,0,1))
plot(NMDS.UWUF.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = TimeColours[factor(MF.morphkeep$Time)]
     , sub = round(NMDS.UWUF.morphonly$stress/100,2)
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



# TIME-- ONE at a time


# 20 minutes
NMDS.UWUF.20.only <- NMDS.UWUF.morphonly$points[grep("-20-", rownames(NMDS.UWUF.morphonly$points)),]

NMDS.UWUF.20.only.CR <- NMDS.UWUF.20.only[grep("CR", rownames(NMDS.UWUF.20.only)),]
NMDS.UWUF.20.only.CR.chull <- chull(NMDS.UWUF.20.only.CR)
NMDS.UWUF.20.only.CR.chull <- c(NMDS.UWUF.20.only.CR.chull, NMDS.UWUF.20.only.CR.chull[1])

NMDS.UWUF.20.only.BL <- NMDS.UWUF.20.only[grep("BL", rownames(NMDS.UWUF.20.only)),]
NMDS.UWUF.20.only.BL.chull <- chull(NMDS.UWUF.20.only.BL)
NMDS.UWUF.20.only.BL.chull <- c(NMDS.UWUF.20.only.BL.chull, NMDS.UWUF.20.only.BL.chull[1])

NMDS.UWUF.20.only.FB <- NMDS.UWUF.20.only[grep("FB", rownames(NMDS.UWUF.20.only)),]
NMDS.UWUF.20.only.FB.chull <- chull(NMDS.UWUF.20.only.FB)
NMDS.UWUF.20.only.FB.chull <- c(NMDS.UWUF.20.only.FB.chull, NMDS.UWUF.20.only.FB.chull[1])

# 60 minutes
NMDS.UWUF.60.only <- NMDS.UWUF.morphonly$points[grep("-60-", rownames(NMDS.UWUF.morphonly$points)),]


NMDS.UWUF.60.only.CR <- NMDS.UWUF.60.only[grep("CR", rownames(NMDS.UWUF.60.only)),]
NMDS.UWUF.60.only.CR.chull <- chull(NMDS.UWUF.60.only.CR)
NMDS.UWUF.60.only.CR.chull <- c(NMDS.UWUF.60.only.CR.chull, NMDS.UWUF.60.only.CR.chull[1])

NMDS.UWUF.60.only.BL <- NMDS.UWUF.60.only[grep("BL", rownames(NMDS.UWUF.60.only)),]
NMDS.UWUF.60.only.BL.chull <- chull(NMDS.UWUF.60.only.BL)
NMDS.UWUF.60.only.BL.chull <- c(NMDS.UWUF.60.only.BL.chull, NMDS.UWUF.60.only.BL.chull[1])

NMDS.UWUF.60.only.FB <- NMDS.UWUF.60.only[grep("FB", rownames(NMDS.UWUF.60.only)),]
NMDS.UWUF.60.only.FB.chull <- chull(NMDS.UWUF.60.only.FB)
NMDS.UWUF.60.only.FB.chull <- c(NMDS.UWUF.60.only.FB.chull, NMDS.UWUF.60.only.FB.chull[1])



# 360 minutes
NMDS.UWUF.360.only <- NMDS.UWUF.morphonly$points[grep("-360-", rownames(NMDS.UWUF.morphonly$points)),]


NMDS.UWUF.360.only.CR <- NMDS.UWUF.360.only[grep("CR", rownames(NMDS.UWUF.360.only)),]
NMDS.UWUF.360.only.CR.chull <- chull(NMDS.UWUF.360.only.CR)
NMDS.UWUF.360.only.CR.chull <- c(NMDS.UWUF.360.only.CR.chull, NMDS.UWUF.360.only.CR.chull[1])

NMDS.UWUF.360.only.BL <- NMDS.UWUF.360.only[grep("BL", rownames(NMDS.UWUF.360.only)),]
NMDS.UWUF.360.only.BL.chull <- chull(NMDS.UWUF.360.only.BL)
NMDS.UWUF.360.only.BL.chull <- c(NMDS.UWUF.360.only.BL.chull, NMDS.UWUF.360.only.BL.chull[1])

NMDS.UWUF.360.only.FB <- NMDS.UWUF.360.only[grep("FB", rownames(NMDS.UWUF.360.only)),]
NMDS.UWUF.360.only.FB.chull <- chull(NMDS.UWUF.360.only.FB)
NMDS.UWUF.360.only.FB.chull <- c(NMDS.UWUF.360.only.FB.chull, NMDS.UWUF.360.only.FB.chull[1])


# 720 minutes
NMDS.UWUF.720.only <- NMDS.UWUF.morphonly$points[grep("-720-", rownames(NMDS.UWUF.morphonly$points)),]


NMDS.UWUF.720.only.CR <- NMDS.UWUF.720.only[grep("CR", rownames(NMDS.UWUF.720.only)),]
NMDS.UWUF.720.only.CR.chull <- chull(NMDS.UWUF.720.only.CR)
NMDS.UWUF.720.only.CR.chull <- c(NMDS.UWUF.720.only.CR.chull, NMDS.UWUF.720.only.CR.chull[1])

NMDS.UWUF.720.only.BL <- NMDS.UWUF.720.only[grep("BL", rownames(NMDS.UWUF.720.only)),]
NMDS.UWUF.720.only.BL.chull <- chull(NMDS.UWUF.720.only.BL)
NMDS.UWUF.720.only.BL.chull <- c(NMDS.UWUF.720.only.BL.chull, NMDS.UWUF.720.only.BL.chull[1])

NMDS.UWUF.720.only.FB <- NMDS.UWUF.720.only[grep("FB", rownames(NMDS.UWUF.720.only)),]
NMDS.UWUF.720.only.FB.chull <- chull(NMDS.UWUF.720.only.FB)
NMDS.UWUF.720.only.FB.chull <- c(NMDS.UWUF.720.only.FB.chull, NMDS.UWUF.720.only.FB.chull[1])


jpeg(paste0("BETAPLOTS/NMDS_",metric,"_TimebyMorph.jpeg"), width = 2000, height = 800, pointsize = 24)
par(mfrow= c(1,4), oma = c(4,4,4,6))
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.20.only
     , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.20.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.20.only.CR[NMDS.UWUF.20.only.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.20.only.BL[NMDS.UWUF.20.only.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.20.only.FB[NMDS.UWUF.20.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.60.only
     , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.UWUF.60.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.60.only.CR[NMDS.UWUF.60.only.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.60.only.BL[NMDS.UWUF.60.only.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.60.only.FB[NMDS.UWUF.60.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.360.only
     , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.360.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.360.only.CR[NMDS.UWUF.360.only.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.360.only.BL[NMDS.UWUF.360.only.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.360.only.FB[NMDS.UWUF.360.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.720.only
     , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.720.only.CR[NMDS.UWUF.720.only.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.720.only.BL[NMDS.UWUF.720.only.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.720.only.FB[NMDS.UWUF.720.only.FB.chull,]
      , col = "blue")
#STOP
dev.off()

############ COMBO DISP BETA ################
# Disp and beta through time combined
ylimits <- c(0.45,0.7)
# xvalues <- log(FB.BL.UWUF.ALL.mean[,1])
xvalues <- as.character(FB.BL.UWUF.ALL.mean[,1])

jpeg(paste0("BETAPLOTS/COMBO_dispbeta_",metric,".jpeg"), width = 1000, height = 800, pointsize = 16)
par(fig = c(0.5,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1
     , at = c(1,2,3,4)
     , labels = c("20 min","1 h","6 h","12 h")
     , las = 2)
points(FB.BL.UWUF.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*0.99
       , x1 = c(1,2,3,4)*0.99
       , y0 = c(FB.BL.UWUF.ALL.mean[,2] - FB.BL.UWUF.ALL.sd[,2]/2)
       , y1 = c(FB.BL.UWUF.ALL.mean[,2] + FB.BL.UWUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.UWUF.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1
       , x1 = c(1,2,3,4)*1
       , y0 = c(FB.CR.UWUF.ALL.mean[,2] - FB.CR.UWUF.ALL.sd[,2]/2)
       , y1 = c(FB.CR.UWUF.ALL.mean[,2] + FB.CR.UWUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.UWUF.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1.01
       , x1 = c(1,2,3,4)*1.01
       , y0 = c(CR.BL.UWUF.ALL.mean[,2] - CR.BL.UWUF.ALL.sd[,2]/2)
       , y1 = c(CR.BL.UWUF.ALL.mean[,2] + CR.BL.UWUF.ALL.sd[,2]/2)
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
       , col = c("darkgreen","purple","grey"))
par(fig = c(0,0.5,0.8,1), mfrow= c(1,4), oma = c(4,4,4,6), new = TRUE)
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.20.only
     # , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.20.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.20.only.CR[NMDS.UWUF.20.only.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.20.only.BL[NMDS.UWUF.20.only.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.20.only.FB[NMDS.UWUF.20.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.60.only
     # , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.UWUF.60.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.60.only.CR[NMDS.UWUF.60.only.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.60.only.BL[NMDS.UWUF.60.only.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.60.only.FB[NMDS.UWUF.60.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.360.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.360.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.360.only.CR[NMDS.UWUF.360.only.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.360.only.BL[NMDS.UWUF.360.only.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.360.only.FB[NMDS.UWUF.360.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.720.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.720.only.CR[NMDS.UWUF.720.only.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.720.only.BL[NMDS.UWUF.720.only.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.720.only.FB[NMDS.UWUF.720.only.FB.chull,]
      , col = "blue")
#STOP
dev.off()

############ PLOT 5760 ################
dm.UWUF.5760<- dm.UWUF.inclWater[grep("(CR|FB|BL)-5760", rownames(dm.UWUF.inclWater)),grep("(CR|FB|BL)-5760", colnames(dm.UWUF.inclWater))]

NMDS.UWUF.5760 <- isoMDS(as.matrix(dm.UWUF.5760), y = cmdscale(as.matrix(dm.UWUF.5760), 2))

# NMDS.UWUF.morphallTimePoints <- NMDS.UWUF.morphAllTime$points[grep("(CR|BL|FB)-5760", rownames(NMDS.UWUF.morphAllTime$points)),]
MF.5760 <- MF[sapply(rownames(NMDS.UWUF.5760$points), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
}),]
MF.5760$Morph <- factor(MF.5760$Morph, levels = c("CR","BL","FB"))
MorphColours <- c("red","magenta","blue")

NMDS.UWUF.5760.CR <- NMDS.UWUF.5760$points[grep("CR", rownames(NMDS.UWUF.5760$points)),]
NMDS.UWUF.5760.CR.chull <- chull(NMDS.UWUF.5760.CR)
NMDS.UWUF.5760.CR.chull <- c(NMDS.UWUF.5760.CR.chull, NMDS.UWUF.5760.CR.chull[1])

NMDS.UWUF.5760.BL <- NMDS.UWUF.5760$points[grep("BL", rownames(NMDS.UWUF.5760$points)),]
NMDS.UWUF.5760.BL.chull <- chull(NMDS.UWUF.5760.BL)
NMDS.UWUF.5760.BL.chull <- c(NMDS.UWUF.5760.BL.chull, NMDS.UWUF.5760.BL.chull[1])

NMDS.UWUF.5760.FB <- NMDS.UWUF.5760$points[grep("FB", rownames(NMDS.UWUF.5760$points)),]
NMDS.UWUF.5760.FB.chull <- chull(NMDS.UWUF.5760.FB)
NMDS.UWUF.5760.FB.chull <- c(NMDS.UWUF.5760.FB.chull, NMDS.UWUF.5760.FB.chull[1])

ANOVA.UWUF.5760 <- adonis(dm.UWUF.5760 ~ Morph, data = MF.5760)
ANOSIM.UWUF.5760 <- anosim(dat = dm.UWUF.5760, grouping = MF.5760$Morph)
capture.output(ANOVA.UWUF.5760, file = paste0("BETAPLOTS/anova_",metric,"_5760only.txt"))

jpeg(paste0("BETAPLOTS/NMDS_",metric,"_5760Only.jpeg"), width = 1000, height = 800, pointsize = 24)
par(fig = c(0,0.7,0,1))
plot(NMDS.UWUF.5760$points
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.5760$Morph)]
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , sub = paste0("Stress: ", round(NMDS.UWUF.5760$stress,2)/100)
     , cex = 1.5
     )
lines(NMDS.UWUF.5760.CR[NMDS.UWUF.5760.CR.chull,]
      , col = "red")
lines(NMDS.UWUF.5760.BL[NMDS.UWUF.5760.BL.chull,]
      , col = "magenta")
lines(NMDS.UWUF.5760.FB[NMDS.UWUF.5760.FB.chull,]
      , col = "blue")
par(fig = c(0.7,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , xlab = ""
     , ylab = ""
     , xaxt = "n"
     , yaxt = "n"
     , bty = "n"
     , pch = "")
legend("center"
       , legend = c("Crust","Blade","Finely Br.")
       , pch = 21
       , col= "black"
       , pt.bg = c("red","magenta","blue"))
dev.off()


####### WUF #############
metric <- "WUF"

### NMDS #####

NMDS.WUF.morphonly <- isoMDS(as.matrix(dm.WUF.morphonly), y = cmdscale(as.matrix(dm.WUF.morphonly), 2))
NMDS.WUF.all <- isoMDS(as.matrix(dm.WUF.inclWater), y = cmdscale(as.matrix(dm.WUF.inclWater), 2))

###### STATS ##########
MF.morphkeep <- MF.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
dm.WUF.morphonly.H <- dm.WUF.morphonly[grep(".", rownames(dm.WUF.morphonly), fixed = TRUE),grep(".", colnames(dm.WUF.morphonly), fixed = TRUE)]
dm.WUF.morphonly.P <- dm.WUF.morphonly[grep("-", rownames(dm.WUF.morphonly), fixed = TRUE),grep("-", colnames(dm.WUF.morphonly), fixed = TRUE)]

ANOVA.WUF.morphtime <- adonis(dm.WUF.morphonly ~ Time*Morph, data = MF.morphkeep, by = "margin")
capture.output(ANOVA.WUF.morphtime, file = paste0("BETAPLOTS/adonis_", metric,"_Hakai.txt"))
ANOSIM.WUF.morphtime <- anosim(dm.WUF.morphonly, grouping = MF.morphkeep$Morph)

# Dispersion across time and between morphs
dist.WUF.morphonly <- as.dist(dm.WUF.morphonly)
betadisp.WUF.time <- betadisper(d = dist.WUF.morphonly, group = MF.morphkeep$Time)
betadisp.WUF.morph <- betadisper(d = dist.WUF.morphonly, group = MF.morphkeep$Morph)

capture.output(betadisp.WUF.time, file = paste0("BETAPLOTS/betadispTime_", metric, "_Hakai.txt"))
capture.output(betadisp.WUF.morph, file = paste0("BETAPLOTS/betadispMorph_", metric, "_Hakai.txt"))


# Grep values for each timepoint and morph, and then plot dispersion within each
TP <- levels(factor(MF.morphkeep$Time))
MorphTypes <- levels(factor(MF.morphkeep$Morph))
# First Timepoint
FirstTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[1],"$"), MF.morphkeep$Time)]
# Between morphs
dm.WUF.FirstTP <- dm.WUF.morphonly[FirstTPNames,FirstTPNames]
FB.BL.Firstvalues <- as.vector(as.dist(dm.WUF.FirstTP[grep("FB", rownames(dm.WUF.FirstTP)), grep("BL", colnames(dm.WUF.FirstTP))]))
FB.CR.Firstvalues <- as.vector(as.dist(dm.WUF.FirstTP[grep("FB", rownames(dm.WUF.FirstTP)), grep("CR", colnames(dm.WUF.FirstTP))]))
CR.BL.Firstvalues <- as.vector(as.dist(dm.WUF.FirstTP[grep("CR", rownames(dm.WUF.FirstTP)), grep("BL", colnames(dm.WUF.FirstTP))]))
# ANOVA
dm.WUF.20 <- dm.WUF.morphonly[grep("-20-", rownames(dm.WUF.morphonly)), grep("-20-", colnames(dm.WUF.morphonly))]
MF.WUF.20.only <- MF.morphkeep[grep("-20-", rownames(MF.morphkeep)),]
ANOVA.WUF.20.only <- adonis(dm.WUF.20 ~ Morph, data = MF.WUF.20.only, by = "margin")

# Second Timepoine
SecondTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[2],"$"), MF.morphkeep$Time)]
# Between morphs
dm.WUF.SecondTP <- dm.WUF.morphonly[SecondTPNames,SecondTPNames]
FB.BL.Secondvalues <- as.vector(as.dist(dm.WUF.SecondTP[grep("FB", rownames(dm.WUF.SecondTP)), grep("BL", colnames(dm.WUF.SecondTP))]))
FB.CR.Secondvalues <- as.vector(as.dist(dm.WUF.SecondTP[grep("FB", rownames(dm.WUF.SecondTP)), grep("CR", colnames(dm.WUF.SecondTP))]))
CR.BL.Secondvalues <- as.vector(as.dist(dm.WUF.SecondTP[grep("CR", rownames(dm.WUF.SecondTP)), grep("BL", colnames(dm.WUF.SecondTP))]))
# ANOVA
dm.WUF.60 <- dm.WUF.morphonly[grep("-60-", rownames(dm.WUF.morphonly)), grep("-60-", colnames(dm.WUF.morphonly))]
MF.WUF.60.only <- MF.morphkeep[grep("-60-", rownames(MF.morphkeep)),]
ANOVA.WUF.60.only <- adonis(dm.WUF.60 ~ Morph, data = MF.WUF.60.only, by = "margin")

# Third Timepoint
ThirdTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[3],"$"), MF.morphkeep$Time)]
# Between morphs
dm.WUF.ThirdTP <- dm.WUF.morphonly[ThirdTPNames,ThirdTPNames]
FB.BL.Thirdvalues <- as.vector(as.dist(dm.WUF.ThirdTP[grep("FB", rownames(dm.WUF.ThirdTP)), grep("BL", colnames(dm.WUF.ThirdTP))]))
FB.CR.Thirdvalues <- as.vector(as.dist(dm.WUF.ThirdTP[grep("FB", rownames(dm.WUF.ThirdTP)), grep("CR", colnames(dm.WUF.ThirdTP))]))
CR.BL.Thirdvalues <- as.vector(as.dist(dm.WUF.ThirdTP[grep("CR", rownames(dm.WUF.ThirdTP)), grep("BL", colnames(dm.WUF.ThirdTP))]))
# ANOVA
dm.WUF.360 <- dm.WUF.morphonly[grep("-360-", rownames(dm.WUF.morphonly)), grep("-360-", colnames(dm.WUF.morphonly))]
MF.WUF.360.only <- MF.morphkeep[grep("-360-", rownames(MF.morphkeep)),]
ANOVA.WUF.360.only <- adonis(dm.WUF.360 ~ Morph, data = MF.WUF.360.only, by = "margin")


# Fourth Timepoint
FourthTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[4],"$"), MF.morphkeep$Time)]
# Between morphs
dm.WUF.FourthTP <- dm.WUF.morphonly[FourthTPNames,FourthTPNames]
FB.BL.Fourthvalues <- as.vector(as.dist(dm.WUF.FourthTP[grep("FB", rownames(dm.WUF.FourthTP)), grep("BL", colnames(dm.WUF.FourthTP))]))
FB.CR.Fourthvalues <- as.vector(as.dist(dm.WUF.FourthTP[grep("FB", rownames(dm.WUF.FourthTP)), grep("CR", colnames(dm.WUF.FourthTP))]))
CR.BL.Fourthvalues <- as.vector(as.dist(dm.WUF.FourthTP[grep("CR", rownames(dm.WUF.FourthTP)), grep("BL", colnames(dm.WUF.FourthTP))]))
# ANOVA
dm.WUF.720 <- dm.WUF.morphonly[grep("-720-", rownames(dm.WUF.morphonly)), grep("-720-", colnames(dm.WUF.morphonly))]
MF.WUF.720.only <- MF.morphkeep[grep("-720-", rownames(MF.morphkeep)),]
ANOVA.WUF.720.only <- adonis(dm.WUF.720 ~ Morph, data = MF.WUF.720.only, by = "margin")


# Combine into single tables

FB.BL.WUF.ALL <- cbind(as.numeric(c(FB.BL.Firstvalues
                                     , FB.BL.Secondvalues
                                     , FB.BL.Thirdvalues
                                     , FB.BL.Fourthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.BL.Firstvalues))
                                       , rep(TP[2], length(FB.BL.Secondvalues))
                                       , rep(TP[3], length(FB.BL.Thirdvalues))
                                       , rep(TP[4], length(FB.BL.Fourthvalues)))))
FB.CR.WUF.ALL <- cbind(as.numeric(c(FB.CR.Firstvalues
                                     , FB.CR.Secondvalues
                                     , FB.CR.Thirdvalues
                                     , FB.CR.Fourthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.CR.Firstvalues))
                                       , rep(TP[2], length(FB.CR.Secondvalues))
                                       , rep(TP[3], length(FB.CR.Thirdvalues))
                                       , rep(TP[4], length(FB.CR.Fourthvalues)))))
CR.BL.WUF.ALL <- cbind(as.numeric(c(CR.BL.Firstvalues
                                     , CR.BL.Secondvalues
                                     , CR.BL.Thirdvalues
                                     , CR.BL.Fourthvalues))
                        , as.numeric(c(rep(TP[1], length(CR.BL.Firstvalues))
                                       , rep(TP[2], length(CR.BL.Secondvalues))
                                       , rep(TP[3], length(CR.BL.Thirdvalues))
                                       , rep(TP[4], length(CR.BL.Fourthvalues)))))

FB.BL.WUF.lm <- lm(FB.BL.WUF.ALL[,1] ~ log(FB.BL.WUF.ALL[,2]))
summary(FB.BL.WUF.lm)

CR.BL.WUF.lm <- lm(CR.BL.WUF.ALL[,1] ~ log(CR.BL.WUF.ALL[,2]))
summary(CR.BL.WUF.lm)

FB.CR.WUF.lm <- lm(FB.CR.WUF.ALL[,1] ~ log(FB.CR.WUF.ALL[,2]))
summary(FB.CR.WUF.lm)

# Plot the distances

FB.BL.WUF.ALL.mean <- aggregate(FB.BL.WUF.ALL, by = list(FB.BL.WUF.ALL[,2]), mean)
FB.BL.WUF.ALL.sd <- aggregate(FB.BL.WUF.ALL, by = list(FB.BL.WUF.ALL[,2]), sd)

CR.BL.WUF.ALL.mean <- aggregate(CR.BL.WUF.ALL, by = list(CR.BL.WUF.ALL[,2]), mean)
CR.BL.WUF.ALL.sd <- aggregate(CR.BL.WUF.ALL, by = list(CR.BL.WUF.ALL[,2]), sd)

FB.CR.WUF.ALL.mean <- aggregate(FB.CR.WUF.ALL, by = list(FB.CR.WUF.ALL[,2]), mean)
FB.CR.WUF.ALL.sd <- aggregate(FB.CR.WUF.ALL, by = list(FB.CR.WUF.ALL[,2]), sd)

######## PLOT DISP ###############
# ylimits <- c(min(FB.BL.WUF.ALL[,1],FB.CR.WUF.ALL[,1],CR.BL.WUF.ALL[,1]), max(FB.BL.WUF.ALL[,1],FB.CR.WUF.ALL[,1],CR.BL.WUF.ALL[,1]))
ylimits <- c(0.1,0.4)

# xvalues <- log(FB.BL.WUF.ALL.mean[,1])
xvalues <- as.character(FB.BL.WUF.ALL.mean[,1])
jpeg(paste0("BETAPLOTS/DispOverTime_",metric,".jpeg"))
par(fig = c(0,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1, at = c(1,2,3,4), labels = xvalues)
points(FB.BL.WUF.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*0.99
       , x1 = c(1,2,3,4)*0.99
       , y0 = c(FB.BL.WUF.ALL.mean[,2] - FB.BL.WUF.ALL.sd[,2]/2)
       , y1 = c(FB.BL.WUF.ALL.mean[,2] + FB.BL.WUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.WUF.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1
       , x1 = c(1,2,3,4)*1
       , y0 = c(FB.CR.WUF.ALL.mean[,2] - FB.CR.WUF.ALL.sd[,2]/2)
       , y1 = c(FB.CR.WUF.ALL.mean[,2] + FB.CR.WUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.WUF.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1.01
       , x1 = c(1,2,3,4)*1.01
       , y0 = c(CR.BL.WUF.ALL.mean[,2] - CR.BL.WUF.ALL.sd[,2]/2)
       , y1 = c(CR.BL.WUF.ALL.mean[,2] + CR.BL.WUF.ALL.sd[,2]/2)
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
       , col = c("darkgreen","purple","grey"))
dev.off()

####### PLOT MORPH ############
# NOCHLP WUF
MF.morphkeep$Morph <- factor(MF.morphkeep$Morph, levels = c('CR','BL','FB'))
MF.morphkeep$Time <- factor(MF.morphkeep$Time, levels = c('20','60','180','360','720','1440'))
MF.morphkeep$Type <- factor(MF.morphkeep$Type, levels = c('P','H'))
MorphColours <- c("red","magenta","blue")


# MAKE POLYGONS for plotting
NMDS.WUF.CR <- NMDS.WUF.morphonly$points[grep("CR", rownames(NMDS.WUF.morphonly$points)),]
NMDS.WUF.BL <- NMDS.WUF.morphonly$points[grep("BL", rownames(NMDS.WUF.morphonly$points)),]
NMDS.WUF.FB <- NMDS.WUF.morphonly$points[grep("FB", rownames(NMDS.WUF.morphonly$points)),]

NMDS.WUF.CR.chull <- chull(NMDS.WUF.CR)
NMDS.WUF.CR.chull <- c(NMDS.WUF.CR.chull, NMDS.WUF.CR.chull[1])

NMDS.WUF.BL.chull <- chull(NMDS.WUF.BL)
NMDS.WUF.BL.chull <- c(NMDS.WUF.BL.chull, NMDS.WUF.BL.chull[1])

NMDS.WUF.FB.chull <- chull(NMDS.WUF.FB)
NMDS.WUF.FB.chull <- c(NMDS.WUF.FB.chull, NMDS.WUF.FB.chull[1])


jpeg(paste0("BETAPLOTS/NMDS_",metric,"_Morph.jpeg"))
par(fig = c(0,0.8,0,1))
plot(NMDS.WUF.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.morphkeep$Morph)]
     , sub = round(NMDS.WUF.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
)
lines(NMDS.WUF.CR[NMDS.WUF.CR.chull,]
      , col = "red")
lines(NMDS.WUF.BL[NMDS.WUF.BL.chull,]
      , col = "magenta")
lines(NMDS.WUF.FB[NMDS.WUF.FB.chull,]
      , col = "blue")
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
)
dev.off()

####### PLOT TIME ############
TimeColours <- c("grey","lightblue","blue","darkblue")

jpeg(paste0("BETAPLOTS/NMDS_",metric,"_Time.jpeg"))
par(fig = c(0,0.8,0,1))
plot(NMDS.WUF.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = TimeColours[factor(MF.morphkeep$Time)]
     , sub = round(NMDS.WUF.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
)
# lines(NMDS.WUF.CR[NMDS.WUF.CR.chull,]
#       , col = "red")
# lines(NMDS.WUF.BL[NMDS.WUF.BL.chull,]
#       , col = "magenta")
# lines(NMDS.WUF.FB[NMDS.WUF.FB.chull,]
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



# TIME-- ONE at a time


# 20 minutes
NMDS.WUF.20.only <- NMDS.WUF.morphonly$points[grep("-20-", rownames(NMDS.WUF.morphonly$points)),]

NMDS.WUF.20.only.CR <- NMDS.WUF.20.only[grep("CR", rownames(NMDS.WUF.20.only)),]
NMDS.WUF.20.only.CR.chull <- chull(NMDS.WUF.20.only.CR)
NMDS.WUF.20.only.CR.chull <- c(NMDS.WUF.20.only.CR.chull, NMDS.WUF.20.only.CR.chull[1])

NMDS.WUF.20.only.BL <- NMDS.WUF.20.only[grep("BL", rownames(NMDS.WUF.20.only)),]
NMDS.WUF.20.only.BL.chull <- chull(NMDS.WUF.20.only.BL)
NMDS.WUF.20.only.BL.chull <- c(NMDS.WUF.20.only.BL.chull, NMDS.WUF.20.only.BL.chull[1])

NMDS.WUF.20.only.FB <- NMDS.WUF.20.only[grep("FB", rownames(NMDS.WUF.20.only)),]
NMDS.WUF.20.only.FB.chull <- chull(NMDS.WUF.20.only.FB)
NMDS.WUF.20.only.FB.chull <- c(NMDS.WUF.20.only.FB.chull, NMDS.WUF.20.only.FB.chull[1])

# 60 minutes
NMDS.WUF.60.only <- NMDS.WUF.morphonly$points[grep("-60-", rownames(NMDS.WUF.morphonly$points)),]


NMDS.WUF.60.only.CR <- NMDS.WUF.60.only[grep("CR", rownames(NMDS.WUF.60.only)),]
NMDS.WUF.60.only.CR.chull <- chull(NMDS.WUF.60.only.CR)
NMDS.WUF.60.only.CR.chull <- c(NMDS.WUF.60.only.CR.chull, NMDS.WUF.60.only.CR.chull[1])

NMDS.WUF.60.only.BL <- NMDS.WUF.60.only[grep("BL", rownames(NMDS.WUF.60.only)),]
NMDS.WUF.60.only.BL.chull <- chull(NMDS.WUF.60.only.BL)
NMDS.WUF.60.only.BL.chull <- c(NMDS.WUF.60.only.BL.chull, NMDS.WUF.60.only.BL.chull[1])

NMDS.WUF.60.only.FB <- NMDS.WUF.60.only[grep("FB", rownames(NMDS.WUF.60.only)),]
NMDS.WUF.60.only.FB.chull <- chull(NMDS.WUF.60.only.FB)
NMDS.WUF.60.only.FB.chull <- c(NMDS.WUF.60.only.FB.chull, NMDS.WUF.60.only.FB.chull[1])



# 360 minutes
NMDS.WUF.360.only <- NMDS.WUF.morphonly$points[grep("-360-", rownames(NMDS.WUF.morphonly$points)),]


NMDS.WUF.360.only.CR <- NMDS.WUF.360.only[grep("CR", rownames(NMDS.WUF.360.only)),]
NMDS.WUF.360.only.CR.chull <- chull(NMDS.WUF.360.only.CR)
NMDS.WUF.360.only.CR.chull <- c(NMDS.WUF.360.only.CR.chull, NMDS.WUF.360.only.CR.chull[1])

NMDS.WUF.360.only.BL <- NMDS.WUF.360.only[grep("BL", rownames(NMDS.WUF.360.only)),]
NMDS.WUF.360.only.BL.chull <- chull(NMDS.WUF.360.only.BL)
NMDS.WUF.360.only.BL.chull <- c(NMDS.WUF.360.only.BL.chull, NMDS.WUF.360.only.BL.chull[1])

NMDS.WUF.360.only.FB <- NMDS.WUF.360.only[grep("FB", rownames(NMDS.WUF.360.only)),]
NMDS.WUF.360.only.FB.chull <- chull(NMDS.WUF.360.only.FB)
NMDS.WUF.360.only.FB.chull <- c(NMDS.WUF.360.only.FB.chull, NMDS.WUF.360.only.FB.chull[1])


# 720 minutes
NMDS.WUF.720.only <- NMDS.WUF.morphonly$points[grep("-720-", rownames(NMDS.WUF.morphonly$points)),]


NMDS.WUF.720.only.CR <- NMDS.WUF.720.only[grep("CR", rownames(NMDS.WUF.720.only)),]
NMDS.WUF.720.only.CR.chull <- chull(NMDS.WUF.720.only.CR)
NMDS.WUF.720.only.CR.chull <- c(NMDS.WUF.720.only.CR.chull, NMDS.WUF.720.only.CR.chull[1])

NMDS.WUF.720.only.BL <- NMDS.WUF.720.only[grep("BL", rownames(NMDS.WUF.720.only)),]
NMDS.WUF.720.only.BL.chull <- chull(NMDS.WUF.720.only.BL)
NMDS.WUF.720.only.BL.chull <- c(NMDS.WUF.720.only.BL.chull, NMDS.WUF.720.only.BL.chull[1])

NMDS.WUF.720.only.FB <- NMDS.WUF.720.only[grep("FB", rownames(NMDS.WUF.720.only)),]
NMDS.WUF.720.only.FB.chull <- chull(NMDS.WUF.720.only.FB)
NMDS.WUF.720.only.FB.chull <- c(NMDS.WUF.720.only.FB.chull, NMDS.WUF.720.only.FB.chull[1])


jpeg(paste0("BETAPLOTS/NMDS_",metric,"_TimebyMorph.jpeg"), width = 2000, height = 800, pointsize = 24)
par(mfrow= c(1,4), oma = c(4,4,4,6))
par(mar = c(4,0,4,0))
plot(NMDS.WUF.20.only
     , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.WUF.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.20.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.20.only.CR[NMDS.WUF.20.only.CR.chull,]
      , col = "red")
lines(NMDS.WUF.20.only.BL[NMDS.WUF.20.only.BL.chull,]
      , col = "magenta")
lines(NMDS.WUF.20.only.FB[NMDS.WUF.20.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.WUF.60.only
     , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.WUF.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.WUF.60.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.60.only.CR[NMDS.WUF.60.only.CR.chull,]
      , col = "red")
lines(NMDS.WUF.60.only.BL[NMDS.WUF.60.only.BL.chull,]
      , col = "magenta")
lines(NMDS.WUF.60.only.FB[NMDS.WUF.60.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.WUF.360.only
     , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.WUF.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.360.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.360.only.CR[NMDS.WUF.360.only.CR.chull,]
      , col = "red")
lines(NMDS.WUF.360.only.BL[NMDS.WUF.360.only.BL.chull,]
      , col = "magenta")
lines(NMDS.WUF.360.only.FB[NMDS.WUF.360.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.WUF.720.only
     , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.WUF.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.720.only.CR[NMDS.WUF.720.only.CR.chull,]
      , col = "red")
lines(NMDS.WUF.720.only.BL[NMDS.WUF.720.only.BL.chull,]
      , col = "magenta")
lines(NMDS.WUF.720.only.FB[NMDS.WUF.720.only.FB.chull,]
      , col = "blue")
#STOP
dev.off()


############ PLOT 5760 ################
dm.WUF.5760<- dm.WUF.inclWater[grep("(CR|FB|BL)-5760", rownames(dm.WUF.inclWater)),grep("(CR|FB|BL)-5760", colnames(dm.WUF.inclWater))]

NMDS.WUF.5760 <- isoMDS(as.matrix(dm.WUF.5760), y = cmdscale(as.matrix(dm.WUF.5760), 2))

# NMDS.WUF.morphallTimePoints <- NMDS.WUF.morphAllTime$points[grep("(CR|BL|FB)-5760", rownames(NMDS.WUF.morphAllTime$points)),]
MF.5760 <- MF[sapply(rownames(NMDS.WUF.5760$points), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
}),]
MF.5760$Morph <- factor(MF.5760$Morph, levels = c("CR","BL","FB"))
MorphColours <- c("red","magenta","blue")

NMDS.WUF.5760.CR <- NMDS.WUF.5760$points[grep("CR", rownames(NMDS.WUF.5760$points)),]
NMDS.WUF.5760.CR.chull <- chull(NMDS.WUF.5760.CR)
NMDS.WUF.5760.CR.chull <- c(NMDS.WUF.5760.CR.chull, NMDS.WUF.5760.CR.chull[1])

NMDS.WUF.5760.BL <- NMDS.WUF.5760$points[grep("BL", rownames(NMDS.WUF.5760$points)),]
NMDS.WUF.5760.BL.chull <- chull(NMDS.WUF.5760.BL)
NMDS.WUF.5760.BL.chull <- c(NMDS.WUF.5760.BL.chull, NMDS.WUF.5760.BL.chull[1])

NMDS.WUF.5760.FB <- NMDS.WUF.5760$points[grep("FB", rownames(NMDS.WUF.5760$points)),]
NMDS.WUF.5760.FB.chull <- chull(NMDS.WUF.5760.FB)
NMDS.WUF.5760.FB.chull <- c(NMDS.WUF.5760.FB.chull, NMDS.WUF.5760.FB.chull[1])

ANOVA.WUF.5760 <- adonis(dm.WUF.5760 ~ Morph, data = MF.5760)
ANOSIM.WUF.5760 <- anosim(dat = dm.WUF.5760, grouping = MF.5760$Morph)
capture.output(ANOVA.WUF.5760, file = paste0("BETAPLOTS/anova_",metric,"_5760only.txt"))

jpeg(paste0("BETAPLOTS/NMDS_",metric,"_5760Only.jpeg"), width = 1000, height = 800, pointsize = 24)
par(fig = c(0,0.7,0,1))
plot(NMDS.WUF.5760$points
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.5760$Morph)]
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , sub = paste0("Stress: ", round(NMDS.WUF.5760$stress,2)/100)
     , cex = 1.5
)
lines(NMDS.WUF.5760.CR[NMDS.WUF.5760.CR.chull,]
      , col = "red")
lines(NMDS.WUF.5760.BL[NMDS.WUF.5760.BL.chull,]
      , col = "magenta")
lines(NMDS.WUF.5760.FB[NMDS.WUF.5760.FB.chull,]
      , col = "blue")
par(fig = c(0.7,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , xlab = ""
     , ylab = ""
     , xaxt = "n"
     , yaxt = "n"
     , bty = "n"
     , pch = "")
legend("center"
       , legend = c("Crust","Blade","Finely Br.")
       , pch = 21
       , col= "black"
       , pt.bg = c("red","magenta","blue"))
dev.off()




####### BC #############
metric <- "BC"

### NMDS #####

NMDS.BC.morphonly <- isoMDS(as.matrix(dm.BC.morphonly), y = cmdscale(as.matrix(dm.BC.morphonly), 2))
NMDS.BC.all <- isoMDS(as.matrix(dm.BC.inclWater), y = cmdscale(as.matrix(dm.BC.inclWater), 2))

###### STATS ##########
MF.morphkeep <- MF.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
dm.BC.morphonly.H <- dm.BC.morphonly[grep(".", rownames(dm.BC.morphonly), fixed = TRUE),grep(".", colnames(dm.BC.morphonly), fixed = TRUE)]
dm.BC.morphonly.P <- dm.BC.morphonly[grep("-", rownames(dm.BC.morphonly), fixed = TRUE),grep("-", colnames(dm.BC.morphonly), fixed = TRUE)]

ANOVA.BC.morphtime <- adonis(dm.BC.morphonly ~ Time*Morph, data = MF.morphkeep, by = "margin")
capture.output(ANOVA.BC.morphtime, file = paste0("BETAPLOTS/adonis_", metric,"_Hakai.txt"))
ANOSIM.BC.morphtime <- anosim(dm.BC.morphonly, grouping = MF.morphkeep$Morph)

# Dispersion across time and between morphs
dist.BC.morphonly <- as.dist(dm.BC.morphonly)
betadisp.BC.time <- betadisper(d = dist.BC.morphonly, group = MF.morphkeep$Time)
betadisp.BC.morph <- betadisper(d = dist.BC.morphonly, group = MF.morphkeep$Morph)

capture.output(betadisp.BC.time, file = paste0("BETAPLOTS/betadispTime_", metric, "_Hakai.txt"))
capture.output(betadisp.BC.morph, file = paste0("BETAPLOTS/betadispMorph_", metric, "_Hakai.txt"))


# Grep values for each timepoint and morph, and then plot dispersion within each
TP <- levels(factor(MF.morphkeep$Time))
MorphTypes <- levels(factor(MF.morphkeep$Morph))
# First Timepoint
FirstTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[1],"$"), MF.morphkeep$Time)]
# Between morphs
dm.BC.FirstTP <- dm.BC.morphonly[FirstTPNames,FirstTPNames]
FB.BL.Firstvalues <- as.vector(as.dist(dm.BC.FirstTP[grep("FB", rownames(dm.BC.FirstTP)), grep("BL", colnames(dm.BC.FirstTP))]))
FB.CR.Firstvalues <- as.vector(as.dist(dm.BC.FirstTP[grep("FB", rownames(dm.BC.FirstTP)), grep("CR", colnames(dm.BC.FirstTP))]))
CR.BL.Firstvalues <- as.vector(as.dist(dm.BC.FirstTP[grep("CR", rownames(dm.BC.FirstTP)), grep("BL", colnames(dm.BC.FirstTP))]))
# ANOVA
dm.BC.20 <- dm.BC.morphonly[grep("-20-", rownames(dm.BC.morphonly)), grep("-20-", colnames(dm.BC.morphonly))]
MF.BC.20.only <- MF.morphkeep[grep("-20-", rownames(MF.morphkeep)),]
ANOVA.BC.20.only <- adonis(dm.BC.20 ~ Morph, data = MF.BC.20.only, by = "margin")

# Second Timepoine
SecondTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[2],"$"), MF.morphkeep$Time)]
# Between morphs
dm.BC.SecondTP <- dm.BC.morphonly[SecondTPNames,SecondTPNames]
FB.BL.Secondvalues <- as.vector(as.dist(dm.BC.SecondTP[grep("FB", rownames(dm.BC.SecondTP)), grep("BL", colnames(dm.BC.SecondTP))]))
FB.CR.Secondvalues <- as.vector(as.dist(dm.BC.SecondTP[grep("FB", rownames(dm.BC.SecondTP)), grep("CR", colnames(dm.BC.SecondTP))]))
CR.BL.Secondvalues <- as.vector(as.dist(dm.BC.SecondTP[grep("CR", rownames(dm.BC.SecondTP)), grep("BL", colnames(dm.BC.SecondTP))]))
# ANOVA
dm.BC.60 <- dm.BC.morphonly[grep("-60-", rownames(dm.BC.morphonly)), grep("-60-", colnames(dm.BC.morphonly))]
MF.BC.60.only <- MF.morphkeep[grep("-60-", rownames(MF.morphkeep)),]
ANOVA.BC.60.only <- adonis(dm.BC.60 ~ Morph, data = MF.BC.60.only, by = "margin")

# Third Timepoint
ThirdTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[3],"$"), MF.morphkeep$Time)]
# Between morphs
dm.BC.ThirdTP <- dm.BC.morphonly[ThirdTPNames,ThirdTPNames]
FB.BL.Thirdvalues <- as.vector(as.dist(dm.BC.ThirdTP[grep("FB", rownames(dm.BC.ThirdTP)), grep("BL", colnames(dm.BC.ThirdTP))]))
FB.CR.Thirdvalues <- as.vector(as.dist(dm.BC.ThirdTP[grep("FB", rownames(dm.BC.ThirdTP)), grep("CR", colnames(dm.BC.ThirdTP))]))
CR.BL.Thirdvalues <- as.vector(as.dist(dm.BC.ThirdTP[grep("CR", rownames(dm.BC.ThirdTP)), grep("BL", colnames(dm.BC.ThirdTP))]))
# ANOVA
dm.BC.360 <- dm.BC.morphonly[grep("-360-", rownames(dm.BC.morphonly)), grep("-360-", colnames(dm.BC.morphonly))]
MF.BC.360.only <- MF.morphkeep[grep("-360-", rownames(MF.morphkeep)),]
ANOVA.BC.360.only <- adonis(dm.BC.360 ~ Morph, data = MF.BC.360.only, by = "margin")


# Fourth Timepoint
FourthTPNames <- rownames(MF.morphkeep)[grep(paste0("^",TP[4],"$"), MF.morphkeep$Time)]
# Between morphs
dm.BC.FourthTP <- dm.BC.morphonly[FourthTPNames,FourthTPNames]
FB.BL.Fourthvalues <- as.vector(as.dist(dm.BC.FourthTP[grep("FB", rownames(dm.BC.FourthTP)), grep("BL", colnames(dm.BC.FourthTP))]))
FB.CR.Fourthvalues <- as.vector(as.dist(dm.BC.FourthTP[grep("FB", rownames(dm.BC.FourthTP)), grep("CR", colnames(dm.BC.FourthTP))]))
CR.BL.Fourthvalues <- as.vector(as.dist(dm.BC.FourthTP[grep("CR", rownames(dm.BC.FourthTP)), grep("BL", colnames(dm.BC.FourthTP))]))
# ANOVA
dm.BC.720 <- dm.BC.morphonly[grep("-720-", rownames(dm.BC.morphonly)), grep("-720-", colnames(dm.BC.morphonly))]
MF.BC.720.only <- MF.morphkeep[grep("-720-", rownames(MF.morphkeep)),]
ANOVA.BC.720.only <- adonis(dm.BC.720 ~ Morph, data = MF.BC.720.only, by = "margin")


# Combine into single tables

FB.BL.BC.ALL <- cbind(as.numeric(c(FB.BL.Firstvalues
                                     , FB.BL.Secondvalues
                                     , FB.BL.Thirdvalues
                                     , FB.BL.Fourthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.BL.Firstvalues))
                                       , rep(TP[2], length(FB.BL.Secondvalues))
                                       , rep(TP[3], length(FB.BL.Thirdvalues))
                                       , rep(TP[4], length(FB.BL.Fourthvalues)))))
FB.CR.BC.ALL <- cbind(as.numeric(c(FB.CR.Firstvalues
                                     , FB.CR.Secondvalues
                                     , FB.CR.Thirdvalues
                                     , FB.CR.Fourthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.CR.Firstvalues))
                                       , rep(TP[2], length(FB.CR.Secondvalues))
                                       , rep(TP[3], length(FB.CR.Thirdvalues))
                                       , rep(TP[4], length(FB.CR.Fourthvalues)))))
CR.BL.BC.ALL <- cbind(as.numeric(c(CR.BL.Firstvalues
                                     , CR.BL.Secondvalues
                                     , CR.BL.Thirdvalues
                                     , CR.BL.Fourthvalues))
                        , as.numeric(c(rep(TP[1], length(CR.BL.Firstvalues))
                                       , rep(TP[2], length(CR.BL.Secondvalues))
                                       , rep(TP[3], length(CR.BL.Thirdvalues))
                                       , rep(TP[4], length(CR.BL.Fourthvalues)))))

FB.BL.BC.lm <- lm(FB.BL.BC.ALL[,1] ~ log(FB.BL.BC.ALL[,2]))
summary(FB.BL.BC.lm)

CR.BL.BC.lm <- lm(CR.BL.BC.ALL[,1] ~ log(CR.BL.BC.ALL[,2]))
summary(CR.BL.BC.lm)

FB.CR.BC.lm <- lm(FB.CR.BC.ALL[,1] ~ log(FB.CR.BC.ALL[,2]))
summary(FB.CR.BC.lm)

# Plot the distances

FB.BL.BC.ALL.mean <- aggregate(FB.BL.BC.ALL, by = list(FB.BL.BC.ALL[,2]), mean)
FB.BL.BC.ALL.sd <- aggregate(FB.BL.BC.ALL, by = list(FB.BL.BC.ALL[,2]), sd)

CR.BL.BC.ALL.mean <- aggregate(CR.BL.BC.ALL, by = list(CR.BL.BC.ALL[,2]), mean)
CR.BL.BC.ALL.sd <- aggregate(CR.BL.BC.ALL, by = list(CR.BL.BC.ALL[,2]), sd)

FB.CR.BC.ALL.mean <- aggregate(FB.CR.BC.ALL, by = list(FB.CR.BC.ALL[,2]), mean)
FB.CR.BC.ALL.sd <- aggregate(FB.CR.BC.ALL, by = list(FB.CR.BC.ALL[,2]), sd)

######## PLOT DISP ###############
# ylimits <- c(min(FB.BL.BC.ALL[,1],FB.CR.BC.ALL[,1],CR.BL.BC.ALL[,1]), max(FB.BL.BC.ALL[,1],FB.CR.BC.ALL[,1],CR.BL.BC.ALL[,1]))
ylimits <- c(0.45,0.9)

# xvalues <- log(FB.BL.BC.ALL.mean[,1])
xvalues <- as.character(FB.BL.BC.ALL.mean[,1])
jpeg(paste0("BETAPLOTS/DispOverTime_",metric,".jpeg"))
par(fig = c(0,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1, at = c(1,2,3,4), labels = xvalues)
points(FB.BL.BC.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*0.99
       , x1 = c(1,2,3,4)*0.99
       , y0 = c(FB.BL.BC.ALL.mean[,2] - FB.BL.BC.ALL.sd[,2]/2)
       , y1 = c(FB.BL.BC.ALL.mean[,2] + FB.BL.BC.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.BC.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1
       , x1 = c(1,2,3,4)*1
       , y0 = c(FB.CR.BC.ALL.mean[,2] - FB.CR.BC.ALL.sd[,2]/2)
       , y1 = c(FB.CR.BC.ALL.mean[,2] + FB.CR.BC.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.BC.ALL.mean[,2] ~ c(1,2,3,4)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1.01
       , x1 = c(1,2,3,4)*1.01
       , y0 = c(CR.BL.BC.ALL.mean[,2] - CR.BL.BC.ALL.sd[,2]/2)
       , y1 = c(CR.BL.BC.ALL.mean[,2] + CR.BL.BC.ALL.sd[,2]/2)
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
       , col = c("darkgreen","purple","grey"))
dev.off()

####### PLOT MORPH ############
# NOCHLP BC
MF.morphkeep$Morph <- factor(MF.morphkeep$Morph, levels = c('CR','BL','FB'))
MF.morphkeep$Time <- factor(MF.morphkeep$Time, levels = c('20','60','180','360','720','1440'))
MF.morphkeep$Type <- factor(MF.morphkeep$Type, levels = c('P','H'))
MorphColours <- c("red","magenta","blue")


# MAKE POLYGONS for plotting
NMDS.BC.CR <- NMDS.BC.morphonly$points[grep("CR", rownames(NMDS.BC.morphonly$points)),]
NMDS.BC.BL <- NMDS.BC.morphonly$points[grep("BL", rownames(NMDS.BC.morphonly$points)),]
NMDS.BC.FB <- NMDS.BC.morphonly$points[grep("FB", rownames(NMDS.BC.morphonly$points)),]

NMDS.BC.CR.chull <- chull(NMDS.BC.CR)
NMDS.BC.CR.chull <- c(NMDS.BC.CR.chull, NMDS.BC.CR.chull[1])

NMDS.BC.BL.chull <- chull(NMDS.BC.BL)
NMDS.BC.BL.chull <- c(NMDS.BC.BL.chull, NMDS.BC.BL.chull[1])

NMDS.BC.FB.chull <- chull(NMDS.BC.FB)
NMDS.BC.FB.chull <- c(NMDS.BC.FB.chull, NMDS.BC.FB.chull[1])


jpeg(paste0("BETAPLOTS/NMDS_",metric,"_Morph.jpeg"))
par(fig = c(0,0.8,0,1))
plot(NMDS.BC.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.morphkeep$Morph)]
     , sub = round(NMDS.BC.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
)
lines(NMDS.BC.CR[NMDS.BC.CR.chull,]
      , col = "red")
lines(NMDS.BC.BL[NMDS.BC.BL.chull,]
      , col = "magenta")
lines(NMDS.BC.FB[NMDS.BC.FB.chull,]
      , col = "blue")
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
)
dev.off()

####### PLOT TIME ############
TimeColours <- c("grey","lightblue","blue","darkblue")

jpeg(paste0("BETAPLOTS/NMDS_",metric,"_Time.jpeg"))
par(fig = c(0,0.8,0,1))
plot(NMDS.BC.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = TimeColours[factor(MF.morphkeep$Time)]
     , sub = round(NMDS.BC.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
)
# lines(NMDS.BC.CR[NMDS.BC.CR.chull,]
#       , col = "red")
# lines(NMDS.BC.BL[NMDS.BC.BL.chull,]
#       , col = "magenta")
# lines(NMDS.BC.FB[NMDS.BC.FB.chull,]
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



# TIME-- ONE at a time


# 20 minutes
NMDS.BC.20.only <- NMDS.BC.morphonly$points[grep("-20-", rownames(NMDS.BC.morphonly$points)),]

NMDS.BC.20.only.CR <- NMDS.BC.20.only[grep("CR", rownames(NMDS.BC.20.only)),]
NMDS.BC.20.only.CR.chull <- chull(NMDS.BC.20.only.CR)
NMDS.BC.20.only.CR.chull <- c(NMDS.BC.20.only.CR.chull, NMDS.BC.20.only.CR.chull[1])

NMDS.BC.20.only.BL <- NMDS.BC.20.only[grep("BL", rownames(NMDS.BC.20.only)),]
NMDS.BC.20.only.BL.chull <- chull(NMDS.BC.20.only.BL)
NMDS.BC.20.only.BL.chull <- c(NMDS.BC.20.only.BL.chull, NMDS.BC.20.only.BL.chull[1])

NMDS.BC.20.only.FB <- NMDS.BC.20.only[grep("FB", rownames(NMDS.BC.20.only)),]
NMDS.BC.20.only.FB.chull <- chull(NMDS.BC.20.only.FB)
NMDS.BC.20.only.FB.chull <- c(NMDS.BC.20.only.FB.chull, NMDS.BC.20.only.FB.chull[1])

# 60 minutes
NMDS.BC.60.only <- NMDS.BC.morphonly$points[grep("-60-", rownames(NMDS.BC.morphonly$points)),]


NMDS.BC.60.only.CR <- NMDS.BC.60.only[grep("CR", rownames(NMDS.BC.60.only)),]
NMDS.BC.60.only.CR.chull <- chull(NMDS.BC.60.only.CR)
NMDS.BC.60.only.CR.chull <- c(NMDS.BC.60.only.CR.chull, NMDS.BC.60.only.CR.chull[1])

NMDS.BC.60.only.BL <- NMDS.BC.60.only[grep("BL", rownames(NMDS.BC.60.only)),]
NMDS.BC.60.only.BL.chull <- chull(NMDS.BC.60.only.BL)
NMDS.BC.60.only.BL.chull <- c(NMDS.BC.60.only.BL.chull, NMDS.BC.60.only.BL.chull[1])

NMDS.BC.60.only.FB <- NMDS.BC.60.only[grep("FB", rownames(NMDS.BC.60.only)),]
NMDS.BC.60.only.FB.chull <- chull(NMDS.BC.60.only.FB)
NMDS.BC.60.only.FB.chull <- c(NMDS.BC.60.only.FB.chull, NMDS.BC.60.only.FB.chull[1])



# 360 minutes
NMDS.BC.360.only <- NMDS.BC.morphonly$points[grep("-360-", rownames(NMDS.BC.morphonly$points)),]


NMDS.BC.360.only.CR <- NMDS.BC.360.only[grep("CR", rownames(NMDS.BC.360.only)),]
NMDS.BC.360.only.CR.chull <- chull(NMDS.BC.360.only.CR)
NMDS.BC.360.only.CR.chull <- c(NMDS.BC.360.only.CR.chull, NMDS.BC.360.only.CR.chull[1])

NMDS.BC.360.only.BL <- NMDS.BC.360.only[grep("BL", rownames(NMDS.BC.360.only)),]
NMDS.BC.360.only.BL.chull <- chull(NMDS.BC.360.only.BL)
NMDS.BC.360.only.BL.chull <- c(NMDS.BC.360.only.BL.chull, NMDS.BC.360.only.BL.chull[1])

NMDS.BC.360.only.FB <- NMDS.BC.360.only[grep("FB", rownames(NMDS.BC.360.only)),]
NMDS.BC.360.only.FB.chull <- chull(NMDS.BC.360.only.FB)
NMDS.BC.360.only.FB.chull <- c(NMDS.BC.360.only.FB.chull, NMDS.BC.360.only.FB.chull[1])


# 720 minutes
NMDS.BC.720.only <- NMDS.BC.morphonly$points[grep("-720-", rownames(NMDS.BC.morphonly$points)),]


NMDS.BC.720.only.CR <- NMDS.BC.720.only[grep("CR", rownames(NMDS.BC.720.only)),]
NMDS.BC.720.only.CR.chull <- chull(NMDS.BC.720.only.CR)
NMDS.BC.720.only.CR.chull <- c(NMDS.BC.720.only.CR.chull, NMDS.BC.720.only.CR.chull[1])

NMDS.BC.720.only.BL <- NMDS.BC.720.only[grep("BL", rownames(NMDS.BC.720.only)),]
NMDS.BC.720.only.BL.chull <- chull(NMDS.BC.720.only.BL)
NMDS.BC.720.only.BL.chull <- c(NMDS.BC.720.only.BL.chull, NMDS.BC.720.only.BL.chull[1])

NMDS.BC.720.only.FB <- NMDS.BC.720.only[grep("FB", rownames(NMDS.BC.720.only)),]
NMDS.BC.720.only.FB.chull <- chull(NMDS.BC.720.only.FB)
NMDS.BC.720.only.FB.chull <- c(NMDS.BC.720.only.FB.chull, NMDS.BC.720.only.FB.chull[1])


jpeg(paste0("BETAPLOTS/NMDS_",metric,"_TimebyMorph.jpeg"), width = 2000, height = 800, pointsize = 24)
par(mfrow= c(1,4), oma = c(4,4,4,6))
par(mar = c(4,0,4,0))
plot(NMDS.BC.20.only
     , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.BC.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.20.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.20.only.CR[NMDS.BC.20.only.CR.chull,]
      , col = "red")
lines(NMDS.BC.20.only.BL[NMDS.BC.20.only.BL.chull,]
      , col = "magenta")
lines(NMDS.BC.20.only.FB[NMDS.BC.20.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.BC.60.only
     , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.BC.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.BC.60.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.60.only.CR[NMDS.BC.60.only.CR.chull,]
      , col = "red")
lines(NMDS.BC.60.only.BL[NMDS.BC.60.only.BL.chull,]
      , col = "magenta")
lines(NMDS.BC.60.only.FB[NMDS.BC.60.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.BC.360.only
     , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.BC.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.360.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.360.only.CR[NMDS.BC.360.only.CR.chull,]
      , col = "red")
lines(NMDS.BC.360.only.BL[NMDS.BC.360.only.BL.chull,]
      , col = "magenta")
lines(NMDS.BC.360.only.FB[NMDS.BC.360.only.FB.chull,]
      , col = "blue")
par(mar = c(4,0,4,0))
plot(NMDS.BC.720.only
     , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.BC.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.720.only.CR[NMDS.BC.720.only.CR.chull,]
      , col = "red")
lines(NMDS.BC.720.only.BL[NMDS.BC.720.only.BL.chull,]
      , col = "magenta")
lines(NMDS.BC.720.only.FB[NMDS.BC.720.only.FB.chull,]
      , col = "blue")
#STOP
dev.off()


############ PLOT 5760 ################
dm.BC.5760<- dm.BC.inclWater[grep("(CR|FB|BL)-5760", rownames(dm.BC.inclWater)),grep("(CR|FB|BL)-5760", colnames(dm.BC.inclWater))]

NMDS.BC.5760 <- isoMDS(as.matrix(dm.BC.5760), y = cmdscale(as.matrix(dm.BC.5760), 2))

# NMDS.BC.morphallTimePoints <- NMDS.BC.morphAllTime$points[grep("(CR|BL|FB)-5760", rownames(NMDS.BC.morphAllTime$points)),]
MF.5760 <- MF[sapply(rownames(NMDS.BC.5760$points), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
}),]
MF.5760$Morph <- factor(MF.5760$Morph, levels = c("CR","BL","FB"))
MorphColours <- c("red","magenta","blue")

NMDS.BC.5760.CR <- NMDS.BC.5760$points[grep("CR", rownames(NMDS.BC.5760$points)),]
NMDS.BC.5760.CR.chull <- chull(NMDS.BC.5760.CR)
NMDS.BC.5760.CR.chull <- c(NMDS.BC.5760.CR.chull, NMDS.BC.5760.CR.chull[1])

NMDS.BC.5760.BL <- NMDS.BC.5760$points[grep("BL", rownames(NMDS.BC.5760$points)),]
NMDS.BC.5760.BL.chull <- chull(NMDS.BC.5760.BL)
NMDS.BC.5760.BL.chull <- c(NMDS.BC.5760.BL.chull, NMDS.BC.5760.BL.chull[1])

NMDS.BC.5760.FB <- NMDS.BC.5760$points[grep("FB", rownames(NMDS.BC.5760$points)),]
NMDS.BC.5760.FB.chull <- chull(NMDS.BC.5760.FB)
NMDS.BC.5760.FB.chull <- c(NMDS.BC.5760.FB.chull, NMDS.BC.5760.FB.chull[1])

ANOVA.BC.5760 <- adonis(dm.BC.5760 ~ Morph, data = MF.5760)
ANOSIM.BC.5760 <- anosim(dat = dm.BC.5760, grouping = MF.5760$Morph)
capture.output(ANOVA.BC.5760, file = paste0("BETAPLOTS/anova_",metric,"_5760only.txt"))

jpeg(paste0("BETAPLOTS/NMDS_",metric,"_5760Only.jpeg"), width = 1000, height = 800, pointsize = 24)
par(fig = c(0,0.7,0,1))
plot(NMDS.BC.5760$points
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.5760$Morph)]
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , sub = paste0("Stress: ", round(NMDS.BC.5760$stress,2)/100)
     , cex = 1.5
)
lines(NMDS.BC.5760.CR[NMDS.BC.5760.CR.chull,]
      , col = "red")
lines(NMDS.BC.5760.BL[NMDS.BC.5760.BL.chull,]
      , col = "magenta")
lines(NMDS.BC.5760.FB[NMDS.BC.5760.FB.chull,]
      , col = "blue")
par(fig = c(0.7,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , xlab = ""
     , ylab = ""
     , xaxt = "n"
     , yaxt = "n"
     , bty = "n"
     , pch = "")
legend("center"
       , legend = c("Crust","Blade","Finely Br.")
       , pch = 21
       , col= "black"
       , pt.bg = c("red","magenta","blue"))
dev.off()


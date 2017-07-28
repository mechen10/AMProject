#!/bin/RScript
set.seed(3)
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


#*** NOTE: ONLY WORKING ON UWUF RIGHT NOW; DO OTHERS AFTER
####### *********UWUF********* ############# 
metric <- "UWUF"

######## --HAKAI-- ###########
### NMDS #####

NMDS.UWUF.morphonly <- isoMDS(as.matrix(dm.UWUF.morphonly), y = cmdscale(as.matrix(dm.UWUF.morphonly), 2))
NMDS.UWUF.all <- isoMDS(as.matrix(dm.UWUF.inclWater), y = cmdscale(as.matrix(dm.UWUF.inclWater), 2))

###### STATS ##########
MF.morphkeep <- MF.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
MF.morphkeep$Morph <- factor(MF.morphkeep$Morph, levels = c('CR','BL','FB'))
MF.morphkeep$Time <- factor(MF.morphkeep$Time, levels = c('20','60','180','360','720','1440','5760'))
MF.morphkeep$Type <- factor(MF.morphkeep$Type, levels = c('P','H'))


# dm.UWUF.morphonly.H <- dm.UWUF.morphonly[grep(".", rownames(dm.UWUF.morphonly), fixed = TRUE),grep(".", colnames(dm.UWUF.morphonly), fixed = TRUE)]
# dm.UWUF.morphonly.P <- dm.UWUF.morphonly[grep("-", rownames(dm.UWUF.morphonly), fixed = TRUE),grep("-", colnames(dm.UWUF.morphonly), fixed = TRUE)]

ANOVA.UWUF.morphtime <- adonis(dm.UWUF.morphonly ~ Time*Morph, data = MF.morphkeep, by = "margin")
capture.output(ANOVA.UWUF.morphtime, file = paste0("BETAPLOTS_H/adonis_", metric,"_Hakai.txt"))
# ANOSIM.UWUF.morphtime <- anosim(dm.UWUF.morphonly, grouping = MF.morphkeep$Morph)

# Dispersion across time and between morphs
dist.UWUF.morphonly <- as.dist(dm.UWUF.inclWater[-grep("W", rownames(dm.UWUF.inclWater)), -grep("W", colnames(dm.UWUF.inclWater))])
MF.incl5760 <- MF.inclWater[-grep("W", rownames(MF.inclWater)),]
betadisp.UWUF.time <- betadisper(d = dist.UWUF.morphonly, group = MF.incl5760$Time)
betadisp.UWUF.morph <- betadisper(d = dist.UWUF.morphonly, group = MF.incl5760$Morph)

capture.output(betadisp.UWUF.time, file = paste0("BETAPLOTS_H/betadispTime_", metric, "_Hakai.txt"))
capture.output(betadisp.UWUF.morph, file = paste0("BETAPLOTS_H/betadispMorph_", metric, "_Hakai.txt"))


# Grep values for each timepoint and morph, and then plot dispersion within each
TP <- levels(factor(MF.inclWater$Time))
MorphTypes <- levels(factor(MF.morphkeep$Morph))
# First Timepoint
FirstTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[1],"$"), MF.inclWater$Time)]
toremove <- grep("W", FirstTPNames)
if (length(toremove) > 0) {
  FirstTPNames <- FirstTPNames[-toremove]
}
# Between morphs
dm.UWUF.FirstTP <- dm.UWUF.inclWater[FirstTPNames,FirstTPNames]
FB.BL.Firstvalues <- as.vector(as.dist(dm.UWUF.FirstTP[grep("FB", rownames(dm.UWUF.FirstTP)), grep("BL", colnames(dm.UWUF.FirstTP))]))
FB.CR.Firstvalues <- as.vector(as.dist(dm.UWUF.FirstTP[grep("FB", rownames(dm.UWUF.FirstTP)), grep("CR", colnames(dm.UWUF.FirstTP))]))
CR.BL.Firstvalues <- as.vector(as.dist(dm.UWUF.FirstTP[grep("CR", rownames(dm.UWUF.FirstTP)), grep("BL", colnames(dm.UWUF.FirstTP))]))

# ANOVA
dm.UWUF.20 <- dm.UWUF.inclWater[grep("(CR|BL|FB)-20-", rownames(dm.UWUF.inclWater)), grep("-20-", colnames(dm.UWUF.inclWater))]
MF.UWUF.20.only <- MF.inclWater[grep("(CR|BL|FB)-20-", rownames(MF.inclWater)),]
ANOVA.UWUF.20.only <- adonis(dm.UWUF.20 ~ Morph, data = MF.UWUF.20.only, by = "margin")

# Second Timepoine
SecondTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[2],"$"), MF.inclWater$Time)]
toremove <- grep("W", SecondTPNames)
if (length(toremove) > 0) {
  SecondTPNames <- SecondTPNames[-toremove]
}
# Between morphs
dm.UWUF.SecondTP <- dm.UWUF.inclWater[SecondTPNames,SecondTPNames]
FB.BL.Secondvalues <- as.vector(as.dist(dm.UWUF.SecondTP[grep("FB", rownames(dm.UWUF.SecondTP)), grep("BL", colnames(dm.UWUF.SecondTP))]))
FB.CR.Secondvalues <- as.vector(as.dist(dm.UWUF.SecondTP[grep("FB", rownames(dm.UWUF.SecondTP)), grep("CR", colnames(dm.UWUF.SecondTP))]))
CR.BL.Secondvalues <- as.vector(as.dist(dm.UWUF.SecondTP[grep("CR", rownames(dm.UWUF.SecondTP)), grep("BL", colnames(dm.UWUF.SecondTP))]))
# ANOVA
dm.UWUF.60 <- dm.UWUF.inclWater[grep("(BL|CR|FB)-60-", rownames(dm.UWUF.inclWater)), grep("-60-", colnames(dm.UWUF.morphonly))]
MF.UWUF.60.only <- MF.inclWater[grep("(BL|CR|FB)-60-", rownames(MF.inclWater)),]
ANOVA.UWUF.60.only <- adonis(dm.UWUF.60 ~ Morph, data = MF.UWUF.60.only, by = "margin")

# Third Timepoint
ThirdTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[3],"$"), MF.inclWater$Time)]
toremove <- grep("W", ThirdTPNames)
if (length(toremove) > 0) {
  ThirdTPNames <- ThirdTPNames[-toremove]
}
# Between morphs
dm.UWUF.ThirdTP <- dm.UWUF.inclWater[ThirdTPNames,ThirdTPNames]
FB.BL.Thirdvalues <- as.vector(as.dist(dm.UWUF.ThirdTP[grep("FB", rownames(dm.UWUF.ThirdTP)), grep("BL", colnames(dm.UWUF.ThirdTP))]))
FB.CR.Thirdvalues <- as.vector(as.dist(dm.UWUF.ThirdTP[grep("FB", rownames(dm.UWUF.ThirdTP)), grep("CR", colnames(dm.UWUF.ThirdTP))]))
CR.BL.Thirdvalues <- as.vector(as.dist(dm.UWUF.ThirdTP[grep("CR", rownames(dm.UWUF.ThirdTP)), grep("BL", colnames(dm.UWUF.ThirdTP))]))
# ANOVA
dm.UWUF.360 <- dm.UWUF.inclWater[grep("(BL|CR|FB)-360-", rownames(dm.UWUF.inclWater)), grep("-360-", colnames(dm.UWUF.morphonly))]
MF.UWUF.360.only <- MF.inclWater[grep("(BL|CR|FB)-360-", rownames(MF.inclWater)),]
ANOVA.UWUF.360.only <- adonis(dm.UWUF.360 ~ Morph, data = MF.UWUF.360.only, by = "margin")


# Fourth Timepoint
FourthTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[4],"$"), MF.inclWater$Time)]
toremove <- grep("W", FourthTPNames)
if (length(toremove) > 0) {
  FourthTPNames <- FourthTPNames[-toremove]
}
# Between morphs
dm.UWUF.FourthTP <- dm.UWUF.inclWater[FourthTPNames,FourthTPNames]
FB.BL.Fourthvalues <- as.vector(as.dist(dm.UWUF.FourthTP[grep("FB", rownames(dm.UWUF.FourthTP)), grep("BL", colnames(dm.UWUF.FourthTP))]))
FB.CR.Fourthvalues <- as.vector(as.dist(dm.UWUF.FourthTP[grep("FB", rownames(dm.UWUF.FourthTP)), grep("CR", colnames(dm.UWUF.FourthTP))]))
CR.BL.Fourthvalues <- as.vector(as.dist(dm.UWUF.FourthTP[grep("CR", rownames(dm.UWUF.FourthTP)), grep("BL", colnames(dm.UWUF.FourthTP))]))
# ANOVA
dm.UWUF.720 <- dm.UWUF.inclWater[grep("(BL|CR|FB)-720-", rownames(dm.UWUF.inclWater)), grep("-720-", colnames(dm.UWUF.morphonly))]
MF.UWUF.720.only <- MF.inclWater[grep("(BL|CR|FB)-720-", rownames(MF.inclWater)),]
ANOVA.UWUF.720.only <- adonis(dm.UWUF.720 ~ Morph, data = MF.UWUF.720.only, by = "margin")

# Fifth Timepoint
FifthTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[5],"$"), MF.inclWater$Time)]
toremove <- grep("W", FifthTPNames)
if (length(toremove) > 0) {
  FifthTPNames <- FifthTPNames[-toremove]
}
# Between morphs
dm.UWUF.FifthTP <- dm.UWUF.inclWater[FifthTPNames,FifthTPNames]
FB.BL.Fifthvalues <- as.vector(as.dist(dm.UWUF.FifthTP[grep("FB", rownames(dm.UWUF.FifthTP)), grep("BL", colnames(dm.UWUF.FifthTP))]))
FB.CR.Fifthvalues <- as.vector(as.dist(dm.UWUF.FifthTP[grep("FB", rownames(dm.UWUF.FifthTP)), grep("CR", colnames(dm.UWUF.FifthTP))]))
CR.BL.Fifthvalues <- as.vector(as.dist(dm.UWUF.FifthTP[grep("CR", rownames(dm.UWUF.FifthTP)), grep("BL", colnames(dm.UWUF.FifthTP))]))
# ANOVA
dm.UWUF.5760 <- dm.UWUF.inclWater[grep("(BL|CR|FB)-5760-", rownames(dm.UWUF.inclWater)), grep("-5760-", colnames(dm.UWUF.inclWater))]
MF.UWUF.5760.only <- MF.inclWater[grep("(BL|CR|FB)-5760-", rownames(MF.inclWater)),]
ANOVA.UWUF.5760.only <- adonis(dm.UWUF.5760 ~ Morph, data = MF.UWUF.5760.only, by = "margin")


# Combine into single tables

FB.BL.UWUF.ALL <- cbind(as.numeric(c(FB.BL.Firstvalues
                                     , FB.BL.Secondvalues
                                     , FB.BL.Thirdvalues
                                     , FB.BL.Fourthvalues
                                     , FB.BL.Fifthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.BL.Firstvalues))
                                       , rep(TP[2], length(FB.BL.Secondvalues))
                                       , rep(TP[3], length(FB.BL.Thirdvalues))
                                       , rep(TP[4], length(FB.BL.Fourthvalues))
                                       , rep(TP[5], length(FB.BL.Fifthvalues)))))
FB.CR.UWUF.ALL <- cbind(as.numeric(c(FB.CR.Firstvalues
                                     , FB.CR.Secondvalues
                                     , FB.CR.Thirdvalues
                                     , FB.CR.Fourthvalues
                                     , FB.CR.Fifthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.CR.Firstvalues))
                                       , rep(TP[2], length(FB.CR.Secondvalues))
                                       , rep(TP[3], length(FB.CR.Thirdvalues))
                                       , rep(TP[4], length(FB.CR.Fourthvalues))
                                       , rep(TP[5], length(FB.CR.Fifthvalues)))))
CR.BL.UWUF.ALL <- cbind(as.numeric(c(CR.BL.Firstvalues
                                     , CR.BL.Secondvalues
                                     , CR.BL.Thirdvalues
                                     , CR.BL.Fourthvalues
                                     , CR.BL.Fifthvalues))
                        , as.numeric(c(rep(TP[1], length(CR.BL.Firstvalues))
                                       , rep(TP[2], length(CR.BL.Secondvalues))
                                       , rep(TP[3], length(CR.BL.Thirdvalues))
                                       , rep(TP[4], length(CR.BL.Fourthvalues))
                                       , rep(TP[5], length(CR.BL.Fifthvalues)))))

# FB.BL.UWUF.lm <- lm(FB.BL.UWUF.ALL[,1] ~ log(FB.BL.UWUF.ALL[,2]))
# summary(FB.BL.UWUF.lm)
# 
# CR.BL.UWUF.lm <- lm(CR.BL.UWUF.ALL[,1] ~ log(CR.BL.UWUF.ALL[,2]))
# summary(CR.BL.UWUF.lm)
# 
# FB.CR.UWUF.lm <- lm(FB.CR.UWUF.ALL[,1] ~ log(FB.CR.UWUF.ALL[,2]))
# summary(FB.CR.UWUF.lm)

# Plot the distances

FB.BL.UWUF.ALL.mean <- aggregate(FB.BL.UWUF.ALL, by = list(FB.BL.UWUF.ALL[,2]), mean)
FB.BL.UWUF.ALL.sd <- aggregate(FB.BL.UWUF.ALL, by = list(FB.BL.UWUF.ALL[,2]), sd)

CR.BL.UWUF.ALL.mean <- aggregate(CR.BL.UWUF.ALL, by = list(CR.BL.UWUF.ALL[,2]), mean)
CR.BL.UWUF.ALL.sd <- aggregate(CR.BL.UWUF.ALL, by = list(CR.BL.UWUF.ALL[,2]), sd)

FB.CR.UWUF.ALL.mean <- aggregate(FB.CR.UWUF.ALL, by = list(FB.CR.UWUF.ALL[,2]), mean)
FB.CR.UWUF.ALL.sd <- aggregate(FB.CR.UWUF.ALL, by = list(FB.CR.UWUF.ALL[,2]), sd)


######## PLOT DISP ###############
# ylimits <- c(min(FB.BL.UWUF.ALL[,1],FB.CR.UWUF.ALL[,1],CR.BL.UWUF.ALL[,1]), max(FB.BL.UWUF.ALL[,1],FB.CR.UWUF.ALL[,1],CR.BL.UWUF.ALL[,1]))
ylimits <- c(0.4,0.7)

# xvalues <- log(FB.BL.UWUF.ALL.mean[,1])
xvalues <- as.character(FB.BL.UWUF.ALL.mean[,1])
pdf(paste0("BETAPLOTS_H/DispOverTime_H_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1, at = c(1,2,3,4,5), labels = c('20 min','1 h', '6 h', '12 h','4 d'))
points(FB.BL.UWUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*0.99
       , x1 = c(1,2,3,4,5)*0.99
       , y0 = c(FB.BL.UWUF.ALL.mean[,2] - FB.BL.UWUF.ALL.sd[,2]/2)
       , y1 = c(FB.BL.UWUF.ALL.mean[,2] + FB.BL.UWUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.UWUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1
       , x1 = c(1,2,3,4,5)*1
       , y0 = c(FB.CR.UWUF.ALL.mean[,2] - FB.CR.UWUF.ALL.sd[,2]/2)
       , y1 = c(FB.CR.UWUF.ALL.mean[,2] + FB.CR.UWUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.UWUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1.01
       , x1 = c(1,2,3,4,5)*1.01
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
       , col = c("darkgreen","purple","grey")
       , lwd = 2)
dev.off()

####### PLOT MORPH ############
# NO 5760 UWUF

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


MorphColours <- c("darkorchid4","dodgerblue","salmon") 

pdf(paste0("BETAPLOTS_H/NMDS_H_",metric,"_Morph.pdf"), pointsize = 14)
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
      , col = MorphColours[1])
lines(NMDS.UWUF.BL[NMDS.UWUF.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.FB[NMDS.UWUF.FB.chull,]
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

# 5760 minutes
NMDS.UWUF.5760.only <- NMDS.UWUF.all$points[grep("(BL|CR|FB)-5760-", rownames(NMDS.UWUF.all$points)),]

NMDS.UWUF.5760.only.CR <- NMDS.UWUF.5760.only[grep("CR", rownames(NMDS.UWUF.5760.only)),]
NMDS.UWUF.5760.only.CR.chull <- chull(NMDS.UWUF.5760.only.CR)
NMDS.UWUF.5760.only.CR.chull <- c(NMDS.UWUF.5760.only.CR.chull, NMDS.UWUF.5760.only.CR.chull[1])

NMDS.UWUF.5760.only.BL <- NMDS.UWUF.5760.only[grep("BL", rownames(NMDS.UWUF.5760.only)),]
NMDS.UWUF.5760.only.BL.chull <- chull(NMDS.UWUF.5760.only.BL)
NMDS.UWUF.5760.only.BL.chull <- c(NMDS.UWUF.5760.only.BL.chull, NMDS.UWUF.5760.only.BL.chull[1])

NMDS.UWUF.5760.only.FB <- NMDS.UWUF.5760.only[grep("FB", rownames(NMDS.UWUF.5760.only)),]
NMDS.UWUF.5760.only.FB.chull <- chull(NMDS.UWUF.5760.only.FB)
NMDS.UWUF.5760.only.FB.chull <- c(NMDS.UWUF.5760.only.FB.chull, NMDS.UWUF.5760.only.FB.chull[1])


# pdf(paste0("BETAPLOTS_H/NMDS_",metric,"_TimebyMorph.pdf"),pointsize = 14, width = 7, height = 4)
# par(mfrow= c(1,4), oma = c(4,4,4,6))
# par(mar = c(4,0,4,0))
# plot(NMDS.UWUF.20.only
#      , main = "20 Minutes"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.UWUF.20.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p = ", ANOVA.UWUF.20.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.UWUF.20.only.CR[NMDS.UWUF.20.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.UWUF.20.only.BL[NMDS.UWUF.20.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.UWUF.20.only.FB[NMDS.UWUF.20.only.FB.chull,]
#       , col = MorphColours[3])
# par(mar = c(4,0,4,0))
# plot(NMDS.UWUF.60.only
#      , main = "1 Hour"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.UWUF.60.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p =  ",ANOVA.UWUF.60.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.UWUF.60.only.CR[NMDS.UWUF.60.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.UWUF.60.only.BL[NMDS.UWUF.60.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.UWUF.60.only.FB[NMDS.UWUF.60.only.FB.chull,]
#       , col = MorphColours[3])
# par(mar = c(4,0,4,0))
# plot(NMDS.UWUF.360.only
#      , main = "6 Hours"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.UWUF.360.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p = ", ANOVA.UWUF.360.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.UWUF.360.only.CR[NMDS.UWUF.360.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.UWUF.360.only.BL[NMDS.UWUF.360.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.UWUF.360.only.FB[NMDS.UWUF.360.only.FB.chull,]
#       , col = MorphColours[3])
# par(mar = c(4,0,4,0))
# plot(NMDS.UWUF.720.only
#      , main = "12 Hours"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.UWUF.720.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p = ", ANOVA.UWUF.720.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.UWUF.720.only.CR[NMDS.UWUF.720.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.UWUF.720.only.BL[NMDS.UWUF.720.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.UWUF.720.only.FB[NMDS.UWUF.720.only.FB.chull,]
#       , col = MorphColours[3])
# #STOP
# dev.off()

############ COMBO DISP BETA ################
# Disp and beta through time combined
# xvalues <- log(FB.BL.UWUF.ALL.mean[,1])
xvalues <- as.character(FB.BL.UWUF.ALL.mean[,1])

# Change factor of individual plots
MF.UWUF.20.only$Morph <- factor(MF.UWUF.20.only$Morph, levels = c("CR","BL","FB"))
MF.UWUF.60.only$Morph <- factor(MF.UWUF.60.only$Morph, levels = c("CR","BL","FB"))
MF.UWUF.360.only$Morph <- factor(MF.UWUF.360.only$Morph, levels = c("CR","BL","FB"))
MF.UWUF.720.only$Morph <- factor(MF.UWUF.720.only$Morph, levels = c("CR","BL","FB"))
MF.UWUF.5760.only$Morph <- factor(MF.UWUF.5760.only$Morph, levels = c("CR","BL","FB"))


pdf(paste0("BETAPLOTS_H/COMBO_H_dispbeta_",metric,".pdf"), pointsize = 14, width = 10, height = 7)
par(fig = c(0,0.8,0.3,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = ''
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1
     , at = c(1,2,3,4,5)
     , labels = c("20 min","1 h","6 h","12 h","4 d")
     , las = 2)
points(FB.BL.UWUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*0.99
       , x1 = c(1,2,3,4,5)*0.99
       , y0 = c(FB.BL.UWUF.ALL.mean[,2] - FB.BL.UWUF.ALL.sd[,2])
       , y1 = c(FB.BL.UWUF.ALL.mean[,2] + FB.BL.UWUF.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.UWUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1
       , x1 = c(1,2,3,4,5)*1
       , y0 = c(FB.CR.UWUF.ALL.mean[,2] - FB.CR.UWUF.ALL.sd[,2])
       , y1 = c(FB.CR.UWUF.ALL.mean[,2] + FB.CR.UWUF.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.UWUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1.01
       , x1 = c(1,2,3,4,5)*1.01
       , y0 = c(CR.BL.UWUF.ALL.mean[,2] - CR.BL.UWUF.ALL.sd[,2])
       , y1 = c(CR.BL.UWUF.ALL.mean[,2] + CR.BL.UWUF.ALL.sd[,2])
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
       , legend = c("FB:CR", "FB:BL", "CR:BL")
       , lty = 1
       , col = c("darkgreen","purple","grey")
       , lwd = 2 
       )
# EACH GETS 0.15 SPACE TOTAL;
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.049,0.1932,0,0.4), new = TRUE)
plot(NMDS.UWUF.20.only
     # , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.20.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.20.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.20.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.20.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)

lines(NMDS.UWUF.20.only.CR[NMDS.UWUF.20.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.20.only.BL[NMDS.UWUF.20.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.20.only.FB[NMDS.UWUF.20.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.2032,0.3474,0,0.4), new = TRUE)
plot(NMDS.UWUF.60.only
     # , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.60.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.60.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.60.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.60.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.60.only.CR[NMDS.UWUF.60.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.60.only.BL[NMDS.UWUF.60.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.60.only.FB[NMDS.UWUF.60.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.3574,0.5016,0,0.4), new = TRUE)
plot(NMDS.UWUF.360.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.360.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.360.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.360.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.360.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.360.only.CR[NMDS.UWUF.360.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.360.only.BL[NMDS.UWUF.360.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.360.only.FB[NMDS.UWUF.360.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.5116,0.6558,0,0.4), new = TRUE)
plot(NMDS.UWUF.720.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.720.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.720.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.720.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.720.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.720.only.CR[NMDS.UWUF.720.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.720.only.BL[NMDS.UWUF.720.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.720.only.FB[NMDS.UWUF.720.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.6658,0.81,0,0.4), new = TRUE)
plot(NMDS.UWUF.5760.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.UWUF.5760.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.5760.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.5760.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.5760.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.5760.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.5760.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.5760.only.CR[NMDS.UWUF.5760.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.5760.only.BL[NMDS.UWUF.5760.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.5760.only.FB[NMDS.UWUF.5760.only.FB.chull,]
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

betadisp.H.UWUF.FB <- betadisp.UWUF.time$distances[grep("FB", names(betadisp.UWUF.time$distances))]
betadisp.H.UWUF.BL <- betadisp.UWUF.time$distances[grep("BL", names(betadisp.UWUF.time$distances))]
betadisp.H.UWUF.CR <- betadisp.UWUF.time$distances[grep("CR", names(betadisp.UWUF.time$distances))]

MF.H.FB <- MF.incl5760[grep("FB", rownames(MF.incl5760)),]
MF.H.BL <- MF.incl5760[grep("BL", rownames(MF.incl5760)),]
MF.H.CR <- MF.incl5760[grep("CR", rownames(MF.incl5760)),]

betadisp.FB.H.agg <- aggregate(betadisp.H.UWUF.FB, by = list(MF.H.FB$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.BL.H.agg <- aggregate(betadisp.H.UWUF.BL, by = list(MF.H.BL$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.CR.H.agg <- aggregate(betadisp.H.UWUF.CR, by = list(MF.H.CR$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )

betadisp.UWUF.time.forstat <- cbind(betadisp.UWUF.time$distances, MF.incl5760[unlist(lapply(names(betadisp.UWUF.time$distances), function(x) {grep(paste0("^",x,"$"), rownames(MF.incl5760))})),])
colnames(betadisp.UWUF.time.forstat)[1] <- c("Distance")
ANOVA.betadisp.UWUF <- anova(lm(Distance ~ Time*Morph, data = betadisp.UWUF.time.forstat))
capture.output(ANOVA.betadisp.UWUF, file = "./BETAPLOTS_H/ANOVA.betadisp.UWUF.txt")

xvalues <- c("20","60","360","720","5760")
ylimits <- c(0.3,0.5)
pdf(paste0("./BETAPLOTS_H/BetaDisp_H_",metric,"_eachmorph.pdf"),pointsize = 14)
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = 'Time'
     , ylab = "Distance (Unweighted Unifrac)"
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
dm.UWUF.5760<- dm.UWUF.inclWater[grep("(CR|FB|BL)-5760", rownames(dm.UWUF.inclWater)),grep("(CR|FB|BL)-5760", colnames(dm.UWUF.inclWater))]

NMDS.UWUF.5760 <- isoMDS(as.matrix(dm.UWUF.5760), y = cmdscale(as.matrix(dm.UWUF.5760), 2))

# NMDS.UWUF.morphallTimePoints <- NMDS.UWUF.morphAllTime$points[grep("(CR|BL|FB)-5760", rownames(NMDS.UWUF.morphAllTime$points)),]
MF.5760 <- MF[sapply(rownames(NMDS.UWUF.5760$points), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
}),]
MF.5760$Morph <- factor(MF.5760$Morph, levels = c("CR","BL","FB"))

NMDS.UWUF.5760.CR <- NMDS.UWUF.5760$points[grep("CR", rownames(NMDS.UWUF.5760$points)),]
NMDS.UWUF.5760.CR.chull <- chull(NMDS.UWUF.5760.CR)
NMDS.UWUF.5760.CR.chull <- c(NMDS.UWUF.5760.CR.chull, NMDS.UWUF.5760.CR.chull[1])

NMDS.UWUF.5760.BL <- NMDS.UWUF.5760$points[grep("BL", rownames(NMDS.UWUF.5760$points)),]
NMDS.UWUF.5760.BL.chull <- chull(NMDS.UWUF.5760.BL)
NMDS.UWUF.5760.BL.chull <- c(NMDS.UWUF.5760.BL.chull, NMDS.UWUF.5760.BL.chull[1])

NMDS.UWUF.5760.FB <- NMDS.UWUF.5760$points[grep("FB", rownames(NMDS.UWUF.5760$points)),]
NMDS.UWUF.5760.FB.chull <- chull(NMDS.UWUF.5760.FB)
NMDS.UWUF.5760.FB.chull <- c(NMDS.UWUF.5760.FB.chull, NMDS.UWUF.5760.FB.chull[1])

# ANOVA.UWUF.5760 <- adonis(dm.UWUF.5760 ~ Morph, data = MF.5760)
# # ANOSIM.UWUF.5760 <- anosim(dat = dm.UWUF.5760, grouping = MF.5760$Morph)
# capture.output(ANOVA.UWUF.5760, file = paste0("BETAPLOTS_H/individualtests/anova_",metric,"_5760only.txt"))

pdf(paste0("BETAPLOTS_H/NMDS_",metric,"_5760Only.pdf"), pointsize = 14)
par(fig = c(0,0.75,0,1))
plot(NMDS.UWUF.5760$points
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.5760$Morph)]
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , sub = paste0("Stress: ", round(NMDS.UWUF.5760$stress,2)/100)
     , cex = 1.5
     , main = "NMDS of community composition (4 days)")
lines(NMDS.UWUF.5760.CR[NMDS.UWUF.5760.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.5760.BL[NMDS.UWUF.5760.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.5760.FB[NMDS.UWUF.5760.FB.chull,]
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
rownames(dm.UWUF.P.morphonly) <- gsub(".","-", rownames(dm.UWUF.P.morphonly), fixed = TRUE)
colnames(dm.UWUF.P.morphonly) <- gsub(".","-", colnames(dm.UWUF.P.morphonly), fixed = TRUE)

rownames(MF.P.morphkeep) <- gsub(".","-", rownames(MF.P.morphkeep), fixed = TRUE)
rownames(MF.P.inclWater) <- gsub(".","-", rownames(MF.P.inclWater), fixed = TRUE)



NMDS.UWUF.P.morphonly <- isoMDS(as.matrix(dm.UWUF.P.morphonly), y = cmdscale(as.matrix(dm.UWUF.P.morphonly), 2))
NMDS.UWUF.P.all <- isoMDS(as.matrix(dm.UWUF.P.inclWater), y = cmdscale(as.matrix(dm.UWUF.P.inclWater), 2))

###### STATS ##########
MF.P.morphkeep <- MF.P.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
MF.P.morphkeep$Morph <- factor(MF.P.morphkeep$Morph, levels = c('CR','BL','FB'))
MF.P.morphkeep$Time <- factor(MF.P.morphkeep$Time, levels = c('20','60','180','360','720','1440'))
MF.P.morphkeep$Type <- factor(MF.P.morphkeep$Type, levels = c('P','H'))


# dm.UWUF.P.morphonly.H <- dm.UWUF.P.morphonly[grep(".", rownames(dm.UWUF.P.morphonly), fixed = TRUE),grep(".", colnames(dm.UWUF.P.morphonly), fixed = TRUE)]
# dm.UWUF.P.morphonly.P <- dm.UWUF.P.morphonly[grep("-", rownames(dm.UWUF.P.morphonly), fixed = TRUE),grep("-", colnames(dm.UWUF.P.morphonly), fixed = TRUE)]

ANOVA.UWUF.P.morphtime <- adonis(dm.UWUF.P.morphonly ~ Time*Morph, data = MF.P.morphkeep, by = "margin")
capture.output(ANOVA.UWUF.P.morphtime, file = paste0("BETAPLOTS_P/adonis_", metric,"_PM.txt"))
# ANOSIM.UWUF.P.morphtime <- anosim(dm.UWUF.P.morphonly, grouping = MF.P.morphkeep$Morph)

# Dispersion across time and between morphs
dist.UWUF.P.morphonly <- as.dist(dm.UWUF.P.morphonly)
betadisp.UWUF.P.time <- betadisper(d = dist.UWUF.P.morphonly, group = MF.P.morphkeep$Time)
betadisp.UWUF.P.morph <- betadisper(d = dist.UWUF.P.morphonly, group = MF.P.morphkeep$Morph)

capture.output(betadisp.UWUF.P.time, file = paste0("BETAPLOTS_P/betadispTime_", metric, "_PM.txt"))
capture.output(betadisp.UWUF.P.morph, file = paste0("BETAPLOTS_P/betadispMorph_", metric, "_PM.txt"))


# Grep values for each timepoint and morph, and then plot dispersion within each
TP <- levels(factor(MF.P.morphkeep$Time))
MorphTypes <- levels(factor(MF.P.morphkeep$Morph))

# First Timepoint
FirstTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[1],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.UWUF.P.FirstTP <- dm.UWUF.P.morphonly[FirstTPNames,FirstTPNames]
FB.BL.Firstvalues <- as.vector(as.dist(dm.UWUF.P.FirstTP[grep("FB", rownames(dm.UWUF.P.FirstTP)), grep("BL", colnames(dm.UWUF.P.FirstTP))]))
FB.CR.Firstvalues <- as.vector(as.dist(dm.UWUF.P.FirstTP[grep("FB", rownames(dm.UWUF.P.FirstTP)), grep("CR", colnames(dm.UWUF.P.FirstTP))]))
CR.BL.Firstvalues <- as.vector(as.dist(dm.UWUF.P.FirstTP[grep("CR", rownames(dm.UWUF.P.FirstTP)), grep("BL", colnames(dm.UWUF.P.FirstTP))]))
# ANOVA
dm.UWUF.P.20 <- dm.UWUF.P.morphonly[grep("^20-", rownames(dm.UWUF.P.morphonly)), grep("^20-", colnames(dm.UWUF.P.morphonly))]
MF.P.UWUF.P.20.only <- MF.P.morphkeep[grep("^20-", rownames(MF.P.morphkeep)),]
ANOVA.UWUF.P.20.only <- adonis(dm.UWUF.P.20 ~ Morph, data = MF.P.UWUF.P.20.only, by = "margin")

# Second Timepoint
SecondTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[2],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.UWUF.P.SecondTP <- dm.UWUF.P.morphonly[SecondTPNames,SecondTPNames]
FB.BL.Secondvalues <- as.vector(as.dist(dm.UWUF.P.SecondTP[grep("FB", rownames(dm.UWUF.P.SecondTP)), grep("BL", colnames(dm.UWUF.P.SecondTP))]))
FB.CR.Secondvalues <- as.vector(as.dist(dm.UWUF.P.SecondTP[grep("FB", rownames(dm.UWUF.P.SecondTP)), grep("CR", colnames(dm.UWUF.P.SecondTP))]))
CR.BL.Secondvalues <- as.vector(as.dist(dm.UWUF.P.SecondTP[grep("CR", rownames(dm.UWUF.P.SecondTP)), grep("BL", colnames(dm.UWUF.P.SecondTP))]))
# ANOVA
dm.UWUF.P.60 <- dm.UWUF.P.morphonly[grep("^60-", rownames(dm.UWUF.P.morphonly)), grep("^60-", colnames(dm.UWUF.P.morphonly))]
MF.P.UWUF.P.60.only <- MF.P.morphkeep[grep("^60-", rownames(MF.P.morphkeep)),]
ANOVA.UWUF.P.60.only <- adonis(dm.UWUF.P.60 ~ Morph, data = MF.P.UWUF.P.60.only, by = "margin")

# Third Timepoint
ThirdTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[3],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.UWUF.P.ThirdTP <- dm.UWUF.P.morphonly[ThirdTPNames,ThirdTPNames]
FB.BL.Thirdvalues <- as.vector(as.dist(dm.UWUF.P.ThirdTP[grep("FB", rownames(dm.UWUF.P.ThirdTP)), grep("BL", colnames(dm.UWUF.P.ThirdTP))]))
FB.CR.Thirdvalues <- as.vector(as.dist(dm.UWUF.P.ThirdTP[grep("FB", rownames(dm.UWUF.P.ThirdTP)), grep("CR", colnames(dm.UWUF.P.ThirdTP))]))
CR.BL.Thirdvalues <- as.vector(as.dist(dm.UWUF.P.ThirdTP[grep("CR", rownames(dm.UWUF.P.ThirdTP)), grep("BL", colnames(dm.UWUF.P.ThirdTP))]))
# ANOVA
dm.UWUF.P.180 <- dm.UWUF.P.morphonly[grep("^180-", rownames(dm.UWUF.P.morphonly)), grep("^180-", colnames(dm.UWUF.P.morphonly))]
MF.P.UWUF.P.180.only <- MF.P.morphkeep[grep("^180-", rownames(MF.P.morphkeep)),]
ANOVA.UWUF.P.180.only <- adonis(dm.UWUF.P.180 ~ Morph, data = MF.P.UWUF.P.180.only, by = "margin")


# Fourth Timepoint
FourthTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[4],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.UWUF.P.FourthTP <- dm.UWUF.P.morphonly[FourthTPNames,FourthTPNames]
FB.BL.Fourthvalues <- as.vector(as.dist(dm.UWUF.P.FourthTP[grep("FB", rownames(dm.UWUF.P.FourthTP)), grep("BL", colnames(dm.UWUF.P.FourthTP))]))
FB.CR.Fourthvalues <- as.vector(as.dist(dm.UWUF.P.FourthTP[grep("FB", rownames(dm.UWUF.P.FourthTP)), grep("CR", colnames(dm.UWUF.P.FourthTP))]))
CR.BL.Fourthvalues <- as.vector(as.dist(dm.UWUF.P.FourthTP[grep("CR", rownames(dm.UWUF.P.FourthTP)), grep("BL", colnames(dm.UWUF.P.FourthTP))]))
# ANOVA
dm.UWUF.P.360 <- dm.UWUF.P.morphonly[grep("^360-", rownames(dm.UWUF.P.morphonly)), grep("^360-", colnames(dm.UWUF.P.morphonly))]
MF.P.UWUF.P.360.only <- MF.P.morphkeep[grep("^360-", rownames(MF.P.morphkeep)),]
ANOVA.UWUF.P.360.only <- adonis(dm.UWUF.P.360 ~ Morph, data = MF.P.UWUF.P.360.only, by = "margin")

# Fifth Timepoint
FifthTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[5],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.UWUF.P.FifthTP <- dm.UWUF.P.morphonly[FifthTPNames,FifthTPNames]
FB.BL.Fifthvalues <- as.vector(as.dist(dm.UWUF.P.FifthTP[grep("FB", rownames(dm.UWUF.P.FifthTP)), grep("BL", colnames(dm.UWUF.P.FifthTP))]))
FB.CR.Fifthvalues <- as.vector(as.dist(dm.UWUF.P.FifthTP[grep("FB", rownames(dm.UWUF.P.FifthTP)), grep("CR", colnames(dm.UWUF.P.FifthTP))]))
CR.BL.Fifthvalues <- as.vector(as.dist(dm.UWUF.P.FifthTP[grep("CR", rownames(dm.UWUF.P.FifthTP)), grep("BL", colnames(dm.UWUF.P.FifthTP))]))
# ANOVA
dm.UWUF.P.720 <- dm.UWUF.P.morphonly[grep("^720-", rownames(dm.UWUF.P.morphonly)), grep("^720-", colnames(dm.UWUF.P.morphonly))]
MF.P.UWUF.P.720.only <- MF.P.morphkeep[grep("^720-", rownames(MF.P.morphkeep)),]
ANOVA.UWUF.P.720.only <- adonis(dm.UWUF.P.720 ~ Morph, data = MF.P.UWUF.P.720.only, by = "margin")

# Sixth Timepoint
SixthTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[6],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.UWUF.P.SixthTP <- dm.UWUF.P.morphonly[SixthTPNames,SixthTPNames]
FB.BL.Sixthvalues <- as.vector(as.dist(dm.UWUF.P.SixthTP[grep("FB", rownames(dm.UWUF.P.SixthTP)), grep("BL", colnames(dm.UWUF.P.SixthTP))]))
FB.CR.Sixthvalues <- as.vector(as.dist(dm.UWUF.P.SixthTP[grep("FB", rownames(dm.UWUF.P.SixthTP)), grep("CR", colnames(dm.UWUF.P.SixthTP))]))
CR.BL.Sixthvalues <- as.vector(as.dist(dm.UWUF.P.SixthTP[grep("CR", rownames(dm.UWUF.P.SixthTP)), grep("BL", colnames(dm.UWUF.P.SixthTP))]))
# ANOVA
dm.UWUF.P.1440<- dm.UWUF.P.morphonly[grep("^1440-", rownames(dm.UWUF.P.morphonly)), grep("^1440-", colnames(dm.UWUF.P.morphonly))]
MF.P.UWUF.P.1440.only <- MF.P.morphkeep[grep("^1440-", rownames(MF.P.morphkeep)),]
ANOVA.UWUF.P.1440.only <- adonis(dm.UWUF.P.1440 ~ Morph, data = MF.P.UWUF.P.1440.only, by = "margin")



# Combine into single tables

FB.BL.UWUF.P.ALL <- cbind(as.numeric(c(FB.BL.Firstvalues
                                     , FB.BL.Secondvalues
                                     , FB.BL.Thirdvalues
                                     , FB.BL.Fourthvalues
                                     , FB.BL.Fifthvalues
                                     , FB.BL.Sixthvalues
                                     ))
                        , as.numeric(c(rep(TP[1], length(FB.BL.Firstvalues))
                                       , rep(TP[2], length(FB.BL.Secondvalues))
                                       , rep(TP[3], length(FB.BL.Thirdvalues))
                                       , rep(TP[4], length(FB.BL.Fourthvalues))
                                       , rep(TP[5], length(FB.BL.Fifthvalues))
                                       , rep(TP[6], length(FB.BL.Sixthvalues))
                                       )))
FB.CR.UWUF.P.ALL <- cbind(as.numeric(c(FB.CR.Firstvalues
                                     , FB.CR.Secondvalues
                                     , FB.CR.Thirdvalues
                                     , FB.CR.Fourthvalues
                                     , FB.CR.Fifthvalues
                                     , FB.CR.Sixthvalues
                                     ))
                        , as.numeric(c(rep(TP[1], length(FB.CR.Firstvalues))
                                       , rep(TP[2], length(FB.CR.Secondvalues))
                                       , rep(TP[3], length(FB.CR.Thirdvalues))
                                       , rep(TP[4], length(FB.CR.Fourthvalues))
                                       , rep(TP[5], length(FB.CR.Fifthvalues))
                                       , rep(TP[6], length(FB.CR.Sixthvalues))
                                       )))
CR.BL.UWUF.P.ALL <- cbind(as.numeric(c(CR.BL.Firstvalues
                                     , CR.BL.Secondvalues
                                     , CR.BL.Thirdvalues
                                     , CR.BL.Fourthvalues
                                     , CR.BL.Fifthvalues
                                     , CR.BL.Sixthvalues
                                     ))
                        , as.numeric(c(rep(TP[1], length(CR.BL.Firstvalues))
                                       , rep(TP[2], length(CR.BL.Secondvalues))
                                       , rep(TP[3], length(CR.BL.Thirdvalues))
                                       , rep(TP[4], length(CR.BL.Fourthvalues))
                                       , rep(TP[5], length(CR.BL.Fifthvalues))
                                       , rep(TP[6], length(CR.BL.Sixthvalues))
                                       )))

# FB.BL.UWUF.P.lm <- lm(FB.BL.UWUF.P.ALL[,1] ~ log(FB.BL.UWUF.P.ALL[,2]))
# summary(FB.BL.UWUF.P.lm)
# 
# CR.BL.UWUF.P.lm <- lm(CR.BL.UWUF.P.ALL[,1] ~ log(CR.BL.UWUF.P.ALL[,2]))
# summary(CR.BL.UWUF.P.lm)
# 
# FB.CR.UWUF.P.lm <- lm(FB.CR.UWUF.P.ALL[,1] ~ log(FB.CR.UWUF.P.ALL[,2]))
# summary(FB.CR.UWUF.P.lm)

# Plot the distances

FB.BL.UWUF.P.ALL.mean <- aggregate(FB.BL.UWUF.P.ALL, by = list(FB.BL.UWUF.P.ALL[,2]), mean)
FB.BL.UWUF.P.ALL.sd <- aggregate(FB.BL.UWUF.P.ALL, by = list(FB.BL.UWUF.P.ALL[,2]), sd)

CR.BL.UWUF.P.ALL.mean <- aggregate(CR.BL.UWUF.P.ALL, by = list(CR.BL.UWUF.P.ALL[,2]), mean)
CR.BL.UWUF.P.ALL.sd <- aggregate(CR.BL.UWUF.P.ALL, by = list(CR.BL.UWUF.P.ALL[,2]), sd)

FB.CR.UWUF.P.ALL.mean <- aggregate(FB.CR.UWUF.P.ALL, by = list(FB.CR.UWUF.P.ALL[,2]), mean)
FB.CR.UWUF.P.ALL.sd <- aggregate(FB.CR.UWUF.P.ALL, by = list(FB.CR.UWUF.P.ALL[,2]), sd)

######## PLOT DISP ###############
# ylimits <- c(min(FB.BL.UWUF.P.ALL[,1],FB.CR.UWUF.P.ALL[,1],CR.BL.UWUF.P.ALL[,1]), max(FB.BL.UWUF.P.ALL[,1],FB.CR.UWUF.P.ALL[,1],CR.BL.UWUF.P.ALL[,1]))
ylimits <- c(0.3,0.7)

# xvalues <- log(FB.BL.UWUF.P.ALL.mean[,1])
xvalues <- as.character(FB.BL.UWUF.P.ALL.mean[,1])
pdf(paste0("BETAPLOTS_P/DispOverTime_P_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1, at = c(1,2,3,4,5,6), labels = c('20 min','1 h','3 h','6 h', '12 h', '24 h'))
points(FB.BL.UWUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*0.99
       , x1 = c(1,2,3,4,5,6)*0.99
       , y0 = c(FB.BL.UWUF.P.ALL.mean[,2] - FB.BL.UWUF.P.ALL.sd[,2]/2)
       , y1 = c(FB.BL.UWUF.P.ALL.mean[,2] + FB.BL.UWUF.P.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.UWUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*1
       , x1 = c(1,2,3,4,5,6)*1
       , y0 = c(FB.CR.UWUF.P.ALL.mean[,2] - FB.CR.UWUF.P.ALL.sd[,2]/2)
       , y1 = c(FB.CR.UWUF.P.ALL.mean[,2] + FB.CR.UWUF.P.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.UWUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*1.01
       , x1 = c(1,2,3,4,5,6)*1.01
       , y0 = c(CR.BL.UWUF.P.ALL.mean[,2] - CR.BL.UWUF.P.ALL.sd[,2]/2)
       , y1 = c(CR.BL.UWUF.P.ALL.mean[,2] + CR.BL.UWUF.P.ALL.sd[,2]/2)
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
# NOCHLP UWUF.P

# SWITCH AROUND AXIS
NMDS.UWUF.P.morphonly$points[,1] <- -(NMDS.UWUF.P.morphonly$points[,1])


# MAKE POLYGONS for plotting
NMDS.UWUF.P.CR <- NMDS.UWUF.P.morphonly$points[grep("CR", rownames(NMDS.UWUF.P.morphonly$points)),]
NMDS.UWUF.P.BL <- NMDS.UWUF.P.morphonly$points[grep("BL", rownames(NMDS.UWUF.P.morphonly$points)),]
NMDS.UWUF.P.FB <- NMDS.UWUF.P.morphonly$points[grep("FB", rownames(NMDS.UWUF.P.morphonly$points)),]

NMDS.UWUF.P.CR.chull <- chull(NMDS.UWUF.P.CR)
NMDS.UWUF.P.CR.chull <- c(NMDS.UWUF.P.CR.chull, NMDS.UWUF.P.CR.chull[1])

NMDS.UWUF.P.BL.chull <- chull(NMDS.UWUF.P.BL)
NMDS.UWUF.P.BL.chull <- c(NMDS.UWUF.P.BL.chull, NMDS.UWUF.P.BL.chull[1])

NMDS.UWUF.P.FB.chull <- chull(NMDS.UWUF.P.FB)
NMDS.UWUF.P.FB.chull <- c(NMDS.UWUF.P.FB.chull, NMDS.UWUF.P.FB.chull[1])


MorphColours <- c("darkorchid4","dodgerblue","salmon") 

pdf(paste0("BETAPLOTS_P/NMDS_P_",metric,"_Morph.pdf"), pointsize = 14)
par(fig = c(0,0.8,0,1))
plot(NMDS.UWUF.P.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.morphkeep$Morph)]
     , sub = round(NMDS.UWUF.P.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
)
lines(NMDS.UWUF.P.CR[NMDS.UWUF.P.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.BL[NMDS.UWUF.P.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.FB[NMDS.UWUF.P.FB.chull,]
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
plot(NMDS.UWUF.P.morphonly$points[,1]
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = TimeColours[factor(MF.P.morphkeep$Time)]
     , sub = round(NMDS.UWUF.P.morphonly$stress/100,2)
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



# TIME-- ONE at a time


# 20 minutes
NMDS.UWUF.P.20.only <- NMDS.UWUF.P.morphonly$points[grep("^20-", rownames(NMDS.UWUF.P.morphonly$points)),]

NMDS.UWUF.P.20.only.CR <- NMDS.UWUF.P.20.only[grep("CR", rownames(NMDS.UWUF.P.20.only)),]
NMDS.UWUF.P.20.only.CR.chull <- chull(NMDS.UWUF.P.20.only.CR)
NMDS.UWUF.P.20.only.CR.chull <- c(NMDS.UWUF.P.20.only.CR.chull, NMDS.UWUF.P.20.only.CR.chull[1])

NMDS.UWUF.P.20.only.BL <- NMDS.UWUF.P.20.only[grep("BL", rownames(NMDS.UWUF.P.20.only)),]
NMDS.UWUF.P.20.only.BL.chull <- chull(NMDS.UWUF.P.20.only.BL)
NMDS.UWUF.P.20.only.BL.chull <- c(NMDS.UWUF.P.20.only.BL.chull, NMDS.UWUF.P.20.only.BL.chull[1])

NMDS.UWUF.P.20.only.FB <- NMDS.UWUF.P.20.only[grep("FB", rownames(NMDS.UWUF.P.20.only)),]
NMDS.UWUF.P.20.only.FB.chull <- chull(NMDS.UWUF.P.20.only.FB)
NMDS.UWUF.P.20.only.FB.chull <- c(NMDS.UWUF.P.20.only.FB.chull, NMDS.UWUF.P.20.only.FB.chull[1])

# 60 minutes
NMDS.UWUF.P.60.only <- NMDS.UWUF.P.morphonly$points[grep("^60-", rownames(NMDS.UWUF.P.morphonly$points)),]

NMDS.UWUF.P.60.only.CR <- NMDS.UWUF.P.60.only[grep("CR", rownames(NMDS.UWUF.P.60.only)),]
NMDS.UWUF.P.60.only.CR.chull <- chull(NMDS.UWUF.P.60.only.CR)
NMDS.UWUF.P.60.only.CR.chull <- c(NMDS.UWUF.P.60.only.CR.chull, NMDS.UWUF.P.60.only.CR.chull[1])

NMDS.UWUF.P.60.only.BL <- NMDS.UWUF.P.60.only[grep("BL", rownames(NMDS.UWUF.P.60.only)),]
NMDS.UWUF.P.60.only.BL.chull <- chull(NMDS.UWUF.P.60.only.BL)
NMDS.UWUF.P.60.only.BL.chull <- c(NMDS.UWUF.P.60.only.BL.chull, NMDS.UWUF.P.60.only.BL.chull[1])

NMDS.UWUF.P.60.only.FB <- NMDS.UWUF.P.60.only[grep("FB", rownames(NMDS.UWUF.P.60.only)),]
NMDS.UWUF.P.60.only.FB.chull <- chull(NMDS.UWUF.P.60.only.FB)
NMDS.UWUF.P.60.only.FB.chull <- c(NMDS.UWUF.P.60.only.FB.chull, NMDS.UWUF.P.60.only.FB.chull[1])

# 180 minutes
NMDS.UWUF.P.180.only <- NMDS.UWUF.P.morphonly$points[grep("^180-", rownames(NMDS.UWUF.P.morphonly$points)),]

NMDS.UWUF.P.180.only.CR <- NMDS.UWUF.P.180.only[grep("CR", rownames(NMDS.UWUF.P.180.only)),]
NMDS.UWUF.P.180.only.CR.chull <- chull(NMDS.UWUF.P.180.only.CR)
NMDS.UWUF.P.180.only.CR.chull <- c(NMDS.UWUF.P.180.only.CR.chull, NMDS.UWUF.P.180.only.CR.chull[1])

NMDS.UWUF.P.180.only.BL <- NMDS.UWUF.P.180.only[grep("BL", rownames(NMDS.UWUF.P.180.only)),]
NMDS.UWUF.P.180.only.BL.chull <- chull(NMDS.UWUF.P.180.only.BL)
NMDS.UWUF.P.180.only.BL.chull <- c(NMDS.UWUF.P.180.only.BL.chull, NMDS.UWUF.P.180.only.BL.chull[1])

NMDS.UWUF.P.180.only.FB <- NMDS.UWUF.P.180.only[grep("FB", rownames(NMDS.UWUF.P.180.only)),]
NMDS.UWUF.P.180.only.FB.chull <- chull(NMDS.UWUF.P.180.only.FB)
NMDS.UWUF.P.180.only.FB.chull <- c(NMDS.UWUF.P.180.only.FB.chull, NMDS.UWUF.P.180.only.FB.chull[1])



# 360 minutes
NMDS.UWUF.P.360.only <- NMDS.UWUF.P.morphonly$points[grep("^360-", rownames(NMDS.UWUF.P.morphonly$points)),]


NMDS.UWUF.P.360.only.CR <- NMDS.UWUF.P.360.only[grep("CR", rownames(NMDS.UWUF.P.360.only)),]
NMDS.UWUF.P.360.only.CR.chull <- chull(NMDS.UWUF.P.360.only.CR)
NMDS.UWUF.P.360.only.CR.chull <- c(NMDS.UWUF.P.360.only.CR.chull, NMDS.UWUF.P.360.only.CR.chull[1])

NMDS.UWUF.P.360.only.BL <- NMDS.UWUF.P.360.only[grep("BL", rownames(NMDS.UWUF.P.360.only)),]
NMDS.UWUF.P.360.only.BL.chull <- chull(NMDS.UWUF.P.360.only.BL)
NMDS.UWUF.P.360.only.BL.chull <- c(NMDS.UWUF.P.360.only.BL.chull, NMDS.UWUF.P.360.only.BL.chull[1])

NMDS.UWUF.P.360.only.FB <- NMDS.UWUF.P.360.only[grep("FB", rownames(NMDS.UWUF.P.360.only)),]
NMDS.UWUF.P.360.only.FB.chull <- chull(NMDS.UWUF.P.360.only.FB)
NMDS.UWUF.P.360.only.FB.chull <- c(NMDS.UWUF.P.360.only.FB.chull, NMDS.UWUF.P.360.only.FB.chull[1])


# 720 minutes
NMDS.UWUF.P.720.only <- NMDS.UWUF.P.morphonly$points[grep("^720-", rownames(NMDS.UWUF.P.morphonly$points)),]


NMDS.UWUF.P.720.only.CR <- NMDS.UWUF.P.720.only[grep("CR", rownames(NMDS.UWUF.P.720.only)),]
NMDS.UWUF.P.720.only.CR.chull <- chull(NMDS.UWUF.P.720.only.CR)
NMDS.UWUF.P.720.only.CR.chull <- c(NMDS.UWUF.P.720.only.CR.chull, NMDS.UWUF.P.720.only.CR.chull[1])

NMDS.UWUF.P.720.only.BL <- NMDS.UWUF.P.720.only[grep("BL", rownames(NMDS.UWUF.P.720.only)),]
NMDS.UWUF.P.720.only.BL.chull <- chull(NMDS.UWUF.P.720.only.BL)
NMDS.UWUF.P.720.only.BL.chull <- c(NMDS.UWUF.P.720.only.BL.chull, NMDS.UWUF.P.720.only.BL.chull[1])

NMDS.UWUF.P.720.only.FB <- NMDS.UWUF.P.720.only[grep("FB", rownames(NMDS.UWUF.P.720.only)),]
NMDS.UWUF.P.720.only.FB.chull <- chull(NMDS.UWUF.P.720.only.FB)
NMDS.UWUF.P.720.only.FB.chull <- c(NMDS.UWUF.P.720.only.FB.chull, NMDS.UWUF.P.720.only.FB.chull[1])


# 1440 minutes
NMDS.UWUF.P.1440.only <- NMDS.UWUF.P.morphonly$points[grep("^1440-", rownames(NMDS.UWUF.P.morphonly$points)),]


NMDS.UWUF.P.1440.only.CR <- NMDS.UWUF.P.1440.only[grep("CR", rownames(NMDS.UWUF.P.1440.only)),]
NMDS.UWUF.P.1440.only.CR.chull <- chull(NMDS.UWUF.P.1440.only.CR)
NMDS.UWUF.P.1440.only.CR.chull <- c(NMDS.UWUF.P.1440.only.CR.chull, NMDS.UWUF.P.1440.only.CR.chull[1])

NMDS.UWUF.P.1440.only.BL <- NMDS.UWUF.P.1440.only[grep("BL", rownames(NMDS.UWUF.P.1440.only)),]
NMDS.UWUF.P.1440.only.BL.chull <- chull(NMDS.UWUF.P.1440.only.BL)
NMDS.UWUF.P.1440.only.BL.chull <- c(NMDS.UWUF.P.1440.only.BL.chull, NMDS.UWUF.P.1440.only.BL.chull[1])

NMDS.UWUF.P.1440.only.FB <- NMDS.UWUF.P.1440.only[grep("FB", rownames(NMDS.UWUF.P.1440.only)),]
NMDS.UWUF.P.1440.only.FB.chull <- chull(NMDS.UWUF.P.1440.only.FB)
NMDS.UWUF.P.1440.only.FB.chull <- c(NMDS.UWUF.P.1440.only.FB.chull, NMDS.UWUF.P.1440.only.FB.chull[1])


pdf(paste0("BETAPLOTS_P/NMDS_",metric,"_TimebyMorph.pdf"),pointsize = 14, width = 7, height = 4)
par(mfrow= c(1,4), oma = c(4,4,4,6))
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.P.20.only
     , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.P.20.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.P.20.only.CR[NMDS.UWUF.P.20.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.20.only.BL[NMDS.UWUF.P.20.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.20.only.FB[NMDS.UWUF.P.20.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.P.60.only
     , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.UWUF.P.60.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.P.60.only.CR[NMDS.UWUF.P.60.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.60.only.BL[NMDS.UWUF.P.60.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.60.only.FB[NMDS.UWUF.P.60.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.P.180.only
     , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.180.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.UWUF.P.180.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.P.180.only.CR[NMDS.UWUF.P.180.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.180.only.BL[NMDS.UWUF.P.180.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.180.only.FB[NMDS.UWUF.P.180.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.P.360.only
     , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.P.360.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.P.360.only.CR[NMDS.UWUF.P.360.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.360.only.BL[NMDS.UWUF.P.360.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.360.only.FB[NMDS.UWUF.P.360.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.P.720.only
     , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.P.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.P.720.only.CR[NMDS.UWUF.P.720.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.720.only.BL[NMDS.UWUF.P.720.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.720.only.FB[NMDS.UWUF.P.720.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.UWUF.P.1440.only
     , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.1440.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.UWUF.P.1440.only.CR[NMDS.UWUF.P.1440.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.1440.only.BL[NMDS.UWUF.P.1440.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.1440.only.FB[NMDS.UWUF.P.1440.only.FB.chull,]
      , col = MorphColours[3])
#STOP
dev.off()


############ COMBO DISP BETA ################
# Disp and beta through time combined
# xvalues <- log(FB.BL.UWUF.P.ALL.mean[,1])
xvalues <- as.character(FB.BL.UWUF.P.ALL.mean[,1])

pdf(paste0("BETAPLOTS_P/COMBO_P_dispbeta_",metric,".pdf"), pointsize = 14, width = 14, height = 8)
par(fig = c(0,0.8,0.23,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = ''
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1
     , at = c(1,2,3,4,5,6)
     , labels = c("20 min","1 h","3 h","6 h","12 h","24 h")
     , las = 2)
points(FB.BL.UWUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*0.99
       , x1 = c(1,2,3,4,5,6)*0.99
       , y0 = c(FB.BL.UWUF.P.ALL.mean[,2] - FB.BL.UWUF.P.ALL.sd[,2])
       , y1 = c(FB.BL.UWUF.P.ALL.mean[,2] + FB.BL.UWUF.P.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.UWUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*1
       , x1 = c(1,2,3,4,5,6)*1
       , y0 = c(FB.CR.UWUF.P.ALL.mean[,2] - FB.CR.UWUF.P.ALL.sd[,2])
       , y1 = c(FB.CR.UWUF.P.ALL.mean[,2] + FB.CR.UWUF.P.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.UWUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1.01
       , x1 = c(1,2,3,4)*1.01
       , y0 = c(CR.BL.UWUF.P.ALL.mean[,2] - CR.BL.UWUF.P.ALL.sd[,2])
       , y1 = c(CR.BL.UWUF.P.ALL.mean[,2] + CR.BL.UWUF.P.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
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
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.05,0.158333,0,0.3), new = TRUE)
plot(NMDS.UWUF.P.20.only
     # , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.P.20.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.20.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.20.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.P.20.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.P.20.only.CR[NMDS.UWUF.P.20.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.20.only.BL[NMDS.UWUF.P.20.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.20.only.FB[NMDS.UWUF.P.20.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.178333,0.28666,0,0.3), new = TRUE)
plot(NMDS.UWUF.P.60.only
     # , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.P.60.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.60.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.60.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.P.60.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.P.60.only.CR[NMDS.UWUF.P.60.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.60.only.BL[NMDS.UWUF.P.60.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.60.only.FB[NMDS.UWUF.P.60.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.30666,0.414999,0,0.3), new = TRUE)
plot(NMDS.UWUF.P.180.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.180.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.P.180.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.180.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.180.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.P.180.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.P.180.only.CR[NMDS.UWUF.P.180.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.180.only.BL[NMDS.UWUF.P.180.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.180.only.FB[NMDS.UWUF.P.180.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.434999,0.543333,0,0.3), new = TRUE)
plot(NMDS.UWUF.P.360.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.P.360.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.360.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.360.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.P.360.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.P.360.only.CR[NMDS.UWUF.P.360.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.360.only.BL[NMDS.UWUF.P.360.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.360.only.FB[NMDS.UWUF.P.360.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.565555,0.671666,0,0.3), new = TRUE)
plot(NMDS.UWUF.P.720.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.P.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.P.720.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.720.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.720.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.P.720.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.P.720.only.CR[NMDS.UWUF.P.720.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.720.only.BL[NMDS.UWUF.P.720.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.720.only.FB[NMDS.UWUF.P.720.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.691666,0.79999,0,0.3), new = TRUE)
plot(NMDS.UWUF.P.1440.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.UWUF.P.1440.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.UWUF.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.UWUF.P.1440.only$aov.tab[1]$Df[1],",",ANOVA.UWUF.P.1440.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.UWUF.P.1440.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.UWUF.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.UWUF.P.1440.only.CR[NMDS.UWUF.P.1440.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.UWUF.P.1440.only.BL[NMDS.UWUF.P.1440.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.UWUF.P.1440.only.FB[NMDS.UWUF.P.1440.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(2,0,2,0), fig = c(0.775,1,0,0.3), new = TRUE)
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

allindividualTests <- c("ANOVA.UWUF.P.20.only"
                        ,"ANOVA.UWUF.P.60.only"
                        ,"ANOVA.UWUF.P.180.only"
                        ,"ANOVA.UWUF.P.360.only"
                        ,"ANOVA.UWUF.P.720.only"
                        ,"ANOVA.UWUF.P.1440.only"
                        
)
for (i in allindividualTests) {
  # print(get(i))
  capture.output(get(i), file = paste0("./BETAPLOTS_P/individualtests/",i,".txt"))
}




########### BETADISP#############
betadisp.UWUF
betadisp.P.UWUF.FB <- betadisp.UWUF.P.time$distances[grep("FB", names(betadisp.UWUF.P.time$distances))]
betadisp.P.UWUF.BL <- betadisp.UWUF.P.time$distances[grep("BL", names(betadisp.UWUF.P.time$distances))]
betadisp.P.UWUF.CR <- betadisp.UWUF.P.time$distances[grep("CR", names(betadisp.UWUF.P.time$distances))]

MF.P.FB <- MF.P.morphkeep[grep("FB", rownames(MF.P.morphkeep)),]
MF.P.BL <- MF.P.morphkeep[grep("BL", rownames(MF.P.morphkeep)),]
MF.P.CR <- MF.P.morphkeep[grep("CR", rownames(MF.P.morphkeep)),]

betadisp.FB.P.agg <- aggregate(betadisp.P.UWUF.FB, by = list(MF.P.FB$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.BL.P.agg <- aggregate(betadisp.P.UWUF.BL, by = list(MF.P.BL$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.CR.P.agg <- aggregate(betadisp.P.UWUF.CR, by = list(MF.P.CR$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )

betadisp.UWUF.time.forstat <- cbind(betadisp.UWUF.time$distances, MF.morphkeep[unlist(lapply(names(betadisp.UWUF.time$distances), function(x) {grep(paste0("^",x,"$"), rownames(MF.morphkeep))})),])
colnames(betadisp.UWUF.time.forstat)[1] <- c("Distance")
ANOVA.betadisp.UWUF <- anova(lm(Distance ~ Time*Morph, data = betadisp.UWUF.time.forstat))
capture.output(ANOVA.betadisp.UWUF, file = "./BETAPLOTS_P/ANOVA.betadisp.UWUF.txt")

ylimits <- c(0.25,0.45)
pdf(paste0("./BETAPLOTS_P/BetaDisp_P_",metric,"_eachmorph.pdf"),pointsize = 14)
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = 'Time'
     , ylab = "Distance (Unweighted Unifrac)"
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
anova.betadisp.UWUF.morph <- anova(betadisp.UWUF.morph)
anova.betadisp.UWUF.time <- anova(betadisp.UWUF.time)
anova.betadisp.P.UWUF.morph <- anova(betadisp.UWUF.P.morph)
anova.betadisp.P.UWUF.time <- anova(betadisp.UWUF.P.time)

# This is PERMADISP-- don't need a table, will quote in text
capture.output(anova.betadisp.UWUF.morph, file = "BETAPLOTS_H/anova.betadisp.UWUF.morph.txt")
capture.output(anova.betadisp.UWUF.time, file = "BETAPLOTS_H/anova.betadisp.UWUF.time.txt")
capture.output(anova.betadisp.P.UWUF.morph, file = "BETAPLOTS_P/anova.betadisp.UWUF.P.morph.txt")
capture.output(anova.betadisp.P.UWUF.time, file = "BETAPLOTS_P/anova.betadisp.UWUF.P.time.txt")

### NOW DO EACH TIME POINT ##
## This is a table with both Hakai and PM for the supplementary figures
permdisp.morphology.across.time <- matrix(ncol = 8, nrow = 2)
rownames(permdisp.morphology.across.time) <- c("P","H")
colnames(permdisp.morphology.across.time) <- c("20","60","180","360","720","1440","5760", "Overall")
for (t in c("20","60","180","360","720","1440")) {
  assign(paste0("betadisp.P.",t), betadisper(dist(get(paste0("dm.UWUF.P.",t))), group = get(paste0("MF.P.UWUF.P.",t,".only"))$Morph))
  assign(paste0("anova.betadisp.P.",t), anova(get(paste0("betadisp.P.",t))))
  ptemp <- get(paste0("anova.betadisp.P.",t))$`Pr(>F)`[1]
  ftemp <- get(paste0("anova.betadisp.P.",t))$`F value`[1]
  dftemp <- get(paste0("anova.betadisp.P.",t))$Df[1]
  
  toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
  permdisp.morphology.across.time["P", paste0(t)] <- toPaste
 }
for (t in c("20","60","360","720","5760")) {
  assign(paste0("betadisp.",t), betadisper(dist(get(paste0("dm.UWUF.",t))), group = get(paste0("MF.UWUF.",t,".only"))$Morph))
  assign(paste0("anova.betadisp.",t), anova(get(paste0("betadisp.",t))))
  ptemp <- get(paste0("anova.betadisp.",t))$`Pr(>F)`[1]
  ftemp <- get(paste0("anova.betadisp.",t))$`F value`[1]
  dftemp <- get(paste0("anova.betadisp.",t))$Df[1]
  
  toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
  permdisp.morphology.across.time["H", paste0(t)] <- toPaste
  
}
# Get overall
ptemp <- anova.betadisp.P.UWUF.morph$`Pr(>F)`[1]
ftemp <- anova.betadisp.P.UWUF.morph$`F value`[1]
dftemp <- anova.betadisp.P.UWUF.morph$Df[1]
toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
permdisp.morphology.across.time["P","Overall"] <- toPaste

ptemp <- anova.betadisp.UWUF.morph$`Pr(>F)`[1]
ftemp <- anova.betadisp.UWUF.morph$`F value`[1]
dftemp <- anova.betadisp.UWUF.morph$Df[1]
toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
permdisp.morphology.across.time["H","Overall"] <- toPaste

colnames(permdisp.morphology.across.time) <- c("20 minutes","1 hour","3 hours","6 hours","12 hours","1 day","4 days", "Overall")

tempfile1 <- permdisp.morphology.across.time
for (r in 1:nrow(tempfile1)) {
  for (c in 1:ncol(tempfile1)) {
    if (is.na(tempfile1[r,c])) {
      permdisp.morphology.across.time[r,c] <- "-"
    }}}


# Make double header
permdisp.morphology.across.time.UWUF <- permdisp.morphology.across.time
rownames(permdisp.morphology.across.time.UWUF) <- c("Reed Point","Hakai")

capture.output(xtable(permdisp.morphology.across.time.UWUF, digits = NULL), file = paste0("BETAPLOTS_LATEX/permdisp.morph.across.time.",metric,".txt"))

### MAKE BETA DIV TABLES-- metrics separately but H and P together; extras will go in supp

anova.morphology.across.time <- matrix(ncol = 8, nrow = 2)
colnames(anova.morphology.across.time) <- c("20","60","180","360","720","1440","5760", "Overall")
rownames(anova.morphology.across.time) <- c("P","H")
for (t in c("20","60","180","360","720","1440")) {
  ptemp <- get(paste0("ANOVA.UWUF.P.",t,".only"))$aov.tab$`Pr(>F)`[1]
  rtemp <- get(paste0("ANOVA.UWUF.P.",t,".only"))$aov.tab$`R2`[1]
  dftemp <- get(paste0("ANOVA.UWUF.P.",t,".only"))$aov.tab$`Df`[1]
  toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
  
  anova.morphology.across.time["P",paste0(t)] <- toPaste
}
for (t in c("20","60","360","720","5760")) {
  ptemp <- get(paste0("ANOVA.UWUF.",t,".only"))$aov.tab$`Pr(>F)`[1]
  rtemp <- get(paste0("ANOVA.UWUF.",t,".only"))$aov.tab$`R2`[1]
  dftemp <- get(paste0("ANOVA.UWUF.",t,".only"))$aov.tab$`Df`[1]
  toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
  
  anova.morphology.across.time["H",paste0(t)] <- toPaste
 }
# Do overall P
ptemp <- ANOVA.UWUF.P.morphtime$aov.tab$`Pr(>F)`[2]
rtemp <- ANOVA.UWUF.P.morphtime$aov.tab$`R2`[2]
dftemp <- ANOVA.UWUF.P.morphtime$aov.tab$`Df`[2]
toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
anova.morphology.across.time["P","Overall"] <- toPaste

# Do overall H
ptemp <- ANOVA.UWUF.morphtime$aov.tab$`Pr(>F)`[2]
rtemp <- ANOVA.UWUF.morphtime$aov.tab$`R2`[2]
dftemp <- ANOVA.UWUF.morphtime$aov.tab$`Df`[2]
toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
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
anova.morphology.across.time.UWUF <- anova.morphology.across.time
rownames(anova.morphology.across.time.UWUF) <- c("Reed Point","Hakai")

capture.output(xtable(anova.morphology.across.time.UWUF, digits = NULL), file = paste0("BETAPLOTS_LATEX/anova.morph.across.time.",metric,".txt"))

### DO FB/CR/BL TEST FOR EACH
# Make table
pairwiseAdonis.all <- matrix(ncol = 7,nrow = 6)
colnames(pairwiseAdonis.all) <- c("20","60","180","360","720","1440", "5760")
rownames(pairwiseAdonis.all) <- c("pFB:BL","pFB:CR","pBL:CR","hFB:BL","hFB:CR","hBL:CR" )
listMorphs <- c("FB","BL","CR")
for (m in 1:(length(listMorphs)-1)) {
  for (n in (m+1):length(listMorphs)) {
    for (t in c("20","60","180","360","720","1440")) {
      tempMF <- get(paste0("MF.P.UWUF.P.",t,".only"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), get(paste0("MF.P.UWUF.P.",t,".only"))$Morph),]
      tempDM <- get(paste0("dm.UWUF.P.",t))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.UWUF.P.",t)))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.UWUF.P.",t))))]
      tempAdonis <- adonis(tempDM ~ Morph, data = tempMF, by = "marginal")
      toPaste <- paste0( tempAdonis$aov.tab$`Pr(>F)`[1]
                         ," (R^2 = ", round(tempAdonis$aov.tab$R2[1],digits = 2)
                         ,", Df = ", tempAdonis$aov.tab$Df[1] 
                         , ")")
      pairwiseAdonis.all[paste0("p",listMorphs[m],":",listMorphs[n]),paste0(t)] <- toPaste
    }
    for (t in c("20","60","360","720","5760")) {
      tempMF <- get(paste0("MF.UWUF.",t,".only"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), get(paste0("MF.UWUF.",t,".only"))$Morph),]
      tempDM <- get(paste0("dm.UWUF.",t))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.UWUF.",t)))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.UWUF.",t))))]
      tempAdonis <- adonis(tempDM ~ Morph, data = tempMF, by = "marginal")
      toPaste <- paste0( tempAdonis$aov.tab$`Pr(>F)`[1]
                         ," (R^2 = ", round(tempAdonis$aov.tab$R2[1],digits = 2)
                         ,", Df = ", tempAdonis$aov.tab$Df[1] 
                         , ")")
      pairwiseAdonis.all[paste0("h",listMorphs[m],":",listMorphs[n]),paste0(t)] <- toPaste
    }
  }
}

# Get rid of NAs
for (r in 1:nrow(pairwiseAdonis.all)) {
  for (c in 1:ncol(pairwiseAdonis.all)) {
    if (is.na(pairwiseAdonis.all[r,c])) {
      pairwiseAdonis.all[r,c] <- "-"
    }
  }
}

# Change Rownames
pairwiseAdonis.all.UWUF <- cbind(c("FB:BL","FB:CR","BL:CR","FB:BL","FB:CR","BL:CR"), pairwiseAdonis.all)
rownames(pairwiseAdonis.all.UWUF) <- c("Reed Point", " ","  ","Hakai","   ","    ")

capture.output(xtable(pairwiseAdonis.all.UWUF), file = paste0("BETAPLOTS_LATEX/pairwiseAdonis.",metric,".txt"))

####### *********WUF********* ############# 
metric <- "WUF"

######## --HAKAI-- ###########
### NMDS #####

NMDS.WUF.morphonly <- isoMDS(as.matrix(dm.WUF.morphonly), y = cmdscale(as.matrix(dm.WUF.morphonly), 2))
NMDS.WUF.all <- isoMDS(as.matrix(dm.WUF.inclWater), y = cmdscale(as.matrix(dm.WUF.inclWater), 2))

###### STATS ##########
MF.morphkeep <- MF.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
MF.morphkeep$Morph <- factor(MF.morphkeep$Morph, levels = c('CR','BL','FB'))
MF.morphkeep$Time <- factor(MF.morphkeep$Time, levels = c('20','60','180','360','720','1440','5760'))
MF.morphkeep$Type <- factor(MF.morphkeep$Type, levels = c('P','H'))


# dm.WUF.morphonly.H <- dm.WUF.morphonly[grep(".", rownames(dm.WUF.morphonly), fixed = TRUE),grep(".", colnames(dm.WUF.morphonly), fixed = TRUE)]
# dm.WUF.morphonly.P <- dm.WUF.morphonly[grep("-", rownames(dm.WUF.morphonly), fixed = TRUE),grep("-", colnames(dm.WUF.morphonly), fixed = TRUE)]

ANOVA.WUF.morphtime <- adonis(dm.WUF.morphonly ~ Time*Morph, data = MF.morphkeep, by = "margin")
capture.output(ANOVA.WUF.morphtime, file = paste0("BETAPLOTS_H/adonis_", metric,"_Hakai.txt"))
# ANOSIM.WUF.morphtime <- anosim(dm.WUF.morphonly, grouping = MF.morphkeep$Morph)

# Dispersion across time and between morphs
dist.WUF.morphonly <- as.dist(dm.WUF.inclWater[-grep("W", rownames(dm.WUF.inclWater)), -grep("W", colnames(dm.WUF.inclWater))])
MF.incl5760 <- MF.inclWater[-grep("W", rownames(MF.inclWater)),]
betadisp.WUF.time <- betadisper(d = dist.WUF.morphonly, group = MF.incl5760$Time)
betadisp.WUF.morph <- betadisper(d = dist.WUF.morphonly, group = MF.incl5760$Morph)

capture.output(betadisp.WUF.time, file = paste0("BETAPLOTS_H/betadispTime_", metric, "_Hakai.txt"))
capture.output(betadisp.WUF.morph, file = paste0("BETAPLOTS_H/betadispMorph_", metric, "_Hakai.txt"))


# Grep values for each timepoint and morph, and then plot dispersion within each
TP <- levels(factor(MF.inclWater$Time))
MorphTypes <- levels(factor(MF.morphkeep$Morph))
# First Timepoint
FirstTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[1],"$"), MF.inclWater$Time)]
toremove <- grep("W", FirstTPNames)
if (length(toremove) > 0) {
  FirstTPNames <- FirstTPNames[-toremove]
}
# Between morphs
dm.WUF.FirstTP <- dm.WUF.inclWater[FirstTPNames,FirstTPNames]
FB.BL.Firstvalues <- as.vector(as.dist(dm.WUF.FirstTP[grep("FB", rownames(dm.WUF.FirstTP)), grep("BL", colnames(dm.WUF.FirstTP))]))
FB.CR.Firstvalues <- as.vector(as.dist(dm.WUF.FirstTP[grep("FB", rownames(dm.WUF.FirstTP)), grep("CR", colnames(dm.WUF.FirstTP))]))
CR.BL.Firstvalues <- as.vector(as.dist(dm.WUF.FirstTP[grep("CR", rownames(dm.WUF.FirstTP)), grep("BL", colnames(dm.WUF.FirstTP))]))
# ANOVA
dm.WUF.20 <- dm.WUF.inclWater[grep("(CR|BL|FB)-20-", rownames(dm.WUF.inclWater)), grep("-20-", colnames(dm.WUF.inclWater))]
MF.WUF.20.only <- MF.inclWater[grep("(CR|BL|FB)-20-", rownames(MF.inclWater)),]
ANOVA.WUF.20.only <- adonis(dm.WUF.20 ~ Morph, data = MF.WUF.20.only, by = "margin")

# Second Timepoine
SecondTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[2],"$"), MF.inclWater$Time)]
toremove <- grep("W", SecondTPNames)
if (length(toremove) > 0) {
  SecondTPNames <- SecondTPNames[-toremove]
}
# Between morphs
dm.WUF.SecondTP <- dm.WUF.inclWater[SecondTPNames,SecondTPNames]
FB.BL.Secondvalues <- as.vector(as.dist(dm.WUF.SecondTP[grep("FB", rownames(dm.WUF.SecondTP)), grep("BL", colnames(dm.WUF.SecondTP))]))
FB.CR.Secondvalues <- as.vector(as.dist(dm.WUF.SecondTP[grep("FB", rownames(dm.WUF.SecondTP)), grep("CR", colnames(dm.WUF.SecondTP))]))
CR.BL.Secondvalues <- as.vector(as.dist(dm.WUF.SecondTP[grep("CR", rownames(dm.WUF.SecondTP)), grep("BL", colnames(dm.WUF.SecondTP))]))
# ANOVA
dm.WUF.60 <- dm.WUF.inclWater[grep("(BL|CR|FB)-60-", rownames(dm.WUF.inclWater)), grep("-60-", colnames(dm.WUF.morphonly))]
MF.WUF.60.only <- MF.inclWater[grep("(BL|CR|FB)-60-", rownames(MF.inclWater)),]
ANOVA.WUF.60.only <- adonis(dm.WUF.60 ~ Morph, data = MF.WUF.60.only, by = "margin")

# Third Timepoint
ThirdTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[3],"$"), MF.inclWater$Time)]
toremove <- grep("W", ThirdTPNames)
if (length(toremove) > 0) {
  ThirdTPNames <- ThirdTPNames[-toremove]
}
# Between morphs
dm.WUF.ThirdTP <- dm.WUF.inclWater[ThirdTPNames,ThirdTPNames]
FB.BL.Thirdvalues <- as.vector(as.dist(dm.WUF.ThirdTP[grep("FB", rownames(dm.WUF.ThirdTP)), grep("BL", colnames(dm.WUF.ThirdTP))]))
FB.CR.Thirdvalues <- as.vector(as.dist(dm.WUF.ThirdTP[grep("FB", rownames(dm.WUF.ThirdTP)), grep("CR", colnames(dm.WUF.ThirdTP))]))
CR.BL.Thirdvalues <- as.vector(as.dist(dm.WUF.ThirdTP[grep("CR", rownames(dm.WUF.ThirdTP)), grep("BL", colnames(dm.WUF.ThirdTP))]))
# ANOVA
dm.WUF.360 <- dm.WUF.inclWater[grep("(BL|CR|FB)-360-", rownames(dm.WUF.inclWater)), grep("-360-", colnames(dm.WUF.morphonly))]
MF.WUF.360.only <- MF.inclWater[grep("(BL|CR|FB)-360-", rownames(MF.inclWater)),]
ANOVA.WUF.360.only <- adonis(dm.WUF.360 ~ Morph, data = MF.WUF.360.only, by = "margin")


# Fourth Timepoint
FourthTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[4],"$"), MF.inclWater$Time)]
toremove <- grep("W", FourthTPNames)
if (length(toremove) > 0) {
  FourthTPNames <- FourthTPNames[-toremove]
}
# Between morphs
dm.WUF.FourthTP <- dm.WUF.inclWater[FourthTPNames,FourthTPNames]
FB.BL.Fourthvalues <- as.vector(as.dist(dm.WUF.FourthTP[grep("FB", rownames(dm.WUF.FourthTP)), grep("BL", colnames(dm.WUF.FourthTP))]))
FB.CR.Fourthvalues <- as.vector(as.dist(dm.WUF.FourthTP[grep("FB", rownames(dm.WUF.FourthTP)), grep("CR", colnames(dm.WUF.FourthTP))]))
CR.BL.Fourthvalues <- as.vector(as.dist(dm.WUF.FourthTP[grep("CR", rownames(dm.WUF.FourthTP)), grep("BL", colnames(dm.WUF.FourthTP))]))
# ANOVA
dm.WUF.720 <- dm.WUF.inclWater[grep("(BL|CR|FB)-720-", rownames(dm.WUF.inclWater)), grep("-720-", colnames(dm.WUF.morphonly))]
MF.WUF.720.only <- MF.inclWater[grep("(BL|CR|FB)-720-", rownames(MF.inclWater)),]
ANOVA.WUF.720.only <- adonis(dm.WUF.720 ~ Morph, data = MF.WUF.720.only, by = "margin")

# Fifth Timepoint
FifthTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[5],"$"), MF.inclWater$Time)]
toremove <- grep("W", FifthTPNames)
if (length(toremove) > 0) {
  FifthTPNames <- FifthTPNames[-toremove]
}
# Between morphs
dm.WUF.FifthTP <- dm.WUF.inclWater[FifthTPNames,FifthTPNames]
FB.BL.Fifthvalues <- as.vector(as.dist(dm.WUF.FifthTP[grep("FB", rownames(dm.WUF.FifthTP)), grep("BL", colnames(dm.WUF.FifthTP))]))
FB.CR.Fifthvalues <- as.vector(as.dist(dm.WUF.FifthTP[grep("FB", rownames(dm.WUF.FifthTP)), grep("CR", colnames(dm.WUF.FifthTP))]))
CR.BL.Fifthvalues <- as.vector(as.dist(dm.WUF.FifthTP[grep("CR", rownames(dm.WUF.FifthTP)), grep("BL", colnames(dm.WUF.FifthTP))]))
# ANOVA
dm.WUF.5760 <- dm.WUF.inclWater[grep("(BL|CR|FB)-5760-", rownames(dm.WUF.inclWater)), grep("-5760-", colnames(dm.WUF.inclWater))]
MF.WUF.5760.only <- MF.inclWater[grep("(BL|CR|FB)-5760-", rownames(MF.inclWater)),]
ANOVA.WUF.5760.only <- adonis(dm.WUF.5760 ~ Morph, data = MF.WUF.5760.only, by = "margin")


# Combine into single tables

FB.BL.WUF.ALL <- cbind(as.numeric(c(FB.BL.Firstvalues
                                     , FB.BL.Secondvalues
                                     , FB.BL.Thirdvalues
                                     , FB.BL.Fourthvalues
                                     , FB.BL.Fifthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.BL.Firstvalues))
                                       , rep(TP[2], length(FB.BL.Secondvalues))
                                       , rep(TP[3], length(FB.BL.Thirdvalues))
                                       , rep(TP[4], length(FB.BL.Fourthvalues))
                                       , rep(TP[5], length(FB.BL.Fifthvalues)))))
FB.CR.WUF.ALL <- cbind(as.numeric(c(FB.CR.Firstvalues
                                     , FB.CR.Secondvalues
                                     , FB.CR.Thirdvalues
                                     , FB.CR.Fourthvalues
                                     , FB.CR.Fifthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.CR.Firstvalues))
                                       , rep(TP[2], length(FB.CR.Secondvalues))
                                       , rep(TP[3], length(FB.CR.Thirdvalues))
                                       , rep(TP[4], length(FB.CR.Fourthvalues))
                                       , rep(TP[5], length(FB.CR.Fifthvalues)))))
CR.BL.WUF.ALL <- cbind(as.numeric(c(CR.BL.Firstvalues
                                     , CR.BL.Secondvalues
                                     , CR.BL.Thirdvalues
                                     , CR.BL.Fourthvalues
                                     , CR.BL.Fifthvalues))
                        , as.numeric(c(rep(TP[1], length(CR.BL.Firstvalues))
                                       , rep(TP[2], length(CR.BL.Secondvalues))
                                       , rep(TP[3], length(CR.BL.Thirdvalues))
                                       , rep(TP[4], length(CR.BL.Fourthvalues))
                                       , rep(TP[5], length(CR.BL.Fifthvalues)))))

# FB.BL.WUF.lm <- lm(FB.BL.WUF.ALL[,1] ~ log(FB.BL.WUF.ALL[,2]))
# summary(FB.BL.WUF.lm)
# 
# CR.BL.WUF.lm <- lm(CR.BL.WUF.ALL[,1] ~ log(CR.BL.WUF.ALL[,2]))
# summary(CR.BL.WUF.lm)
# 
# FB.CR.WUF.lm <- lm(FB.CR.WUF.ALL[,1] ~ log(FB.CR.WUF.ALL[,2]))
# summary(FB.CR.WUF.lm)

# Plot the distances

FB.BL.WUF.ALL.mean <- aggregate(FB.BL.WUF.ALL, by = list(FB.BL.WUF.ALL[,2]), mean)
FB.BL.WUF.ALL.sd <- aggregate(FB.BL.WUF.ALL, by = list(FB.BL.WUF.ALL[,2]), sd)

CR.BL.WUF.ALL.mean <- aggregate(CR.BL.WUF.ALL, by = list(CR.BL.WUF.ALL[,2]), mean)
CR.BL.WUF.ALL.sd <- aggregate(CR.BL.WUF.ALL, by = list(CR.BL.WUF.ALL[,2]), sd)

FB.CR.WUF.ALL.mean <- aggregate(FB.CR.WUF.ALL, by = list(FB.CR.WUF.ALL[,2]), mean)
FB.CR.WUF.ALL.sd <- aggregate(FB.CR.WUF.ALL, by = list(FB.CR.WUF.ALL[,2]), sd)


######## PLOT DISP ###############
# ylimits <- c(min(FB.BL.WUF.ALL[,1],FB.CR.WUF.ALL[,1],CR.BL.WUF.ALL[,1]), max(FB.BL.WUF.ALL[,1],FB.CR.WUF.ALL[,1],CR.BL.WUF.ALL[,1]))
ylimits <- c(0,0.4)

# xvalues <- log(FB.BL.WUF.ALL.mean[,1])
xvalues <- as.character(FB.BL.WUF.ALL.mean[,1])
pdf(paste0("BETAPLOTS_H/DispOverTime_H_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1, at = c(1,2,3,4,5), labels = c('20 min','1 h', '6 h', '12 h','4 d'))
points(FB.BL.WUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*0.99
       , x1 = c(1,2,3,4,5)*0.99
       , y0 = c(FB.BL.WUF.ALL.mean[,2] - FB.BL.WUF.ALL.sd[,2]/2)
       , y1 = c(FB.BL.WUF.ALL.mean[,2] + FB.BL.WUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.WUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1
       , x1 = c(1,2,3,4,5)*1
       , y0 = c(FB.CR.WUF.ALL.mean[,2] - FB.CR.WUF.ALL.sd[,2]/2)
       , y1 = c(FB.CR.WUF.ALL.mean[,2] + FB.CR.WUF.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.WUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1.01
       , x1 = c(1,2,3,4,5)*1.01
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
       , col = c("darkgreen","purple","grey")
       , lwd = 2)
dev.off()

####### PLOT MORPH ############
# NO 5760 WUF

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


MorphColours <- c("darkorchid4","dodgerblue","salmon") 

pdf(paste0("BETAPLOTS_H/NMDS_H_",metric,"_Morph.pdf"), pointsize = 14)
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
      , col = MorphColours[1])
lines(NMDS.WUF.BL[NMDS.WUF.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.FB[NMDS.WUF.FB.chull,]
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

# 5760 minutes
NMDS.WUF.5760.only <- NMDS.WUF.all$points[grep("(BL|CR|FB)-5760-", rownames(NMDS.WUF.all$points)),]

NMDS.WUF.5760.only.CR <- NMDS.WUF.5760.only[grep("CR", rownames(NMDS.WUF.5760.only)),]
NMDS.WUF.5760.only.CR.chull <- chull(NMDS.WUF.5760.only.CR)
NMDS.WUF.5760.only.CR.chull <- c(NMDS.WUF.5760.only.CR.chull, NMDS.WUF.5760.only.CR.chull[1])

NMDS.WUF.5760.only.BL <- NMDS.WUF.5760.only[grep("BL", rownames(NMDS.WUF.5760.only)),]
NMDS.WUF.5760.only.BL.chull <- chull(NMDS.WUF.5760.only.BL)
NMDS.WUF.5760.only.BL.chull <- c(NMDS.WUF.5760.only.BL.chull, NMDS.WUF.5760.only.BL.chull[1])

NMDS.WUF.5760.only.FB <- NMDS.WUF.5760.only[grep("FB", rownames(NMDS.WUF.5760.only)),]
NMDS.WUF.5760.only.FB.chull <- chull(NMDS.WUF.5760.only.FB)
NMDS.WUF.5760.only.FB.chull <- c(NMDS.WUF.5760.only.FB.chull, NMDS.WUF.5760.only.FB.chull[1])


# pdf(paste0("BETAPLOTS_H/NMDS_",metric,"_TimebyMorph.pdf"),pointsize = 14, width = 7, height = 4)
# par(mfrow= c(1,4), oma = c(4,4,4,6))
# par(mar = c(4,0,4,0))
# plot(NMDS.WUF.20.only
#      , main = "20 Minutes"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.WUF.20.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p = ", ANOVA.WUF.20.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.WUF.20.only.CR[NMDS.WUF.20.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.WUF.20.only.BL[NMDS.WUF.20.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.WUF.20.only.FB[NMDS.WUF.20.only.FB.chull,]
#       , col = MorphColours[3])
# par(mar = c(4,0,4,0))
# plot(NMDS.WUF.60.only
#      , main = "1 Hour"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.WUF.60.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p =  ",ANOVA.WUF.60.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.WUF.60.only.CR[NMDS.WUF.60.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.WUF.60.only.BL[NMDS.WUF.60.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.WUF.60.only.FB[NMDS.WUF.60.only.FB.chull,]
#       , col = MorphColours[3])
# par(mar = c(4,0,4,0))
# plot(NMDS.WUF.360.only
#      , main = "6 Hours"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.WUF.360.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p = ", ANOVA.WUF.360.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.WUF.360.only.CR[NMDS.WUF.360.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.WUF.360.only.BL[NMDS.WUF.360.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.WUF.360.only.FB[NMDS.WUF.360.only.FB.chull,]
#       , col = MorphColours[3])
# par(mar = c(4,0,4,0))
# plot(NMDS.WUF.720.only
#      , main = "12 Hours"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.WUF.720.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p = ", ANOVA.WUF.720.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.WUF.720.only.CR[NMDS.WUF.720.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.WUF.720.only.BL[NMDS.WUF.720.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.WUF.720.only.FB[NMDS.WUF.720.only.FB.chull,]
#       , col = MorphColours[3])
# #STOP
# dev.off()

############ COMBO DISP BETA ################
# Disp and beta through time combined
# xvalues <- log(FB.BL.WUF.ALL.mean[,1])
xvalues <- as.character(FB.BL.WUF.ALL.mean[,1])

# Change factor of individual plots
MF.WUF.20.only$Morph <- factor(MF.WUF.20.only$Morph, levels = c("CR","BL","FB"))
MF.WUF.60.only$Morph <- factor(MF.WUF.60.only$Morph, levels = c("CR","BL","FB"))
MF.WUF.360.only$Morph <- factor(MF.WUF.360.only$Morph, levels = c("CR","BL","FB"))
MF.WUF.720.only$Morph <- factor(MF.WUF.720.only$Morph, levels = c("CR","BL","FB"))
MF.WUF.5760.only$Morph <- factor(MF.WUF.5760.only$Morph, levels = c("CR","BL","FB"))


pdf(paste0("BETAPLOTS_H/COMBO_H_dispbeta_",metric,".pdf"), pointsize = 14, width = 10, height = 7)
par(fig = c(0,0.8,0.3,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = ''
     , ylab = "Distance (Unweighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1
     , at = c(1,2,3,4,5)
     , labels = c("20 min","1 h","6 h","12 h","4 d")
     , las = 2)
points(FB.BL.WUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*0.99
       , x1 = c(1,2,3,4,5)*0.99
       , y0 = c(FB.BL.WUF.ALL.mean[,2] - FB.BL.WUF.ALL.sd[,2])
       , y1 = c(FB.BL.WUF.ALL.mean[,2] + FB.BL.WUF.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.WUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1
       , x1 = c(1,2,3,4,5)*1
       , y0 = c(FB.CR.WUF.ALL.mean[,2] - FB.CR.WUF.ALL.sd[,2])
       , y1 = c(FB.CR.WUF.ALL.mean[,2] + FB.CR.WUF.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.WUF.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1.01
       , x1 = c(1,2,3,4,5)*1.01
       , y0 = c(CR.BL.WUF.ALL.mean[,2] - CR.BL.WUF.ALL.sd[,2])
       , y1 = c(CR.BL.WUF.ALL.mean[,2] + CR.BL.WUF.ALL.sd[,2])
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
       , legend = c("FB:CR", "FB:BL", "CR:BL")
       , lty = 1
       , col = c("darkgreen","purple","grey")
       , lwd = 2 
)
# EACH GETS 0.15 SPACE TOTAL;
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.049,0.1932,0,0.4), new = TRUE)
plot(NMDS.WUF.20.only
     # , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.WUF.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.20.only$aov.tab[1]$Df[1],",",ANOVA.WUF.20.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.20.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.20.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.20.only.CR[NMDS.WUF.20.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.20.only.BL[NMDS.WUF.20.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.20.only.FB[NMDS.WUF.20.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.2032,0.3474,0,0.4), new = TRUE)
plot(NMDS.WUF.60.only
     # , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.WUF.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.60.only$aov.tab[1]$Df[1],",",ANOVA.WUF.60.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.60.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.60.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.60.only.CR[NMDS.WUF.60.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.60.only.BL[NMDS.WUF.60.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.60.only.FB[NMDS.WUF.60.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.3574,0.5016,0,0.4), new = TRUE)
plot(NMDS.WUF.360.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.WUF.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.360.only$aov.tab[1]$Df[1],",",ANOVA.WUF.360.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.360.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.360.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.360.only.CR[NMDS.WUF.360.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.360.only.BL[NMDS.WUF.360.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.360.only.FB[NMDS.WUF.360.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.5116,0.6558,0,0.4), new = TRUE)
plot(NMDS.WUF.720.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.WUF.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.720.only$aov.tab[1]$Df[1],",",ANOVA.WUF.720.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.720.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.720.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.720.only.CR[NMDS.WUF.720.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.720.only.BL[NMDS.WUF.720.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.720.only.FB[NMDS.WUF.720.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.6658,0.81,0,0.4), new = TRUE)
plot(NMDS.WUF.5760.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.WUF.5760.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.5760.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.5760.only$aov.tab[1]$Df[1],",",ANOVA.WUF.5760.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.5760.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.5760.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.5760.only.CR[NMDS.WUF.5760.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.5760.only.BL[NMDS.WUF.5760.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.5760.only.FB[NMDS.WUF.5760.only.FB.chull,]
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

betadisp.H.WUF.FB <- betadisp.WUF.time$distances[grep("FB", names(betadisp.WUF.time$distances))]
betadisp.H.WUF.BL <- betadisp.WUF.time$distances[grep("BL", names(betadisp.WUF.time$distances))]
betadisp.H.WUF.CR <- betadisp.WUF.time$distances[grep("CR", names(betadisp.WUF.time$distances))]

MF.H.FB <- MF.incl5760[grep("FB", rownames(MF.incl5760)),]
MF.H.BL <- MF.incl5760[grep("BL", rownames(MF.incl5760)),]
MF.H.CR <- MF.incl5760[grep("CR", rownames(MF.incl5760)),]

betadisp.FB.H.agg <- aggregate(betadisp.H.WUF.FB, by = list(MF.H.FB$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.BL.H.agg <- aggregate(betadisp.H.WUF.BL, by = list(MF.H.BL$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.CR.H.agg <- aggregate(betadisp.H.WUF.CR, by = list(MF.H.CR$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )

betadisp.WUF.time.forstat <- cbind(betadisp.WUF.time$distances, MF.incl5760[unlist(lapply(names(betadisp.WUF.time$distances), function(x) {grep(paste0("^",x,"$"), rownames(MF.incl5760))})),])
colnames(betadisp.WUF.time.forstat)[1] <- c("Distance")
ANOVA.betadisp.WUF <- anova(lm(Distance ~ Time*Morph, data = betadisp.WUF.time.forstat))
capture.output(ANOVA.betadisp.WUF, file = "./BETAPLOTS_H/ANOVA.betadisp.WUF.txt")

xvalues <- c("20","60","360","720","5760")
ylimits <- c(0,0.35)
pdf(paste0("./BETAPLOTS_H/BetaDisp_H_",metric,"_eachmorph.pdf"),pointsize = 14)
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = 'Time'
     , ylab = "Distance (Weighted Unifrac)"
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
dm.WUF.5760<- dm.WUF.inclWater[grep("(CR|FB|BL)-5760", rownames(dm.WUF.inclWater)),grep("(CR|FB|BL)-5760", colnames(dm.WUF.inclWater))]

NMDS.WUF.5760 <- isoMDS(as.matrix(dm.WUF.5760), y = cmdscale(as.matrix(dm.WUF.5760), 2))

# NMDS.WUF.morphallTimePoints <- NMDS.WUF.morphAllTime$points[grep("(CR|BL|FB)-5760", rownames(NMDS.WUF.morphAllTime$points)),]
MF.5760 <- MF[sapply(rownames(NMDS.WUF.5760$points), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
}),]
MF.5760$Morph <- factor(MF.5760$Morph, levels = c("CR","BL","FB"))

NMDS.WUF.5760.CR <- NMDS.WUF.5760$points[grep("CR", rownames(NMDS.WUF.5760$points)),]
NMDS.WUF.5760.CR.chull <- chull(NMDS.WUF.5760.CR)
NMDS.WUF.5760.CR.chull <- c(NMDS.WUF.5760.CR.chull, NMDS.WUF.5760.CR.chull[1])

NMDS.WUF.5760.BL <- NMDS.WUF.5760$points[grep("BL", rownames(NMDS.WUF.5760$points)),]
NMDS.WUF.5760.BL.chull <- chull(NMDS.WUF.5760.BL)
NMDS.WUF.5760.BL.chull <- c(NMDS.WUF.5760.BL.chull, NMDS.WUF.5760.BL.chull[1])

NMDS.WUF.5760.FB <- NMDS.WUF.5760$points[grep("FB", rownames(NMDS.WUF.5760$points)),]
NMDS.WUF.5760.FB.chull <- chull(NMDS.WUF.5760.FB)
NMDS.WUF.5760.FB.chull <- c(NMDS.WUF.5760.FB.chull, NMDS.WUF.5760.FB.chull[1])

# ANOVA.WUF.5760 <- adonis(dm.WUF.5760 ~ Morph, data = MF.5760)
# # ANOSIM.WUF.5760 <- anosim(dat = dm.WUF.5760, grouping = MF.5760$Morph)
# capture.output(ANOVA.WUF.5760, file = paste0("BETAPLOTS_H/anova_",metric,"_5760only.txt"))

pdf(paste0("BETAPLOTS_H/NMDS_",metric,"_5760Only.pdf"), pointsize = 14)
par(fig = c(0,0.75,0,1))
plot(NMDS.WUF.5760$points
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.5760$Morph)]
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , sub = paste0("Stress: ", round(NMDS.WUF.5760$stress,2)/100)
     , cex = 1.5
     , main = "NMDS of community composition (4 days)")
lines(NMDS.WUF.5760.CR[NMDS.WUF.5760.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.5760.BL[NMDS.WUF.5760.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.5760.FB[NMDS.WUF.5760.FB.chull,]
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


# Combine individual stats and print out

allindividualTests <- c("ANOVA.WUF.20.only"
                        ,"ANOVA.WUF.60.only"
                        ,"ANOVA.WUF.360.only"
                        ,"ANOVA.WUF.720.only"
                        ,"ANOVA.WUF.5760.only"
)
for (i in allindividualTests) {
  # print(get(i))
  capture.output(get(i), file = paste0("./BETAPLOTS_H/individualtests/",i,".txt"))
}


######## --PM-- ###########
### NMDS #####
rownames(dm.WUF.P.morphonly) <- gsub(".","-", rownames(dm.WUF.P.morphonly), fixed = TRUE)
colnames(dm.WUF.P.morphonly) <- gsub(".","-", colnames(dm.WUF.P.morphonly), fixed = TRUE)

rownames(MF.P.morphkeep) <- gsub(".","-", rownames(MF.P.morphkeep), fixed = TRUE)
rownames(MF.P.inclWater) <- gsub(".","-", rownames(MF.P.inclWater), fixed = TRUE)



NMDS.WUF.P.morphonly <- isoMDS(as.matrix(dm.WUF.P.morphonly), y = cmdscale(as.matrix(dm.WUF.P.morphonly), 2))
NMDS.WUF.P.all <- isoMDS(as.matrix(dm.WUF.P.inclWater), y = cmdscale(as.matrix(dm.WUF.P.inclWater), 2))

###### STATS ##########
MF.P.morphkeep <- MF.P.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
MF.P.morphkeep$Morph <- factor(MF.P.morphkeep$Morph, levels = c('CR','BL','FB'))
MF.P.morphkeep$Time <- factor(MF.P.morphkeep$Time, levels = c('20','60','180','360','720','1440'))
MF.P.morphkeep$Type <- factor(MF.P.morphkeep$Type, levels = c('P','H'))


# dm.WUF.P.morphonly.H <- dm.WUF.P.morphonly[grep(".", rownames(dm.WUF.P.morphonly), fixed = TRUE),grep(".", colnames(dm.WUF.P.morphonly), fixed = TRUE)]
# dm.WUF.P.morphonly.P <- dm.WUF.P.morphonly[grep("-", rownames(dm.WUF.P.morphonly), fixed = TRUE),grep("-", colnames(dm.WUF.P.morphonly), fixed = TRUE)]

ANOVA.WUF.P.morphtime <- adonis(dm.WUF.P.morphonly ~ Time*Morph, data = MF.P.morphkeep, by = "margin")
capture.output(ANOVA.WUF.P.morphtime, file = paste0("BETAPLOTS_P/adonis_", metric,"_PM.txt"))
# ANOSIM.WUF.P.morphtime <- anosim(dm.WUF.P.morphonly, grouping = MF.P.morphkeep$Morph)

# Dispersion across time and between morphs
dist.WUF.P.morphonly <- as.dist(dm.WUF.P.morphonly)
betadisp.WUF.P.time <- betadisper(d = dist.WUF.P.morphonly, group = MF.P.morphkeep$Time)
betadisp.WUF.P.morph <- betadisper(d = dist.WUF.P.morphonly, group = MF.P.morphkeep$Morph)

capture.output(betadisp.WUF.P.time, file = paste0("BETAPLOTS_P/betadispTime_", metric, "_PM.txt"))
capture.output(betadisp.WUF.P.morph, file = paste0("BETAPLOTS_P/betadispMorph_", metric, "_PM.txt"))


# Grep values for each timepoint and morph, and then plot dispersion within each
TP <- levels(factor(MF.P.morphkeep$Time))
MorphTypes <- levels(factor(MF.P.morphkeep$Morph))

# First Timepoint
FirstTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[1],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.WUF.P.FirstTP <- dm.WUF.P.morphonly[FirstTPNames,FirstTPNames]
FB.BL.Firstvalues <- as.vector(as.dist(dm.WUF.P.FirstTP[grep("FB", rownames(dm.WUF.P.FirstTP)), grep("BL", colnames(dm.WUF.P.FirstTP))]))
FB.CR.Firstvalues <- as.vector(as.dist(dm.WUF.P.FirstTP[grep("FB", rownames(dm.WUF.P.FirstTP)), grep("CR", colnames(dm.WUF.P.FirstTP))]))
CR.BL.Firstvalues <- as.vector(as.dist(dm.WUF.P.FirstTP[grep("CR", rownames(dm.WUF.P.FirstTP)), grep("BL", colnames(dm.WUF.P.FirstTP))]))
# ANOVA
dm.WUF.P.20 <- dm.WUF.P.morphonly[grep("^20-", rownames(dm.WUF.P.morphonly)), grep("^20-", colnames(dm.WUF.P.morphonly))]
MF.P.WUF.P.20.only <- MF.P.morphkeep[grep("^20-", rownames(MF.P.morphkeep)),]
ANOVA.WUF.P.20.only <- adonis(dm.WUF.P.20 ~ Morph, data = MF.P.WUF.P.20.only, by = "margin")

# Second Timepoint
SecondTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[2],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.WUF.P.SecondTP <- dm.WUF.P.morphonly[SecondTPNames,SecondTPNames]
FB.BL.Secondvalues <- as.vector(as.dist(dm.WUF.P.SecondTP[grep("FB", rownames(dm.WUF.P.SecondTP)), grep("BL", colnames(dm.WUF.P.SecondTP))]))
FB.CR.Secondvalues <- as.vector(as.dist(dm.WUF.P.SecondTP[grep("FB", rownames(dm.WUF.P.SecondTP)), grep("CR", colnames(dm.WUF.P.SecondTP))]))
CR.BL.Secondvalues <- as.vector(as.dist(dm.WUF.P.SecondTP[grep("CR", rownames(dm.WUF.P.SecondTP)), grep("BL", colnames(dm.WUF.P.SecondTP))]))
# ANOVA
dm.WUF.P.60 <- dm.WUF.P.morphonly[grep("^60-", rownames(dm.WUF.P.morphonly)), grep("^60-", colnames(dm.WUF.P.morphonly))]
MF.P.WUF.P.60.only <- MF.P.morphkeep[grep("^60-", rownames(MF.P.morphkeep)),]
ANOVA.WUF.P.60.only <- adonis(dm.WUF.P.60 ~ Morph, data = MF.P.WUF.P.60.only, by = "margin")

# Third Timepoint
ThirdTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[3],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.WUF.P.ThirdTP <- dm.WUF.P.morphonly[ThirdTPNames,ThirdTPNames]
FB.BL.Thirdvalues <- as.vector(as.dist(dm.WUF.P.ThirdTP[grep("FB", rownames(dm.WUF.P.ThirdTP)), grep("BL", colnames(dm.WUF.P.ThirdTP))]))
FB.CR.Thirdvalues <- as.vector(as.dist(dm.WUF.P.ThirdTP[grep("FB", rownames(dm.WUF.P.ThirdTP)), grep("CR", colnames(dm.WUF.P.ThirdTP))]))
CR.BL.Thirdvalues <- as.vector(as.dist(dm.WUF.P.ThirdTP[grep("CR", rownames(dm.WUF.P.ThirdTP)), grep("BL", colnames(dm.WUF.P.ThirdTP))]))
# ANOVA
dm.WUF.P.180 <- dm.WUF.P.morphonly[grep("^180-", rownames(dm.WUF.P.morphonly)), grep("^180-", colnames(dm.WUF.P.morphonly))]
MF.P.WUF.P.180.only <- MF.P.morphkeep[grep("^180-", rownames(MF.P.morphkeep)),]
ANOVA.WUF.P.180.only <- adonis(dm.WUF.P.180 ~ Morph, data = MF.P.WUF.P.180.only, by = "margin")


# Fourth Timepoint
FourthTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[4],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.WUF.P.FourthTP <- dm.WUF.P.morphonly[FourthTPNames,FourthTPNames]
FB.BL.Fourthvalues <- as.vector(as.dist(dm.WUF.P.FourthTP[grep("FB", rownames(dm.WUF.P.FourthTP)), grep("BL", colnames(dm.WUF.P.FourthTP))]))
FB.CR.Fourthvalues <- as.vector(as.dist(dm.WUF.P.FourthTP[grep("FB", rownames(dm.WUF.P.FourthTP)), grep("CR", colnames(dm.WUF.P.FourthTP))]))
CR.BL.Fourthvalues <- as.vector(as.dist(dm.WUF.P.FourthTP[grep("CR", rownames(dm.WUF.P.FourthTP)), grep("BL", colnames(dm.WUF.P.FourthTP))]))
# ANOVA
dm.WUF.P.360 <- dm.WUF.P.morphonly[grep("^360-", rownames(dm.WUF.P.morphonly)), grep("^360-", colnames(dm.WUF.P.morphonly))]
MF.P.WUF.P.360.only <- MF.P.morphkeep[grep("^360-", rownames(MF.P.morphkeep)),]
ANOVA.WUF.P.360.only <- adonis(dm.WUF.P.360 ~ Morph, data = MF.P.WUF.P.360.only, by = "margin")

# Fifth Timepoint
FifthTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[5],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.WUF.P.FifthTP <- dm.WUF.P.morphonly[FifthTPNames,FifthTPNames]
FB.BL.Fifthvalues <- as.vector(as.dist(dm.WUF.P.FifthTP[grep("FB", rownames(dm.WUF.P.FifthTP)), grep("BL", colnames(dm.WUF.P.FifthTP))]))
FB.CR.Fifthvalues <- as.vector(as.dist(dm.WUF.P.FifthTP[grep("FB", rownames(dm.WUF.P.FifthTP)), grep("CR", colnames(dm.WUF.P.FifthTP))]))
CR.BL.Fifthvalues <- as.vector(as.dist(dm.WUF.P.FifthTP[grep("CR", rownames(dm.WUF.P.FifthTP)), grep("BL", colnames(dm.WUF.P.FifthTP))]))
# ANOVA
dm.WUF.P.720 <- dm.WUF.P.morphonly[grep("^720-", rownames(dm.WUF.P.morphonly)), grep("^720-", colnames(dm.WUF.P.morphonly))]
MF.P.WUF.P.720.only <- MF.P.morphkeep[grep("^720-", rownames(MF.P.morphkeep)),]
ANOVA.WUF.P.720.only <- adonis(dm.WUF.P.720 ~ Morph, data = MF.P.WUF.P.720.only, by = "margin")

# Sixth Timepoint
SixthTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[6],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.WUF.P.SixthTP <- dm.WUF.P.morphonly[SixthTPNames,SixthTPNames]
FB.BL.Sixthvalues <- as.vector(as.dist(dm.WUF.P.SixthTP[grep("FB", rownames(dm.WUF.P.SixthTP)), grep("BL", colnames(dm.WUF.P.SixthTP))]))
FB.CR.Sixthvalues <- as.vector(as.dist(dm.WUF.P.SixthTP[grep("FB", rownames(dm.WUF.P.SixthTP)), grep("CR", colnames(dm.WUF.P.SixthTP))]))
CR.BL.Sixthvalues <- as.vector(as.dist(dm.WUF.P.SixthTP[grep("CR", rownames(dm.WUF.P.SixthTP)), grep("BL", colnames(dm.WUF.P.SixthTP))]))
# ANOVA
dm.WUF.P.1440<- dm.WUF.P.morphonly[grep("^1440-", rownames(dm.WUF.P.morphonly)), grep("^1440-", colnames(dm.WUF.P.morphonly))]
MF.P.WUF.P.1440.only <- MF.P.morphkeep[grep("^1440-", rownames(MF.P.morphkeep)),]
ANOVA.WUF.P.1440.only <- adonis(dm.WUF.P.1440 ~ Morph, data = MF.P.WUF.P.1440.only, by = "margin")



# Combine into single tables

FB.BL.WUF.P.ALL <- cbind(as.numeric(c(FB.BL.Firstvalues
                                       , FB.BL.Secondvalues
                                       , FB.BL.Thirdvalues
                                       , FB.BL.Fourthvalues
                                       , FB.BL.Fifthvalues
                                       , FB.BL.Sixthvalues
))
, as.numeric(c(rep(TP[1], length(FB.BL.Firstvalues))
               , rep(TP[2], length(FB.BL.Secondvalues))
               , rep(TP[3], length(FB.BL.Thirdvalues))
               , rep(TP[4], length(FB.BL.Fourthvalues))
               , rep(TP[5], length(FB.BL.Fifthvalues))
               , rep(TP[6], length(FB.BL.Sixthvalues))
)))
FB.CR.WUF.P.ALL <- cbind(as.numeric(c(FB.CR.Firstvalues
                                       , FB.CR.Secondvalues
                                       , FB.CR.Thirdvalues
                                       , FB.CR.Fourthvalues
                                       , FB.CR.Fifthvalues
                                       , FB.CR.Sixthvalues
))
, as.numeric(c(rep(TP[1], length(FB.CR.Firstvalues))
               , rep(TP[2], length(FB.CR.Secondvalues))
               , rep(TP[3], length(FB.CR.Thirdvalues))
               , rep(TP[4], length(FB.CR.Fourthvalues))
               , rep(TP[5], length(FB.CR.Fifthvalues))
               , rep(TP[6], length(FB.CR.Sixthvalues))
)))
CR.BL.WUF.P.ALL <- cbind(as.numeric(c(CR.BL.Firstvalues
                                       , CR.BL.Secondvalues
                                       , CR.BL.Thirdvalues
                                       , CR.BL.Fourthvalues
                                       , CR.BL.Fifthvalues
                                       , CR.BL.Sixthvalues
))
, as.numeric(c(rep(TP[1], length(CR.BL.Firstvalues))
               , rep(TP[2], length(CR.BL.Secondvalues))
               , rep(TP[3], length(CR.BL.Thirdvalues))
               , rep(TP[4], length(CR.BL.Fourthvalues))
               , rep(TP[5], length(CR.BL.Fifthvalues))
               , rep(TP[6], length(CR.BL.Sixthvalues))
)))

# FB.BL.WUF.P.lm <- lm(FB.BL.WUF.P.ALL[,1] ~ log(FB.BL.WUF.P.ALL[,2]))
# summary(FB.BL.WUF.P.lm)
# 
# CR.BL.WUF.P.lm <- lm(CR.BL.WUF.P.ALL[,1] ~ log(CR.BL.WUF.P.ALL[,2]))
# summary(CR.BL.WUF.P.lm)
# 
# FB.CR.WUF.P.lm <- lm(FB.CR.WUF.P.ALL[,1] ~ log(FB.CR.WUF.P.ALL[,2]))
# summary(FB.CR.WUF.P.lm)

# Plot the distances

FB.BL.WUF.P.ALL.mean <- aggregate(FB.BL.WUF.P.ALL, by = list(FB.BL.WUF.P.ALL[,2]), mean)
FB.BL.WUF.P.ALL.sd <- aggregate(FB.BL.WUF.P.ALL, by = list(FB.BL.WUF.P.ALL[,2]), sd)

CR.BL.WUF.P.ALL.mean <- aggregate(CR.BL.WUF.P.ALL, by = list(CR.BL.WUF.P.ALL[,2]), mean)
CR.BL.WUF.P.ALL.sd <- aggregate(CR.BL.WUF.P.ALL, by = list(CR.BL.WUF.P.ALL[,2]), sd)

FB.CR.WUF.P.ALL.mean <- aggregate(FB.CR.WUF.P.ALL, by = list(FB.CR.WUF.P.ALL[,2]), mean)
FB.CR.WUF.P.ALL.sd <- aggregate(FB.CR.WUF.P.ALL, by = list(FB.CR.WUF.P.ALL[,2]), sd)

######## PLOT DISP ###############
# ylimits <- c(min(FB.BL.WUF.P.ALL[,1],FB.CR.WUF.P.ALL[,1],CR.BL.WUF.P.ALL[,1]), max(FB.BL.WUF.P.ALL[,1],FB.CR.WUF.P.ALL[,1],CR.BL.WUF.P.ALL[,1]))
ylimits <- c(0,0.65)

# xvalues <- log(FB.BL.WUF.P.ALL.mean[,1])
xvalues <- as.character(FB.BL.WUF.P.ALL.mean[,1])
pdf(paste0("BETAPLOTS_P/DispOverTime_P_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Weighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1, at = c(1,2,3,4,5,6), labels = c('20 min','1 h','3 h','6 h', '12 h', '24 h'))
points(FB.BL.WUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*0.99
       , x1 = c(1,2,3,4,5,6)*0.99
       , y0 = c(FB.BL.WUF.P.ALL.mean[,2] - FB.BL.WUF.P.ALL.sd[,2]/2)
       , y1 = c(FB.BL.WUF.P.ALL.mean[,2] + FB.BL.WUF.P.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.WUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*1
       , x1 = c(1,2,3,4,5,6)*1
       , y0 = c(FB.CR.WUF.P.ALL.mean[,2] - FB.CR.WUF.P.ALL.sd[,2]/2)
       , y1 = c(FB.CR.WUF.P.ALL.mean[,2] + FB.CR.WUF.P.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.WUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*1.01
       , x1 = c(1,2,3,4,5,6)*1.01
       , y0 = c(CR.BL.WUF.P.ALL.mean[,2] - CR.BL.WUF.P.ALL.sd[,2]/2)
       , y1 = c(CR.BL.WUF.P.ALL.mean[,2] + CR.BL.WUF.P.ALL.sd[,2]/2)
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
# NOCHLP WUF.P

# SWITCH AROUND AXIS
NMDS.WUF.P.morphonly$points[,1] <- -(NMDS.WUF.P.morphonly$points[,1])

# MAKE POLYGONS for plotting
NMDS.WUF.P.CR <- NMDS.WUF.P.morphonly$points[grep("CR", rownames(NMDS.WUF.P.morphonly$points)),]
NMDS.WUF.P.BL <- NMDS.WUF.P.morphonly$points[grep("BL", rownames(NMDS.WUF.P.morphonly$points)),]
NMDS.WUF.P.FB <- NMDS.WUF.P.morphonly$points[grep("FB", rownames(NMDS.WUF.P.morphonly$points)),]

NMDS.WUF.P.CR.chull <- chull(NMDS.WUF.P.CR)
NMDS.WUF.P.CR.chull <- c(NMDS.WUF.P.CR.chull, NMDS.WUF.P.CR.chull[1])

NMDS.WUF.P.BL.chull <- chull(NMDS.WUF.P.BL)
NMDS.WUF.P.BL.chull <- c(NMDS.WUF.P.BL.chull, NMDS.WUF.P.BL.chull[1])

NMDS.WUF.P.FB.chull <- chull(NMDS.WUF.P.FB)
NMDS.WUF.P.FB.chull <- c(NMDS.WUF.P.FB.chull, NMDS.WUF.P.FB.chull[1])


MorphColours <- c("darkorchid4","dodgerblue","salmon") 

pdf(paste0("BETAPLOTS_P/NMDS_P_",metric,"_Morph.pdf"), pointsize = 14)
par(fig = c(0,0.8,0,1))
plot(NMDS.WUF.P.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.morphkeep$Morph)]
     , sub = round(NMDS.WUF.P.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
)
lines(NMDS.WUF.P.CR[NMDS.WUF.P.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.BL[NMDS.WUF.P.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.FB[NMDS.WUF.P.FB.chull,]
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
plot(NMDS.WUF.P.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = TimeColours[factor(MF.P.morphkeep$Time)]
     , sub = round(NMDS.WUF.P.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
)
# lines(NMDS.WUF.P.CR[NMDS.WUF.P.CR.chull,]
#       , col = "red")
# lines(NMDS.WUF.P.BL[NMDS.WUF.P.BL.chull,]
#       , col = "magenta")
# lines(NMDS.WUF.P.FB[NMDS.WUF.P.FB.chull,]
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



# TIME-- ONE at a time


# 20 minutes
NMDS.WUF.P.20.only <- NMDS.WUF.P.morphonly$points[grep("^20-", rownames(NMDS.WUF.P.morphonly$points)),]

NMDS.WUF.P.20.only.CR <- NMDS.WUF.P.20.only[grep("CR", rownames(NMDS.WUF.P.20.only)),]
NMDS.WUF.P.20.only.CR.chull <- chull(NMDS.WUF.P.20.only.CR)
NMDS.WUF.P.20.only.CR.chull <- c(NMDS.WUF.P.20.only.CR.chull, NMDS.WUF.P.20.only.CR.chull[1])

NMDS.WUF.P.20.only.BL <- NMDS.WUF.P.20.only[grep("BL", rownames(NMDS.WUF.P.20.only)),]
NMDS.WUF.P.20.only.BL.chull <- chull(NMDS.WUF.P.20.only.BL)
NMDS.WUF.P.20.only.BL.chull <- c(NMDS.WUF.P.20.only.BL.chull, NMDS.WUF.P.20.only.BL.chull[1])

NMDS.WUF.P.20.only.FB <- NMDS.WUF.P.20.only[grep("FB", rownames(NMDS.WUF.P.20.only)),]
NMDS.WUF.P.20.only.FB.chull <- chull(NMDS.WUF.P.20.only.FB)
NMDS.WUF.P.20.only.FB.chull <- c(NMDS.WUF.P.20.only.FB.chull, NMDS.WUF.P.20.only.FB.chull[1])

# 60 minutes
NMDS.WUF.P.60.only <- NMDS.WUF.P.morphonly$points[grep("^60-", rownames(NMDS.WUF.P.morphonly$points)),]

NMDS.WUF.P.60.only.CR <- NMDS.WUF.P.60.only[grep("CR", rownames(NMDS.WUF.P.60.only)),]
NMDS.WUF.P.60.only.CR.chull <- chull(NMDS.WUF.P.60.only.CR)
NMDS.WUF.P.60.only.CR.chull <- c(NMDS.WUF.P.60.only.CR.chull, NMDS.WUF.P.60.only.CR.chull[1])

NMDS.WUF.P.60.only.BL <- NMDS.WUF.P.60.only[grep("BL", rownames(NMDS.WUF.P.60.only)),]
NMDS.WUF.P.60.only.BL.chull <- chull(NMDS.WUF.P.60.only.BL)
NMDS.WUF.P.60.only.BL.chull <- c(NMDS.WUF.P.60.only.BL.chull, NMDS.WUF.P.60.only.BL.chull[1])

NMDS.WUF.P.60.only.FB <- NMDS.WUF.P.60.only[grep("FB", rownames(NMDS.WUF.P.60.only)),]
NMDS.WUF.P.60.only.FB.chull <- chull(NMDS.WUF.P.60.only.FB)
NMDS.WUF.P.60.only.FB.chull <- c(NMDS.WUF.P.60.only.FB.chull, NMDS.WUF.P.60.only.FB.chull[1])

# 180 minutes
NMDS.WUF.P.180.only <- NMDS.WUF.P.morphonly$points[grep("^180-", rownames(NMDS.WUF.P.morphonly$points)),]

NMDS.WUF.P.180.only.CR <- NMDS.WUF.P.180.only[grep("CR", rownames(NMDS.WUF.P.180.only)),]
NMDS.WUF.P.180.only.CR.chull <- chull(NMDS.WUF.P.180.only.CR)
NMDS.WUF.P.180.only.CR.chull <- c(NMDS.WUF.P.180.only.CR.chull, NMDS.WUF.P.180.only.CR.chull[1])

NMDS.WUF.P.180.only.BL <- NMDS.WUF.P.180.only[grep("BL", rownames(NMDS.WUF.P.180.only)),]
NMDS.WUF.P.180.only.BL.chull <- chull(NMDS.WUF.P.180.only.BL)
NMDS.WUF.P.180.only.BL.chull <- c(NMDS.WUF.P.180.only.BL.chull, NMDS.WUF.P.180.only.BL.chull[1])

NMDS.WUF.P.180.only.FB <- NMDS.WUF.P.180.only[grep("FB", rownames(NMDS.WUF.P.180.only)),]
NMDS.WUF.P.180.only.FB.chull <- chull(NMDS.WUF.P.180.only.FB)
NMDS.WUF.P.180.only.FB.chull <- c(NMDS.WUF.P.180.only.FB.chull, NMDS.WUF.P.180.only.FB.chull[1])



# 360 minutes
NMDS.WUF.P.360.only <- NMDS.WUF.P.morphonly$points[grep("^360-", rownames(NMDS.WUF.P.morphonly$points)),]


NMDS.WUF.P.360.only.CR <- NMDS.WUF.P.360.only[grep("CR", rownames(NMDS.WUF.P.360.only)),]
NMDS.WUF.P.360.only.CR.chull <- chull(NMDS.WUF.P.360.only.CR)
NMDS.WUF.P.360.only.CR.chull <- c(NMDS.WUF.P.360.only.CR.chull, NMDS.WUF.P.360.only.CR.chull[1])

NMDS.WUF.P.360.only.BL <- NMDS.WUF.P.360.only[grep("BL", rownames(NMDS.WUF.P.360.only)),]
NMDS.WUF.P.360.only.BL.chull <- chull(NMDS.WUF.P.360.only.BL)
NMDS.WUF.P.360.only.BL.chull <- c(NMDS.WUF.P.360.only.BL.chull, NMDS.WUF.P.360.only.BL.chull[1])

NMDS.WUF.P.360.only.FB <- NMDS.WUF.P.360.only[grep("FB", rownames(NMDS.WUF.P.360.only)),]
NMDS.WUF.P.360.only.FB.chull <- chull(NMDS.WUF.P.360.only.FB)
NMDS.WUF.P.360.only.FB.chull <- c(NMDS.WUF.P.360.only.FB.chull, NMDS.WUF.P.360.only.FB.chull[1])


# 720 minutes
NMDS.WUF.P.720.only <- NMDS.WUF.P.morphonly$points[grep("^720-", rownames(NMDS.WUF.P.morphonly$points)),]


NMDS.WUF.P.720.only.CR <- NMDS.WUF.P.720.only[grep("CR", rownames(NMDS.WUF.P.720.only)),]
NMDS.WUF.P.720.only.CR.chull <- chull(NMDS.WUF.P.720.only.CR)
NMDS.WUF.P.720.only.CR.chull <- c(NMDS.WUF.P.720.only.CR.chull, NMDS.WUF.P.720.only.CR.chull[1])

NMDS.WUF.P.720.only.BL <- NMDS.WUF.P.720.only[grep("BL", rownames(NMDS.WUF.P.720.only)),]
NMDS.WUF.P.720.only.BL.chull <- chull(NMDS.WUF.P.720.only.BL)
NMDS.WUF.P.720.only.BL.chull <- c(NMDS.WUF.P.720.only.BL.chull, NMDS.WUF.P.720.only.BL.chull[1])

NMDS.WUF.P.720.only.FB <- NMDS.WUF.P.720.only[grep("FB", rownames(NMDS.WUF.P.720.only)),]
NMDS.WUF.P.720.only.FB.chull <- chull(NMDS.WUF.P.720.only.FB)
NMDS.WUF.P.720.only.FB.chull <- c(NMDS.WUF.P.720.only.FB.chull, NMDS.WUF.P.720.only.FB.chull[1])


# 1440 minutes
NMDS.WUF.P.1440.only <- NMDS.WUF.P.morphonly$points[grep("^1440-", rownames(NMDS.WUF.P.morphonly$points)),]


NMDS.WUF.P.1440.only.CR <- NMDS.WUF.P.1440.only[grep("CR", rownames(NMDS.WUF.P.1440.only)),]
NMDS.WUF.P.1440.only.CR.chull <- chull(NMDS.WUF.P.1440.only.CR)
NMDS.WUF.P.1440.only.CR.chull <- c(NMDS.WUF.P.1440.only.CR.chull, NMDS.WUF.P.1440.only.CR.chull[1])

NMDS.WUF.P.1440.only.BL <- NMDS.WUF.P.1440.only[grep("BL", rownames(NMDS.WUF.P.1440.only)),]
NMDS.WUF.P.1440.only.BL.chull <- chull(NMDS.WUF.P.1440.only.BL)
NMDS.WUF.P.1440.only.BL.chull <- c(NMDS.WUF.P.1440.only.BL.chull, NMDS.WUF.P.1440.only.BL.chull[1])

NMDS.WUF.P.1440.only.FB <- NMDS.WUF.P.1440.only[grep("FB", rownames(NMDS.WUF.P.1440.only)),]
NMDS.WUF.P.1440.only.FB.chull <- chull(NMDS.WUF.P.1440.only.FB)
NMDS.WUF.P.1440.only.FB.chull <- c(NMDS.WUF.P.1440.only.FB.chull, NMDS.WUF.P.1440.only.FB.chull[1])


pdf(paste0("BETAPLOTS_P/NMDS_",metric,"_TimebyMorph.pdf"),pointsize = 14, width = 7, height = 4)
par(mfrow= c(1,4), oma = c(4,4,4,6))
par(mar = c(4,0,4,0))
plot(NMDS.WUF.P.20.only
     , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.P.20.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.P.20.only.CR[NMDS.WUF.P.20.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.20.only.BL[NMDS.WUF.P.20.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.20.only.FB[NMDS.WUF.P.20.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.WUF.P.60.only
     , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.WUF.P.60.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.P.60.only.CR[NMDS.WUF.P.60.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.60.only.BL[NMDS.WUF.P.60.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.60.only.FB[NMDS.WUF.P.60.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.WUF.P.180.only
     , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.180.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.WUF.P.180.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.P.180.only.CR[NMDS.WUF.P.180.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.180.only.BL[NMDS.WUF.P.180.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.180.only.FB[NMDS.WUF.P.180.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.WUF.P.360.only
     , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.P.360.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.P.360.only.CR[NMDS.WUF.P.360.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.360.only.BL[NMDS.WUF.P.360.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.360.only.FB[NMDS.WUF.P.360.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.WUF.P.720.only
     , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.P.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.P.720.only.CR[NMDS.WUF.P.720.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.720.only.BL[NMDS.WUF.P.720.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.720.only.FB[NMDS.WUF.P.720.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.WUF.P.1440.only
     , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.1440.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.WUF.P.1440.only.CR[NMDS.WUF.P.1440.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.1440.only.BL[NMDS.WUF.P.1440.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.1440.only.FB[NMDS.WUF.P.1440.only.FB.chull,]
      , col = MorphColours[3])
#STOP
dev.off()

############ COMBO DISP BETA ################
# Disp and beta through time combined

# xvalues <- log(FB.BL.WUF.P.ALL.mean[,1])
xvalues <- as.character(FB.BL.WUF.P.ALL.mean[,1])

pdf(paste0("BETAPLOTS_P/COMBO_P_dispbeta_",metric,".pdf"), pointsize = 14, width = 14, height = 8)
par(fig = c(0,0.8,0.23,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = ''
     , ylab = "Distance (Weighted Unifrac)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1
     , at = c(1,2,3,4,5,6)
     , labels = c("20 min","1 h","3 h","6 h","12 h","24 h")
     , las = 2)
points(FB.BL.WUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*0.99
       , x1 = c(1,2,3,4,5,6)*0.99
       , y0 = c(FB.BL.WUF.P.ALL.mean[,2] - FB.BL.WUF.P.ALL.sd[,2])
       , y1 = c(FB.BL.WUF.P.ALL.mean[,2] + FB.BL.WUF.P.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.WUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*1
       , x1 = c(1,2,3,4,5,6)*1
       , y0 = c(FB.CR.WUF.P.ALL.mean[,2] - FB.CR.WUF.P.ALL.sd[,2])
       , y1 = c(FB.CR.WUF.P.ALL.mean[,2] + FB.CR.WUF.P.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.WUF.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1.01
       , x1 = c(1,2,3,4)*1.01
       , y0 = c(CR.BL.WUF.P.ALL.mean[,2] - CR.BL.WUF.P.ALL.sd[,2])
       , y1 = c(CR.BL.WUF.P.ALL.mean[,2] + CR.BL.WUF.P.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
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
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.05,0.158333,0,0.3), new = TRUE)
plot(NMDS.WUF.P.20.only
     # , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.P.20.only$aov.tab[1]$Df[1],",",ANOVA.WUF.P.20.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.P.20.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.P.20.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.P.20.only.CR[NMDS.WUF.P.20.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.20.only.BL[NMDS.WUF.P.20.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.20.only.FB[NMDS.WUF.P.20.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.178333,0.28666,0,0.3), new = TRUE)
plot(NMDS.WUF.P.60.only
     # , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.P.60.only$aov.tab[1]$Df[1],",",ANOVA.WUF.P.60.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.P.60.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.P.60.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.P.60.only.CR[NMDS.WUF.P.60.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.60.only.BL[NMDS.WUF.P.60.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.60.only.FB[NMDS.WUF.P.60.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.30666,0.414999,0,0.3), new = TRUE)
plot(NMDS.WUF.P.180.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.180.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.P.180.only$aov.tab[1]$Df[1],",",ANOVA.WUF.P.180.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.P.180.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.P.180.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.P.180.only.CR[NMDS.WUF.P.180.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.180.only.BL[NMDS.WUF.P.180.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.180.only.FB[NMDS.WUF.P.180.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.434999,0.543333,0,0.3), new = TRUE)
plot(NMDS.WUF.P.360.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.P.360.only$aov.tab[1]$Df[1],",",ANOVA.WUF.P.360.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.P.360.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.P.360.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.P.360.only.CR[NMDS.WUF.P.360.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.360.only.BL[NMDS.WUF.P.360.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.360.only.FB[NMDS.WUF.P.360.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.565555,0.671666,0,0.3), new = TRUE)
plot(NMDS.WUF.P.720.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.P.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.P.720.only$aov.tab[1]$Df[1],",",ANOVA.WUF.P.720.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.P.720.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.P.720.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.P.720.only.CR[NMDS.WUF.P.720.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.720.only.BL[NMDS.WUF.P.720.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.720.only.FB[NMDS.WUF.P.720.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.691666,0.79999,0,0.3), new = TRUE)
plot(NMDS.WUF.P.1440.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.WUF.P.1440.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.WUF.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.WUF.P.1440.only$aov.tab[1]$Df[1],",",ANOVA.WUF.P.1440.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.WUF.P.1440.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.WUF.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.WUF.P.1440.only.CR[NMDS.WUF.P.1440.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.WUF.P.1440.only.BL[NMDS.WUF.P.1440.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.WUF.P.1440.only.FB[NMDS.WUF.P.1440.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(2,0,2,0), fig = c(0.775,1,0,0.3), new = TRUE)
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

allindividualTests <- c("ANOVA.WUF.P.20.only"
                        ,"ANOVA.WUF.P.60.only"
                        ,"ANOVA.WUF.P.180.only"
                        ,"ANOVA.WUF.P.360.only"
                        ,"ANOVA.WUF.P.720.only"
                        ,"ANOVA.WUF.P.1440.only"
                        
)
for (i in allindividualTests) {
  # print(get(i))
  capture.output(get(i), file = paste0("./BETAPLOTS_P/individualtests/",i,".txt"))
}

########### BETADISP#############
betadisp.P.WUF.FB <- betadisp.WUF.P.time$distances[grep("FB", names(betadisp.WUF.P.time$distances))]
betadisp.P.WUF.BL <- betadisp.WUF.P.time$distances[grep("BL", names(betadisp.WUF.P.time$distances))]
betadisp.P.WUF.CR <- betadisp.WUF.P.time$distances[grep("CR", names(betadisp.WUF.P.time$distances))]

MF.P.FB <- MF.P.morphkeep[grep("FB", rownames(MF.P.morphkeep)),]
MF.P.BL <- MF.P.morphkeep[grep("BL", rownames(MF.P.morphkeep)),]
MF.P.CR <- MF.P.morphkeep[grep("CR", rownames(MF.P.morphkeep)),]

betadisp.FB.P.agg <- aggregate(betadisp.P.WUF.FB, by = list(MF.P.FB$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.BL.P.agg <- aggregate(betadisp.P.WUF.BL, by = list(MF.P.BL$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.CR.P.agg <- aggregate(betadisp.P.WUF.CR, by = list(MF.P.CR$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )

betadisp.WUF.time.forstat <- cbind(betadisp.WUF.time$distances, MF.morphkeep[unlist(lapply(names(betadisp.WUF.time$distances), function(x) {grep(paste0("^",x,"$"), rownames(MF.morphkeep))})),])
colnames(betadisp.WUF.time.forstat)[1] <- c("Distance")
ANOVA.betadisp.WUF <- anova(lm(Distance ~ Time*Morph, data = betadisp.WUF.time.forstat))
capture.output(ANOVA.betadisp.WUF, file = "./BETAPLOTS_P/ANOVA.betadisp.WUF.txt")

ylimits <- c(0,0.45)
pdf(paste0("./BETAPLOTS_P/BetaDisp_P_",metric,"_eachmorph.pdf"),pointsize = 14)
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = 'Time'
     , ylab = "Distance (Weighted Unifrac)"
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
anova.betadisp.WUF.morph <- anova(betadisp.WUF.morph)
anova.betadisp.WUF.time <- anova(betadisp.WUF.time)
anova.betadisp.P.WUF.morph <- anova(betadisp.WUF.P.morph)
anova.betadisp.P.WUF.time <- anova(betadisp.WUF.P.time)

# This is PERMADISP-- don't need a table, will quote in text
capture.output(anova.betadisp.WUF.morph, file = "BETAPLOTS_H/anova.betadisp.WUF.morph.txt")
capture.output(anova.betadisp.WUF.time, file = "BETAPLOTS_H/anova.betadisp.WUF.time.txt")
capture.output(anova.betadisp.P.WUF.morph, file = "BETAPLOTS_P/anova.betadisp.WUF.P.morph.txt")
capture.output(anova.betadisp.P.WUF.time, file = "BETAPLOTS_P/anova.betadisp.WUF.P.time.txt")

### NOW DO EACH TIME POINT ##
## This is a table with both Hakai and PM for the supplementary figures
permdisp.morphology.across.time <- matrix(ncol = 8, nrow = 2)
rownames(permdisp.morphology.across.time) <- c("P","H")
colnames(permdisp.morphology.across.time) <- c("20","60","180","360","720","1440","5760", "Overall")
for (t in c("20","60","180","360","720","1440")) {
  assign(paste0("betadisp.P.",t), betadisper(dist(get(paste0("dm.WUF.P.",t))), group = get(paste0("MF.P.WUF.P.",t,".only"))$Morph))
  assign(paste0("anova.betadisp.P.",t), anova(get(paste0("betadisp.P.",t))))
  ptemp <- get(paste0("anova.betadisp.P.",t))$`Pr(>F)`[1]
  ftemp <- get(paste0("anova.betadisp.P.",t))$`F value`[1]
  dftemp <- get(paste0("anova.betadisp.P.",t))$Df[1]
  
  toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
  permdisp.morphology.across.time["P", paste0(t)] <- toPaste
}
for (t in c("20","60","360","720","5760")) {
  assign(paste0("betadisp.",t), betadisper(dist(get(paste0("dm.WUF.",t))), group = get(paste0("MF.WUF.",t,".only"))$Morph))
  assign(paste0("anova.betadisp.",t), anova(get(paste0("betadisp.",t))))
  ptemp <- get(paste0("anova.betadisp.",t))$`Pr(>F)`[1]
  ftemp <- get(paste0("anova.betadisp.",t))$`F value`[1]
  dftemp <- get(paste0("anova.betadisp.",t))$Df[1]
  
  toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
  permdisp.morphology.across.time["H", paste0(t)] <- toPaste
  
}
# Get overall
ptemp <- anova.betadisp.P.WUF.morph$`Pr(>F)`[1]
ftemp <- anova.betadisp.P.WUF.morph$`F value`[1]
dftemp <- anova.betadisp.P.WUF.morph$Df[1]
toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
permdisp.morphology.across.time["P","Overall"] <- toPaste

ptemp <- anova.betadisp.WUF.morph$`Pr(>F)`[1]
ftemp <- anova.betadisp.WUF.morph$`F value`[1]
dftemp <- anova.betadisp.WUF.morph$Df[1]
toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
permdisp.morphology.across.time["H","Overall"] <- toPaste

colnames(permdisp.morphology.across.time) <- c("20 minutes","1 hour","3 hours","6 hours","12 hours","1 day","4 days", "Overall")

tempfile1 <- permdisp.morphology.across.time
for (r in 1:nrow(tempfile1)) {
  for (c in 1:ncol(tempfile1)) {
    if (is.na(tempfile1[r,c])) {
      permdisp.morphology.across.time[r,c] <- "-"
    }}}


# Make double header
permdisp.morphology.across.time.WUF <- permdisp.morphology.across.time
rownames(permdisp.morphology.across.time.WUF) <- c("Reed Point","Hakai")

capture.output(xtable(permdisp.morphology.across.time.WUF, digits = NULL), file = paste0("BETAPLOTS_LATEX/permdisp.morph.across.time.",metric,".txt"))

### MAKE BETA DIV TABLES-- metrics separately but H and P together; extras will go in supp

anova.morphology.across.time <- matrix(ncol = 8, nrow = 2)
colnames(anova.morphology.across.time) <- c("20","60","180","360","720","1440","5760", "Overall")
rownames(anova.morphology.across.time) <- c("P","H")
for (t in c("20","60","180","360","720","1440")) {
  ptemp <- get(paste0("ANOVA.WUF.P.",t,".only"))$aov.tab$`Pr(>F)`[1]
  rtemp <- get(paste0("ANOVA.WUF.P.",t,".only"))$aov.tab$`R2`[1]
  dftemp <- get(paste0("ANOVA.WUF.P.",t,".only"))$aov.tab$`Df`[1]
  toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
  
  anova.morphology.across.time["P",paste0(t)] <- toPaste
}
for (t in c("20","60","360","720","5760")) {
  ptemp <- get(paste0("ANOVA.WUF.",t,".only"))$aov.tab$`Pr(>F)`[1]
  rtemp <- get(paste0("ANOVA.WUF.",t,".only"))$aov.tab$`R2`[1]
  dftemp <- get(paste0("ANOVA.WUF.",t,".only"))$aov.tab$`Df`[1]
  toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
  
  anova.morphology.across.time["H",paste0(t)] <- toPaste
}
# Do overall P
ptemp <- ANOVA.WUF.P.morphtime$aov.tab$`Pr(>F)`[2]
rtemp <- ANOVA.WUF.P.morphtime$aov.tab$`R2`[2]
dftemp <- ANOVA.WUF.P.morphtime$aov.tab$`Df`[2]
toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
anova.morphology.across.time["P","Overall"] <- toPaste

# Do overall H
ptemp <- ANOVA.WUF.morphtime$aov.tab$`Pr(>F)`[2]
rtemp <- ANOVA.WUF.morphtime$aov.tab$`R2`[2]
dftemp <- ANOVA.WUF.morphtime$aov.tab$`Df`[2]
toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
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
anova.morphology.across.time.WUF <- anova.morphology.across.time
rownames(anova.morphology.across.time.WUF) <- c("Reed Point","Hakai")

capture.output(xtable(anova.morphology.across.time.WUF, digits = NULL), file = paste0("BETAPLOTS_LATEX/anova.morph.across.time.",metric,".txt"))

### DO FB/CR/BL TEST FOR EACH
# Make table
pairwiseAdonis.all <- matrix(ncol = 7,nrow = 6)
colnames(pairwiseAdonis.all) <- c("20","60","180","360","720","1440", "5760")
rownames(pairwiseAdonis.all) <- c("pFB:BL","pFB:CR","pBL:CR","hFB:BL","hFB:CR","hBL:CR" )
listMorphs <- c("FB","BL","CR")
for (m in 1:(length(listMorphs)-1)) {
  for (n in (m+1):length(listMorphs)) {
    for (t in c("20","60","180","360","720","1440")) {
      tempMF <- get(paste0("MF.P.WUF.P.",t,".only"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), get(paste0("MF.P.WUF.P.",t,".only"))$Morph),]
      tempDM <- get(paste0("dm.WUF.P.",t))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.WUF.P.",t)))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.WUF.P.",t))))]
      tempAdonis <- adonis(tempDM ~ Morph, data = tempMF, by = "marginal")
      toPaste <- paste0( tempAdonis$aov.tab$`Pr(>F)`[1]
                         ," (R^2 = ", round(tempAdonis$aov.tab$R2[1],digits = 2)
                         ,", Df = ", tempAdonis$aov.tab$Df[1] 
                         , ")")
      pairwiseAdonis.all[paste0("p",listMorphs[m],":",listMorphs[n]),paste0(t)] <- toPaste
    }
    for (t in c("20","60","360","720","5760")) {
      tempMF <- get(paste0("MF.WUF.",t,".only"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), get(paste0("MF.WUF.",t,".only"))$Morph),]
      tempDM <- get(paste0("dm.WUF.",t))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.WUF.",t)))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.WUF.",t))))]
      tempAdonis <- adonis(tempDM ~ Morph, data = tempMF, by = "marginal")
      toPaste <- paste0( tempAdonis$aov.tab$`Pr(>F)`[1]
                         ," (R^2 = ", round(tempAdonis$aov.tab$R2[1],digits = 2)
                         ,", Df = ", tempAdonis$aov.tab$Df[1] 
                         , ")")
      pairwiseAdonis.all[paste0("h",listMorphs[m],":",listMorphs[n]),paste0(t)] <- toPaste
    }
  }
}

# Get rid of NAs
for (r in 1:nrow(pairwiseAdonis.all)) {
  for (c in 1:ncol(pairwiseAdonis.all)) {
    if (is.na(pairwiseAdonis.all[r,c])) {
      pairwiseAdonis.all[r,c] <- "-"
    }
  }
}

# Change Rownames
pairwiseAdonis.all.WUF <- cbind(c("FB:BL","FB:CR","BL:CR","FB:BL","FB:CR","BL:CR"), pairwiseAdonis.all)
rownames(pairwiseAdonis.all.WUF) <- c("Reed Point", " ","  ","Hakai","   ","    ")

capture.output(xtable(pairwiseAdonis.all.WUF), file = paste0("BETAPLOTS_LATEX/pairwiseAdonis.",metric,".txt"))

####### *********BC********* ############# 
metric <- "BC"

######## --HAKAI-- ###########
### NMDS #####

NMDS.BC.morphonly <- isoMDS(as.matrix(dm.BC.morphonly), y = cmdscale(as.matrix(dm.BC.morphonly), 2))
NMDS.BC.all <- isoMDS(as.matrix(dm.BC.inclWater), y = cmdscale(as.matrix(dm.BC.inclWater), 2))

###### STATS ##########
MF.morphkeep <- MF.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
MF.morphkeep$Morph <- factor(MF.morphkeep$Morph, levels = c('CR','BL','FB'))
MF.morphkeep$Time <- factor(MF.morphkeep$Time, levels = c('20','60','180','360','720','1440','5760'))
MF.morphkeep$Type <- factor(MF.morphkeep$Type, levels = c('P','H'))


# dm.BC.morphonly.H <- dm.BC.morphonly[grep(".", rownames(dm.BC.morphonly), fixed = TRUE),grep(".", colnames(dm.BC.morphonly), fixed = TRUE)]
# dm.BC.morphonly.P <- dm.BC.morphonly[grep("-", rownames(dm.BC.morphonly), fixed = TRUE),grep("-", colnames(dm.BC.morphonly), fixed = TRUE)]

ANOVA.BC.morphtime <- adonis(dm.BC.morphonly ~ Time*Morph, data = MF.morphkeep, by = "margin")
capture.output(ANOVA.BC.morphtime, file = paste0("BETAPLOTS_H/adonis_", metric,"_Hakai.txt"))
# ANOSIM.BC.morphtime <- anosim(dm.BC.morphonly, grouping = MF.morphkeep$Morph)

# Dispersion across time and between morphs
dist.BC.morphonly <- as.dist(dm.BC.inclWater[-grep("W", rownames(dm.BC.inclWater)), -grep("W", colnames(dm.BC.inclWater))])
MF.incl5760 <- MF.inclWater[-grep("W", rownames(MF.inclWater)),]
betadisp.BC.time <- betadisper(d = dist.BC.morphonly, group = MF.incl5760$Time)
betadisp.BC.morph <- betadisper(d = dist.BC.morphonly, group = MF.incl5760$Morph)

capture.output(betadisp.BC.time, file = paste0("BETAPLOTS_H/betadispTime_", metric, "_Hakai.txt"))
capture.output(betadisp.BC.morph, file = paste0("BETAPLOTS_H/betadispMorph_", metric, "_Hakai.txt"))


# Grep values for each timepoint and morph, and then plot dispersion within each
TP <- levels(factor(MF.inclWater$Time))
MorphTypes <- levels(factor(MF.morphkeep$Morph))
# First Timepoint
FirstTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[1],"$"), MF.inclWater$Time)]
toremove <- grep("W", FirstTPNames)
if (length(toremove) > 0) {
  FirstTPNames <- FirstTPNames[-toremove]
}
# Between morphs
dm.BC.FirstTP <- dm.BC.inclWater[FirstTPNames,FirstTPNames]
FB.BL.Firstvalues <- as.vector(as.dist(dm.BC.FirstTP[grep("FB", rownames(dm.BC.FirstTP)), grep("BL", colnames(dm.BC.FirstTP))]))
FB.CR.Firstvalues <- as.vector(as.dist(dm.BC.FirstTP[grep("FB", rownames(dm.BC.FirstTP)), grep("CR", colnames(dm.BC.FirstTP))]))
CR.BL.Firstvalues <- as.vector(as.dist(dm.BC.FirstTP[grep("CR", rownames(dm.BC.FirstTP)), grep("BL", colnames(dm.BC.FirstTP))]))
# ANOVA
dm.BC.20 <- dm.BC.inclWater[grep("(CR|BL|FB)-20-", rownames(dm.BC.inclWater)), grep("-20-", colnames(dm.BC.inclWater))]
MF.BC.20.only <- MF.inclWater[grep("(CR|BL|FB)-20-", rownames(MF.inclWater)),]
ANOVA.BC.20.only <- adonis(dm.BC.20 ~ Morph, data = MF.BC.20.only, by = "margin")

# Second Timepoine
SecondTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[2],"$"), MF.inclWater$Time)]
toremove <- grep("W", SecondTPNames)
if (length(toremove) > 0) {
  SecondTPNames <- SecondTPNames[-toremove]
}
# Between morphs
dm.BC.SecondTP <- dm.BC.inclWater[SecondTPNames,SecondTPNames]
FB.BL.Secondvalues <- as.vector(as.dist(dm.BC.SecondTP[grep("FB", rownames(dm.BC.SecondTP)), grep("BL", colnames(dm.BC.SecondTP))]))
FB.CR.Secondvalues <- as.vector(as.dist(dm.BC.SecondTP[grep("FB", rownames(dm.BC.SecondTP)), grep("CR", colnames(dm.BC.SecondTP))]))
CR.BL.Secondvalues <- as.vector(as.dist(dm.BC.SecondTP[grep("CR", rownames(dm.BC.SecondTP)), grep("BL", colnames(dm.BC.SecondTP))]))
# ANOVA
dm.BC.60 <- dm.BC.inclWater[grep("(BL|CR|FB)-60-", rownames(dm.BC.inclWater)), grep("-60-", colnames(dm.BC.morphonly))]
MF.BC.60.only <- MF.inclWater[grep("(BL|CR|FB)-60-", rownames(MF.inclWater)),]
ANOVA.BC.60.only <- adonis(dm.BC.60 ~ Morph, data = MF.BC.60.only, by = "margin")

# Third Timepoint
ThirdTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[3],"$"), MF.inclWater$Time)]
toremove <- grep("W", ThirdTPNames)
if (length(toremove) > 0) {
  ThirdTPNames <- ThirdTPNames[-toremove]
}
# Between morphs
dm.BC.ThirdTP <- dm.BC.inclWater[ThirdTPNames,ThirdTPNames]
FB.BL.Thirdvalues <- as.vector(as.dist(dm.BC.ThirdTP[grep("FB", rownames(dm.BC.ThirdTP)), grep("BL", colnames(dm.BC.ThirdTP))]))
FB.CR.Thirdvalues <- as.vector(as.dist(dm.BC.ThirdTP[grep("FB", rownames(dm.BC.ThirdTP)), grep("CR", colnames(dm.BC.ThirdTP))]))
CR.BL.Thirdvalues <- as.vector(as.dist(dm.BC.ThirdTP[grep("CR", rownames(dm.BC.ThirdTP)), grep("BL", colnames(dm.BC.ThirdTP))]))
# ANOVA
dm.BC.360 <- dm.BC.inclWater[grep("(BL|CR|FB)-360-", rownames(dm.BC.inclWater)), grep("-360-", colnames(dm.BC.morphonly))]
MF.BC.360.only <- MF.inclWater[grep("(BL|CR|FB)-360-", rownames(MF.inclWater)),]
ANOVA.BC.360.only <- adonis(dm.BC.360 ~ Morph, data = MF.BC.360.only, by = "margin")


# Fourth Timepoint
FourthTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[4],"$"), MF.inclWater$Time)]
toremove <- grep("W", FourthTPNames)
if (length(toremove) > 0) {
  FourthTPNames <- FourthTPNames[-toremove]
}
# Between morphs
dm.BC.FourthTP <- dm.BC.inclWater[FourthTPNames,FourthTPNames]
FB.BL.Fourthvalues <- as.vector(as.dist(dm.BC.FourthTP[grep("FB", rownames(dm.BC.FourthTP)), grep("BL", colnames(dm.BC.FourthTP))]))
FB.CR.Fourthvalues <- as.vector(as.dist(dm.BC.FourthTP[grep("FB", rownames(dm.BC.FourthTP)), grep("CR", colnames(dm.BC.FourthTP))]))
CR.BL.Fourthvalues <- as.vector(as.dist(dm.BC.FourthTP[grep("CR", rownames(dm.BC.FourthTP)), grep("BL", colnames(dm.BC.FourthTP))]))
# ANOVA
dm.BC.720 <- dm.BC.inclWater[grep("(BL|CR|FB)-720-", rownames(dm.BC.inclWater)), grep("-720-", colnames(dm.BC.morphonly))]
MF.BC.720.only <- MF.inclWater[grep("(BL|CR|FB)-720-", rownames(MF.inclWater)),]
ANOVA.BC.720.only <- adonis(dm.BC.720 ~ Morph, data = MF.BC.720.only, by = "margin")

# Fifth Timepoint
FifthTPNames <- rownames(MF.inclWater)[grep(paste0("^",TP[5],"$"), MF.inclWater$Time)]
toremove <- grep("W", FifthTPNames)
if (length(toremove) > 0) {
  FifthTPNames <- FifthTPNames[-toremove]
}
# Between morphs
dm.BC.FifthTP <- dm.BC.inclWater[FifthTPNames,FifthTPNames]
FB.BL.Fifthvalues <- as.vector(as.dist(dm.BC.FifthTP[grep("FB", rownames(dm.BC.FifthTP)), grep("BL", colnames(dm.BC.FifthTP))]))
FB.CR.Fifthvalues <- as.vector(as.dist(dm.BC.FifthTP[grep("FB", rownames(dm.BC.FifthTP)), grep("CR", colnames(dm.BC.FifthTP))]))
CR.BL.Fifthvalues <- as.vector(as.dist(dm.BC.FifthTP[grep("CR", rownames(dm.BC.FifthTP)), grep("BL", colnames(dm.BC.FifthTP))]))
# ANOVA
dm.BC.5760 <- dm.BC.inclWater[grep("(BL|CR|FB)-5760-", rownames(dm.BC.inclWater)), grep("-5760-", colnames(dm.BC.inclWater))]
MF.BC.5760.only <- MF.inclWater[grep("(BL|CR|FB)-5760-", rownames(MF.inclWater)),]
ANOVA.BC.5760.only <- adonis(dm.BC.5760 ~ Morph, data = MF.BC.5760.only, by = "margin")


# Combine into single tables

FB.BL.BC.ALL <- cbind(as.numeric(c(FB.BL.Firstvalues
                                     , FB.BL.Secondvalues
                                     , FB.BL.Thirdvalues
                                     , FB.BL.Fourthvalues
                                     , FB.BL.Fifthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.BL.Firstvalues))
                                       , rep(TP[2], length(FB.BL.Secondvalues))
                                       , rep(TP[3], length(FB.BL.Thirdvalues))
                                       , rep(TP[4], length(FB.BL.Fourthvalues))
                                       , rep(TP[5], length(FB.BL.Fifthvalues)))))
FB.CR.BC.ALL <- cbind(as.numeric(c(FB.CR.Firstvalues
                                     , FB.CR.Secondvalues
                                     , FB.CR.Thirdvalues
                                     , FB.CR.Fourthvalues
                                     , FB.CR.Fifthvalues))
                        , as.numeric(c(rep(TP[1], length(FB.CR.Firstvalues))
                                       , rep(TP[2], length(FB.CR.Secondvalues))
                                       , rep(TP[3], length(FB.CR.Thirdvalues))
                                       , rep(TP[4], length(FB.CR.Fourthvalues))
                                       , rep(TP[5], length(FB.CR.Fifthvalues)))))
CR.BL.BC.ALL <- cbind(as.numeric(c(CR.BL.Firstvalues
                                     , CR.BL.Secondvalues
                                     , CR.BL.Thirdvalues
                                     , CR.BL.Fourthvalues
                                     , CR.BL.Fifthvalues))
                        , as.numeric(c(rep(TP[1], length(CR.BL.Firstvalues))
                                       , rep(TP[2], length(CR.BL.Secondvalues))
                                       , rep(TP[3], length(CR.BL.Thirdvalues))
                                       , rep(TP[4], length(CR.BL.Fourthvalues))
                                       , rep(TP[5], length(CR.BL.Fifthvalues)))))

# FB.BL.BC.lm <- lm(FB.BL.BC.ALL[,1] ~ log(FB.BL.BC.ALL[,2]))
# summary(FB.BL.BC.lm)
# 
# CR.BL.BC.lm <- lm(CR.BL.BC.ALL[,1] ~ log(CR.BL.BC.ALL[,2]))
# summary(CR.BL.BC.lm)
# 
# FB.CR.BC.lm <- lm(FB.CR.BC.ALL[,1] ~ log(FB.CR.BC.ALL[,2]))
# summary(FB.CR.BC.lm)

# Plot the distances

FB.BL.BC.ALL.mean <- aggregate(FB.BL.BC.ALL, by = list(FB.BL.BC.ALL[,2]), mean)
FB.BL.BC.ALL.sd <- aggregate(FB.BL.BC.ALL, by = list(FB.BL.BC.ALL[,2]), sd)

CR.BL.BC.ALL.mean <- aggregate(CR.BL.BC.ALL, by = list(CR.BL.BC.ALL[,2]), mean)
CR.BL.BC.ALL.sd <- aggregate(CR.BL.BC.ALL, by = list(CR.BL.BC.ALL[,2]), sd)

FB.CR.BC.ALL.mean <- aggregate(FB.CR.BC.ALL, by = list(FB.CR.BC.ALL[,2]), mean)
FB.CR.BC.ALL.sd <- aggregate(FB.CR.BC.ALL, by = list(FB.CR.BC.ALL[,2]), sd)


######## PLOT DISP ###############
# ylimits <- c(min(FB.BL.BC.ALL[,1],FB.CR.BC.ALL[,1],CR.BL.BC.ALL[,1]), max(FB.BL.BC.ALL[,1],FB.CR.BC.ALL[,1],CR.BL.BC.ALL[,1]))
ylimits <- c(0.1,0.9)

# xvalues <- log(FB.BL.BC.ALL.mean[,1])
xvalues <- as.character(FB.BL.BC.ALL.mean[,1])
pdf(paste0("BETAPLOTS_H/DispOverTime_H_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Bray-Curtis)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1, at = c(1,2,3,4,5), labels = c('20 min','1 h', '6 h', '12 h','4 d'))
points(FB.BL.BC.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*0.99
       , x1 = c(1,2,3,4,5)*0.99
       , y0 = c(FB.BL.BC.ALL.mean[,2] - FB.BL.BC.ALL.sd[,2]/2)
       , y1 = c(FB.BL.BC.ALL.mean[,2] + FB.BL.BC.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.BC.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1
       , x1 = c(1,2,3,4,5)*1
       , y0 = c(FB.CR.BC.ALL.mean[,2] - FB.CR.BC.ALL.sd[,2]/2)
       , y1 = c(FB.CR.BC.ALL.mean[,2] + FB.CR.BC.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.BC.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1.01
       , x1 = c(1,2,3,4,5)*1.01
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
       , col = c("darkgreen","purple","grey")
       , lwd = 2)
dev.off()

####### PLOT MORPH ############
# NO 5760 BC

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


MorphColours <- c("darkorchid4","dodgerblue","salmon") 

pdf(paste0("BETAPLOTS_H/NMDS_H_",metric,"_Morph.pdf"), pointsize = 14)
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
      , col = MorphColours[1])
lines(NMDS.BC.BL[NMDS.BC.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.FB[NMDS.BC.FB.chull,]
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

# 5760 minutes
NMDS.BC.5760.only <- NMDS.BC.all$points[grep("(BL|CR|FB)-5760-", rownames(NMDS.BC.all$points)),]

NMDS.BC.5760.only.CR <- NMDS.BC.5760.only[grep("CR", rownames(NMDS.BC.5760.only)),]
NMDS.BC.5760.only.CR.chull <- chull(NMDS.BC.5760.only.CR)
NMDS.BC.5760.only.CR.chull <- c(NMDS.BC.5760.only.CR.chull, NMDS.BC.5760.only.CR.chull[1])

NMDS.BC.5760.only.BL <- NMDS.BC.5760.only[grep("BL", rownames(NMDS.BC.5760.only)),]
NMDS.BC.5760.only.BL.chull <- chull(NMDS.BC.5760.only.BL)
NMDS.BC.5760.only.BL.chull <- c(NMDS.BC.5760.only.BL.chull, NMDS.BC.5760.only.BL.chull[1])

NMDS.BC.5760.only.FB <- NMDS.BC.5760.only[grep("FB", rownames(NMDS.BC.5760.only)),]
NMDS.BC.5760.only.FB.chull <- chull(NMDS.BC.5760.only.FB)
NMDS.BC.5760.only.FB.chull <- c(NMDS.BC.5760.only.FB.chull, NMDS.BC.5760.only.FB.chull[1])


# pdf(paste0("BETAPLOTS_H/NMDS_",metric,"_TimebyMorph.pdf"),pointsize = 14, width = 7, height = 4)
# par(mfrow= c(1,4), oma = c(4,4,4,6))
# par(mar = c(4,0,4,0))
# plot(NMDS.BC.20.only
#      , main = "20 Minutes"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.BC.20.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p = ", ANOVA.BC.20.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.BC.20.only.CR[NMDS.BC.20.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.BC.20.only.BL[NMDS.BC.20.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.BC.20.only.FB[NMDS.BC.20.only.FB.chull,]
#       , col = MorphColours[3])
# par(mar = c(4,0,4,0))
# plot(NMDS.BC.60.only
#      , main = "1 Hour"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.BC.60.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p =  ",ANOVA.BC.60.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.BC.60.only.CR[NMDS.BC.60.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.BC.60.only.BL[NMDS.BC.60.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.BC.60.only.FB[NMDS.BC.60.only.FB.chull,]
#       , col = MorphColours[3])
# par(mar = c(4,0,4,0))
# plot(NMDS.BC.360.only
#      , main = "6 Hours"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.BC.360.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p = ", ANOVA.BC.360.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.BC.360.only.CR[NMDS.BC.360.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.BC.360.only.BL[NMDS.BC.360.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.BC.360.only.FB[NMDS.BC.360.only.FB.chull,]
#       , col = MorphColours[3])
# par(mar = c(4,0,4,0))
# plot(NMDS.BC.720.only
#      , main = "12 Hours"
#      , pch = 21
#      , col = "black"
#      , bg = MorphColours[factor(MF.BC.720.only$Morph)]
#      , xaxt = 'n'
#      , yaxt = 'n'
#      , xlab = paste0("p = ", ANOVA.BC.720.only$aov.tab[6]$`Pr(>F)`[1])
#      , ylab = ''
#      , cex = 2
#      , cex.lab = 2
#      , cex.main = 2
# )
# lines(NMDS.BC.720.only.CR[NMDS.BC.720.only.CR.chull,]
#       , col = MorphColours[1])
# lines(NMDS.BC.720.only.BL[NMDS.BC.720.only.BL.chull,]
#       , col = MorphColours[2])
# lines(NMDS.BC.720.only.FB[NMDS.BC.720.only.FB.chull,]
#       , col = MorphColours[3])
# #STOP
# dev.off()

############ COMBO DISP BETA ################
# Disp and beta through time combined
# xvalues <- log(FB.BL.BC.ALL.mean[,1])
xvalues <- as.character(FB.BL.BC.ALL.mean[,1])

# Change factor of individual plots
MF.BC.20.only$Morph <- factor(MF.BC.20.only$Morph, levels = c("CR","BL","FB"))
MF.BC.60.only$Morph <- factor(MF.BC.60.only$Morph, levels = c("CR","BL","FB"))
MF.BC.360.only$Morph <- factor(MF.BC.360.only$Morph, levels = c("CR","BL","FB"))
MF.BC.720.only$Morph <- factor(MF.BC.720.only$Morph, levels = c("CR","BL","FB"))
MF.BC.5760.only$Morph <- factor(MF.BC.5760.only$Morph, levels = c("CR","BL","FB"))


pdf(paste0("BETAPLOTS_H/COMBO_H_dispbeta_",metric,".pdf"), pointsize = 14, width = 10, height = 7)
par(fig = c(0,0.8,0.3,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = ''
     , ylab = "Distance (Bray-Curtis)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1
     , at = c(1,2,3,4,5)
     , labels = c("20 min","1 h","6 h","12 h","4 d")
     , las = 2)
points(FB.BL.BC.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*0.99
       , x1 = c(1,2,3,4,5)*0.99
       , y0 = c(FB.BL.BC.ALL.mean[,2] - FB.BL.BC.ALL.sd[,2])
       , y1 = c(FB.BL.BC.ALL.mean[,2] + FB.BL.BC.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.BC.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1
       , x1 = c(1,2,3,4,5)*1
       , y0 = c(FB.CR.BC.ALL.mean[,2] - FB.CR.BC.ALL.sd[,2])
       , y1 = c(FB.CR.BC.ALL.mean[,2] + FB.CR.BC.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.BC.ALL.mean[,2] ~ c(1,2,3,4,5)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5)*1.01
       , x1 = c(1,2,3,4,5)*1.01
       , y0 = c(CR.BL.BC.ALL.mean[,2] - CR.BL.BC.ALL.sd[,2])
       , y1 = c(CR.BL.BC.ALL.mean[,2] + CR.BL.BC.ALL.sd[,2])
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
       , legend = c("FB:CR", "FB:BL", "CR:BL")
       , lty = 1
       , col = c("darkgreen","purple","grey")
       , lwd = 2 
)
# EACH GETS 0.15 SPACE TOTAL;
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.049,0.1932,0,0.4), new = TRUE)
plot(NMDS.BC.20.only
     # , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.BC.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.20.only$aov.tab[1]$Df[1],",",ANOVA.BC.20.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.20.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.20.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.20.only.CR[NMDS.BC.20.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.20.only.BL[NMDS.BC.20.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.20.only.FB[NMDS.BC.20.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.2032,0.3474,0,0.4), new = TRUE)
plot(NMDS.BC.60.only
     # , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.BC.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.60.only$aov.tab[1]$Df[1],",",ANOVA.BC.60.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.60.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.60.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.60.only.CR[NMDS.BC.60.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.60.only.BL[NMDS.BC.60.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.60.only.FB[NMDS.BC.60.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.3574,0.5016,0,0.4), new = TRUE)
plot(NMDS.BC.360.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.BC.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.360.only$aov.tab[1]$Df[1],",",ANOVA.BC.360.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.360.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.360.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.360.only.CR[NMDS.BC.360.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.360.only.BL[NMDS.BC.360.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.360.only.FB[NMDS.BC.360.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.5116,0.6558,0,0.4), new = TRUE)
plot(NMDS.BC.720.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.BC.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.720.only$aov.tab[1]$Df[1],",",ANOVA.BC.720.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.720.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.720.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.720.only.CR[NMDS.BC.720.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.720.only.BL[NMDS.BC.720.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.720.only.FB[NMDS.BC.720.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.6658,0.81,0,0.4), new = TRUE)
plot(NMDS.BC.5760.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.BC.5760.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.5760.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.5760.only$aov.tab[1]$Df[1],",",ANOVA.BC.5760.only$aov.tab[1]$Df[3]) 
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.5760.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.5760.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.5760.only.CR[NMDS.BC.5760.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.5760.only.BL[NMDS.BC.5760.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.5760.only.FB[NMDS.BC.5760.only.FB.chull,]
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

betadisp.H.BC.FB <- betadisp.BC.time$distances[grep("FB", names(betadisp.BC.time$distances))]
betadisp.H.BC.BL <- betadisp.BC.time$distances[grep("BL", names(betadisp.BC.time$distances))]
betadisp.H.BC.CR <- betadisp.BC.time$distances[grep("CR", names(betadisp.BC.time$distances))]

MF.H.FB <- MF.incl5760[grep("FB", rownames(MF.incl5760)),]
MF.H.BL <- MF.incl5760[grep("BL", rownames(MF.incl5760)),]
MF.H.CR <- MF.incl5760[grep("CR", rownames(MF.incl5760)),]

betadisp.FB.H.agg <- aggregate(betadisp.H.BC.FB, by = list(MF.H.FB$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.BL.H.agg <- aggregate(betadisp.H.BC.BL, by = list(MF.H.BL$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.CR.H.agg <- aggregate(betadisp.H.BC.CR, by = list(MF.H.CR$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )

betadisp.BC.time.forstat <- cbind(betadisp.BC.time$distances, MF.incl5760[unlist(lapply(names(betadisp.BC.time$distances), function(x) {grep(paste0("^",x,"$"), rownames(MF.incl5760))})),])
colnames(betadisp.BC.time.forstat)[1] <- c("Distance")
ANOVA.betadisp.BC <- anova(lm(Distance ~ Time*Morph, data = betadisp.BC.time.forstat))
capture.output(ANOVA.betadisp.BC, file = "./BETAPLOTS_H/ANOVA.betadisp.BC.txt")

xvalues <- c("20","60","360","720","5760")
ylimits <- c(0,0.6)
pdf(paste0("./BETAPLOTS_H/BetaDisp_H_",metric,"_eachmorph.pdf"),pointsize = 14)
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = 'Time'
     , ylab = "Distance (Bray-Curtis)"
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
dm.BC.5760<- dm.BC.inclWater[grep("(CR|FB|BL)-5760", rownames(dm.BC.inclWater)),grep("(CR|FB|BL)-5760", colnames(dm.BC.inclWater))]

NMDS.BC.5760 <- isoMDS(as.matrix(dm.BC.5760), y = cmdscale(as.matrix(dm.BC.5760), 2))

# NMDS.BC.morphallTimePoints <- NMDS.BC.morphAllTime$points[grep("(CR|BL|FB)-5760", rownames(NMDS.BC.morphAllTime$points)),]
MF.5760 <- MF[sapply(rownames(NMDS.BC.5760$points), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
}),]
MF.5760$Morph <- factor(MF.5760$Morph, levels = c("CR","BL","FB"))

NMDS.BC.5760.CR <- NMDS.BC.5760$points[grep("CR", rownames(NMDS.BC.5760$points)),]
NMDS.BC.5760.CR.chull <- chull(NMDS.BC.5760.CR)
NMDS.BC.5760.CR.chull <- c(NMDS.BC.5760.CR.chull, NMDS.BC.5760.CR.chull[1])

NMDS.BC.5760.BL <- NMDS.BC.5760$points[grep("BL", rownames(NMDS.BC.5760$points)),]
NMDS.BC.5760.BL.chull <- chull(NMDS.BC.5760.BL)
NMDS.BC.5760.BL.chull <- c(NMDS.BC.5760.BL.chull, NMDS.BC.5760.BL.chull[1])

NMDS.BC.5760.FB <- NMDS.BC.5760$points[grep("FB", rownames(NMDS.BC.5760$points)),]
NMDS.BC.5760.FB.chull <- chull(NMDS.BC.5760.FB)
NMDS.BC.5760.FB.chull <- c(NMDS.BC.5760.FB.chull, NMDS.BC.5760.FB.chull[1])

# ANOVA.BC.5760 <- adonis(dm.BC.5760 ~ Morph, data = MF.5760)
# # ANOSIM.BC.5760 <- anosim(dat = dm.BC.5760, grouping = MF.5760$Morph)
# capture.output(ANOVA.BC.5760, file = paste0("BETAPLOTS_H/anova_",metric,"_5760only.txt"))

pdf(paste0("BETAPLOTS_H/NMDS_",metric,"_5760Only.pdf"), pointsize = 14)
par(fig = c(0,0.75,0,1))
plot(NMDS.BC.5760$points
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.5760$Morph)]
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , sub = paste0("Stress: ", round(NMDS.BC.5760$stress,2)/100)
     , cex = 1.5
     , main = "NMDS of community composition (4 days)")
lines(NMDS.BC.5760.CR[NMDS.BC.5760.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.5760.BL[NMDS.BC.5760.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.5760.FB[NMDS.BC.5760.FB.chull,]
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


# Combine individual stats and print out

allindividualTests <- c("ANOVA.BC.20.only"
                        ,"ANOVA.BC.60.only"
                        ,"ANOVA.BC.360.only"
                        ,"ANOVA.BC.720.only"
                        ,"ANOVA.BC.5760.only"
)
for (i in allindividualTests) {
  # print(get(i))
  capture.output(get(i), file = paste0("./BETAPLOTS_H/individualtests/",i,".txt"))
}

######## --PM-- ###########
### NMDS #####
rownames(dm.BC.P.morphonly) <- gsub(".","-", rownames(dm.BC.P.morphonly), fixed = TRUE)
colnames(dm.BC.P.morphonly) <- gsub(".","-", colnames(dm.BC.P.morphonly), fixed = TRUE)

rownames(MF.P.morphkeep) <- gsub(".","-", rownames(MF.P.morphkeep), fixed = TRUE)
rownames(MF.P.inclWater) <- gsub(".","-", rownames(MF.P.inclWater), fixed = TRUE)



NMDS.BC.P.morphonly <- isoMDS(as.matrix(dm.BC.P.morphonly), y = cmdscale(as.matrix(dm.BC.P.morphonly), 2))
NMDS.BC.P.all <- isoMDS(as.matrix(dm.BC.P.inclWater), y = cmdscale(as.matrix(dm.BC.P.inclWater), 2))

###### STATS ##########
MF.P.morphkeep <- MF.P.morphkeep[,c('Morph','Time','Type','TypeMorphTime')]
MF.P.morphkeep$Morph <- factor(MF.P.morphkeep$Morph, levels = c('CR','BL','FB'))
MF.P.morphkeep$Time <- factor(MF.P.morphkeep$Time, levels = c('20','60','180','360','720','1440'))
MF.P.morphkeep$Type <- factor(MF.P.morphkeep$Type, levels = c('P','H'))


# dm.BC.P.morphonly.H <- dm.BC.P.morphonly[grep(".", rownames(dm.BC.P.morphonly), fixed = TRUE),grep(".", colnames(dm.BC.P.morphonly), fixed = TRUE)]
# dm.BC.P.morphonly.P <- dm.BC.P.morphonly[grep("-", rownames(dm.BC.P.morphonly), fixed = TRUE),grep("-", colnames(dm.BC.P.morphonly), fixed = TRUE)]

ANOVA.BC.P.morphtime <- adonis(dm.BC.P.morphonly ~ Time*Morph, data = MF.P.morphkeep, by = "margin")
capture.output(ANOVA.BC.P.morphtime, file = paste0("BETAPLOTS_P/adonis_", metric,"_PM.txt"))
# ANOSIM.BC.P.morphtime <- anosim(dm.BC.P.morphonly, grouping = MF.P.morphkeep$Morph)

# Dispersion across time and between morphs
dist.BC.P.morphonly <- as.dist(dm.BC.P.morphonly)
betadisp.BC.P.time <- betadisper(d = dist.BC.P.morphonly, group = MF.P.morphkeep$Time)
betadisp.BC.P.morph <- betadisper(d = dist.BC.P.morphonly, group = MF.P.morphkeep$Morph)

capture.output(betadisp.BC.P.time, file = paste0("BETAPLOTS_P/betadispTime_", metric, "_PM.txt"))
capture.output(betadisp.BC.P.morph, file = paste0("BETAPLOTS_P/betadispMorph_", metric, "_PM.txt"))


# Grep values for each timepoint and morph, and then plot dispersion within each
TP <- levels(factor(MF.P.morphkeep$Time))
MorphTypes <- levels(factor(MF.P.morphkeep$Morph))

# First Timepoint
FirstTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[1],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.BC.P.FirstTP <- dm.BC.P.morphonly[FirstTPNames,FirstTPNames]
FB.BL.Firstvalues <- as.vector(as.dist(dm.BC.P.FirstTP[grep("FB", rownames(dm.BC.P.FirstTP)), grep("BL", colnames(dm.BC.P.FirstTP))]))
FB.CR.Firstvalues <- as.vector(as.dist(dm.BC.P.FirstTP[grep("FB", rownames(dm.BC.P.FirstTP)), grep("CR", colnames(dm.BC.P.FirstTP))]))
CR.BL.Firstvalues <- as.vector(as.dist(dm.BC.P.FirstTP[grep("CR", rownames(dm.BC.P.FirstTP)), grep("BL", colnames(dm.BC.P.FirstTP))]))
# ANOVA
dm.BC.P.20 <- dm.BC.P.morphonly[grep("^20-", rownames(dm.BC.P.morphonly)), grep("^20-", colnames(dm.BC.P.morphonly))]
MF.P.BC.P.20.only <- MF.P.morphkeep[grep("^20-", rownames(MF.P.morphkeep)),]
ANOVA.BC.P.20.only <- adonis(dm.BC.P.20 ~ Morph, data = MF.P.BC.P.20.only, by = "margin")

# Second Timepoint
SecondTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[2],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.BC.P.SecondTP <- dm.BC.P.morphonly[SecondTPNames,SecondTPNames]
FB.BL.Secondvalues <- as.vector(as.dist(dm.BC.P.SecondTP[grep("FB", rownames(dm.BC.P.SecondTP)), grep("BL", colnames(dm.BC.P.SecondTP))]))
FB.CR.Secondvalues <- as.vector(as.dist(dm.BC.P.SecondTP[grep("FB", rownames(dm.BC.P.SecondTP)), grep("CR", colnames(dm.BC.P.SecondTP))]))
CR.BL.Secondvalues <- as.vector(as.dist(dm.BC.P.SecondTP[grep("CR", rownames(dm.BC.P.SecondTP)), grep("BL", colnames(dm.BC.P.SecondTP))]))
# ANOVA
dm.BC.P.60 <- dm.BC.P.morphonly[grep("^60-", rownames(dm.BC.P.morphonly)), grep("^60-", colnames(dm.BC.P.morphonly))]
MF.P.BC.P.60.only <- MF.P.morphkeep[grep("^60-", rownames(MF.P.morphkeep)),]
ANOVA.BC.P.60.only <- adonis(dm.BC.P.60 ~ Morph, data = MF.P.BC.P.60.only, by = "margin")

# Third Timepoint
ThirdTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[3],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.BC.P.ThirdTP <- dm.BC.P.morphonly[ThirdTPNames,ThirdTPNames]
FB.BL.Thirdvalues <- as.vector(as.dist(dm.BC.P.ThirdTP[grep("FB", rownames(dm.BC.P.ThirdTP)), grep("BL", colnames(dm.BC.P.ThirdTP))]))
FB.CR.Thirdvalues <- as.vector(as.dist(dm.BC.P.ThirdTP[grep("FB", rownames(dm.BC.P.ThirdTP)), grep("CR", colnames(dm.BC.P.ThirdTP))]))
CR.BL.Thirdvalues <- as.vector(as.dist(dm.BC.P.ThirdTP[grep("CR", rownames(dm.BC.P.ThirdTP)), grep("BL", colnames(dm.BC.P.ThirdTP))]))
# ANOVA
dm.BC.P.180 <- dm.BC.P.morphonly[grep("^180-", rownames(dm.BC.P.morphonly)), grep("^180-", colnames(dm.BC.P.morphonly))]
MF.P.BC.P.180.only <- MF.P.morphkeep[grep("^180-", rownames(MF.P.morphkeep)),]
ANOVA.BC.P.180.only <- adonis(dm.BC.P.180 ~ Morph, data = MF.P.BC.P.180.only, by = "margin")


# Fourth Timepoint
FourthTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[4],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.BC.P.FourthTP <- dm.BC.P.morphonly[FourthTPNames,FourthTPNames]
FB.BL.Fourthvalues <- as.vector(as.dist(dm.BC.P.FourthTP[grep("FB", rownames(dm.BC.P.FourthTP)), grep("BL", colnames(dm.BC.P.FourthTP))]))
FB.CR.Fourthvalues <- as.vector(as.dist(dm.BC.P.FourthTP[grep("FB", rownames(dm.BC.P.FourthTP)), grep("CR", colnames(dm.BC.P.FourthTP))]))
CR.BL.Fourthvalues <- as.vector(as.dist(dm.BC.P.FourthTP[grep("CR", rownames(dm.BC.P.FourthTP)), grep("BL", colnames(dm.BC.P.FourthTP))]))
# ANOVA
dm.BC.P.360 <- dm.BC.P.morphonly[grep("^360-", rownames(dm.BC.P.morphonly)), grep("^360-", colnames(dm.BC.P.morphonly))]
MF.P.BC.P.360.only <- MF.P.morphkeep[grep("^360-", rownames(MF.P.morphkeep)),]
ANOVA.BC.P.360.only <- adonis(dm.BC.P.360 ~ Morph, data = MF.P.BC.P.360.only, by = "margin")

# Fifth Timepoint
FifthTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[5],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.BC.P.FifthTP <- dm.BC.P.morphonly[FifthTPNames,FifthTPNames]
FB.BL.Fifthvalues <- as.vector(as.dist(dm.BC.P.FifthTP[grep("FB", rownames(dm.BC.P.FifthTP)), grep("BL", colnames(dm.BC.P.FifthTP))]))
FB.CR.Fifthvalues <- as.vector(as.dist(dm.BC.P.FifthTP[grep("FB", rownames(dm.BC.P.FifthTP)), grep("CR", colnames(dm.BC.P.FifthTP))]))
CR.BL.Fifthvalues <- as.vector(as.dist(dm.BC.P.FifthTP[grep("CR", rownames(dm.BC.P.FifthTP)), grep("BL", colnames(dm.BC.P.FifthTP))]))
# ANOVA
dm.BC.P.720 <- dm.BC.P.morphonly[grep("^720-", rownames(dm.BC.P.morphonly)), grep("^720-", colnames(dm.BC.P.morphonly))]
MF.P.BC.P.720.only <- MF.P.morphkeep[grep("^720-", rownames(MF.P.morphkeep)),]
ANOVA.BC.P.720.only <- adonis(dm.BC.P.720 ~ Morph, data = MF.P.BC.P.720.only, by = "margin")

# Sixth Timepoint
SixthTPNames <- rownames(MF.P.morphkeep)[grep(paste0("^",TP[6],"$"), MF.P.morphkeep$Time)]
# Between morphs
dm.BC.P.SixthTP <- dm.BC.P.morphonly[SixthTPNames,SixthTPNames]
FB.BL.Sixthvalues <- as.vector(as.dist(dm.BC.P.SixthTP[grep("FB", rownames(dm.BC.P.SixthTP)), grep("BL", colnames(dm.BC.P.SixthTP))]))
FB.CR.Sixthvalues <- as.vector(as.dist(dm.BC.P.SixthTP[grep("FB", rownames(dm.BC.P.SixthTP)), grep("CR", colnames(dm.BC.P.SixthTP))]))
CR.BL.Sixthvalues <- as.vector(as.dist(dm.BC.P.SixthTP[grep("CR", rownames(dm.BC.P.SixthTP)), grep("BL", colnames(dm.BC.P.SixthTP))]))
# ANOVA
dm.BC.P.1440<- dm.BC.P.morphonly[grep("^1440-", rownames(dm.BC.P.morphonly)), grep("^1440-", colnames(dm.BC.P.morphonly))]
MF.P.BC.P.1440.only <- MF.P.morphkeep[grep("^1440-", rownames(MF.P.morphkeep)),]
ANOVA.BC.P.1440.only <- adonis(dm.BC.P.1440 ~ Morph, data = MF.P.BC.P.1440.only, by = "margin")



# Combine into single tables

FB.BL.BC.P.ALL <- cbind(as.numeric(c(FB.BL.Firstvalues
                                      , FB.BL.Secondvalues
                                      , FB.BL.Thirdvalues
                                      , FB.BL.Fourthvalues
                                      , FB.BL.Fifthvalues
                                      , FB.BL.Sixthvalues
))
, as.numeric(c(rep(TP[1], length(FB.BL.Firstvalues))
               , rep(TP[2], length(FB.BL.Secondvalues))
               , rep(TP[3], length(FB.BL.Thirdvalues))
               , rep(TP[4], length(FB.BL.Fourthvalues))
               , rep(TP[5], length(FB.BL.Fifthvalues))
               , rep(TP[6], length(FB.BL.Sixthvalues))
)))
FB.CR.BC.P.ALL <- cbind(as.numeric(c(FB.CR.Firstvalues
                                      , FB.CR.Secondvalues
                                      , FB.CR.Thirdvalues
                                      , FB.CR.Fourthvalues
                                      , FB.CR.Fifthvalues
                                      , FB.CR.Sixthvalues
))
, as.numeric(c(rep(TP[1], length(FB.CR.Firstvalues))
               , rep(TP[2], length(FB.CR.Secondvalues))
               , rep(TP[3], length(FB.CR.Thirdvalues))
               , rep(TP[4], length(FB.CR.Fourthvalues))
               , rep(TP[5], length(FB.CR.Fifthvalues))
               , rep(TP[6], length(FB.CR.Sixthvalues))
)))
CR.BL.BC.P.ALL <- cbind(as.numeric(c(CR.BL.Firstvalues
                                      , CR.BL.Secondvalues
                                      , CR.BL.Thirdvalues
                                      , CR.BL.Fourthvalues
                                      , CR.BL.Fifthvalues
                                      , CR.BL.Sixthvalues
))
, as.numeric(c(rep(TP[1], length(CR.BL.Firstvalues))
               , rep(TP[2], length(CR.BL.Secondvalues))
               , rep(TP[3], length(CR.BL.Thirdvalues))
               , rep(TP[4], length(CR.BL.Fourthvalues))
               , rep(TP[5], length(CR.BL.Fifthvalues))
               , rep(TP[6], length(CR.BL.Sixthvalues))
)))

# FB.BL.BC.P.lm <- lm(FB.BL.BC.P.ALL[,1] ~ log(FB.BL.BC.P.ALL[,2]))
# summary(FB.BL.BC.P.lm)
# 
# CR.BL.BC.P.lm <- lm(CR.BL.BC.P.ALL[,1] ~ log(CR.BL.BC.P.ALL[,2]))
# summary(CR.BL.BC.P.lm)
# 
# FB.CR.BC.P.lm <- lm(FB.CR.BC.P.ALL[,1] ~ log(FB.CR.BC.P.ALL[,2]))
# summary(FB.CR.BC.P.lm)

# Plot the distances

FB.BL.BC.P.ALL.mean <- aggregate(FB.BL.BC.P.ALL, by = list(FB.BL.BC.P.ALL[,2]), mean)
FB.BL.BC.P.ALL.sd <- aggregate(FB.BL.BC.P.ALL, by = list(FB.BL.BC.P.ALL[,2]), sd)

CR.BL.BC.P.ALL.mean <- aggregate(CR.BL.BC.P.ALL, by = list(CR.BL.BC.P.ALL[,2]), mean)
CR.BL.BC.P.ALL.sd <- aggregate(CR.BL.BC.P.ALL, by = list(CR.BL.BC.P.ALL[,2]), sd)

FB.CR.BC.P.ALL.mean <- aggregate(FB.CR.BC.P.ALL, by = list(FB.CR.BC.P.ALL[,2]), mean)
FB.CR.BC.P.ALL.sd <- aggregate(FB.CR.BC.P.ALL, by = list(FB.CR.BC.P.ALL[,2]), sd)

######## PLOT DISP ###############
# ylimits <- c(min(FB.BL.BC.P.ALL[,1],FB.CR.BC.P.ALL[,1],CR.BL.BC.P.ALL[,1]), max(FB.BL.BC.P.ALL[,1],FB.CR.BC.P.ALL[,1],CR.BL.BC.P.ALL[,1]))
ylimits <- c(0.3,1)

# xvalues <- log(FB.BL.BC.P.ALL.mean[,1])
xvalues <- as.character(FB.BL.BC.P.ALL.mean[,1])
pdf(paste0("BETAPLOTS_P/DispOverTime_P_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.8,0,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = "Time"
     , ylab = "Distance (Bray-Curtis)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1, at = c(1,2,3,4,5,6), labels = c('20 min','1 h','3 h','6 h', '12 h', '24 h'))
points(FB.BL.BC.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*0.99
       , x1 = c(1,2,3,4,5,6)*0.99
       , y0 = c(FB.BL.BC.P.ALL.mean[,2] - FB.BL.BC.P.ALL.sd[,2]/2)
       , y1 = c(FB.BL.BC.P.ALL.mean[,2] + FB.BL.BC.P.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.BC.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*1
       , x1 = c(1,2,3,4,5,6)*1
       , y0 = c(FB.CR.BC.P.ALL.mean[,2] - FB.CR.BC.P.ALL.sd[,2]/2)
       , y1 = c(FB.CR.BC.P.ALL.mean[,2] + FB.CR.BC.P.ALL.sd[,2]/2)
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.BC.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*1.01
       , x1 = c(1,2,3,4,5,6)*1.01
       , y0 = c(CR.BL.BC.P.ALL.mean[,2] - CR.BL.BC.P.ALL.sd[,2]/2)
       , y1 = c(CR.BL.BC.P.ALL.mean[,2] + CR.BL.BC.P.ALL.sd[,2]/2)
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
# NOCHLP BC.P

NMDS.BC.P.morphonly$points <- NMDS.BC.P.morphonly$points[, c(2,1)]

# MAKE POLYGONS for plotting
NMDS.BC.P.CR <- NMDS.BC.P.morphonly$points[grep("CR", rownames(NMDS.BC.P.morphonly$points)),]
NMDS.BC.P.BL <- NMDS.BC.P.morphonly$points[grep("BL", rownames(NMDS.BC.P.morphonly$points)),]
NMDS.BC.P.FB <- NMDS.BC.P.morphonly$points[grep("FB", rownames(NMDS.BC.P.morphonly$points)),]

NMDS.BC.P.CR.chull <- chull(NMDS.BC.P.CR)
NMDS.BC.P.CR.chull <- c(NMDS.BC.P.CR.chull, NMDS.BC.P.CR.chull[1])

NMDS.BC.P.BL.chull <- chull(NMDS.BC.P.BL)
NMDS.BC.P.BL.chull <- c(NMDS.BC.P.BL.chull, NMDS.BC.P.BL.chull[1])

NMDS.BC.P.FB.chull <- chull(NMDS.BC.P.FB)
NMDS.BC.P.FB.chull <- c(NMDS.BC.P.FB.chull, NMDS.BC.P.FB.chull[1])


MorphColours <- c("darkorchid4","dodgerblue","salmon") 

pdf(paste0("BETAPLOTS_P/NMDS_P_",metric,"_Morph.pdf"), pointsize = 14)
par(fig = c(0,0.8,0,1))
plot(NMDS.BC.P.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.morphkeep$Morph)]
     , sub = round(NMDS.BC.P.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
)
lines(NMDS.BC.P.CR[NMDS.BC.P.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.BL[NMDS.BC.P.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.FB[NMDS.BC.P.FB.chull,]
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
plot(NMDS.BC.P.morphonly$points
     , main = "NMDS of Artificial Seaweed Shapes"
     , pch = 21
     , col = "black"
     , bg = TimeColours[factor(MF.P.morphkeep$Time)]
     , sub = round(NMDS.BC.P.morphonly$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2
)
# lines(NMDS.BC.P.CR[NMDS.BC.P.CR.chull,]
#       , col = "red")
# lines(NMDS.BC.P.BL[NMDS.BC.P.BL.chull,]
#       , col = "magenta")
# lines(NMDS.BC.P.FB[NMDS.BC.P.FB.chull,]
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



# TIME-- ONE at a time


# 20 minutes
NMDS.BC.P.20.only <- NMDS.BC.P.morphonly$points[grep("^20-", rownames(NMDS.BC.P.morphonly$points)),]

NMDS.BC.P.20.only.CR <- NMDS.BC.P.20.only[grep("CR", rownames(NMDS.BC.P.20.only)),]
NMDS.BC.P.20.only.CR.chull <- chull(NMDS.BC.P.20.only.CR)
NMDS.BC.P.20.only.CR.chull <- c(NMDS.BC.P.20.only.CR.chull, NMDS.BC.P.20.only.CR.chull[1])

NMDS.BC.P.20.only.BL <- NMDS.BC.P.20.only[grep("BL", rownames(NMDS.BC.P.20.only)),]
NMDS.BC.P.20.only.BL.chull <- chull(NMDS.BC.P.20.only.BL)
NMDS.BC.P.20.only.BL.chull <- c(NMDS.BC.P.20.only.BL.chull, NMDS.BC.P.20.only.BL.chull[1])

NMDS.BC.P.20.only.FB <- NMDS.BC.P.20.only[grep("FB", rownames(NMDS.BC.P.20.only)),]
NMDS.BC.P.20.only.FB.chull <- chull(NMDS.BC.P.20.only.FB)
NMDS.BC.P.20.only.FB.chull <- c(NMDS.BC.P.20.only.FB.chull, NMDS.BC.P.20.only.FB.chull[1])

# 60 minutes
NMDS.BC.P.60.only <- NMDS.BC.P.morphonly$points[grep("^60-", rownames(NMDS.BC.P.morphonly$points)),]

NMDS.BC.P.60.only.CR <- NMDS.BC.P.60.only[grep("CR", rownames(NMDS.BC.P.60.only)),]
NMDS.BC.P.60.only.CR.chull <- chull(NMDS.BC.P.60.only.CR)
NMDS.BC.P.60.only.CR.chull <- c(NMDS.BC.P.60.only.CR.chull, NMDS.BC.P.60.only.CR.chull[1])

NMDS.BC.P.60.only.BL <- NMDS.BC.P.60.only[grep("BL", rownames(NMDS.BC.P.60.only)),]
NMDS.BC.P.60.only.BL.chull <- chull(NMDS.BC.P.60.only.BL)
NMDS.BC.P.60.only.BL.chull <- c(NMDS.BC.P.60.only.BL.chull, NMDS.BC.P.60.only.BL.chull[1])

NMDS.BC.P.60.only.FB <- NMDS.BC.P.60.only[grep("FB", rownames(NMDS.BC.P.60.only)),]
NMDS.BC.P.60.only.FB.chull <- chull(NMDS.BC.P.60.only.FB)
NMDS.BC.P.60.only.FB.chull <- c(NMDS.BC.P.60.only.FB.chull, NMDS.BC.P.60.only.FB.chull[1])

# 180 minutes
NMDS.BC.P.180.only <- NMDS.BC.P.morphonly$points[grep("^180-", rownames(NMDS.BC.P.morphonly$points)),]

NMDS.BC.P.180.only.CR <- NMDS.BC.P.180.only[grep("CR", rownames(NMDS.BC.P.180.only)),]
NMDS.BC.P.180.only.CR.chull <- chull(NMDS.BC.P.180.only.CR)
NMDS.BC.P.180.only.CR.chull <- c(NMDS.BC.P.180.only.CR.chull, NMDS.BC.P.180.only.CR.chull[1])

NMDS.BC.P.180.only.BL <- NMDS.BC.P.180.only[grep("BL", rownames(NMDS.BC.P.180.only)),]
NMDS.BC.P.180.only.BL.chull <- chull(NMDS.BC.P.180.only.BL)
NMDS.BC.P.180.only.BL.chull <- c(NMDS.BC.P.180.only.BL.chull, NMDS.BC.P.180.only.BL.chull[1])

NMDS.BC.P.180.only.FB <- NMDS.BC.P.180.only[grep("FB", rownames(NMDS.BC.P.180.only)),]
NMDS.BC.P.180.only.FB.chull <- chull(NMDS.BC.P.180.only.FB)
NMDS.BC.P.180.only.FB.chull <- c(NMDS.BC.P.180.only.FB.chull, NMDS.BC.P.180.only.FB.chull[1])



# 360 minutes
NMDS.BC.P.360.only <- NMDS.BC.P.morphonly$points[grep("^360-", rownames(NMDS.BC.P.morphonly$points)),]


NMDS.BC.P.360.only.CR <- NMDS.BC.P.360.only[grep("CR", rownames(NMDS.BC.P.360.only)),]
NMDS.BC.P.360.only.CR.chull <- chull(NMDS.BC.P.360.only.CR)
NMDS.BC.P.360.only.CR.chull <- c(NMDS.BC.P.360.only.CR.chull, NMDS.BC.P.360.only.CR.chull[1])

NMDS.BC.P.360.only.BL <- NMDS.BC.P.360.only[grep("BL", rownames(NMDS.BC.P.360.only)),]
NMDS.BC.P.360.only.BL.chull <- chull(NMDS.BC.P.360.only.BL)
NMDS.BC.P.360.only.BL.chull <- c(NMDS.BC.P.360.only.BL.chull, NMDS.BC.P.360.only.BL.chull[1])

NMDS.BC.P.360.only.FB <- NMDS.BC.P.360.only[grep("FB", rownames(NMDS.BC.P.360.only)),]
NMDS.BC.P.360.only.FB.chull <- chull(NMDS.BC.P.360.only.FB)
NMDS.BC.P.360.only.FB.chull <- c(NMDS.BC.P.360.only.FB.chull, NMDS.BC.P.360.only.FB.chull[1])


# 720 minutes
NMDS.BC.P.720.only <- NMDS.BC.P.morphonly$points[grep("^720-", rownames(NMDS.BC.P.morphonly$points)),]


NMDS.BC.P.720.only.CR <- NMDS.BC.P.720.only[grep("CR", rownames(NMDS.BC.P.720.only)),]
NMDS.BC.P.720.only.CR.chull <- chull(NMDS.BC.P.720.only.CR)
NMDS.BC.P.720.only.CR.chull <- c(NMDS.BC.P.720.only.CR.chull, NMDS.BC.P.720.only.CR.chull[1])

NMDS.BC.P.720.only.BL <- NMDS.BC.P.720.only[grep("BL", rownames(NMDS.BC.P.720.only)),]
NMDS.BC.P.720.only.BL.chull <- chull(NMDS.BC.P.720.only.BL)
NMDS.BC.P.720.only.BL.chull <- c(NMDS.BC.P.720.only.BL.chull, NMDS.BC.P.720.only.BL.chull[1])

NMDS.BC.P.720.only.FB <- NMDS.BC.P.720.only[grep("FB", rownames(NMDS.BC.P.720.only)),]
NMDS.BC.P.720.only.FB.chull <- chull(NMDS.BC.P.720.only.FB)
NMDS.BC.P.720.only.FB.chull <- c(NMDS.BC.P.720.only.FB.chull, NMDS.BC.P.720.only.FB.chull[1])


# 1440 minutes
NMDS.BC.P.1440.only <- NMDS.BC.P.morphonly$points[grep("^1440-", rownames(NMDS.BC.P.morphonly$points)),]


NMDS.BC.P.1440.only.CR <- NMDS.BC.P.1440.only[grep("CR", rownames(NMDS.BC.P.1440.only)),]
NMDS.BC.P.1440.only.CR.chull <- chull(NMDS.BC.P.1440.only.CR)
NMDS.BC.P.1440.only.CR.chull <- c(NMDS.BC.P.1440.only.CR.chull, NMDS.BC.P.1440.only.CR.chull[1])

NMDS.BC.P.1440.only.BL <- NMDS.BC.P.1440.only[grep("BL", rownames(NMDS.BC.P.1440.only)),]
NMDS.BC.P.1440.only.BL.chull <- chull(NMDS.BC.P.1440.only.BL)
NMDS.BC.P.1440.only.BL.chull <- c(NMDS.BC.P.1440.only.BL.chull, NMDS.BC.P.1440.only.BL.chull[1])

NMDS.BC.P.1440.only.FB <- NMDS.BC.P.1440.only[grep("FB", rownames(NMDS.BC.P.1440.only)),]
NMDS.BC.P.1440.only.FB.chull <- chull(NMDS.BC.P.1440.only.FB)
NMDS.BC.P.1440.only.FB.chull <- c(NMDS.BC.P.1440.only.FB.chull, NMDS.BC.P.1440.only.FB.chull[1])


pdf(paste0("BETAPLOTS_P/NMDS_",metric,"_TimebyMorph.pdf"),pointsize = 14, width = 7, height = 4)
par(mfrow= c(1,4), oma = c(4,4,4,6))
par(mar = c(4,0,4,0))
plot(NMDS.BC.P.20.only
     , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.P.20.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.P.20.only.CR[NMDS.BC.P.20.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.20.only.BL[NMDS.BC.P.20.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.20.only.FB[NMDS.BC.P.20.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.BC.P.60.only
     , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.BC.P.60.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.P.60.only.CR[NMDS.BC.P.60.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.60.only.BL[NMDS.BC.P.60.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.60.only.FB[NMDS.BC.P.60.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.BC.P.180.only
     , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.180.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p =  ",ANOVA.BC.P.180.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.P.180.only.CR[NMDS.BC.P.180.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.180.only.BL[NMDS.BC.P.180.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.180.only.FB[NMDS.BC.P.180.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.BC.P.360.only
     , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.P.360.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.P.360.only.CR[NMDS.BC.P.360.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.360.only.BL[NMDS.BC.P.360.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.360.only.FB[NMDS.BC.P.360.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.BC.P.720.only
     , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.P.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.P.720.only.CR[NMDS.BC.P.720.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.720.only.BL[NMDS.BC.P.720.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.720.only.FB[NMDS.BC.P.720.only.FB.chull,]
      , col = MorphColours[3])
par(mar = c(4,0,4,0))
plot(NMDS.BC.P.1440.only
     , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.1440.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
     , cex = 2
     , cex.lab = 2
     , cex.main = 2
)
lines(NMDS.BC.P.1440.only.CR[NMDS.BC.P.1440.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.1440.only.BL[NMDS.BC.P.1440.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.1440.only.FB[NMDS.BC.P.1440.only.FB.chull,]
      , col = MorphColours[3])
#STOP
dev.off()

############ COMBO DISP BETA ################
# Disp and beta through time combined

# xvalues <- log(FB.BL.BC.P.ALL.mean[,1])
xvalues <- as.character(FB.BL.BC.P.ALL.mean[,1])

pdf(paste0("BETAPLOTS_P/COMBO_P_dispbeta_",metric,".pdf"), pointsize = 14, width = 14, height = 8)
par(fig = c(0,0.8,0.23,1))
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = ''
     , ylab = "Distance (Bray-Curtis)"
     , ylim = ylimits
     , xaxt = 'n')
axis(side = 1
     , at = c(1,2,3,4,5,6)
     , labels = c("20 min","1 h","3 h","6 h","12 h","24 h")
     , las = 2)
points(FB.BL.BC.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "purple"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*0.99
       , x1 = c(1,2,3,4,5,6)*0.99
       , y0 = c(FB.BL.BC.P.ALL.mean[,2] - FB.BL.BC.P.ALL.sd[,2])
       , y1 = c(FB.BL.BC.P.ALL.mean[,2] + FB.BL.BC.P.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(FB.CR.BC.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "darkgreen"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4,5,6)*1
       , x1 = c(1,2,3,4,5,6)*1
       , y0 = c(FB.CR.BC.P.ALL.mean[,2] - FB.CR.BC.P.ALL.sd[,2])
       , y1 = c(FB.CR.BC.P.ALL.mean[,2] + FB.CR.BC.P.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
points(CR.BL.BC.P.ALL.mean[,2] ~ c(1,2,3,4,5,6)
       , type = 'l'
       , col = "grey"
       , lwd = 2
       , lty = 1)
arrows(x0 = c(1,2,3,4)*1.01
       , x1 = c(1,2,3,4)*1.01
       , y0 = c(CR.BL.BC.P.ALL.mean[,2] - CR.BL.BC.P.ALL.sd[,2])
       , y1 = c(CR.BL.BC.P.ALL.mean[,2] + CR.BL.BC.P.ALL.sd[,2])
       , angle = 90
       , code = 3
       , length = 0.03)
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
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.05,0.158333,0,0.3), new = TRUE)
plot(NMDS.BC.P.20.only
     # , main = "20 Minutes"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.20.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.P.20.only$aov.tab[1]$Df[1],",",ANOVA.BC.P.20.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.P.20.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.P.20.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.P.20.only.CR[NMDS.BC.P.20.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.20.only.BL[NMDS.BC.P.20.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.20.only.FB[NMDS.BC.P.20.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.178333,0.28666,0,0.3), new = TRUE)
plot(NMDS.BC.P.60.only
     # , main = "1 Hour"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.60.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.P.60.only$aov.tab[1]$Df[1],",",ANOVA.BC.P.60.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.P.60.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.P.60.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.P.60.only.CR[NMDS.BC.P.60.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.60.only.BL[NMDS.BC.P.60.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.60.only.FB[NMDS.BC.P.60.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.30666,0.414999,0,0.3), new = TRUE)
plot(NMDS.BC.P.180.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.180.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.P.180.only$aov.tab[1]$Df[1],",",ANOVA.BC.P.180.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.P.180.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.P.180.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.P.180.only.CR[NMDS.BC.P.180.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.180.only.BL[NMDS.BC.P.180.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.180.only.FB[NMDS.BC.P.180.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.434999,0.543333,0,0.3), new = TRUE)
plot(NMDS.BC.P.360.only
     # , main = "6 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.360.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = 'n'
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.P.360.only$aov.tab[1]$Df[1],",",ANOVA.BC.P.360.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.P.360.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.P.360.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.P.360.only.CR[NMDS.BC.P.360.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.360.only.BL[NMDS.BC.P.360.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.360.only.FB[NMDS.BC.P.360.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.565555,0.671666,0,0.3), new = TRUE)
plot(NMDS.BC.P.720.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.720.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.P.720.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.P.720.only$aov.tab[1]$Df[1],",",ANOVA.BC.P.720.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.P.720.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.P.720.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.P.720.only.CR[NMDS.BC.P.720.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.720.only.BL[NMDS.BC.P.720.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.720.only.FB[NMDS.BC.P.720.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(3,0,2,0), fig = c(0.691666,0.79999,0,0.3), new = TRUE)
plot(NMDS.BC.P.1440.only
     # , main = "12 Hours"
     , pch = 21
     , col = "black"
     , bg = MorphColours[factor(MF.P.BC.P.1440.only$Morph)]
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = paste0("p = ", ANOVA.BC.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
     , ylab = ''
)
title(xlab = paste0("Df = ", ANOVA.BC.P.1440.only$aov.tab[1]$Df[1],",",ANOVA.BC.P.1440.only$aov.tab[1]$Df[3])
      , line = 0.25
      , cex.axis = 0.5)
title(xlab = substitute(paste(R^2 == X), list(X = format(ANOVA.BC.P.1440.only$aov.tab[5]$R2[1], digits = 2))) 
      , line= 1
      , cex.axis = 0.5)
title(xlab = paste0("p = ", ANOVA.BC.P.1440.only$aov.tab[6]$`Pr(>F)`[1])
      , line = 1.75
      , cex.axis = 0.5)
lines(NMDS.BC.P.1440.only.CR[NMDS.BC.P.1440.only.CR.chull,]
      , col = MorphColours[1])
lines(NMDS.BC.P.1440.only.BL[NMDS.BC.P.1440.only.BL.chull,]
      , col = MorphColours[2])
lines(NMDS.BC.P.1440.only.FB[NMDS.BC.P.1440.only.FB.chull,]
      , col = MorphColours[3])
par(oma = c(1,0.1,1,0.1), mar = c(2,0,2,0), fig = c(0.775,1,0,0.3), new = TRUE)
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

allindividualTests <- c("ANOVA.BC.P.20.only"
                        ,"ANOVA.BC.P.60.only"
                        ,"ANOVA.BC.P.180.only"
                        ,"ANOVA.BC.P.360.only"
                        ,"ANOVA.BC.P.720.only"
                        ,"ANOVA.BC.P.1440.only"
                        
)
for (i in allindividualTests) {
  # print(get(i))
  capture.output(get(i), file = paste0("./BETAPLOTS_P/individualtests/",i,".txt"))
}

########### BETADISP#############
betadisp.P.BC.FB <- betadisp.BC.P.time$distances[grep("FB", names(betadisp.BC.P.time$distances))]
betadisp.P.BC.BL <- betadisp.BC.P.time$distances[grep("BL", names(betadisp.BC.P.time$distances))]
betadisp.P.BC.CR <- betadisp.BC.P.time$distances[grep("CR", names(betadisp.BC.P.time$distances))]

MF.P.FB <- MF.P.morphkeep[grep("FB", rownames(MF.P.morphkeep)),]
MF.P.BL <- MF.P.morphkeep[grep("BL", rownames(MF.P.morphkeep)),]
MF.P.CR <- MF.P.morphkeep[grep("CR", rownames(MF.P.morphkeep)),]

betadisp.FB.P.agg <- aggregate(betadisp.P.BC.FB, by = list(MF.P.FB$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.BL.P.agg <- aggregate(betadisp.P.BC.BL, by = list(MF.P.BL$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )
betadisp.CR.P.agg <- aggregate(betadisp.P.BC.CR, by = list(MF.P.CR$Time), FUN = function(x) {c(mean = mean(x),sd = sd(x))} )

betadisp.BC.time.forstat <- cbind(betadisp.BC.time$distances, MF.morphkeep[unlist(lapply(names(betadisp.BC.time$distances), function(x) {grep(paste0("^",x,"$"), rownames(MF.morphkeep))})),])
colnames(betadisp.BC.time.forstat)[1] <- c("Distance")
ANOVA.betadisp.BC <- anova(lm(Distance ~ Time*Morph, data = betadisp.BC.time.forstat))
capture.output(ANOVA.betadisp.BC, file = "./BETAPLOTS_P/ANOVA.betadisp.BC.txt")

ylimits <- c(0.15,0.7)
pdf(paste0("./BETAPLOTS_P/BetaDisp_P_",metric,"_eachmorph.pdf"),pointsize = 14)
plot(xvalues, NULL
     , main = "Dispersion of morphologies across time"
     , xlab = 'Time'
     , ylab = "Distance (Bray-Curtis)"
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

anova.betadisp.BC.morph <- anova(betadisp.BC.morph)
anova.betadisp.BC.time <- anova(betadisp.BC.time)
anova.betadisp.P.BC.morph <- anova(betadisp.BC.P.morph)
anova.betadisp.P.BC.time <- anova(betadisp.BC.P.time)

# This is PERMADISP-- don't need a table, will quote in text
capture.output(anova.betadisp.BC.morph, file = "BETAPLOTS_H/anova.betadisp.BC.morph.txt")
capture.output(anova.betadisp.BC.time, file = "BETAPLOTS_H/anova.betadisp.BC.time.txt")
capture.output(anova.betadisp.P.BC.morph, file = "BETAPLOTS_P/anova.betadisp.BC.P.morph.txt")
capture.output(anova.betadisp.P.BC.time, file = "BETAPLOTS_P/anova.betadisp.BC.P.time.txt")

### NOW DO EACH TIME POINT ##
## This is a table with both Hakai and PM for the supplementary figures
permdisp.morphology.across.time <- matrix(ncol = 8, nrow = 2)
rownames(permdisp.morphology.across.time) <- c("P","H")
colnames(permdisp.morphology.across.time) <- c("20","60","180","360","720","1440","5760", "Overall")
for (t in c("20","60","180","360","720","1440")) {
  assign(paste0("betadisp.P.",t), betadisper(dist(get(paste0("dm.BC.P.",t))), group = get(paste0("MF.P.BC.P.",t,".only"))$Morph))
  assign(paste0("anova.betadisp.P.",t), anova(get(paste0("betadisp.P.",t))))
  ptemp <- get(paste0("anova.betadisp.P.",t))$`Pr(>F)`[1]
  ftemp <- get(paste0("anova.betadisp.P.",t))$`F value`[1]
  dftemp <- paste0(get(paste0("anova.betadisp.P.",t))$Df[1],",",get(paste0("anova.betadisp.P.",t))$Df[3])
  
  toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
  permdisp.morphology.across.time["P", paste0(t)] <- toPaste
}
for (t in c("20","60","360","720","5760")) {
  assign(paste0("betadisp.",t), betadisper(dist(get(paste0("dm.BC.",t))), group = get(paste0("MF.BC.",t,".only"))$Morph))
  assign(paste0("anova.betadisp.",t), anova(get(paste0("betadisp.",t))))
  ptemp <- get(paste0("anova.betadisp.",t))$`Pr(>F)`[1]
  ftemp <- get(paste0("anova.betadisp.",t))$`F value`[1]
  dftemp <- paste0(get(paste0("anova.betadisp.",t))$Df[1],",",get(paste0("anova.betadisp.",t))$Df[3])
  
  toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
  permdisp.morphology.across.time["H", paste0(t)] <- toPaste
  
}
# Get overall
ptemp <- anova.betadisp.P.BC.morph$`Pr(>F)`[1]
ftemp <- anova.betadisp.P.BC.morph$`F value`[1]
dftemp <- paste0(anova.betadisp.P.BC.morph$Df[1],",",anova.betadisp.P.BC.morph$Df[12])
toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
permdisp.morphology.across.time["P","Overall"] <- toPaste

ptemp <- anova.betadisp.BC.morph$`Pr(>F)`[1]
ftemp <- anova.betadisp.BC.morph$`F value`[1]
dftemp <- paste0(anova.betadisp.BC.morph$Df[1],",",anova.betadisp.BC.morph$Df[2])
toPaste <- paste0(signif(ptemp,2), " (f=",round(ftemp,2),",Df=",dftemp,")")
permdisp.morphology.across.time["H","Overall"] <- toPaste

colnames(permdisp.morphology.across.time) <- c("20 minutes","1 hour","3 hours","6 hours","12 hours","1 day","4 days", "Overall")

tempfile1 <- permdisp.morphology.across.time
for (r in 1:nrow(tempfile1)) {
  for (c in 1:ncol(tempfile1)) {
    if (is.na(tempfile1[r,c])) {
      permdisp.morphology.across.time[r,c] <- "-"
    }}}


# Make double header
permdisp.morphology.across.time.BC <- permdisp.morphology.across.time
rownames(permdisp.morphology.across.time.BC) <- c("Reed Point","Hakai")

capture.output(xtable(permdisp.morphology.across.time.BC, digits = NULL), file = paste0("BETAPLOTS_LATEX/permdisp.morph.across.time.",metric,".txt"))

### MAKE BETA DIV TABLES-- metrics separately but H and P together; extras will go in supp

anova.morphology.across.time <- matrix(ncol = 8, nrow = 2)
colnames(anova.morphology.across.time) <- c("20","60","180","360","720","1440","5760", "Overall")
rownames(anova.morphology.across.time) <- c("P","H")
for (t in c("20","60","180","360","720","1440")) {
  ptemp <- get(paste0("ANOVA.BC.P.",t,".only"))$aov.tab$`Pr(>F)`[1]
  rtemp <- get(paste0("ANOVA.BC.P.",t,".only"))$aov.tab$`R2`[1]
  dftemp <- paste0(get(paste0("ANOVA.BC.P.",t,".only"))$aov.tab$`Df`[1],",",get(paste0("ANOVA.BC.P.",t,".only"))$aov.tab$`Df`[3])
  toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
  
  anova.morphology.across.time["P",paste0(t)] <- toPaste
}
for (t in c("20","60","360","720","5760")) {
  ptemp <- get(paste0("ANOVA.BC.",t,".only"))$aov.tab$`Pr(>F)`[1]
  rtemp <- get(paste0("ANOVA.BC.",t,".only"))$aov.tab$`R2`[1]
  dftemp <- paste0(get(paste0("ANOVA.BC.",t,".only"))$aov.tab$`Df`[1],",",get(paste0("ANOVA.BC.",t,".only"))$aov.tab$`Df`[3])
  toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
  
  anova.morphology.across.time["H",paste0(t)] <- toPaste
}
# Do overall P
ptemp <- ANOVA.BC.P.morphtime$aov.tab$`Pr(>F)`[2]
rtemp <- ANOVA.BC.P.morphtime$aov.tab$`R2`[2]
dftemp <- paste0(ANOVA.BC.P.morphtime$aov.tab$`Df`[2],",",ANOVA.BC.P.morphtime$aov.tab$`Df`[5])
toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
anova.morphology.across.time["P","Overall"] <- toPaste

# Do overall H
ptemp <- ANOVA.BC.morphtime$aov.tab$`Pr(>F)`[2]
rtemp <- ANOVA.BC.morphtime$aov.tab$`R2`[2]
dftemp <- paste0(ANOVA.BC.morphtime$aov.tab$`Df`[2],",",ANOVA.BC.morphtime$aov.tab$`Df`[5])
toPaste <- paste0(signif(ptemp,2), " (R^2=",round(rtemp,2),",Df=",dftemp,")")
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
anova.morphology.across.time.BC <- anova.morphology.across.time
rownames(anova.morphology.across.time.BC) <- c("Reed Point","Hakai")

capture.output(xtable(anova.morphology.across.time.BC, digits = NULL), file = paste0("BETAPLOTS_LATEX/anova.morph.across.time.",metric,".txt"))

### DO FB/CR/BL TEST FOR EACH
# Make table
pairwiseAdonis.all <- matrix(ncol = 7,nrow = 6)
colnames(pairwiseAdonis.all) <- c("20","60","180","360","720","1440", "5760")
rownames(pairwiseAdonis.all) <- c("pFB:BL","pFB:CR","pBL:CR","hFB:BL","hFB:CR","hBL:CR" )
listMorphs <- c("FB","BL","CR")
for (m in 1:(length(listMorphs)-1)) {
  for (n in (m+1):length(listMorphs)) {
    for (t in c("20","60","180","360","720","1440")) {
      tempMF <- get(paste0("MF.P.BC.P.",t,".only"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), get(paste0("MF.P.BC.P.",t,".only"))$Morph),]
      tempDM <- get(paste0("dm.BC.P.",t))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.BC.P.",t)))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.BC.P.",t))))]
      tempAdonis <- adonis(tempDM ~ Morph, data = tempMF, by = "marginal")
      newP <- p.adjust(tempAdonis$aov.tab$`Pr(>F)`[1], method = "fdr", n = 3)
      toPaste <- paste0( newP
                         ," (R^2 = ", round(tempAdonis$aov.tab$R2[1],digits = 2)
                         ,", Df = ", paste0(tempAdonis$aov.tab$Df[1],",",tempAdonis$aov.tab$Df[3]) 
                         , ")")
      pairwiseAdonis.all[paste0("p",listMorphs[m],":",listMorphs[n]),paste0(t)] <- toPaste
    }
    for (t in c("20","60","360","720","5760")) {
      tempMF <- get(paste0("MF.BC.",t,".only"))[grep(paste0(listMorphs[m],"|",listMorphs[n]), get(paste0("MF.BC.",t,".only"))$Morph),]
      tempDM <- get(paste0("dm.BC.",t))[grep(paste0(listMorphs[m],"|",listMorphs[n]), rownames(get(paste0("dm.BC.",t)))),grep(paste0(listMorphs[m],"|",listMorphs[n]), colnames(get(paste0("dm.BC.",t))))]
      tempAdonis <- adonis(tempDM ~ Morph, data = tempMF, by = "marginal")
      toPaste <- paste0( tempAdonis$aov.tab$`Pr(>F)`[1]
                         ," (R^2 = ", round(tempAdonis$aov.tab$R2[1],digits = 2)
                         ,", Df = ", paste0(tempAdonis$aov.tab$Df[1],",",tempAdonis$aov.tab$Df[3])
                         , ")")
      pairwiseAdonis.all[paste0("h",listMorphs[m],":",listMorphs[n]),paste0(t)] <- toPaste
    }
  }
}

# Get rid of NAs
for (r in 1:nrow(pairwiseAdonis.all)) {
  for (c in 1:ncol(pairwiseAdonis.all)) {
    if (is.na(pairwiseAdonis.all[r,c])) {
      pairwiseAdonis.all[r,c] <- "-"
    }
  }
}

# Change Rownames
pairwiseAdonis.all.BC <- cbind(c("FB:BL","FB:CR","BL:CR","FB:BL","FB:CR","BL:CR"), pairwiseAdonis.all)
rownames(pairwiseAdonis.all.BC) <- c("Reed Point", " ","  ","Hakai","   ","    ")

capture.output(xtable(pairwiseAdonis.all.BC), file = paste0("BETAPLOTS_LATEX/pairwiseAdonis.",metric,".txt"))


############ PRINTING LATEX-DEP STATS ############
# COMBINING BETA DISP
permdisp.morphology.across.time.MASTER <- matrix(ncol = 8, nrow = 6)
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

capture.output(xtable(permdisp.morphology.across.time.MASTER.edit), file = paste0("BETAPLOTS_LATEX/permdisp.morphology.across.time.MASTER.txt"))





# COMBINING BETA ANOVAs
anova.morphology.across.time.MASTER <- matrix(ncol = 8, nrow = 6)
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

capture.output(xtable(anova.morphology.across.time.MASTER.edit), file = paste0("BETAPLOTS_LATEX/anova.morphology.across.time.MASTER.txt"))


# Combine all pairwise Adonises too? No. Too many














# "\newcommand{\newCommandName}{text to insert}"
# BETA ANOVA RESULTS
# 
# BetaDivResultsToPrint <- matrix(ncol = 1)
# BetaDivResultsToPrint <- matrix()
# for (i in c("UWUF","WUF","BC")) {
#   # Assign all names
#   assign(paste0("aov.",i,".H.Time.p"), format(get(paste0("ANOVA.",i,".morphtime"))$aov.tab$`Pr(>F)`[1], digits= 2))
#   assign(paste0("aov.",i,".H.Time.r"), format(get(paste0("ANOVA.",i,".morphtime"))$aov.tab$R2[1], digits = 2))
#   assign(paste0("aov.",i,".H.Morph.p"), format(get(paste0("ANOVA.",i,".morphtime"))$aov.tab$`Pr(>F)`[2], digits = 2))
#   assign(paste0("aov.",i,".H.Morph.r"), format(get(paste0("ANOVA.",i,".morphtime"))$aov.tab$R2[2], digits = 2))
#   assign(paste0("aov.",i,".H.TimeMorph.p"), format(get(paste0("ANOVA.",i,".morphtime"))$aov.tab$`Pr(>F)`[3], digits = 2))
#   assign(paste0("aov.",i,".H.TimeMorph.r"), format(get(paste0("ANOVA.",i,".morphtime"))$aov.tab$R2[3], digits = 2))
#   
#   assign(paste0("aov.",i,".P.Time.p"), format(get(paste0("ANOVA.",i,".P.morphtime"))$aov.tab$`Pr(>F)`[1], digits= 2))
#   assign(paste0("aov.",i,".P.Time.r"), format(get(paste0("ANOVA.",i,".P.morphtime"))$aov.tab$R2[1], digits = 2))
#   assign(paste0("aov.",i,".P.Morph.p"), format(get(paste0("ANOVA.",i,".P.morphtime"))$aov.tab$`Pr(>F)`[2], digits = 2))
#   assign(paste0("aov.",i,".P.Morph.r"), format(get(paste0("ANOVA.",i,".P.morphtime"))$aov.tab$R2[2], digits = 2))
#   assign(paste0("aov.",i,".P.TimeMorph.p"), format(get(paste0("ANOVA.",i,".P.morphtime"))$aov.tab$`Pr(>F)`[3], digits = 2))
#   assign(paste0("aov.",i,".P.TimeMorph.r"), format(get(paste0("ANOVA.",i,".P.morphtime"))$aov.tab$R2[3], digits = 2))
#   
#   tempFile <- c(  paste0("aov.",i,".H.Time.p")
#                 , paste0("aov.",i,".H.Time.r")
#                 , paste0("aov.",i,".H.Morph.p")
#                 , paste0("aov.",i,".H.Morph.r")
#                 , paste0("aov.",i,".H.TimeMorph.p")
#                 , paste0("aov.",i,".H.TimeMorph.r")
# 
#                 , paste0("aov.",i,".P.Time.p")
#                 , paste0("aov.",i,".P.Time.r")
#                 , paste0("aov.",i,".P.Morph.p")
#                 , paste0("aov.",i,".P.Morph.r")
#                 , paste0("aov.",i,".P.TimeMorph.p")
#                 , paste0("aov.",i,".P.TimeMorph.r")
#   )
#   
# 
#   for (i in tempFile) {
#     tempName <- gsub("[.]","", i)
#     #"\newcommand{\newCommandName}{text to insert}"
#     tempInsert <- paste0("\\","newcommand{","\\",tempName,"}{",get(i),"}")
#      BetaDivResultsToPrint <- rbind(BetaDivResultsToPrint, paste0("\\newcommand{\\",tempName,"}{",get(i),"}"))
#   }
# }
# write.matrix(BetaDivResultsToPrint,file = "./BETAPLOTS_LATEX/BetaDivStats.txt")
# 



# BETA INDIVIDUAL RESULTS
# 
# ANOVA.BC.5760.only <- ANOVA.BC.5760
# ANOVA.UWUF.5760.only <- ANOVA.UWUF.5760
# ANOVA.WUF.5760.only <- ANOVA.WUF.5760
# 
# BetaDivSepToPrint <- matrix(ncol = 1)
# translated <- cbind(c("20","60","180","360","720","1440","5760"), c("A","B","C","D","E","F","G"))
# for (i in c("UWUF","WUF","BC")) {
#   tempFile <- c()
#   for (j in c("20","60","360","720","5760")) {
#     TRANS <- translated[grep(paste0("^",j,"$"), translated[,1]),2]
#     # Assign all names
#     assign(paste0("aov.",i,".H.",TRANS,".p"), format(get(paste0("ANOVA.",i,".",j,".only"))$aov.tab$`Pr(>F)`[1], digits= 2))
#     assign(paste0("aov.",i,".H.",TRANS,".r"), format(get(paste0("ANOVA.",i,".",j,".only"))$aov.tab$R2[1], digits = 2))
#     tempFile <- c(  tempFile
#                     , paste0("aov.",i,".H.",TRANS,".p")
#                     , paste0("aov.",i,".H.",TRANS,".r")
#     )
#     
#   }
#   for (l in c("20","60","180","360","720","1440")) {
#     TRANS <- translated[grep(paste0("^",l,"$"), translated[,1]),2]
#     # Assign all names
#     assign(paste0("aov.",i,".P.",TRANS,".p"), format(get(paste0("ANOVA.",i,".P.",l,".only"))$aov.tab$`Pr(>F)`[1], digits= 2))
#     assign(paste0("aov.",i,".P.",TRANS,".r"), format(get(paste0("ANOVA.",i,".P.",l,".only"))$aov.tab$R2[1], digits = 2))
#     tempFile <- c(  tempFile
#                     , paste0("aov.",i,".P.",TRANS,".p")
#                     , paste0("aov.",i,".P.",TRANS,".r"))
#   }
#   
# 
#   for (k in tempFile) {
#     tempName <- gsub("[.]","", k)
#     #"\newcommand{\newCommandName}{text to insert}"
#     BetaDivSepToPrint <- rbind(BetaDivSepToPrint, paste0("\\newcommand{\\",tempName,"}{",get(k),"}"))
#   }
#   
# }
# write.matrix(BetaDivSepToPrint,file = "./BETAPLOTS_LATEX/BetaSepStats.txt")




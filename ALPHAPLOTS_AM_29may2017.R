#!/bin/Rscript

# setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/ALPHA_R")
library("car")
library("nlme")
library("optparse")
library("xtable")
# library("gridExtra")
########################### OPT PARSE #################################

option_list = list(
  make_option(c("-m", "--mappingFP"), type="character",
              help="Mapping file with alpha data added"),
  make_option(c("-a", "--alphaNames"), type="character",
              help="Comma separated list of alpha metrics used")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

MPFP = opt$mappingFP
alphaNames = opt$alphaNames
alphaList <- unlist(strsplit(alphaNames, ","))
########################### FOR TESTING #################################

setwd('/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis')
MPFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/TEMP_frombotclust/MF_withalpha.txt'
alphaNames = 'chao1_even_4000_normalized_alpha,PD_whole_tree_even_4000_normalized_alpha,observed_otus_even_4000_normalized_alpha'
alphaList <- unlist(strsplit(alphaNames, ","))

# setwd('/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis')
# MPFP <- '/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_AM/1_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha.txt'
# alphaNames <-  'chao1_even_4000_normalized_alpha,PD_whole_tree_even_4000_normalized_alpha,observed_otus_even_4000_normalized_alpha'
# alphaList <- unlist(strsplit(alphaNames, ","))

########################### LOAD AND FILTER #################################
# Load files 
MF.Alpha.all <- read.delim(paste0(MPFP)
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE)

# Delete samples that don't have alpha value

MF.Alpha.all <- MF.Alpha.all[-grep("N/A", MF.Alpha.all$PD_whole_tree_even_4000_alpha),]



# Grep out only samples
storeOnlySeaweed <- c()
storeAllTypes <- c()
for (i in 1:nrow(MF.Alpha.all)) {
  if (MF.Alpha.all[i,'Morph'] %in% c("CR","BL","FB")) {
    storeOnlySeaweed <- c(storeOnlySeaweed, i)
  } 
  if (MF.Alpha.all[i,'Morph'] %in% c("CR","BL","FB","W","H2O")) {
    storeAllTypes <- c(storeAllTypes, i)
  }
}

MF.Alpha.morphonly <- MF.Alpha.all[storeOnlySeaweed,]
MF.Alpha.morphwater <- MF.Alpha.all[storeAllTypes,]

# Make Hakai and Port moody
storeOnlyH <- c()
storeOnlyP <- c()
for (i in 1:nrow(MF.Alpha.morphonly)) {
  if (MF.Alpha.morphonly[i,'Type'] %in% c("H")) {
    storeOnlyH <- c(storeOnlyH, i)
  } 
  if (MF.Alpha.morphonly[i,'Type'] %in% c("P")) {
    storeOnlyP <- c(storeOnlyP, i)
  }
}

MF.Alpha.morphonly.H <- MF.Alpha.morphonly[storeOnlyH,]
MF.Alpha.morphonly.P <- MF.Alpha.morphonly[storeOnlyP,]


# Replace headers to simple format

colnames(MF.Alpha.morphonly.H) <- gsub("chao1_even_4000_alpha", "chao1", colnames(MF.Alpha.morphonly.H))
colnames(MF.Alpha.morphonly.P) <- gsub("chao1_even_4000_alpha", "chao1", colnames(MF.Alpha.morphonly.P))

colnames(MF.Alpha.morphonly.H) <- gsub("PD_whole_tree_even_4000_alpha", "PD_whole_tree", colnames(MF.Alpha.morphonly.H))
colnames(MF.Alpha.morphonly.P) <- gsub("PD_whole_tree_even_4000_alpha", "PD_whole_tree", colnames(MF.Alpha.morphonly.P))

colnames(MF.Alpha.morphonly.H) <- gsub("observed_otus_even_4000_alpha", "observed_otus", colnames(MF.Alpha.morphonly.H))
colnames(MF.Alpha.morphonly.P) <- gsub("observed_otus_even_4000_alpha", "observed_otus", colnames(MF.Alpha.morphonly.P))

# Makes all numeric


MF.Alpha.morphonly.H[,'chao1'] <- as.numeric(MF.Alpha.morphonly.H[,'chao1'])
MF.Alpha.morphonly.H[,'PD_whole_tree'] <- as.numeric(MF.Alpha.morphonly.H[,'PD_whole_tree'])
MF.Alpha.morphonly.H[,'observed_otus'] <- as.numeric(MF.Alpha.morphonly.H[,'observed_otus'])


MF.Alpha.morphonly.P[,'chao1'] <- as.numeric(MF.Alpha.morphonly.P[,'chao1'])
MF.Alpha.morphonly.P[,'PD_whole_tree'] <- as.numeric(MF.Alpha.morphonly.P[,'PD_whole_tree'])
MF.Alpha.morphonly.P[,'observed_otus'] <- as.numeric(MF.Alpha.morphonly.P[,'observed_otus'])


########################### MKDIR #################################
system('mkdir ALPHAPLOTS_H')
system('mkdir ALPHAPLOTS_P')
system('mkdir ALPHAPLOTS_latex')
system('mkdir ALPHA_compare')

##################### REP COUNT ######################
# Count replicates of each treatment

H.Replicates <- table(MF.Alpha.morphonly.H$ColRep)
P.Replicates <- table(MF.Alpha.morphonly.P$ColRep)

replicateTable <- matrix(ncol = 7, nrow = 6)
colnames(replicateTable) <- c("20","60","180","360","720","1440","5760")
rownames(replicateTable) <- c("PFB","PBL","PCR","HFB","HBL","HCR")
for (m in c("FB","BL","CR")) {
  for (t in c("20","60","180","360","720","1440")) {
    replicateTable[paste0("P",m), paste0(t)] <-  P.Replicates[grep(paste0(m,t), names(P.Replicates))]
  }
  for (t in c("20","60","360","720","5760")) {
    replicateTable[paste0("H",m), paste0(t)] <-  H.Replicates[grep(paste0(m,t), names(H.Replicates))]
  }
  
}

for (r in 1:nrow(replicateTable)) {
  for (c in 1:ncol(replicateTable)) {
    if (is.na(replicateTable[r,c])) {
      replicateTable[r,c] <- "-"
    }
  }
}

newRepTable <- cbind(c("Reed Point","","","Hakai","",""),c("FB","BL","CR","FB","BL","CR"), replicateTable)
colnames(newRepTable) <- c("Site","Morphology","20 m","1 h","3 h","6 h","12 h","1 d","4 d")

capture.output(print(xtable(newRepTable), include.rownames=FALSE), file = "ALPHAPLOTS_latex/Replicates.txt")

########### ****CHAO1**** ##############
metric <- 'chao1'

########### Stats ##############

H.morph.lm <- lm(chao1 ~ Time*Morph, data = MF.Alpha.morphonly.H)
anova.H.morph.lm <- Anova(H.morph.lm, type = 3)
capture.output(xtable(anova.H.morph.lm), file = paste0("ALPHAPLOTS_latex/anova_Hakai_", metric, ".txt"))
P.morph.lm <- lm(chao1 ~ Time*Morph, data = MF.Alpha.morphonly.P)
anova.P.morph.lm <- Anova(P.morph.lm, type = 3)
capture.output(xtable(anova.P.morph.lm), file = paste0("ALPHAPLOTS_latex/anova_PM_", metric, ".txt"))

# Separated out Stats- Hakai
overall.anova.morphology.through.time.H <- matrix(nrow = 5, ncol = 3)
colnames(overall.anova.morphology.through.time.H) <- c("p-value","F-value","Df")
rownames(overall.anova.morphology.through.time.H) <- c("20","60","360","720","5760")
pairwise.ttest.morphology.through.time.H <- matrix(nrow = 5, ncol = 3)
colnames(pairwise.ttest.morphology.through.time.H) <- c("FB:BL","FB:CR","BL:CR")
rownames(pairwise.ttest.morphology.through.time.H) <- c("20","60","360","720","5760")

for (t in c("20","60","360","720","5760")) {
  
  # Get overall anova at each time point and make into a table; print this table as latex 
  assign(paste0("MF.Alpha.H.",t), MF.Alpha.morphonly.H[grep(paste0("^",t,"$"), MF.Alpha.morphonly.H$Time),])
  assign(paste0("H.morph.",t,".lm"), lm(chao1 ~ Morph, data = get(paste0("MF.Alpha.H.",t))))
  assign(paste0("anova.H.morph.",t,".lm"), anova(get(paste0("H.morph.",t,".lm"))))
  overall.anova.morphology.through.time.H[paste0(t),"p-value"] <- get(paste0("anova.H.morph.",t,".lm"))$`Pr(>F)`[1]
  overall.anova.morphology.through.time.H[paste0(t), "F-value"] <- get(paste0("anova.H.morph.",t,".lm"))$`F value`[1]
  overall.anova.morphology.through.time.H[paste0(t), "Df"] <- get(paste0("anova.H.morph.",t,".lm"))$Df[1]
 
  # Now, do pairwise tests
  assign(paste0("pairwisettest.H.morph.",t) , pairwise.t.test(x = get(paste0("MF.Alpha.H.",t))[,paste0(metric)]
                                              , g = get(paste0("MF.Alpha.H.",t))[,"Morph"]
                                              , p.adjust.method = "BH"
  ))
  pairwise.ttest.morphology.through.time.H[paste0(t),"FB:BL"] <- get(paste0("pairwisettest.H.morph.",t))$p.value[2,1]
  pairwise.ttest.morphology.through.time.H[paste0(t),"FB:CR"] <- get(paste0("pairwisettest.H.morph.",t))$p.value[2,2]
  pairwise.ttest.morphology.through.time.H[paste0(t),"BL:CR"] <- get(paste0("pairwisettest.H.morph.",t))$p.value[1,1]
  
}
# Change rownames to something nice
rownames(overall.anova.morphology.through.time.H) <- c("20 min","1 h","6 h","12 h","4 d")
rownames(pairwise.ttest.morphology.through.time.H) <- c("20 min","1 h","6 h","12 h","4 d")

# Round the numbers to make sure all are appropriate length
tempfile <- overall.anova.morphology.through.time.H
for (r in 1:nrow(tempfile)) {
  for (c in 1:ncol(tempfile)) {
    if (tempfile[r,c] < 0.001) {
      overall.anova.morphology.through.time.H[r,c] <- "<0.001"
    } else {
      overall.anova.morphology.through.time.H[r,c] <- round(tempfile[r,c], digits = 3)
    }
  }
}

tempfile2 <- pairwise.ttest.morphology.through.time.H
for (r in 1:nrow(tempfile2)) {
  for (c in 1:ncol(tempfile2)) {
    if (tempfile2[r,c] < 0.001) {
      pairwise.ttest.morphology.through.time.H[r,c] <- "<0.001"
    } else {
      pairwise.ttest.morphology.through.time.H[r,c] <- round(tempfile2[r,c], digits = 3)
    }
  }
}

capture.output(xtable(overall.anova.morphology.through.time.H), file = paste0("ALPHAPLOTS_latex/overall_anova_morphology_through_time_H",metric,".txt"))
capture.output(xtable(pairwise.ttest.morphology.through.time.H), file = paste0("ALPHAPLOTS_latex/pairwise_ttest_morphology_through_time_H",metric,".txt"))




# Separated out Stats- PM
overall.anova.morphology.through.time.P <- matrix(nrow = 6, ncol = 3)
colnames(overall.anova.morphology.through.time.P) <- c("p-value","F-value","Df")
rownames(overall.anova.morphology.through.time.P) <- c("20","60","180","360","720","1440")
pairwise.ttest.morphology.through.time.P <- matrix(nrow = 6, ncol = 3)
colnames(pairwise.ttest.morphology.through.time.P) <- c("FB:BL","FB:CR","BL:CR")
rownames(pairwise.ttest.morphology.through.time.P) <- c("20","60","180","360","720","1440")

for (t in c("20","60","180","360","720","1440")) {
  
  # Get overall anova at each time point and make into a table; print this table as latex 
  assign(paste0("MF.Alpha.P.",t), MF.Alpha.morphonly.P[grep(paste0("^",t,"$"), MF.Alpha.morphonly.P$Time),])
  assign(paste0("P.morph.",t,".lm"), lm(chao1 ~ Morph, data = get(paste0("MF.Alpha.P.",t))))
  assign(paste0("anova.P.morph.",t,".lm"), anova(get(paste0("P.morph.",t,".lm"))))
  overall.anova.morphology.through.time.P[paste0(t),"p-value"] <- get(paste0("anova.P.morph.",t,".lm"))$`Pr(>F)`[1]
  overall.anova.morphology.through.time.P[paste0(t), "F-value"] <- get(paste0("anova.P.morph.",t,".lm"))$`F value`[1]
  overall.anova.morphology.through.time.P[paste0(t), "Df"] <- get(paste0("anova.P.morph.",t,".lm"))$Df[1]
  
  # Now, do pairwise tests
  assign(paste0("pairwisettest.P.morph.",t) , pairwise.t.test(x = get(paste0("MF.Alpha.P.",t))[,paste0(metric)]
                                                              , g = get(paste0("MF.Alpha.P.",t))[,"Morph"]
                                                              , p.adjust.method = "BH"
  ))
  pairwise.ttest.morphology.through.time.P[paste0(t),"FB:BL"] <- get(paste0("pairwisettest.P.morph.",t))$p.value[2,1]
  pairwise.ttest.morphology.through.time.P[paste0(t),"FB:CR"] <- get(paste0("pairwisettest.P.morph.",t))$p.value[2,2]
  pairwise.ttest.morphology.through.time.P[paste0(t),"BL:CR"] <- get(paste0("pairwisettest.P.morph.",t))$p.value[1,1]
  
}
# Change rownames to something nice
rownames(overall.anova.morphology.through.time.P) <- c("20 min","1 h","3 h","6 h","12 h","1 d")
rownames(pairwise.ttest.morphology.through.time.P) <- c("20 min","1 h","3 h","6 h","12 h","1 d")

# Round the numbers to make sure all are appropriate length
tempfile <- overall.anova.morphology.through.time.P
for (r in 1:nrow(tempfile)) {
  for (c in 1:ncol(tempfile)) {
    if (tempfile[r,c] < 0.001) {
      overall.anova.morphology.through.time.P[r,c] <- "<0.001"
    } else {
      overall.anova.morphology.through.time.P[r,c] <- round(tempfile[r,c], digits = 3)
    }
  }
}

tempfile2 <- pairwise.ttest.morphology.through.time.P
for (r in 1:nrow(tempfile2)) {
  for (c in 1:ncol(tempfile2)) {
    if (tempfile2[r,c] < 0.001) {
      pairwise.ttest.morphology.through.time.P[r,c] <- "<0.001"
    } else {
      pairwise.ttest.morphology.through.time.P[r,c] <- round(tempfile2[r,c], digits = 3)
    }
  }
}

capture.output(xtable(overall.anova.morphology.through.time.P), file = paste0("ALPHAPLOTS_latex/overall_anova_morphology_through_time_P",metric,".txt"))
capture.output(xtable(pairwise.ttest.morphology.through.time.P), file = paste0("ALPHAPLOTS_latex/pairwise_ttest_morphology_through_time_P",metric,".txt"))

########### Plotting Hakai ##############

H.CR.Alpha <- MF.Alpha.morphonly.H[grep("CR", MF.Alpha.morphonly.H$Morph),]
H.BL.Alpha <- MF.Alpha.morphonly.H[grep("BL", MF.Alpha.morphonly.H$Morph),]
H.FB.Alpha <- MF.Alpha.morphonly.H[grep("FB", MF.Alpha.morphonly.H$Morph),]

H.CR.mean <- aggregate(H.CR.Alpha[,paste0(metric)], by=list(Time=H.CR.Alpha$Time), FUN=mean)
H.CR.sd <- aggregate(H.CR.Alpha$chao1, by=list(Time=H.CR.Alpha$Time), FUN=sd)
H.BL.mean <- aggregate(H.BL.Alpha$chao1, by=list(Time=H.BL.Alpha$Time), FUN=mean)
H.BL.sd <- aggregate(H.BL.Alpha$chao1, by=list(Time=H.BL.Alpha$Time), FUN=sd)
H.FB.mean <- aggregate(H.FB.Alpha$chao1, by=list(Time=H.FB.Alpha$Time), FUN=mean)
H.FB.sd <- aggregate(H.FB.Alpha$chao1, by=list(Time=H.FB.Alpha$Time), FUN=sd)

TimeChar <- H.CR.mean$Time[1:4]

pdf(paste0("ALPHAPLOTS_H/Alpha_diversity_Hakai_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0) ~ TimeChar[1:4]
     , main = paste0(metric)
     , ylim = c(min(MF.Alpha.morphonly.H$chao1), max(MF.Alpha.morphonly.H$chao1))
     , xlab = "Time"
     , ylab = "Richness (chao1)"
     , xaxt = 'n'
     )
axis(1
      , at = TimeChar[1:4]
     , labels = c("20m","1h","6h","12h")
     , las = 2)
points(TimeChar[1:4], H.CR.mean$x[1:4]
       , lty = 3
       , type = 'l')
arrows(x0 = TimeChar*0.99
         , x1 = TimeChar*0.99
         , y0 = c(H.CR.mean$x - H.CR.sd$x/2)[1:4]
         , y1 = c(H.CR.mean$x + H.CR.sd$x/2)[1:4]
         , angle = 90
       , code = 3
       , length = 0.02
       )
points(TimeChar[1:4], H.BL.mean$x[1:4]
       , lty = 2
       , type = 'l')
arrows(x0 = TimeChar*1.0
       , x1 = TimeChar*1.0
       , y0 = c(H.BL.mean$x - H.BL.sd$x/2)[1:4]
       , y1 = c(H.BL.mean$x + H.BL.sd$x/2)[1:4]
       , angle = 90
       , code = 3
       , length = 0.02)
points(TimeChar[1:4], H.FB.mean$x[1:4]
       , lty = 1
       , type = 'l')
arrows(x0 = TimeChar*1.01
       , x1 = TimeChar*1.01
       , y0 = c(H.FB.mean$x - H.FB.sd$x/2)[1:4]
       , y1 = c(H.FB.mean$x + H.FB.sd$x/2)[1:4]
       , angle = 90
       , code = 3
       , length = 0.02)
text(TimeChar, 1.05*c(H.FB.mean$x + H.FB.sd$x/2)[1:4]
       , labels = c('','**','***','')
       , cex = 2
       )
par(fig = c(0.63,1,0,1), mar= c(5,0,5,0), new = TRUE)
plot(0,0
     , bty = 'n'
     , xlab = ""
     , ylab = ""
     , xaxt = "n"
     , yaxt = "n"
     , pch = '')
legend("top"
       , legend = c("Finely Branched","Bladed","Crustose")
       , lty = c(1,2,3)
       , bty = 'n')
text(-1,0
     , labels = c("ANOVA p-Values:")
     , pos= 4)
text(-1,-0.15
     , labels = c(paste0("20 min : "))
     , pos= 4)
text(-1,-0.3
     , labels = c(paste0("1 h: "))
      , pos = 4)
text(-1,-0.45
     , labels = c(paste0("6 h: "))
     , pos = 4)
text(-1,-0.6
     , labels = c(paste0("12 h: "))
     , pos = 4)
text(-1,-0.8
     , labels = c(paste0("TimexMorph: "))
     , pos = 4)

text(0.15,-0.15
     , labels = c(paste0(format(anova.H.morph.20.lm$`Pr(>F)`[1],digits = 2), ''))
     , pos= 4)
text(0.15,-0.3
     , labels = c(paste0(format(anova.H.morph.60.lm$`Pr(>F)`[1], digits = 2), '**'))
     , pos = 4)
text(0.15,-0.45
     , labels = c(paste0(format(anova.H.morph.360.lm$`Pr(>F)`[1], digits = 2), '***'))
     , pos = 4)
text(0.15,-0.6
     , labels = c(paste0(format(anova.H.morph.720.lm$`Pr(>F)`[1], digits = 2), ''))
     , pos = 4)
text(0.15,-0.8
     , labels = c(paste0(format(anova.H.morph.lm$`Pr(>F)`[4], digits = 2)))
     , pos = 4)
dev.off()

########### Plotting PM ##############

P.CR.Alpha <- MF.Alpha.morphonly.P[grep("CR", MF.Alpha.morphonly.P$Morph),]
P.BL.Alpha <- MF.Alpha.morphonly.P[grep("BL", MF.Alpha.morphonly.P$Morph),]
P.FB.Alpha <- MF.Alpha.morphonly.P[grep("FB", MF.Alpha.morphonly.P$Morph),]

P.CR.mean <- aggregate(P.CR.Alpha[,paste0(metric)], by=list(Time=P.CR.Alpha$Time), FUN=mean)
P.CR.sd <- aggregate(P.CR.Alpha$chao1, by=list(Time=P.CR.Alpha$Time), FUN=sd)
P.BL.mean <- aggregate(P.BL.Alpha$chao1, by=list(Time=P.BL.Alpha$Time), FUN=mean)
P.BL.sd <- aggregate(P.BL.Alpha$chao1, by=list(Time=P.BL.Alpha$Time), FUN=sd)
P.FB.mean <- aggregate(P.FB.Alpha$chao1, by=list(Time=P.FB.Alpha$Time), FUN=mean)
P.FB.sd <- aggregate(P.FB.Alpha$chao1, by=list(Time=P.FB.Alpha$Time), FUN=sd)

TimeChar <- P.CR.mean$Time[1:6]

pdf(paste0("ALPHAPLOTS_P/Alpha_diversity_PM_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0,0,0) ~ TimeChar[1:6]
     , main = paste0(metric)
     , ylim = c(min(MF.Alpha.morphonly.P$chao1), max(MF.Alpha.morphonly.P$chao1))
     , xlab = "Time"
     , ylab = "Richness (chao1)"
     , xaxt = 'n'
)
axis(1
     , at = c(-20,10,0,0,0,0)+TimeChar[1:6]
     , labels = c("20m","1h","3h","6h","12h","24h")
     , las = 2
     , cex.axis = 1)
points(TimeChar[1:6], P.CR.mean$x[1:6]
       , lty = 3
       , type = 'l')
arrows(x0 = TimeChar*0.99
       , x1 = TimeChar*0.99
       , y0 = c(P.CR.mean$x - P.CR.sd$x/2)[1:6]
       , y1 = c(P.CR.mean$x + P.CR.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02
)
points(TimeChar[1:6], P.BL.mean$x[1:6]
       , lty = 2
       , type = 'l')
arrows(x0 = TimeChar*1.0
       , x1 = TimeChar*1.0
       , y0 = c(P.BL.mean$x - P.BL.sd$x/2)[1:6]
       , y1 = c(P.BL.mean$x + P.BL.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02)
points(TimeChar[1:6], P.FB.mean$x[1:6]
       , lty = 1
       , type = 'l')
arrows(x0 = TimeChar*1.01
       , x1 = TimeChar*1.01
       , y0 = c(P.FB.mean$x - P.FB.sd$x/2)[1:6]
       , y1 = c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02)
text(TimeChar, 1.05*c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
       , labels = c('*','**','*','','*','')
       , cex = 2)
par(fig = c(0.63,1,0,1), mar= c(5,0,5,0), new = TRUE)
plot(0,0
     , bty = 'n'
     , xlab = ""
     , ylab = ""
     , xaxt = "n"
     , yaxt = "n"
     , pch = '')
legend("top"
       , legend = c("Finely Branched","Bladed","Crustose")
       , lty = c(1,2,3)
       , bty = 'n')
text(-1,0.15
     , labels = c("ANOVA p-Values:")
     , pos= 4)
text(-1,0
     , labels = c(paste0("20 min: "))
     , pos= 4)
text(-1,-0.15
     , labels = c(paste0("1 h: "))
     , pos = 4)
text(-1,-0.3
     , labels = c(paste0("3 h: "))
     , pos = 4)
text(-1,-0.45
     , labels = c(paste0("6 h: "))
     , pos = 4)
text(-1,-0.6
     , labels = c(paste0("12 h: "))
     , pos = 4)
text(-1,-0.75
     , labels = c(paste0("24 h: "))
     , pos = 4)
text(-1,-0.95
     , labels = c(paste0("TimexMorph: "))
     , pos = 4)

text(0.15,0
     , labels = c(paste0(format(anova.P.morph.20.lm$`Pr(>F)`[1],digits = 2), "*"))
     , pos= 4)
text(0.15,-0.15
     , labels = c(paste0(format(anova.P.morph.60.lm$`Pr(>F)`[1], digits = 2), "**"))
     , pos = 4)
text(0.15,-0.3
     , labels = c(paste0(format(anova.P.morph.180.lm$`Pr(>F)`[1], digits = 2), "*"))
     , pos = 4)
text(0.15,-0.45
     , labels = c(paste0(format(anova.P.morph.360.lm$`Pr(>F)`[1], digits = 2), ""))
     , pos = 4)
text(0.15,-0.6
     , labels = c(paste0(format(anova.P.morph.720.lm$`Pr(>F)`[1], digits = 2), "*"))
     , pos = 4)
text(0.15,-0.75
     , labels = c(paste0(format(anova.P.morph.1440.lm$`Pr(>F)`[1], digits = 2), ""))
     , pos = 4)
text(0.15,-0.95
     , labels = c(paste0(format(anova.P.morph.lm$`Pr(>F)`[4], digits = 2)))
     , pos = 4)

dev.off()


########### ****PD_WHOLE_TREE**** ##############
metric <- 'PD_whole_tree'


########### Stats ##############

H.morph.lm <- lm(chao1 ~ Time*Morph, data = MF.Alpha.morphonly.H)
anova.H.morph.lm <- Anova(H.morph.lm, type = 3)
capture.output(xtable(anova.H.morph.lm), file = paste0("ALPHAPLOTS_latex/anova_Hakai_", metric, ".txt"))
P.morph.lm <- lm(chao1 ~ Time*Morph, data = MF.Alpha.morphonly.P)
anova.P.morph.lm <- Anova(P.morph.lm, type = 3)
capture.output(xtable(anova.P.morph.lm), file = paste0("ALPHAPLOTS_latex/anova_PM_", metric, ".txt"))

# Separated out Stats- Hakai
overall.anova.morphology.through.time.H <- matrix(nrow = 5, ncol = 3)
colnames(overall.anova.morphology.through.time.H) <- c("p-value","F-value","Df")
rownames(overall.anova.morphology.through.time.H) <- c("20","60","360","720","5760")
pairwise.ttest.morphology.through.time.H <- matrix(nrow = 5, ncol = 3)
colnames(pairwise.ttest.morphology.through.time.H) <- c("FB:BL","FB:CR","BL:CR")
rownames(pairwise.ttest.morphology.through.time.H) <- c("20","60","360","720","5760")

for (t in c("20","60","360","720","5760")) {
  
  # Get overall anova at each time point and make into a table; print this table as latex 
  assign(paste0("MF.Alpha.H.",t), MF.Alpha.morphonly.H[grep(paste0("^",t,"$"), MF.Alpha.morphonly.H$Time),])
  assign(paste0("H.morph.",t,".lm"), lm(chao1 ~ Morph, data = get(paste0("MF.Alpha.H.",t))))
  assign(paste0("anova.H.morph.",t,".lm"), anova(get(paste0("H.morph.",t,".lm"))))
  overall.anova.morphology.through.time.H[paste0(t),"p-value"] <- get(paste0("anova.H.morph.",t,".lm"))$`Pr(>F)`[1]
  overall.anova.morphology.through.time.H[paste0(t), "F-value"] <- get(paste0("anova.H.morph.",t,".lm"))$`F value`[1]
  overall.anova.morphology.through.time.H[paste0(t), "Df"] <- get(paste0("anova.H.morph.",t,".lm"))$Df[1]
  
  # Now, do pairwise tests
  assign(paste0("pairwisettest.H.morph.",t) , pairwise.t.test(x = get(paste0("MF.Alpha.H.",t))[,paste0(metric)]
                                                              , g = get(paste0("MF.Alpha.H.",t))[,"Morph"]
                                                              , p.adjust.method = "BH"
  ))
  pairwise.ttest.morphology.through.time.H[paste0(t),"FB:BL"] <- get(paste0("pairwisettest.H.morph.",t))$p.value[2,1]
  pairwise.ttest.morphology.through.time.H[paste0(t),"FB:CR"] <- get(paste0("pairwisettest.H.morph.",t))$p.value[2,2]
  pairwise.ttest.morphology.through.time.H[paste0(t),"BL:CR"] <- get(paste0("pairwisettest.H.morph.",t))$p.value[1,1]
  
}
# Change rownames to something nice
rownames(overall.anova.morphology.through.time.H) <- c("20 min","1 h","6 h","12 h","4 d")
rownames(pairwise.ttest.morphology.through.time.H) <- c("20 min","1 h","6 h","12 h","4 d")

# Round the numbers to make sure all are appropriate length
tempfile <- overall.anova.morphology.through.time.H
for (r in 1:nrow(tempfile)) {
  for (c in 1:ncol(tempfile)) {
    if (tempfile[r,c] < 0.001) {
      overall.anova.morphology.through.time.H[r,c] <- "<0.001"
    } else {
      overall.anova.morphology.through.time.H[r,c] <- round(tempfile[r,c], digits = 3)
    }
  }
}

tempfile2 <- pairwise.ttest.morphology.through.time.H
for (r in 1:nrow(tempfile2)) {
  for (c in 1:ncol(tempfile2)) {
    if (tempfile2[r,c] < 0.001) {
      pairwise.ttest.morphology.through.time.H[r,c] <- "<0.001"
    } else {
      pairwise.ttest.morphology.through.time.H[r,c] <- round(tempfile2[r,c], digits = 3)
    }
  }
}

capture.output(xtable(overall.anova.morphology.through.time.H), file = paste0("ALPHAPLOTS_latex/overall_anova_morphology_through_time_H",metric,".txt"))
capture.output(xtable(pairwise.ttest.morphology.through.time.H), file = paste0("ALPHAPLOTS_latex/pairwise_ttest_morphology_through_time_H",metric,".txt"))




# Separated out Stats- PM
overall.anova.morphology.through.time.P <- matrix(nrow = 6, ncol = 3)
colnames(overall.anova.morphology.through.time.P) <- c("p-value","F-value","Df")
rownames(overall.anova.morphology.through.time.P) <- c("20","60","180","360","720","1440")
pairwise.ttest.morphology.through.time.P <- matrix(nrow = 6, ncol = 3)
colnames(pairwise.ttest.morphology.through.time.P) <- c("FB:BL","FB:CR","BL:CR")
rownames(pairwise.ttest.morphology.through.time.P) <- c("20","60","180","360","720","1440")

for (t in c("20","60","180","360","720","1440")) {
  
  # Get overall anova at each time point and make into a table; print this table as latex 
  assign(paste0("MF.Alpha.P.",t), MF.Alpha.morphonly.P[grep(paste0("^",t,"$"), MF.Alpha.morphonly.P$Time),])
  assign(paste0("P.morph.",t,".lm"), lm(chao1 ~ Morph, data = get(paste0("MF.Alpha.P.",t))))
  assign(paste0("anova.P.morph.",t,".lm"), anova(get(paste0("P.morph.",t,".lm"))))
  overall.anova.morphology.through.time.P[paste0(t),"p-value"] <- get(paste0("anova.P.morph.",t,".lm"))$`Pr(>F)`[1]
  overall.anova.morphology.through.time.P[paste0(t), "F-value"] <- get(paste0("anova.P.morph.",t,".lm"))$`F value`[1]
  overall.anova.morphology.through.time.P[paste0(t), "Df"] <- get(paste0("anova.P.morph.",t,".lm"))$Df[1]
  
  # Now, do pairwise tests
  assign(paste0("pairwisettest.P.morph.",t) , pairwise.t.test(x = get(paste0("MF.Alpha.P.",t))[,paste0(metric)]
                                                              , g = get(paste0("MF.Alpha.P.",t))[,"Morph"]
                                                              , p.adjust.method = "BH"
  ))
  pairwise.ttest.morphology.through.time.P[paste0(t),"FB:BL"] <- get(paste0("pairwisettest.P.morph.",t))$p.value[2,1]
  pairwise.ttest.morphology.through.time.P[paste0(t),"FB:CR"] <- get(paste0("pairwisettest.P.morph.",t))$p.value[2,2]
  pairwise.ttest.morphology.through.time.P[paste0(t),"BL:CR"] <- get(paste0("pairwisettest.P.morph.",t))$p.value[1,1]
  
}
# Change rownames to something nice
rownames(overall.anova.morphology.through.time.P) <- c("20 min","1 h","3 h","6 h","12 h","1 d")
rownames(pairwise.ttest.morphology.through.time.P) <- c("20 min","1 h","3 h","6 h","12 h","1 d")

# Round the numbers to make sure all are appropriate length
tempfile <- overall.anova.morphology.through.time.P
for (r in 1:nrow(tempfile)) {
  for (c in 1:ncol(tempfile)) {
    if (tempfile[r,c] < 0.001) {
      overall.anova.morphology.through.time.P[r,c] <- "<0.001"
    } else {
      overall.anova.morphology.through.time.P[r,c] <- round(tempfile[r,c], digits = 3)
    }
  }
}

tempfile2 <- pairwise.ttest.morphology.through.time.P
for (r in 1:nrow(tempfile2)) {
  for (c in 1:ncol(tempfile2)) {
    if (tempfile2[r,c] < 0.001) {
      pairwise.ttest.morphology.through.time.P[r,c] <- "<0.001"
    } else {
      pairwise.ttest.morphology.through.time.P[r,c] <- round(tempfile2[r,c], digits = 3)
    }
  }
}

capture.output(xtable(overall.anova.morphology.through.time.P), file = paste0("ALPHAPLOTS_latex/overall_anova_morphology_through_time_P",metric,".txt"))
capture.output(xtable(pairwise.ttest.morphology.through.time.P), file = paste0("ALPHAPLOTS_latex/pairwise_ttest_morphology_through_time_P",metric,".txt"))

########### Plotting Hakai ##############

H.CR.Alpha <- MF.Alpha.morphonly.H[grep("CR", MF.Alpha.morphonly.H$Morph),]
H.BL.Alpha <- MF.Alpha.morphonly.H[grep("BL", MF.Alpha.morphonly.H$Morph),]
H.FB.Alpha <- MF.Alpha.morphonly.H[grep("FB", MF.Alpha.morphonly.H$Morph),]

H.CR.mean <- aggregate(H.CR.Alpha[,paste0(metric)], by=list(Time=H.CR.Alpha$Time), FUN=mean)
H.CR.sd <- aggregate(H.CR.Alpha$chao1, by=list(Time=H.CR.Alpha$Time), FUN=sd)
H.BL.mean <- aggregate(H.BL.Alpha$chao1, by=list(Time=H.BL.Alpha$Time), FUN=mean)
H.BL.sd <- aggregate(H.BL.Alpha$chao1, by=list(Time=H.BL.Alpha$Time), FUN=sd)
H.FB.mean <- aggregate(H.FB.Alpha$chao1, by=list(Time=H.FB.Alpha$Time), FUN=mean)
H.FB.sd <- aggregate(H.FB.Alpha$chao1, by=list(Time=H.FB.Alpha$Time), FUN=sd)

TimeChar <- H.CR.mean$Time[1:4]

pdf(paste0("ALPHAPLOTS_H/Alpha_diversity_Hakai_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0) ~ TimeChar[1:4]
     , main = paste0(metric)
     , ylim = c(min(MF.Alpha.morphonly.H$chao1), max(MF.Alpha.morphonly.H$chao1))
     , xlab = "Time"
     , ylab = "Richness (PD_whole_tree)"
     , xaxt = 'n'
)
axis(1
     , at = TimeChar[1:4]
     , labels = c("20m","1h","6h","12h")
     , las = 2)
points(TimeChar[1:4], H.CR.mean$x[1:4]
       , lty = 3
       , type = 'l')
arrows(x0 = TimeChar*0.99
       , x1 = TimeChar*0.99
       , y0 = c(H.CR.mean$x - H.CR.sd$x/2)[1:4]
       , y1 = c(H.CR.mean$x + H.CR.sd$x/2)[1:4]
       , angle = 90
       , code = 3
       , length = 0.02
)
points(TimeChar[1:4], H.BL.mean$x[1:4]
       , lty = 2
       , type = 'l')
arrows(x0 = TimeChar*1.0
       , x1 = TimeChar*1.0
       , y0 = c(H.BL.mean$x - H.BL.sd$x/2)[1:4]
       , y1 = c(H.BL.mean$x + H.BL.sd$x/2)[1:4]
       , angle = 90
       , code = 3
       , length = 0.02)
points(TimeChar[1:4], H.FB.mean$x[1:4]
       , lty = 1
       , type = 'l')
arrows(x0 = TimeChar*1.01
       , x1 = TimeChar*1.01
       , y0 = c(H.FB.mean$x - H.FB.sd$x/2)[1:4]
       , y1 = c(H.FB.mean$x + H.FB.sd$x/2)[1:4]
       , angle = 90
       , code = 3
       , length = 0.02)
text(TimeChar, 1.05*c(H.FB.mean$x + H.FB.sd$x/2)[1:4]
     , labels = c('','**','***','')
     , cex = 2
)
par(fig = c(0.63,1,0,1), mar= c(5,0,5,0), new = TRUE)
plot(0,0
     , bty = 'n'
     , xlab = ""
     , ylab = ""
     , xaxt = "n"
     , yaxt = "n"
     , pch = '')
legend("top"
       , legend = c("Finely Branched","Bladed","Crustose")
       , lty = c(1,2,3)
       , bty = 'n')
text(-1,0
     , labels = c("ANOVA p-Values:")
     , pos= 4)
text(-1,-0.15
     , labels = c(paste0("20 min : "))
     , pos= 4)
text(-1,-0.3
     , labels = c(paste0("1 h: "))
     , pos = 4)
text(-1,-0.45
     , labels = c(paste0("6 h: "))
     , pos = 4)
text(-1,-0.6
     , labels = c(paste0("12 h: "))
     , pos = 4)
text(-1,-0.8
     , labels = c(paste0("TimexMorph: "))
     , pos = 4)

text(0.15,-0.15
     , labels = c(paste0(format(anova.H.morph.20.lm$`Pr(>F)`[1],digits = 2), ''))
     , pos= 4)
text(0.15,-0.3
     , labels = c(paste0(format(anova.H.morph.60.lm$`Pr(>F)`[1], digits = 2), '**'))
     , pos = 4)
text(0.15,-0.45
     , labels = c(paste0(format(anova.H.morph.360.lm$`Pr(>F)`[1], digits = 2), '***'))
     , pos = 4)
text(0.15,-0.6
     , labels = c(paste0(format(anova.H.morph.720.lm$`Pr(>F)`[1], digits = 2), ''))
     , pos = 4)
text(0.15,-0.8
     , labels = c(paste0(format(anova.H.morph.lm$`Pr(>F)`[4], digits = 2)))
     , pos = 4)
dev.off()

########### Plotting PM ##############

P.CR.Alpha <- MF.Alpha.morphonly.P[grep("CR", MF.Alpha.morphonly.P$Morph),]
P.BL.Alpha <- MF.Alpha.morphonly.P[grep("BL", MF.Alpha.morphonly.P$Morph),]
P.FB.Alpha <- MF.Alpha.morphonly.P[grep("FB", MF.Alpha.morphonly.P$Morph),]

P.CR.mean <- aggregate(P.CR.Alpha[,paste0(metric)], by=list(Time=P.CR.Alpha$Time), FUN=mean)
P.CR.sd <- aggregate(P.CR.Alpha$chao1, by=list(Time=P.CR.Alpha$Time), FUN=sd)
P.BL.mean <- aggregate(P.BL.Alpha$chao1, by=list(Time=P.BL.Alpha$Time), FUN=mean)
P.BL.sd <- aggregate(P.BL.Alpha$chao1, by=list(Time=P.BL.Alpha$Time), FUN=sd)
P.FB.mean <- aggregate(P.FB.Alpha$chao1, by=list(Time=P.FB.Alpha$Time), FUN=mean)
P.FB.sd <- aggregate(P.FB.Alpha$chao1, by=list(Time=P.FB.Alpha$Time), FUN=sd)

TimeChar <- P.CR.mean$Time[1:6]

pdf(paste0("ALPHAPLOTS_P/Alpha_diversity_PM_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0,0,0) ~ TimeChar[1:6]
     , main = paste0(metric)
     , ylim = c(min(MF.Alpha.morphonly.P$chao1), max(MF.Alpha.morphonly.P$chao1))
     , xlab = "Time"
     , ylab = "Richness (PD_whole_tree)"
     , xaxt = 'n'
)
axis(1
     , at = c(-20,10,0,0,0,0)+TimeChar[1:6]
     , labels = c("20m","1h","3h","6h","12h","24h")
     , las = 2
     , cex.axis = 1)
points(TimeChar[1:6], P.CR.mean$x[1:6]
       , lty = 3
       , type = 'l')
arrows(x0 = TimeChar*0.99
       , x1 = TimeChar*0.99
       , y0 = c(P.CR.mean$x - P.CR.sd$x/2)[1:6]
       , y1 = c(P.CR.mean$x + P.CR.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02
)
points(TimeChar[1:6], P.BL.mean$x[1:6]
       , lty = 2
       , type = 'l')
arrows(x0 = TimeChar*1.0
       , x1 = TimeChar*1.0
       , y0 = c(P.BL.mean$x - P.BL.sd$x/2)[1:6]
       , y1 = c(P.BL.mean$x + P.BL.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02)
points(TimeChar[1:6], P.FB.mean$x[1:6]
       , lty = 1
       , type = 'l')
arrows(x0 = TimeChar*1.01
       , x1 = TimeChar*1.01
       , y0 = c(P.FB.mean$x - P.FB.sd$x/2)[1:6]
       , y1 = c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02)
text(TimeChar, 1.05*c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
     , labels = c('*','**','*','','*','')
     , cex = 2)
par(fig = c(0.63,1,0,1), mar= c(5,0,5,0), new = TRUE)
plot(0,0
     , bty = 'n'
     , xlab = ""
     , ylab = ""
     , xaxt = "n"
     , yaxt = "n"
     , pch = '')
legend("top"
       , legend = c("Finely Branched","Bladed","Crustose")
       , lty = c(1,2,3)
       , bty = 'n')
text(-1,0.15
     , labels = c("ANOVA p-Values:")
     , pos= 4)
text(-1,0
     , labels = c(paste0("20 min: "))
     , pos= 4)
text(-1,-0.15
     , labels = c(paste0("1 h: "))
     , pos = 4)
text(-1,-0.3
     , labels = c(paste0("3 h: "))
     , pos = 4)
text(-1,-0.45
     , labels = c(paste0("6 h: "))
     , pos = 4)
text(-1,-0.6
     , labels = c(paste0("12 h: "))
     , pos = 4)
text(-1,-0.75
     , labels = c(paste0("24 h: "))
     , pos = 4)
text(-1,-0.95
     , labels = c(paste0("TimexMorph: "))
     , pos = 4)

text(0.15,0
     , labels = c(paste0(format(anova.P.morph.20.lm$`Pr(>F)`[1],digits = 2), "*"))
     , pos= 4)
text(0.15,-0.15
     , labels = c(paste0(format(anova.P.morph.60.lm$`Pr(>F)`[1], digits = 2), "**"))
     , pos = 4)
text(0.15,-0.3
     , labels = c(paste0(format(anova.P.morph.180.lm$`Pr(>F)`[1], digits = 2), "*"))
     , pos = 4)
text(0.15,-0.45
     , labels = c(paste0(format(anova.P.morph.360.lm$`Pr(>F)`[1], digits = 2), ""))
     , pos = 4)
text(0.15,-0.6
     , labels = c(paste0(format(anova.P.morph.720.lm$`Pr(>F)`[1], digits = 2), "*"))
     , pos = 4)
text(0.15,-0.75
     , labels = c(paste0(format(anova.P.morph.1440.lm$`Pr(>F)`[1], digits = 2), ""))
     , pos = 4)
text(0.15,-0.95
     , labels = c(paste0(format(anova.P.morph.lm$`Pr(>F)`[4], digits = 2)))
     , pos = 4)

dev.off()



########### ****OBSERVED OTUs**** ##############
metric <- 'observed_otus'

########### Stats ##############

H.morph.lm <- lm(chao1 ~ Time*Morph, data = MF.Alpha.morphonly.H)
anova.H.morph.lm <- Anova(H.morph.lm, type = 3)
capture.output(xtable(anova.H.morph.lm), file = paste0("ALPHAPLOTS_latex/anova_Hakai_", metric, ".txt"))
P.morph.lm <- lm(chao1 ~ Time*Morph, data = MF.Alpha.morphonly.P)
anova.P.morph.lm <- Anova(P.morph.lm, type = 3)
capture.output(xtable(anova.P.morph.lm), file = paste0("ALPHAPLOTS_latex/anova_PM_", metric, ".txt"))

# Separated out Stats- Hakai
overall.anova.morphology.through.time.H <- matrix(nrow = 5, ncol = 3)
colnames(overall.anova.morphology.through.time.H) <- c("p-value","F-value","Df")
rownames(overall.anova.morphology.through.time.H) <- c("20","60","360","720","5760")
pairwise.ttest.morphology.through.time.H <- matrix(nrow = 5, ncol = 3)
colnames(pairwise.ttest.morphology.through.time.H) <- c("FB:BL","FB:CR","BL:CR")
rownames(pairwise.ttest.morphology.through.time.H) <- c("20","60","360","720","5760")

for (t in c("20","60","360","720","5760")) {
  
  # Get overall anova at each time point and make into a table; print this table as latex 
  assign(paste0("MF.Alpha.H.",t), MF.Alpha.morphonly.H[grep(paste0("^",t,"$"), MF.Alpha.morphonly.H$Time),])
  assign(paste0("H.morph.",t,".lm"), lm(chao1 ~ Morph, data = get(paste0("MF.Alpha.H.",t))))
  assign(paste0("anova.H.morph.",t,".lm"), anova(get(paste0("H.morph.",t,".lm"))))
  overall.anova.morphology.through.time.H[paste0(t),"p-value"] <- get(paste0("anova.H.morph.",t,".lm"))$`Pr(>F)`[1]
  overall.anova.morphology.through.time.H[paste0(t), "F-value"] <- get(paste0("anova.H.morph.",t,".lm"))$`F value`[1]
  overall.anova.morphology.through.time.H[paste0(t), "Df"] <- get(paste0("anova.H.morph.",t,".lm"))$Df[1]
  
  # Now, do pairwise tests
  assign(paste0("pairwisettest.H.morph.",t) , pairwise.t.test(x = get(paste0("MF.Alpha.H.",t))[,paste0(metric)]
                                                              , g = get(paste0("MF.Alpha.H.",t))[,"Morph"]
                                                              , p.adjust.method = "BH"
  ))
  pairwise.ttest.morphology.through.time.H[paste0(t),"FB:BL"] <- get(paste0("pairwisettest.H.morph.",t))$p.value[2,1]
  pairwise.ttest.morphology.through.time.H[paste0(t),"FB:CR"] <- get(paste0("pairwisettest.H.morph.",t))$p.value[2,2]
  pairwise.ttest.morphology.through.time.H[paste0(t),"BL:CR"] <- get(paste0("pairwisettest.H.morph.",t))$p.value[1,1]
  
}
# Change rownames to something nice
rownames(overall.anova.morphology.through.time.H) <- c("20 min","1 h","6 h","12 h","4 d")
rownames(pairwise.ttest.morphology.through.time.H) <- c("20 min","1 h","6 h","12 h","4 d")

# Round the numbers to make sure all are appropriate length
tempfile <- overall.anova.morphology.through.time.H
for (r in 1:nrow(tempfile)) {
  for (c in 1:ncol(tempfile)) {
    if (tempfile[r,c] < 0.001) {
      overall.anova.morphology.through.time.H[r,c] <- "<0.001"
    } else {
      overall.anova.morphology.through.time.H[r,c] <- round(tempfile[r,c], digits = 3)
    }
  }
}

tempfile2 <- pairwise.ttest.morphology.through.time.H
for (r in 1:nrow(tempfile2)) {
  for (c in 1:ncol(tempfile2)) {
    if (tempfile2[r,c] < 0.001) {
      pairwise.ttest.morphology.through.time.H[r,c] <- "<0.001"
    } else {
      pairwise.ttest.morphology.through.time.H[r,c] <- round(tempfile2[r,c], digits = 3)
    }
  }
}

capture.output(xtable(overall.anova.morphology.through.time.H), file = paste0("ALPHAPLOTS_latex/overall_anova_morphology_through_time_H",metric,".txt"))
capture.output(xtable(pairwise.ttest.morphology.through.time.H), file = paste0("ALPHAPLOTS_latex/pairwise_ttest_morphology_through_time_H",metric,".txt"))




# Separated out Stats- PM
overall.anova.morphology.through.time.P <- matrix(nrow = 6, ncol = 3)
colnames(overall.anova.morphology.through.time.P) <- c("p-value","F-value","Df")
rownames(overall.anova.morphology.through.time.P) <- c("20","60","180","360","720","1440")
pairwise.ttest.morphology.through.time.P <- matrix(nrow = 6, ncol = 3)
colnames(pairwise.ttest.morphology.through.time.P) <- c("FB:BL","FB:CR","BL:CR")
rownames(pairwise.ttest.morphology.through.time.P) <- c("20","60","180","360","720","1440")

for (t in c("20","60","180","360","720","1440")) {
  
  # Get overall anova at each time point and make into a table; print this table as latex 
  assign(paste0("MF.Alpha.P.",t), MF.Alpha.morphonly.P[grep(paste0("^",t,"$"), MF.Alpha.morphonly.P$Time),])
  assign(paste0("P.morph.",t,".lm"), lm(chao1 ~ Morph, data = get(paste0("MF.Alpha.P.",t))))
  assign(paste0("anova.P.morph.",t,".lm"), anova(get(paste0("P.morph.",t,".lm"))))
  overall.anova.morphology.through.time.P[paste0(t),"p-value"] <- get(paste0("anova.P.morph.",t,".lm"))$`Pr(>F)`[1]
  overall.anova.morphology.through.time.P[paste0(t), "F-value"] <- get(paste0("anova.P.morph.",t,".lm"))$`F value`[1]
  overall.anova.morphology.through.time.P[paste0(t), "Df"] <- get(paste0("anova.P.morph.",t,".lm"))$Df[1]
  
  # Now, do pairwise tests
  assign(paste0("pairwisettest.P.morph.",t) , pairwise.t.test(x = get(paste0("MF.Alpha.P.",t))[,paste0(metric)]
                                                              , g = get(paste0("MF.Alpha.P.",t))[,"Morph"]
                                                              , p.adjust.method = "BH"
  ))
  pairwise.ttest.morphology.through.time.P[paste0(t),"FB:BL"] <- get(paste0("pairwisettest.P.morph.",t))$p.value[2,1]
  pairwise.ttest.morphology.through.time.P[paste0(t),"FB:CR"] <- get(paste0("pairwisettest.P.morph.",t))$p.value[2,2]
  pairwise.ttest.morphology.through.time.P[paste0(t),"BL:CR"] <- get(paste0("pairwisettest.P.morph.",t))$p.value[1,1]
  
}
# Change rownames to something nice
rownames(overall.anova.morphology.through.time.P) <- c("20 min","1 h","3 h","6 h","12 h","1 d")
rownames(pairwise.ttest.morphology.through.time.P) <- c("20 min","1 h","3 h","6 h","12 h","1 d")

# Round the numbers to make sure all are appropriate length
tempfile <- overall.anova.morphology.through.time.P
for (r in 1:nrow(tempfile)) {
  for (c in 1:ncol(tempfile)) {
    if (tempfile[r,c] < 0.001) {
      overall.anova.morphology.through.time.P[r,c] <- "<0.001"
    } else {
      overall.anova.morphology.through.time.P[r,c] <- round(tempfile[r,c], digits = 3)
    }
  }
}

tempfile2 <- pairwise.ttest.morphology.through.time.P
for (r in 1:nrow(tempfile2)) {
  for (c in 1:ncol(tempfile2)) {
    if (tempfile2[r,c] < 0.001) {
      pairwise.ttest.morphology.through.time.P[r,c] <- "<0.001"
    } else {
      pairwise.ttest.morphology.through.time.P[r,c] <- round(tempfile2[r,c], digits = 3)
    }
  }
}

capture.output(xtable(overall.anova.morphology.through.time.P), file = paste0("ALPHAPLOTS_latex/overall_anova_morphology_through_time_P",metric,".txt"))
capture.output(xtable(pairwise.ttest.morphology.through.time.P), file = paste0("ALPHAPLOTS_latex/pairwise_ttest_morphology_through_time_P",metric,".txt"))

########### Plotting Hakai ##############

H.CR.Alpha <- MF.Alpha.morphonly.H[grep("CR", MF.Alpha.morphonly.H$Morph),]
H.BL.Alpha <- MF.Alpha.morphonly.H[grep("BL", MF.Alpha.morphonly.H$Morph),]
H.FB.Alpha <- MF.Alpha.morphonly.H[grep("FB", MF.Alpha.morphonly.H$Morph),]

H.CR.mean <- aggregate(H.CR.Alpha[,paste0(metric)], by=list(Time=H.CR.Alpha$Time), FUN=mean)
H.CR.sd <- aggregate(H.CR.Alpha$chao1, by=list(Time=H.CR.Alpha$Time), FUN=sd)
H.BL.mean <- aggregate(H.BL.Alpha$chao1, by=list(Time=H.BL.Alpha$Time), FUN=mean)
H.BL.sd <- aggregate(H.BL.Alpha$chao1, by=list(Time=H.BL.Alpha$Time), FUN=sd)
H.FB.mean <- aggregate(H.FB.Alpha$chao1, by=list(Time=H.FB.Alpha$Time), FUN=mean)
H.FB.sd <- aggregate(H.FB.Alpha$chao1, by=list(Time=H.FB.Alpha$Time), FUN=sd)

TimeChar <- H.CR.mean$Time[1:4]

pdf(paste0("ALPHAPLOTS_H/Alpha_diversity_Hakai_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0) ~ TimeChar[1:4]
     , main = paste0(metric)
     , ylim = c(min(MF.Alpha.morphonly.H$chao1), max(MF.Alpha.morphonly.H$chao1))
     , xlab = "Time"
     , ylab = "Richness (observed_otus)"
     , xaxt = 'n'
)
axis(1
     , at = TimeChar[1:4]
     , labels = c("20m","1h","6h","12h")
     , las = 2)
points(TimeChar[1:4], H.CR.mean$x[1:4]
       , lty = 3
       , type = 'l')
arrows(x0 = TimeChar*0.99
       , x1 = TimeChar*0.99
       , y0 = c(H.CR.mean$x - H.CR.sd$x/2)[1:4]
       , y1 = c(H.CR.mean$x + H.CR.sd$x/2)[1:4]
       , angle = 90
       , code = 3
       , length = 0.02
)
points(TimeChar[1:4], H.BL.mean$x[1:4]
       , lty = 2
       , type = 'l')
arrows(x0 = TimeChar*1.0
       , x1 = TimeChar*1.0
       , y0 = c(H.BL.mean$x - H.BL.sd$x/2)[1:4]
       , y1 = c(H.BL.mean$x + H.BL.sd$x/2)[1:4]
       , angle = 90
       , code = 3
       , length = 0.02)
points(TimeChar[1:4], H.FB.mean$x[1:4]
       , lty = 1
       , type = 'l')
arrows(x0 = TimeChar*1.01
       , x1 = TimeChar*1.01
       , y0 = c(H.FB.mean$x - H.FB.sd$x/2)[1:4]
       , y1 = c(H.FB.mean$x + H.FB.sd$x/2)[1:4]
       , angle = 90
       , code = 3
       , length = 0.02)
text(TimeChar, 1.05*c(H.FB.mean$x + H.FB.sd$x/2)[1:4]
     , labels = c('','**','***','')
     , cex = 2
)
par(fig = c(0.63,1,0,1), mar= c(5,0,5,0), new = TRUE)
plot(0,0
     , bty = 'n'
     , xlab = ""
     , ylab = ""
     , xaxt = "n"
     , yaxt = "n"
     , pch = '')
legend("top"
       , legend = c("Finely Branched","Bladed","Crustose")
       , lty = c(1,2,3)
       , bty = 'n')
text(-1,0
     , labels = c("ANOVA p-Values:")
     , pos= 4)
text(-1,-0.15
     , labels = c(paste0("20 min : "))
     , pos= 4)
text(-1,-0.3
     , labels = c(paste0("1 h: "))
     , pos = 4)
text(-1,-0.45
     , labels = c(paste0("6 h: "))
     , pos = 4)
text(-1,-0.6
     , labels = c(paste0("12 h: "))
     , pos = 4)
text(-1,-0.8
     , labels = c(paste0("TimexMorph: "))
     , pos = 4)

text(0.15,-0.15
     , labels = c(paste0(format(anova.H.morph.20.lm$`Pr(>F)`[1],digits = 2), ''))
     , pos= 4)
text(0.15,-0.3
     , labels = c(paste0(format(anova.H.morph.60.lm$`Pr(>F)`[1], digits = 2), '**'))
     , pos = 4)
text(0.15,-0.45
     , labels = c(paste0(format(anova.H.morph.360.lm$`Pr(>F)`[1], digits = 2), '***'))
     , pos = 4)
text(0.15,-0.6
     , labels = c(paste0(format(anova.H.morph.720.lm$`Pr(>F)`[1], digits = 2), ''))
     , pos = 4)
text(0.15,-0.8
     , labels = c(paste0(format(anova.H.morph.lm$`Pr(>F)`[4], digits = 2)))
     , pos = 4)
dev.off()

########### Plotting PM ##############

P.CR.Alpha <- MF.Alpha.morphonly.P[grep("CR", MF.Alpha.morphonly.P$Morph),]
P.BL.Alpha <- MF.Alpha.morphonly.P[grep("BL", MF.Alpha.morphonly.P$Morph),]
P.FB.Alpha <- MF.Alpha.morphonly.P[grep("FB", MF.Alpha.morphonly.P$Morph),]

P.CR.mean <- aggregate(P.CR.Alpha[,paste0(metric)], by=list(Time=P.CR.Alpha$Time), FUN=mean)
P.CR.sd <- aggregate(P.CR.Alpha$chao1, by=list(Time=P.CR.Alpha$Time), FUN=sd)
P.BL.mean <- aggregate(P.BL.Alpha$chao1, by=list(Time=P.BL.Alpha$Time), FUN=mean)
P.BL.sd <- aggregate(P.BL.Alpha$chao1, by=list(Time=P.BL.Alpha$Time), FUN=sd)
P.FB.mean <- aggregate(P.FB.Alpha$chao1, by=list(Time=P.FB.Alpha$Time), FUN=mean)
P.FB.sd <- aggregate(P.FB.Alpha$chao1, by=list(Time=P.FB.Alpha$Time), FUN=sd)

TimeChar <- P.CR.mean$Time[1:6]

pdf(paste0("ALPHAPLOTS_P/Alpha_diversity_PM_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0,0,0) ~ TimeChar[1:6]
     , main = paste0(metric)
     , ylim = c(min(MF.Alpha.morphonly.P$chao1), max(MF.Alpha.morphonly.P$chao1))
     , xlab = "Time"
     , ylab = "Richness (observed_otus)"
     , xaxt = 'n'
)
axis(1
     , at = c(-20,10,0,0,0,0)+TimeChar[1:6]
     , labels = c("20m","1h","3h","6h","12h","24h")
     , las = 2
     , cex.axis = 1)
points(TimeChar[1:6], P.CR.mean$x[1:6]
       , lty = 3
       , type = 'l')
arrows(x0 = TimeChar*0.99
       , x1 = TimeChar*0.99
       , y0 = c(P.CR.mean$x - P.CR.sd$x/2)[1:6]
       , y1 = c(P.CR.mean$x + P.CR.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02
)
points(TimeChar[1:6], P.BL.mean$x[1:6]
       , lty = 2
       , type = 'l')
arrows(x0 = TimeChar*1.0
       , x1 = TimeChar*1.0
       , y0 = c(P.BL.mean$x - P.BL.sd$x/2)[1:6]
       , y1 = c(P.BL.mean$x + P.BL.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02)
points(TimeChar[1:6], P.FB.mean$x[1:6]
       , lty = 1
       , type = 'l')
arrows(x0 = TimeChar*1.01
       , x1 = TimeChar*1.01
       , y0 = c(P.FB.mean$x - P.FB.sd$x/2)[1:6]
       , y1 = c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02)
text(TimeChar, 1.05*c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
     , labels = c('*','**','*','','*','')
     , cex = 2)
par(fig = c(0.63,1,0,1), mar= c(5,0,5,0), new = TRUE)
plot(0,0
     , bty = 'n'
     , xlab = ""
     , ylab = ""
     , xaxt = "n"
     , yaxt = "n"
     , pch = '')
legend("top"
       , legend = c("Finely Branched","Bladed","Crustose")
       , lty = c(1,2,3)
       , bty = 'n')
text(-1,0.15
     , labels = c("ANOVA p-Values:")
     , pos= 4)
text(-1,0
     , labels = c(paste0("20 min: "))
     , pos= 4)
text(-1,-0.15
     , labels = c(paste0("1 h: "))
     , pos = 4)
text(-1,-0.3
     , labels = c(paste0("3 h: "))
     , pos = 4)
text(-1,-0.45
     , labels = c(paste0("6 h: "))
     , pos = 4)
text(-1,-0.6
     , labels = c(paste0("12 h: "))
     , pos = 4)
text(-1,-0.75
     , labels = c(paste0("24 h: "))
     , pos = 4)
text(-1,-0.95
     , labels = c(paste0("TimexMorph: "))
     , pos = 4)

text(0.15,0
     , labels = c(paste0(format(anova.P.morph.20.lm$`Pr(>F)`[1],digits = 2), "*"))
     , pos= 4)
text(0.15,-0.15
     , labels = c(paste0(format(anova.P.morph.60.lm$`Pr(>F)`[1], digits = 2), "**"))
     , pos = 4)
text(0.15,-0.3
     , labels = c(paste0(format(anova.P.morph.180.lm$`Pr(>F)`[1], digits = 2), "*"))
     , pos = 4)
text(0.15,-0.45
     , labels = c(paste0(format(anova.P.morph.360.lm$`Pr(>F)`[1], digits = 2), ""))
     , pos = 4)
text(0.15,-0.6
     , labels = c(paste0(format(anova.P.morph.720.lm$`Pr(>F)`[1], digits = 2), "*"))
     , pos = 4)
text(0.15,-0.75
     , labels = c(paste0(format(anova.P.morph.1440.lm$`Pr(>F)`[1], digits = 2), ""))
     , pos = 4)
text(0.15,-0.95
     , labels = c(paste0(format(anova.P.morph.lm$`Pr(>F)`[4], digits = 2)))
     , pos = 4)

dev.off()


############# ******WATER COMPARISONS****** ##############
allwaters <- MF.Alpha.morphwater[grep("(H_W)|(P_H2)", MF.Alpha.morphwater$TypeMorphTime),]
for (i in alphaList) {
  temp.lm <- lm(as.numeric(allwaters[,i]) ~ allwaters$Type+allwaters$Time)
  anova.allwaters <- Anova(temp.lm, type = "III")
  capture.output(xtable(anova.allwaters, digits = 3), file = paste0("./ALPHA_compare/wateronly_",i,".txt"))

}

############# ******PER MORPH COMPARISONS****** ##############

#### HAKAI vs PM DIVERSITY ####

for (metric in c("chao1","PD_whole_tree","observed_otus")) {
  # Empty matrix for results; name the headers
  results.HPM.adivcompare <- matrix(nrow = 5, ncol = 3)
  colnames(results.HPM.adivcompare) <- c("p-value","t-statistic","Df")
  rownames(results.HPM.adivcompare) <- c("20","60","360","720", "Overall")
  # Iterate through all shared timepoints
  for (i in c("20","60","360","720")) {
    assign(paste0("ttest.HvsP.",i), t.test(get(paste0("MF.Alpha.H.",i))[,metric], get(paste0("MF.Alpha.P.",i))[,metric]))
    results.HPM.adivcompare[paste0(i), "p-value"] <- get(paste0("ttest.HvsP.",i))$p.value
    results.HPM.adivcompare[paste0(i), "t-statistic"] <- get(paste0("ttest.HvsP.",i))$statistic
    results.HPM.adivcompare[paste0(i), "Df"] <- get(paste0("ttest.HvsP.",i))$parameter
  }
  # Add an overall test too
  ttest.HvsP.all <- t.test(MF.Alpha.morphonly.H[,metric], MF.Alpha.morphonly.P[,metric])
  results.HPM.adivcompare["Overall","p-value"] <- ttest.HvsP.all$p.value
  results.HPM.adivcompare["Overall","t-statistic"] <- ttest.HvsP.all$statistic
  results.HPM.adivcompare["Overall","Df"] <- ttest.HvsP.all$parameter
  rownames(results.HPM.adivcompare) <- c("20 min","1 hour","6 hours","12 hours","Overall")
  # change p < 0.001 to character so formatting makes sense
  temp.pvalues <- results.HPM.adivcompare
  for (r in 1:nrow(temp.pvalues)) {
    for (c in 1:ncol(temp.pvalues)) {
      if (temp.pvalues[r,c] < 0.001) {
        results.HPM.adivcompare[r,c] <- "<0.001"
      } else (
        results.HPM.adivcompare[r,c] <- round(temp.pvalues[r,c], digits = 3)
      )
    }
  }
  

  capture.output(xtable(results.HPM.adivcompare), file = paste0("ALPHAPLOTS_latex/LATEX_HPM.adivcompare.",metric,".txt"))
  
  pdf(paste0("./ALPHA_compare/HandPMcompare_",metric,".pdf"),pointsize = 14)
  plot(NULL, NULL
       , xlim = c(1,7)
       , ylim = c(min(as.numeric(MF.Alpha.all[,paste0(metric,"_even_4000_alpha")])), max(as.numeric(MF.Alpha.all[,paste0(metric,"_even_4000_alpha")]))*1.1)
       , ylab = paste0("Diversity (", metric, ")")
       , xaxt = "n"
       , xlab = "Time"
       , main = paste0("Diversity (",metric,") of Hakai vs Reed Point")
  )
  axis(side = 1
       , labels = c("20 m", "1 h","3 h","6 h","12 h","1 d","4 d")
       , at = c(1,2,3,4,5,6,7)
       , las = 2
  )
  points(MF.Alpha.morphonly.H[,metric] ~ factor(MF.Alpha.morphonly.H[,"Time"], levels = c(20,60,180,360,720,1440,5760))
         , col = "red"
  )
  points(MF.Alpha.morphonly.P[,metric] ~ factor(MF.Alpha.morphonly.P[,"Time"], levels = c(20,60,180,360,720,1440,5760))
         , col = "blue")
  legend("topright"
         , legend = c("Hakai","Reed Point")
         , col = c("red","blue")
         , pch = 21)
  axis( side = 1
        , line = 5
        , las = 2
        , at = c(1,2,4,5)
        , labels = )
  
  dev.off()
}




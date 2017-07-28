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
alphaNames = 'chao1_even_4000_alpha,PD_whole_tree_even_4000_alpha,observed_otus_even_4000_alpha'
alphaList <- unlist(strsplit(alphaNames, ","))
metrics <- gsub(pattern = "_even_4000_alpha", "", alphaList)

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
types <- c("H","P")
storeOnlyH <- c()
storeOnlyP <- c()
for (i in 1:nrow(MF.Alpha.morphonly)) {
  if (MF.Alpha.morphonly[i,'Type'] %in% c(paste0(types[1]))) {
    storeOnlyH <- c(storeOnlyH, i)
  } 
  if (MF.Alpha.morphonly[i,'Type'] %in% c(paste0(types[2]))) {
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
############## ***** ALL METRICS ****** ############
for (metric in metrics) {
  ########### Stats ##############
  
  H.morph.lm <- lm(get(paste0(metric)) ~ Time*Morph, data = MF.Alpha.morphonly.H)
  assign(paste0("anova.H.morph.lm.",metric), Anova(H.morph.lm, type = 3))
  capture.output(xtable(get(paste0("anova.H.morph.lm.",metric))), file = paste0("ALPHAPLOTS_latex/anova_Hakai_", metric, ".txt"))
  P.morph.lm <- lm(get(paste0(metric)) ~ Time*Morph, data = MF.Alpha.morphonly.P)
  assign(paste0("anova.P.morph.lm.", metric), Anova(P.morph.lm, type = 3))
  capture.output(xtable(get(paste0("anova.P.morph.lm.", metric))), file = paste0("ALPHAPLOTS_latex/anova_PM_", metric, ".txt"))
  
  assign(paste0("pairwise.ttest.overall.H.",metric), pairwise.t.test(MF.Alpha.morphonly.H[,paste0(metric)], MF.Alpha.morphonly.H[,"Morph"], p.adjust.method = "BH"))
  assign(paste0("pairwise.ttest.overall.P.",metric), pairwise.t.test(MF.Alpha.morphonly.P[,paste0(metric)], MF.Alpha.morphonly.P[,"Morph"], p.adjust.method = "BH"))

  # Separated out Stats- Hakai
  overall.anova.morphology.through.time <- matrix(ncol = 8, nrow = 2)
  rownames(overall.anova.morphology.through.time) <- c("P","H")
  colnames(overall.anova.morphology.through.time) <- c("20","60","180","360","720","1440","5760","Overall")
  pairwise.ttest.morphology.through.time <- matrix(ncol = 8, nrow = 6)
  rownames(pairwise.ttest.morphology.through.time) <- c("pFB:BL","pFB:CR","pBL:CR","hFB:BL","hFB:CR","hBL:CR")
  colnames(pairwise.ttest.morphology.through.time) <- c("20","60","180","360","720","1440","5760","Overall")
  
  #PM
  for (t in c("20","60","180","360","720","1440")) {
    # Get overall anova at each time point and make into a table; print this table as latex 
    assign(paste0("MF.Alpha.P.",t), MF.Alpha.morphonly.P[grep(paste0("^",t,"$"), MF.Alpha.morphonly.P$Time),])
    assign(paste0("P.morph.",t,".lm"), lm(get(metric) ~ Morph, data = get(paste0("MF.Alpha.P.",t))))
    assign(paste0("anova.P.morph.",t,".lm"), anova(get(paste0("P.morph.",t,".lm"))))
    ptemp <- get(paste0("anova.P.morph.",t,".lm"))$`Pr(>F)`[1]
    if (ptemp < 0.001) {
      ptemp <- paste0("\\bf{<0.001}")
    } else if (ptemp <= 0.05) {
      ptemp <- paste0("\\bf{",signif(ptemp,2),"}")
    } else {
      ptemp <- signif(ptemp,2)
    }
    
    ftemp <- get(paste0("anova.P.morph.",t,".lm"))$`F value`[1]
    dftemp <- paste0(get(paste0("anova.P.morph.",t,".lm"))$Df[1],",",get(paste0("anova.P.morph.",t,".lm"))$Df[2])
    toPaste <- paste0(ptemp, " ($F_{",dftemp,"}$","=",round(ftemp,2),")")
    
    overall.anova.morphology.through.time["P",paste0(t)] <- toPaste
    
    # Now, do pairwise tests
    assign(paste0("pairwisettest.P.morph.",t) , pairwise.t.test(x = get(paste0("MF.Alpha.P.",t))[,paste0(metric)]
                                                                , g = get(paste0("MF.Alpha.P.",t))[,"Morph"]
                                                                , p.adjust.method = "BH"
    ))
    fbbltemp <- get(paste0("pairwisettest.P.morph.",t))$p.value[2,1]
    fbcrtemp <- get(paste0("pairwisettest.P.morph.",t))$p.value[2,2]
    blcrtemp <- get(paste0("pairwisettest.P.morph.",t))$p.value[1,1]
    
    ### MAKE BF
    for (i in c("fbbltemp","fbcrtemp","blcrtemp")) {
      if (get(paste0(i)) < 0.001) {
        assign(paste0(i), "\\bf{<0.001}")
      } else if (get(paste0(i)) <= 0.05) {
        assign(paste0(i), paste0("\\bf{",signif(get(paste0(i)),2),"}"))
      } else {
        assign(paste0(i), paste0(signif(get(paste0(i)),2)))
      }
    }
    
    pairwise.ttest.morphology.through.time["pFB:BL", paste0(t)] <- fbbltemp
    pairwise.ttest.morphology.through.time["pFB:CR", paste0(t)] <- fbcrtemp
    pairwise.ttest.morphology.through.time["pBL:CR", paste0(t)] <- blcrtemp
 
  }
  # Hakai
  for (t in c("20","60","360","720","5760")) {
    # Get overall anova at each time point and make into a table; print this table as latex 
    assign(paste0("MF.Alpha.H.",t), MF.Alpha.morphonly.H[grep(paste0("^",t,"$"), MF.Alpha.morphonly.H$Time),])
    assign(paste0("H.morph.",t,".lm"), lm(get(metric) ~ Morph, data = get(paste0("MF.Alpha.H.",t))))
    assign(paste0("anova.H.morph.",t,".lm"), anova(get(paste0("H.morph.",t,".lm"))))
    ptemp <- get(paste0("anova.H.morph.",t,".lm"))$`Pr(>F)`[1]
    if (ptemp < 0.001) {
      ptemp <- paste0("\\bf{<0.001}")
    } else if (ptemp <= 0.05) {
      ptemp <- paste0("\\bf{",signif(ptemp,2),"}")
    } else {
      ptemp <- signif(ptemp,2)
    }
    ftemp <- get(paste0("anova.H.morph.",t,".lm"))$`F value`[1]
    dftemp <- paste0(get(paste0("anova.H.morph.",t,".lm"))$Df[1],",",get(paste0("anova.H.morph.",t,".lm"))$Df[2])
    toPaste <- paste0(ptemp, " ($F_{",dftemp,"}$=",round(ftemp,2), ")")
    
    overall.anova.morphology.through.time["H",paste0(t)] <- toPaste
    
    # Now, do pairwise tests
    assign(paste0("pairwisettest.H.morph.",t) , pairwise.t.test(x = get(paste0("MF.Alpha.H.",t))[,paste0(metric)]
                                                                , g = get(paste0("MF.Alpha.H.",t))[,"Morph"]
                                                                , p.adjust.method = "BH"
    ))
    
    fbbltemp <- get(paste0("pairwisettest.H.morph.",t))$p.value[2,1]
    fbcrtemp <- get(paste0("pairwisettest.H.morph.",t))$p.value[2,2]
    blcrtemp <- get(paste0("pairwisettest.H.morph.",t))$p.value[1,1]
    
    ### MAKE BF
    for (i in c("fbbltemp","fbcrtemp","blcrtemp")) {
      if (get(paste0(i)) < 0.001) {
        assign(paste0(i), "\\bf{<0.001}")
      } else if (get(paste0(i)) <= 0.05) {
        assign(paste0(i), paste0("\\bf{",signif(get(paste0(i)),2),"}"))
      } else {
        assign(paste0(i), paste0(signif(get(paste0(i)),2)))
      }
    }
    
    pairwise.ttest.morphology.through.time["hFB:BL", paste0(t)] <- fbbltemp
    pairwise.ttest.morphology.through.time["hFB:CR", paste0(t)] <- fbcrtemp
    pairwise.ttest.morphology.through.time["hBL:CR", paste0(t)] <- blcrtemp
    
  }
  
  
  # Add overall ones
  ptemp <- get(paste0("anova.P.morph.lm.",metric))$`Pr(>F)`[3]
    if (ptemp < 0.001) {
      assign(paste0("ptemp"), "\\bf{<0.001}")
    } else if (ptemp <= 0.05) {
      assign(paste0("ptemp"), paste0("\\bf{",signif(ptemp,2),"}"))
    } else {
      assign(paste0("ptemp"), paste0(signif(ptemp,2)))
    }
  ftemp <- get(paste0("anova.P.morph.lm.",metric))$`F value`[3]
  dftemp <- paste0(get(paste0("anova.P.morph.lm.",metric))$Df[3],",",get(paste0("anova.P.morph.lm.",metric))$Df[5])
  toPaste <- paste0(ptemp, " ($F_{",dftemp,"}$=",round(ftemp,2),")")
  overall.anova.morphology.through.time["P","Overall"] <- toPaste
  
  ptemp <- get(paste0("anova.H.morph.lm.",metric))$`Pr(>F)`[3]
  if (ptemp < 0.001) {
    assign(paste0("ptemp"), "\\bf{<0.001}")
  } else if (ptemp <= 0.05) {
    assign(paste0("ptemp"), paste0("\\bf{",signif(ptemp,2),"}"))
  } else {
    assign(paste0("ptemp"), paste0(signif(ptemp,2)))
  }
  ftemp <- get(paste0("anova.H.morph.lm.",metric))$`F value`[3]
  dftemp <- paste0(get(paste0("anova.H.morph.lm.",metric))$Df[3],",",get(paste0("anova.H.morph.lm.",metric))$Df[5])
  toPaste <- paste0(ptemp, " ($F_{",dftemp,"}$=",round(ftemp,2),")")
  overall.anova.morphology.through.time["H","Overall"] <- toPaste
  
  #PAIRWISE TURN
  fbbltemp <- get(paste0("pairwise.ttest.overall.P.",metric))$p.value[2,1]
  fbcrtemp <- get(paste0("pairwise.ttest.overall.P.",metric))$p.value[2,2]
  blcrtemp <- get(paste0("pairwise.ttest.overall.P.",metric))$p.value[1,1]
  ### MAKE BF
  for (i in c("fbbltemp","fbcrtemp","blcrtemp")) {
    if (get(paste0(i)) < 0.001) {
      assign(paste0(i), "\\bf{<0.001}")
    } else if (get(paste0(i)) <= 0.05) {
      assign(paste0(i), paste0("\\bf{",signif(get(paste0(i)),2),"}"))
    } else {
      assign(paste0(i), paste0(signif(get(paste0(i)),2)))
    }
  }
  
  pairwise.ttest.morphology.through.time["pFB:BL","Overall"] <- fbbltemp
  pairwise.ttest.morphology.through.time["pFB:CR","Overall"] <- fbcrtemp
  pairwise.ttest.morphology.through.time["pBL:CR","Overall"] <- blcrtemp
  
  fbbltemp <- get(paste0("pairwise.ttest.overall.H.",metric))$p.value[2,1]
  fbcrtemp <- get(paste0("pairwise.ttest.overall.H.",metric))$p.value[2,2]
  blcrtemp <- get(paste0("pairwise.ttest.overall.H.",metric))$p.value[1,1]
  ### MAKE BF
  for (i in c("fbbltemp","fbcrtemp","blcrtemp")) {
    if (get(paste0(i)) < 0.001) {
      assign(paste0(i), "\\bf{<0.001}")
    } else if (get(paste0(i)) <= 0.05) {
      assign(paste0(i), paste0("\\bf{",signif(get(paste0(i)),2),"}"))
    } else {
      assign(paste0(i), paste0(signif(get(paste0(i)),2)))
    }
  }
  pairwise.ttest.morphology.through.time["hFB:BL","Overall"] <- fbbltemp
  pairwise.ttest.morphology.through.time["hFB:CR","Overall"] <- fbcrtemp
  pairwise.ttest.morphology.through.time["hBL:CR","Overall"] <- blcrtemp
  
  
  # Round the numbers to make sure all are appropriate length
  tempfile <- overall.anova.morphology.through.time
  for (r in 1:nrow(tempfile)) {
    for (c in 1:ncol(tempfile)) {
      if (is.na(tempfile[r,c])) {
        overall.anova.morphology.through.time[r,c] <- "-"
      }}}
  
  tempfile2 <- pairwise.ttest.morphology.through.time
  for (r in 1:nrow(tempfile2)) {
    for (c in 1:ncol(tempfile2)) {
      if (is.na(tempfile2[r,c])) {
        pairwise.ttest.morphology.through.time[r,c] <- "-"
      } 
    }
  }
  
  # Change colnames to something nice
  colnames(overall.anova.morphology.through.time) <- c("20 minutes","1 hour","3 hours","6 hours","12 hours","1 day","4 days", "Overall")
  colnames(pairwise.ttest.morphology.through.time) <- c("20 minutes","1 hour","3 hours","6 hours","12 hours","1 day","4 days", "Overall")
  
  # Change rownames to something nice
  rownames(overall.anova.morphology.through.time) <- c("Reed Point","Hakai")
  assign(paste0("pairwise.ttest.morphology.through.time.",metric),cbind(c("Reed Point",""," ","Hakai","  ","   "),Factors = rep(c("FB:BL","FB:CR","BL:CR"),2),pairwise.ttest.morphology.through.time ))
  assign(paste0("pairwise.ttest.morphology.through.time.",metric), data.frame(get(paste0("pairwise.ttest.morphology.through.time.",metric)), row.names = 1, check.names = FALSE))
  assign(paste0("overall.anova.morphology.through.time.",metric), overall.anova.morphology.through.time)

  capture.output(print(xtable(get(paste0("overall.anova.morphology.through.time.",metric)))
                 , file = paste0("ALPHAPLOTS_latex/overall_anova_morphology_through_time_H",metric,".txt")
                 , math.style.negative = TRUE
                 , math.style.exponents = TRUE
                 , sanitize.text.function = function(x) {x}))
  capture.output(print(xtable(get(paste0("pairwise.ttest.morphology.through.time.",metric)))
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x})
                 , file = paste0("ALPHAPLOTS_latex/pairwise_ttest_morphology_through_time_H",metric,".txt"))
  
  
  ########### Plotting Hakai ##############
  
  H.CR.Alpha <- MF.Alpha.morphonly.H[grep("CR", MF.Alpha.morphonly.H$Morph),]
  H.BL.Alpha <- MF.Alpha.morphonly.H[grep("BL", MF.Alpha.morphonly.H$Morph),]
  H.FB.Alpha <- MF.Alpha.morphonly.H[grep("FB", MF.Alpha.morphonly.H$Morph),]
  
  H.CR.mean <- aggregate(H.CR.Alpha[,paste0(metric)], by=list(Time=H.CR.Alpha$Time), FUN=mean)
  H.CR.sd <- aggregate(H.CR.Alpha[,paste0(metric)], by=list(Time=H.CR.Alpha$Time), FUN=sd)
  H.BL.mean <- aggregate(H.BL.Alpha[,paste0(metric)], by=list(Time=H.BL.Alpha$Time), FUN=mean)
  H.BL.sd <- aggregate(H.BL.Alpha[,paste0(metric)], by=list(Time=H.BL.Alpha$Time), FUN=sd)
  H.FB.mean <- aggregate(H.FB.Alpha[,paste0(metric)], by=list(Time=H.FB.Alpha$Time), FUN=mean)
  H.FB.sd <- aggregate(H.FB.Alpha[,paste0(metric)], by=list(Time=H.FB.Alpha$Time), FUN=sd)
  
  TimeChar <- H.CR.mean$Time[1:4]
  
  pdf(paste0("ALPHAPLOTS_H/Alpha_diversity_Hakai_",metric,".pdf"), pointsize = 14)
  par(fig = c(0,0.7,0,1))
  plot(c(0,0,0,0) ~ TimeChar[1:4]
       , main = paste0(metric)
       , ylim = c(min(MF.Alpha.morphonly.H[,paste0(metric)]), max(MF.Alpha.morphonly.H[,paste0(metric)]))
       , xlab = "Time"
       , ylab = paste0("Richness (",metric,")")
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
  dev.off()
  
  ########### Plotting PM ##############
  
  P.CR.Alpha <- MF.Alpha.morphonly.P[grep("CR", MF.Alpha.morphonly.P$Morph),]
  P.BL.Alpha <- MF.Alpha.morphonly.P[grep("BL", MF.Alpha.morphonly.P$Morph),]
  P.FB.Alpha <- MF.Alpha.morphonly.P[grep("FB", MF.Alpha.morphonly.P$Morph),]
  
  P.CR.mean <- aggregate(P.CR.Alpha[,paste0(metric)], by=list(Time=P.CR.Alpha$Time), FUN=mean)
  P.CR.sd <- aggregate(P.CR.Alpha[,paste0(metric)], by=list(Time=P.CR.Alpha$Time), FUN=sd)
  P.BL.mean <- aggregate(P.BL.Alpha[,paste0(metric)], by=list(Time=P.BL.Alpha$Time), FUN=mean)
  P.BL.sd <- aggregate(P.BL.Alpha[,paste0(metric)], by=list(Time=P.BL.Alpha$Time), FUN=sd)
  P.FB.mean <- aggregate(P.FB.Alpha[,paste0(metric)], by=list(Time=P.FB.Alpha$Time), FUN=mean)
  P.FB.sd <- aggregate(P.FB.Alpha[,paste0(metric)], by=list(Time=P.FB.Alpha$Time), FUN=sd)
  
  TimeChar <- P.CR.mean$Time[1:6]
  
  pdf(paste0("ALPHAPLOTS_P/Alpha_diversity_PM_",metric,".pdf"), pointsize = 14)
  par(fig = c(0,0.7,0,1))
  plot(c(0,0,0,0,0,0) ~ TimeChar[1:6]
       , main = paste0(metric)
       , ylim = c(min(MF.Alpha.morphonly.P[,paste0(metric)]), max(MF.Alpha.morphonly.P[,paste0(metric)]))
       , xlab = "Time"
       , ylab = paste0("Richness (",metric,")")
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
  
  dev.off()
  
}

############# ******WATER COMPARISONS****** ##############
allwaters <- MF.Alpha.morphwater[grep("(H_W)|(P_H2)", MF.Alpha.morphwater$TypeMorphTime),]
for (i in alphaList) {
  temp.lm <- lm(as.numeric(allwaters[,i]) ~ allwaters$Type)
  anova.allwaters <- Anova(temp.lm, type = "III")
  capture.output(xtable(anova.allwaters, digits = 3), file = paste0("./ALPHA_compare/wateronly_",i,".txt"))

}

############ *********Water vs substrate ********############
# PM

MF.Alpha.P.mw <- MF.Alpha.morphwater[MF.Alpha.morphwater$Type == "P",]

for (i in alphaList) {
  met <- gsub("_even_4000_alpha","", i)
  water.alpha.P <- as.numeric(MF.Alpha.P.mw[grep("P_H", MF.Alpha.P.mw$TypeMorphTime), i])
  AM.alpha.P <- MF.Alpha.P.mw[-grep("P_H", MF.Alpha.P.mw$TypeMorphTime), ]
  alphaCompare.P <- matrix(ncol = 6, nrow = 4)
  rownames(alphaCompare.P) <- c("mean","sd","ttest","ttest2")
  colnames(alphaCompare.P) <- c("20","60","180","360","720","1440")
  for (t in c("20","60","180","360","720","1440")) {
    alphaTemp <- as.numeric(AM.alpha.P[AM.alpha.P$Time == t,i])
    alphaCompare.P["mean",t] <- mean(alphaTemp)
    alphaCompare.P["sd", t] <- sd(alphaTemp)
    
    ttestTemp <-  t.test(as.numeric(water.alpha.P), as.numeric(alphaTemp))
    ptemp <- signif(ttestTemp$p.value,2)
    ttemp <- round(ttestTemp$statistic,2)
    dftemp <- round(ttestTemp$parameter,2)
    toPaste <- paste0(" (t=",ttemp,", df=",dftemp,")")
    
    alphaCompare.P["ttest",t] <- paste0("p=",ptemp)
    alphaCompare.P["ttest2",t] <- toPaste
  }
  
  minD <- min(as.numeric(AM.alpha.P[,i]))
  maxD <- max(as.numeric(AM.alpha.P[,i]))
  pdf(paste0("ALPHA_compare/watervssubstrate_P_",met,".pdf"), pointsize = 14)
  par(mar = c(5.5,4,4,4))
  plot(alphaCompare.P["mean",]
       , pch = 19
       , col = "black"
       , xlab = ""
       , ylab = paste0("Diversity (",met,")")
       , xaxt = "n"
       , ylim = c(minD,maxD)
  )
  title(xlab = "Time"
        , line = 4.5)
  axis(side = 1
       , at = c(1,2,3,4,5,6)
       , labels = c("20 m","1 h","3 h", "6 h","12 h","1 d"))
  arrows(x0 = c(1,2,3,4,5,6)
         , y0 = as.numeric(alphaCompare.P["mean",])- as.numeric(alphaCompare.P["sd",])
         , y1 = as.numeric(alphaCompare.P["mean",])+ as.numeric(alphaCompare.P["sd",])
         , angle = 90
         , code = 3)
  abline(h = mean(water.alpha.P) + c(sd(water.alpha.P),0,-sd(water.alpha.P))
         , lty = c("dotted","solid","dotted")
         , col = "blue")
  axis(side = 1
       , at = c(1,2,3,4,5,6)
       , labels = alphaCompare.P["ttest",]
       , cex.axis = 0.7
       , line = 1
       , lwd = 0)
  axis(side = 1
       , at = c(1,2,3,4,5,6)
       , labels = alphaCompare.P["ttest2",]
       , cex.axis = 0.5
       , line = 2
       , lwd = 0)
  dev.off()
}

## HAKAI

MF.Alpha.H.mw <- MF.Alpha.morphwater[MF.Alpha.morphwater$Type == "H",]

for (i in alphaList) {
  met <- gsub("_even_4000_alpha","", i)
  water.alpha.H <- as.numeric(MF.Alpha.H.mw[grep("H_W", MF.Alpha.H.mw$TypeMorphTime), i])
  AM.alpha.H <- MF.Alpha.H.mw[-grep("H_W", MF.Alpha.H.mw$TypeMorphTime), ]
  alphaCompare.H <- matrix(ncol = 5, nrow = 4)
  rownames(alphaCompare.H) <- c("mean","sd","ttest","ttest2")
  colnames(alphaCompare.H) <- c("20","60","360","720","5760")
  for (t in c("20","60","360","720","5760")) {
    alphaTemp <- as.numeric(AM.alpha.H[AM.alpha.H$Time == t,i])
    alphaCompare.H["mean",t] <- mean(alphaTemp)
    alphaCompare.H["sd", t] <- sd(alphaTemp)
    
    ttestTemp <-  t.test(as.numeric(water.alpha.H), as.numeric(alphaTemp))
    ptemp <- signif(ttestTemp$p.value,2)
    ttemp <- round(ttestTemp$statistic,2)
    dftemp <- round(ttestTemp$parameter,2)
    toPaste <- paste0(" (t=",ttemp,", df=",dftemp,")")
    
    alphaCompare.H["ttest",t] <- paste0("p=",ptemp)
    alphaCompare.H["ttest2",t] <- toPaste
  }
  
  minD <- min(as.numeric(AM.alpha.H[,i]))
  maxD <- max(as.numeric(AM.alpha.H[,i]))
  pdf(paste0("ALPHA_compare/watervssubstrate_H_",met,".pdf"), pointsize = 14)
  par(mar = c(5.5,4,4,4))
  plot(alphaCompare.H["mean",]
       , pch = 19
       , col = "black"
       , xlab = ""
       , ylab = paste0("Diversity (",met,")")
       , xaxt = "n"
       , ylim = c(minD,maxD)
  )
  title(xlab = "Time"
        , line = 4.5)
  axis(side = 1
       , at = c(1,2,3,4,5)
       , labels = c("20 m","1 h", "6 h","12 h","4 d"))
  arrows(x0 = c(1,2,3,4,5)
         , y0 = as.numeric(alphaCompare.H["mean",])- as.numeric(alphaCompare.H["sd",])
         , y1 = as.numeric(alphaCompare.H["mean",])+ as.numeric(alphaCompare.H["sd",])
         , angle = 90
         , code = 3)
  abline(h = mean(water.alpha.H) + c(sd(water.alpha.H),0,-sd(water.alpha.H))
         , lty = c("dotted","solid","dotted")
         , col = "blue")
  axis(side = 1
       , at = c(1,2,3,4,5)
       , labels = alphaCompare.H["ttest",]
       , cex.axis = 1
       , line = 1
       , lwd = 0)
  axis(side = 1
       , at = c(1,2,3,4,5)
       , labels = alphaCompare.H["ttest2",]
       , cex.axis = 0.5
       , line = 2
       , lwd = 0)
  dev.off()
}

  
  

############# ******PER MORPH COMPARISONS****** ##############

#### HAKAI vs PM DIVERSITY ####

for (metric in c("chao1","PD_whole_tree","observed_otus")) {
  # Empty matrix for results; name the headers
  results.HPM.adivcompare <- matrix(ncol = 5, nrow = 3)
  rownames(results.HPM.adivcompare) <- c("p-value","t-statistic","Df")
  colnames(results.HPM.adivcompare) <- c("20","60","360","720", "Overall")
  # Iterate through all shared timepoints
  for (i in c("20","60","360","720")) {
    assign(paste0("ttest.HvsP.",i), t.test(get(paste0("MF.Alpha.H.",i))[,metric], get(paste0("MF.Alpha.P.",i))[,metric]))
    results.HPM.adivcompare[ "p-value", paste0(i)] <- get(paste0("ttest.HvsP.",i))$p.value
    results.HPM.adivcompare["t-statistic",paste0(i)] <- get(paste0("ttest.HvsP.",i))$statistic
    results.HPM.adivcompare["Df",paste0(i)] <- get(paste0("ttest.HvsP.",i))$parameter
  }
  # Add an overall test too
  ttest.HvsP.all <- t.test(MF.Alpha.morphonly.H[,metric], MF.Alpha.morphonly.P[,metric])
  results.HPM.adivcompare["p-value", "Overall"] <- ttest.HvsP.all$p.value
  results.HPM.adivcompare["t-statistic", "Overall"] <- ttest.HvsP.all$statistic
  results.HPM.adivcompare["Df","Overall"] <- ttest.HvsP.all$parameter
  colnames(results.HPM.adivcompare) <- c("20 min","1 hour","6 hours","12 hours","Overall")
  # change p < 0.001 to character so formatting makes sense
  temp.pvalues <- results.HPM.adivcompare
  for (r in 1:nrow(temp.pvalues)) {
    for (c in 1:ncol(temp.pvalues)) {
      if (temp.pvalues[r,c] < 0.001) {
        results.HPM.adivcompare[r,c] <- "<0.001"
      } else (
        results.HPM.adivcompare[r,c] <- as.character(round(temp.pvalues[r,c], digits = 3))
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

############# ********INDIVIDUAL COMBINED************ ################

# Combine into master
overall.anova.morph.through.time.MASTER <- matrix(ncol = 8, nrow = 6)
colnames(overall.anova.morph.through.time.MASTER) <- colnames(overall.anova.morphology.through.time.chao1)
rownames(overall.anova.morph.through.time.MASTER) <- c("pchao1","pPD_whole_tree","pobserved_otus","hchao1","hPD_whole_tree","hobserved_otus")

for (met in metrics) {
  PMTemp <- get(paste0("overall.anova.morphology.through.time.",met))[1,]
  HKTemp <- get(paste0("overall.anova.morphology.through.time.",met))[2,]
  
  overall.anova.morph.through.time.MASTER[paste0("p",met),] <- PMTemp
  overall.anova.morph.through.time.MASTER[paste0("h",met),] <- HKTemp
}

# Shift names

overall.anova.morph.through.time.MASTER.edit <- cbind(rep(c("Chao1","PD_whole_tree","Observed_OTUs"),2), overall.anova.morph.through.time.MASTER)
rownames(overall.anova.morph.through.time.MASTER.edit) <- c("Reed Point",""," ","Hakai","  ","   ")

capture.output(xtable(overall.anova.morph.through.time.MASTER.edit), file = "ALPHAPLOTS_latex/overall.anova.morph.through.time.MASTER.txt")

### ONLY Chao1 into full table

# Add extra column to overall
overallTemp <- cbind(Factors = c("Morph","Morph"), overall.anova.morphology.through.time.chao1)

# Now combine...
LATEX_overall_pairwise_combo_chao1 <- rbind(overallTemp[1,]
                                            ,  apply(pairwise.ttest.morphology.through.time.chao1[1:3,], MARGIN = 2, FUN = factor)
                                            , overallTemp[2,]
                                            , apply(pairwise.ttest.morphology.through.time.chao1[4:6,], MARGIN = 2, FUN = factor)
                                            )
LATEX_overall_pairwise_combo_chao1 <- cbind("Sites" = c("Reed Point",""," ","   ", "Hakai","","    ","     ")
                                            , "Test type" = c("ANOVA","","Pairwise t-tests","","ANOVA","","Pairwise t-tests","")
                                            , LATEX_overall_pairwise_combo_chao1)

LATEX_overall_pairwise_combo_chao1 <- gsub("^","\\\\makecell{", LATEX_overall_pairwise_combo_chao1)
LATEX_overall_pairwise_combo_chao1 <- gsub("$","}", LATEX_overall_pairwise_combo_chao1)
LATEX_overall_pairwise_combo_chao1 <- gsub("(","\\\\(", LATEX_overall_pairwise_combo_chao1, fixed = TRUE)



capture.output(print(xtable(LATEX_overall_pairwise_combo_chao1)
                     , include.rownames = FALSE
                     , math.style.negative = TRUE
                     , math.style.exponents = TRUE
                     , sanitize.text.function = function(x) {x})
               , file = "ALPHAPLOTS_latex/LATEX_overall_pairwise_combo_chao1.txt")

                                            

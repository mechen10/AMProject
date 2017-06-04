#!/bin/Rscript

# setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/ALPHA_R")
library("car")
library("nlme")
library("optparse")
library("gridExtra")
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

# setwd('/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis')
# MPFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/TEMP_frombotclust/MF_withalpha.txt'
# alphaNames = 'chao1_even_4000_normalized_alpha,PD_whole_tree_even_4000_normalized_alpha,observed_otus_even_4000_normalized_alpha'
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
  if (MF.Alpha.all[i,'Morph'] %in% c("CR","BL","FB","W")) {
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

########################### STATS #################################
system('mkdir ALPHAPLOTS_H')
system('mkdir ALPHAPLOTS_P')

########### ****CHAO1**** ##############
metric <- 'chao1'

########### Stats ##############

H.morph.lm <- lm(chao1 ~ Time*Morph, data = MF.Alpha.morphonly.H)
anova.H.morph.lm <- Anova(H.morph.lm, type = 3)
capture.output(anova.H.morph.lm, file = paste0("ALPHAPLOTS_H/anova_Hakai_", metric, ".txt"))
P.morph.lm <- lm(chao1 ~ Time*Morph, data = MF.Alpha.morphonly.P)
anova.P.morph.lm <- Anova(P.morph.lm, type = 3)
capture.output(anova.P.morph.lm, file = paste0("ALPHAPLOTS_P/anova_PM_", metric, ".txt"))

# Separated out Stats- Hakai
MF.Alpha.H.20 <- MF.Alpha.morphonly.H[grep("^20$", MF.Alpha.morphonly.H$Time),]
MF.Alpha.H.60 <- MF.Alpha.morphonly.H[grep("^60$", MF.Alpha.morphonly.H$Time),]
MF.Alpha.H.360 <- MF.Alpha.morphonly.H[grep("^360$", MF.Alpha.morphonly.H$Time),]
MF.Alpha.H.720 <- MF.Alpha.morphonly.H[grep("^720$", MF.Alpha.morphonly.H$Time),]

# Anova and capture output
H.morph.20.lm <- lm(chao1 ~ Morph, data = MF.Alpha.H.20)
anova.H.morph.20.lm <- anova(H.morph.20.lm)
capture.output(anova.H.morph.20.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_20.txt"))

H.morph.60.lm <- lm(chao1 ~ Morph, data = MF.Alpha.H.60)
anova.H.morph.60.lm <- anova(H.morph.60.lm)
capture.output(anova.H.morph.60.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_60.txt"))

H.morph.360.lm <- lm(chao1 ~ Morph, data = MF.Alpha.H.360)
anova.H.morph.360.lm <- anova(H.morph.360.lm)
capture.output(anova.H.morph.360.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_360.txt"))

H.morph.720.lm <- lm(chao1 ~ Morph, data = MF.Alpha.H.720)
anova.H.morph.720.lm <- anova(H.morph.720.lm)
capture.output(anova.H.morph.720.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_720.txt"))

# Separated out Stats- PM
MF.Alpha.P.20 <- MF.Alpha.morphonly.P[grep("^20$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.60 <- MF.Alpha.morphonly.P[grep("^60$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.180 <- MF.Alpha.morphonly.P[grep("^180$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.360 <- MF.Alpha.morphonly.P[grep("^360$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.720 <- MF.Alpha.morphonly.P[grep("^720$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.1440 <- MF.Alpha.morphonly.P[grep("^1440$", MF.Alpha.morphonly.P$Time),]


# Anova and capture output
P.morph.20.lm <- lm(chao1 ~ Morph, data = MF.Alpha.P.20)
anova.P.morph.20.lm <- anova(P.morph.20.lm)
capture.output(anova.P.morph.20.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_20.txt"))

P.morph.60.lm <- lm(chao1 ~ Morph, data = MF.Alpha.P.60)
anova.P.morph.60.lm <- anova(P.morph.60.lm)
capture.output(anova.P.morph.60.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_60.txt"))

P.morph.180.lm <- lm(chao1 ~ Morph, data = MF.Alpha.P.180)
anova.P.morph.180.lm <- anova(P.morph.180.lm)
capture.output(anova.P.morph.180.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_180.txt"))

P.morph.360.lm <- lm(chao1 ~ Morph, data = MF.Alpha.P.360)
anova.P.morph.360.lm <- anova(P.morph.360.lm)
capture.output(anova.P.morph.360.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_360.txt"))

P.morph.720.lm <- lm(chao1 ~ Morph, data = MF.Alpha.P.720)
anova.P.morph.720.lm <- anova(P.morph.720.lm)
capture.output(anova.P.morph.720.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_1440.txt"))

P.morph.1440.lm <- lm(chao1 ~ Morph, data = MF.Alpha.P.1440)
anova.P.morph.1440.lm <- anova(P.morph.1440.lm)
capture.output(anova.P.morph.1440.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_1440.txt"))

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

# TimeChar <- as.character(CR$Time)[1:4]
# TimeChar <- log(H.CR.mean$Time[1:4],10)
TimeChar <- H.CR.mean$Time[1:4]

# # Make table for stats
# p <- c(format(anova.H.morph.20.lm$`Pr(>F)`[1],digits = 2)
#   , format(anova.H.morph.60.lm$`Pr(>F)`[1],digits = 2)
#   , format(anova.H.morph.360.lm$`Pr(>F)`[1],digits = 2)
#   , format(anova.H.morph.720.lm$`Pr(>F)`[1],digits = 2)
#   , format(anova.H.morph.lm$`Pr(>F)`[4], digits = 2)
#   )
# f <- c(format(anova.H.morph.20.lm$`F value`[1],digits = 2)
#              , format(anova.H.morph.60.lm$`F value`[1],digits = 2)
#              , format(anova.H.morph.360.lm$`F value`[1],digits = 2)
#              , format(anova.H.morph.720.lm$`F value`[1],digits = 2)
#              , format(anova.H.morph.lm$`F value`[4], digits = 2)
# )
# statsTable <- cbind(p,f)
# rownames(statsTable) <- c("20 min"
#                           , "1 h"
#                           , "6 h"
#                           , "12 h"
#                           , "Time:Morph")

# Min and max values
# TimeChar <- factor(TimeChar, levels = c("20","60","360","720"))
pdf(paste0("ALPHAPLOTS_H/Alpha_diversity_Hakai_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0) ~ TimeChar[1:4]
     , main = "Richness through time"
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
       # , col = 'darkgreen'
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
       # , col = "green"
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
       # , col = "yellow"
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
     , labels = c(paste0("20 min: "))
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

# TimeChar <- as.character(CR$Time)[1:4]
# TimeChar <- log(P.CR.mean$Time[1:4],10)
TimeChar <- P.CR.mean$Time[1:6]


# Min and max values
# TimeChar <- factor(TimeChar, levels = c("20","60","360","720"))
pdf(paste0("ALPHAPLOTS_P/Alpha_diversity_PM_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0,0,0) ~ TimeChar[1:6]
     , main = "Richness through time"
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
       # , col = 'darkgreen'
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
       # , col = "green"
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
       # , col = "yellow"
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

H.morph.lm <- lm(PD_whole_tree ~ Time*Morph, data = MF.Alpha.morphonly.H)
anova.H.morph.lm <- Anova(H.morph.lm, type = 3)
capture.output(anova.H.morph.lm, file = paste0("ALPHAPLOTS_H/anova_Hakai_", metric, ".txt"))
P.morph.lm <- lm(PD_whole_tree ~ Time*Morph, data = MF.Alpha.morphonly.P)
anova.P.morph.lm <- Anova(P.morph.lm, type = 3)
capture.output(anova.P.morph.lm, file = paste0("ALPHAPLOTS_P/anova_PM_", metric, ".txt"))

# Separated out Stats- Hakai
MF.Alpha.H.20 <- MF.Alpha.morphonly.H[grep("^20$", MF.Alpha.morphonly.H$Time),]
MF.Alpha.H.60 <- MF.Alpha.morphonly.H[grep("^60$", MF.Alpha.morphonly.H$Time),]
MF.Alpha.H.360 <- MF.Alpha.morphonly.H[grep("^360$", MF.Alpha.morphonly.H$Time),]
MF.Alpha.H.720 <- MF.Alpha.morphonly.H[grep("^720$", MF.Alpha.morphonly.H$Time),]

# Anova and capture output
H.morph.20.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.H.20)
anova.H.morph.20.lm <- anova(H.morph.20.lm)
capture.output(anova.H.morph.20.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_20.txt"))

H.morph.60.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.H.60)
anova.H.morph.60.lm <- anova(H.morph.60.lm)
capture.output(anova.H.morph.60.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_60.txt"))

H.morph.360.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.H.360)
anova.H.morph.360.lm <- anova(H.morph.360.lm)
capture.output(anova.H.morph.360.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_360.txt"))

H.morph.720.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.H.720)
anova.H.morph.720.lm <- anova(H.morph.720.lm)
capture.output(anova.H.morph.720.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_720.txt"))

# Separated out Stats- PM
MF.Alpha.P.20 <- MF.Alpha.morphonly.P[grep("^20$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.60 <- MF.Alpha.morphonly.P[grep("^60$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.180 <- MF.Alpha.morphonly.P[grep("^180$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.360 <- MF.Alpha.morphonly.P[grep("^360$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.720 <- MF.Alpha.morphonly.P[grep("^720$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.1440 <- MF.Alpha.morphonly.P[grep("^1440$", MF.Alpha.morphonly.P$Time),]


# Anova and capture output
P.morph.20.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.P.20)
anova.P.morph.20.lm <- anova(P.morph.20.lm)
capture.output(anova.P.morph.20.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_20.txt"))

P.morph.60.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.P.60)
anova.P.morph.60.lm <- anova(P.morph.60.lm)
capture.output(anova.P.morph.60.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_60.txt"))

P.morph.180.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.P.180)
anova.P.morph.180.lm <- anova(P.morph.180.lm)
capture.output(anova.P.morph.180.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_180.txt"))

P.morph.360.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.P.360)
anova.P.morph.360.lm <- anova(P.morph.360.lm)
capture.output(anova.P.morph.360.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_360.txt"))

P.morph.720.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.P.720)
anova.P.morph.720.lm <- anova(P.morph.720.lm)
capture.output(anova.P.morph.720.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_1440.txt"))

P.morph.1440.lm <- lm(PD_whole_tree ~ Morph, data = MF.Alpha.P.1440)
anova.P.morph.1440.lm <- anova(P.morph.1440.lm)
capture.output(anova.P.morph.1440.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_1440.txt"))

########### Plotting Hakai ##############

H.CR.Alpha <- MF.Alpha.morphonly.H[grep("CR", MF.Alpha.morphonly.H$Morph),]
H.BL.Alpha <- MF.Alpha.morphonly.H[grep("BL", MF.Alpha.morphonly.H$Morph),]
H.FB.Alpha <- MF.Alpha.morphonly.H[grep("FB", MF.Alpha.morphonly.H$Morph),]

H.CR.mean <- aggregate(H.CR.Alpha[,paste0(metric)], by=list(Time=H.CR.Alpha$Time), FUN=mean)
H.CR.sd <- aggregate(H.CR.Alpha$PD_whole_tree, by=list(Time=H.CR.Alpha$Time), FUN=sd)
H.BL.mean <- aggregate(H.BL.Alpha$PD_whole_tree, by=list(Time=H.BL.Alpha$Time), FUN=mean)
H.BL.sd <- aggregate(H.BL.Alpha$PD_whole_tree, by=list(Time=H.BL.Alpha$Time), FUN=sd)
H.FB.mean <- aggregate(H.FB.Alpha$PD_whole_tree, by=list(Time=H.FB.Alpha$Time), FUN=mean)
H.FB.sd <- aggregate(H.FB.Alpha$PD_whole_tree, by=list(Time=H.FB.Alpha$Time), FUN=sd)

# TimeChar <- as.character(CR$Time)[1:4]
# TimeChar <- log(H.CR.mean$Time[1:4],10)
TimeChar <- H.CR.mean$Time[1:4]

# Min and max values
# TimeChar <- factor(TimeChar, levels = c("20","60","360","720"))
pdf(paste0("ALPHAPLOTS_H/Alpha_diversity_Hakai_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0) ~ TimeChar[1:4]
     , main = "Richness through time"
     , ylim = c(min(MF.Alpha.morphonly.H$PD_whole_tree), max(MF.Alpha.morphonly.H$PD_whole_tree))
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
       # , col = 'darkgreen'
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
       # , col = "green"
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
       # , col = "yellow"
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
     , labels = c(paste0("20 min: "))
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
     , labels = c(paste0(format(anova.H.morph.720.lm$`Pr(>F)`[1], digits = 2), '*'))
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
P.CR.sd <- aggregate(P.CR.Alpha$PD_whole_tree, by=list(Time=P.CR.Alpha$Time), FUN=sd)
P.BL.mean <- aggregate(P.BL.Alpha$PD_whole_tree, by=list(Time=P.BL.Alpha$Time), FUN=mean)
P.BL.sd <- aggregate(P.BL.Alpha$PD_whole_tree, by=list(Time=P.BL.Alpha$Time), FUN=sd)
P.FB.mean <- aggregate(P.FB.Alpha$PD_whole_tree, by=list(Time=P.FB.Alpha$Time), FUN=mean)
P.FB.sd <- aggregate(P.FB.Alpha$PD_whole_tree, by=list(Time=P.FB.Alpha$Time), FUN=sd)

# TimeChar <- as.character(CR$Time)[1:4]
# TimeChar <- log(P.CR.mean$Time[1:4],10)
TimeChar <- P.CR.mean$Time[1:6]


# Min and max values
# TimeChar <- factor(TimeChar, levels = c("20","60","360","720"))
pdf(paste0("ALPHAPLOTS_P/Alpha_diversity_PM_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0,0,0) ~ TimeChar[1:6]
     , main = "Richness through time"
     , ylim = c(min(MF.Alpha.morphonly.P$PD_whole_tree), max(MF.Alpha.morphonly.P$PD_whole_tree))
     , xlab = "Time"
     , ylab = "Richness (PD_whole_tree)"
     , xaxt = 'n'
)
axis(1
     , at = c(-20,10,0,0,0,0) + TimeChar[1:6]
     , labels = c("20m","1h","3h","6h","12h","24h")
     , las = 2
     , cex.axis = 1)
points(TimeChar[1:6], P.CR.mean$x[1:6]
       , lty = 3
       # , col = 'darkgreen'
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
       # , col = "green"
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
       # , col = "yellow"
       , type = 'l')
arrows(x0 = TimeChar*1.01
       , x1 = TimeChar*1.01
       , y0 = c(P.FB.mean$x - P.FB.sd$x/2)[1:6]
       , y1 = c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02)
text(TimeChar, 1.05*c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
     , labels = c('','**','*','','*','')
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
     , labels = c(paste0(format(anova.P.morph.20.lm$`Pr(>F)`[1],digits = 2), ""))
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

H.morph.lm <- lm(observed_otus ~ Time*Morph, data = MF.Alpha.morphonly.H)
anova.H.morph.lm <- Anova(H.morph.lm, type = 3)
capture.output(anova.H.morph.lm, file = paste0("ALPHAPLOTS_H/anova_Hakai_", metric, ".txt"))
P.morph.lm <- lm(observed_otus ~ Time*Morph, data = MF.Alpha.morphonly.P)
anova.P.morph.lm <- Anova(P.morph.lm, type = 3)
capture.output(anova.P.morph.lm, file = paste0("ALPHAPLOTS_P/anova_PM_", metric, ".txt"))

# Separated out Stats- Hakai
MF.Alpha.H.20 <- MF.Alpha.morphonly.H[grep("^20$", MF.Alpha.morphonly.H$Time),]
MF.Alpha.H.60 <- MF.Alpha.morphonly.H[grep("^60$", MF.Alpha.morphonly.H$Time),]
MF.Alpha.H.360 <- MF.Alpha.morphonly.H[grep("^360$", MF.Alpha.morphonly.H$Time),]
MF.Alpha.H.720 <- MF.Alpha.morphonly.H[grep("^720$", MF.Alpha.morphonly.H$Time),]

# Anova and capture output
H.morph.20.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.H.20)
anova.H.morph.20.lm <- anova(H.morph.20.lm)
capture.output(anova.H.morph.20.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_20.txt"))

H.morph.60.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.H.60)
anova.H.morph.60.lm <- anova(H.morph.60.lm)
capture.output(anova.H.morph.60.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_60.txt"))

H.morph.360.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.H.360)
anova.H.morph.360.lm <- anova(H.morph.360.lm)
capture.output(anova.H.morph.360.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_360.txt"))

H.morph.720.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.H.720)
anova.H.morph.720.lm <- anova(H.morph.720.lm)
capture.output(anova.H.morph.720.lm, file = paste0("ALPHAPLOTS_H/anova_sep_H_",metric,"_720.txt"))

# Separated out Stats- PM
MF.Alpha.P.20 <- MF.Alpha.morphonly.P[grep("^20$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.60 <- MF.Alpha.morphonly.P[grep("^60$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.180 <- MF.Alpha.morphonly.P[grep("^180$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.360 <- MF.Alpha.morphonly.P[grep("^360$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.720 <- MF.Alpha.morphonly.P[grep("^720$", MF.Alpha.morphonly.P$Time),]
MF.Alpha.P.1440 <- MF.Alpha.morphonly.P[grep("^1440$", MF.Alpha.morphonly.P$Time),]


# Anova and capture output
P.morph.20.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.P.20)
anova.P.morph.20.lm <- anova(P.morph.20.lm)
capture.output(anova.P.morph.20.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_20.txt"))

P.morph.60.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.P.60)
anova.P.morph.60.lm <- anova(P.morph.60.lm)
capture.output(anova.P.morph.60.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_60.txt"))

P.morph.180.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.P.180)
anova.P.morph.180.lm <- anova(P.morph.180.lm)
capture.output(anova.P.morph.180.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_180.txt"))

P.morph.360.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.P.360)
anova.P.morph.360.lm <- anova(P.morph.360.lm)
capture.output(anova.P.morph.360.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_360.txt"))

P.morph.720.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.P.720)
anova.P.morph.720.lm <- anova(P.morph.720.lm)
capture.output(anova.P.morph.720.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_1440.txt"))

P.morph.1440.lm <- lm(observed_otus ~ Morph, data = MF.Alpha.P.1440)
anova.P.morph.1440.lm <- anova(P.morph.1440.lm)
capture.output(anova.P.morph.1440.lm, file = paste0("ALPHAPLOTS_P/anova_sep_P_",metric,"_1440.txt"))

########### Plotting Hakai ##############

H.CR.Alpha <- MF.Alpha.morphonly.H[grep("CR", MF.Alpha.morphonly.H$Morph),]
H.BL.Alpha <- MF.Alpha.morphonly.H[grep("BL", MF.Alpha.morphonly.H$Morph),]
H.FB.Alpha <- MF.Alpha.morphonly.H[grep("FB", MF.Alpha.morphonly.H$Morph),]

H.CR.mean <- aggregate(H.CR.Alpha[,paste0(metric)], by=list(Time=H.CR.Alpha$Time), FUN=mean)
H.CR.sd <- aggregate(H.CR.Alpha$observed_otus, by=list(Time=H.CR.Alpha$Time), FUN=sd)
H.BL.mean <- aggregate(H.BL.Alpha$observed_otus, by=list(Time=H.BL.Alpha$Time), FUN=mean)
H.BL.sd <- aggregate(H.BL.Alpha$observed_otus, by=list(Time=H.BL.Alpha$Time), FUN=sd)
H.FB.mean <- aggregate(H.FB.Alpha$observed_otus, by=list(Time=H.FB.Alpha$Time), FUN=mean)
H.FB.sd <- aggregate(H.FB.Alpha$observed_otus, by=list(Time=H.FB.Alpha$Time), FUN=sd)

# TimeChar <- as.character(CR$Time)[1:4]
# TimeChar <- log(H.CR.mean$Time[1:4],10)
TimeChar <- H.CR.mean$Time[1:4]

# Min and max values
# TimeChar <- factor(TimeChar, levels = c("20","60","360","720"))
pdf(paste0("ALPHAPLOTS_H/Alpha_diversity_Hakai_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0) ~ TimeChar[1:4]
     , main = "Richness through time"
     , ylim = c(min(MF.Alpha.morphonly.H$observed_otus), max(MF.Alpha.morphonly.H$observed_otus))
     , xlab = "Time"
     , ylab = "Richness (observed_otus)"
     ,xaxt = 'n'
)
axis(1
     , at = TimeChar[1:4]
     , labels = c("20m","1h","6h","12h")
     , las = 2)
points(TimeChar[1:4], H.CR.mean$x[1:4]
       , lty = 3
       # , col = 'darkgreen'
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
       # , col = "green"
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
       # , col = "yellow"
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
     , labels = c(paste0("20 min: "))
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
P.CR.sd <- aggregate(P.CR.Alpha$observed_otus, by=list(Time=P.CR.Alpha$Time), FUN=sd)
P.BL.mean <- aggregate(P.BL.Alpha$observed_otus, by=list(Time=P.BL.Alpha$Time), FUN=mean)
P.BL.sd <- aggregate(P.BL.Alpha$observed_otus, by=list(Time=P.BL.Alpha$Time), FUN=sd)
P.FB.mean <- aggregate(P.FB.Alpha$observed_otus, by=list(Time=P.FB.Alpha$Time), FUN=mean)
P.FB.sd <- aggregate(P.FB.Alpha$observed_otus, by=list(Time=P.FB.Alpha$Time), FUN=sd)

# TimeChar <- as.character(CR$Time)[1:4]
# TimeChar <- log(P.CR.mean$Time[1:4],10)
TimeChar <- P.CR.mean$Time[1:6]


# Min and max values
# TimeChar <- factor(TimeChar, levels = c("20","60","360","720"))
pdf(paste0("ALPHAPLOTS_P/Alpha_diversity_PM_",metric,".pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(c(0,0,0,0,0,0) ~ TimeChar[1:6]
     , main = "Richness through time"
     , ylim = c(min(MF.Alpha.morphonly.P$observed_otus), max(MF.Alpha.morphonly.P$observed_otus))
     , xlab = "Time"
     , ylab = "Richness (observed_otus)"
     , xaxt = 'n'
)
axis(1
     , at = c(-20,10,0,0,0,0) + TimeChar[1:6]
     , labels = c("20m","1h","3h","6h","12h","24h")
     , las = 2
     , cex.axis = 1)
points(TimeChar[1:6], P.CR.mean$x[1:6]
       , lty = 3
       # , col = 'darkgreen'
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
       # , col = "green"
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
       # , col = "yellow"
       , type = 'l')
arrows(x0 = TimeChar*1.01
       , x1 = TimeChar*1.01
       , y0 = c(P.FB.mean$x - P.FB.sd$x/2)[1:6]
       , y1 = c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
       , angle = 90
       , code = 3
       , length = 0.02)
text(TimeChar, 1.05*c(P.FB.mean$x + P.FB.sd$x/2)[1:6]
     , labels = c('','**','*','','*','')
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
     , labels = c(paste0(format(anova.P.morph.20.lm$`Pr(>F)`[1],digits = 2), ""))
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


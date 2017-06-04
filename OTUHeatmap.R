#!/bin/Rscript

library("optparse")
########################### OPT PARSE #################################
option_list = list(
  make_option(c("-i", "--inputFolder"), type="character",
              help="Folder with reduced OTU tables that are only core OTUs for each treatment type")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

inputFolder = opt$inputFolder

########################### FOR TESTING #################################
setwd('/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis')
inputFolder <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_ArtificialMacroalgae/1_analysis/ALLMORPH_cores'

########################### LOAD DATA #################################

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

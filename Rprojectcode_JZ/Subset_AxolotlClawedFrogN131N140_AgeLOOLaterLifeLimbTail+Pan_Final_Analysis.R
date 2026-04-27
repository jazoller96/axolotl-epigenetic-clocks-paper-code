rm(list=ls())
options(stringAsFactors=F)
library(tidyverse)
library(glmnet)
library(WGCNA)
setwd("~/Dropbox/MyResearchFiles/Horvath_mammalian_meth")
library(devtools)
library(MammalMethylClock)

load("AllNormalizedDataCBUAndInfo.RData")
## anAgeUpdated  datAllSamp  infoAllSamp  infoCBUAllSamp  probe_mappability_table  probe_amin_table_*
infoN131 <- read.csv("N131.ET0087.SalamanderMaxYun/SampleSheetAgeN131final.csv")
infoN140 <- read.csv("N140.2021-9212FrogChristofNiehrs/SampleSheetAgeN140final.csv")

###############################################################################
rm(datAllSamp,infoAllSamp,infoCBUAllSamp)
infoAllSamp1.csv=c('N131.ET0087.SalamanderMaxYun/SampleSheetAgeN131final.csv')
infoAllSamp2.csv=c('N140.2021-9212FrogChristofNiehrs/SampleSheetAgeN140final.csv')
datAllSamp_tp1.rdata=c('N131.ET0087.SalamanderMaxYun/NormalizedData/all_probes_sesame_normalized.Rdata')
datAllSamp_tp2.rdata=c('N140.2021-9212FrogChristofNiehrs/NormalizedData/all_probes_sesame_normalized.Rdata')

###############################################################################
### LOADING ALL DATA ###
###############################################################################
#load sample sheet data
infoAllSamp1=read.csv(infoAllSamp1.csv, as.is=T) %>%
  dplyr::select(Basename,SpeciesLatinName,OriginalOrderInBatch,Age,ConfidenceInAgeEstimate,
                CanBeUsedForAgingStudies,Tissue,Female,SpeciesCommonName,ExternalSampleID,Folder,
                Experiment,RegenExperimentGroup,AnimalID,AnimalName)
infoAllSamp2=read.csv(infoAllSamp2.csv, as.is=T) %>%
  dplyr::select(Basename,SpeciesLatinName,OriginalOrderInBatch,Age,ConfidenceInAgeEstimate,
                CanBeUsedForAgingStudies,Tissue,Female,SpeciesCommonName,ExternalSampleID,Folder) %>%
  dplyr::mutate(Experiment=NA,RegenExperimentGroup=NA,AnimalID=NA,AnimalName=NA)
infoAllSamp <- rbind(infoAllSamp1,infoAllSamp2)
rm(infoAllSamp1.csv,infoAllSamp2.csv,infoAllSamp1,infoAllSamp2)
if ("Female" %in% colnames(infoAllSamp)) {
  infoAllSamp$Female[which(is.na(infoAllSamp$Female))] <- "NA"
  infoAllSamp$Female <- factor(infoAllSamp$Female, levels=c(0,1,"NA"))
  levels(infoAllSamp$Female) <- c("Male","Female","NA")
}

#load reformatted DNA methylation data
datAllSamp1 <- transpose_dat(loadRData(datAllSamp_tp1.rdata) %>% as.data.frame(), "Basename")
datAllSamp2 <- transpose_dat(loadRData(datAllSamp_tp2.rdata) %>% as.data.frame(), "Basename")
datAllSamp <- rbind(datAllSamp1,datAllSamp2)
rm(datAllSamp_tp1.rdata,datAllSamp_tp2.rdata,datAllSamp1,datAllSamp2)

###############################################################################
### APPENDING anAge DATA TO SAMPLE SHEET ###
###############################################################################
anAgeUpdated.temp <- anAgeUpdated %>%
  dplyr::select(SpeciesLatinName,Female.maturity..days.,Male.maturity..days.,
                averagedMaturity.yrs,maxAgeCaesar,gestationYears)
infoAllSamp$idx <- 1:nrow(infoAllSamp)
infoAllSamp <- base::merge(infoAllSamp, anAgeUpdated.temp, by="SpeciesLatinName", all.x=T, sort=F)
infoAllSamp <- infoAllSamp[order(infoAllSamp$idx),]
infoAllSamp <- select(infoAllSamp, -idx)
rm(anAgeUpdated.temp)

###############################################################################
### REFINING DATA ###
###############################################################################
infoCBUAllSamp <- infoAllSamp %>% dplyr::filter(CanBeUsedForAgingStudies == "yes") %>%
  dplyr::filter(ConfidenceInAgeEstimate >= 90) %>%
  dplyr::filter(!is.na(Age)) %>%
  dplyr::filter(Basename %in% datAllSamp$Basename)

###############################################################################

### Sample filtering and partitioning
# anAgeUpdated <- dplyr::select(anAgeUpdated,SpeciesLatinName,Female.maturity..days.,
#                               Male.maturity..days.,averagedMaturity.yrs,maxAgeCaesar)
anAgeUpdatedaxolotlclawedfrogn131n140 <- anAgeUpdated %>% dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum","Xenopus laevis","Xenopus tropicalis"))
infoAllaxolotlclawedfrogn131n140 <- infoAllSamp %>%
  dplyr::filter(Folder %in% c("N131.ET0087.SalamanderMaxYun","N140.2021-9212FrogChristofNiehrs")) %>%
  dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum","Xenopus laevis","Xenopus tropicalis")) %>%
  dplyr::mutate(RelAge = fun_relative.trans(Age, maxAgeCaesar))
infoCBUaxolotlclawedfrogn131n140 <- infoCBUAllSamp %>%
  dplyr::filter(Folder %in% c("N131.ET0087.SalamanderMaxYun","N140.2021-9212FrogChristofNiehrs")) %>%
  dplyr::filter(!Tissue %in% c("Embryo")) %>%
  dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum","Xenopus laevis","Xenopus tropicalis")) %>%
  dplyr::mutate(RelAge = fun_relative.trans(Age, maxAgeCaesar)) %>%
  dplyr::filter(Experiment %in% c("AxolotlClock") | is.na(Experiment)) %>%
  dplyr::filter(Age > 4.0 | !SpeciesLatinName %in% c("Ambystoma mexicanum")) %>%
  dplyr::filter(Age > 2.0 | !SpeciesLatinName %in% c("Xenopus laevis","Xenopus tropicalis")) %>%
  dplyr::filter(Tissue %in% c("Limb","Tail") | !SpeciesLatinName %in% c("Ambystoma mexicanum"))
infoCBUaxolotlclawedfrogn131n140$SpeciesLatinName <- factor(infoCBUaxolotlclawedfrogn131n140$SpeciesLatinName)
infoCBUaxolotlclawedfrogn131n140$SpeciesCommonName <- factor(infoCBUaxolotlclawedfrogn131n140$SpeciesCommonName)
infoCBUaxolotlclawedfrogn131n140$Tissue <- factor(infoCBUaxolotlclawedfrogn131n140$Tissue)

infoNBUaxolotlclawedfrogn131n140 <- infoAllaxolotlclawedfrogn131n140 %>% dplyr::filter(!Basename %in% infoCBUaxolotlclawedfrogn131n140$Basename) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes") %>%
  dplyr::filter(!Tissue %in% c("Embryo")) %>%
  dplyr::filter(Age > 4.0 | !SpeciesLatinName %in% c("Ambystoma mexicanum")) %>%
  dplyr::filter(Age > 2.0 | !SpeciesLatinName %in% c("Xenopus laevis","Xenopus tropicalis")) %>%
  dplyr::filter(Tissue %in% c("Limb","Tail") | !SpeciesLatinName %in% c("Ambystoma mexicanum"))
infoNBUaxolotlclawedfrogn131n140$SpeciesLatinName <- factor(infoNBUaxolotlclawedfrogn131n140$SpeciesLatinName)
infoNBUaxolotlclawedfrogn131n140$SpeciesCommonName <- factor(infoNBUaxolotlclawedfrogn131n140$SpeciesCommonName)
infoNBUaxolotlclawedfrogn131n140$Tissue <- factor(infoNBUaxolotlclawedfrogn131n140$Tissue)



### Data refinement of CpG Sites, based on shared probe mappings
# probe_mappability_axolotlclawedfrogn131n140 <- na.omit(probe_mappability_table[,which(colnames(probe_mappability_table) %in% c("probeID","AmbystomaMexicanum","XenopusLaevis","XenopusTropicalis"))])
# #probe_mappability_table is a proper subset of datAllSamp
# datAllSamp_subCPGaxolotlclawedfrogn131n140 <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_mappability_axolotlclawedfrogn131n140$probeID))]
axolotlclawedfrogn131n140_SLNvecC <- c("CGid",colnames(probe_amin_table_amphibian)[unlist(sapply(c("Ambystoma_mexicanum","Xenopus_laevis","Xenopus_tropicalis"),grep,colnames(probe_amin_table_amphibian)))])
probe_amin_axolotlclawedfrogn131n140 <- probe_amin_table_amphibian[,which(colnames(probe_amin_table_amphibian) %in% axolotlclawedfrogn131n140_SLNvecC)]
probe_amin_axolotlclawedfrogn131n140$num_shared <- rowSums(sapply(as.data.frame(!is.na(probe_amin_axolotlclawedfrogn131n140)), as.numeric))-1
probe_amin_axolotlclawedfrogn131n140 <- dplyr::filter(probe_amin_axolotlclawedfrogn131n140, num_shared >= 1)
#probe_amin_table_amphibian is a proper subset of datAllSamp
datAllSamp_subCPGaxolotlclawedfrogn131n140 <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_amin_axolotlclawedfrogn131n140$CGid))]
probe_joseph_axolotln131 <- read_tsv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_probe_joseph_subCPGcombinationmiddlefilter.tsv')
probe_joseph_axolotln131 <- na.omit(probe_joseph_axolotln131)
#probe_joseph_table is a proper subset of datAllSamp
datAllSamp_subCPGaxonewtcombinationmiddlefilter <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_joseph_axolotln131$CGid))]



table(as.character(infoCBUaxolotlclawedfrogn131n140$Tissue))
latin2common_axolotlclawedfrogn131n140 <- unique(dplyr::select(infoCBUaxolotlclawedfrogn131n140,SpeciesLatinName,SpeciesCommonName))

###############################################################################

### Pre-Analysis, using full clocks
set.seed(1236)
yxs.list <- alignDatToInfo(infoCBUaxolotlclawedfrogn131n140,datAllSamp,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)

### Analysis, using full clocks
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticAge'
PREDVAR="DNAmAgeLOO"
RESVAR="AgeAccelLOO"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="SpeciesLatinName"
ALPHA=0.5
NFOLD=10
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticSqrtAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticSqrtAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticSqrtAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticSqrtAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_sqrt.trans,fun_inv=fun_sqrt.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLogAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLogAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLogAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLogAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log.trans,fun_inv=fun_log.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLog2Age.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLog2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLog2Age.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLog2Age'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log2.trans,fun_inv=fun_log2.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLLinAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLLinAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLLinAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLLinAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin.trans,fun_inv=fun_llin.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
# OUTVAR="Age"
# fun_VAR1="averagedMaturity.yrs"
# fun_VAR2="gestationYears"
# out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLLin2Age.RData'
# output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLLin2Age_PredictedValues.csv'
# out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLLin2Age.png'
# out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticLLin2Age'
# ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin2.trans,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
# ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
#                          "Basename", all=T, sort=F)
# ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
#                          "Basename", all=T, sort=F)
# write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="RelAge"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticRelativeAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticRelativeAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticRelativeAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_EpigeneticRelativeAge'
PREDVAR="DNAmRelAgeLOO"
RESVAR="RelAgeAccelLOO"
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

### Pre-Analysis, using subsetted clocks
set.seed(1236)
yxs.list <- alignDatToInfo(infoCBUaxolotlclawedfrogn131n140,datAllSamp_subCPGaxolotlclawedfrogn131n140,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)

### Analysis, using subsetted clocks
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticAge'
PREDVAR="DNAmAgeLOO"
RESVAR="AgeAccelLOO"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="SpeciesLatinName"
ALPHA=0.5
NFOLD=10
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticSqrtAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticSqrtAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticSqrtAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticSqrtAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_sqrt.trans,fun_inv=fun_sqrt.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLogAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLogAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLogAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLogAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log.trans,fun_inv=fun_log.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLog2Age.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLog2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLog2Age.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLog2Age'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log2.trans,fun_inv=fun_log2.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLLinAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLLinAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLLinAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLLinAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin.trans,fun_inv=fun_llin.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
# OUTVAR="Age"
# fun_VAR1="averagedMaturity.yrs"
# fun_VAR2="gestationYears"
# out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLLin2Age.RData'
# output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLLin2Age_PredictedValues.csv'
# out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLLin2Age.png'
# out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLLin2Age'
# ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin2.trans,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
# ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
#                          "Basename", all=T, sort=F)
# ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
#                          "Basename", all=T, sort=F)
# write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="RelAge"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticRelativeAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticRelativeAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticRelativeAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticRelativeAge'
PREDVAR="DNAmRelAgeLOO"
RESVAR="RelAgeAccelLOO"
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

### Pre-Analysis, using middle filter-subsetted clocks
set.seed(1236)
yxs.list <- alignDatToInfo(infoCBUaxolotlclawedfrogn131n140,datAllSamp_subCPGaxonewtcombinationmiddlefilter,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)

### Analysis, using middle filter-subsetted clocks
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticAge'
PREDVAR="DNAmAgeLOO"
RESVAR="AgeAccelLOO"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="SpeciesLatinName"
ALPHA=0.5
NFOLD=10
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticSqrtAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticSqrtAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticSqrtAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticSqrtAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_sqrt.trans,fun_inv=fun_sqrt.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLogAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLogAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLogAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLogAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log.trans,fun_inv=fun_log.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLog2Age.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLog2Age.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLog2Age'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log2.trans,fun_inv=fun_log2.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLLinAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLLinAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLLinAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLLinAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin.trans,fun_inv=fun_llin.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
# OUTVAR="Age"
# fun_VAR1="averagedMaturity.yrs"
# fun_VAR2="gestationYears"
# out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLLin2Age.RData'
# output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLLin2Age_PredictedValues.csv'
# out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLLin2Age.png'
# out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticLLin2Age'
# ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin2.trans,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
# ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
#                          "Basename", all=T, sort=F)
# ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
#                          "Basename", all=T, sort=F)
# write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="RelAge"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticRelativeAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticRelativeAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticRelativeAge.png'
out.png.title='Subset_AxolotlClawedFrogN131N140_LOOLaterLifeLimbTail+Pan_Final_subCPGaxonewtcombinationmiddlefilter_EpigeneticRelativeAge'
PREDVAR="DNAmRelAgeLOO"
RESVAR="RelAgeAccelLOO"
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotlclawedfrogn131n140, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)



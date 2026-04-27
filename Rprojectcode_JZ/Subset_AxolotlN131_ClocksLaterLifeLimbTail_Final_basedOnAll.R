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

###############################################################################
rm(datAllSamp,infoAllSamp,infoCBUAllSamp)
infoAllSamp.csv=c('N131.ET0087.SalamanderMaxYun/SampleSheetAgeN131final.csv')
datAllSamp_tp.rdata=c('N131.ET0087.SalamanderMaxYun/NormalizedData/all_probes_sesame_normalized.Rdata')

###############################################################################
### LOADING ALL DATA ###
###############################################################################
#load sample sheet data
infoAllSamp=read.csv(infoAllSamp.csv, as.is=T)
infoAllSamp <- infoAllSamp %>%
  dplyr::select(Basename,SpeciesLatinName,OriginalOrderInBatch,Age,ConfidenceInAgeEstimate,
                CanBeUsedForAgingStudies,Tissue,Female,SpeciesCommonName,ExternalSampleID,Folder,
                Experiment,RegenExperimentGroup,AnimalID,AnimalName)
rm(infoAllSamp.csv)
if ("Female" %in% colnames(infoAllSamp)) {
  infoAllSamp$Female[which(is.na(infoAllSamp$Female))] <- "NA"
  infoAllSamp$Female <- factor(infoAllSamp$Female, levels=c(0,1,"NA"))
  levels(infoAllSamp$Female) <- c("Male","Female","NA")
}

#load reformatted DNA methylation data
datAllSamp <- transpose_dat(loadRData(datAllSamp_tp.rdata) %>% as.data.frame(), "Basename")
rm(datAllSamp_tp.rdata)

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
anAgeUpdatedaxolotln131 <- anAgeUpdated %>% dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
infoAllaxolotln131 <- infoAllSamp %>%
  dplyr::filter(Folder %in% c("N131.ET0087.SalamanderMaxYun")) %>%
  dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
infoCBUaxolotln131 <- infoCBUAllSamp %>%
  dplyr::filter(Folder %in% c("N131.ET0087.SalamanderMaxYun")) %>%
  dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum")) %>%
  dplyr::filter(Experiment %in% c("AxolotlClock")) %>%
  dplyr::filter(Age > 4.0) %>%
  dplyr::filter(Tissue %in% c("Limb","Tail"))
infoCBUaxolotln131$SpeciesLatinName <- factor(infoCBUaxolotln131$SpeciesLatinName)
infoCBUaxolotln131$SpeciesCommonName <- factor(infoCBUaxolotln131$SpeciesCommonName)
infoCBUaxolotln131$Tissue <- factor(infoCBUaxolotln131$Tissue)

infoNBUaxolotln131 <- infoAllaxolotln131 %>% dplyr::filter(!Basename %in% infoCBUaxolotln131$Basename) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes") %>%
  dplyr::filter(Age > 4.0) %>%
  dplyr::filter(Tissue %in% c("Limb","Tail"))
infoNBUaxolotln131$SpeciesLatinName <- factor(infoNBUaxolotln131$SpeciesLatinName)
infoNBUaxolotln131$SpeciesCommonName <- factor(infoNBUaxolotln131$SpeciesCommonName)
infoNBUaxolotln131$Tissue <- factor(infoNBUaxolotln131$Tissue)



### Data refinement of CpG Sites, based on shared probe mappings
# probe_mappability_axolotln131 <- na.omit(probe_mappability_table[,which(colnames(probe_mappability_table) %in% c("probeID","AmbystomaMexicanum"))])
# #probe_mappability_table is a proper subset of datAllSamp
# datAllSamp_subCPGaxolotln131 <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_mappability_axolotln131$probeID))]
axolotln131_SLNvecC <- c("CGid",colnames(probe_amin_table_amphibian)[unlist(sapply(c("Ambystoma_mexicanum"),grep,colnames(probe_amin_table_amphibian)))])
probe_amin_axolotln131 <- probe_amin_table_amphibian[,which(colnames(probe_amin_table_amphibian) %in% axolotln131_SLNvecC)]
probe_amin_axolotln131 <- na.omit(probe_amin_axolotln131)
#probe_amin_table_amphibian is a proper subset of datAllSamp
datAllSamp_subCPGaxolotln131 <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_amin_axolotln131$CGid))]
probe_joseph_axolotln131 <- read_tsv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_probe_joseph_subCPGcombinationmiddlefilter.tsv')
probe_joseph_axolotln131 <- na.omit(probe_joseph_axolotln131)
#probe_joseph_table is a proper subset of datAllSamp
datAllSamp_subCPGcombinationmiddlefilter <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_joseph_axolotln131$CGid))]



table(as.character(infoCBUaxolotln131$Tissue))
latin2common_axolotln131 <- unique(dplyr::select(infoCBUaxolotln131,SpeciesLatinName,SpeciesCommonName))

###############################################################################

### Pre-Analysis, Training from All
set.seed(1236)
yxs.train.list <- alignDatToInfo(infoCBUaxolotln131,datAllSamp,"Basename","Basename")
yxs.other.list <- alignDatToInfo(infoNBUaxolotln131,datAllSamp,"Basename","Basename")
ys.train <- yxs.train.list[[1]]
xs.train <- yxs.train.list[[2]]
ys.other <- yxs.other.list[[1]]
xs.other <- yxs.other.list[[2]]
rm(yxs.train.list,yxs.other.list)

### Analysis, Training from All
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Square-Root Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticSqrtAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticSqrtAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticSqrtAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticSqrtAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_sqrt.trans,fun_inv=fun_sqrt.inv)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Logarithm Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLogAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLogAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLogAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLogAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_log.trans,fun_inv=fun_log.inv)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Logarithm Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLog2Age.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLog2Age.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLog2Age'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_log2.trans,fun_inv=fun_log2.inv)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)


###############################################################################

### Pre-Analysis, Training from All, using Special Transformation
#summary(infoCBUaxolotln131$averagedMaturity.yrs)
#summary(infoCBUaxolotln131$maxAgeCaesar)
infoCBUaxolotln1312 <- infoCBUaxolotln131 %>%
  dplyr::filter(!is.na(averagedMaturity.yrs)) # no change

set.seed(1236)
yxs.train.list <- alignDatToInfo(infoCBUaxolotln131,datAllSamp,"Basename","Basename")
yxs.other.list <- alignDatToInfo(infoNBUaxolotln131,datAllSamp,"Basename","Basename")
ys.train <- yxs.train.list[[1]]
xs.train <- yxs.train.list[[2]]
ys.other <- yxs.other.list[[1]]
xs.other <- yxs.other.list[[2]]
rm(yxs.train.list,yxs.other.list)

## Log+Linear Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLLinAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLLinAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLLinAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLLinAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_llin.trans,fun_inv=fun_llin.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLLin2Age.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLLin2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLLin2Age.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_basedOnAll_EpigeneticLLin2Age'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_llin2.trans,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

### Pre-Analysis, using subsetted Clock, Training from All
set.seed(1236)
yxs.train.list <- alignDatToInfo(infoCBUaxolotln131,datAllSamp_subCPGaxolotln131,"Basename","Basename")
yxs.other.list <- alignDatToInfo(infoNBUaxolotln131,datAllSamp_subCPGaxolotln131,"Basename","Basename")
ys.train <- yxs.train.list[[1]]
xs.train <- yxs.train.list[[2]]
ys.other <- yxs.other.list[[1]]
xs.other <- yxs.other.list[[2]]
rm(yxs.train.list,yxs.other.list)

### Analysis, using subsetted Clock, Training from All
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Square-Root Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticSqrtAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticSqrtAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticSqrtAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticSqrtAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_sqrt.trans,fun_inv=fun_sqrt.inv)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Logarithm Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLogAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLogAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLogAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLogAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_log.trans,fun_inv=fun_log.inv)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Logarithm Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_log2.trans,fun_inv=fun_log2.inv)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)


###############################################################################

### Pre-Analysis, using subsetted Clock, Training from All, using Special Transformation
#summary(infoCBUaxolotln131$averagedMaturity.yrs)
#summary(infoCBUaxolotln131$maxAgeCaesar)
infoCBUaxolotln1312 <- infoCBUaxolotln131 %>%
  dplyr::filter(!is.na(averagedMaturity.yrs)) # no change

set.seed(1236)
yxs.train.list <- alignDatToInfo(infoCBUaxolotln131,datAllSamp_subCPGaxolotln131,"Basename","Basename")
yxs.other.list <- alignDatToInfo(infoNBUaxolotln131,datAllSamp_subCPGaxolotln131,"Basename","Basename")
ys.train <- yxs.train.list[[1]]
xs.train <- yxs.train.list[[2]]
ys.other <- yxs.other.list[[1]]
xs.other <- yxs.other.list[[2]]
rm(yxs.train.list,yxs.other.list)

## Log+Linear Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLinAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLinAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLinAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLinAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_llin.trans,fun_inv=fun_llin.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_llin2.trans,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

### Pre-Analysis, using middle filter-subsetted Clock, Training from All
set.seed(1236)
yxs.train.list <- alignDatToInfo(infoCBUaxolotln131,datAllSamp_subCPGcombinationmiddlefilter,"Basename","Basename")
yxs.other.list <- alignDatToInfo(infoNBUaxolotln131,datAllSamp_subCPGcombinationmiddlefilter,"Basename","Basename")
ys.train <- yxs.train.list[[1]]
xs.train <- yxs.train.list[[2]]
ys.other <- yxs.other.list[[1]]
xs.other <- yxs.other.list[[2]]
rm(yxs.train.list,yxs.other.list)

### Analysis, using middle filter-subsetted Clock, Training from All
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Square-Root Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticSqrtAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticSqrtAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticSqrtAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticSqrtAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_sqrt.trans,fun_inv=fun_sqrt.inv)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Logarithm Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLogAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLogAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLogAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLogAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_log.trans,fun_inv=fun_log.inv)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Logarithm Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_log2.trans,fun_inv=fun_log2.inv)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)


###############################################################################

### Pre-Analysis, using middle filter-subsetted Clock, Training from All, using Special Transformation
#summary(infoCBUaxolotln131$averagedMaturity.yrs)
#summary(infoCBUaxolotln131$maxAgeCaesar)
infoCBUaxolotln1312 <- infoCBUaxolotln131 %>%
  dplyr::filter(!is.na(averagedMaturity.yrs)) # no change

set.seed(1236)
yxs.train.list <- alignDatToInfo(infoCBUaxolotln131,datAllSamp_subCPGcombinationmiddlefilter,"Basename","Basename")
yxs.other.list <- alignDatToInfo(infoNBUaxolotln131,datAllSamp_subCPGcombinationmiddlefilter,"Basename","Basename")
ys.train <- yxs.train.list[[1]]
xs.train <- yxs.train.list[[2]]
ys.other <- yxs.other.list[[1]]
xs.other <- yxs.other.list[[2]]
rm(yxs.train.list,yxs.other.list)

## Log+Linear Transformed
OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLinAge.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLinAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLinAge.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLinAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_llin.trans,fun_inv=fun_llin.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

OUTVAR="Age"
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age.csv'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age.png'
out.png.title='Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinOtherVAR="AgeAccelinOtherbasedOnAll"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveBuildClock(xs.train,ys.train,xs.other,ys.other,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinOtherVAR,ALPHA,NFOLD,fun_trans=fun_llin2.trans,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinOtherVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)



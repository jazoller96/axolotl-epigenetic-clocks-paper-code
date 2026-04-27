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
  dplyr::filter(Age <= 4.0)
infoCBUaxolotln131$SpeciesLatinName <- factor(infoCBUaxolotln131$SpeciesLatinName)
infoCBUaxolotln131$SpeciesCommonName <- factor(infoCBUaxolotln131$SpeciesCommonName)
infoCBUaxolotln131$Tissue <- factor(infoCBUaxolotln131$Tissue)

infoNBUaxolotln131 <- infoAllaxolotln131 %>% dplyr::filter(!Basename %in% infoCBUaxolotln131$Basename) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes") %>%
  dplyr::filter(Age <= 4.0)
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

### Pre-Analysis, using full clocks
set.seed(1236)
yxs.list <- alignDatToInfo(infoCBUaxolotln131,datAllSamp,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)

### Analysis, using full clocks
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticAge'
PREDVAR="DNAmAgeLOO"
RESVAR="AgeAccelLOO"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticSqrtAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticSqrtAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticSqrtAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticSqrtAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_sqrt.trans,fun_inv=fun_sqrt.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLogAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLogAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLogAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLogAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log.trans,fun_inv=fun_log.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLog2Age.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLog2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLog2Age.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLog2Age'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log2.trans,fun_inv=fun_log2.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLLinAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLLinAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLLinAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLLinAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin.trans,fun_inv=fun_llin.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLLin2Age.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLLin2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLLin2Age.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_EpigeneticLLin2Age'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin2.trans,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

### Pre-Analysis, using subsetted clocks
set.seed(1236)
yxs.list <- alignDatToInfo(infoCBUaxolotln131,datAllSamp_subCPGaxolotln131,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)

### Analysis, using subsetted clocks
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticAge'
PREDVAR="DNAmAgeLOO"
RESVAR="AgeAccelLOO"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticSqrtAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticSqrtAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticSqrtAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticSqrtAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_sqrt.trans,fun_inv=fun_sqrt.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLogAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLogAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLogAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLogAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log.trans,fun_inv=fun_log.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log2.trans,fun_inv=fun_log2.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLLinAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLLinAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLLinAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLLinAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin.trans,fun_inv=fun_llin.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin2Age.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin2Age.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin2Age'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin2.trans,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

### Pre-Analysis, using middle filter-subsetted clocks
set.seed(1236)
yxs.list <- alignDatToInfo(infoCBUaxolotln131,datAllSamp_subCPGcombinationmiddlefilter,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)

### Analysis, using middle filter-subsetted clocks
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticAge'
PREDVAR="DNAmAgeLOO"
RESVAR="AgeAccelLOO"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticSqrtAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticSqrtAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticSqrtAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticSqrtAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_sqrt.trans,fun_inv=fun_sqrt.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLogAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLogAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLogAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLogAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log.trans,fun_inv=fun_log.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log2.trans,fun_inv=fun_log2.inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLinAge.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLinAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLinAge.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLinAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin.trans,fun_inv=fun_llin.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
OUTVAR="Age"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
out.rdata='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin2Age.RData'
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin2Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin2Age.png'
out.png.title='Subset_AxolotlN131_LOOEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin2Age'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin2.trans,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)



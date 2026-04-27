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
  dplyr::filter(Experiment %in% c("AxolotlClock"))
infoCBUaxolotln131$SpeciesLatinName <- factor(infoCBUaxolotln131$SpeciesLatinName)
infoCBUaxolotln131$SpeciesCommonName <- factor(infoCBUaxolotln131$SpeciesCommonName)
infoCBUaxolotln131$Tissue <- factor(infoCBUaxolotln131$Tissue)

infoNBUaxolotln131 <- infoAllaxolotln131 %>% dplyr::filter(!Basename %in% infoCBUaxolotln131$Basename) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
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



table(as.character(infoCBUaxolotln131$Tissue))
latin2common_axolotln131 <- unique(dplyr::select(infoCBUaxolotln131,SpeciesLatinName,SpeciesCommonName))

###############################################################################
infoAllaxolotln131expmtlimbregenall <- infoAllaxolotln131 %>%
  dplyr::filter(Experiment %in% c("LimbRegeneration","LimbRegenerationSupplementary","LimbRegenerationSupplementary2",
                                  "LimbRegenerationExpanded"))
infoAllaxolotln131expmtlimbregenall$SpeciesLatinName <- factor(infoAllaxolotln131expmtlimbregenall$SpeciesLatinName)
infoAllaxolotln131expmtlimbregenall$SpeciesCommonName <- factor(infoAllaxolotln131expmtlimbregenall$SpeciesCommonName)
infoAllaxolotln131expmtlimbregenall$Tissue <- factor(infoAllaxolotln131expmtlimbregenall$Tissue)
infoAllaxolotln131expmtlimbregenall$RegenExperimentGroup <- factor(infoAllaxolotln131expmtlimbregenall$RegenExperimentGroup)
# infoACBUaxolotln131expmtlimbregenall <- infoAllaxolotln131expmtlimbregenall %>%
#   dplyr::filter(CanBeUsedForAgingStudies == "yes")
# infoACBUaxolotln131expmtlimbregenall$SpeciesLatinName <- factor(infoACBUaxolotln131expmtlimbregenall$SpeciesLatinName)
# infoACBUaxolotln131expmtlimbregenall$SpeciesCommonName <- factor(infoACBUaxolotln131expmtlimbregenall$SpeciesCommonName)
# infoACBUaxolotln131expmtlimbregenall$Tissue <- factor(infoACBUaxolotln131expmtlimbregenall$Tissue)
# infoACBUaxolotln131expmtlimbregenall$RegenExperimentGroup <- factor(infoACBUaxolotln131expmtlimbregenall$RegenExperimentGroup)
yxs.other.list <- alignDatToInfo(infoAllaxolotln131expmtlimbregenall,datAllSamp,"Basename","Basename")
ys.other <- yxs.other.list[[1]]
xs.other <- yxs.other.list[[2]]
rm(yxs.other.list)

###### Fitting Final Axolotl Early Life Clock to AxolotlN131ExpmtLimbRegen samples ######
OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="RegenExperimentGroup"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=10,square.axes=F)
ys.output <- base::merge(infoAllaxolotln131expmtlimbregenall, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131expmtlimbregenall, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

###### Fitting Final Axolotl Early Life LimbTail Clock to AxolotlN131ExpmtLimbRegen samples ######
OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen'
PREDVAR="DNAmAgebasedOnAllLimbTail"
RESVAR="AgeAccelbasedOnAllLimbTail"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="RegenExperimentGroup"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=10,square.axes=F)
ys.output <- base::merge(infoAllaxolotln131expmtlimbregenall, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131expmtlimbregenall, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

###### Fitting Final Axolotl Early Life Limb Clock to AxolotlN131ExpmtLimbRegen samples ######
OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen'
PREDVAR="DNAmAgebasedOnAllLimb"
RESVAR="AgeAccelbasedOnAllLimb"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="RegenExperimentGroup"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=10,square.axes=F)
ys.output <- base::merge(infoAllaxolotln131expmtlimbregenall, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131expmtlimbregenall, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

###### Fitting Final Axolotl Early Life Tail Clock to AxolotlN131ExpmtLimbRegen samples ######
OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen'
PREDVAR="DNAmAgebasedOnAllTail"
RESVAR="AgeAccelbasedOnAllTail"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="RegenExperimentGroup"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=10,square.axes=F)
ys.output <- base::merge(infoAllaxolotln131expmtlimbregenall, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln131expmtlimbregenall, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################



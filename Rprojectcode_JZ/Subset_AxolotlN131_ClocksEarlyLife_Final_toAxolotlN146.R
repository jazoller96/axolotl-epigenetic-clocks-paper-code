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
infoN146 <- read.csv("N146.ET0131AxolotlBohanGladyshev/SampleSheetAgeN146final.csv")

###############################################################################
rm(datAllSamp,infoAllSamp,infoCBUAllSamp)
infoAllSamp.csv=c('N146.ET0131AxolotlBohanGladyshev/SampleSheetAgeN146final.csv')
datAllSamp_tp.rdata=c('N146.ET0131AxolotlBohanGladyshev/NormalizedData/all_probes_sesame_normalized.Rdata')

###############################################################################
### LOADING ALL DATA ###
###############################################################################
#load sample sheet data
infoAllSamp=read.csv(infoAllSamp.csv, as.is=T)
infoAllSamp <- infoAllSamp %>%
  dplyr::select(Basename,SpeciesLatinName,OriginalOrderInBatch,Age,ConfidenceInAgeEstimate,
                CanBeUsedForAgingStudies,Tissue,Female,SpeciesCommonName,ExternalSampleID,Folder)
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
anAgeUpdatedaxolotln146 <- anAgeUpdated %>% dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
infoAllaxolotln146 <- infoAllSamp %>%
  dplyr::filter(Folder %in% c("N146.ET0131AxolotlBohanGladyshev"))
infoCBUaxolotln146 <- infoCBUAllSamp %>%
  dplyr::filter(Folder %in% c("N146.ET0131AxolotlBohanGladyshev"))
infoCBUaxolotln146$SpeciesLatinName <- factor(infoCBUaxolotln146$SpeciesLatinName)
infoCBUaxolotln146$SpeciesCommonName <- factor(infoCBUaxolotln146$SpeciesCommonName)
infoCBUaxolotln146$Tissue <- factor(infoCBUaxolotln146$Tissue)

infoNBUaxolotln146 <- infoAllaxolotln146 %>% dplyr::filter(!Basename %in% infoCBUaxolotln146$Basename) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
infoNBUaxolotln146$SpeciesLatinName <- factor(infoNBUaxolotln146$SpeciesLatinName)
infoNBUaxolotln146$SpeciesCommonName <- factor(infoNBUaxolotln146$SpeciesCommonName)
infoNBUaxolotln146$Tissue <- factor(infoNBUaxolotln146$Tissue)



### Data refinement of CpG Sites, based on shared probe mappings
# probe_mappability_axolotln146 <- na.omit(probe_mappability_table[,which(colnames(probe_mappability_table) %in% c("probeID","AmbystomaMexicanum"))])
# #probe_mappability_table is a proper subset of datAllSamp
# datAllSamp_subCPGaxolotln146 <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_mappability_axolotln146$probeID))]
axolotln146_SLNvecC <- c("CGid",colnames(probe_amin_table_amphibian)[unlist(sapply(c("Ambystoma_mexicanum"),grep,colnames(probe_amin_table_amphibian)))])
probe_amin_axolotln146 <- probe_amin_table_amphibian[,which(colnames(probe_amin_table_amphibian) %in% axolotln146_SLNvecC)]
probe_amin_axolotln146 <- na.omit(probe_amin_axolotln146)
#probe_amin_table_amphibian is a proper subset of datAllSamp
datAllSamp_subCPGaxolotln146 <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_amin_axolotln146$CGid))]



table(as.character(infoCBUaxolotln146$Tissue))
latin2common_axolotln146 <- unique(dplyr::select(infoCBUaxolotln146,SpeciesLatinName,SpeciesCommonName))

###############################################################################
infoACBUaxolotln146 <- infoAllaxolotln146 %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
infoACBUaxolotln146$SpeciesLatinName <- factor(infoACBUaxolotln146$SpeciesLatinName)
infoACBUaxolotln146$SpeciesCommonName <- factor(infoACBUaxolotln146$SpeciesCommonName)
infoACBUaxolotln146$Tissue <- factor(infoACBUaxolotln146$Tissue)
yxs.other.list <- alignDatToInfo(infoACBUaxolotln146,datAllSamp,"Basename","Basename")
ys.other <- yxs.other.list[[1]]
xs.other <- yxs.other.list[[2]]
rm(yxs.other.list)

###### Fitting Final Axolotl Early Life Clock to AxolotlN146 Early Life samples ######
OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticAge.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying Log2 as well
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying LLin2 as well
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)

OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying Log2 as well
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying LLin2 as well
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLife_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)

###############################################################################

###### Fitting Final Axolotl Early Life LimbTail Clock to AxolotlN146 Early Life samples ######
OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146'
PREDVAR="DNAmAgebasedOnAllLimbTail"
RESVAR="AgeAccelbasedOnAllLimbTail"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying Log2 as well
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying LLin2 as well
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)

OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146'
PREDVAR="DNAmAgebasedOnAllLimbTail"
RESVAR="AgeAccelbasedOnAllLimbTail"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying Log2 as well
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying LLin2 as well
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)

###############################################################################

###### Fitting Final Axolotl Early Life Limb Clock to AxolotlN146 Early Life samples ######
OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticAge.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146'
PREDVAR="DNAmAgebasedOnAllLimb"
RESVAR="AgeAccelbasedOnAllLimb"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying Log2 as well
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying LLin2 as well
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)

OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146'
PREDVAR="DNAmAgebasedOnAllLimb"
RESVAR="AgeAccelbasedOnAllLimb"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying Log2 as well
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying LLin2 as well
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)

###############################################################################

###### Fitting Final Axolotl Early Life Tail Clock to AxolotlN146 Early Life samples ######
OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticAge.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146'
PREDVAR="DNAmAgebasedOnAllTail"
RESVAR="AgeAccelbasedOnAllTail"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying Log2 as well
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying LLin2 as well
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)

OUTVAR="Age"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146'
PREDVAR="DNAmAgebasedOnAllTail"
RESVAR="AgeAccelbasedOnAllTail"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticAge_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying Log2 as well
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_log2.inv,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)
## Trying LLin2 as well
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="gestationYears"
in.valbeta <- read.csv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age.csv')
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146'
ys.output <- saveApplyClock(in.valbeta,xs.other,ys.other,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_llin2.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,oma.right=3,square.axes=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllaxolotln146, ys.output[, c("Basename", PREDVAR, RESVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife.png'
out.png.title='Subset_AxolotlN131_ClockEarlyLifeTail_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLLin2Age_toAxolotlN146EarlyLife'
ys.output <- ys.output %>% dplyr::filter(Age <= 4.0)
ys.output$Tissue <- factor(ys.output$Tissue)
saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=3)

###############################################################################



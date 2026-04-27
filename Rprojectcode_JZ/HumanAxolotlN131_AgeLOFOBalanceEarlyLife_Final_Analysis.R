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
infoAllSamp_human <- infoAllSamp %>% dplyr::filter(SpeciesLatinName %in% c("Homo sapiens"))
datAllSamp_human <- datAllSamp %>% dplyr::filter(Basename %in% infoAllSamp_human$Basename)
rm(datAllSamp,infoAllSamp,infoCBUAllSamp)
infoAllSamp.csv=c('N131.ET0087.SalamanderMaxYun/SampleSheetAgeN131final.csv')
datAllSamp_tp.rdata=c('N131.ET0087.SalamanderMaxYun/NormalizedData/all_probes_sesame_normalized.Rdata')

###############################################################################
### LOADING ALL DATA ###
###############################################################################
#load sample sheet data
infoAllSamp=read.csv(infoAllSamp.csv, as.is=T) %>%
  dplyr::select(Basename,SpeciesLatinName,OriginalOrderInBatch,Age,ConfidenceInAgeEstimate,
                CanBeUsedForAgingStudies,Tissue,Female,SpeciesCommonName,ExternalSampleID,Folder,
                Experiment,RegenExperimentGroup,AnimalID,AnimalName)
infoAllSamp_human <- infoAllSamp_human %>%
  dplyr::select(Basename,SpeciesLatinName,OriginalOrderInBatch,Age,ConfidenceInAgeEstimate,
                CanBeUsedForAgingStudies,Tissue,Female,SpeciesCommonName,ExternalSampleID,Folder) %>%
  dplyr::mutate(Experiment=NA,RegenExperimentGroup=NA,AnimalID=NA,AnimalName=NA)
infoAllSamp <- rbind(infoAllSamp,infoAllSamp_human)
rm(infoAllSamp.csv)
rm(infoAllSamp_human)
if ("Female" %in% colnames(infoAllSamp)) {
  infoAllSamp$Female[which(is.na(infoAllSamp$Female))] <- "NA"
  infoAllSamp$Female <- factor(infoAllSamp$Female, levels=c(0,1,"NA"))
  levels(infoAllSamp$Female) <- c("Male","Female","NA")
}

#load reformatted DNA methylation data
datAllSamp <- transpose_dat(loadRData(datAllSamp_tp.rdata) %>% as.data.frame(), "Basename")
datAllSamp <- rbind(datAllSamp,datAllSamp_human)
rm(datAllSamp_tp.rdata)
rm(datAllSamp_human)

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
anAgeUpdatedhumanaxolotln131 <- anAgeUpdated %>% dplyr::filter(SpeciesLatinName %in% c("Homo sapiens","Ambystoma mexicanum"))
infoAllhumanaxolotln131 <- infoAllSamp %>%
  dplyr::filter(SpeciesLatinName %in% c("Homo sapiens","Ambystoma mexicanum")) %>%
  dplyr::mutate(RelAge = fun_relative.trans(Age, maxAgeCaesar))
infoCBUhumanaxolotln131 <- infoCBUAllSamp %>%
  dplyr::filter(SpeciesLatinName %in% c("Homo sapiens","Ambystoma mexicanum")) %>%
  dplyr::filter(!Folder %in% c("N33.2019-9152HumanBloodAnil")) %>%
  dplyr::mutate(RelAge = fun_relative.trans(Age, maxAgeCaesar)) %>%
  dplyr::filter(Experiment %in% c("AxolotlClock") | is.na(Experiment)) %>%
  dplyr::filter(Age <= 4.0 | !SpeciesLatinName %in% c("Ambystoma mexicanum")) %>%
  dplyr::filter(Age <= 25.0 | !SpeciesLatinName %in% c("Homo sapiens"))
infoCBUhumanaxolotln131$SpeciesLatinName <- factor(infoCBUhumanaxolotln131$SpeciesLatinName)
infoCBUhumanaxolotln131$SpeciesCommonName <- factor(infoCBUhumanaxolotln131$SpeciesCommonName)
infoCBUhumanaxolotln131$Tissue <- factor(infoCBUhumanaxolotln131$Tissue)

infoNBUhumanaxolotln131 <- infoAllhumanaxolotln131 %>% dplyr::filter(!Basename %in% infoCBUhumanaxolotln131$Basename) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes") %>%
  dplyr::filter(Age <= 4.0 | !SpeciesLatinName %in% c("Ambystoma mexicanum")) %>%
  dplyr::filter(Age <= 25.0 | !SpeciesLatinName %in% c("Homo sapiens"))
infoNBUhumanaxolotln131$SpeciesLatinName <- factor(infoNBUhumanaxolotln131$SpeciesLatinName)
infoNBUhumanaxolotln131$SpeciesCommonName <- factor(infoNBUhumanaxolotln131$SpeciesCommonName)
infoNBUhumanaxolotln131$Tissue <- factor(infoNBUhumanaxolotln131$Tissue)

#summary(infoCBUhumanaxolotln131$averagedMaturity.yrs)
#summary(infoCBUhumanaxolotln131$maxAgeCaesar)
infoCBUhumanaxolotln1312 <- infoCBUhumanaxolotln131 %>%
  dplyr::filter(!is.na(averagedMaturity.yrs)) # TODO check



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



table(as.character(infoCBUhumanaxolotln131$Tissue))
table(as.character(infoCBUhumanaxolotln131$SpeciesLatinName))
latin2common_humanaxolotln131 <- unique(dplyr::select(infoCBUhumanaxolotln131,SpeciesLatinName,SpeciesCommonName))

###############################################################################

### Pre-Analysis
set.seed(123456)
yxs.list <- alignDatToInfo(infoCBUhumanaxolotln131,datAllSamp,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)
#Defining balanced folds
SPECVAR="FoldNumber"
ys[,SPECVAR] <- NA
for (spec in levels(ys[,"SpeciesLatinName"])) {
  idx <- which(ys[,"SpeciesLatinName"] %in% c(spec))
  ys[idx,SPECVAR] <- sample.int(length(idx))%%10+1
}
rm(spec, idx)
ys[,SPECVAR] <- factor(ys[,SPECVAR])

## Log+Linear Transformed
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticLLin3Age.RData'
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticLLin3Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticLLin3Age.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticLLin3Age'
PREDVAR="DNAmAgeLOFO10Balance"
RESVAR="AgeAccelLOFO10Balance"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="SpeciesLatinName"
ALPHA=0.5
NFOLD=10
loglambda.seq <- seq(-6,-1,length.out=100)
ys.output <- saveLOSOEstimation(xs,ys,OUTVAR,SPECVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin3.trans,fun_inv=fun_llin3.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,loglambda.seq=loglambda.seq)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticLLin3Age_AllAxolotlN131PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticLLin3Age_AllAxolotlN131.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticLLin3Age_AllAxolotlN131'
COLVAR="Tissue"
temp <- ys.output %>% dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
write.table(temp,output.csv,sep=',',row.names=F,quote=F)
temp <- temp[!is.na(temp[,PREDVAR]),]
temp <- temp %>% dplyr::mutate(Tissue = factor(Tissue))
saveValidationPlot(temp,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=7)
rm(temp)

## Relative Age
OUTVAR="RelAge"
out.rdata='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticRelativeAge.RData'
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticRelativeAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticRelativeAge.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticRelativeAge'
PREDVAR="DNAmRelAgeLOFO10Balance"
RESVAR="RelAgeAccelLOFO10Balance"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="SpeciesLatinName"
ALPHA=0.5
NFOLD=10
loglambda.seq <- seq(-8,-3,length.out=100)
ys.output <- saveLOSOEstimation(xs,ys,OUTVAR,SPECVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR,loglambda.seq=loglambda.seq)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticRelativeAge_AllAxolotlN131PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticRelativeAge_AllAxolotlN131.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_EpigeneticRelativeAge_AllAxolotlN131'
COLVAR="Tissue"
temp <- ys.output %>% dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
write.table(temp,output.csv,sep=',',row.names=F,quote=F)
temp <- temp[!is.na(temp[,PREDVAR]),]
temp <- temp %>% dplyr::mutate(Tissue = factor(Tissue))
saveValidationPlot(temp,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=7)
rm(temp)

###############################################################################

### Pre-Analysis, using Axolotl-subsetted clocks
set.seed(123456)
yxs.list <- alignDatToInfo(infoCBUhumanaxolotln131,datAllSamp_subCPGaxolotln131,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)
#Defining balanced folds
SPECVAR="FoldNumber"
ys[,SPECVAR] <- NA
for (spec in levels(ys[,"SpeciesLatinName"])) {
  idx <- which(ys[,"SpeciesLatinName"] %in% c(spec))
  ys[idx,SPECVAR] <- sample.int(length(idx))%%10+1
}
rm(spec, idx)
ys[,SPECVAR] <- factor(ys[,SPECVAR])

## Log+Linear Transformed
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age.RData'
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age'
PREDVAR="DNAmAgeLOFO10Balance"
RESVAR="AgeAccelLOFO10Balance"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="SpeciesLatinName"
ALPHA=0.5
NFOLD=10
loglambda.seq <- seq(-6,-1,length.out=100)
ys.output <- saveLOSOEstimation(xs,ys,OUTVAR,SPECVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin3.trans,fun_inv=fun_llin3.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,loglambda.seq=loglambda.seq)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age_AllAxolotlN131PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age_AllAxolotlN131.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age_AllAxolotlN131'
COLVAR="Tissue"
temp <- ys.output %>% dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
write.table(temp,output.csv,sep=',',row.names=F,quote=F)
temp <- temp[!is.na(temp[,PREDVAR]),]
temp <- temp %>% dplyr::mutate(Tissue = factor(Tissue))
saveValidationPlot(temp,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=7)
rm(temp)

## Relative Age
OUTVAR="RelAge"
out.rdata='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticRelativeAge.RData'
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticRelativeAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticRelativeAge.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticRelativeAge'
PREDVAR="DNAmRelAgeLOFO10Balance"
RESVAR="RelAgeAccelLOFO10Balance"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="SpeciesLatinName"
ALPHA=0.5
NFOLD=10
loglambda.seq <- seq(-8,-3,length.out=100)
ys.output <- saveLOSOEstimation(xs,ys,OUTVAR,SPECVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR,loglambda.seq=loglambda.seq)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticRelativeAge_AllAxolotlN131PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticRelativeAge_AllAxolotlN131.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticRelativeAge_AllAxolotlN131'
COLVAR="Tissue"
temp <- ys.output %>% dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
write.table(temp,output.csv,sep=',',row.names=F,quote=F)
temp <- temp[!is.na(temp[,PREDVAR]),]
temp <- temp %>% dplyr::mutate(Tissue = factor(Tissue))
saveValidationPlot(temp,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=7)
rm(temp)

###############################################################################

### Pre-Analysis, using Axolotl-middle filter-subsetted clocks
set.seed(123456)
yxs.list <- alignDatToInfo(infoCBUhumanaxolotln131,datAllSamp_subCPGcombinationmiddlefilter,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)
#Defining balanced folds
SPECVAR="FoldNumber"
ys[,SPECVAR] <- NA
for (spec in levels(ys[,"SpeciesLatinName"])) {
  idx <- which(ys[,"SpeciesLatinName"] %in% c(spec))
  ys[idx,SPECVAR] <- sample.int(length(idx))%%10+1
}
rm(spec, idx)
ys[,SPECVAR] <- factor(ys[,SPECVAR])

## Log+Linear Transformed
OUTVAR="Age"
out.rdata='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin3Age.RData'
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin3Age_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin3Age.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin3Age'
PREDVAR="DNAmAgeLOFO10Balance"
RESVAR="AgeAccelLOFO10Balance"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="SpeciesLatinName"
ALPHA=0.5
NFOLD=10
loglambda.seq <- seq(-6,-1,length.out=100)
ys.output <- saveLOSOEstimation(xs,ys,OUTVAR,SPECVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin3.trans,fun_inv=fun_llin3.inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR,loglambda.seq=loglambda.seq)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin3Age_AllAxolotlN131PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin3Age_AllAxolotlN131.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticLLin3Age_AllAxolotlN131'
COLVAR="Tissue"
temp <- ys.output %>% dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
write.table(temp,output.csv,sep=',',row.names=F,quote=F)
temp <- temp[!is.na(temp[,PREDVAR]),]
temp <- temp %>% dplyr::mutate(Tissue = factor(Tissue))
saveValidationPlot(temp,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=7)
rm(temp)

## Relative Age
OUTVAR="RelAge"
out.rdata='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticRelativeAge.RData'
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticRelativeAge_PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticRelativeAge.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticRelativeAge'
PREDVAR="DNAmRelAgeLOFO10Balance"
RESVAR="RelAgeAccelLOFO10Balance"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="SpeciesLatinName"
ALPHA=0.5
NFOLD=10
loglambda.seq <- seq(-8,-3,length.out=100)
ys.output <- saveLOSOEstimation(xs,ys,OUTVAR,SPECVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR,loglambda.seq=loglambda.seq)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllhumanaxolotln131, ys.output[, c("Basename", PREDVAR, RESVAR, "count_probes_selected", "log_lambda_hat")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
output.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticRelativeAge_AllAxolotlN131PredictedValues.csv'
out.png='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticRelativeAge_AllAxolotlN131.png'
out.png.title='HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGcombinationmiddlefilter_EpigeneticRelativeAge_AllAxolotlN131'
COLVAR="Tissue"
temp <- ys.output %>% dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
write.table(temp,output.csv,sep=',',row.names=F,quote=F)
temp <- temp[!is.na(temp[,PREDVAR]),]
temp <- temp %>% dplyr::mutate(Tissue = factor(Tissue))
saveValidationPlot(temp,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str=paste0(out.png.title,'\n'),width=5,height=6,oma.right=7)
rm(temp)



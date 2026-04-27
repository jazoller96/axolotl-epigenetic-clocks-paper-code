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
                Experiment,PassageNumber,Comparison4.DecitabineConcentration,Comparison5.DecitabineDurationTreatment.Days,AnimalID,AnimalName)
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
probe_joseph_axolotln131 <- read_tsv('SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_probe_joseph_subCPGcombinationmiddlefilter.tsv')
probe_joseph_axolotln131 <- na.omit(probe_joseph_axolotln131)
#probe_joseph_table is a proper subset of datAllSamp
datAllSamp_subCPGcombinationmiddlefilter <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_joseph_axolotln131$CGid))]



table(as.character(infoCBUaxolotln131$Tissue))
latin2common_axolotln131 <- unique(dplyr::select(infoCBUaxolotln131,SpeciesLatinName,SpeciesCommonName))

###############################################################################
infoAllaxolotln131expmtdemethylation <- infoAllaxolotln131 %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock","DemethylationTreatmentCells_Revision","DemethylationTreatmentCells","DemethylationTreatmentLarva"))
infoAllaxolotln131expmtdemethylation$SpeciesLatinName <- factor(infoAllaxolotln131expmtdemethylation$SpeciesLatinName)
infoAllaxolotln131expmtdemethylation$SpeciesCommonName <- factor(infoAllaxolotln131expmtdemethylation$SpeciesCommonName)
infoAllaxolotln131expmtdemethylation$Tissue <- factor(infoAllaxolotln131expmtdemethylation$Tissue)
infoAllaxolotln131expmtdemethylation$Experiment <- factor(infoAllaxolotln131expmtdemethylation$Experiment)
# infoACBUaxolotln131expmtdemethylation <- infoAllaxolotln131expmtdemethylation %>%
#   dplyr::filter(CanBeUsedForAgingStudies == "yes")
# infoACBUaxolotln131expmtdemethylation$SpeciesLatinName <- factor(infoACBUaxolotln131expmtdemethylation$SpeciesLatinName)
# infoACBUaxolotln131expmtdemethylation$SpeciesCommonName <- factor(infoACBUaxolotln131expmtdemethylation$SpeciesCommonName)
# infoACBUaxolotln131expmtdemethylation$Tissue <- factor(infoACBUaxolotln131expmtdemethylation$Tissue)
# infoACBUaxolotln131expmtdemethylation$Experiment <- factor(infoACBUaxolotln131expmtdemethylation$Experiment)
yxs.other.list <- alignDatToInfo(infoAllaxolotln131expmtdemethylation,datAllSamp_subCPGcombinationmiddlefilter,"Basename","Basename")
ys.other <- yxs.other.list[[1]]
xs.other <- yxs.other.list[[2]]
rm(yxs.other.list)
yxs.reference.list <- alignDatToInfo(infoAllaxolotln131,datAllSamp_subCPGcombinationmiddlefilter,"Basename","Basename")
ys.reference <- yxs.reference.list[[1]]
xs.reference <- yxs.reference.list[[2]]
rm(yxs.reference.list)

cpg_names_high <- names(which(colMeans(xs.reference, na.rm=T) > 0.50))
cpg_names_low <- names(which(colMeans(xs.reference, na.rm=T) < 0.50))

###### Calculating Overall Mean Methylation Values for AxolotlN131ExpmtDemethylation samples ######
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylation_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
PREDVAR="MeanMethylation"
ys.output <- infoAllaxolotln131expmtdemethylation
ys.output[,PREDVAR]=rowMeans(xs.other)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

###### Calculating High Mean Methylation Values for AxolotlN131ExpmtDemethylation samples ######
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationHigh_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
PREDVAR="MeanMethylationHigh"
ys.output <- infoAllaxolotln131expmtdemethylation
ys.output[,PREDVAR]=rowMeans(xs.other[,cpg_names_high])
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################

###### Calculating Low Mean Methylation Values for AxolotlN131ExpmtDemethylation samples ######
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationLow_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
PREDVAR="MeanMethylationLow"
ys.output <- infoAllaxolotln131expmtdemethylation
ys.output[,PREDVAR]=rowMeans(xs.other[,cpg_names_low])
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

###############################################################################



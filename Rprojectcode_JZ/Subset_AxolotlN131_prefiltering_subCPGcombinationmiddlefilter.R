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

infoAllnewtn131 <- infoAllSamp %>%
  dplyr::filter(Folder %in% c("N131.ET0087.SalamanderMaxYun")) %>%
  dplyr::filter(SpeciesLatinName %in% c("Pleurodeles waltl"))



### Data refinement of CpG Sites, based on shared probe mappings
# probe_mappability_axolotln131 <- na.omit(probe_mappability_table[,which(colnames(probe_mappability_table) %in% c("probeID","AmbystomaMexicanum"))])
# #probe_mappability_table is a proper subset of datAllSamp
# datAllSamp_subCPGaxolotln131 <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_mappability_axolotln131$probeID))]
axolotln131_SLNvecC <- c("CGid",colnames(probe_amin_table_amphibian)[unlist(sapply(c("Ambystoma_mexicanum"),grep,colnames(probe_amin_table_amphibian)))])
probe_amin_axolotln131 <- probe_amin_table_amphibian[,which(colnames(probe_amin_table_amphibian) %in% axolotln131_SLNvecC)]
probe_amin_axolotln131 <- na.omit(probe_amin_axolotln131)

###############################################################################

### Data- and Biologically-driven refinement of CpG Sites, based on combination middle filter + probe mappability filter
infoACBUtrainaxolotln131 <- infoAllaxolotln131 %>%
  dplyr::filter(Experiment %in% c("AxolotlClock")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
infoACBUtrainaxolotln131$SpeciesLatinName <- factor(infoACBUtrainaxolotln131$SpeciesLatinName)
infoACBUtrainaxolotln131$SpeciesCommonName <- factor(infoACBUtrainaxolotln131$SpeciesCommonName)
infoACBUtrainaxolotln131$Tissue <- factor(infoACBUtrainaxolotln131$Tissue)
yxs.list <- alignDatToInfo(infoACBUtrainaxolotln131,datAllSamp,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)

infoACBUtrainnewtn131 <- infoAllnewtn131 %>%
  dplyr::filter(Experiment %in% c("Newt")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
infoACBUtrainnewtn131$SpeciesLatinName <- factor(infoACBUtrainnewtn131$SpeciesLatinName)
infoACBUtrainnewtn131$SpeciesCommonName <- factor(infoACBUtrainnewtn131$SpeciesCommonName)
infoACBUtrainnewtn131$Tissue <- factor(infoACBUtrainnewtn131$Tissue)
yxs2.list <- alignDatToInfo(infoACBUtrainnewtn131,datAllSamp,"Basename","Basename")
ys2 <- yxs2.list[[1]]
xs2 <- yxs2.list[[2]]
rm(yxs2.list)

## Middle filter (in Axolotl data)
cpg_names_middlefilter.axolotl <- names(which(abs(colMeans(xs, na.rm=T) - 0.5) > 0.03))
## Middle filter (in Newt data)
cpg_names_middlefilter.newt <- names(which(abs(colMeans(xs2, na.rm=T) - 0.5) > 0.03))
## Probe mappability filter (to Axolotl genome)
cpg_names_probe.axolotl <- probe_amin_axolotln131$CGid
## Combination filters
cpg_names_combinationmiddlefilter <- intersect(cpg_names_middlefilter.axolotl,cpg_names_middlefilter.newt)
cpg_names_combinationmiddlefilter <- intersect(cpg_names_combinationmiddlefilter,cpg_names_probe.axolotl)
probe_joseph_axolotln131_combinationmiddlefilter <- base::merge(data.frame(CGid=probe_amin_table$CGid,stringsAsFactors=F),
                                                                data.frame(CGid=cpg_names_combinationmiddlefilter,"yes",stringsAsFactors=F),
                                                                "CGid", all=T, sort=F)
probe_joseph_axolotln131_combinationmiddlefilter <- base::merge(data.frame(CGid=probe_amin_table$CGid,stringsAsFactors=F),
                                                                probe_joseph_axolotln131_combinationmiddlefilter,
                                                                "CGid", all=T, sort=F)
colnames(probe_joseph_axolotln131_combinationmiddlefilter)=c("CGid","AmbystomaMexicanum")

probe_joseph_axolotln131_combinationmiddlefilter.tsv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_probe_joseph_subCPGcombinationmiddlefilter.tsv'
write.table(probe_joseph_axolotln131_combinationmiddlefilter,probe_joseph_axolotln131_combinationmiddlefilter.tsv,sep='\t',row.names=F,quote=F)



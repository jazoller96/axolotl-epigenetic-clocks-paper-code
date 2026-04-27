rm(list=ls())
options(stringAsFactors=F)
library(tidyverse)
library(glmnet)
library(WGCNA)
setwd("~/Dropbox/MyResearchFiles/Horvath_mammalian_meth")
library(devtools)
library(MammalMethylClock)

### Plotting separately by tissue
OUTVAR="Age"
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYears_Final_EpigeneticLog2Age_PredictedValues.csv'
PREDVAR="DNAmAgeLOO"
PANELVAR="Tissue"
input.info_pred <- dplyr::filter(read.csv(output.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::select(Basename, Age, DNAmAgeLOO, Tissue) %>%
  dplyr::mutate(Tissue = factor(Tissue))
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYears_Final_EpigeneticLog2Age_PANEL.png'
out.png.title='DNAmAgeLOO for Axolotl, by Tissue'
saveValidationPanelPlot(input.info_pred,OUTVAR,PREDVAR,PANELVAR,out.png,paste0(out.png.title,'\n'),mfrow=c(3,3),width=13,height=14)
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYears_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
input.info_pred <- dplyr::filter(read.csv(output.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::select(Basename, Age, DNAmAgeLOO, Tissue) %>%
  dplyr::mutate(Tissue = factor(Tissue))
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYears_Final_subCPGaxolotln131_EpigeneticLog2Age_PANEL.png'
out.png.title='DNAmAgeLOO for Axolotl, by Tissue'
saveValidationPanelPlot(input.info_pred,OUTVAR,PREDVAR,PANELVAR,out.png,paste0(out.png.title,'\n'),mfrow=c(3,3),width=13,height=14)
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYears_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
input.info_pred <- dplyr::filter(read.csv(output.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::select(Basename, Age, DNAmAgeLOO, Tissue) %>%
  dplyr::mutate(Tissue = factor(Tissue))
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYears_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PANEL.png'
out.png.title='DNAmAgeLOO for Axolotl, by Tissue'
saveValidationPanelPlot(input.info_pred,OUTVAR,PREDVAR,PANELVAR,out.png,paste0(out.png.title,'\n'),mfrow=c(3,3),width=13,height=14)



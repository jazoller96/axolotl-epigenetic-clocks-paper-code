rm(list=ls())
options(stringAsFactors=F)
library(tidyverse)
library(glmnet)
library(WGCNA)
setwd("~/Dropbox/MyResearchFiles/Horvath_mammalian_meth")
library(devtools)
library(MammalMethylClock)

### (-1) Saving Data Deposition Sample Sheet and Normalized Data (for Publication Journal)
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/info_DataDeposition_SampleSheet.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/info_DataDeposition_SampleSheet.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/dat_DataDeposition_NormalizedData.csv'
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/dat_DataDeposition_NormalizedData.csv'
file.copy(in.csv, out.csv, overwrite=T)
rm(in.csv, out.csv)

### (0) Saving Middle Filter Probe Mapping File
in.tsv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_probe_amin_subCPGaxolotl.tsv'
in.valbeta <- read_tsv(in.tsv)
out.tsv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/probe_amin_axolotln131_probemapping.tsv'
write.table(in.valbeta,out.tsv,sep='\t',row.names=F,quote=F)
rm(in.tsv, out.tsv)
in.tsv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_probe_joseph_subCPGcombinationmiddlefilter.tsv'
in.valbeta <- read_tsv(in.tsv)
out.tsv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/probe_joseph_axolotln131_combinationmiddlefilter.tsv'
write.table(in.valbeta,out.tsv,sep='\t',row.names=F,quote=F)
rm(in.tsv, out.tsv)

### (1) Combining and Saving Tables of Coefficients
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_Clock_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalLimbTailJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalLimbJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalTailJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalLaterLifeJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalLaterLifeLimbTailJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalLaterLifeLimbJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalLaterLifeTailJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalEarlyLifeJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalEarlyLifeLimbTailJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalEarlyLifeLimbJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datCoefN131FinalEarlyLifeTailJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)

### (2) Combining and Saving Tables of Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_Clock_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOO_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalLimbTailJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalLimbJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalTailJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)

## Later Life Clocks
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalLaterLifeJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalLaterLifeLimbTailJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalLaterLifeLimbJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalLaterLifeTailJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)

## Early Life Clocks
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalEarlyLifeJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalEarlyLifeLimbTailJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalEarlyLifeLimbJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgeLOO","AgeAccelLOO")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAgeN131FinalEarlyLifeTailJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, output.csv, ys.input1, ys.input2)

## HumanAxolotl Early Life LOFO Results (Absolute Age LOFO only --> for Training Data description purposes)
input1.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_LOFOOnly_HumanAxolotlN131_FinalEarlyLifeJosephZoller_LLin3.csv'
ys.output <- ys.input1 %>%
  dplyr::filter(!is.na(DNAmAgeLOFO10Balance)) # REMOVING NON-TRAINING SAMPLES
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, output.csv, ys.input1)
input1.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLifeLimbTail+Pan_Final_subCPGaxolotln131_EpigeneticLLin3Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_LOFOOnly_HumanAxolotlN131_FinalEarlyLifeLimbTail+PanJosephZoller_LLin3.csv'
ys.output <- ys.input1 %>%
  dplyr::filter(!is.na(DNAmAgeLOFO10Balance)) # REMOVING NON-TRAINING SAMPLES
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, output.csv, ys.input1)

## AxolotlClawedFrog Early Life LOO Results (Absolute Age LOO only --> for Training Data description purposes)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOEarlyLife_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_LOOOnly_AxolotlClawedFrogN131N140_FinalEarlyLifeJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  dplyr::filter(!is.na(DNAmAgeLOO)) # REMOVING NON-TRAINING SAMPLES
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, output.csv, ys.input1)
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOEarlyLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLog2Age_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_LOOOnly_AxolotlClawedFrogN131N140_FinalEarlyLifeLimbTail+PanJosephZoller_Log2.csv'
ys.output <- ys.input1 %>%
  dplyr::filter(!is.na(DNAmAgeLOO)) # REMOVING NON-TRAINING SAMPLES
ys.output <- dplyr::arrange(ys.output, SpeciesLatinName, Folder, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, output.csv, ys.input1)

## ExpmtLimbRepRegen Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRepRegen_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRepRegen_PredictedValues.csv'
input3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRepRegen_PredictedValues.csv'
input4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRepRegen_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)[, -1*grep("AgeAccel", colnames(read.csv(input1.csv, as.is=T)))]
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimbTail")]
ys.input3 <- read.csv(input3.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimb")]
ys.input4 <- read.csv(input4.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllTail")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_AxolotlN131EarlyLifetoAxolotlN131ExpmtLimbRepRegen_FinalCOMBINEDJosephZoller.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input3, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input4, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, input3.csv, input4.csv, output.csv, ys.input1, ys.input2, ys.input3, ys.input4)

## ExpmtTailRegen Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
input3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
input4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)[, -1*grep("AgeAccel", colnames(read.csv(input1.csv, as.is=T)))]
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimbTail")]
ys.input3 <- read.csv(input3.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimb")]
ys.input4 <- read.csv(input4.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllTail")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_AxolotlN131EarlyLifetoAxolotlN131ExpmtTailRegen_FinalCOMBINEDJosephZoller.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input3, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input4, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, input3.csv, input4.csv, output.csv, ys.input1, ys.input2, ys.input3, ys.input4)

## ExpmtDemethylation Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
input3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
input4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)#[, -1*grep("AgeAccel", colnames(read.csv(input1.csv, as.is=T)))]
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimbTail")]
ys.input3 <- read.csv(input3.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimb")]
ys.input4 <- read.csv(input4.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllTail")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_AxolotlN131EarlyLifetoAxolotlN131ExpmtDemethylation_FinalCOMBINEDJosephZoller.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input3, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input4, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, input3.csv, input4.csv, output.csv, ys.input1, ys.input2, ys.input3, ys.input4)

## AxolotlN131Branch2 Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131Branch2_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131Branch2_PredictedValues.csv'
input3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131Branch2_PredictedValues.csv'
input4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131Branch2_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)[, -1*grep("AgeAccel", colnames(read.csv(input1.csv, as.is=T)))]
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimbTail")]
ys.input3 <- read.csv(input3.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimb")]
ys.input4 <- read.csv(input4.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllTail")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_AxolotlN131EarlyLifetoAxolotlN131Branch2_FinalCOMBINEDJosephZoller.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input3, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input4, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, input3.csv, input4.csv, output.csv, ys.input1, ys.input2, ys.input3, ys.input4)

## AxolotlN146Early Life Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
input3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
input4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)[, -1*grep("AgeAccel", colnames(read.csv(input1.csv, as.is=T)))]
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimbTail")]
ys.input3 <- read.csv(input3.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimb")]
ys.input4 <- read.csv(input4.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllTail")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_AxolotlN131EarlyLifetoAxolotlN146_FinalCOMBINEDJosephZoller.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input3, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input4, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, input3.csv, input4.csv, output.csv, ys.input1, ys.input2, ys.input3, ys.input4)

## ExpmtTERTTailRegen Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTERTTailRegen_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTERTTailRegen_PredictedValues.csv'
input3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTERTTailRegen_PredictedValues.csv'
input4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTERTTailRegen_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)[, -1*grep("AgeAccel", colnames(read.csv(input1.csv, as.is=T)))]
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimbTail")]
ys.input3 <- read.csv(input3.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimb")]
ys.input4 <- read.csv(input4.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllTail")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_AxolotlN131EarlyLifetoAxolotlN131ExpmtTERTTailRegen_FinalCOMBINEDJosephZoller.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input3, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input4, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, input3.csv, input4.csv, output.csv, ys.input1, ys.input2, ys.input3, ys.input4)

## AxolotlN131SkinAndTestSet Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131SkinAndTestSet_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131SkinAndTestSet_PredictedValues.csv'
input3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131SkinAndTestSet_PredictedValues.csv'
input4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131SkinAndTestSet_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)[, -1*grep("AgeAccel", colnames(read.csv(input1.csv, as.is=T)))]
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimbTail")]
ys.input3 <- read.csv(input3.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimb")]
ys.input4 <- read.csv(input4.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllTail")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_AxolotlN131EarlyLifetoAxolotlN131SkinAndTestSet_FinalCOMBINEDJosephZoller.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input3, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input4, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, input3.csv, input4.csv, output.csv, ys.input1, ys.input2, ys.input3, ys.input4)

## ExpmtLimbPersistRegen Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
input3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
input4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)[, -1*grep("AgeAccel", colnames(read.csv(input1.csv, as.is=T)))]
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimbTail")]
ys.input3 <- read.csv(input3.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimb")]
ys.input4 <- read.csv(input4.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllTail")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_AxolotlN131EarlyLifetoAxolotlN131ExpmtLimbPersistRegen_FinalCOMBINEDJosephZoller.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input3, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input4, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, input3.csv, input4.csv, output.csv, ys.input1, ys.input2, ys.input3, ys.input4)

## ExpmtLimbRegen Predictions
input1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
input2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
input3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
input4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
ys.input1 <- read.csv(input1.csv, as.is=T)[, -1*grep("AgeAccel", colnames(read.csv(input1.csv, as.is=T)))]
ys.input2 <- read.csv(input2.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimbTail")]
ys.input3 <- read.csv(input3.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllLimb")]
ys.input4 <- read.csv(input4.csv, as.is=T)[, c("Basename","DNAmAgebasedOnAllTail")]
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datPredictedAge_AxolotlN131EarlyLifetoAxolotlN131ExpmtLimbRegen_FinalCOMBINEDJosephZoller.csv'
ys.output <- ys.input1 %>%
  base::merge(ys.input2, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input3, "Basename", all.x=T, sort=F) %>%
  base::merge(ys.input4, "Basename", all.x=T, sort=F)
ys.output <- dplyr::arrange(ys.output, OriginalOrderInBatch)
attr(ys.output,"row.names") <- seq.int(nrow(ys.output))
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
rm(input1.csv, input2.csv, input3.csv, input4.csv, output.csv, ys.input1, ys.input2, ys.input3, ys.input4)

### (3) Copying Images of Plots
## Figure 1
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOOCOMBINED_Final_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Fig1_LOOCOMBINED_Final_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Figure 2
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOFO_HumanAxolotlN131_COMBINED_Final_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Fig2_LOFO_HumanAxolotlN131_COMBINED_Final_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Figure 2
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOO_AxolotlClawedFrogN131N140_COMBINED_Final_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Fig2_LOO_AxolotlClawedFrogN131N140_COMBINED_Final_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Figure 6
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRepRegen_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Fig6_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRepRegen_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Figure 6 WHISKERPLOTS
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRepRegen_PANEL-WHISKERPLOTS.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Fig6_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRepRegen_PANEL-WHISKERPLOTS.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Table 6
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRepRegen.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Tab6_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRepRegen.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Figure 6
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Fig6_ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Figure 6 WHISKERPLOTS
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen_PANEL-WHISKERPLOTS.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Fig6_ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen_PANEL-WHISKERPLOTS.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Table 6
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Tab6_ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## (ALTERNATE) Figure 6
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen-ALL_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Fig6_ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen-ALL_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Figure 6 WHISKERPLOTS
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen-ALL_PANEL-WHISKERPLOTS.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Fig6_ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen-ALL_PANEL-WHISKERPLOTS.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Table 6
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen-ALL.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Tab6_ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen-ALL.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Extended Data Figure 4
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOOEarlyLife_AxolotlN131_Final_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig4_LOOEarlyLife_AxolotlN131_Final_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Extended Data Figure 5
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOO+LOFO10EarlyLifeCOMBINED_AxolotlN131_Final_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig5_LOO+LOFO10EarlyLifeCOMBINED_AxolotlN131_Final_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Extended Data Figure 6
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtCONTROLS_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig6_ClockCOMBINED_AxolotlN131_Final_toExpmtCONTROLS_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Extended Data Figure 6A
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig6A_ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Extended Data Table 6A
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataTab6A_ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## (ALTERNATE) Extended Data Figure 6A
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation-ALL_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig6A_ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation-ALL_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Extended Data Table 6A
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation-ALL.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataTab6A_ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation-ALL.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Extended Data Figure 6B
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig6B_ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Extended Data Table 6B
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataTab6B_ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## (ALTERNATE) Extended Data Figure 6B
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION-ALL_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig6B_ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION-ALL_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Extended Data Table 6B
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION-ALL.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataTab6B_ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION-ALL.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Extended Data Figure 6A
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig6A_MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Extended Data Table 6A
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataTab6A_MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## (ALTERNATE) Extended Data Figure 6A
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation-ALL_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig6A_MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation-ALL_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Extended Data Table 6A
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation-ALL.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataTab6A_MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation-ALL.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Extended Data Figure 6B
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig6B_MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Extended Data Table 6B
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataTab6B_MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## (ALTERNATE) Extended Data Figure 6B
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION-ALL_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig6B_MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION-ALL_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Extended Data Table 6B
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION-ALL.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataTab6B_MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION-ALL.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Extended Data Figure 12
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOOBeforeAfterTwoYearsCOMBINED_AxolotlN131_Final_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/ExtendDataFig12_LOOBeforeAfterTwoYearsCOMBINED_AxolotlN131_Final_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toAxolotlN131Branch2EarlyLife_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toAxolotlN131Branch2EarlyLife_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toAxolotlN146EarlyLife_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toAxolotlN146EarlyLife_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTERTTailRegen_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtTERTTailRegen_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Table
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTERTTailRegen.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalTab_ClockCOMBINED_AxolotlN131_Final_toExpmtTERTTailRegen.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toSkinAndTestSet_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toSkinAndTestSet_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_SkinAndTestSet_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_MeanMethylationCOMBINED_AxolotlN131_SkinAndTestSet_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Figure WHISKERPLOTS
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen_PANEL-WHISKERPLOTS.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen_PANEL-WHISKERPLOTS.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Table
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalTab_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## (ALTERNATE) Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen-ALL_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen-ALL_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Internal Figure WHISKERPLOTS
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen-ALL_PANEL-WHISKERPLOTS.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen-ALL_PANEL-WHISKERPLOTS.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Internal Table
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen-ALL.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalTab_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen-ALL.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Figure WHISKERPLOTS
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen_PANEL-WHISKERPLOTS.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen_PANEL-WHISKERPLOTS.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Table
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalTab_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## (ALTERNATE) Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen-ALL_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen-ALL_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Internal Figure WHISKERPLOTS
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen-ALL_PANEL-WHISKERPLOTS.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen-ALL_PANEL-WHISKERPLOTS.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Internal Table
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen-ALL.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalTab_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen-ALL.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Figure WHISKERPLOTS
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded_PANEL-WHISKERPLOTS.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded_PANEL-WHISKERPLOTS.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## Internal Table
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalTab_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## (ALTERNATE) Internal Figure
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded-ALL_PANEL.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded-ALL_PANEL.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Internal Figure WHISKERPLOTS
in.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded-ALL_PANEL-WHISKERPLOTS.pdf'
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalFig_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded-ALL_PANEL-WHISKERPLOTS.pdf'
file.copy(in.pdf, out.pdf, overwrite=T)
rm(in.pdf, out.pdf)
## (ALTERNATE) Internal Table
in.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded-ALL.xlsx'
out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/InternalTab_ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded-ALL.xlsx'
file.copy(in.xlsx, out.xlsx, overwrite=T)
rm(in.xlsx, out.xlsx)
## Unused Figures
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_Clock_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_Clock_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockLaterLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockLaterLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockLaterLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockLaterLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockLaterLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOO_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOO_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOO_Final_subCPGaxolotln131_EpigeneticLog2Age_PANEL.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOO_Final_subCPGaxolotln131_EpigeneticLog2Age_PANEL.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLife_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOLaterLife_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PANEL.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOLaterLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PANEL.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOLaterLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOLaterLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOLaterLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
## IDENTICAL to Extended Data Figure 4
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PANEL.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PANEL.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOEarlyLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOEarlyLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_LOOEarlyLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)

### (4) Saving EWAS Tables
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeEWAS_Final_Analysis/Subset_AxolotlN131_EWASxTissue_Final_subCPGaxolotln131_Log2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datEWASxTissueN131FinalJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeEWAS_Final_Analysis/Subset_AxolotlN131_EWASxTissue_Final_subCPGaxolotln131_Log2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_EWASxTissue_Final_subCPGaxolotln131_Log2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeEWAS_Final_Analysis/Subset_AxolotlN131_EWASLaterLifexTissue_Final_subCPGaxolotln131_Log2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datEWASxTissueN131FinalLaterLifeJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeEWAS_Final_Analysis/Subset_AxolotlN131_EWASLaterLifexTissue_Final_subCPGaxolotln131_Log2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_EWASLaterLifexTissue_Final_subCPGaxolotln131_Log2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)
in.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeEWAS_Final_Analysis/Subset_AxolotlN131_EWASEarlyLifexTissue_Final_subCPGaxolotln131_Log2Age.csv'
in.valbeta <- read.csv(in.csv)
out.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/datEWASxTissueN131FinalEarlyLifeJosephZoller_Log2.csv'
write.table(in.valbeta,out.csv,sep=',',row.names=F,quote=F)
rm(in.csv, out.csv)
in.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeEWAS_Final_Analysis/Subset_AxolotlN131_EWASEarlyLifexTissue_Final_subCPGaxolotln131_Log2Age.png'
out.png='SpeciesSubsetAnalyses/Subset_AxolotlN131_PredictionResultsJosephZoller/Subset_AxolotlN131_EWASEarlyLifexTissue_Final_subCPGaxolotln131_Log2Age.png'
file.copy(in.png, out.png, overwrite=T)
rm(in.png, out.png)



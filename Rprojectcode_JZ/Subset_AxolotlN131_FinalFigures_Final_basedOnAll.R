rm(list=ls())
options(stringAsFactors=F)
library(tidyverse)
library(glmnet)
library(WGCNA)
setwd("~/Dropbox/MyResearchFiles/Horvath_mammalian_meth")
library(devtools)
library(MammalMethylClock)

### FOR DATA DEPOSITION: We must report which samples were used for the paper
Basename.list_data_depo <- list()
#############################################################################

###############################################################################
### Figure 1: Plotting all pan-tissue and tissue-specific LOO clocks together
###############################################################################
panel.mains=c('Axolotl PanTissue','Axolotl LimbTail','Axolotl Limb','Axolotl Tail',
              'Axolotl Early Life PanTissue','Axolotl Early Life LimbTail','Axolotl Early Life Limb','Axolotl Early Life Tail')#,
              #'Axolotl Later Life PanTissue','Axolotl Later Life LimbTail','Axolotl Later Life Limb','Axolotl Later Life Tail')
y.axis.labs=c('DNAmAgeLOO Pan Tissue','DNAmAgeLOO LimbTail','DNAmAgeLOO Limb','DNAmAgeLOO Tail',
              'DNAmAgeLOO Early Life Pan Tissue','DNAmAgeLOO Early Life LimbTail','DNAmAgeLOO Early Life Limb','DNAmAgeLOO Early Life Tail')#,
              #'DNAmAgeLOO Later Life Pan Tissue','DNAmAgeLOO Later Life LimbTail','DNAmAgeLOO Later Life Limb','DNAmAgeLOO Later Life Tail')
x.axis.labs=c('Age','Age','Age','Age',
              'Age','Age','Age','Age')#,
              #'Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOO_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output5.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output6.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output7.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
output8.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
# output9.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
# output10.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
# output11.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
# output12.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOLaterLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgeLOO)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgeLOO)
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgeLOO)
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgeLOO)
input5.info_pred <- dplyr::filter(read.csv(output5.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=5, Outcome=Age, Prediction=DNAmAgeLOO)
input6.info_pred <- dplyr::filter(read.csv(output6.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=6, Outcome=Age, Prediction=DNAmAgeLOO)
input7.info_pred <- dplyr::filter(read.csv(output7.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=7, Outcome=Age, Prediction=DNAmAgeLOO)
input8.info_pred <- dplyr::filter(read.csv(output8.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=8, Outcome=Age, Prediction=DNAmAgeLOO)
# input9.info_pred <- dplyr::filter(read.csv(output9.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
#   dplyr::mutate(PanelName=9, Outcome=Age, Prediction=DNAmAgeLOO)
# input10.info_pred <- dplyr::filter(read.csv(output10.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
#   dplyr::mutate(PanelName=10, Outcome=Age, Prediction=DNAmAgeLOO)
# input11.info_pred <- dplyr::filter(read.csv(output11.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
#   dplyr::mutate(PanelName=11, Outcome=Age, Prediction=DNAmAgeLOO)
# input12.info_pred <- dplyr::filter(read.csv(output12.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
#   dplyr::mutate(PanelName=12, Outcome=Age, Prediction=DNAmAgeLOO)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred,
                             input5.info_pred, input6.info_pred, input7.info_pred, input8.info_pred) %>%#,
                             #input9.info_pred, input10.info_pred, input11.info_pred, input12.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("Tissue","Tissue","Tissue","Tissue",
             "Tissue","Tissue","Tissue","Tissue")#,
             #"Tissue","Tissue","Tissue","Tissue")#c("SpeciesTissue","SpeciesTissue","SpeciesTissue","SpeciesTissue")
NUMVAR.VEC=c(NA,NA,NA,NA,
             NA,NA,NA,NA)#,
             #NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOOCOMBINED_Final_PANEL.pdf'
out.pdf.title='Leave-One-Out Analysis of All Final Epigenetic Clocks'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,4)#c(3,4)
width=17
height=10#14
panel.labs=letters
oma.right=0
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2) {
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  ## NEW, CUSTOM b/c reviewer request
  MeanAE <- mean(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MeanAE_str <- paste0("MeanAE=",signif(MeanAE,3))
  MAE_str <- paste0(MAE_str,', ',MeanAE_str)
  ##
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 1) {
    legend('topleft',
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["LOOCOMBINED_Final_PANEL"]] <- unique(input.info_pred$Basename)
#######################################
#######################################

###############################################################################
### Figure 2: Plotting all HumanAxolotl LOFO clocks together
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Human-Axolotl Early Life PanTissue','Human-Axolotl Early Life PanTissue',
              'Human-Axolotl Early Life LimbTail+Pan','Human-Axolotl Early Life LimbTail+Pan')
y.axis.labs=c('DNAmAgeLOO Early Life Pan Tissue','DNAmRelativeAgeLOO Early Life Pan Tissue',
              'DNAmAgeLOO Early Life LimbTail+PanTissue','DNAmRelativeAgeLOO Early Life LimbTail+PanTissue')
x.axis.labs=c('Age','RelativeAge',
              'Age','RelativeAge')
output1.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticRelativeAge_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLifeLimbTail+Pan_Final_subCPGaxolotln131_EpigeneticLLin3Age_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLifeLimbTail+Pan_Final_subCPGaxolotln131_EpigeneticRelativeAge_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgeLOFO10Balance)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgeLOFO10Balance)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmRelAgeLOFO10Balance)) %>%
  dplyr::mutate(PanelName=2, Outcome=RelAge, Prediction=DNAmRelAgeLOFO10Balance)
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgeLOFO10Balance)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgeLOFO10Balance)
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmRelAgeLOFO10Balance)) %>%
  dplyr::mutate(PanelName=4, Outcome=RelAge, Prediction=DNAmRelAgeLOFO10Balance)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred,
                             input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum','Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("SpeciesTissue","SpeciesTissue",
             "SpeciesTissue","SpeciesTissue")
NUMVAR.VEC=c(NA,NA,
             NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOFO_HumanAxolotlN131_COMBINED_Final_PANEL.pdf'
out.pdf.title='Leave-One-Fold-Out Analysis of All Final Human-Axolotl Epigenetic Clocks'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=9+2.3 #wider for exterior legend
height=10
panel.labs=letters
oma.right=17
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2) {
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  ## NEW, CUSTOM b/c reviewer request
  MeanAE <- mean(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MeanAE_str <- paste0("MeanAE=",signif(MeanAE,3))
  MAE_str <- paste0(MAE_str,', ',MeanAE_str)
  ##
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 2) {
    legend('right',inset=-1.00,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["LOFO_HumanAxolotlN131_COMBINED_Final_PANEL"]] <- setdiff(unique(input.info_pred$Basename),
                                                                                   Basename.list_data_depo[["LOOCOMBINED_Final_PANEL"]])
#######################################
#######################################

###############################################################################
### Figure 2: Plotting all AxolotlClawedFrog LOO clocks together
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl-Frog Early Life PanTissue','Axolotl-Frog Early Life PanTissue',
              'Axolotl-Frog Early Life LimbTail+Pan','Axolotl-Frog Early Life LimbTail+Pan')
y.axis.labs=c('DNAmAgeLOO Early Life Pan Tissue','DNAmAgeLOO Early Life Pan Tissue',
              'DNAmAgeLOO Early Life LimbTail+PanTissue','DNAmAgeLOO Early Life LimbTail+PanTissue')
x.axis.labs=c('Age','RelativeAge',
              'Age','RelativeAge')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOEarlyLife_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLog2Age_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOEarlyLife_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticRelativeAge_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOEarlyLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLog2Age_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOEarlyLifeLimbTail+Pan_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticRelativeAge_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgeLOO)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmRelAgeLOO)) %>%
  dplyr::mutate(PanelName=2, Outcome=RelAge, Prediction=DNAmRelAgeLOO)
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgeLOO)
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmRelAgeLOO)) %>%
  dplyr::mutate(PanelName=4, Outcome=RelAge, Prediction=DNAmRelAgeLOO)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred,
                             input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum','Xenopus laevis','Xenopus tropicalis'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue,
                                                             ignore=c("Xenopus laevis","Xenopus tropicalis")))

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("SpeciesTissue","SpeciesTissue",
             "SpeciesTissue","SpeciesTissue")
NUMVAR.VEC=c(NA,NA,
             NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOO_AxolotlClawedFrogN131N140_COMBINED_Final_PANEL.pdf'
out.pdf.title='Leave-One-Out Analysis of All Final Axolotl-Clawed Frog Epigenetic Clocks'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=9+2.3 #wider for exterior legend
height=10
panel.labs=letters
oma.right=17
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2) {
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  ## NEW, CUSTOM b/c reviewer request
  MeanAE <- mean(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MeanAE_str <- paste0("MeanAE=",signif(MeanAE,3))
  MAE_str <- paste0(MAE_str,', ',MeanAE_str)
  ##
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 2) {
    legend('right',inset=-1.00,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["LOO_AxolotlClawedFrogN131N140_COMBINED_Final_PANEL"]] <- setdiff(unique(input.info_pred$Basename),
                                                                                           Basename.list_data_depo[["LOOCOMBINED_Final_PANEL"]])
#######################################
#######################################

###############################################################################
### Figure 6 + Table 6: Plotting all AxolotlN131 clocks applied to AxolotlN131 Repetitive Limb Regen
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Repetitive Limb Regen','Axolotl Repetitive Limb Regen','Axolotl Repetitive Limb Regen','Axolotl Repetitive Limb Regen')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRepRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRepRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRepRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRepRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup,
                                               levels=c("MatureLeftLimb","RegeneratedLeftLimbRound1",
                                                        "RegeneratedLeftLimbRound2","RegeneratedLeftLimbRound3",
                                                        "MatureRightLimb(AgedControlForThriceRegenerated)"))
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
NUMVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRepRegen_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Repetitive Limb Regen Experiment Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=13.0 #wider for exterior legend
height=10
panel.labs=letters
oma.right=28
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRepRegen.xlsx'

#######################################
#######################################
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, AnimalID, RegenExperimentGroup, Prediction) %>%
  dplyr::group_by(Clock, AnimalID) %>%
  tidyr::spread(RegenExperimentGroup, Prediction)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  #reporting t-test p-value for comparing treatment vs. control at final time step
  rows_ij_ttest <- which(as.numeric(ys.panelfactor)==i &
                           ys.colfactor %in% c("RegeneratedLeftLimbRound3",
                                               "MatureRightLimb(AgedControlForThriceRegenerated)"))
  pval <- t.test(ys.prediction[rows_ij_ttest] ~ ys.colfactor[rows_ij_ttest], na.rm=T)$p.value
  pval_str <- paste0("Round3 p-value=",signif(pval,2))
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  title(pval_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.65,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRepRegen_PANEL"]] <- unique(input.info_pred$Basename)
#######################################
#######################################

########  WHISKER PLOTS (ALT)  ########
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(sub("(.*)(\\.)", "\\1-WHISKERPLOTS\\2", out.pdf),
    width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  #reporting t-test p-value for comparing treatment vs. control at final time step
  rows_ij_ttest <- which(as.numeric(ys.panelfactor)==i &
                           ys.colfactor %in% c("RegeneratedLeftLimbRound3",
                                               "MatureRightLimb(AgedControlForThriceRegenerated)"))
  pval <- t.test(ys.prediction[rows_ij_ttest] ~ ys.colfactor[rows_ij_ttest], na.rm=T)$p.value
  pval_str <- paste0("Round3 p-value=",signif(pval,2))
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]

  df_box <- data.frame(y=ys.prediction[rows_i], g=ys.colfactor[rows_i], Age=ys.outcome[rows_i]) %>%
    dplyr::arrange(g, Age) %>% dplyr::group_by(g) %>%
    dplyr::summarise(x=mean(range(Age)),
                     median=median(y),
                     lower=median(y)-IQR(y),
                     upper=median(y)+IQR(y))
  plot(median~x,data=df_box,
       ylim=range(c(lower,upper),na.rm=T),
       xlim=axis_limits(x),
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab="Experiment Group",
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors,pch=19)
  arrows(df_box$x,df_box$lower,df_box$x,df_box$upper,
         length=0.05,angle=90,code=3,
         col=ys.colors)
  # plot(pch='',y=range(df_box$y),x=axis_limits(df_box$Age),
  #      main=paste0(PANEL_str,' (N=',N,')'),
  #      ylab=ylab,xlab="Experiment Group",
  #      cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf)
  # boxplot(add=T,y~as.numeric(g),data=df_box,at=c(0.4794521,0.9863014,0.9863014,1.2005),
  #         cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #         axes=F,boxwex=0.05,col=ys.colors)

  title(pval_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.65,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Figure 6 + Table 6: Plotting all AxolotlN131 clocks applied to AxolotlN131 Tail Regen
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Tail Regen All','Axolotl Tail Regen All','Axolotl Tail Regen All','Axolotl Tail Regen All')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, Experiment, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$Experiment <- factor(input.info_pred$Experiment)
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup)
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, Experiment, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("Experiment","Experiment","Experiment","Experiment")
NUMVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
ys.colors_original <- NA
ys.numbers_original <- c("1","2","3","4","5","6","7","8")
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Tail Regen Experiment Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=13.0 #wider for exterior legend
height=10
panel.labs=letters
oma.right=28
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, Experiment, AnimalID, RegenExperimentGroup, Prediction)
## handling replicates: giving distinct labels
table.input_pred <- table.input_pred %>% dplyr::mutate(RegenExperimentGroup = as.character(RegenExperimentGroup)) %>%
  dplyr::group_by(Clock, AnimalID, RegenExperimentGroup, Experiment) %>%
  dplyr::mutate(RegenExperimentGroup = ifelse(row_number()==1, RegenExperimentGroup, paste(RegenExperimentGroup, row_number(), sep="-"))) %>%
  dplyr::ungroup()
##
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Clock, Experiment, AnimalID) %>%
  tidyr::spread(RegenExperimentGroup, Prediction) %>%
  dplyr::arrange(Clock, Experiment, AnimalID)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8,
         font=2)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.60,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen_PANEL"]] <- unique(input.info_pred$Basename)
#######################################
#######################################

########  WHISKER PLOTS (ALT)  ########
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(sub("(.*)(\\.)", "\\1-WHISKERPLOTS\\2", out.pdf),
    width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]

  ## Handling two layers of grouping: repeating boxes and colors ##
  df_box <- data.frame(y=ys.prediction[rows_i], g1=ys.colfactor[rows_i], g2=ys.numfactor[rows_i], Age=ys.outcome[rows_i]) %>%
    dplyr::arrange(g1, g2, Age) %>% dplyr::mutate(g2 = as.numeric(g2))
  #### separating Tail controls by age ####
  df_box$g2[which(df_box$g1 %in% c("TailRegenerationAgedControl") & df_box$Age < 0.9)] <- 1
  df_box$g2[which(df_box$g1 %in% c("TailRegenerationAgedControl") & df_box$Age > 0.9)] <- 2
  ########
  df_box$g <- paste(df_box$g1, df_box$g2, sep="-")
  df_box$g <- factor(df_box$g, levels=unique(df_box$g))
  ####
  df_box <- df_box %>% dplyr::group_by(g1,g) %>%
    dplyr::summarise(x=mean(range(Age)),
                     median=median(y),
                     lower=median(y)-IQR(y),
                     upper=median(y)+IQR(y))
  plot(median~x,data=df_box,
       ylim=range(c(lower,upper),na.rm=T),
       xlim=axis_limits(x),
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab="Experiment Group",
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors[df_box$g1],pch=19)
  arrows(df_box$x,df_box$lower,df_box$x,df_box$upper,
         length=0.05,angle=90,code=3,
         col=ys.colors[df_box$g1])
  #Creating spaghetti plots instead of scatter plots
  for (j in 2:length(levels(df_box$g1))) {
    expmt <- levels(df_box$g1)[j]
    rows_i_expmt <- which(df_box$g1==expmt)
    lines(y=df_box$median[rows_i_expmt],x=df_box$x[rows_i_expmt],
          lty='dashed',lwd=0.4,col=ys.colors[j])
  }
  # boxplot(y~g, data=df_box,
  #         main=paste0(PANEL_str,' (N=',N,')'),
  #         ylab=ylab,xlab="Experiment Group",names=NA,
  #         cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #         col=rep(ys.colors, times=table(unique(dplyr::select(df_box, g1, g2))$g1)))

  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.60,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### (ALTERNATE) Figure 6 + Table 6: Plotting all AxolotlN131 clocks applied to AxolotlN131 Tail Regen (ALL SAMPLES)
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Tail Regen All','Axolotl Tail Regen All','Axolotl Tail Regen All','Axolotl Tail Regen All')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTailRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail)
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb)
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, Experiment, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$Experiment <- factor(input.info_pred$Experiment)
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup)
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, Experiment, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("Experiment","Experiment","Experiment","Experiment")
NUMVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
ys.colors_original <- NA
ys.numbers_original <- c("1","2","3","4","5","6","7","8")
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen-ALL_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Tail Regen Experiment Data (ALL SAMPLES)'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=13.0 #wider for exterior legend
height=10
panel.labs=letters
oma.right=28
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTailRegen-ALL.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, Experiment, AnimalID, RegenExperimentGroup, Prediction)
## handling replicates: giving distinct labels
table.input_pred <- table.input_pred %>% dplyr::mutate(RegenExperimentGroup = as.character(RegenExperimentGroup)) %>%
  dplyr::group_by(Clock, AnimalID, RegenExperimentGroup, Experiment) %>%
  dplyr::mutate(RegenExperimentGroup = ifelse(row_number()==1, RegenExperimentGroup, paste(RegenExperimentGroup, row_number(), sep="-"))) %>%
  dplyr::ungroup()
##
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Clock, Experiment, AnimalID) %>%
  tidyr::spread(RegenExperimentGroup, Prediction) %>%
  dplyr::arrange(Clock, Experiment, AnimalID)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8,
         font=2)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.60,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

########  WHISKER PLOTS (ALT)  ########
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(sub("(.*)(\\.)", "\\1-WHISKERPLOTS\\2", out.pdf),
    width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]

  ## Handling two layers of grouping: repeating boxes and colors ##
  df_box <- data.frame(y=ys.prediction[rows_i], g1=ys.colfactor[rows_i], g2=ys.numfactor[rows_i], Age=ys.outcome[rows_i]) %>%
    dplyr::arrange(g1, g2, Age) %>% dplyr::mutate(g2 = as.numeric(g2))
  #### separating Tail controls by age ####
  df_box$g2[which(df_box$g1 %in% c("TailRegenerationAgedControl") & df_box$Age < 0.9)] <- 1
  df_box$g2[which(df_box$g1 %in% c("TailRegenerationAgedControl") & df_box$Age > 0.9)] <- 2
  ########
  df_box$g <- paste(df_box$g1, df_box$g2, sep="-")
  df_box$g <- factor(df_box$g, levels=unique(df_box$g))
  ####
  df_box <- df_box %>% dplyr::group_by(g1,g) %>%
    dplyr::summarise(x=mean(range(Age)),
                     median=median(y),
                     lower=median(y)-IQR(y),
                     upper=median(y)+IQR(y))
  plot(median~x,data=df_box,
       ylim=range(c(lower,upper),na.rm=T),
       xlim=axis_limits(x),
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab="Experiment Group",
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors[df_box$g1],pch=19)
  arrows(df_box$x,df_box$lower,df_box$x,df_box$upper,
         length=0.05,angle=90,code=3,
         col=ys.colors[df_box$g1])
  #Creating spaghetti plots instead of scatter plots
  for (j in 2:length(levels(df_box$g1))) {
    expmt <- levels(df_box$g1)[j]
    rows_i_expmt <- which(df_box$g1==expmt)
    lines(y=df_box$median[rows_i_expmt],x=df_box$x[rows_i_expmt],
          lty='dashed',lwd=0.4,col=ys.colors[j])
  }
  # boxplot(y~g, data=df_box,
  #         main=paste0(PANEL_str,' (N=',N,')'),
  #         ylab=ylab,xlab="Experiment Group",names=NA,
  #         cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #         col=rep(ys.colors, times=table(unique(dplyr::select(df_box, g1, g2))$g1)))

  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.60,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Extended Data Figure 4: Plotting pan-tissue Early Life LOO clock by Tissue (RE-PLOTTING FOR PDF)
###############################################################################

## NOTE: IDENTICAL to 'SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PANEL.png'

rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
y.axis.labs=rep('DNAmAgeLOO',100)
x.axis.labs=rep('Age',100)
output.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
input0.info_pred <- dplyr::filter(read.csv(output.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::select(Basename, Age, DNAmAgeLOO, Tissue) %>%
  dplyr::mutate(PanelName = "All")
input1.info_pred <- dplyr::filter(read.csv(output.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::select(Basename, Age, DNAmAgeLOO, Tissue) %>%
  dplyr::mutate(PanelName = Tissue)
input.info_pred <- bind_rows(input0.info_pred, input1.info_pred)
input.info_pred$PanelName <- factor(input.info_pred$PanelName) #exploiting the fact that "All" comes before all tissue names alphabetically
input.info_pred$Tissue <- factor(input.info_pred$Tissue)

panel.mains=levels(input.info_pred$PanelName)

OUTVAR="Age"
PREDVAR="DNAmAgeLOO"
PANELVAR="PanelName"
COLVAR.VEC=rep("Tissue",100)
NUMVAR.VEC=rep(NA,100)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOOEarlyLife_AxolotlN131_Final_PANEL.pdf'
out.pdf.title='DNAmAgeLOO for Axolotl, by Tissue'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(3,3)
width=13
height=14
panel.labs=letters
oma.right=0
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2) {
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  ## NEW, CUSTOM b/c reviewer request
  MeanAE <- mean(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MeanAE_str <- paste0("MeanAE=",signif(MeanAE,3))
  MAE_str <- paste0(MAE_str,', ',MeanAE_str)
  ##
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["LOOEarlyLife_AxolotlN131_Final_PANEL"]] <- setdiff(unique(input.info_pred$Basename),
                                                                             Basename.list_data_depo[["LOOCOMBINED_Final_PANEL"]])
#######################################
#######################################

###############################################################################
### Extended Data Figure 5: Plotting Early Life LOO and LOFO clocks together
###############################################################################
panel.mains=c('Axolotl Early Life PanTissue','Axolotl Early Life LimbTail','Axolotl Early Life Limb','Axolotl Early Life Tail',
              'Axolotl Early Life PanTissue','Axolotl Early Life LimbTail','Axolotl Early Life Limb','Axolotl Early Life Tail')
y.axis.labs=c('DNAmAgeLOO Early Life Pan Tissue','DNAmAgeLOO Early Life LimbTail','DNAmAgeLOO Early Life Limb','DNAmAgeLOO Early Life Tail',
              'DNAmAgeLOFO10 Early Life Pan Tissue','DNAmAgeLOFO10 Early Life LimbTail','DNAmAgeLOFO10 Early Life Limb','DNAmAgeLOFO10 Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age',
              'Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOEarlyLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output5.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOFO_Final_Analysis/Subset_AxolotlN131_LOFO10EarlyLife_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output6.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOFO_Final_Analysis/Subset_AxolotlN131_LOFO10EarlyLifeLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output7.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOFO_Final_Analysis/Subset_AxolotlN131_LOFO10EarlyLifeLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
output8.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOFO_Final_Analysis/Subset_AxolotlN131_LOFO10EarlyLifeTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgeLOO)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgeLOO)
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgeLOO)
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgeLOO)
input5.info_pred <- dplyr::filter(read.csv(output5.csv, as.is=T), !is.na(DNAmAgeLOFO10)) %>%
  dplyr::mutate(PanelName=5, Outcome=Age, Prediction=DNAmAgeLOFO10)
input6.info_pred <- dplyr::filter(read.csv(output6.csv, as.is=T), !is.na(DNAmAgeLOFO10)) %>%
  dplyr::mutate(PanelName=6, Outcome=Age, Prediction=DNAmAgeLOFO10)
input7.info_pred <- dplyr::filter(read.csv(output7.csv, as.is=T), !is.na(DNAmAgeLOFO10)) %>%
  dplyr::mutate(PanelName=7, Outcome=Age, Prediction=DNAmAgeLOFO10)
input8.info_pred <- dplyr::filter(read.csv(output8.csv, as.is=T), !is.na(DNAmAgeLOFO10)) %>%
  dplyr::mutate(PanelName=8, Outcome=Age, Prediction=DNAmAgeLOFO10)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred,
                             input5.info_pred, input6.info_pred, input7.info_pred, input8.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("Tissue","Tissue","Tissue","Tissue",
             "Tissue","Tissue","Tissue","Tissue")#c("SpeciesTissue","SpeciesTissue","SpeciesTissue","SpeciesTissue")
NUMVAR.VEC=c(NA,NA,NA,NA,
             NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOO+LOFO10EarlyLifeCOMBINED_AxolotlN131_Final_PANEL.pdf'
out.pdf.title='Leave-One-Out and Leave-One-Fold-Out (10-Fold) Analyses of Early Life Final Epigenetic Clocks'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,4)
width=17
height=10
panel.labs=letters
oma.right=0
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2) {
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  ## NEW, CUSTOM b/c reviewer request
  MeanAE <- mean(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MeanAE_str <- paste0("MeanAE=",signif(MeanAE,3))
  MAE_str <- paste0(MAE_str,', ',MeanAE_str)
  ##
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 1) {
    legend('topleft',
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["LOO+LOFOEarlyLifeCOMBINED_AxolotlN131_Final_PANEL"]] <- unique(input.info_pred$Basename)
#######################################
#######################################

###############################################################################
### Extended Data Figure 6: Plotting all AxolotlN131 clocks applied to AxolotlN131 Select Expmt Controls
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Expmt Controls','Axolotl Expmt Controls')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail')
x.axis.labs=c('Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("PersistenceOfLimbRegenerativeRejuvenation","TailRegenerationAgedControl")) %>%
  dplyr::filter(RegenExperimentGroup %in% c("MatureLeftLimb(AgedControl)","OriginalTail")) %>%
  dplyr::arrange(Experiment, RegenExperimentGroup)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("PersistenceOfLimbRegenerativeRejuvenation","TailRegenerationAgedControl")) %>%
  dplyr::filter(RegenExperimentGroup %in% c("MatureLeftLimb(AgedControl)","OriginalTail")) %>%
  dplyr::arrange(Experiment, RegenExperimentGroup)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, Experiment, RegenExperimentGroup, ExternalSampleID, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup)

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup")#c("SpeciesTissue","SpeciesTissue")
NUMVAR.VEC=c(NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtCONTROLS_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Experiment Data Control Samples'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(1,2)
width=9+1.8 #wider for exterior legend
height=5.5
panel.labs=letters
oma.right=13
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2) {
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  lim <- axis_square_limits(c(0.6,2.0),c(0.6,2.0)) #axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i]) # axes where requested to expanded manually
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  ## NEW, CUSTOM b/c reviewer request
  MeanAE <- mean(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MeanAE_str <- paste0("MeanAE=",signif(MeanAE,3))
  MAE_str <- paste0(MAE_str,', ',MeanAE_str)
  ##
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 2) {
    legend('right',inset=-1.10,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["ClockCOMBINED_AxolotlN131_Final_toExpmtCONTROLS_PANEL"]] <- unique(input.info_pred$Basename)
#######################################
#######################################

###############################################################################
### Extended Data Figure 6A + Extended Data Table 6A: Plotting all AxolotlN131 clocks applied to AxolotlN131 Demethylation vs. Concentration
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains_list=list(c('Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells'),
                      c('Axolotl Demethylation Larva','Axolotl Demethylation Larva','Axolotl Demethylation Larva','Axolotl Demethylation Larva'))
y.axis.labs_list=list(c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail'),
                      c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail'))
x.axis.labs_list=list(c('Concentration','Concentration','Concentration','Concentration'),
                      c('Concentration','Concentration','Concentration','Concentration'))
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input5.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=5, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input6.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=6, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input7.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=7, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input8.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=8, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
#Creating 2 separate data frames, because "Concentration" has different factor levels in each experimental set
inputCell.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Prediction, PanelName, Experiment, Concentration, Passage, ExternalSampleID) %>%
  dplyr::mutate(Concentration = ifelse(is.na(Concentration), "None", Concentration))
inputCell.info_pred$PanelName <- factor(inputCell.info_pred$PanelName)
inputCell.info_pred$Experiment <- factor(inputCell.info_pred$Experiment)
inputCell.info_pred$Concentration <- relevel(factor(inputCell.info_pred$Concentration), ref="None")
inputCell.info_pred$Passage <- factor(inputCell.info_pred$Passage)
inputCell.info_pred$ExternalSampleID <- factor(inputCell.info_pred$ExternalSampleID)
inputLarva.info_pred <- bind_rows(input5.info_pred, input6.info_pred, input7.info_pred, input8.info_pred) %>%
  dplyr::select(Basename, Prediction, PanelName, Experiment, Concentration, Duration, ExternalSampleID) %>%
  dplyr::filter(Duration %in% c(12,14))
inputLarva.info_pred$PanelName <- factor(inputLarva.info_pred$PanelName)
inputLarva.info_pred$Experiment <- factor(inputLarva.info_pred$Experiment)
inputLarva.info_pred$Concentration <- factor(inputLarva.info_pred$Concentration)
inputLarva.info_pred$Duration <- factor(inputLarva.info_pred$Duration)
inputLarva.info_pred$ExternalSampleID <- factor(inputLarva.info_pred$ExternalSampleID)

inputCell.info_pred <- dplyr::arrange(inputCell.info_pred, PanelName, Experiment, Concentration, Passage, ExternalSampleID) #ensure correct order
inputLarva.info_pred <- dplyr::arrange(inputLarva.info_pred, PanelName, Experiment, Concentration, Duration, ExternalSampleID) #ensure correct order
GROUPVAR="Concentration"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC_list=list(c(NA,NA,NA,NA),
                     c(NA,NA,NA,NA))
ys.colors_original_list <- list(NA,NA)
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Demethylation Experiment Data'
ys.output_list <- list(inputCell.info_pred, inputLarva.info_pred)
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,4)
width=17
height=10
panel.labs_list=list(letters[1:4],
                     letters[5:8])#letters
oma.right=0
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred_list = list(length(ys.output_list))
for (j in 1:length(ys.output_list)) {
  input.info_pred <- ys.output_list[[j]]

  if (j == 1) {
    table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, Concentration, Prediction)
    table.input_pred <- table.input_pred %>%
      dplyr::group_by(Clock, Concentration) %>%
      dplyr::mutate(Num = row_number()) %>%
      tidyr::spread(Num, Prediction)
  }
  if (j == 2) {
    table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, Concentration, Prediction)
    table.input_pred <- table.input_pred %>%
      dplyr::group_by(Clock, Concentration) %>%
      dplyr::mutate(Num = row_number()) %>%
      tidyr::spread(Num, Prediction)
  }
  levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
  table.input_pred_list[[j]] <- table.input_pred

  rm(input.info_pred,table.input_pred)
}
write.xlsx(table.input_pred_list,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold"),
           sheetName=c("DemethylationTreatmentCells","DemethylationTreatmentLarva")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (j in 1:length(ys.output_list)) {
  #Extracting vectors from vector-lists, for each experiment separately
  ys.output <- ys.output_list[[j]]
  panel.mains <- panel.mains_list[[j]]
  y.axis.labs <- y.axis.labs_list[[j]]
  x.axis.labs <- x.axis.labs_list[[j]]
  COLVAR.VEC <- COLVAR.VEC_list[[j]]
  ys.colors_original <- ys.colors_original_list[[j]]
  panel.labs <- panel.labs_list[[j]]

  ys.group <- ys.output[,GROUPVAR]
  ys.prediction <- ys.output[,PREDVAR]
  ys.panelfactor <- ys.output[,PANELVAR]
  for (i in 1:length(levels(ys.panelfactor))) {
    COLVAR=COLVAR.VEC[i]
    ys.colors <- ys.colors_original
    if (!is.na(COLVAR)) {
      ys.colfactor <- ys.output[,COLVAR]
      # ys.colors contains the palette of distinct colors
      # ys.colors_vec contains the color assignment for every row in the data set
      if (is.na(ys.colors[1])) {
        if (length(levels(ys.colfactor)) <= 2){
          ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
        } else if (length(levels(ys.colfactor)) <= 8) {
          ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
        } else {
          ys.colors <- rainbow(length(levels(ys.colfactor)))
        }
      }
    } else {
      ys.colfactor <- rep(1, nrow(ys.output))
      ys.colors <- 'black'
    }
    ys.colors_vec <- ys.colors[ys.colfactor]

    rows_i <- which(as.numeric(ys.panelfactor)==i)
    PVAL <- kruskal.test(ys.prediction[rows_i], ys.group[rows_i])$p.value
    PVAL_str <- paste0("p=",signif(PVAL,2))
    PANEL_str <- panel.mains[i]
    N <- length(which(!is.na(ys.group[rows_i]) & !is.na(ys.prediction[rows_i])))
    ylab=y.axis.labs[i]
    xlab=x.axis.labs[i]
    plab=panel.labs[i]
    boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
            main=paste0(PANEL_str,' (N=',N,')'),
            ylab=ylab,xlab=xlab,
            cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
            addScatterplot=T,las=1,notch=F)
    points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
           pch=16,cex=1.5,
           col=ys.colors_vec[rows_i])
    # if (j == 1) {
    #   boxp <- boxplot(ys.prediction[rows_i] ~ ys.colfactor[rows_i] + ys.group[rows_i],
    #                   main=paste0(PANEL_str,' (N=',N,')'),
    #                   ylab=ylab,xlab=xlab,
    #                   cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
    #                   col=ys.colors,las=1,notch=F,
    #                   at=setdiff(1:25, 5*1:5),boxwex=1.0,xaxt='n')
    #   axis(side=1,at=seq(2.5,by=5,length.out=5),
    #        labels=levels(ys.group),
    #        cex.axis=1/redf)
    #   #for normal plot: jitter(rep(1:ncol(boxp$stats), boxp$n))
    #   points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
    #          pch=16,cex=1.5,
    #          col='black')
    #   # points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
    #   #        pch=21,cex=1.5,
    #   #        col='blue',bg=ys.colors_vec[rows_i],lwd=3)
    # } else {
    #   boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
    #           main=paste0(PANEL_str,' (N=',N,')'),
    #           ylab=ylab,xlab=xlab,
    #           cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
    #           addScatterplot=T,las=1,notch=F)
    #   points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
    #          pch=16,cex=1.5,
    #          col=ys.colors_vec[rows_i])
    # }
    title(PVAL_str,outer=F,line=0.4,cex.main=1/redf)
    mtext(plab,adj=0,font=2,cex=1.4)
    # if (j == 1 & i == 4) {
    #   legend('right',inset=-0.30,
    #          legend=levels(ys.colfactor),
    #          col=ys.colors,
    #          pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
    # }
  }
  rm(ys.output)
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### (ALTERNATE) Extended Data Figure 6A + Extended Data Table 6A: Plotting all AxolotlN131 clocks applied to AxolotlN131 Demethylation (ALL SAMPLES) vs. Concentration
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains_list=list(c('Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells'),
                      c('Axolotl Demethylation Larva','Axolotl Demethylation Larva','Axolotl Demethylation Larva','Axolotl Demethylation Larva'))
y.axis.labs_list=list(c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail'),
                      c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail'))
x.axis.labs_list=list(c('Concentration','Concentration','Concentration','Concentration'),
                      c('Concentration','Concentration','Concentration','Concentration'))
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock"))
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock"))
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock"))
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock"))
input5.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=5, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva"))
input6.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=6, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva"))
input7.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=7, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva"))
input8.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=8, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva"))
#Creating 2 separate data frames, because "Concentration" has different factor levels in each experimental set
inputCell.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Prediction, PanelName, Experiment, Concentration, Passage, ExternalSampleID) %>%
  dplyr::mutate(Concentration = ifelse(is.na(Concentration), "None", Concentration))
inputCell.info_pred$PanelName <- factor(inputCell.info_pred$PanelName)
inputCell.info_pred$Experiment <- factor(inputCell.info_pred$Experiment)
inputCell.info_pred$Concentration <- relevel(factor(inputCell.info_pred$Concentration), ref="None")
inputCell.info_pred$Passage <- factor(inputCell.info_pred$Passage)
inputCell.info_pred$ExternalSampleID <- factor(inputCell.info_pred$ExternalSampleID)
inputLarva.info_pred <- bind_rows(input5.info_pred, input6.info_pred, input7.info_pred, input8.info_pred) %>%
  dplyr::select(Basename, Prediction, PanelName, Experiment, Concentration, Duration, ExternalSampleID) %>%
  dplyr::filter(Duration %in% c(12,14))
inputLarva.info_pred$PanelName <- factor(inputLarva.info_pred$PanelName)
inputLarva.info_pred$Experiment <- factor(inputLarva.info_pred$Experiment)
inputLarva.info_pred$Concentration <- factor(inputLarva.info_pred$Concentration)
inputLarva.info_pred$Duration <- factor(inputLarva.info_pred$Duration)
inputLarva.info_pred$ExternalSampleID <- factor(inputLarva.info_pred$ExternalSampleID)

inputCell.info_pred <- dplyr::arrange(inputCell.info_pred, PanelName, Experiment, Concentration, Passage, ExternalSampleID) #ensure correct order
inputLarva.info_pred <- dplyr::arrange(inputLarva.info_pred, PanelName, Experiment, Concentration, Duration, ExternalSampleID) #ensure correct order
GROUPVAR="Concentration"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC_list=list(c(NA,NA,NA,NA),
                     c(NA,NA,NA,NA))
ys.colors_original_list <- list(NA,NA)
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation-ALL_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Demethylation Experiment Data (ALL SAMPLES)'
ys.output_list <- list(inputCell.info_pred, inputLarva.info_pred)
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,4)
width=17
height=10
panel.labs_list=list(letters[1:4],
                     letters[5:8])#letters
oma.right=0
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylation-ALL.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred_list = list(length(ys.output_list))
for (j in 1:length(ys.output_list)) {
  input.info_pred <- ys.output_list[[j]]

  if (j == 1) {
    table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, Concentration, Prediction)
    table.input_pred <- table.input_pred %>%
      dplyr::group_by(Clock, Concentration) %>%
      dplyr::mutate(Num = row_number()) %>%
      tidyr::spread(Num, Prediction)
  }
  if (j == 2) {
    table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, Concentration, Prediction)
    table.input_pred <- table.input_pred %>%
      dplyr::group_by(Clock, Concentration) %>%
      dplyr::mutate(Num = row_number()) %>%
      tidyr::spread(Num, Prediction)
  }
  levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
  table.input_pred_list[[j]] <- table.input_pred

  rm(input.info_pred,table.input_pred)
}
write.xlsx(table.input_pred_list,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold"),
           sheetName=c("DemethylationTreatmentCells","DemethylationTreatmentLarva")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (j in 1:length(ys.output_list)) {
  #Extracting vectors from vector-lists, for each experiment separately
  ys.output <- ys.output_list[[j]]
  panel.mains <- panel.mains_list[[j]]
  y.axis.labs <- y.axis.labs_list[[j]]
  x.axis.labs <- x.axis.labs_list[[j]]
  COLVAR.VEC <- COLVAR.VEC_list[[j]]
  ys.colors_original <- ys.colors_original_list[[j]]
  panel.labs <- panel.labs_list[[j]]

  ys.group <- ys.output[,GROUPVAR]
  ys.prediction <- ys.output[,PREDVAR]
  ys.panelfactor <- ys.output[,PANELVAR]
  for (i in 1:length(levels(ys.panelfactor))) {
    COLVAR=COLVAR.VEC[i]
    ys.colors <- ys.colors_original
    if (!is.na(COLVAR)) {
      ys.colfactor <- ys.output[,COLVAR]
      # ys.colors contains the palette of distinct colors
      # ys.colors_vec contains the color assignment for every row in the data set
      if (is.na(ys.colors[1])) {
        if (length(levels(ys.colfactor)) <= 2){
          ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
        } else if (length(levels(ys.colfactor)) <= 8) {
          ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
        } else {
          ys.colors <- rainbow(length(levels(ys.colfactor)))
        }
      }
    } else {
      ys.colfactor <- rep(1, nrow(ys.output))
      ys.colors <- 'black'
    }
    ys.colors_vec <- ys.colors[ys.colfactor]

    rows_i <- which(as.numeric(ys.panelfactor)==i)
    PVAL <- kruskal.test(ys.prediction[rows_i], ys.group[rows_i])$p.value
    PVAL_str <- paste0("p=",signif(PVAL,2))
    PANEL_str <- panel.mains[i]
    N <- length(which(!is.na(ys.group[rows_i]) & !is.na(ys.prediction[rows_i])))
    ylab=y.axis.labs[i]
    xlab=x.axis.labs[i]
    plab=panel.labs[i]
    boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
            main=paste0(PANEL_str,' (N=',N,')'),
            ylab=ylab,xlab=xlab,
            cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
            addScatterplot=T,las=1,notch=F)
    points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
           pch=16,cex=1.5,
           col=ys.colors_vec[rows_i])
    # if (j == 1) {
    #   boxp <- boxplot(ys.prediction[rows_i] ~ ys.colfactor[rows_i] + ys.group[rows_i],
    #                   main=paste0(PANEL_str,' (N=',N,')'),
    #                   ylab=ylab,xlab=xlab,
    #                   cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
    #                   col=ys.colors,las=1,notch=F,
    #                   at=setdiff(1:25, 5*1:5),boxwex=1.0,xaxt='n')
    #   axis(side=1,at=seq(2.5,by=5,length.out=5),
    #        labels=levels(ys.group),
    #        cex.axis=1/redf)
    #   #for normal plot: jitter(rep(1:ncol(boxp$stats), boxp$n))
    #   points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
    #          pch=16,cex=1.5,
    #          col='black')
    #   # points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
    #   #        pch=21,cex=1.5,
    #   #        col='blue',bg=ys.colors_vec[rows_i],lwd=3)
    # } else {
    #   boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
    #           main=paste0(PANEL_str,' (N=',N,')'),
    #           ylab=ylab,xlab=xlab,
    #           cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
    #           addScatterplot=T,las=1,notch=F)
    #   points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
    #          pch=16,cex=1.5,
    #          col=ys.colors_vec[rows_i])
    # }
    title(PVAL_str,outer=F,line=0.4,cex.main=1/redf)
    mtext(plab,adj=0,font=2,cex=1.4)
    # if (j == 1 & i == 4) {
    #   legend('right',inset=-0.30,
    #          legend=levels(ys.colfactor),
    #          col=ys.colors,
    #          pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
    # }
  }
  rm(ys.output)
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Extended Data Figure 6B + Extended Data Table 6B: Plotting all AxolotlN131 clocks applied to AxolotlN131 Demethylation vs. Concentration
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Concentration','Concentration','Concentration','Concentration')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Prediction, PanelName, Experiment, Concentration, Passage, ExternalSampleID) %>%
  dplyr::mutate(Concentration = ifelse(is.na(Concentration), "None", Concentration))
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Experiment <- factor(input.info_pred$Experiment)
input.info_pred$Concentration <- relevel(factor(input.info_pred$Concentration), ref="None")
input.info_pred$Passage <- factor(input.info_pred$Passage)
input.info_pred$ExternalSampleID <- factor(input.info_pred$ExternalSampleID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, Experiment, Concentration, Passage, ExternalSampleID) #ensure correct order
GROUPVAR="Concentration"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Demethylation Experiment REVISION Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(1,4)
width=17
height=6
panel.labs=letters
oma.right=0
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, Concentration, Prediction)
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Clock, Concentration) %>%
  dplyr::mutate(Num = row_number()) %>%
  tidyr::spread(Num, Prediction)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.group <- ys.output[,GROUPVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PVAL <- kruskal.test(ys.prediction[rows_i], ys.group[rows_i])$p.value
  PVAL_str <- paste0("p=",signif(PVAL,2))
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.group[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
          main=paste0(PANEL_str,' (N=',N,')'),
          ylab=ylab,xlab=xlab,
          cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
          addScatterplot=T,las=1,notch=F)
  points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
         pch=16,cex=1.5,
         col=ys.colors_vec[rows_i])
  # if (j == 1) {
  #   boxp <- boxplot(ys.prediction[rows_i] ~ ys.colfactor[rows_i] + ys.group[rows_i],
  #                   main=paste0(PANEL_str,' (N=',N,')'),
  #                   ylab=ylab,xlab=xlab,
  #                   cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #                   col=ys.colors,las=1,notch=F,
  #                   at=setdiff(1:25, 5*1:5),boxwex=1.0,xaxt='n')
  #   axis(side=1,at=seq(2.5,by=5,length.out=5),
  #        labels=levels(ys.group),
  #        cex.axis=1/redf)
  #   #for normal plot: jitter(rep(1:ncol(boxp$stats), boxp$n))
  #   points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
  #          pch=16,cex=1.5,
  #          col='black')
  #   # points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
  #   #        pch=21,cex=1.5,
  #   #        col='blue',bg=ys.colors_vec[rows_i],lwd=3)
  # } else {
  #   boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
  #           main=paste0(PANEL_str,' (N=',N,')'),
  #           ylab=ylab,xlab=xlab,
  #           cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #           addScatterplot=T,las=1,notch=F)
  #   points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
  #          pch=16,cex=1.5,
  #          col=ys.colors_vec[rows_i])
  # }
  title(PVAL_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,adj=0,font=2,cex=1.4)
  # if (j == 1 & i == 4) {
  #   legend('right',inset=-0.30,
  #          legend=levels(ys.colfactor),
  #          col=ys.colors,
  #          pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  # }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION_PANEL"]] <- unique(input.info_pred$Basename)
#######################################
#######################################

###############################################################################
### (ALTERNATE) Extended Data Figure 6B + Extended Data Table 6B: Plotting all AxolotlN131 clocks applied to AxolotlN131 Demethylation (ALL SAMPLES) vs. Concentration
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Concentration','Concentration','Concentration','Concentration')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtDemethylation_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision"))
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision"))
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision"))
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision"))
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Prediction, PanelName, Experiment, Concentration, Passage, ExternalSampleID) %>%
  dplyr::mutate(Concentration = ifelse(is.na(Concentration), "None", Concentration))
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Experiment <- factor(input.info_pred$Experiment)
input.info_pred$Concentration <- relevel(factor(input.info_pred$Concentration), ref="None")
input.info_pred$Passage <- factor(input.info_pred$Passage)
input.info_pred$ExternalSampleID <- factor(input.info_pred$ExternalSampleID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, Experiment, Concentration, Passage, ExternalSampleID) #ensure correct order
GROUPVAR="Concentration"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION-ALL_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Demethylation Experiment REVISION Data (ALL SAMPLES)'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(1,4)
width=17
height=6
panel.labs=letters
oma.right=0
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtDemethylationREVISION-ALL.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, Concentration, Prediction)
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Clock, Concentration) %>%
  dplyr::mutate(Num = row_number()) %>%
  tidyr::spread(Num, Prediction)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.group <- ys.output[,GROUPVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PVAL <- kruskal.test(ys.prediction[rows_i], ys.group[rows_i])$p.value
  PVAL_str <- paste0("p=",signif(PVAL,2))
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.group[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
          main=paste0(PANEL_str,' (N=',N,')'),
          ylab=ylab,xlab=xlab,
          cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
          addScatterplot=T,las=1,notch=F)
  points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
         pch=16,cex=1.5,
         col=ys.colors_vec[rows_i])
  # if (j == 1) {
  #   boxp <- boxplot(ys.prediction[rows_i] ~ ys.colfactor[rows_i] + ys.group[rows_i],
  #                   main=paste0(PANEL_str,' (N=',N,')'),
  #                   ylab=ylab,xlab=xlab,
  #                   cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #                   col=ys.colors,las=1,notch=F,
  #                   at=setdiff(1:25, 5*1:5),boxwex=1.0,xaxt='n')
  #   axis(side=1,at=seq(2.5,by=5,length.out=5),
  #        labels=levels(ys.group),
  #        cex.axis=1/redf)
  #   #for normal plot: jitter(rep(1:ncol(boxp$stats), boxp$n))
  #   points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
  #          pch=16,cex=1.5,
  #          col='black')
  #   # points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
  #   #        pch=21,cex=1.5,
  #   #        col='blue',bg=ys.colors_vec[rows_i],lwd=3)
  # } else {
  #   boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
  #           main=paste0(PANEL_str,' (N=',N,')'),
  #           ylab=ylab,xlab=xlab,
  #           cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #           addScatterplot=T,las=1,notch=F)
  #   points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
  #          pch=16,cex=1.5,
  #          col=ys.colors_vec[rows_i])
  # }
  title(PVAL_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,adj=0,font=2,cex=1.4)
  # if (j == 1 & i == 4) {
  #   legend('right',inset=-0.30,
  #          legend=levels(ys.colfactor),
  #          col=ys.colors,
  #          pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  # }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Extended Data Figure 6A + Extended Data Table 6A: Plotting overall mean methylation in AxolotlN131 Demethylation vs. Concentration
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains_list=list(c('Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells'),
                      c('Axolotl Demethylation Larva','Axolotl Demethylation Larva','Axolotl Demethylation Larva'))
y.axis.labs_list=list(c('Overall Mean Methylation','High Mean Methylation','Low Mean Methylation'),
                      c('Overall Mean Methylation','High Mean Methylation','Low Mean Methylation'))
x.axis.labs_list=list(c('Concentration','Concentration','Concentration'),
                      c('Concentration','Concentration','Concentration'))
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylation_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationHigh_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationLow_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(MeanMethylation)) %>%
  dplyr::mutate(PanelName=1, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylation) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock")) %>% #never use "DemethylationTreatmentCells"
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(MeanMethylationHigh)) %>%
  dplyr::mutate(PanelName=2, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylationHigh) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(MeanMethylationLow)) %>%
  dplyr::mutate(PanelName=3, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylationLow) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(MeanMethylation)) %>%
  dplyr::mutate(PanelName=4, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Mean=MeanMethylation) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input5.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(MeanMethylationHigh)) %>%
  dplyr::mutate(PanelName=5, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Mean=MeanMethylationHigh) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input6.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(MeanMethylationLow)) %>%
  dplyr::mutate(PanelName=6, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Mean=MeanMethylationLow) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
#Creating 2 separate data frames, because "Concentration" has different factor levels in each experimental set
inputCell.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred) %>%
  dplyr::select(Basename, Mean, PanelName, Experiment, Concentration, Passage, ExternalSampleID) %>%
  dplyr::mutate(Concentration = ifelse(is.na(Concentration), "None", Concentration))
inputCell.info_pred$PanelName <- factor(inputCell.info_pred$PanelName)
inputCell.info_pred$Experiment <- factor(inputCell.info_pred$Experiment)
inputCell.info_pred$Concentration <- relevel(factor(inputCell.info_pred$Concentration), ref="None")
inputCell.info_pred$Passage <- factor(inputCell.info_pred$Passage)
inputCell.info_pred$ExternalSampleID <- factor(inputCell.info_pred$ExternalSampleID)
inputLarva.info_pred <- bind_rows(input4.info_pred, input5.info_pred, input6.info_pred) %>%
  dplyr::select(Basename, Mean, PanelName, Experiment, Concentration, Duration, ExternalSampleID) %>%
  dplyr::filter(Duration %in% c(12,14))
inputLarva.info_pred$PanelName <- factor(inputLarva.info_pred$PanelName)
inputLarva.info_pred$Experiment <- factor(inputLarva.info_pred$Experiment)
inputLarva.info_pred$Concentration <- factor(inputLarva.info_pred$Concentration)
inputLarva.info_pred$Duration <- factor(inputLarva.info_pred$Duration)
inputLarva.info_pred$ExternalSampleID <- factor(inputLarva.info_pred$ExternalSampleID)

inputCell.info_pred <- dplyr::arrange(inputCell.info_pred, PanelName, Experiment, Concentration, Passage, ExternalSampleID) #ensure correct order
inputLarva.info_pred <- dplyr::arrange(inputLarva.info_pred, PanelName, Experiment, Concentration, Duration, ExternalSampleID) #ensure correct order
GROUPVAR="Concentration"
PREDVAR="Mean"
PANELVAR="PanelName"
COLVAR.VEC_list=list(c(NA,NA,NA),
                     c(NA,NA,NA))
ys.colors_original_list <- list(NA,NA)
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation_PANEL.pdf'
out.pdf.title='Mean Methylation of Demethylation Experiment Data'
ys.output_list <- list(inputCell.info_pred, inputLarva.info_pred)
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,3)
width=13
height=10
panel.labs_list=list(letters[1:3],
                     letters[4:6])#letters
oma.right=0
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred_list = list(length(ys.output_list))
for (j in 1:length(ys.output_list)) {
  input.info_pred <- ys.output_list[[j]]

  if (j == 1) {
    table.input_pred <- input.info_pred %>% dplyr::select(Probes=PanelName, Concentration, Mean)
    table.input_pred <- table.input_pred %>%
      dplyr::group_by(Probes, Concentration) %>%
      dplyr::mutate(Num = row_number()) %>%
      tidyr::spread(Num, Mean)
  }
  if (j == 2) {
    table.input_pred <- input.info_pred %>% dplyr::select(Probes=PanelName, Concentration, Mean)
    table.input_pred <- table.input_pred %>%
      dplyr::group_by(Probes, Concentration) %>%
      dplyr::mutate(Num = row_number()) %>%
      tidyr::spread(Num, Mean)
  }
  levels(table.input_pred$Probes) <- c("Overall","High","Low")
  table.input_pred_list[[j]] <- table.input_pred

  rm(input.info_pred,table.input_pred)
}
write.xlsx(table.input_pred_list,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold"),
           sheetName=c("DemethylationTreatmentCells","DemethylationTreatmentLarva")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
par(mgp=c(3.5,0.8,0)) #further axis titles for decimal values
for (j in 1:length(ys.output_list)) {
  #Extracting vectors from vector-lists, for each experiment separately
  ys.output <- ys.output_list[[j]]
  panel.mains <- panel.mains_list[[j]]
  y.axis.labs <- y.axis.labs_list[[j]]
  x.axis.labs <- x.axis.labs_list[[j]]
  COLVAR.VEC <- COLVAR.VEC_list[[j]]
  ys.colors_original <- ys.colors_original_list[[j]]
  panel.labs <- panel.labs_list[[j]]

  ys.group <- ys.output[,GROUPVAR]
  ys.prediction <- ys.output[,PREDVAR]
  ys.panelfactor <- ys.output[,PANELVAR]
  for (i in 1:length(levels(ys.panelfactor))) {
    COLVAR=COLVAR.VEC[i]
    ys.colors <- ys.colors_original
    if (!is.na(COLVAR)) {
      ys.colfactor <- ys.output[,COLVAR]
      # ys.colors contains the palette of distinct colors
      # ys.colors_vec contains the color assignment for every row in the data set
      if (is.na(ys.colors[1])) {
        if (length(levels(ys.colfactor)) <= 2){
          ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
        } else if (length(levels(ys.colfactor)) <= 8) {
          ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
        } else {
          ys.colors <- rainbow(length(levels(ys.colfactor)))
        }
      }
    } else {
      ys.colfactor <- rep(1, nrow(ys.output))
      ys.colors <- 'black'
    }
    ys.colors_vec <- ys.colors[ys.colfactor]

    rows_i <- which(as.numeric(ys.panelfactor)==i)
    PVAL <- kruskal.test(ys.prediction[rows_i], ys.group[rows_i])$p.value
    PVAL_str <- paste0("p=",signif(PVAL,2))
    PANEL_str <- panel.mains[i]
    N <- length(which(!is.na(ys.group[rows_i]) & !is.na(ys.prediction[rows_i])))
    ylab=y.axis.labs[i]
    xlab=x.axis.labs[i]
    plab=panel.labs[i]
    boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
            main=paste0(PANEL_str,' (N=',N,')'),
            ylab=ylab,xlab=xlab,
            cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
            addScatterplot=T,las=1,notch=F)
    points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
           pch=16,cex=1.5,
           col=ys.colors_vec[rows_i])
    # if (j == 1) {
    #   boxp <- boxplot(ys.prediction[rows_i] ~ ys.colfactor[rows_i] + ys.group[rows_i],
    #                   main=paste0(PANEL_str,' (N=',N,')'),
    #                   ylab=ylab,xlab=xlab,
    #                   cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
    #                   col=ys.colors,las=1,notch=F,
    #                   at=setdiff(1:20, 5*1:4),boxwex=1.0,xaxt='n')
    #   axis(side=1,at=seq(2.5,by=5,length.out=4),
    #        labels=levels(ys.group),
    #        cex.axis=1/redf)
    #   #for normal plot: jitter(rep(1:ncol(boxp$stats), boxp$n))
    #   points(jitter(rep(setdiff(1:20, 5*1:4), boxp$n)),ys.prediction[rows_i],
    #          pch=16,cex=1.5,
    #          col='black')
    #   # points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
    #   #        pch=21,cex=1.5,
    #   #        col='blue',bg=ys.colors_vec[rows_i],lwd=3)
    # } else {
    #   boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
    #           main=paste0(PANEL_str,' (N=',N,')'),
    #           ylab=ylab,xlab=xlab,
    #           cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
    #           addScatterplot=T,las=1,notch=F)
    #   points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
    #          pch=16,cex=1.5,
    #          col=ys.colors_vec[rows_i])
    # }
    title(PVAL_str,outer=F,line=0.4,cex.main=1/redf)
    mtext(plab,adj=0,font=2,cex=1.4)
    # if (j == 1 & i == 3) {
    #   legend('right',inset=-0.30,
    #          legend=levels(ys.colfactor),
    #          col=ys.colors,
    #          pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
    # }
  }
  rm(ys.output)
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### (ALTERNATE) Extended Data Figure 6A + Extended Data Table 6A: Plotting overall mean methylation in AxolotlN131 Demethylation (ALL SAMPLES) vs. Concentration
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains_list=list(c('Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells'),
                      c('Axolotl Demethylation Larva','Axolotl Demethylation Larva','Axolotl Demethylation Larva'))
y.axis.labs_list=list(c('Overall Mean Methylation','High Mean Methylation','Low Mean Methylation'),
                      c('Overall Mean Methylation','High Mean Methylation','Low Mean Methylation'))
x.axis.labs_list=list(c('Concentration','Concentration','Concentration'),
                      c('Concentration','Concentration','Concentration'))
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylation_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationHigh_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationLow_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(MeanMethylation)) %>%
  dplyr::mutate(PanelName=1, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylation) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock")) #never use "DemethylationTreatmentCells"
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(MeanMethylationHigh)) %>%
  dplyr::mutate(PanelName=2, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylationHigh) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock"))
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(MeanMethylationLow)) %>%
  dplyr::mutate(PanelName=3, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylationLow) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Clock"))
input4.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(MeanMethylation)) %>%
  dplyr::mutate(PanelName=4, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Mean=MeanMethylation) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva"))
input5.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(MeanMethylationHigh)) %>%
  dplyr::mutate(PanelName=5, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Mean=MeanMethylationHigh) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva"))
input6.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(MeanMethylationLow)) %>%
  dplyr::mutate(PanelName=6, Concentration=Comparison4.DecitabineConcentration, Duration=Comparison5.DecitabineDurationTreatment.Days,
                Mean=MeanMethylationLow) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentLarva"))
#Creating 2 separate data frames, because "Concentration" has different factor levels in each experimental set
inputCell.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred) %>%
  dplyr::select(Basename, Mean, PanelName, Experiment, Concentration, Passage, ExternalSampleID) %>%
  dplyr::mutate(Concentration = ifelse(is.na(Concentration), "None", Concentration))
inputCell.info_pred$PanelName <- factor(inputCell.info_pred$PanelName)
inputCell.info_pred$Experiment <- factor(inputCell.info_pred$Experiment)
inputCell.info_pred$Concentration <- relevel(factor(inputCell.info_pred$Concentration), ref="None")
inputCell.info_pred$Passage <- factor(inputCell.info_pred$Passage)
inputCell.info_pred$ExternalSampleID <- factor(inputCell.info_pred$ExternalSampleID)
inputLarva.info_pred <- bind_rows(input4.info_pred, input5.info_pred, input6.info_pred) %>%
  dplyr::select(Basename, Mean, PanelName, Experiment, Concentration, Duration, ExternalSampleID) %>%
  dplyr::filter(Duration %in% c(12,14))
inputLarva.info_pred$PanelName <- factor(inputLarva.info_pred$PanelName)
inputLarva.info_pred$Experiment <- factor(inputLarva.info_pred$Experiment)
inputLarva.info_pred$Concentration <- factor(inputLarva.info_pred$Concentration)
inputLarva.info_pred$Duration <- factor(inputLarva.info_pred$Duration)
inputLarva.info_pred$ExternalSampleID <- factor(inputLarva.info_pred$ExternalSampleID)

inputCell.info_pred <- dplyr::arrange(inputCell.info_pred, PanelName, Experiment, Concentration, Passage, ExternalSampleID) #ensure correct order
inputLarva.info_pred <- dplyr::arrange(inputLarva.info_pred, PanelName, Experiment, Concentration, Duration, ExternalSampleID) #ensure correct order
GROUPVAR="Concentration"
PREDVAR="Mean"
PANELVAR="PanelName"
COLVAR.VEC_list=list(c(NA,NA,NA),
                     c(NA,NA,NA))
ys.colors_original_list <- list(NA,NA)
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation-ALL_PANEL.pdf'
out.pdf.title='Mean Methylation of Demethylation Experiment Data (ALL SAMPLES)'
ys.output_list <- list(inputCell.info_pred, inputLarva.info_pred)
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,3)
width=13
height=10
panel.labs_list=list(letters[1:3],
                     letters[4:6])#letters
oma.right=0
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylation-ALL.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred_list = list(length(ys.output_list))
for (j in 1:length(ys.output_list)) {
  input.info_pred <- ys.output_list[[j]]

  if (j == 1) {
    table.input_pred <- input.info_pred %>% dplyr::select(Probes=PanelName, Concentration, Mean)
    table.input_pred <- table.input_pred %>%
      dplyr::group_by(Probes, Concentration) %>%
      dplyr::mutate(Num = row_number()) %>%
      tidyr::spread(Num, Mean)
  }
  if (j == 2) {
    table.input_pred <- input.info_pred %>% dplyr::select(Probes=PanelName, Concentration, Mean)
    table.input_pred <- table.input_pred %>%
      dplyr::group_by(Probes, Concentration) %>%
      dplyr::mutate(Num = row_number()) %>%
      tidyr::spread(Num, Mean)
  }
  levels(table.input_pred$Probes) <- c("Overall","High","Low")
  table.input_pred_list[[j]] <- table.input_pred

  rm(input.info_pred,table.input_pred)
}
write.xlsx(table.input_pred_list,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold"),
           sheetName=c("DemethylationTreatmentCells","DemethylationTreatmentLarva")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
par(mgp=c(3.5,0.8,0)) #further axis titles for decimal values
for (j in 1:length(ys.output_list)) {
  #Extracting vectors from vector-lists, for each experiment separately
  ys.output <- ys.output_list[[j]]
  panel.mains <- panel.mains_list[[j]]
  y.axis.labs <- y.axis.labs_list[[j]]
  x.axis.labs <- x.axis.labs_list[[j]]
  COLVAR.VEC <- COLVAR.VEC_list[[j]]
  ys.colors_original <- ys.colors_original_list[[j]]
  panel.labs <- panel.labs_list[[j]]

  ys.group <- ys.output[,GROUPVAR]
  ys.prediction <- ys.output[,PREDVAR]
  ys.panelfactor <- ys.output[,PANELVAR]
  for (i in 1:length(levels(ys.panelfactor))) {
    COLVAR=COLVAR.VEC[i]
    ys.colors <- ys.colors_original
    if (!is.na(COLVAR)) {
      ys.colfactor <- ys.output[,COLVAR]
      # ys.colors contains the palette of distinct colors
      # ys.colors_vec contains the color assignment for every row in the data set
      if (is.na(ys.colors[1])) {
        if (length(levels(ys.colfactor)) <= 2){
          ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
        } else if (length(levels(ys.colfactor)) <= 8) {
          ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
        } else {
          ys.colors <- rainbow(length(levels(ys.colfactor)))
        }
      }
    } else {
      ys.colfactor <- rep(1, nrow(ys.output))
      ys.colors <- 'black'
    }
    ys.colors_vec <- ys.colors[ys.colfactor]

    rows_i <- which(as.numeric(ys.panelfactor)==i)
    PVAL <- kruskal.test(ys.prediction[rows_i], ys.group[rows_i])$p.value
    PVAL_str <- paste0("p=",signif(PVAL,2))
    PANEL_str <- panel.mains[i]
    N <- length(which(!is.na(ys.group[rows_i]) & !is.na(ys.prediction[rows_i])))
    ylab=y.axis.labs[i]
    xlab=x.axis.labs[i]
    plab=panel.labs[i]
    boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
            main=paste0(PANEL_str,' (N=',N,')'),
            ylab=ylab,xlab=xlab,
            cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
            addScatterplot=T,las=1,notch=F)
    points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
           pch=16,cex=1.5,
           col=ys.colors_vec[rows_i])
    # if (j == 1) {
    #   boxp <- boxplot(ys.prediction[rows_i] ~ ys.colfactor[rows_i] + ys.group[rows_i],
    #                   main=paste0(PANEL_str,' (N=',N,')'),
    #                   ylab=ylab,xlab=xlab,
    #                   cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
    #                   col=ys.colors,las=1,notch=F,
    #                   at=setdiff(1:20, 5*1:4),boxwex=1.0,xaxt='n')
    #   axis(side=1,at=seq(2.5,by=5,length.out=4),
    #        labels=levels(ys.group),
    #        cex.axis=1/redf)
    #   #for normal plot: jitter(rep(1:ncol(boxp$stats), boxp$n))
    #   points(jitter(rep(setdiff(1:20, 5*1:4), boxp$n)),ys.prediction[rows_i],
    #          pch=16,cex=1.5,
    #          col='black')
    #   # points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
    #   #        pch=21,cex=1.5,
    #   #        col='blue',bg=ys.colors_vec[rows_i],lwd=3)
    # } else {
    #   boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
    #           main=paste0(PANEL_str,' (N=',N,')'),
    #           ylab=ylab,xlab=xlab,
    #           cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
    #           addScatterplot=T,las=1,notch=F)
    #   points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
    #          pch=16,cex=1.5,
    #          col=ys.colors_vec[rows_i])
    # }
    title(PVAL_str,outer=F,line=0.4,cex.main=1/redf)
    mtext(plab,adj=0,font=2,cex=1.4)
    # if (j == 1 & i == 3) {
    #   legend('right',inset=-0.30,
    #          legend=levels(ys.colfactor),
    #          col=ys.colors,
    #          pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
    # }
  }
  rm(ys.output)
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Extended Data Figure 6B + Extended Data Table 6B: Plotting overall mean methylation in AxolotlN131 Demethylation vs. Concentration
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells')
y.axis.labs=c('Overall Mean Methylation','High Mean Methylation','Low Mean Methylation')
x.axis.labs=c('Concentration','Concentration','Concentration')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylation_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationHigh_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationLow_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(MeanMethylation)) %>%
  dplyr::mutate(PanelName=1, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylation) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(MeanMethylationHigh)) %>%
  dplyr::mutate(PanelName=2, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylationHigh) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(MeanMethylationLow)) %>%
  dplyr::mutate(PanelName=3, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylationLow) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred) %>%
  dplyr::select(Basename, Mean, PanelName, Experiment, Concentration, Passage, ExternalSampleID) %>%
  dplyr::mutate(Concentration = ifelse(is.na(Concentration), "None", Concentration))
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Experiment <- factor(input.info_pred$Experiment)
input.info_pred$Concentration <- relevel(factor(input.info_pred$Concentration), ref="None")
input.info_pred$Passage <- factor(input.info_pred$Passage)
input.info_pred$ExternalSampleID <- factor(input.info_pred$ExternalSampleID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, Experiment, Concentration, Passage, ExternalSampleID) #ensure correct order
GROUPVAR="Concentration"
PREDVAR="Mean"
PANELVAR="PanelName"
COLVAR.VEC=c(NA,NA,NA)
ys.colors_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION_PANEL.pdf'
out.pdf.title='Mean Methylation of Demethylation Experiment REVISION Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(1,3)
width=13
height=6
panel.labs=letters
oma.right=0
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Probes=PanelName, Concentration, Mean)
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Probes, Concentration) %>%
  dplyr::mutate(Num = row_number()) %>%
  tidyr::spread(Num, Mean)
levels(table.input_pred$Probes) <- c("Overall","High","Low")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.group <- ys.output[,GROUPVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
par(mgp=c(3.5,0.8,0)) #further axis titles for decimal values

for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PVAL <- kruskal.test(ys.prediction[rows_i], ys.group[rows_i])$p.value
  PVAL_str <- paste0("p=",signif(PVAL,2))
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.group[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
          main=paste0(PANEL_str,' (N=',N,')'),
          ylab=ylab,xlab=xlab,
          cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
          addScatterplot=T,las=1,notch=F)
  points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
         pch=16,cex=1.5,
         col=ys.colors_vec[rows_i])
  # if (j == 1) {
  #   boxp <- boxplot(ys.prediction[rows_i] ~ ys.colfactor[rows_i] + ys.group[rows_i],
  #                   main=paste0(PANEL_str,' (N=',N,')'),
  #                   ylab=ylab,xlab=xlab,
  #                   cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #                   col=ys.colors,las=1,notch=F,
  #                   at=setdiff(1:20, 5*1:4),boxwex=1.0,xaxt='n')
  #   axis(side=1,at=seq(2.5,by=5,length.out=4),
  #        labels=levels(ys.group),
  #        cex.axis=1/redf)
  #   #for normal plot: jitter(rep(1:ncol(boxp$stats), boxp$n))
  #   points(jitter(rep(setdiff(1:20, 5*1:4), boxp$n)),ys.prediction[rows_i],
  #          pch=16,cex=1.5,
  #          col='black')
  #   # points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
  #   #        pch=21,cex=1.5,
  #   #        col='blue',bg=ys.colors_vec[rows_i],lwd=3)
  # } else {
  #   boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
  #           main=paste0(PANEL_str,' (N=',N,')'),
  #           ylab=ylab,xlab=xlab,
  #           cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #           addScatterplot=T,las=1,notch=F)
  #   points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
  #          pch=16,cex=1.5,
  #          col=ys.colors_vec[rows_i])
  # }
  title(PVAL_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,adj=0,font=2,cex=1.4)
  # if (j == 1 & i == 3) {
  #   legend('right',inset=-0.30,
  #          legend=levels(ys.colfactor),
  #          col=ys.colors,
  #          pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  # }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
# ### FOR DATA DEPOSITION
# Basename.list_data_depo[["MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION_PANEL"]] <- unique(input.info_pred$Basename)
#######################################
#######################################

###############################################################################
### (ALTERNATE) Extended Data Figure 6B + Extended Data Table 6B: Plotting overall mean methylation in AxolotlN131 Demethylation (ALL SAMPLES) vs. Concentration
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Demethylation Cells','Axolotl Demethylation Cells','Axolotl Demethylation Cells')
y.axis.labs=c('Overall Mean Methylation','High Mean Methylation','Low Mean Methylation')
x.axis.labs=c('Concentration','Concentration','Concentration')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylation_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationHigh_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationLow_subCPGaxolotln131_AxolotlN131ExpmtDemethylation.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(MeanMethylation)) %>%
  dplyr::mutate(PanelName=1, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylation) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision"))
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(MeanMethylationHigh)) %>%
  dplyr::mutate(PanelName=2, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylationHigh) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision"))
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(MeanMethylationLow)) %>%
  dplyr::mutate(PanelName=3, Concentration=Comparison4.DecitabineConcentration, Passage=PassageNumber,
                Mean=MeanMethylationLow) %>%
  dplyr::filter(Experiment %in% c("DemethylationTreatmentCells_Revision"))
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred) %>%
  dplyr::select(Basename, Mean, PanelName, Experiment, Concentration, Passage, ExternalSampleID) %>%
  dplyr::mutate(Concentration = ifelse(is.na(Concentration), "None", Concentration))
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Experiment <- factor(input.info_pred$Experiment)
input.info_pred$Concentration <- relevel(factor(input.info_pred$Concentration), ref="None")
input.info_pred$Passage <- factor(input.info_pred$Passage)
input.info_pred$ExternalSampleID <- factor(input.info_pred$ExternalSampleID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, Experiment, Concentration, Passage, ExternalSampleID) #ensure correct order
GROUPVAR="Concentration"
PREDVAR="Mean"
PANELVAR="PanelName"
COLVAR.VEC=c(NA,NA,NA)
ys.colors_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION-ALL_PANEL.pdf'
out.pdf.title='Mean Methylation of Demethylation Experiment REVISION Data (ALL SAMPLES)'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(1,3)
width=13
height=6
panel.labs=letters
oma.right=0
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_ExpmtDemethylationREVISION-ALL.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Probes=PanelName, Concentration, Mean)
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Probes, Concentration) %>%
  dplyr::mutate(Num = row_number()) %>%
  tidyr::spread(Num, Mean)
levels(table.input_pred$Probes) <- c("Overall","High","Low")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.group <- ys.output[,GROUPVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
par(mgp=c(3.5,0.8,0)) #further axis titles for decimal values

for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PVAL <- kruskal.test(ys.prediction[rows_i], ys.group[rows_i])$p.value
  PVAL_str <- paste0("p=",signif(PVAL,2))
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.group[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
          main=paste0(PANEL_str,' (N=',N,')'),
          ylab=ylab,xlab=xlab,
          cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
          addScatterplot=T,las=1,notch=F)
  points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
         pch=16,cex=1.5,
         col=ys.colors_vec[rows_i])
  # if (j == 1) {
  #   boxp <- boxplot(ys.prediction[rows_i] ~ ys.colfactor[rows_i] + ys.group[rows_i],
  #                   main=paste0(PANEL_str,' (N=',N,')'),
  #                   ylab=ylab,xlab=xlab,
  #                   cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #                   col=ys.colors,las=1,notch=F,
  #                   at=setdiff(1:20, 5*1:4),boxwex=1.0,xaxt='n')
  #   axis(side=1,at=seq(2.5,by=5,length.out=4),
  #        labels=levels(ys.group),
  #        cex.axis=1/redf)
  #   #for normal plot: jitter(rep(1:ncol(boxp$stats), boxp$n))
  #   points(jitter(rep(setdiff(1:20, 5*1:4), boxp$n)),ys.prediction[rows_i],
  #          pch=16,cex=1.5,
  #          col='black')
  #   # points(jitter(rep(setdiff(1:25, 5*1:5), boxp$n)),ys.prediction[rows_i],
  #   #        pch=21,cex=1.5,
  #   #        col='blue',bg=ys.colors_vec[rows_i],lwd=3)
  # } else {
  #   boxplot(ys.prediction[rows_i] ~ ys.group[rows_i],
  #           main=paste0(PANEL_str,' (N=',N,')'),
  #           ylab=ylab,xlab=xlab,
  #           cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #           addScatterplot=T,las=1,notch=F)
  #   points(jitter(as.numeric(ys.group[rows_i])),ys.prediction[rows_i],
  #          pch=16,cex=1.5,
  #          col=ys.colors_vec[rows_i])
  # }
  title(PVAL_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,adj=0,font=2,cex=1.4)
  # if (j == 1 & i == 3) {
  #   legend('right',inset=-0.30,
  #          legend=levels(ys.colfactor),
  #          col=ys.colors,
  #          pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  # }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Extended Data Figure 12: Plotting Before/After-Two-Years pan-tissue and tissue-specific LOO clocks together
###############################################################################
panel.mains=c('Axolotl Before Two Years PanTissue','Axolotl Before Two Years LimbTail','Axolotl Before Two Years Limb','Axolotl Before Two Years Tail',
              'Axolotl After Two Years PanTissue','Axolotl After Two Years LimbTail','Axolotl After Two Years Limb','Axolotl After Two Years Tail')
y.axis.labs=c('DNAmAgeLOO Before Two Years Pan Tissue','DNAmAgeLOO Before Two Years LimbTail','DNAmAgeLOO Before Two Years Limb','DNAmAgeLOO Before Two Years Tail',
              'DNAmAgeLOO After Two Years Pan Tissue','DNAmAgeLOO After Two Years LimbTail','DNAmAgeLOO After Two Years Limb','DNAmAgeLOO After Two Years Tail')
x.axis.labs=c('Age','Age','Age','Age',
              'Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYears_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYearsLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYearsLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOBeforeTwoYearsTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output5.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOAfterTwoYears_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output6.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOAfterTwoYearsLimbTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
output7.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOAfterTwoYearsLimb_Final_subCPGcombinationmiddlefilter_EpigeneticLog2Age_PredictedValues.csv'
output8.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_AgeLOO_Final_Analysis/Subset_AxolotlN131_LOOAfterTwoYearsTail_Final_subCPGaxolotln131_EpigeneticLog2Age_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgeLOO)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgeLOO)
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgeLOO)
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgeLOO)
input5.info_pred <- dplyr::filter(read.csv(output5.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=5, Outcome=Age, Prediction=DNAmAgeLOO)
input6.info_pred <- dplyr::filter(read.csv(output6.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=6, Outcome=Age, Prediction=DNAmAgeLOO)
input7.info_pred <- dplyr::filter(read.csv(output7.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=7, Outcome=Age, Prediction=DNAmAgeLOO)
input8.info_pred <- dplyr::filter(read.csv(output8.csv, as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::mutate(PanelName=8, Outcome=Age, Prediction=DNAmAgeLOO)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred,
                             input5.info_pred, input6.info_pred, input7.info_pred, input8.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("Tissue","Tissue","Tissue","Tissue",
             "Tissue","Tissue","Tissue","Tissue")#c("SpeciesTissue","SpeciesTissue","SpeciesTissue","SpeciesTissue")
NUMVAR.VEC=c(NA,NA,NA,NA,
             NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/LOOBeforeAfterTwoYearsCOMBINED_AxolotlN131_Final_PANEL.pdf'
out.pdf.title='Leave-One-Out Analysis of Before/After-Two-Years Final Epigenetic Clocks'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,4)
width=17
height=10
panel.labs=letters
oma.right=0
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2) {
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  ## NEW, CUSTOM b/c reviewer request
  MeanAE <- mean(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MeanAE_str <- paste0("MeanAE=",signif(MeanAE,3))
  MAE_str <- paste0(MAE_str,', ',MeanAE_str)
  ##
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 1) {
    legend('topleft',
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
### FOR DATA DEPOSITION
Basename.list_data_depo[["LOOBeforeAfterTwoYearsCOMBINED_AxolotlN131_Final_PANEL"]] <- unique(input.info_pred$Basename)
#######################################
#######################################

###############################################################################
### Internal Figure: Plotting all AxolotlN131 clocks applied to AxolotlN131 Branch 2 EarlyLife
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl N131 Branch 2','Axolotl N131 Branch 2','Axolotl N131 Branch 2','Axolotl N131 Branch 2')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131Branch2_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131Branch2_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131Branch2_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131Branch2_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail)
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb)
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName) %>%
  dplyr::filter(Outcome <= 4.0)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("Tissue","Tissue","Tissue","Tissue")
NUMVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toAxolotlN131Branch2EarlyLife_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Axolotl N131 Branch 2 Early Life Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=9.7 #wider for exterior legend
height=10
panel.labs=letters
oma.right=5
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.30,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Internal Figure: Plotting all AxolotlN131 clocks applied to AxolotlN146 EarlyLife
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl N146 EarlyLife','Axolotl N146 EarlyLife','Axolotl N146 EarlyLife','Axolotl N146 EarlyLife')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN146_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, ExternalSampleID) %>%
  dplyr::filter(Outcome <= 4.0)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("Tissue","Tissue","Tissue","Tissue")
NUMVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toAxolotlN146EarlyLife_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Axolotl N146 Early Life Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=9.7 #wider for exterior legend
height=10
panel.labs=letters
oma.right=5
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.30,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Internal Figure + Internal Table: Plotting all AxolotlN131 clocks applied to AxolotlN131 TERT Tail Regen
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl TERT Tail Regen','Axolotl TERT Tail Regen','Axolotl TERT Tail Regen','Axolotl TERT Tail Regen')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTERTTailRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTERTTailRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTERTTailRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtTERTTailRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, Experiment, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup)
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
NUMVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTERTTailRegen_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to TERT Tail Regen Experiment Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=11.3 #wider for exterior legend
height=10
panel.labs=letters
oma.right=16
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtTERTTailRegen.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, AnimalID, RegenExperimentGroup, Prediction)
## handling replicates: giving distinct labels
table.input_pred <- table.input_pred %>% dplyr::mutate(RegenExperimentGroup = as.character(RegenExperimentGroup)) %>%
  dplyr::group_by(Clock, AnimalID, RegenExperimentGroup) %>%
  dplyr::mutate(RegenExperimentGroup = ifelse(row_number()==1, RegenExperimentGroup, paste(RegenExperimentGroup, row_number(), sep="-"))) %>%
  dplyr::ungroup()
##
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Clock, AnimalID) %>%
  tidyr::spread(RegenExperimentGroup, Prediction) %>%
  dplyr::arrange(Clock, AnimalID)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8,
         font=2)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.95,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Internal Figure: Plotting all AxolotlN131 clocks applied to AxolotlN131 Skin and Test Set
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Skin and Test Set','Axolotl Skin and Test Set','Axolotl Skin and Test Set')
y.axis.labs=c('DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131SkinAndTestSet_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131SkinAndTestSet_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131SkinAndTestSet_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Age >= 4.1 & Age <= 5.0) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Age >= 4.1 & Age <= 5.0) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Age >= 4.1 & Age <= 5.0) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, Experiment, AnimalName)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$Experiment <- factor(input.info_pred$Experiment)
input.info_pred$Group <- as.character(factor(input.info_pred$AnimalName, levels=c('EED','EZh')))
input.info_pred$Group[is.na(input.info_pred$Group)] <- 'Skin'
input.info_pred$Group <- factor(input.info_pred$Group, levels=c('Skin','EED','EZh'))
input.info_pred$AnimalName <- NULL

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, Experiment, Group) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("Group","Group","Group")
NUMVAR.VEC=c(NA,NA,NA)
ys.colors_original <- c('black','dodgerblue','darkorange')
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toSkinAndTestSet_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Skin and Test Set Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(1,3)
width=14.0 #wider for exterior legend
height=5
panel.labs=letters
oma.right=7
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
    COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
    abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 3) {
    legend('right',inset=-0.30,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Internal Figure + Internal Table: Plotting overall mean methylation in AxolotlN131 Skin and Test Set vs. Age
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Skin and Test Set','Axolotl Skin and Test Set','Axolotl Skin and Test Set')
y.axis.labs=c('Overall Mean Methylation','High Mean Methylation','Low Mean Methylation')
x.axis.labs=c('Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylation_subCPGaxolotln131_AxolotlN131SkinAndTestSet.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationHigh_subCPGaxolotln131_AxolotlN131SkinAndTestSet.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_MeanMethylation/Subset_AxolotlN131_MeanMethylationLow_subCPGaxolotln131_AxolotlN131SkinAndTestSet.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(MeanMethylation)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=MeanMethylation) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes") %>%
  dplyr::filter(ConfidenceInAgeEstimate >= 90)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(MeanMethylationHigh)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=MeanMethylationHigh) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes") %>%
  dplyr::filter(ConfidenceInAgeEstimate >= 90)
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(MeanMethylationLow)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=MeanMethylationLow) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes") %>%
  dplyr::filter(ConfidenceInAgeEstimate >= 90)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, Experiment, AnimalName)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$Experiment <- factor(input.info_pred$Experiment)
input.info_pred$Group <- as.character(factor(input.info_pred$AnimalName, levels=c('EED','EZh')))
input.info_pred$Group[is.na(input.info_pred$Group)] <- 'Skin'
input.info_pred$Group <- factor(input.info_pred$Group, levels=c('Skin','EED','EZh'))
input.info_pred$AnimalName <- NULL

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, Experiment, Group) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("Group","Group","Group")
NUMVAR.VEC=c(NA,NA,NA)
ys.colors_original <- c('black','dodgerblue','darkorange')
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/MeanMethylationCOMBINED_AxolotlN131_SkinAndTestSet_PANEL.pdf'
out.pdf.title='Mean Methylation of Skin and Test Set Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(1,3)
width=14.0 #wider for exterior legend
height=5
panel.labs=letters
oma.right=7
pointsize=12

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  # if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
  #   COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
  #              use='pairwise.complete.obs')
  #   COR_str <- paste0("cor=",signif(COR,2))
  # }
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  abline(h=0.5,lty="dashed")
  # title(MAE_str,outer=F,line=0.4,cex.main=1/redf)
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 3) {
    legend('right',inset=-0.30,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Internal Figure + Internal Table: Plotting all AxolotlN131 clocks applied to AxolotlN131 Persistence of Limb Regen
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Persistence of Limb Regen','Axolotl Persistence of Limb Regen','Axolotl Persistence of Limb Regen','Axolotl Persistence of Limb Regen')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup,
                                               levels=c("MatureLeftLimb","RegeneratedLeftLimbRound1",
                                                        "MatureRightLimb(AgedControlForOnceRegenerated)","MatureLeftLimb(AgedControl)"))
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
NUMVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Persistence of Limb Regen Experiment Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=13.0 #wider for exterior legend
height=10
panel.labs=letters
oma.right=28
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, AnimalID, RegenExperimentGroup, Prediction)
## handling replicates: giving distinct labels
table.input_pred <- table.input_pred %>% dplyr::mutate(RegenExperimentGroup = as.character(RegenExperimentGroup)) %>%
  dplyr::group_by(Clock, AnimalID, RegenExperimentGroup) %>%
  dplyr::mutate(RegenExperimentGroup = ifelse(row_number()==1, RegenExperimentGroup, paste(RegenExperimentGroup, row_number(), sep="-"))) %>%
  dplyr::ungroup()
##
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Clock, AnimalID) %>%
  tidyr::spread(RegenExperimentGroup, Prediction) %>%
  dplyr::arrange(Clock, AnimalID)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.65,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

########  WHISKER PLOTS (ALT)  ########
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(sub("(.*)(\\.)", "\\1-WHISKERPLOTS\\2", out.pdf),
    width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  
  df_box <- data.frame(y=ys.prediction[rows_i], g=ys.colfactor[rows_i], Age=ys.outcome[rows_i]) %>%
    dplyr::arrange(g, Age) %>% dplyr::group_by(g) %>%
    dplyr::summarise(x=mean(range(Age)),
                     median=median(y),
                     lower=median(y)-IQR(y),
                     upper=median(y)+IQR(y))
  plot(median~x,data=df_box,
       ylim=range(c(lower,upper),na.rm=T),
       xlim=axis_limits(x),
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab="Experiment Group",
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors,pch=19)
  arrows(df_box$x,df_box$lower,df_box$x,df_box$upper,
         length=0.05,angle=90,code=3,
         col=ys.colors)
  # plot(pch='',y=range(df_box$y),x=axis_limits(df_box$Age),
  #      main=paste0(PANEL_str,' (N=',N,')'),
  #      ylab=ylab,xlab="Experiment Group",
  #      cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf)
  # boxplot(add=T,y~as.numeric(g),data=df_box,at=c(0.4794521,0.9863014,0.9863014,1.2005),
  #         cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #         axes=F,boxwex=0.05,col=ys.colors)
  
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.65,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### (ALTERNATE) Internal Figure + Internal Table: Plotting all AxolotlN131 clocks applied to AxolotlN131 Persistence of Limb Regen (ALL SAMPLES)
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Persistence of Limb Regen','Axolotl Persistence of Limb Regen','Axolotl Persistence of Limb Regen','Axolotl Persistence of Limb Regen')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbPersistRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll)
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail)
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb)
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail)
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup,
                                               levels=c("MatureLeftLimb","RegeneratedLeftLimbRound1",
                                                        "MatureRightLimb(AgedControlForOnceRegenerated)","MatureLeftLimb(AgedControl)"))
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
NUMVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen-ALL_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Persistence of Limb Regen Experiment Data (ALL SAMPLES)'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=13.0 #wider for exterior legend
height=10
panel.labs=letters
oma.right=28
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbPersistRegen-ALL.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, AnimalID, RegenExperimentGroup, Prediction)
## handling replicates: giving distinct labels
table.input_pred <- table.input_pred %>% dplyr::mutate(RegenExperimentGroup = as.character(RegenExperimentGroup)) %>%
  dplyr::group_by(Clock, AnimalID, RegenExperimentGroup) %>%
  dplyr::mutate(RegenExperimentGroup = ifelse(row_number()==1, RegenExperimentGroup, paste(RegenExperimentGroup, row_number(), sep="-"))) %>%
  dplyr::ungroup()
##
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Clock, AnimalID) %>%
  tidyr::spread(RegenExperimentGroup, Prediction) %>%
  dplyr::arrange(Clock, AnimalID)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.65,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

########  WHISKER PLOTS (ALT)  ########
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(sub("(.*)(\\.)", "\\1-WHISKERPLOTS\\2", out.pdf),
    width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  
  df_box <- data.frame(y=ys.prediction[rows_i], g=ys.colfactor[rows_i], Age=ys.outcome[rows_i]) %>%
    dplyr::arrange(g, Age) %>% dplyr::group_by(g) %>%
    dplyr::summarise(x=mean(range(Age)),
                     median=median(y),
                     lower=median(y)-IQR(y),
                     upper=median(y)+IQR(y))
  plot(median~x,data=df_box,
       ylim=range(c(lower,upper),na.rm=T),
       xlim=axis_limits(x),
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab="Experiment Group",
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors,pch=19)
  arrows(df_box$x,df_box$lower,df_box$x,df_box$upper,
         length=0.05,angle=90,code=3,
         col=ys.colors)
  # plot(pch='',y=range(df_box$y),x=axis_limits(df_box$Age),
  #      main=paste0(PANEL_str,' (N=',N,')'),
  #      ylab=ylab,xlab="Experiment Group",
  #      cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf)
  # boxplot(add=T,y~as.numeric(g),data=df_box,at=c(0.4794521,0.9863014,0.9863014,1.2005),
  #         cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #         axes=F,boxwex=0.05,col=ys.colors)
  
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-1.65,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Internal Figure + Internal Table: Plotting all AxolotlN131 clocks applied to AxolotlN131 Limb Regen
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Limb Regen','Axolotl Limb Regen','Axolotl Limb Regen','Axolotl Limb Regen',
              #'Axolotl Limb Regen Supp.','Axolotl Limb Regen Supp.','Axolotl Limb Regen Supp.','Axolotl Limb Regen Supp.',
              'Axolotl Limb Regen Supp. 2','Axolotl Limb Regen Supp. 2','Axolotl Limb Regen Supp. 2','Axolotl Limb Regen Supp. 2')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail',
              #'DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail',
              'DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age',
              #'Age','Age','Age','Age',
              'Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("LimbRegeneration")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegeneration")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("LimbRegeneration")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegeneration")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
# input5.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
#   dplyr::mutate(PanelName=5, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
#   dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary")) %>%
#   dplyr::filter(CanBeUsedForAgingStudies == "yes")
# input6.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
#   dplyr::mutate(PanelName=6, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
#   dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary")) %>%
#   dplyr::filter(CanBeUsedForAgingStudies == "yes")
# input7.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
#   dplyr::mutate(PanelName=7, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
#   dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary")) %>%
#   dplyr::filter(CanBeUsedForAgingStudies == "yes")
# input8.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
#   dplyr::mutate(PanelName=8, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
#   dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary")) %>%
#   dplyr::filter(CanBeUsedForAgingStudies == "yes")
input9.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=9, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary2")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input10.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=10, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary2")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input11.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=11, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary2")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input12.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=12, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary2")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred,
                             #input5.info_pred, input6.info_pred, input7.info_pred, input8.info_pred,
                             input9.info_pred, input10.info_pred, input11.info_pred, input12.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup,
                                               levels=c("OriginalLimb",
                                                        "EarlyBlastema","MidBlastema",
                                                        "StumpMidBlastema","PaletteBlastema",
                                                        "StumpPaletteBlastema"))
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup",
             #"RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup",
             "RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
NUMVAR.VEC=c(NA,NA,NA,NA,
             #NA,NA,NA,NA,
             NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Limb Regen Experiment Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,4)
width=19.6 #wider for exterior legend
height=10
panel.labs=letters
oma.right=18
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
input.info_pred.list <- list(dplyr::filter(input.info_pred, as.numeric(PanelName) <= 4) %>%
                               dplyr::mutate(PanelName = factor(PanelName)),
                             dplyr::filter(input.info_pred, as.numeric(PanelName) >= 5 & as.numeric(PanelName) <= 8) %>%
                               dplyr::mutate(PanelName = factor(PanelName)))#,
#dplyr::filter(input.info_pred, as.numeric(PanelName) >= 9) %>%
#dplyr::mutate(PanelName = factor(PanelName)))
table.input_pred.list <- list(data.frame(), data.frame())
for (j in 1:length(input.info_pred.list)) {
  input.info_pred_j <- input.info_pred.list[[j]]
  table.input_pred <- input.info_pred_j %>% dplyr::select(Clock=PanelName, AnimalID, RegenExperimentGroup, Prediction)
  ## handling replicates: giving distinct labels
  table.input_pred <- table.input_pred %>% dplyr::mutate(RegenExperimentGroup = as.character(RegenExperimentGroup)) %>%
    dplyr::group_by(Clock, AnimalID, RegenExperimentGroup) %>%
    dplyr::mutate(RegenExperimentGroup = ifelse(row_number()==1, RegenExperimentGroup, paste(RegenExperimentGroup, row_number(), sep="-"))) %>%
    dplyr::ungroup()
  ##
  table.input_pred <- table.input_pred %>%
    dplyr::group_by(Clock, AnimalID) %>%
    tidyr::spread(RegenExperimentGroup, Prediction) %>%
    dplyr::arrange(Clock, AnimalID)
  levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
  table.input_pred.list[[j]] <- table.input_pred
}
rm(input.info_pred_j, input.info_pred.list)
write.xlsx(table.input_pred.list,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold"),
           #sheetName=c("LimbRegeneration","LimbRegenerationSupplementary","LimbRegenerationSupplementary2")) #?buildWorkbook
           sheetName=c("LimbRegeneration","LimbRegenerationSupplementary2")) #?buildWorkbook
rm(table.input_pred.list)

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.75,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

########  WHISKER PLOTS (ALT)  ########
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(sub("(.*)(\\.)", "\\1-WHISKERPLOTS\\2", out.pdf),
    width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  
  df_box <- data.frame(y=ys.prediction[rows_i], g=ys.colfactor[rows_i], Age=ys.outcome[rows_i]) %>%
    dplyr::arrange(g, Age) %>% dplyr::group_by(g) %>%
    dplyr::summarise(x=mean(range(Age)),
                     median=median(y),
                     lower=median(y)-IQR(y),
                     upper=median(y)+IQR(y))
  plot(median~x,data=df_box,
       ylim=range(c(lower,upper),na.rm=T),
       xlim=axis_limits(x),
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab="Experiment Group",
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors,pch=19)
  arrows(df_box$x,df_box$lower,df_box$x,df_box$upper,
         length=0.05,angle=90,code=3,
         col=ys.colors)
  # plot(pch='',y=range(df_box$y),x=axis_limits(df_box$Age),
  #      main=paste0(PANEL_str,' (N=',N,')'),
  #      ylab=ylab,xlab="Experiment Group",
  #      cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf)
  # boxplot(add=T,y~as.numeric(g),data=df_box,at=c(0.4794521,0.9863014,0.9863014,1.2005),
  #         cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #         axes=F,boxwex=0.05,col=ys.colors)
  
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.75,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### (ALTERNATE) Internal Figure + Internal Table: Plotting all AxolotlN131 clocks applied to AxolotlN131 Limb Regen (ALL SAMPLES)
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Limb Regen','Axolotl Limb Regen','Axolotl Limb Regen','Axolotl Limb Regen',
              #'Axolotl Limb Regen Supp.','Axolotl Limb Regen Supp.','Axolotl Limb Regen Supp.','Axolotl Limb Regen Supp.',
              'Axolotl Limb Regen Supp. 2','Axolotl Limb Regen Supp. 2','Axolotl Limb Regen Supp. 2','Axolotl Limb Regen Supp. 2')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail',
              #'DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail',
              'DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age',
              #'Age','Age','Age','Age',
              'Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("LimbRegeneration"))
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegeneration"))
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("LimbRegeneration"))
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegeneration"))
# input5.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
#   dplyr::mutate(PanelName=5, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
#   dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary"))
# input6.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
#   dplyr::mutate(PanelName=6, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
#   dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary"))
# input7.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
#   dplyr::mutate(PanelName=7, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
#   dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary"))
# input8.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
#   dplyr::mutate(PanelName=8, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
#   dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary"))
input9.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=9, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary2"))
input10.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=10, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary2"))
input11.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=11, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary2"))
input12.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=12, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationSupplementary2"))
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred,
                             #input5.info_pred, input6.info_pred, input7.info_pred, input8.info_pred,
                             input9.info_pred, input10.info_pred, input11.info_pred, input12.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup,
                                               levels=c("OriginalLimb",
                                                        "EarlyBlastema","MidBlastema",
                                                        "StumpMidBlastema","PaletteBlastema",
                                                        "StumpPaletteBlastema"))
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup",
             #"RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup",
             "RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
NUMVAR.VEC=c(NA,NA,NA,NA,
             #NA,NA,NA,NA,
             NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen-ALL_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Limb Regen Experiment Data (ALL SAMPLES)'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,4)
width=19.6 #wider for exterior legend
height=10
panel.labs=letters
oma.right=18
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegen-ALL.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
input.info_pred.list <- list(dplyr::filter(input.info_pred, as.numeric(PanelName) <= 4) %>%
                               dplyr::mutate(PanelName = factor(PanelName)),
                             dplyr::filter(input.info_pred, as.numeric(PanelName) >= 5 & as.numeric(PanelName) <= 8) %>%
                               dplyr::mutate(PanelName = factor(PanelName)))#,
#dplyr::filter(input.info_pred, as.numeric(PanelName) >= 9) %>%
#dplyr::mutate(PanelName = factor(PanelName)))
table.input_pred.list <- list(data.frame(), data.frame())
for (j in 1:length(input.info_pred.list)) {
  input.info_pred_j <- input.info_pred.list[[j]]
  table.input_pred <- input.info_pred_j %>% dplyr::select(Clock=PanelName, AnimalID, RegenExperimentGroup, Prediction)
  ## handling replicates: giving distinct labels
  table.input_pred <- table.input_pred %>% dplyr::mutate(RegenExperimentGroup = as.character(RegenExperimentGroup)) %>%
    dplyr::group_by(Clock, AnimalID, RegenExperimentGroup) %>%
    dplyr::mutate(RegenExperimentGroup = ifelse(row_number()==1, RegenExperimentGroup, paste(RegenExperimentGroup, row_number(), sep="-"))) %>%
    dplyr::ungroup()
  ##
  table.input_pred <- table.input_pred %>%
    dplyr::group_by(Clock, AnimalID) %>%
    tidyr::spread(RegenExperimentGroup, Prediction) %>%
    dplyr::arrange(Clock, AnimalID)
  levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
  table.input_pred.list[[j]] <- table.input_pred
}
rm(input.info_pred_j, input.info_pred.list)
write.xlsx(table.input_pred.list,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold"),
           #sheetName=c("LimbRegeneration","LimbRegenerationSupplementary","LimbRegenerationSupplementary2")) #?buildWorkbook
           sheetName=c("LimbRegeneration","LimbRegenerationSupplementary2")) #?buildWorkbook
rm(table.input_pred.list)

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.75,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

########  WHISKER PLOTS (ALT)  ########
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(sub("(.*)(\\.)", "\\1-WHISKERPLOTS\\2", out.pdf),
    width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  
  df_box <- data.frame(y=ys.prediction[rows_i], g=ys.colfactor[rows_i], Age=ys.outcome[rows_i]) %>%
    dplyr::arrange(g, Age) %>% dplyr::group_by(g) %>%
    dplyr::summarise(x=mean(range(Age)),
                     median=median(y),
                     lower=median(y)-IQR(y),
                     upper=median(y)+IQR(y))
  plot(median~x,data=df_box,
       ylim=range(c(lower,upper),na.rm=T),
       xlim=axis_limits(x),
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab="Experiment Group",
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors,pch=19)
  arrows(df_box$x,df_box$lower,df_box$x,df_box$upper,
         length=0.05,angle=90,code=3,
         col=ys.colors)
  # plot(pch='',y=range(df_box$y),x=axis_limits(df_box$Age),
  #      main=paste0(PANEL_str,' (N=',N,')'),
  #      ylab=ylab,xlab="Experiment Group",
  #      cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf)
  # boxplot(add=T,y~as.numeric(g),data=df_box,at=c(0.4794521,0.9863014,0.9863014,1.2005),
  #         cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #         axes=F,boxwex=0.05,col=ys.colors)
  
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.75,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### Internal Figure + Internal Table: Plotting all AxolotlN131 clocks applied to AxolotlN131 Limb Regen Expanded
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Limb Regen Expanded','Axolotl Limb Regen Expanded','Axolotl Limb Regen Expanded','Axolotl Limb Regen Expanded')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationExpanded")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationExpanded")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationExpanded")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationExpanded")) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup,
                                               levels=c("OriginalLimb",
                                                        "EarlyBlastema","MidBlastema",
                                                        "StumpMidBlastema","PaletteBlastema",
                                                        "StumpPaletteBlastema"))
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
NUMVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Limb Regen Expanded Experiment Data'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=10.9 #wider for exterior legend
height=10
panel.labs=letters
oma.right=13
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, AnimalID, RegenExperimentGroup, Prediction)
## handling replicates: giving distinct labels
table.input_pred <- table.input_pred %>% dplyr::mutate(RegenExperimentGroup = as.character(RegenExperimentGroup)) %>%
  dplyr::group_by(Clock, AnimalID, RegenExperimentGroup) %>%
  dplyr::mutate(RegenExperimentGroup = ifelse(row_number()==1, RegenExperimentGroup, paste(RegenExperimentGroup, row_number(), sep="-"))) %>%
  dplyr::ungroup()
##
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Clock, AnimalID) %>%
  tidyr::spread(RegenExperimentGroup, Prediction) %>%
  dplyr::arrange(Clock, AnimalID)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.75,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

########  WHISKER PLOTS (ALT)  ########
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(sub("(.*)(\\.)", "\\1-WHISKERPLOTS\\2", out.pdf),
    width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  
  df_box <- data.frame(y=ys.prediction[rows_i], g=ys.colfactor[rows_i], Age=ys.outcome[rows_i]) %>%
    dplyr::arrange(g, Age) %>% dplyr::group_by(g) %>%
    dplyr::summarise(x=mean(range(Age)),
                     median=median(y),
                     lower=median(y)-IQR(y),
                     upper=median(y)+IQR(y))
  plot(median~x,data=df_box,
       ylim=range(c(lower,upper),na.rm=T),
       xlim=axis_limits(x),
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab="Experiment Group",
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors,pch=19)
  arrows(df_box$x,df_box$lower,df_box$x,df_box$upper,
         length=0.05,angle=90,code=3,
         col=ys.colors)
  # plot(pch='',y=range(df_box$y),x=axis_limits(df_box$Age),
  #      main=paste0(PANEL_str,' (N=',N,')'),
  #      ylab=ylab,xlab="Experiment Group",
  #      cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf)
  # boxplot(add=T,y~as.numeric(g),data=df_box,at=c(0.4794521,0.9863014,0.9863014,1.2005),
  #         cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #         axes=F,boxwex=0.05,col=ys.colors)
  
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.75,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### (ALTERNATE) Internal Figure + Internal Table: Plotting all AxolotlN131 clocks applied to AxolotlN131 Limb Regen Expanded (ALL SAMPLES)
###############################################################################
rm(list=ls(pattern="input"))
rm(ys.colors,ys.output,lim,l_lim,u_lim)
panel.mains=c('Axolotl Limb Regen Expanded','Axolotl Limb Regen Expanded','Axolotl Limb Regen Expanded','Axolotl Limb Regen Expanded')
y.axis.labs=c('DNAmAge Early Life Pan Tissue','DNAmAge Early Life LimbTail','DNAmAge Early Life Limb','DNAmAge Early Life Tail')
x.axis.labs=c('Age','Age','Age','Age')
output1.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLife_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output2.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimbTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output3.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeLimb_subCPGcombinationmiddlefilter_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
output4.csv='SpeciesSubsetAnalyses/Subset_AxolotlN131_Clocks_Final/Subset_AxolotlN131_ClockEarlyLifeTail_subCPGaxolotln131_basedOnAll_EpigeneticLog2Age_toAxolotlN131ExpmtLimbRegen_PredictedValues.csv'
input1.info_pred <- dplyr::filter(read.csv(output1.csv, as.is=T), !is.na(DNAmAgebasedOnAll)) %>%
  dplyr::mutate(PanelName=1, Outcome=Age, Prediction=DNAmAgebasedOnAll) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationExpanded"))
input2.info_pred <- dplyr::filter(read.csv(output2.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimbTail)) %>%
  dplyr::mutate(PanelName=2, Outcome=Age, Prediction=DNAmAgebasedOnAllLimbTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationExpanded"))
input3.info_pred <- dplyr::filter(read.csv(output3.csv, as.is=T), !is.na(DNAmAgebasedOnAllLimb)) %>%
  dplyr::mutate(PanelName=3, Outcome=Age, Prediction=DNAmAgebasedOnAllLimb) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationExpanded"))
input4.info_pred <- dplyr::filter(read.csv(output4.csv, as.is=T), !is.na(DNAmAgebasedOnAllTail)) %>%
  dplyr::mutate(PanelName=4, Outcome=Age, Prediction=DNAmAgebasedOnAllTail) %>%
  dplyr::filter(Experiment %in% c("LimbRegenerationExpanded"))
input.info_pred <- bind_rows(input1.info_pred, input2.info_pred, input3.info_pred, input4.info_pred) %>%
  dplyr::select(Basename, Outcome, Prediction, PanelName, Tissue, SpeciesLatinName, RegenExperimentGroup, AnimalID)
input.info_pred$PanelName <- factor(input.info_pred$PanelName)
input.info_pred$Tissue <- factor(input.info_pred$Tissue)
input.info_pred$SpeciesLatinName <- factor(input.info_pred$SpeciesLatinName,
                                           levels=c('Ambystoma mexicanum'))#,'Homo sapiens'))
input.info_pred$SpeciesTissue <- factor(paste.species_tissue(input.info_pred$SpeciesLatinName, input.info_pred$Tissue))
input.info_pred$RegenExperimentGroup <- factor(input.info_pred$RegenExperimentGroup,
                                               levels=c("OriginalLimb",
                                                        "EarlyBlastema","MidBlastema",
                                                        "StumpMidBlastema","PaletteBlastema",
                                                        "StumpPaletteBlastema"))
input.info_pred$AnimalID <- factor(input.info_pred$AnimalID)

input.info_pred <- dplyr::arrange(input.info_pred, PanelName, RegenExperimentGroup, AnimalID) #ensure correct order

OUTVAR="Outcome"
PREDVAR="Prediction"
PANELVAR="PanelName"
COLVAR.VEC=c("RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup","RegenExperimentGroup")
NUMVAR.VEC=c(NA,NA,NA,NA)
ys.colors_original <- NA
ys.numbers_original <- NA
out.pdf='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded-ALL_PANEL.pdf'
out.pdf.title='Final Epigenetic Axolotl Clocks applied to Limb Regen Expanded Experiment Data (ALL SAMPLES)'
ys.output <- input.info_pred
TITLE_str=paste0(out.pdf.title,'\n')
mfrow=c(2,2)
width=10.9 #wider for exterior legend
height=10
panel.labs=letters
oma.right=13
pointsize=12

out.xlsx='SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/ClockCOMBINED_AxolotlN131_Final_toExpmtLimbRegenExpanded-ALL.xlsx'

#######################################
#######################################
library(data.table)
library(openxlsx)
table.input_pred <- input.info_pred %>% dplyr::select(Clock=PanelName, AnimalID, RegenExperimentGroup, Prediction)
## handling replicates: giving distinct labels
table.input_pred <- table.input_pred %>% dplyr::mutate(RegenExperimentGroup = as.character(RegenExperimentGroup)) %>%
  dplyr::group_by(Clock, AnimalID, RegenExperimentGroup) %>%
  dplyr::mutate(RegenExperimentGroup = ifelse(row_number()==1, RegenExperimentGroup, paste(RegenExperimentGroup, row_number(), sep="-"))) %>%
  dplyr::ungroup()
##
table.input_pred <- table.input_pred %>%
  dplyr::group_by(Clock, AnimalID) %>%
  tidyr::spread(RegenExperimentGroup, Prediction) %>%
  dplyr::arrange(Clock, AnimalID)
levels(table.input_pred$Clock) <- c("Pan Tissue","LimbTail","Limb","Tail")
write.xlsx(table.input_pred,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#######################################
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(out.pdf,width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  NUMVAR=NUMVAR.VEC[i]
  ys.colors <- ys.colors_original
  ys.numbers <- ys.numbers_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.na(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.na(ys.numbers[1])) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  # lim <- axis_square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
  # l_lim=lim[1]
  # u_lim=lim[2]
  #DNAmAge on Expmt Samples is highly variable, so we free the axes
  xlim <- axis_limits(ys.outcome[rows_i])#axis_limits(ys.outcome)
  ylim <- axis_limits(ys.prediction[rows_i])#axis_limits(ys.prediction)
  l_lim=xlim[1]
  u_lim=xlim[2]
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
       type=type,main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab=xlab,pch=16,
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors_vec[rows_i],xlim=xlim,ylim=ylim)#xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         col=ys.colors_vec[rows_i],
         labels=ys.numbers_vec[rows_i],cex=0.8)
  }
  #Creating spaghetti plots instead of scatter plots
  for (animal in levels(ys.output[rows_i,"AnimalID"])) {
    rows_i_animal <- which(ys.output[rows_i,"AnimalID"]==animal)
    lines(y=ys.prediction[rows_i][rows_i_animal],x=ys.outcome[rows_i][rows_i_animal],
          lty='dashed',lwd=0.4)
  }
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.75,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

########  WHISKER PLOTS (ALT)  ########
#######################################
require(RColorBrewer)
ys.outcome <- ys.output[,OUTVAR]
ys.prediction <- ys.output[,PREDVAR]
ys.panelfactor <- ys.output[,PANELVAR]

pdf(sub("(.*)(\\.)", "\\1-WHISKERPLOTS\\2", out.pdf),
    width=width,height=height,pointsize=pointsize)#,units='in',res=600,pointsize=pointsize)
par(mfrow=mfrow)
redf=reduct_factor(mfrow)
par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
for (i in 1:length(levels(ys.panelfactor))) {
  COLVAR=COLVAR.VEC[i]
  ys.colors <- ys.colors_original
  if (!is.na(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.na(ys.colors[1])) {
      if (length(levels(ys.colfactor)) <= 2){
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  
  rows_i <- which(as.numeric(ys.panelfactor)==i)
  PANEL_str <- panel.mains[i]
  N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
  ylab=y.axis.labs[i]
  xlab=x.axis.labs[i]
  plab=panel.labs[i]
  
  df_box <- data.frame(y=ys.prediction[rows_i], g=ys.colfactor[rows_i], Age=ys.outcome[rows_i]) %>%
    dplyr::arrange(g, Age) %>% dplyr::group_by(g) %>%
    dplyr::summarise(x=mean(range(Age)),
                     median=median(y),
                     lower=median(y)-IQR(y),
                     upper=median(y)+IQR(y))
  plot(median~x,data=df_box,
       ylim=range(c(lower,upper),na.rm=T),
       xlim=axis_limits(x),
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=ylab,xlab="Experiment Group",
       cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
       col=ys.colors,pch=19)
  arrows(df_box$x,df_box$lower,df_box$x,df_box$upper,
         length=0.05,angle=90,code=3,
         col=ys.colors)
  # plot(pch='',y=range(df_box$y),x=axis_limits(df_box$Age),
  #      main=paste0(PANEL_str,' (N=',N,')'),
  #      ylab=ylab,xlab="Experiment Group",
  #      cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf)
  # boxplot(add=T,y~as.numeric(g),data=df_box,at=c(0.4794521,0.9863014,0.9863014,1.2005),
  #         cex.main=1/redf,cex.lab=1/redf,cex.axis=1/redf,
  #         axes=F,boxwex=0.05,col=ys.colors)
  
  mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  if (i == 4) {
    legend('right',inset=-0.75,
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1/redf,pt.cex=2/redf,xpd=NA)
  }
}
title(TITLE_str,outer=T,line=-1,cex.main=1/redf)
dev.off()
#######################################
#######################################

###############################################################################
### DATA DEPOSITION SAMPLE SHEET
## Creating merged sample info sheet with Figure tracking column
#convert list of character vectors into single data frame
info_data_depo <- do.call(rbind, lapply(names(Basename.list_data_depo), function(Figure) {
  Basenames <- Basename.list_data_depo[[Figure]]
  if (length(Basenames) == 0) return(NULL)  # Skip empty vectors
  data.frame(Figure = Figure, Basename = Basenames, stringsAsFactors = FALSE)
}))
#sapply(Basename.list_data_depo, length) ## CHECK: checking samples sizes with manuscript
#sort(table(info_data_depo$Basename))
info_data_depo <- info_data_depo %>% group_by(Basename) %>%
  summarise(Figure = paste(Figure, collapse = " & "), .groups = "drop") %>% ungroup() %>%
  dplyr::select(Figure, Basename)

#load sample sheet data
infoAllSamp.csv=c('N131.ET0087.SalamanderMaxYun/SampleSheetAgeN131final.csv')
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

infoHuman.csv=c('SpeciesSubsetAnalyses/HumanAxolotlN131_AgeLOFOBalance_Final_Analysis/HumanAxolotlN131_LOFO10BalanceEarlyLife_Final_subCPGaxolotln131_EpigeneticLLin3Age_PredictedValues.csv')
infoHuman=read.csv(infoHuman.csv, as.is=T) %>% dplyr::filter(!SpeciesLatinName %in% c("Ambystoma mexicanum"))
infoHuman <- infoHuman %>%
  dplyr::select(Basename,SpeciesLatinName,OriginalOrderInBatch,Age,ConfidenceInAgeEstimate,
                CanBeUsedForAgingStudies,Tissue,Female,SpeciesCommonName,ExternalSampleID,Folder,
                Experiment,RegenExperimentGroup,AnimalID,AnimalName)
rm(infoHuman.csv)
if ("Female" %in% colnames(infoHuman)) {
  infoHuman$Female[which(is.na(infoHuman$Female))] <- "NA"
  infoHuman$Female <- factor(infoHuman$Female, levels=c(0,1,"NA"))
  levels(infoHuman$Female) <- c("Male","Female","NA")
}

infoFrogs.csv=c('SpeciesSubsetAnalyses/Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis/Subset_AxolotlClawedFrogN131N140_LOOEarlyLife_Final_subCPGaxolotlclawedfrogn131n140_EpigeneticLog2Age_PredictedValues.csv')
infoFrogs=read.csv(infoFrogs.csv, as.is=T) %>% dplyr::filter(!SpeciesLatinName %in% c("Ambystoma mexicanum"))
infoFrogs <- infoFrogs %>%
  dplyr::select(Basename,SpeciesLatinName,OriginalOrderInBatch,Age,ConfidenceInAgeEstimate,
                CanBeUsedForAgingStudies,Tissue,Female,SpeciesCommonName,ExternalSampleID,Folder,
                Experiment,RegenExperimentGroup,AnimalID,AnimalName)
rm(infoFrogs.csv)
if ("Female" %in% colnames(infoFrogs)) {
  infoFrogs$Female[which(is.na(infoFrogs$Female))] <- "NA"
  infoFrogs$Female <- factor(infoFrogs$Female, levels=c(0,1,"NA"))
  levels(infoFrogs$Female) <- c("Male","Female","NA")
}

#merge and save
info_data_depo <- dplyr::left_join(info_data_depo, rbind(infoAllSamp, infoHuman, infoFrogs), by="Basename") %>%
  dplyr::mutate(SpeciesLatinName = factor(SpeciesLatinName, levels=c("Ambystoma mexicanum","Homo sapiens",
                                                                     "Xenopus laevis","Xenopus tropicalis"))) %>%
  dplyr::arrange(SpeciesLatinName, OriginalOrderInBatch) %>%
  dplyr::filter(SpeciesLatinName %in% c("Ambystoma mexicanum"))
out.xlsx="SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/info_DataDeposition_SampleSheet.xlsx"
write.xlsx(info_data_depo,out.xlsx,zoom=160,colWidths=10,
           borders="all",headerStyle=createStyle(border="TopBottomLeftRight",textDecoration="bold")) #?buildWorkbook

#load methylation data and save
datAllSamp_tp.rdata=c('N131.ET0087.SalamanderMaxYun/NormalizedData/all_probes_sesame_normalized.Rdata')
dat_data_depo <- transpose_dat(loadRData(datAllSamp_tp.rdata) %>% as.data.frame(), "Basename") %>%
  dplyr::filter(Basename %in% info_data_depo$Basename)
if (length(unique(dat_data_depo$Basename)) != length(unique(info_data_depo$Basename))) stop("Info / Data Mismatch")
out.csv="SpeciesSubsetAnalyses/Subset_AxolotlN131_FinalFigures_Final/dat_DataDeposition_NormalizedData.csv"
write.table(dat_data_depo,out.csv,sep=',',row.names=F,quote=F)

###############################################################################



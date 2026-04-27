cd ~/Dropbox/MyResearchFiles/Horvath_mammalian_meth/SpeciesSubsetAnalyses/
clear > Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_prefiltering_subCPGaxolotl.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_prefiltering_subCPGcombinationmiddlefilter.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_MeanMethylation_AxolotlN131ExpmtDemethylation.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_MeanMethylation_AxolotlN131SkinAndTestSet.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_Clocks_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksLimbTail_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksLimb_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksTail_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksLaterLife_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksLaterLifeLimbTail_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksLaterLifeLimb_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksLaterLifeTail_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLifeLimbTail_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLifeLimb_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLifeTail_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_toAxolotlN131ExpmtLimbPersistRegen.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_toAxolotlN131ExpmtLimbRepRegen.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_toAxolotlN131ExpmtLimbRegen.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_toAxolotlN131ExpmtTailRegen.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_toAxolotlN131ExpmtDemethylation.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_toAxolotlN131ExpmtTERTTailRegen.R Subset_Axolotl.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_toAxolotlN131SkinAndTestSet.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_toAxolotlN131Branch2.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksEarlyLife_Final_toAxolotlN146.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeEWAS_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeEWASLaterLife_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeEWASEarlyLife_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeEWASExpmtTERTTailRegen_Final_Analysis.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOO_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOLimbTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOLimb_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOLaterLife_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOLaterLifeLimbTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOLaterLifeLimb_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOLaterLifeTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOEarlyLife_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOEarlyLifeLimbTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOEarlyLifeLimb_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOEarlyLifeTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOO_Final_Analysis_2.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOLaterLife_Final_Analysis_2.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOEarlyLife_Final_Analysis_2.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOFOEarlyLife_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOFOEarlyLifeLimbTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOFOEarlyLifeLimb_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOFOEarlyLifeTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOFOEarlyLife_Final_Analysis_2.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksAfterTwoYears_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksAfterTwoYearsLimbTail_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksAfterTwoYearsLimb_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksAfterTwoYearsTail_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksBeforeTwoYears_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksBeforeTwoYearsLimbTail_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksBeforeTwoYearsLimb_Final_basedOnAll.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_ClocksBeforeTwoYearsTail_Final_basedOnAll.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOAfterTwoYears_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOAfterTwoYearsLimbTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOAfterTwoYearsLimb_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOAfterTwoYearsTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOBeforeTwoYears_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOBeforeTwoYearsLimbTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOBeforeTwoYearsLimb_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOBeforeTwoYearsTail_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOAfterTwoYears_Final_Analysis_2.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_AgeLOOBeforeTwoYears_Final_Analysis_2.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet HumanAxolotlN131_AgeLOFOBalanceEarlyLife_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet HumanAxolotlN131_AgeLOFOBalanceEarlyLifeLimbTail+Pan_Final_Analysis.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlClawedFrogN131N140_AgeLOO_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlClawedFrogN131N140_AgeLOOLimbTail+Pan_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlClawedFrogN131N140_AgeLOOLaterLife_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlClawedFrogN131N140_AgeLOOLaterLifeLimbTail+Pan_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlClawedFrogN131N140_AgeLOOEarlyLife_Final_Analysis.R Subset_AxolotlN131.Rout
R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlClawedFrogN131N140_AgeLOOEarlyLifeLimbTail+Pan_Final_Analysis.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_FinalFigures_Final_basedOnAll.R Subset_AxolotlN131.Rout

R CMD BATCH --no-restore --no-save --quiet Subset_AxolotlN131_PredictionResults.R Subset_AxolotlN131.Rout

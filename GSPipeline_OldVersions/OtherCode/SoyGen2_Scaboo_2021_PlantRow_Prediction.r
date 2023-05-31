



##### Get Predictions  

#  Marker Set: 1K Marker Set Agriplex Set
#  Training Set  (NUST lines + SCB lines) as training set  (Complete data after QC) 
#  Target Set (2021_Plant Rows MO) 
#  GP Model - RidgeRegression in 'bWGR' package 

  
  setwd("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data/Miranda data/")
  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
     
  
  infileNUSTSCB_PR1K <- "SCB_PR_Raw_Genotypes_Merged_Raw_1K_KNNI_200_60.vcf"
   
  gt2d_NUST__PR1K <- VCFtoDF(infileNUSTSCB_PR1K)
  
  
  # dim(gt2d_NUST__PR1K)
# [1] 1243 5159
  
  
  
  
  infileVCFPR <- "2021_PlantRow_Mod.vcf"
  gt2d_SCBPR <- VCFtoDF(infileVCFPR)
  testIDs <- colnames(gt2d_SCBPR)[6:ncol(gt2d_SCBPR)]
  
  testIDsMod <- gsub("X","",testIDs)
  testIDsMod2 <- gsub("[-_.]","",testIDsMod)
  TargetIDs <- colnames(gt2d_NUST__PR1K)[which(gsub("[-_.]","",colnames(gt2d_NUST__PR1K)) %in% testIDsMod2)]
  
  IDs <-  gsub("[-_.]","",colnames(gt2d_NUST__PR1K)[6:ncol(gt2d_NUST__PR1K)] )
  TrainIDs <- setdiff(IDs,TargetIDs)
  
 
  
  NUST_BLUEs <- read.csv(infileNUSTBLUEs)
  SCB_BLUEs <- read.csv(infileSCBBLUEs)[,-1]
  colnames(SCB__BLUEs)[1] <- "GermplasmId"
  
  
  BLUEsIndices <-  which(colnames(NUST_BLUEs) %in% colnames(SCB__BLUEs))
  
  SCB_BLUEs <- rbind(NUST_BLUEs[,BLUEsIndices],SCB__BLUEs)
  
#####  
  
  nTraits <- ncol(SCB_BLUEs)-1
  trait <- "YieldBuA"
  SCB_PR_Data_Table_Num_1K <- getMergedData(gt2d_NUST__PR1K,SCB_BLUEs,TargetIDs)
  
  traits <- colnames(SCB_BLUEs)[-1]
  
  outputDF_SCB_PR_1K_List <- list()
  nTr <- 1
  
  for(trait in traits){
  
  
  SCB_PR_Data_Table_Num_1K_Filt <- getProcessedData(SCB_PR_Data_Table_Num_1K,trait)
  
     
  nTraits <- ncol(SCB_BLUEs)-1
  GPModel <- "rrBLUP (bWGR)"
  outputDF_SCB_PR_1K <- getRankedPredictedValues( SCB_PR_Data_Table_Num_1K_Filt,nTraits,trait,GPModel,fixedX=NULL,optTS=NULL)
  outputDF_SCB_PR_1K_List[[nTr]] <- outputDF_SCB_PR_1K
   
   nTr<- nTr+1
  }
  
   nTr <- 1
   outputDF_SCB_PR_1K_Comb <- outputDF_SCB_PR_1K_List[[nTr]]
   
   for(nTr in 2:nTraits){
      outputDF_SCB_PR_1K_Comb <- merge(outputDF_SCB_PR_1K_Comb,outputDF_SCB_PR_1K_List[[nTr]],by="LineID") 
  
   }
   
  Reliability_Indices <- grep("Upper", colnames(outputDF_SCB_PR_1K_Comb))
  outputSCB_PR_1K_Traits1 <- outputDF_SCB_PR_1K_Comb[,-Reliability_Indices[-1]]
  outputSCB_PR_1K_Traits <- cbind(outputSCB_PR_1K_Traits1[,-Reliability_Indices[1]],outputSCB_PR_1K_Traits1[,Reliability_Indices[1]])
 
 write.csv(outputSCB_PR_1K_Traits,"2021_PlantRow_Predictions_All_Traits.csv")
 write.csv(outputDF_SCB_PR_1K,"2021_PlantRow_Predictions.csv")
  
  
  

##########################################################
  
  
  
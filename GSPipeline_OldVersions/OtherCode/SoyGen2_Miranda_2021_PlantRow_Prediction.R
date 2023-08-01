



##### Get Predictions  

#  Marker Set: 1K Marker Set Agriplex Set
#  Training Set  (NUST lines + MRD AYT lines) as training set  (Complete data after QC) 
#  Target Set (2021_Plant Rows NDSU) 
#  GP Model - RidgeRegression in 'bWGR' package 

  
  setwd("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data/Miranda data/")
  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
     
  
  infileNUSTMRD_AYTPR1K <- "NUST_MRD_PR_AYT_Raw_Genotypes_Merged_Raw_1K_KNNI_200_60.vcf"
   
  gt2d_NUST_AYT_PR1K <- VCFtoDF(infileNUSTMRD_AYTPR1K)
  
  
  # dim(gt2d_NUST_AYT_PR1K)
# [1] 1243 5159
  
  
  
  
  infileVCFPR <- "2021_PlantRow_Mod.vcf"
  gt2d_MRDPR <- VCFtoDF(infileVCFPR)
  testIDs <- colnames(gt2d_MRDPR)[6:ncol(gt2d_MRDPR)]
  
  testIDsMod <- gsub("X","",testIDs)
  testIDsMod2 <- gsub("[-_.]","",testIDsMod)
  TargetIDs <- colnames(gt2d_NUST_AYT_PR1K)[which(gsub("[-_.]","",colnames(gt2d_NUST_AYT_PR1K)) %in% testIDsMod2)]
  
  IDs <-  gsub("[-_.]","",colnames(gt2d_NUST_AYT_PR1K)[6:ncol(gt2d_NUST_AYT_PR1K)] )
  TrainIDs <- setdiff(IDs,TargetIDs)
  
 
  infileNUSTBLUEs <- "Pheno.csv"
  infileMRDBLUEs <- "NDSU_AYT_BLUEs_2019_2021.csv"
  
  NUST_BLUEs <- read.csv(infileNUSTBLUEs)
  MRD_AYT_BLUEs <- read.csv(infileMRDBLUEs)[,-1]
  colnames(MRD_AYT_BLUEs)[1] <- "GermplasmId"
  
  
  BLUEsIndices <-  which(colnames(NUST_BLUEs) %in% colnames(MRD_AYT_BLUEs))
  
  NUST_MRD_BLUEs <- rbind(NUST_BLUEs[,BLUEsIndices],MRD_AYT_BLUEs)
  
#####  
  
  nTraits <- ncol(NUST_MRD_BLUEs)-1
  trait <- "YieldBuA"
  NUST_MRD_PR_Data_Table_Num_1K <- getMergedData(gt2d_NUST_AYT_PR1K,NUST_MRD_BLUEs,TargetIDs)
  
  traits <- colnames(NUST_MRD_BLUEs)[-1]
  
  outputDF_NUST_MRD_PR_1K_List <- list()
  nTr <- 1
  
  for(trait in traits){
  
  
  NUST_MRD_PR_Data_Table_Num_1K_Filt <- getProcessedData(NUST_MRD_PR_Data_Table_Num_1K,trait)
  
     
  nTraits <- ncol(NUST_MRD_BLUEs)-1
  GPModel <- "rrBLUP (bWGR)"
  outputDF_NUST_MRD_PR_1K <- getRankedPredictedValues( NUST_MRD_PR_Data_Table_Num_1K_Filt,nTraits,trait,GPModel,fixedX=NULL,optTS=NULL)
  outputDF_NUST_MRD_PR_1K_List[[nTr]] <- outputDF_NUST_MRD_PR_1K
   
   nTr<- nTr+1
  }
  
   nTr <- 1
   outputDF_NUST_MRD_PR_1K_Comb <- outputDF_NUST_MRD_PR_1K_List[[nTr]]
   
   for(nTr in 2:nTraits){
      outputDF_NUST_MRD_PR_1K_Comb <- merge(outputDF_NUST_MRD_PR_1K_Comb,outputDF_NUST_MRD_PR_1K_List[[nTr]],by="LineID") 
  
   }
   
  Reliability_Indices <- grep("Upper", colnames(outputDF_NUST_MRD_PR_1K_Comb))
  outputNUST_MRD_PR_1K_Traits1 <- outputDF_NUST_MRD_PR_1K_Comb[,-Reliability_Indices[-1]]
  outputNUST_MRD_PR_1K_Traits <- cbind(outputNUST_MRD_PR_1K_Traits1[,-Reliability_Indices[1]],outputNUST_MRD_PR_1K_Traits1[,Reliability_Indices[1]])
 
 write.csv(outputNUST_MRD_PR_1K_Traits,"2021_PlantRow_Predictions_All_Traits.csv")
 write.csv(outputDF_NUST_MRD_PR_1K,"2021_PlantRow_Predictions.csv")
  
  
  

##########################################################
  
  
### 6K NUST_MRD_AYT_Raw
  infileNUSTMRD_AYTPR6K <- "NUST_MRD_PR_AYT_Raw_Genotypes_Merged_Raw_6K_KNNI_200_60.vcf"
  gt2d_NUST_AYT_PR6K <- VCFtoDF(infileNUSTMRD_AYTPR6K)
  dim(gt2d_NUST_AYT_PR6K)
  
  BLUEsIndices <-  which(colnames(NUST_BLUEs) %in% colnames(MRD_AYT_BLUEs))
  
  NUST_MRD_BLUEs <- rbind(NUST_BLUEs[,BLUEsIndices],MRD_AYT_BLUEs)
  
  nTraits <- ncol(NUST_MRD_BLUEs)-1
  trait <- "YieldBuA"
  
  NUST_MRD_PR_Data_Table_Num_6K <- getMergedData(gt2d_NUST_AYT_PR6K,NUST_MRD_BLUEs,TargetIDs)
  
  
  
  
  trait <- "YieldBuA" 
  NUST_MRD_PR_Data_Table_Num_6K_Filt <- getProcessedData(NUST_MRD_PR_Data_Table_Num_6K,trait)
  
     
  nTraits <- ncol(NUST_MRD_BLUEs)-1
  GPModel <- "rrBLUP (bWGR)"
  outputDF_NUST_MRD_PR_6K <- getRankedPredictedValues(NUST_MRD_PR_Data_Table_Num_6K_Filt,nTraits,trait,GPModel,fixedX=NULL,optTS=NULL)
   
   
   Data_Table_Num_Filt_List <- NUST_MRD_PR_Data_Table_Num_6K_Filt
  
  write.csv(outputDF_NUST_MRD_PR_6K,"2021_PlantRow_Predictions.csv")

  
  
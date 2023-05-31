
library(vcfR)
library(dplyr)
library(rrBLUP)

### Test 2 for NUST Data 
# ## Set input filenames 
# infileVCF1 <- "2021_03_19_NUST_Geno_Master_Filt_Imp_TASSEL.vcf" 
# infileMetaTable <- "NUST_2004_2021_Line_Program_Metadata_V2.csv"
# infileBLUEs <- "BLUE_GermplasmId_Output_table_wide_duplic_checks__2021_10_21_years28.csv" 
# infileLinkMap <- "GeneticMap_SoySNP_6K_cM.csv"
# infileVCF <- "NUST_geno_filt_imp.vcf"
# infileTestSet <- "Agriplex_2021_091_PRISM.csv"

################### Set working dir  #####

# NUST_Meta_Table <- read.csv("NUST_2004_2021_Line_Program_Metadata_V2.csv")
# NUST_Genotypes_VCF <- read.table("2021_03_19_NUST_Geno_Master_Filt_Imp_TASSEL.vcf")
# NUST_BLUEs <- read.csv("BLUE_GermplasmId_Output_table_wide_duplic_checks__2021_10_21_years28.csv",header=TRUE)
# LinkMap <- read.csv(infileLinkMap)


# WorkDir <- "C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data"
# setwd(WorkDir)

setwd("C:/Users/ivanv/Desktop/UMN_GIT/DataShare_Demo/")
infileBLUEs <- "Pheno.csv" 
infileVCF <- "Geno.vcf"
infileTestSet <- "Target.csv"


NUST_Genotypes_VCF <- read.table(infileVCF)
NUST_BLUEs <- read.csv(infileBLUEs,header=TRUE)
NUST_Test_Data_Table <- read.csv(infileTestSet)


#### Single Trait 
  
 nTraits <- ncol(NUST_BLUEs)-1
 testIDs <- as.character(NUST_Test_Data_Table[,1])
 trait <- "YieldBuA"
  
# trait <- "Oil"
## VCF to DF
 gt2d_NUST <- VCFtoDF(infileVCF)  

 Pheno <- NUST_BLUEs
 NUST_Data_Table_Num <- getMergedData(gt2d_NUST,NUST_BLUEs,testIDs)
 NUST_Data_Table_Num_Filt_1K <- getProcessedData(NUST_Data_Table_Num,trait) 
 
 ####
 
### rTASSEL Test
 
 infileVCF <- "Geno.vcf"
 
    tasGeno <- rTASSEL::readGenotypeTableFromPath(
		path = infileVCFstr
	)
	
	tasSumExp <- rTASSEL::getSumExpFromGenotypeTable(
    tasObj = tasGeno )
	
	
	infileVCF <- "Geno.vcf"
	
	tasGeno <- rTASSEL::readGenotypeTableFromPath(
	  path = infileVCF
	)
	
	
	
	tasGenoDF <- (SummarizedExperiment::assays(tasSumExp)[[1]])
	SummarizedExperiment::colData(tasSumExp)[,"Sample"]
    SummarizedExperiment::rowData(tasSumExp)[,"tasselIndex"] 
	
	
## GP ST Test

  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
  GPModel <- "rrBLUP (bWGR)"

##################  
  
  colSelected <- c("emRR","emBB","emBL")
  emCVR_outTable <- rbind(emCVR_NUST[colSelected],emCVR_SCB[colSelected],emCVR_NUST_SCB[colSelected]) 
 
  
## STPGA 
  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/App/GS_Pipeline_Jan_2022_FnsApp.R")
  noToReduce <-  2237
  nTrainToSelect <- 500 
  testIds <- testIDs 
  optTS <- NULL
  trait <- "YieldBuA" 
  GAParameters <- list(InitPop=NULL,npop=100, nelite=10, mutprob=.5, mutintensity = 1,
         niterations=100,minitbefstop=50, tabu=TRUE,errorstat="PEVMEAN2",tabumemsize = 1,plotiters=FALSE, lambda=1e-6, mc.cores=10)
  
  
  NUST_Data_Table_Num_Filt <- NUST_Data_Table_Num_Filt_1K
  Train_STPGA <- c()
  system.time({
	 if(is.null(optTS)){
	  TS_STPGA <- getOptimalTS(NUST_Data_Table_Num_Filt,trait,nTraits,noToReduce,nTrainToSelect,GAParameters)
	 }
   })
   if(!is.null(optTS)){
    TS_STPGA <- as.character(as.vector(optTS))
   }
   TS_Random <- getRandomTS(NUST_Data_Table_Num_Filt,trait,nTraits,noToReduce,nTrainToSelect) 
   PA_Table <- getTSComparisons(NUST_Data_Table_Num_Filt,TS_STPGA,TS_Random,trait,nTraits,testIDs)
	 
#save(TS_STPGA,TS_Random,PA_Table,file="TS_Demo_PreComputed_Out.RData")  
 
  Train_STPGA <- TS_STPGA
  Train_Random <- TS_Random
  Data_Table_Num_Filt_List <- NUST_Data_Table_Num_Filt
  TS_MT <- getTSComparisonsMT(Data_Table_Num_Filt_List,Train_STPGA,Train_Random,trait,nTraits,testIds,optTS=NULL)


###########################################################################################
## CVR
 
 k <- 5
 nIter <- 5
 emCVR <- getemCVR(NUST_Data_Table_Num_Filt,trait,nTraits,k,nIter)
 
 
  PATable <- rbind(c("RR","BB","BL"),round(emCVR[c("emRR","emBB","emBL")],digits=2)) 
  rownames(PATable)<- c("Prediction Model","Prediction Accuracy")
  colnames(PATable) <- rep("",3)
  PATable
  
 
 
 
    PATableComb <- c()
    for(nSelT in 1:length(trait)){
            i <- nSelT
           PATab <- getemCVR(NUST_Data_Table_Num_Filt,trait[i],nTraits,k,nIter)
           PATable <- round(PATab[c("emRR","emBB","emBL")],digits=2)
           PATable2 <- c(trait[i],PATable)
           PATableComb <- rbind(PATableComb,PATable2)
    }
       PATableComb <- rbind(c("Trait","RR","BB","BL"),PATableComb)
       colnames(PATableComb) <- rep("",ncol(PATableComb))
   }
 
 
 k <- 2
 nIter <- 2
 mtCVR <- getMTCVR(NUST_Data_Table_Num_Filt,trait,nTraits,k,nIter) 
 
### 
################# 
    
  outputDF_NUST_Agp_1K <- getRankedPredictedValues(NUST_Data_Table_Num_Filt_1K,nTraits,trait[1],GPModel,optTS=NULL)
 
 
 
  
####### Multiple Traits
  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
 nTraits <- ncol(NUST_BLUEs)-1
 testIDs <- as.character(NUST_Test_Data_Table[,1])
 trait <- c("YieldBuA","Oil","Protein")
 
 gt2d_NUST <- VCFtoDF(infileVCF)


 Pheno <- NUST_BLUEs
  
 NUST_Data_Table_Num <- getMergedData(gt2d_NUST,NUST_BLUEs,testIDs)
 NUST_Data_Table_Num_Filt <- getProcessedData(NUST_Data_Table_Num,trait)
   
  
### MT Test
  
  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
  
  trait <- c("YieldBuA","Oil")
  Data_Table_Num_Filt_List <- NUST_Data_Table_Num_Filt_1K
  GPModelMT <- "BRR (BGLR)"
  outputDF_BRR <- getRankedPredictedValuesMT(NUST_Data_Table_Num_Filt,nTraits,trait,GPModelMT,optTS=NULL)
 
  GPModelMT <- "SpikeSlab (BGLR)"
  outputDF_SS<- getRankedPredictedValuesMT(NUST_Data_Table_Num_Filt,nTraits,trait,GPModelMT,optTS=NULL)

  GPModelMT <- "RKHS (BGLR)"
  outputDF_RKHS <- getRankedPredictedValuesMT(NUST_Data_Table_Num_Filt,nTraits,trait,GPModelMT,optTS=NULL)
 
    
  GPModelMT <- "RKHS (BGLR)"
  Data_Table_Num_Filt_List <- NUST_Data_Table_Num_Filt

  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
  
  trait <- c("YieldBuA","Oil")
  GPModelMT <- "GBLUP (SOMMER)"
  outputDF_GBLUP <- getRankedPredictedValuesMT(NUST_Data_Table_Num_Filt,nTraits,trait,GPModelMT,optTS=NULL)
 
  

 
####### Imputation Expmnt 1C Compare NUST-NUST AgPlex data 1K vs 6K projected from 1K 
 
 infileBLUEs <- "Pheno.csv" 
 infileVCF <- "NUST_AgPlex_Merged_Genotype_Table_KNNImp_100_60.vcf"
 infileTestSet <- "Target.csv"
 gt2d_NUST_6K <- VCFtoDF(infileVCF)
 NUST_BLUEs <- read.csv(infileBLUEs,header=TRUE)
 
 
#### 
 
 Pheno <- NUST_BLUEs
 testIDs_Full <- as.character(NUST_Test_Data_Table[,1])
 NUST_Data_Table_Num_6K <- getMergedData(gt2d_NUST_6K,NUST_BLUEs,testIDs)
 trait <- "YieldBuA"
 NUST_Data_Table_Num_Filt_6K <- getProcessedData(NUST_Data_Table_Num_6K,trait)
 
 GPModel <- "rrBLUP (bWGR)"
 outputDF_NUST_Agp_1K <- getRankedPredictedValues(NUST_Data_Table_Num_Filt_1K,nTraits,trait[1],GPModel,optTS=NULL)
 outputDF_NUST_Agp_6K <- getRankedPredictedValues(NUST_Data_Table_Num_Filt_6K,nTraits,trait[1],GPModel,optTS=NULL)
 
 outMerged <- merge(outputDF_NUST_Agp_1K,outputDF_NUST_Agp_6K,by="LineID") 
 cor(outMerged[,2],outMerged[,4]) 
 round(cor(outMerged[,2],outMerged[,4]),digits=2)
 #[1] 0.47
 
## BLUEs estimation for Siddhi's evaluations 
## Accuracy between 1K and 6K_Imputed set in predicting target values - 0.47 



### Imputation Expmnt 1A: using only NUST 6K Genotype Table   
### Mask genotypes in contiguous segments using only agriplex markers  in filtered and imputed NUST Raw genotype file 
### in 50% of randomly selected samples 

  WorkDir <- "C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data"
  setwd(WorkDir)
 
  NUST_Geno_Filt_Imp_V2 <- read.table("NUST_geno_filt_imp_v2.vcf",header=TRUE)
  NUST_Geno_Raw_Filt_Imp_V2 <- read.table("NUST_Geno_Raw_Filtered_Imputed_V2.vcf",header=TRUE)
  diff_ssID_Indices <- which(! NUST_Geno_Raw_Filt_Imp_V2[,"ID"] %in% NUST_Geno_Filt_Imp_V2[,"ID"])
  nCol_Geno <- ncol(NUST_Geno_Raw_Filt_Imp_V2)
  nColSel <- sample(c(10:nCol_Geno),0.5*nCol_Geno)
   
#########################################################

  ssID_ColSel <- colnames(NUST_Geno_Raw_Filt_Imp_V2)[nColSel]
  NUST_Geno_Raw_Filt_Imp_V2_Mod <- NUST_Geno_Raw_Filt_Imp_V2
  NUST_Geno_Raw_Filt_Imp_V2_Mod[diff_ssID_Indices,nColSel] <- apply(NUST_Geno_Raw_Filt_Imp_V2[diff_ssID_Indices,nColSel],2,function(x) gsub("[0-1]/[0-1]","./.",x))
  write.table(NUST_Geno_Raw_Filt_Imp_V2_Mod,"NUST_Geno_Raw_Filt_Imputed_V2_Mod.txt",sep="\t",quote=FALSE,row.names=FALSE) 
  NUST_Geno_Raw_Filt_Imp_V3 <- NUST_Geno_Raw_Filt_Imp_V2[,nColSel] 
  
 
#########################################################
### Mask genotypes in contiguous segments using only agriplex markers  in filtered and imputed NUST Raw genotype file 
### in 40% of randomly selected samples 

  nColSel1 <- sample(c(10:nCol_Geno),0.4*nCol_Geno)
  NUST_Geno_Raw_Filt_Imp_V2_Mod <- NUST_Geno_Raw_Filt_Imp_V2
  NUST_Geno_Raw_Filt_Imp_V2_Mod[diff_ssID_Indices,nColSel1] <- apply(NUST_Geno_Raw_Filt_Imp_V2[diff_ssID_Indices,nColSel1],2,function(x) gsub("[0-1]/[0-1]","./.",x))
  write.table(NUST_Geno_Raw_Filt_Imp_V2_Mod,"NUST_Geno_Raw_Filt_Imputed_V2_Mod_0.4.txt",sep="\t",quote=FALSE,row.names=FALSE) 
 
  
##### Estimate % NA imputed and imputation accuracy 
 source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
 
  checkDiff <- function(x,gt2d_NUST_DF,gt2d_NUST_Mod_DF){
  
    NA_Ind <- which(is.na(gt2d_NUST_DF[,x]))
	NA_IndMod <- which(is.na(gt2d_NUST_Mod_DF[,x]))
	NA_Indices <- unique(c(NA_Ind,NA_IndMod))
	length(NA_Indices)
	  
	Diff <- as.numeric(as.character(unlist(gt2d_NUST_DF[-NA_Indices,x]))) - as.numeric(as.character(unlist(gt2d_NUST_Mod_DF[-NA_Indices,x])))
  	Diff_Length <- length(which(Diff !=0))
   	return(list(Diff,Diff_Length))
  }

#################### 

  infileVCF_Raw <- "NUST_Geno_Raw_Filtered_Imputed.vcf"
  gt2d_NUST <- VCFtoDF(infileVCF_Raw)
  
  infileVCF_RawMod <- "NUST_Geno_Raw_Filt_Imputed_V2_Mod.vcf"
  gt2d_NUST_Mod <- VCFtoDF(infileVCF_RawMod)
  
  NA_Raw <- sum(unlist(apply(gt2d_NUST,2,function(x) length(which(is.na(x))))))
  NA_Raw_Mod <- sum(unlist(apply(gt2d_NUST_Mod,2,function(x) length(which(is.na(x))))))

   
#############################################
## Convert VCF to DF 
   
 Accuracy_Imp_List <- list() 
 NA_Imputed_List <- list() 
 
 gt2d_NUST_VCF_List <- list()
 
 L_List <- c(60,70,80,90,100)	
 
 for(nL in 1:length(L_List)){ 
   
  infileVCF_l_List <- paste("NUST_Geno_Raw_Filt_Imputed_V2_Mod_KNNI_",L_List[nL],"_",c(20,30,40,50,60),".vcf",sep="")
  gt2d_NUST_l_List <- list()
 
  # NA_l_List <- rep(0,length(infileVCF_l_List))
  # NA_Imputed <- rep(0,length(infileVCF_l_List))
  # Accuracy_Imp <- rep(0,length(infileVCF_l_List))
  
	  for(nVCF in 1:length(infileVCF_l_List)){
	  
		gt2d_NUST_l_List[[nVCF]] <- VCFtoDF(infileVCF_l_List[nVCF])
		#NA_l_List[nVCF] <- sum(unlist(apply(gt2d_NUST_l_List[[nVCF]],2,function(x) length(which(is.na(x))))))
		
	  } 
	gt2d_NUST_VCF_List[[nL]] <- gt2d_NUST_l_List
 }	  
 
 
 
  save(gt2d_NUST_VCF_List,"gt2d_NUST_VCF_Imp_Params.RData")
##################################### 
##### Imputation Accuracy 
 load("gt2d_NUST_VCF_Imp_Params.RData")
  
 Accuracy_Imp_List <- list() 
 NA_Imputed_List <- list() 
 
 L_List <- c(60,70,80,90,100)	
 
 for(nL in 1:length(L_List)){ 
   
  infileVCF_l_List <- paste("NUST_Geno_Raw_Filt_Imputed_V2_Mod_KNNI_",L_List[nL],"_",c(20,30,40,50,60),".vcf",sep="")
  gt2d_NUST_l_List <- list()
  NA_l_List <- rep(0,length(infileVCF_l_List))
  NA_Imputed <- rep(0,length(infileVCF_l_List))
  Accuracy_Imp <- rep(0,length(infileVCF_l_List))
  
	  for(nVCF in 1:length(infileVCF_l_List)){
	  
		# gt2d_NUST_l_List[[nVCF]] <- VCFtoDF(infileVCF_l_List[nVCF])
		
		gt2d_NUST_l_List[[nVCF]] <- gt2d_NUST_VCF_List[[nL]][[nVCF]]
		NA_l_List[nVCF] <- sum(unlist(apply(gt2d_NUST_l_List[[nVCF]],2,function(x) length(which(is.na(x))))))
		
	## %NA that are imputed 
	 
		NA_Imputed[nVCF] <- 1-(NA_l_List[nVCF]/NA_Raw_Mod)
		nColNUSTGeno <- ncol(gt2d_NUST_l_List[[nVCF]])
		#gt2d_NUST_DF <- as.data.frame(gt2d_NUST)
		gt2d_NUST_l_DF <- as.data.frame(gt2d_NUST_l_List[[nVCF]])
		
		gt2d_NUST_DF <- as.data.frame(gt2d_NUST[diff_ssID_Indices,6:nColNUSTGeno])
	    gt2d_NUST_DF1 <- apply(gt2d_NUST_DF,2,function(x) gsub("AA",1,x))
		gt2d_NUST_DF2 <- apply(gt2d_NUST_DF1,2,function(x) gsub("AB",0,x))
		gt2d_NUST_DF3 <- apply(gt2d_NUST_DF2,2,function(x) gsub("BB",-1,x)) 
  
		gt2d_NUST_Mod_DF <- as.data.frame(gt2d_NUST_l_DF[diff_ssID_Indices,6:nColNUSTGeno])
		gt2d_NUST_Mod_DF1 <- apply(gt2d_NUST_Mod_DF,2,function(x) gsub("AA",1,x))
		gt2d_NUST_Mod_DF2 <- apply(gt2d_NUST_Mod_DF1,2,function(x) gsub("AB",0,x))
		gt2d_NUST_Mod_DF3 <- apply(gt2d_NUST_Mod_DF2,2,function(x) gsub("BB",-1,x)) 
	  
	    nColNUSTGenoDF <- ncol(gt2d_NUST_DF3)
	 
		Diff_Imp <- lapply(c(1:nColNUSTGenoDF),function(x) checkDiff(x,gt2d_NUST_DF3,gt2d_NUST_Mod_DF3))
		
		Diff_Imp_Cnt <- (unlist(lapply(Diff_Imp,function(x)x[[2]])))
		Diff_Imp_Per <- 1- (sum(Diff_Imp_Cnt)/ (nrow(gt2d_NUST_DF3)*ncol(gt2d_NUST_DF3)))
	    Accuracy_Imp[nVCF] <- Diff_Imp_Per
		
	  
	  }
   Accuracy_Imp_List[[nL]] <- Accuracy_Imp
   NA_Imputed_List[[nL]] <- NA_Imputed 
 }

#### 
  

###### Imputation Expmnt 1B: 
#### 1K-6K Imputed Genotype Set, 6K NUST_Genotypes  
 
 WorkDir <- "C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data"
 setwd(WorkDir)
 infileBLUEs <- "Pheno.csv" 
 infileVCF <- "Geno.vcf"
 infileTestSet <- "Target.csv"
 source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
 
 NUST_BLUEs <- read.csv(infileBLUEs,header=TRUE)
 Pheno <- NUST_BLUEs
 
 NUST_Genotypes_VCF <- read.table(infileVCF)
 NUST_BLUEs <- read.csv(infileBLUEs,header=TRUE)
 NUST_Test_Data_Table <- read.csv(infileTestSet)

#### Single Trait 
  
 nTraits <- ncol(NUST_BLUEs)-1
 testIDs <- as.character(NUST_Test_Data_Table[,1])
 trait <- "YieldBuA"
  
## VCF to DF

 gt2d_NUST <- VCFtoDF(infileVCF)  

 Pheno <- NUST_BLUEs
 NUST_Data_Table_Num <- getMergedData(gt2d_NUST,NUST_BLUEs,testIDs)
 NUST_Data_Table_Num_Filt_1K <- getProcessedData(NUST_Data_Table_Num,trait) 
  
 k <- 5
 nIter <- 5
 emCVR_NUST_1K <- getemCVR(NUST_Data_Table_Num_Filt_1K,trait,nTraits,k,nIter)
 
 
### Imputed 6K  #gt2d_NUST_Imputed6K <- gt2d_NUST_VCF_List[[5]][[5]]

 infileVCF_NUST_6K_Imp <- "NUST_Geno_Raw_Filt_Imputed_V2_Mod_KNNI_100_60.vcf"
 gt2d_NUST_Imputed6K <- VCFtoDF(infileVCF_NUST_6K_Imp) 
 trait <- "YieldBuA"
 nTraits <- ncol(NUST_BLUEs)-1
 NUST_Data_Table_Num_Imp_6K <- getMergedData(gt2d_NUST_Imputed6K,NUST_BLUEs,testIDs)
 NUST_Data_Table_Num_Filt_Imp_6K <- getProcessedData(NUST_Data_Table_Num_Imp_6K,trait)
 
 k <- 5
 nIter <- 5
 emCVR_NUST_Imp_6K <- getemCVR(NUST_Data_Table_Num_Filt_Imp_6K,trait,nTraits,k,nIter)
 
### NUST_Raw_6K

 infileVCF_NUST_6K <- "NUST_Geno_Raw_6K_Filtered_Imputed_100_60.vcf"
 gt2d_NUST_6K <- VCFtoDF(infileVCF_NUST_6K) 
 trait <- "YieldBuA"
 nTraits <- ncol(NUST_BLUEs)-1
 testIDs <- NULL
 NUST_Data_Table_Num_6K <- getMergedData(gt2d_NUST_6K,NUST_BLUEs,testIDs)
 NUST_Data_Table_Num_Filt_6K <- getProcessedData(NUST_Data_Table_Num_6K,trait)
 
 k <- 5
 nIter <- 5
 emCVR_NUST_6K <- getemCVR(NUST_Data_Table_Num_Filt_6K,trait,nTraits,k,nIter)
 
#####
# emCVR_NUST_6K[c("emRR","emBB","emBL")]
# emRR   emBB   emBL 
# 0.7588 0.7678 0.7881 
 
### 
# emCVR_NUST_Imp_6K[c("emRR","emBB","emBL")]
# emRR   emBB   emBL 
# 0.7380 0.7407 0.7748 
### 
# emCVR_NUST_1K[c("emRR","emBB","emBL")]
# emRR   emBB   emBL 
# 0.7537 0.7250 0.7532 



 
 
######## Other Params Interval region
      
  # infileVCF_10MB <- "NUST_AgPlex_Comb6K_Imp_10MB.vcf"
  # infileVCF_100MB <- "NUST_AgPlex_Comb6K_Imp_100MB.vcf"
 
  # gt2d_NUST_10MB <- VCFtoDF(infileVCF_10MB)
  # gt2d_NUST_100MB <- VCFtoDF(infileVCF_100MB)

  # gt2d_NUST_10MB_DF <- as.data.frame(gt2d_NUST_10MB)
  # gt2d_NUST_100MB_DF <- as.data.frame(gt2d_NUST_100MB)

  # Diff_NUST_Geno <- lapply(c(1:ncol(gt2d_NUST_10MB_DF)),function(x) setdiff(gt2d_NUST_10MB_DF[,x],gt2d_NUST_100MB_DF[,x]))
  # DiffIndices_10MB_100MB <- which(unlist(lapply(Diff_NUST_Geno,length)) !=0)  

  # NUST_Data_Table_Num_Filt_6K[,nColSel]

#######################################################################################
    
## Test merging unfiltered data 
 # PhenoData has 8352 unique Ids; There are no dups when MG is combined with ID
 # There are 1098 Ids with MG attached to the lineID 
 
  # Pheno_Mod <- Pheno 
  # Pheno_Mod[,1] <- StrainIDMod
  # colnames(Pheno_Mod)[1] <- "StrainID"
 
# When MG is removed from the IDs, there are 435 Ids that are duplicated  
## Among the duplicated Ids, there are 338 unique ids and 97 ids that are duplicated 
  
 # (a) How many lines have data across MGs?
 # (b) How many of those lines in (a) also have genotypic information? 
 # (c) Does the merged table contain all the lines in (b)




####
  # Data_Table1_Num <- merge(GenoTable_Filtered,Pheno_Mod,by="StrainID",all=TRUE)
  # dim(Data_Table1_Num)
  # #[1] 8451 1220
  
  # Data_Table1_Num_X <- merge(GenoTable_Filtered,Pheno_Mod,by="StrainID",all.x=TRUE)
  # dim(Data_Table1_Num_X)
  # # [1] 2475 1220


  # # length(which(duplicated(Data_Table1_Num_X[,1])))
  # #[1] 111 
  
  # length(unique(Data_Table1_Num_X[which(duplicated(Data_Table1_Num_X[,1])),1]))
    
  # Data_Table1_Num_Y <- merge(GenoTable_Filtered,Pheno_Mod,by="StrainID",all.y=TRUE)
  # dim(Data_Table1_Num_Y)
  # # [1] 8352 1220
   
  # dupId_Data <- Data_Table1_Num[which(duplicated(Data_Table1_Num[,"StrainID"])),"StrainID"]
  # dupId_Indices <- lapply(dupId_Data,function(x) which(as.character(Data_Table1_Num[,"StrainID"]) %in% as.character(x)))
  # dupId_Indices_Len <- lapply(dupId_Indices,length)
  
  # table(unlist(dupId_Indices_Len)) 
  
  # which(unlist(dupId_Indices_Len) >6)
  
  # dupStrainID_above6 <- dupId_Data[which(unlist(dupId_Indices_Len) >6)]
  
  # dupId_PhenoData <- Pheno_Mod[which(duplicated(Pheno_Mod[,"StrainID"])),"StrainID"]
  # dupId_Indices_Pheno <- lapply(dupId_PhenoData,function(x) which(Pheno_Mod[,"StrainID"] %in% x))
  # PhenoMod_Dup <- lapply(dupId_Indices_Pheno,function(x)Pheno_Mod[x,"MG"])
  
  # length(which(duplicated(Pheno_Mod[which(duplicated(Pheno_Mod[,1])),1])))
  
####################################################################################################################################
 
 # testIDs
 
 
 # noToReduce <- 500 
 # nTrainToSelect <- 200 
 # testIds <- testIDs
 # testIds <- testIDs 
 # optTS <- NULL 

 # save(NUST_Data_Table_Num_Filt,trait,nTraits,testIDs,noToReduce,nTrainToSelect,file="TS_Demo_PreComputed_In.RData") 
 
 

 # Data_Table_Num_Filt_List <- list(Train_Data_Table_Num_Filt,Test_Genotypes_Table_Mod_Num_Filt) 

 # trait <- "YieldBuA"
 # nTraits <- ncol(NUST_BLUEs)-1
 
 
 system.time(rlang::hash(NUST_Data_Table_Num_Filt))
 

 
 
 
 
 ####
  
  setwd("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data/SNP_Wm82_Data_Versions/")
 
####  
  
  snp6K_a1 <- read.table("SNP6K_Wm82.a1.bed")
  snp6K_a2 <- read.table("SNP6K_Wm82.a2.bed")
  snp6K_a4 <- read.table("SNP6K_Wm82.a4.bed") 
    
  snp50K_a1 <- read.table("SNP50K_Wm82.a1.bed")
  snp50K_a2 <- read.table("SNP50K_Wm82.a2.bed")
  snp50K_a4 <- read.table("SNP50K_Wm82.a4.bed") 
  
#### Merge snp6K
  
  colnames(snp6K_a1) <- snp6K_a4[1,]
  colnames(snp6K_a2) <- snp6K_a4[1,]
  colnames(snp6K_a4) <- snp6K_a4[1,]
  snp6K_Comb1 <- merge(snp6K_a1,snp6K_a2,by="SSID")
  snp6K_Versions <- snp6K_Comb1
  colnames(snp6K_Versions) <- gsub(".x",".V1",colnames(snp6K_Versions))
  colnames(snp6K_Versions) <- gsub(".y",".V2",colnames(snp6K_Versions)) 
  dim(snp6K_Versions)
  #[1] 5867    7
  
 ### Merge snp50K
  colnames(snp50K_a1) <- snp6K_a4[1,]
  colnames(snp50K_a2) <- snp6K_a4[1,]
  colnames(snp50K_a4) <- snp6K_a4[1,]
  snp50K_Comb1 <- merge(snp50K_a1,snp50K_a2,by="SSID")
  snp50K_Versions <- snp50K_Comb1
  colnames(snp50K_Versions) <- gsub(".x",".V1",colnames(snp50K_Versions))
  colnames(snp50K_Versions) <- gsub(".y",".V2",colnames(snp50K_Versions)) 
  dim(snp50K_Versions)  

 # [1] 54326     7 
 
 
######  Marker Position in V1
 
  infileVCF <- "Scaboo_DNA6K_16_19_20_21_QC_PO_Genotypes.vcf"
  gt2d_SCB <- VCFtoDF(infileVCF)
  colnames(snp6K_Versions)[1] <- "SNPID"
  SCB_SNP_Version <- merge(snp6K_Versions,gt2d_SCB,by="SNPID")
  Position_Diff_SCB <- as.numeric(as.character(SCB_SNP_Version[,"START.V1"]))- as.numeric(as.character(SCB_SNP_Version[,"Position"]))
  table(Position_Diff_SCB)
  Position_Diff_SCB
  # 0 
  # 5867
 
### NUST data with marker position in V1
  NUST_SNP_Version <- merge(snp6K_Versions,gt2d_NUST,by="SNPID")
  Position_Diff_NUST <- as.numeric(as.character(NUST_SNP_Version[,"START.V1"]))- as.numeric(as.character(NUST_SNP_Version[,"Position"]))
  table(Position_Diff_NUST)
  # 0 
  # 1178 
  
##  MRD data - marker position in V1 
  
  gt2d_MRD1 <- VCFtoDF(infileVCF1)
  gt2d_MRD2 <- VCFtoDF(infileVCF2) 
  
  MRD1_SNP_Version <- merge(snp6K_Versions,gt2d_MRD1,by="SNPID")
  Position_Diff_MRD1 <- as.numeric(as.character(MRD1_SNP_Version[,"START.V1"]))- as.numeric(as.character(MRD1_SNP_Version[,"Position"]))
  table(Position_Diff_MRD1)
  # 0 
  # 1215
 
  MRD2_SNP_Version <- merge(snp6K_Versions,gt2d_MRD2,by="SNPID")
  Position_Diff_MRD2 <- as.numeric(as.character(MRD2_SNP_Version[,"START.V1"]))- as.numeric(as.character(MRD2_SNP_Version[,"Position"]))
  table(Position_Diff_MRD2)
  # Position_Diff_MRD2
  # 0 
  # 1215 
  
  
 ## BD data Version 2 
 
  BD1_SNP_Version <- merge(snp6K_Versions,gt2d_BD1,by="SNPID")
  Position_Diff_BD1_V1 <- as.numeric(as.character(BD1_SNP_Version[,"START.V1"]))- as.numeric(as.character(BD1_SNP_Version[,"Position"]))

  Position_Diff_BD1_V2 <- as.numeric(as.character(BD1_SNP_Version[,"START.V2"]))- as.numeric(as.character(BD1_SNP_Version[,"Position"]))
  length(which(Position_Diff_BD1_V2==0))
  # [1] 909
  length(Position_Diff_BD1_V2)
  # [1] 958
 
  BD2_SNP_Version <- merge(snp6K_Versions,gt2d_BD2,by="SNPID")
  Position_Diff_BD2_V1 <- as.numeric(as.character(BD2_SNP_Version[,"START.V1"]))- as.numeric(as.character(BD2_SNP_Version[,"Position"]))
  Position_Diff_BD2_V2 <- as.numeric(as.character(BD2_SNP_Version[,"START.V2"]))- as.numeric(as.character(BD2_SNP_Version[,"Position"]))

  
  table(Position_Diff_BD2_V1)
  length(Position_Diff_BD2_V2)
  # [1] 958
  length(which(Position_Diff_BD2_V2==0))
  # [1] 909
###### 

### Compare predictions with evaluations 

Pred_Agp <- read.csv("Predictions_Agriplex.csv",header=TRUE)
Pheno_Agp <- read.csv("Phenotypes for the lines genotyped with Agriplex.csv",header=TRUE)

Agp_Comb <- merge(Pred_Agp,Pheno_Agp,by="")
plot(Agp_Comb[,"Yield_bu_acre_13"] ~ as.factor(Agp_Comb[,"Book.Name"]),data=Agp_Comb) 

### 1K 

colnames(outputDF_NUST_Agp_1K)[2] <- "Pred_YieldBuA"
Agp_Comb_1K <- merge(Pred_Agp,outputDF_NUST_Agp_1K,by="LineID",all=TRUE) 

Agp_Comb_Tb_1K <- as_tibble(Agp_Comb_1K)
Agp_Comb_Lines_1K <- Agp_Comb_Tb_1K %>% group_by(LineID)
Agp_Comb_Lines_Mean_1K <- (Agp_Comb_Lines_1K) %>% summarise( 
                       Pred_YieldBuA = mean(Pred_YieldBuA,na.rm=TRUE),
                       YieldBuA =mean(YieldBuA,na.rm=TRUE))

Agp_Comb_Lines_MeanDF <- as.data.frame(Agp_Comb_Lines_Mean_1K) 

cor(Agp_Comb_Lines_MeanDF[,2],Agp_Comb_Lines_MeanDF[,3]) 
# [1] 0.5175288 


### 6K
colnames(outputDF_NUST_Agp_6K)[2] <- "Pred_YieldBuA"

Agp_Comb_6K <- merge(Pred_Agp,outputDF_NUST_Agp_6K,by="LineID",all=TRUE) 

Agp_Comb_Tb_6K <- as_tibble(Agp_Comb_6K)
Agp_Comb_Lines_6K <- Agp_Comb_Tb_6K %>% group_by(LineID)
Agp_Comb_Lines_Mean_6K <- (Agp_Comb_Lines_6K) %>% summarise( 
                       Pred_YieldBuA = mean(Pred_YieldBuA,na.rm=TRUE),
                       YieldBuA =mean(YieldBuA,na.rm=TRUE)) 

Agp_Comb_Lines_MeanDF_6K <- as.data.frame(Agp_Comb_Lines_Mean_6K) 


NA_Indices <- unique(unlist(apply(Agp_Comb_Lines_MeanDF_6K[,2:3],2,function(x) which(is.na(x)))))

cor(Agp_Comb_Lines_MeanDF_6K[-NA_Indices,2],Agp_Comb_Lines_MeanDF_6K[-NA_Indices,3]) 
# [1] 0.2237749

#######################################################################################################################

## Scaboo Data Trial 1 CV & GP 

######## 



 length(Geno_SCB_StrainID)
 #[1] 3424
 trGenoSCB_ID <-  Pheno_SCB[which(Pheno_SCB[,1] %in% Geno_SCB_StrainID),1]
 length(trGenoSCB_ID)
 #[1] 3184
  targetGenoSCB_ID <- setdiff(Geno_SCB_StrainID,trGenoSCB_ID)
  length(targetGenoSCB_ID)
 # [1] 240
 
  
 SCB_Data_Table_Num <- getMergedData(gt2d_SCB,Pheno_SCB,targetGenoSCB_ID)
 SCB_Data_Table_Num_Filt <- getProcessedData(SCB_Data_Table_Num,trait)
 	 
 gt2d <- gt2d_SCB 
 Pheno <- Pheno_SCB 
 testIDs <- targetGenoSCB_ID

 
 k <- 5
 nIter <- 5
 emCVR_SCB <- getemCVR(SCB_Data_Table_Num_Filt,trait,nTraits,k,nIter) 
  

 GPModel <- "rrBLUP (bWGR)"
 outputDF_SCB <- getRankedPredictedValues(SCB_Data_Table_Num_Filt,nTraits,trait,GPModel,optTS=NULL)
  

########  NUST_SCB Combined Data

### Test Set
 SCB_Test_Data_Table <- SCB_Data_Table_Num_Filt[[2]]
 NUST_Test_Data_Table <- NUST_Data_Table_Num_Filt[[2]]
 
 SCB_Test_Data <- SCB_Test_Data_Table[,-1]
 NUST_Test_Data <- NUST_Test_Data_Table[,-1]
 
 SCB_Test_Data_Tr <- cbind(rownames(t(SCB_Test_Data)),t(SCB_Test_Data))
 NUST_Test_Data_Tr <- cbind(rownames(t(NUST_Test_Data)),t(NUST_Test_Data)) 
 
 colnames(NUST_Test_Data_Tr)[1] <- "SSID"
 colnames(SCB_Test_Data_Tr)[1] <- "SSID" 
 
 NUST_SCB_Test_Data_Tr <- merge(SCB_Test_Data_Tr,NUST_Test_Data_Tr,by="SSID")
 
 dim(NUST_SCB_Test_Data_Tr)
 
### Train Set

 SCB_Train_Data_Table <- SCB_Data_Table_Num_Filt[[1]]
 NUST_Train_Data_Table <- NUST_Data_Table_Num_Filt[[1]]
 
 SCB_Train_Data <- SCB_Train_Data_Table[,-c(1,2)]
 NUST_Train_Data <- NUST_Train_Data_Table[,-c(1,2)]
 
 SCB_Train_Data_Tr <- cbind(rownames(t(SCB_Train_Data)),t(SCB_Train_Data))
 NUST_Train_Data_Tr <- cbind(rownames(t(NUST_Train_Data)),t(NUST_Train_Data)) 
 
 colnames(NUST_Train_Data_Tr) <- c("SSID",NUST_Train_Data_Table[,1])
 colnames(SCB_Train_Data_Tr)<- c("SSID",SCB_Train_Data_Table[,1])
 NUST_SCB_Train_Data_Tr <- merge(SCB_Train_Data_Tr,NUST_Train_Data_Tr,by="SSID")
 
 NUST_SCB_Train_Data <- t(NUST_SCB_Train_Data_Tr[c(4:1155,1:3,1156),])
 NUST_SCB_Train_Data_Table <- cbind(rownames(NUST_SCB_Train_Data), NUST_SCB_Train_Data) 
 
 NUST_SCB_Test_Data <- t(NUST_SCB_Test_Data_Tr) 
 NUST_SCB_Test_Data_Table <- cbind(rownames(NUST_SCB_Test_Data), NUST_SCB_Test_Data) 

 colnames(NUST_SCB_Train_Data_Table) <-  NUST_SCB_Train_Data_Table[1,]  
 colnames(NUST_SCB_Test_Data_Table) <-  NUST_SCB_Test_Data_Table[1,]  
 
### NUST_SCB_Data_List 
 
 NUST_SCB_Data_Table_Num_Filt <- list(NUST_SCB_Train_Data_Table[-1,],NUST_SCB_Test_Data_Table[-1,])
 
### NUST_SCB Data CVR and GP ###
  
 k <- 5
 nIter <- 5
 emCVR <- getemCVR(NUST_SCB_Data_Table_Num_Filt,trait,nTraits,k,nIter) 
 
 
 GPModel <- "rrBLUP (bWGR)"
 nTraits <- 4
 trait <- "YieldBuA" 
 outputDF <- getRankedPredictedValues(NUST_SCB_Data_Table_Num_Filt,nTraits,trait,GPModel,optTS=NULL)
  
  
  
  
    if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install(version = "3.14")
 
  options(repos = BiocManager::repositories())
  library(BiocManager)
   
 
  options(rsconnect.http.trace = TRUE)
  
   if(!require("devtools", quietly = TRUE)){
     install.packages("devtools")
   }
   library(devtools)
   
    if(!require("rTASSEL", quietly = TRUE)){
	 devtools::install_bitbucket(
		repo = "bucklerlab/rTASSEL",
		#host = "bitbucket.org",
		ref = "master",
		build_vignettes = FALSE
		#INSTALL_opts = "--no-multiarch"
	 ) 
  } 

#### 
    genoFile <- infileVCF
    ext <- tools::file_ext(genoFile$datapath)
    req(genoFile)
    validate(need(ext == "vcf", "Please upload a vcf file"))
    
	#noTas <-  getTasObj(genoFile$datapath) 
   
    tasGeno <- rTASSEL::readGenotypeTableFromPath(
		path = genoFile$datapath
	)
   
    siteMinCnt <- 
	MAF <- 
	minNotMissing <- 
	
	tasGenoFilt1 <- rTASSEL::filterGenotypeTableSites(
		tasObj = tasGeno,
		siteMinCount = siteMinCnt,
		siteMinAlleleFreq = MAF,
		siteMaxAlleleFreq = 1.0,
		siteRangeFilterType = "none"
    )
  
  
    tasGenoFilt2 <- rTASSEL::filterGenotypeTableTaxa(
	   tasGeno,
	   minNotMissing = MinNotMissing,
	   minHeterozygous = 0,
	   maxHeterozygous = 1,
	   taxa = NULL
    )
   
	GenoFilt1 <- getFilteredSitesGenoData(genoTas,siteMinCnt,MAF)}) 
    GenoFilt2 <- getFilteredTaxaGenoData(genoTas,minNotMissing)
    
	Geno_DF <- reactive(getGenoTas_to_DF(GenoTas(),Geno()))
    GenoFilt1_DF <- reactive(getGenoTas_to_DF(GenoFilt1(),Geno()))
    GenoFilt2_DF <- reactive(getGenoTas_to_DF(GenoFilt2(),Geno()))
    setTasGenoFilt2 <- reactive(input$setGenoFilt2Tas)
	
  
  
  l <- 30
  k <- 10
  impMethod <- "LDKNNI"
   
  if(impMethod=="LDKNNI"){
   tasGenoImp <- imputeLDKNNi(FiltGeno, highLDSSites = l, knnTaxa = k, maxDistance = 1e+07)
  } 
  if(impMethod=="Numeric"){
   tasGenoImp <- imputeNumeric(FiltGeno,byMean = TRUE,nearestNeighbors = 5, distance = c("Euclidean", "Manhattan", "Cosine")[1])
  }
	

   GenoImp <- getImputedData(FiltGeno(),l,k,impMethod) 
   
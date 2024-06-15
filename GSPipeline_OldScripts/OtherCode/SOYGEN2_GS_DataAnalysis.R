
#### Soygen2 Data Prep 
	
########## GATK VCF to  Variants Table in terminal 

# singularity exec ./gatk_latest.sif gatk  VariantsToTable \
		 # -V ./merged_vcf_SNPS.vcf.gz \
		 # -F CHROM -F POS -GF GT -GF GP \
		 # -O mergedVCF_SNPS_V2.table


# singularity exec ./gatk_latest.sif gatk  VariantsToTable \
		 # -V Scaboo_DNA6K_16_19_20_21_QC_PO_Genotypes.vcf  \
		 # -GF GT \
		 # -O Scaboo_DNA6K_16_19_20_21_QC_PO_Genotypes.table 



########## Diers 211004_10-19-21 VCF Modifications  in terminal 


 # awk 'BEGIN{OFS="\t"} !match($3,"^ss[0-9]*"){$3="ss"NR OFS $3}{print $0}' Diers211004_10-19-21_Final.vcf > Diers211004_10-19-21_Final_Mod.vcf
 # head -25 Diers211004_10-19-21_Final.vcf > Diers211004_10-19-21_Final_Header.vcf    
 # tail -1008 Diers211004_10-19-21_Final_Mod.vcf > Diers211004_10-19-21_Final_Mod2.vcf
 # cat Diers211004_10-19-21_Final_Header.vcf Diers211004_10-19-21_Final_Mod2.vcf > Diers211004_10-19-21_Final_Mod3.vcf
 # head -1032 Diers211004_10-19-21_Final_Mod3.vcf > Diers211004_10-19-21_Final_Mod4.vcf
 

 # singularity exec ./gatk_latest.sif gatk IndexFeatureFile -I ./Diers211004_10-19-21_Final_Mod4.vcf 

 
 # singularity exec ./gatk_latest.sif gatk  VariantsToTable \
		 # -V Diers211004_10-19-21_Final_Mod4.vcf    \
		 # -F CHROM -F POS -F ID -F REF -F ALT -GF GT \
		 # -O Diers211004_10-19-21_Final_Mod4.table 


########## Diers 16to19_11-3-21 VCF Modifications in terminal 

 # awk 'BEGIN{OFS="\t"} !match($3,"^ss[0-9]*"){$3="ss"NR OFS $3}{print $0}' Diers16to19_11-3-21_Final.vcf  > Diers16to19_11-3-21_Final_Mod.vcf 
 # head -25 Diers16to19_11-3-21_Final.vcf > Diers16to19_11-3-21_Final_Header.vcf    
 # tail -1007 Diers16to19_11-3-21_Final_Mod.vcf > Diers16to19_11-3-21_Final_Mod2.vcf
 # cat Diers16to19_11-3-21_Final_Header.vcf Diers16to19_11-3-21_Final_Mod2.vcf > Diers16to19_11-3-21_Final_Mod3.vcf
 # head -1031 Diers16to19_11-3-21_Final_Mod3.vcf > Diers16to19_11-3-21_Final_Mod4.vcf
 
 # singularity exec ./gatk_latest.sif gatk  VariantsToTable \
		 # -V Diers16to19_11-3-21_Final_Mod4.vcf    \
		 # -F CHROM -F POS -F ID -F REF -F ALT -GF GT \
		 # -O Diers16to19_11-3-21_Final_Mod4.table 


##### needs testing 
 # singularity exec ./gatk_latest.sif gatk MergeVcfs -I ./Diers211004_10-19-21_Final_Mod4.vcf -I ./Diers16to19_11-3-21_Final_Mod4.vcf -O Diers_2021_Combined.vcf 
 # singularity exec ./gatk_latest.sif gatk CombineVcfs -I ./Diers211004_10-19-21_Final_Mod4.vcf -I ./Diers16to19_11-3-21_Final_Mod4.vcf -O Diers_2021_Combined.vcf 
 
 
 
 ### Test NUST 
   
 WorkDir <- "C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data"
 setwd(WorkDir)
 source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
 
 infileBLUEs <- "Pheno.csv" 
 infileVCF <- "NUST_Geno_Raw_Filtered_Imputed.vcf"
 infileTestSet <- "Target.csv"
 gt2d_NUST_Raw <- VCFtoDF(infileVCF)
 NUST_BLUEs <- read.csv(infileBLUEs,header=TRUE)
 
 
 infileTargetVCF <- "TargetSet_raw.vcf" 
 gt2d_NUST_Trgt <-  VCFtoDF(infileTargetVCF)  
 

 
 Pheno <- NUST_BLUEs
 testIDs <- NULL
 
 
 NUST_Data_Table_Num <- getMergedData(gt2d_NUST_Raw,NUST_BLUEs,testIDs)
 
 trait <- "YieldBuA"
 NUST_Data_Table_Num_Filt <- getProcessedData(NUST_Data_Table_Num,trait)
 
 
 
 k <- 5
 nIter <- 5
 emCVR_NUST <- getemCVR(NUST_Data_Table_Num_Filt,trait,nTraits,k,nIter)
 
 
 GPModel <- "rrBLUP (bWGR)"
 outputDF <- getRankedPredictedValues(NUST_Data_Table_Num_Filt,nTraits,trait[1],GPModel,optTS=NULL)
 
 
 
 
##### Test Scaboo data 

##########

#Steps in  TASSEL: 

# 1) Sort Genotype File : Input:  "Scaboo_DNA6K_16_19_20_21_QC_PO_Genotypes.vcf"
#                         Output: "Scaboo_DNA6K_16_19_20_21_QC_PO_Genotypes_Sorted.vcf"

# ii) Load "Scaboo_DNA6K_16_19_20_21_QC_PO_Genotypes_Sorted.vcf"
# 
# ## Genotypic master file has 6000 sites 3424 taxa
# 
# iii) Apply Minimum allele Freq 0.02  - 4844 sites -  3424-taxa
# 
# iv) Min sites present in individual - 0.9 - 4844 sites - 3268 taxa 
# 
# After applying these two filters, 4844 sites and 3268 taxa - 
# 
# iv) imputation with LDKNNI 
# save as  "Scaboo_DNA6K_16_19_20_21_QC_PO_Genotypes_Sorted_Filter_Imputed.vcf"
   
  
  setwd("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data/Scaboo data/")
  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
    
  infileVCF_SCB <- "Scaboo_DNA6K_16_19_20_21_QC_PO_Genotypes_Sorted_Filter_Imputed.vcf"
  infileBLUEs <- "Scaboo_pheno_BLUEs.csv"
  SCB_Genotypes_VCF <- read.table(infileVCF)
  SCB_BLUEs <- read.csv(infileBLUEs,header=TRUE)

  nTraits <- ncol(SCB_BLUEs)-1
  trait <- "YieldBuA"

# testIDs <- as.character(Test_Data_Table[,1])
 
## VCF to DF
  gt2d_SCB <- VCFtoDF(infileVCF_SCB)
  Pheno_SCB <- SCB_BLUEs 
  
  Geno_SCB <- t(gt2d_SCB)
  dim(Geno_SCB)
  #[1] 3429 6000

  
  Geno_SCB_StrainID <- rownames(Geno_SCB)[6:nrow(Geno_SCB)]
  length(Geno_SCB_StrainID)
  #[1] 3424
   
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
  
  
   
 k <- 5
 nIter <- 5
 nTraits <- 5
 trait <- "YieldBuA"
 emCVR_NUSTRaw_SCB <- getemCVR(NUST_Raw_SCB_Data_Table_Num_Filt,trait,nTraits,k,nIter) 
  
   
 ####  NUST_Raw_SCB Combined

  
 gt2d_NUST_Raw_SCB <- merge(gt2d_NUST_Raw,gt2d_SCB,by="SNPID")
 
 SCB_Traits_Table <- cbind(SCB_BLUEs,rep("MO",nrow(SCB_BLUEs)))
 NUST_Traits_Table <- cbind(NUST_BLUEs,rep("NUST",nrow(NUST_BLUEs)))
  
 colnames(NUST_Traits_Table) <- c(colnames(NUST_BLUEs),"Trial")
 colnames(SCB_Traits_Table) <- c(colnames(SCB_BLUEs),"Trial")
 
 Trait_Indices <-  which(colnames(NUST_Traits_Table) %in% colnames(SCB_Traits_Table))
 
 NUST_SCB_Traits_Table <- rbind(NUST_Traits_Table[,Trait_Indices],SCB_Traits_Table)
  
 #Pheno_NUST_RAW_SCB <- merge(NUST_BLUEs,SCB_BLUEs,by="GermplasmId",all.x=TRUE)
 
 Pheno_NUST_RAW_SCB <-  NUST_SCB_Traits_Table
 testIDs <- NULL
 NUST_Data_Table_Num <- getMergedData(gt2d_NUST_Raw_SCB,Pheno_NUST_RAW_SCB,testIDs)
 
 
 NUST_Raw_SCB_Data_Table_Num <- NUST_Data_Table_Num
 NUST_Raw_SCB_Data_Table_Num_Filt <- getProcessedData(NUST_Raw_SCB_Data_Table_Num,trait)
 
 
 k <- 5
 nIter <- 5
 nTraits <- 5
 trait <- "YieldBuA"
 emCVR_NUSTRaw_SCB <- getemCVR(NUST_Raw_SCB_Data_Table_Num_Filt,trait,nTraits,k,nIter) 
  
 
 
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
 
### NUST_SCB Data CVR and GP
  
 k <- 5
 nIter <- 5
 emCVR <- getemCVR(NUST_SCB_Data_Table_Num_Filt,trait,nTraits,k,nIter) 
 
 
 GPModel <- "rrBLUP (bWGR)"
 nTraits <- 4
 trait <- "YieldBuA" 
 outputDF <- getRankedPredictedValues(NUST_SCB_Data_Table_Num_Filt,nTraits,trait,GPModel,optTS=NULL)
  
    
   
##### Test Miranda data 

  setwd("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data/Miranda data/")

  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
  
  infileVCF1 <- "2021_NDSU_AYT.vcf"
  infileVCF2 <- "2021_PlantRow.vcf"
  infilePheno1<- "NDSU_AYT_means table.csv"
  infilePheno2 <- "Plant Row Weights.csv"
  
  MRD_Genotypes_VCF1 <- read.table(infileVCF1)
  MRD_Genotypes_VCF2 <- read.table(infileVCF2)
  
  gt2d_MRD1 <- VCFtoDF(infileVCF1)
  gt2d_MRD2 <- VCFtoDF(infileVCF2)
  
  MRD_BLUEs1 <- read.csv(infilePheno1,header=TRUE) 
  MRD_BLUEs2 <- read.csv(infilePheno2,header=TRUE) 
  
### Compare REF/ALT in Agriplex in 2021_NDSU_AYT and Scaboo 6K data  
  
  gt_MRD1_SCB <- merge(gt2d_MRD1,gt2d_SCB,by="SNPID")
  gt_MRD1_SCB_Allele_Comp <- (gt_MRD1_SCB[,c("SNPID","Chrom.x","Position.x","REF.x","REF.y","ALT.x","ALT.y")]) 
  REF_Diff_Indices1 <- which(gt_MRD1_SCB_Allele_Comp[,"REF.x"]!= gt_MRD1_SCB_Allele_Comp[,"REF.y"])
  Ref_Alt_Diff_Indices1 <- which(gt_MRD1_SCB_Allele_Comp[REF_Diff_Indices1,"REF.x"]== gt_MRD1_SCB_Allele_Comp[REF_Diff_Indices1,"ALT.y"]) 
     

  
  dim(gt_MRD1_SCB)
  #[1] 1181 3435
  length(REF_Diff_Indices1)
  # [1] 468
  length(Ref_Alt_Diff_Indices1)
  # [1] 468
  
 ### Compare REF/ALT in Agriplex in 2021_PlantRow and Scaboo 6K data 
   
  gt_MRD2_SCB <- merge(gt2d_MRD2,gt2d_SCB,by="SNPID")
  gt_MRD2_SCB_Allele_Comp <- (gt_MRD2_SCB[,c("SNPID","Chrom.x","Position.x","REF.x","REF.y","ALT.x","ALT.y")]) 
  REF_Diff_Indices2 <- which(gt_MRD2_SCB_Allele_Comp[,"REF.x"]!= gt_MRD2_SCB_Allele_Comp[,"REF.y"])
  Ref_Alt_Diff_Indices2 <- which(gt_MRD2_SCB_Allele_Comp[REF_Diff_Indices2,"REF.x"]== gt_MRD2_SCB_Allele_Comp[REF_Diff_Indices2,"ALT.y"]) 

  dim(gt_MRD2_SCB)
  # [1] 1181 5734
  length(REF_Diff_Indices2)
  #[1] 451
  length(Ref_Alt_Diff_Indices2)
  #[1] 451
  
  
### Agriplex Marker ID 
  
  # AgPlex_Markers_MRD <- read.table("Agriplex markers extracted from 50k results.txt")
  
  AgPlex_Markers_MRD <- read.csv("Agriplex markers extracted from 50k results.csv",header=TRUE)
  length(grep("ss",AgPlex_Markers_MRD[,1]))
  length(AgPlex_Markers_MRD[,"alleles.50K"])
  length(grep("ss",as.character(AgPlex_Markers_MRD[,"X.ssID"])))

  alleles.AgP <-strsplit(AgPlex_Markers_MRD[,"alleles.AP"],"/")
  alleles.50K <- strsplit(AgPlex_Markers_MRD[,"alleles.50K"],"/")
    
### 
  setdiff(as.character(AgPlex_Markers_MRD[1:1243,"ID"]),as.character(AgPlex_Markers_MRD[1:1243,"X.ssID"]))
  # character(0)
  
 # length(which(unlist(lapply(alleles.AgP,length)) ==2))
 #[1] 1243
 # length(which(unlist(lapply(alleles.50K,length)) ==2))
 # [1] 39998
  
 
  
 alleles.AgP.Indices <-  which(unlist(lapply(alleles.AgP,length)) ==2)
 alleles.50K.Indices <- which(unlist(lapply(alleles.50K,length)) ==2)
  
  checkAllele <- function(x) {
  
    if(alleles.AgP[[x]][1]==alleles.50K[[x]][1]){s <- 1} 
    if(alleles.AgP[[x]][1]==alleles.50K[[x]][2]){s <- 2}
	return(s)
  }
  
  
  allele.AgP.Checks <- lapply(alleles.AgP.Indices,function(x) checkAllele(x) )
  
  length(which(unlist(allele.AgP.Checks)==1))
  #[1] 716
  
  length(which(unlist(allele.AgP.Checks)==2))
  #[1] 527
  
##### 

## Step 1 - Get Processed Data Table from input files
# gt2d <- VCFtoDF(infileVCF1)
# NUST_Data_Table_Num_Filt <- getProcessedData(gt2d,NUST_BLUEs,NUST_Meta_Table)
 
#gt2d <- VCFtoDF(infileVCF)

# NUST_Data_Table_Num_Filt <- getProcessedData(gt2d,NUST_BLUEs,testIDs,trait)
## Opt TS
#  Data_Table_Num_Filt_List <- NUST_Data_Table_Num_Filt
 
  
  colSelected <- c("emRR","emBB","emBL")
  emCVR_outTable <- rbind(emCVR_NUST[colSelected],emCVR_SCB[colSelected],emCVR_NUST_SCB[colSelected]) 
 
  setwd("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data/Miranda data/SNP_Wm82_Data_Versions/")
 
 
  
 
 ## Diers Data 
 
 
  setwd("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/Data/Diers Final Data/")

  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
  
  infileVCF1 <- "Diers16to19_11-3-21_Final_Mod4.vcf"
  infileVCF2 <- "Diers211004_10-19-21_Final_Mod4.vcf"
  infilePheno1<- "Diers_pheno_BLUEs.csv"

  
  BD_Genotypes_VCF1 <- read.table(infileVCF1)
  BD_Genotypes_VCF2 <- read.table(infileVCF2)
  
  gt2d_BD1 <- VCFtoDF(infileVCF1)
  gt2d_BD2 <- VCFtoDF(infileVCF2)
  
  BD_BLUEs1 <- read.csv(infilePheno1,header=TRUE) 
  
### Compare REF/ALT in Agriplex in 2021_NDSU_AYT and Scaboo 6K data  
  
  gt_BD1_SCB <- merge(gt2d_BD1,gt2d_SCB,by="SNPID")
  gt_BD1_SCB_Allele_Comp <- (gt_BD1_SCB[,c("SNPID","Chrom.x","Position.x","REF.x","REF.y","ALT.x","ALT.y")]) 
  REF_Diff_Indices1 <- which(gt_BD1_SCB_Allele_Comp[,"REF.x"]!= gt_BD1_SCB_Allele_Comp[,"REF.y"])
  Ref_Alt_Diff_Indices1 <- which(gt_BD1_SCB_Allele_Comp[REF_Diff_Indices1,"REF.x"]== gt_BD1_SCB_Allele_Comp[REF_Diff_Indices1,"ALT.y"]) 
    
  dim(gt_BD1_SCB)
#[1]  968 3761

  length(REF_Diff_Indices1)
  #  [1] 292

  length(Ref_Alt_Diff_Indices1)
  # [1] 243 
 
### BCD2 _SCB 
 
  
  gt_BD2_SCB <- merge(gt2d_BD2,gt2d_SCB,by="SNPID")
  gt_BD2_SCB_Allele_Comp <- (gt_BD2_SCB[,c("SNPID","Chrom.x","Position.x","REF.x","REF.y","ALT.x","ALT.y")]) 
  REF_Diff_Indices1 <- which(gt_BD2_SCB_Allele_Comp[,"REF.x"]!= gt_BD2_SCB_Allele_Comp[,"REF.y"])
  Ref_Alt_Diff_Indices1 <- which(gt_BD2_SCB_Allele_Comp[REF_Diff_Indices1,"REF.x"]== gt_BD2_SCB_Allele_Comp[REF_Diff_Indices1,"ALT.y"]) 
    
  dim(gt_BD2_SCB)
#[1]   968 4843
  length(REF_Diff_Indices1)
  #[1] 293
  length(Ref_Alt_Diff_Indices1) 
  # [1] 244
  Ref_Alt_Diff_Var <- setdiff(REF_Diff_Indices1,REF_Diff_Indices1[Ref_Alt_Diff_Indices1])
  length(Ref_Alt_Diff_Var)
  # [1] 49 
  
#######################################################################################################
 
 
 
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
  
############ ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

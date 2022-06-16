
#### Set dir where you have the data files as working dir.. 
WorkDir <- "C:/Users/ivanv/Desktop/UMN_GIT/DataShare_Demo"
setwd(WorkDir)
infileBLUEs <- "Pheno.csv" 
infileVCF <- "Geno.vcf"
infileTestSet <- "Target.csv"

### Set dir to the one with the app functions file and source App Functions file 
source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")

NUST_Genotypes_VCF <- read.table(infileVCF)
NUST_BLUEs <- read.csv(infileBLUEs,header=TRUE)
NUST_Test_Data_Table <- read.csv(infileTestSet)

#### Single Trait 
  
 nTraits <- ncol(NUST_BLUEs)-1
 testIDs <- as.character(NUST_Test_Data_Table[,1])
  trait <- c("YieldBuA","Oil")
 Pheno <- NUST_BLUEs
   
## VCF to DF
 gt2d_NUST <- VCFtoDF(infileVCF)  

### Process Data

 NUST_Data_Table_Num <- getMergedData(gt2d_NUST,NUST_BLUEs,testIDs)
 NUST_Data_Table_Num_Filt <- getProcessedData(NUST_Data_Table_Num,trait) 

 getPredictionData(NUST_Data_Table_Num_Filt,)
	
## GP ST Test

  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
  GPModel <- "rrBLUP (bWGR)"
 
  nTraits <- ncol(NUST_BLUEs)-1
  outputDF_NUST_Agp <- getRankedPredictedValues(NUST_Data_Table_Num_Filt,nTraits,trait[1],GPModel,optTS=NULL)
 
  
 
### GP MT Test
  
  source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
  
  trait <- c("YieldBuA","Oil")
  
  GPModelMT <- "BRR (BGLR)"
  outputDF_BRR <- getRankedPredictedValuesMT(NUST_Data_Table_Num_Filt,nTraits,trait,GPModelMT,optTS=NULL)
 
  GPModelMT <- "SpikeSlab (BGLR)"
  outputDF_SS<- getRankedPredictedValuesMT(NUST_Data_Table_Num_Filt,nTraits,trait,GPModelMT,optTS=NULL)

  GPModelMT <- "RKHS (BGLR)"
  outputDF_RKHS <- getRankedPredictedValuesMT(NUST_Data_Table_Num_Filt,nTraits,trait,GPModelMT,optTS=NULL)
 
 
#########################################


### SoyNAM Data download from GIGWA in vcf format: 

mkdir ZipTemp
for i in *.zip; do unzip -o $i -d ./ZipTemp; echo $i | awk '{out="";split($0,a,"__");out=a[2]".vcf";print out}' | xargs -I '{}'  mv ./ZipTemp/*.vcf ./'{}'; rm ./ZipTemp/*.vcf ; done


## Combine all VCF files
 
 i <- 1
 infileVCF <- paste("project",i,".vcf",sep="")
 gt2d_NAM_Fam<- VCFtoDF(infileVCF)  

 
  gt2d_NAM   <- t(gt2d_NAM_Fam) 
   
  for(i in 2:39){ 
   
    infileVCF <- paste("project",i,".vcf",sep="")
    gt2d_NAM_Fam<- VCFtoDF(infileVCF)  

    gt2d_NAM <- rbind(gt2d_NAM,t(gt2d_NAM_Fam)[-c(1:5),])
  }
  
  tb_NAM <- as_tibble(t(gt2d_NAM)) 
  
####  
  
  outFile <- 'SoyNAM_Geno.vcf'
  ingeno <- tb_NAM
  outfile <- outFile
  DFToVCF(tb_NAM,outFile) 
  
 ####  
  
  gt2d_NAM_Filt <- (gt2d_NAM[-(grep("IA3023",rownames(gt2d_NAM))[-1]),])
  
  tb_NAM <- as_tibble(t(gt2d_NAM_Filt)) 
  
  outFile <- 'SoyNAM_Geno.vcf'
  DFToVCF_NAM(tb_NAM,outFile) 
  ###
  
  
  
  
  
  vcf <- read.vcfR(infileVCF, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = F,return.alleles=T)
  fix_T <- as_tibble(getFIX(vcf))
  
  RefAlleles <- gt[,grep("IA3023", colnames(gt))]]
  
 ##### SoyNAM Phenotypic data 
  
  install.packages("SoyNAM")
  library(SoyNAM)
  data(soynam)
  dim(data.line)
  SoyNAM_Pheno <- data.line
  head(SoyNAM_Pheno)
  
######  
## ImputeAccuracy_Test 
   
  infileVCF_100_60 <- "NUST_Geno_Raw_Filt_Imputed_V2_Mod_KNNI_100_60.vcf"
  infileVCF_70_30 <- "NUST_Geno_Raw_Filt_Imputed_V2_Mod_KNNI_70_30.vcf"
  
  gt2d_NUST_60_20 <- VCFtoDF(infileVCF_60_20)
  gt2d_NUST_70_30 <- VCFtoDF(infileVCF_70_30)
  
  infileVCF_80_40 <- "NUST_Geno_Raw_Filt_Imputed_V2_Mod_KNNI_80_40.vcf"
  gt2d_NUST_80_40 <- VCFtoDF(infileVCF_80_40)
  
  infileVCF_90_50 <- "NUST_Geno_Raw_Filt_Imputed_V2_Mod_KNNI_90_50.vcf"
  gt2d_NUST_90_50 <- VCFtoDF(infileVCF_90_50)
  
  infileVCF_100_60 <- "NUST_Geno_Raw_Filt_Imputed_V2_Mod_KNNI_100_60.vcf"
  gt2d_NUST_100_60 <- VCFtoDF(infileVCF_100_60)
  
  
  NA_Raw <- sum(unlist(apply(gt2d_NUST,2,function(x) length(which(is.na(x))))))
  NA_Raw_Mod <- sum(unlist(apply(gt2d_NUST_Mod,2,function(x) length(which(is.na(x))))))
 
  NA_60_20 <- sum(unlist(apply(gt2d_NUST_60_20,2,function(x) length(which(is.na(x))))))
  NA_70_30 <- sum(unlist(apply(gt2d_NUST_70_30,2,function(x) length(which(is.na(x))))))
  NA_80_40 <- sum(unlist(apply(gt2d_NUST_80_40,2,function(x) length(which(is.na(x))))))
  NA_90_50 <- sum(unlist(apply(gt2d_NUST_90_50,2,function(x) length(which(is.na(x))))))
  NA_100_60 <- sum(unlist(apply(gt2d_NUST_100_60,2,function(x) length(which(is.na(x)))))) 
  
## %NA that are imputed  
   1-(NA_60_20/NA_Raw_Mod)
  #[1] 0.7602596 
   
  1-(NA_70_30/NA_Raw_Mod)
  #[1] 0.8894492
    
  1-(NA_80_40/NA_Raw_Mod)
  # [1] 0.9607631
  
  1-(NA_90_50/NA_Raw_Mod)
  #[1] 0.986077

  1-(NA_100_60/NA_Raw_Mod)
  #[1] 0.9959506
  
 
  gt2d_NUST_DF <- as.data.frame(gt2d_NUST)
  gt2d_NUST_80_40_DF <- as.data.frame(gt2d_NUST_80_40)
 
  Diff_Imp <- lapply(c(6:nColNUSTGeno),function(x) checkDiff(x,gt2d_NUST_DF,gt2d_NUST_80_40_DF))
  length(which(unlist(lapply(Diff_Imp,function(x)x[[2]])) !=0))
  
  gt2d_NUST_70_30_DF <- as.data.frame(gt2d_NUST_70_30)
  Diff_Imp_70_30 <- lapply(c(6:nColNUSTGeno),function(x) checkDiff(x,gt2d_NUST_DF,gt2d_NUST_70_30_DF))
  length(which(unlist(lapply(Diff_Imp_70_30,function(x)x[[2]])) !=0))
  #gt2d_NUST_Mod <- gt2d_NUST_70_30_DF 
  
### Single Column Test   
  x <- 10
  gt2d_NUST_Mod <- gt2d_NUST_l_List[[nVCF]]
  
  NA_Ind <- which(is.na(gt2d_NUST[,x]))
  NA_IndMod <- which(is.na(gt2d_NUST_Mod[,x]))
  NA_Indices <- unique(c(NA_Ind,NA_IndMod))
  
  gt2d_NUST_DF <- as.data.frame(gt2d_NUST[,6:nColNUSTGeno])
  gt2d_NUST_DF1 <- apply(gt2d_NUST_DF,2,function(x) gsub("AA",1,x))
  gt2d_NUST_DF2 <- apply(gt2d_NUST_DF1,2,function(x) gsub("AB",0,x))
  gt2d_NUST_DF3 <- apply(gt2d_NUST_DF2,2,function(x) gsub("BB",-1,x)) 
  
  gt2d_NUST_Mod_DF <- as.data.frame(gt2d_NUST_Mod[,6:nColNUSTGeno])
  gt2d_NUST_Mod_DF1 <- apply(gt2d_NUST_Mod_DF,2,function(x) gsub("AA",1,x))
  gt2d_NUST_Mod_DF2 <- apply(gt2d_NUST_Mod_DF1,2,function(x) gsub("AB",0,x))
  gt2d_NUST_Mod_DF3 <- apply(gt2d_NUST_Mod_DF2,2,function(x) gsub("BB",-1,x)) 
  
  Diff <- as.numeric(as.character(unlist(gt2d_NUST_DF3[-NA_Indices,x])))- as.numeric(as.character(unlist(gt2d_NUST_Mod_DF3[-NA_Indices,x])))
  
  Diff_Length <- length(which(Diff !=0))
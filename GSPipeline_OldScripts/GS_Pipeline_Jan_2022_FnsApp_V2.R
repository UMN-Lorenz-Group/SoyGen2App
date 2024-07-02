##### SoygenGS App Functions

  library(rrBLUP)
  library(vcfR)
  library(dplyr)  
  library(NAM)
  library(bWGR)
  library(STPGA)
  library(BGLR)
  library(sommer)
  
VCFtoDF <- function(infile){
  
  vcf <- read.vcfR(infile, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = F)
  fix_T <- as_tibble(getFIX(vcf))
  
  gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
  colnames(gt2) <- colnames(gt)
  
  gt2a <- apply(gt,2, function(x) gsub("1/1","BB",x))
  gt2b <- gsub("0[/|]0","AA",gt2a)
  gt2c <- gsub("[10][/|][10]","AB",gt2b)
  gt2d <- gsub("\\.[/|]\\.","NA",gt2c)
  
  # gt2[(gt == "1/1")|(gt == "1|1")] <- 'BB'
  # gt2[(gt == "0/0")|(gt == "0|0")] <- 'AA'
  # gt2[(gt == "0/1")|(gt == "1/0")|(gt == "0|1")|(gt == "1|0")] <- 'AB'
  # gt2[(gt == "\\./\\.")|(gt == "\\.|\\.")] <- NA
  
  gt2 <- as_tibble(gt2d) %>%
    mutate(SNPID = fix_T$ID)
  
  gt.simple <- fix_T %>% dplyr::select("ID", "CHROM", POS, REF, ALT) %>%
    dplyr::rename(SNPID=ID, Chrom=CHROM, Position=POS) %>%
    dplyr::left_join(gt2, by = 'SNPID')
  
  return(gt.simple)
}
 



VCFtoDF_V2<- function(infile){
  
 
  vcf <- read.vcfR(infile, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = F,IDtoRowNames = FALSE)
  fix_T <- as_tibble(getFIX(vcf))
  
  gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
  colnames(gt2) <- colnames(gt)
  
  gt2a <- apply(gt,2, function(x) gsub("1/1","BB",x))
  gt2b <- gsub("0[/|]0","AA",gt2a)
  gt2c <- gsub("[10][/|][10]","AB",gt2b)
  gt2d <- gsub("\\.[/|]\\.","NA",gt2c)
  
  # gt2[(gt == "1/1")|(gt == "1|1")] <- 'BB'
  # gt2[(gt == "0/0")|(gt == "0|0")] <- 'AA'
  # gt2[(gt == "0/1")|(gt == "1/0")|(gt == "0|1")|(gt == "1|0")] <- 'AB'
  # gt2[(gt == "\\./\\.")|(gt == "\\.|\\.")] <- NA
  
  gt2 <- as_tibble(gt2d) %>% mutate(SNPID = fix_T$ID)
   
  # gt.simple <- fix_T %>% dplyr::select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT) %>% 
  # dplyr::rename(SNPID=ID, Chrom=CHROM, Position=POS) %>% dplyr::left_join(gt2, by = 'SNPID')
  gt.simple <- fix_T %>% dplyr::select("ID", "CHROM", POS, REF, ALT) %>%
    dplyr::rename(SNPID=ID, Chrom=CHROM, Position=POS) %>%
    dplyr::left_join(gt2, by = 'SNPID')
  
  return(gt.simple)
}

######


VCFtoDF_V3<- function(infile){
  
  
  vcf <- read.vcfR(infile, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = F,IDtoRowNames = FALSE)
  fix_T <- as_tibble(getFIX(vcf,getINFO=TRUE))
  
  gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
  colnames(gt2) <- colnames(gt)
  
  gt2a <- apply(gt,2, function(x) gsub("1/1","BB",x))
  gt2b <- gsub("0[/|]0","AA",gt2a)
  gt2c <- gsub("[10][/|][10]","AB",gt2b)
  gt2d <- gsub("\\.[/|]\\.","NA",gt2c)
  
  # gt2[(gt == "1/1")|(gt == "1|1")] <- 'BB'
  # gt2[(gt == "0/0")|(gt == "0|0")] <- 'AA'
  # gt2[(gt == "0/1")|(gt == "1/0")|(gt == "0|1")|(gt == "1|0")] <- 'AB'
  # gt2[(gt == "\\./\\.")|(gt == "\\.|\\.")] <- NA
  
  gt2 <- as_tibble(gt2d) %>% mutate(ID = fix_T$ID)
  
  gt.simple <- fix_T %>% dplyr::select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO) %>% 
  dplyr::left_join(gt2, by = 'ID')
  #
  
  return(gt.simple)
}




###                         SimpleToVCF                            ###
###   a function to generate vcf file from simple format           ###
######################################################################
# input:
#    infile: the path and file name of the simple format file (.csv)
#    outfile: the path and file name for the vcf file
# output:
# vcf file
# vcf dataframe
SimpleToVCF <- function(infile, outfile){
  input.geno <- as.data.frame(read.csv(infile))
  # order input.geno
  input.geno <- input.geno[order(input.geno$Chrom, input.geno$Position), ]
  # function to change AA, AB, BB to 1/1, 1/0, 0/0, change NA to ./. for each row
  recode.marker.to.vcf <- function(geno.temp){
    geno2 <- geno.temp[-c(1:5)]
    geno2[geno2 == "AA"] = '0/0'
    geno2[geno2 == "BB"] = '1/1'
    geno2[geno2 == "AB"] = '1/0'
    geno2[is.na(geno2)] = './.' # missing values are coded as 'NN'
    return(geno2)
  }
  convert.geno <- t(apply(input.geno, 1, recode.marker.to.vcf))
  vcf.format <- cbind(input.geno[,c(1:5)], convert.geno) %>%
    as_tibble() %>%
    rename(ID = SNPID, CHROM = Chrom, POS = Position) %>%
    mutate(QUAL = '.', FILTER = 'PASS', INFO = '.', FORMAT = 'GT') %>%
    dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, everything())
  header <- c("##fileformat=VCFv4.0",
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
              "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">",
              "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">",
              "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">",
              "##FORMAT=<ID=PL,Number=.,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">",
              "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">",
              "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
              "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">")
  write_lines(header, outfile, sep = "\n", append = F)
  line1 <- str_c(colnames(vcf.format), collapse = '\t') %>%
    str_c('#', ., "")
  write_lines(line1, outfile, append = T)
  write_tsv(vcf.format, outfile, append = T)
  return(vcf.format)
}

######

getProcessedData_NUST_withFilters <- function(gt2d,NUST_BLUEs,NUST_Meta_Table){

	NAIndices <- apply(gt2d[,6:ncol(gt2d)],2,function(x) which(is.na(x)))
	b <- unlist(lapply(NAIndices,length))
	Non_ZeroIndices <- which(b!=0)
	markerID <- as.vector(unlist(gt2d[,1]))



	NUST_Genotypes_VCF_gt2d <- gt2d[,-c(1:5)]
	NUST_Genotypes_VCF <- NUST_Genotypes_VCF_gt2d[,-Non_ZeroIndices]
	NUST_Meta_Table_StrainID <- NUST_Meta_Table[,1]
	NUST_Genotypes_VCF_ID <- colnames(NUST_Genotypes_VCF)


### Remove special characters from strain ID in meta table and genotype table

	NUST_Meta_Table_StrainID_Mod <- gsub("-","",NUST_Meta_Table_StrainID) 
	NUST_Meta_Table_StrainID_Mod1 <- gsub("\\.","",NUST_Meta_Table_StrainID_Mod) 
	NUST_Meta_Table_StrainID_Mod2 <- gsub("\\(","",NUST_Meta_Table_StrainID_Mod1) 
	NUST_Meta_Table_StrainID_Mod3 <- gsub("\\)","",NUST_Meta_Table_StrainID_Mod2) 
	NUST_Meta_Table_StrainID_Mod4 <- gsub("\\_","",NUST_Meta_Table_StrainID_Mod3) 


	NUST_Genotypes_VCF_ID <- gsub("-","",NUST_Genotypes_VCF_ID) 
	NUST_Genotypes_VCF_ID1 <- gsub("\\.","",NUST_Genotypes_VCF_ID) 
	NUST_Genotypes_VCF_ID2 <- gsub("\\(","",NUST_Genotypes_VCF_ID1) 
	NUST_Genotypes_VCF_ID3 <- gsub("\\)","",NUST_Genotypes_VCF_ID2) 
	NUST_Genotypes_VCF_ID4 <- gsub("\\_","",NUST_Genotypes_VCF_ID3)

### Checks 1

	length(which(NUST_Meta_Table_StrainID_Mod4 %in%  NUST_Genotypes_VCF_ID4))
	length(NUST_Genotypes_VCF_ID4)
	length(unique(NUST_Meta_Table_StrainID_Mod4))
	length(unique(NUST_Genotypes_VCF_ID4))

### Remove duplicated IDs from NUST genotype table

	NUST_Meta_Table_Mod <- NUST_Meta_Table
	NUST_Meta_Table_Mod[,1] <- NUST_Meta_Table_StrainID_Mod4
	NUST_Genotypes_Table_Mod <- cbind(NUST_Genotypes_VCF_ID4,t(NUST_Genotypes_VCF)) 

	colnames(NUST_Genotypes_Table_Mod)[1] <- "Strain" 
	colnames(NUST_Meta_Table_Mod)[1] <- "Strain"

	NUST_Meta_Table_Mod2 <- NUST_Meta_Table_Mod[which(!duplicated(NUST_Meta_Table_StrainID_Mod4)),]

### Checks 2 (Equal length)

	length(NUST_Meta_Table_Mod2[,1]) 
	length(unique(NUST_Meta_Table_Mod2[,1]))

## NUST merged table comprising meta data and genotypes table  


	NUST_Merged_Table <- merge(NUST_Meta_Table_Mod2,NUST_Genotypes_Table_Mod,by="Strain")
	init <- ncol(NUST_Meta_Table_Mod2)+1
	final <- ncol(NUST_Merged_Table)
	colnames(NUST_Merged_Table)[init:final] <- markerID

	########## IDs that are persent in genotypes table and not in the merged table 
	diffIDs <- setdiff(NUST_Genotypes_Table_Mod[,1],NUST_Merged_Table[,1])
	diffIndices <- which(NUST_Genotypes_Table_Mod[,1] %in% diffIDs)

	### Numeric coded genotype table from merged data table 


	NUST_Genotypes_Table_Mod_Merged <- NUST_Merged_Table[,c(1,init:final)]
	NUST_Genotypes_Table_Mod_Num1 <- apply(NUST_Genotypes_Table_Mod_Merged[,-1],2,function(x) gsub("BB","-1",x)) 
	NUST_Genotypes_Table_Mod_Num2 <- apply(NUST_Genotypes_Table_Mod_Num1,2,function(x) gsub("AB","0",x)) 
	NUST_Genotypes_Table_Mod_Num3 <- apply(NUST_Genotypes_Table_Mod_Num2,2,function(x) gsub("AA","1",x)) 
	NUST_Genotypes_Table_Mod_Num <- apply(NUST_Genotypes_Table_Mod_Num3,2,function(x) as.numeric(x)+1)

	### Merged Numeric Table

	NUST_Genotypes_Table_Mod_Num_Comb <- cbind(NUST_Genotypes_Table_Mod_Merged[,1],NUST_Genotypes_Table_Mod_Num)
	colnames(NUST_Genotypes_Table_Mod_Num_Comb)[1] <- "Strain" 

	NUST_Merged_Table_Num <- merge(NUST_Meta_Table_Mod2,NUST_Genotypes_Table_Mod_Num_Comb,by="Strain")
	dim(NUST_Merged_Table_Num)
	colnames(NUST_Merged_Table_Num)
	
	
# Filter Genotype Table based on year and Prog_Breeder_Factor

	dupIndices <- which(duplicated(as.character(NUST_Merged_Table[,1])))
	dupStrain <- as.character(NUST_Merged_Table[dupIndices,1])

	dupIndices_in_Table <- which(as.character(NUST_Merged_Table[,1]) %in% dupStrain)

	NUST_Merged_Table_noDup <- NUST_Merged_Table[-dupIndices,]

	rownames(NUST_Merged_Table_noDup) <- NUST_Merged_Table_noDup[,1]
	year_filt_indices <- (which(as.numeric(as.character(NUST_Merged_Table_noDup[,7])) >= 2010))
	NUST_Merged_Table_Filt1 <- NUST_Merged_Table_noDup[year_filt_indices,] 

	dim(NUST_Merged_Table_noDup)
	dim(NUST_Merged_Table_Filt1)


	dupIndices_Num <- which(duplicated(as.character(NUST_Merged_Table_Num[,1])))
	dupStrain_Num <- as.character(NUST_Merged_Table_Num[dupIndices,1])
	dupIndices_in_Table_Num <- which(as.character(NUST_Merged_Table_Num[,1]) %in% dupStrain_Num)

	NUST_Merged_Table_Num_noDup <- NUST_Merged_Table_Num[-dupIndices_Num,]
	rownames(NUST_Merged_Table_Num_noDup) <- NUST_Merged_Table_Num_noDup[,1]
	year_filt_indices <- (which(as.numeric(as.character(NUST_Merged_Table_Num_noDup[,7])) >= 2010))
	NUST_Merged_Table_Num_Filt1 <- NUST_Merged_Table_Num_noDup[year_filt_indices,] 


### Apply Filter(min program size 50 to reduce the number of programs to 12) 
####

	minProgramSize <- 50
	BrProg_col <- 9
	NUST_Merged_Table_noDup[,BrProg_col] <- gsub("_","-",NUST_Merged_Table_noDup[,BrProg_col])
	BrProg_factor <- factor(NUST_Merged_Table[,BrProg_col]) #,labels=c(1:length(levels(factor(NUST_Merged_Table_Filt1[,9])))))
	BrProg_factor_rm <- which(table(BrProg_factor)< minProgramSize)
	BrProg_filt_name <- names(table(BrProg_factor))[-BrProg_factor_rm]
	BrProg_factor_filt_indices <- which(as.character(BrProg_factor) %in% BrProg_filt_name)

	BrProg_factor_filt <- as.character(BrProg_factor)[BrProg_factor_filt_indices]
	BrProg_factor_filt_table <- table(as.character(BrProg_factor)[BrProg_factor_filt_indices])
	length(BrProg_factor_filt_table)
	BrProg_filt <- BrProg_factor[BrProg_factor_filt_indices]



	MG_col <- 10 
	MG_factor <- factor(NUST_Merged_Table_noDup[,MG_col])
	MG_filt <- MG_factor[BrProg_factor_filt_indices]

### Filtered Genotype Table 

	init <- ncol(NUST_Meta_Table_Mod2)+1
	final <- ncol(NUST_Merged_Table_noDup)
	NUST_Genotypes_Table_Mod_Filt1 <- NUST_Merged_Table_noDup[BrProg_factor_filt_indices,init:final]


    NUST_Genotypes_Table_Mod_Num_Filt1 <- NUST_Merged_Table_Num[BrProg_factor_filt_indices,init:final]

    rownames(NUST_Genotypes_Table_Mod_Filt1) <- rownames(NUST_Merged_Table_noDup[BrProg_factor_filt_indices,init:final])
    rownames(NUST_Genotypes_Table_Mod_Num_Filt1) <- rownames(NUST_Merged_Table_Num_noDup[BrProg_factor_filt_indices,init:final])

### Filter for markers with genetic map information

	# markers_in_GenoTable <- colnames(NUST_Merged_Table_noDup[BrProg_factor_filt_indices,init:final]) 
	# markerFilter <- which(markers_in_GenoTable %in% markers_with_map)

	# length(markers_in_GenoTable[markerFilter])
	# length(which(markers_with_map %in% markers_in_GenoTable))
	# length(which(is.na(markers_with_map[which(markers_with_map %in% markers_in_GenoTable)])))

	NUST_Genotypes_Table_Mod_Filt <- as.matrix(NUST_Genotypes_Table_Mod_Filt1)
	NUST_Genotypes_Table_Mod_Num_Filt <- apply(as.matrix(NUST_Genotypes_Table_Mod_Num_Filt1),2,as.numeric)

	rownames(NUST_Genotypes_Table_Mod_Num_Filt) <- rownames(NUST_Genotypes_Table_Mod_Num_Filt1)
	
	dim(NUST_Genotypes_Table_Mod_Filt)
	dim(NUST_Genotypes_Table_Mod_Num_Filt)

######### Data Prep for GP model training

 StrainID_List <- strsplit(as.character(NUST_BLUEs[,1]),"_")

 StrainID <- unlist(lapply(StrainID_List,function(x) paste(x[[1]],x[[2]],sep="")))
 StrainID1 <- gsub("MGIII","",StrainID) 
 StrainID2 <- gsub("MGII","",StrainID1) 
 StrainID3 <- gsub("MGI","",StrainID2) 
 StrainID4 <- gsub("MG0","",StrainID3) 
 StrainID5 <- gsub("MG0","",StrainID4)
 StrainID6 <- gsub(" ","",StrainID5)
 StrainIDMod <- StrainID6
 #rownames(NUST_BLUEs) <-  StrainIDMod 
 length(which(StrainIDMod %in% rownames(NUST_Genotypes_Table_Mod_Filt)))
 
 commonStrainID <- StrainIDMod[(which(StrainIDMod %in% rownames(NUST_Genotypes_Table_Mod_Filt)))]
 
 Diff_StrainID <- setdiff(rownames(NUST_Genotypes_Table_Mod_Filt),commonStrainID)

 NUST_GenoTable_Filtered <- cbind(rownames(NUST_Genotypes_Table_Mod_Filt),NUST_Genotypes_Table_Mod_Filt)
 colnames(NUST_GenoTable_Filtered)[1] <- "StrainID" 
 
 
 NUST_BLUEs_Filt <- NUST_BLUEs[(which(StrainIDMod %in% rownames(NUST_Genotypes_Table_Mod_Filt))),]
 NUST_BLUEs_Filt[,1] <-  commonStrainID  
 colnames(NUST_BLUEs_Filt)[1] <- "StrainID" 
 
 NUST_Data_Table1 <- merge(NUST_GenoTable_Filtered,NUST_BLUEs_Filt,by="StrainID")
 
 dupIndices_Data <- which(duplicated(NUST_Data_Table1[,"StrainID"]))
 NUST_Data_Table <- NUST_Data_Table1[-dupIndices_Data,]
 
 NUST_GenoTable_Num_Filtered <- cbind(rownames(NUST_Genotypes_Table_Mod_Num_Filt),NUST_Genotypes_Table_Mod_Num_Filt)
 colnames(NUST_GenoTable_Num_Filtered)[1] <- "StrainID" 
 
   
####### Filtered Table in Numeric Format 
 
 NUST_Data_Table1_Num <- merge(NUST_GenoTable_Num_Filtered,NUST_BLUEs_Filt,by="StrainID")
 dupIndices_Data <- which(duplicated(NUST_Data_Table1_Num[,"StrainID"]))
 NUST_Data_Table_Num <- NUST_Data_Table1_Num[-dupIndices_Data,]

 NAIndices <- which(is.na(NUST_Data_Table_Num[,"YieldBuA"]))
 NUST_Data_Table_Num_Filt <- NUST_Data_Table_Num[-NAIndices,]
    
	
return(NUST_Data_Table_Num_Filt)

}




getMergedData <- function(gt2d,Pheno,testIDs){

	NAIndices <- apply(gt2d[,6:ncol(gt2d)],2,function(x) which(is.na(x)))
	b <- unlist(lapply(NAIndices,length))
	Non_ZeroIndices <- which(b!=0)
	markerID <- as.vector(unlist(gt2d[,1]))

    #markerID_Filt <- markerID[-Non_ZeroIndices]
    #Genotypes_VCF_gt2d <- gt2d[,-c(1:5)]
	
	
	Non_ZeroIndices <- NULL
	Genotypes_VCF <-   gt2d[,-c(1:5,Non_ZeroIndices)]
	
	
	Genotypes_VCF_ID <- colnames(Genotypes_VCF)



### Remove special characters from strain ID in meta table and genotype table

	Genotypes_VCF_ID <- gsub("-","",Genotypes_VCF_ID) 
	Genotypes_VCF_ID1 <- gsub("\\.","",Genotypes_VCF_ID) 
	Genotypes_VCF_ID2 <- gsub("\\(","",Genotypes_VCF_ID1) 
	Genotypes_VCF_ID3 <- gsub("\\)","",Genotypes_VCF_ID2) 
	Genotypes_VCF_ID4 <- gsub("\\_","",Genotypes_VCF_ID3)

### Checks 1

  # TestIDs <- testIDs # TestIDsMod <- gsub("[_-]","",TestIDs)

    
	TestIDs <- gsub("[_-]","",testIDs)

	length(Genotypes_VCF_ID4)
	length(unique(Genotypes_VCF_ID4))

### Remove duplicated IDs from genotype table


	Genotypes_Table_Mod <- cbind(Genotypes_VCF_ID4,t(Genotypes_VCF)) 

	colnames(Genotypes_Table_Mod) <- c("Strain",markerID)

### Numeric coded genotype table from merged data table 

	Genotypes_Table_Mod_Merged <- Genotypes_Table_Mod
	Genotypes_Table_Mod_Num1 <- apply(Genotypes_Table_Mod_Merged[,-1],2,function(x) gsub("BB","-1",x)) 
	Genotypes_Table_Mod_Num2 <- apply(Genotypes_Table_Mod_Num1,2,function(x) gsub("AB","0",x)) 
	Genotypes_Table_Mod_Num3 <- apply(Genotypes_Table_Mod_Num2,2,function(x) gsub("AA","1",x)) 
	
	
	if(length(which(rownames(Genotypes_Table_Mod_Num3) %in% "SNPID")) >=1){
	Genotypes_Table_Mod_Num3A <- Genotypes_Table_Mod_Num3[-which(rownames(Genotypes_Table_Mod_Num3) %in% "SNPID"),]
    }
	if(length(which(rownames(Genotypes_Table_Mod_Num3) %in% "SNPID")) <1){
	Genotypes_Table_Mod_Num3A <- Genotypes_Table_Mod_Num3
    }
	
	Genotypes_Table_Mod_Num <- apply(Genotypes_Table_Mod_Num3A,2,function(x) as.numeric(x)+1)


	Genotypes_Table_Mod_Num_Comb <- cbind(Genotypes_Table_Mod_Merged[,1],Genotypes_Table_Mod_Num)
	colnames(Genotypes_Table_Mod_Num_Comb)[1] <- "Strain" 

	
	
	Genotype_Table_Num <- Genotypes_Table_Mod_Num_Comb
	dim(Genotype_Table_Num)
	colnames(Genotype_Table_Num)
	
# Filter Genotype Table and remove duplicate IDs
	
	dupIndices_Num <- which(duplicated(as.character(Genotype_Table_Num[,1])))
	dupStrain_Num <- as.character(Genotype_Table_Num[dupIndices_Num,1])
	dupIndices_in_Table_Num <- which(as.character(Genotype_Table_Num[,1]) %in% dupStrain_Num)

   if(length(dupIndices_Num) >=1) {
	Genotype_Table_Num_noDup <- Genotype_Table_Num[-dupIndices_Num,]
	rownames(Genotype_Table_Num_noDup) <- Genotype_Table_Num_noDup[,1]
	Genotype_Table_Num_Filt1 <- Genotype_Table_Num_noDup
	
	
	Genotypes_Table_Mod_Num_Filt <- apply(as.matrix(Genotype_Table_Num_Filt1[,-1]),2,as.numeric)
	rownames(Genotypes_Table_Mod_Num_Filt) <- rownames(Genotype_Table_Num_Filt1)
	dim(Genotypes_Table_Mod_Num_Filt)
	Genotypes_Table_Mod_Num_Filt <- cbind(rownames(Genotypes_Table_Mod_Num_Filt),Genotypes_Table_Mod_Num_Filt)
	
   }
   if(length(dupIndices_Num) <1) {
	Genotype_Table_Num_noDup <- Genotype_Table_Num
	rownames(Genotype_Table_Num_noDup) <- Genotype_Table_Num_noDup[,1]
	Genotype_Table_Num_Filt1 <- Genotype_Table_Num_noDup
		
	Genotypes_Table_Mod_Num_Filt <- apply(as.matrix(Genotype_Table_Num_Filt1[,-1]),2,as.numeric)
	rownames(Genotypes_Table_Mod_Num_Filt) <- rownames(Genotype_Table_Num_Filt1)
	dim(Genotypes_Table_Mod_Num_Filt)
	Genotypes_Table_Mod_Num_Filt <- cbind(rownames(Genotypes_Table_Mod_Num_Filt),Genotypes_Table_Mod_Num_Filt)
   }

### Filtered Genotype Table 
	
	
#### 
    if(!is.null(TestIDs)){
	  testIndices <- which(as.character(Genotypes_Table_Mod_Num_Filt[,1]) %in% TestIDs)
	  StrainIDs <- as.character(Genotypes_Table_Mod_Num_Filt[,1])
	  TrainIDs <- setdiff(StrainIDs,TestIDs)
	  trainIndices <- which(as.character(Genotypes_Table_Mod_Num_Filt[,1]) %in% TrainIDs)
	  Test_Genotypes_Table_Mod_Num_Filt <- Genotypes_Table_Mod_Num_Filt[testIndices,]
      Train_Genotypes_Table_Mod_Num_Filt <- Genotypes_Table_Mod_Num_Filt[trainIndices,]
    }
  
   if(is.null(TestIDs)){ 
     
	  print("Load Target File")
   }
  
######### Data Prep for GP model training

############# Process IDs
### PhenoIDs  

	StrainID_List <- strsplit(as.character(Pheno[,1]),"[-_]")

	StrainID <- unlist(lapply(StrainID_List,function(x) paste(x,collapse="")))
	 
	if(length(grep("MG",as.character(StrainID))) >1){
		  
		  StrainIDMod <- gsub("MG.*","",as.character(StrainID))
		  length(which(StrainIDMod %in% rownames(Train_Genotypes_Table_Mod_Num_Filt)))
	}

    if(length(grep("MG",as.character(StrainID))) <1){
			StrainIDMod <- as.character(StrainID)
    }		 
	 
 
    # StrainID1 <- gsub("MGIII","",StrainID) 
	# StrainID2 <- gsub("MGII","",StrainID1) 
	# StrainID3 <- gsub("MGI","",StrainID2) 
	# StrainID4 <- gsub("MG0","",StrainID3) 
	# StrainID5 <- gsub("MG00","",StrainID4)
	# StrainID6 <- gsub(" ","",StrainID5)
	# StrainIDMod <- StrainID6

	
### GenoIDs  
  
  trainStrainID <- gsub("[-_]","",(Train_Genotypes_Table_Mod_Num_Filt)[,1])

  commonStrainID <- StrainIDMod[(which(StrainIDMod %in% trainStrainID))]
  Diff_StrainID <- setdiff(trainStrainID,commonStrainID)

  GenoTable_Filtered <- cbind(trainStrainID,Train_Genotypes_Table_Mod_Num_Filt)
  colnames(GenoTable_Filtered)[1] <- "StrainID" 
  
  PhenoTable_Filtered <- Pheno[(which(StrainIDMod %in% trainStrainID)),]
  PhenoTable_Filtered[,1] <-  commonStrainID  
  colnames(PhenoTable_Filtered)[1] <- "StrainID" 
  
### merge geno and pheno tables  
  #Data_Table1_Num <- merge(GenoTable_Filtered,PhenoTable_Filtered,by="StrainID")
  
  Data_Table1_Num <- merge(GenoTable_Filtered,PhenoTable_Filtered,by="StrainID",all=TRUE)
 
 
### Dealing with duplicated data 
 
  dupId_Data <- Data_Table1_Num[which(duplicated(Data_Table1_Num[,"StrainID"])),"StrainID"]
  dupId_Indices <- lapply(dupId_Data,function(x) which(as.character(Data_Table1_Num[,"StrainID"]) %in% as.character(x)))
  dupId_Indices_Len <- lapply(dupId_Indices,length)
  
  table(unlist(dupId_Indices_Len)) 
  dupData_List <- lapply(dupId_Indices,function(x) Data_Table1_Num[as.vector(x),])
  dupData_NA_List <- lapply(dupData_List,function(x) apply(x,1,function(y) length(which(is.na(y)))))
  

 


## Check difference of genotypic scores 
 
  # findDiff <- function(x){
      # d <- matrix(0,nrow=nrow(x),ncol=ncol(x))
	  # for(i in 2:nrow(x)){
	     # ind <- which(!is.na(x[1,1:1205]))
	     # d[i-1,ind] <- (as.numeric(as.character(x[1,ind]))-as.numeric(as.character(x[i,ind])))
	  # }
	  # return(d)
  # }
  
  # dupData_Diff_List <- lapply(dupData_List,function(x) findDiff(x))
 
 ### Simple remove 
  # Data_Table_Num <- Data_Table1_Num[-dupIndices_Data,]
  # rownames(Data_Table_Num) <- Data_Table1_Num[-dupIndices_Data,"StrainID"]
  
  ## Checks 
   # CommonSampleID_Mod <- (which(as.character(Data_Table1_Num[,1]) %in% trainStrainID))
   # length(unique(Data_Table1_Num[CommonSampleID_Mod,1]))
   # length(which(duplicated(Data_Table1_Num[CommonSampleID_Mod,1])))
   # commonGenoMod <- (Data_Table1_Num[(which(duplicated(Data_Table1_Num[CommonSampleID_Mod,1]))),])
   # apply(commonGenoMod,1,function(x) length(which(is.na(as.character(x)))))
   
  
  
  
####### Filtered Table in Numeric Format 
####################################################################################
   #Data_Table_Num <- Data_Table1_Num
   #rownames(Data_Table_Num) <- Data_Table1_Num[,"StrainID"]
  
   Train_Data_Table_Num_Filt <- Data_Table1_Num 
  	
   return(list(Train_Data_Table_Num_Filt,Test_Genotypes_Table_Mod_Num_Filt))

}


getProcessedData <- function(Data_Table_Num_List,trait){

     TrainData_Table_Num <- Data_Table_Num_List[[1]]
	 TestData_Table_Num_Filt<- Data_Table_Num_List[[2]]
	 
	 ### Remove lines with 'NA' for trait values
	 if(length(trait)==1){
	  NAIndices <-  which(is.na(TrainData_Table_Num[,trait]))
	  if(length(NAIndices)>1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num[-NAIndices,]
	  }
	  if(length(NAIndices)<1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num
	  }
	 }
	 
	 if(length(trait)>1){
	   NAIndices <- c(unlist(apply(TrainData_Table_Num[,trait],2,function(x)which(is.na(x)))))
	   if(length(NAIndices)>1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num[-NAIndices,]
	   }
	   if(length(NAIndices)<1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num
	   }
      
	  }
     return(list(Train_Data_Table_Num_Filt,TestData_Table_Num_Filt)) 
	 
} 





getPredictionData <- function(Data_Table_Num_List,noCandidates){

     TrainData_Table_Num_Filt <- Data_Table_Num_List[[1]]
	 TestData_Table_Num_Filt<- Data_Table_Num_List[[2]]
	 set.seed(125)
	 CandidateIndices <- sample(c(1:nrow(TrainData_Table_Num_Filt)),noCandidates)
	 Candidates <- as.character(TrainData_Table_Num_Filt[CandidateIndices,1])
	 
	 if(length(Candidates)==nrow(TrainData_Table_Num_Filt)){
	  
      Candidate_Data_Table_Num_Filt <- TrainData_Table_Num_Filt
	 }
	 if(length(Candidates)< nrow(TrainData_Table_Num_Filt)){
	  
	  Candidate_Data_Table_Num_Filt <- TrainData_Table_Num_Filt[CandidateIndices,]
      
	 }
     
     return(list(Candidate_Data_Table_Num_Filt,TestData_Table_Num_Filt)) 
	 
} 


### Step 2 : Get predicted genetic values usin kin.blup 

 getRankedPredictedValues <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModel,optTS=NULL){ 
    
   fixedX <- NULL
   
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	 
	
	 TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	
	 if(is.null(optTS)){
	    
		 trainPheno0 <- as.numeric(as.character(TrainData_Table_Num_Filt[,trait]))
		 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 }
	 
	 if(!is.null(optTS)){
	     strainID <- as.character(TrainData_Table_Num_Filt[,1])
		 trainSetID <- as.character(as.vector(optTS))
	     trainIndices <- which(strainID %in% trainSetID)
		 trainPheno0 <- as.numeric(as.character(TrainData_Table_Num_Filt[trainIndices,trait]))
		 geno_012 <- apply(TrainData_Table_Num_Filt[trainIndices,c(initCol:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[trainIndices,1])
	 } 
	 
	 
	 if(anyNA(trainGeno)){	 
		trainGeno_Imp0 <- snpQC(trainGeno,impute=TRUE,remove=FALSE)
	 }else{trainGeno_Imp0 <- trainGeno}
	 
	 
	  if(anyNA(TestGenoTable)){ 
	    testGeno_Imp0 <- snpQC(TestGenoTable,impute=TRUE,remove=FALSE)
		testGeno_Imp <- apply(testGeno_Imp0,2,function(x) as.numeric(x-1))
		testNA <- apply(testGeno_Imp,2,function(x) length(which(is.na(x))))
		testNAIndices <-(which(unlist(testNA) !=0))
	 }else{testGeno_Imp <- apply(TestGenoTable,2,function(x) as.numeric(x-1)) 
	    testNAIndices <- NULL
	  }
	  
	  
	 
	 if(anyNA(trainPheno0)){
	   ph_NA_Indices <- which(is.na(trainPheno0))
	   trainPheno <- trainPheno0[-ph_NA_Indices] 
	   trainGeno_Imp <- trainGeno_Imp0[-ph_NA_Indices,]
	   Geno <- as.character(TrainData_Table_Num_Filt[-ph_NA_Indices,1])
	   
	 }
	 if(!anyNA(trainPheno0)){
	  
	   trainPheno <- trainPheno0
	   trainGeno_Imp <- trainGeno_Imp0
	   Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 } 

## Kinship Matrix
	 A <- A.mat(trainGeno_Imp)	
	 colnames(A) <- Geno
	 rownames(A) <- Geno

## Prepare Data Table for GP 
	 Data <- cbind.data.frame(Geno,trainPheno)
	 colnames(Data) <- c("Geno","Pheno")
	 Geno <- "Geno"
	 Pheno <- "Pheno"
	
#### Impute trainGeno and train using mixed.solve
	
	
	 
	if(GPModel == "rrBLUP (rrBLUP)"){ 
	 if(!is.null(fixedX)){
	   pred <- mixed.solve(trainPheno,Z=trainGeno_Imp,X=fixedX,SE=FALSE,return.Hinv =FALSE) 
	   Mean <- fixedX %*% pred$beta
	 }
	 if(is.null(fixedX)){
	     pred <- mixed.solve(trainPheno,Z=trainGeno_Imp,SE=FALSE,return.Hinv =FALSE) 
	     Mean <- as.numeric(pred$beta)
	  }
	 Effects <- pred$u
   PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
   SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
 
	# Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
	 if(!is.null(testNAIndices)){
		 Test_PredictedValues <-  Mean + (testGeno_Imp[,-testNAIndices] %*% Effects[-testNAIndices])
	 }
	 if(is.null(testNAIndices)){
		 Test_PredictedValues <- Mean+ (testGeno_Imp %*% Effects)
	 }
	 
	 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
	 Test_StrainID <- rownames(TestData_Table_Num_Filt)
    
	 }
	 
	 
	 if(GPModel == "rrBLUP (bWGR)"){
	 
		 pred <- emRR(trainPheno,trainGeno_Imp) 
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b
		 
		 
		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
		 
	    if(!is.null(testNAIndices)){
		 Test_PredictedValues <-  Mean + (testGeno_Imp[,-testNAIndices] %*% Effects[-testNAIndices])
		}
		if(is.null(testNAIndices)){
		 Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
		}
		
		 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
		 Test_StrainID <- rownames(TestData_Table_Num_Filt)
    
	 }
	 
	 
	  if(GPModel == "BayesB (bWGR)"){ 
	 
		 pred <- emBB(trainPheno,trainGeno_Imp) 
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b
		 
		 
		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
		 
			
		 # Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
		 
		 if(!is.null(testNAIndices)){
		   Test_PredictedValues <-  Mean + (testGeno_Imp[,-testNAIndices] %*% Effects[-testNAIndices])
		 }
	     if(is.null(testNAIndices)){
		   Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
	     }
		 
		 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
		 Test_StrainID <- rownames(TestData_Table_Num_Filt)
      }
	  
	  if(GPModel == "BayesLASSO (bWGR)"){ 
	 
		 pred <- emBL(trainPheno,trainGeno_Imp) 
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b
		 	 
		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
		 
			
		 #Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
		 
		 if(!is.null(testNAIndices)){
		   Test_PredictedValues <-  Mean + (testGeno_Imp[,-testNAIndices] %*% Effects[-testNAIndices])
		 }
	     if(is.null(testNAIndices)){
		   Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
	     }
		 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
		 Test_StrainID <- rownames(TestData_Table_Num_Filt)
      }
	  
### 
	# pred.kb <- kin.blup(as.data.frame(Data),Geno,Pheno,GAUSS=FALSE,K=A,covariate=NULL,PEV=TRUE,n.core=1,theta.seq=NULL)
		 
### Upper bound of Reliability
    if(length(testNAIndices)>0){
       trainGeno_Imp2 <- apply(trainGeno_Imp[,-testNAIndices],2,function(x) x+1)  
	}
	if(length(testNAIndices)==0){
	    trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) x+1)  
	}
    cleanData <- cleanREP(trainPheno,trainGeno_Imp2)
    M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M
      
    getU <- function(M.Pdt,v){
	   v.hat <- M.Pdt %*% v
	   U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   return(U)
    }
	
	if(length(testNAIndices)>0){
      U <- apply(testGeno_Imp[,-testNAIndices],1,function(x) getU(M.Pdt,x)) 
	}
	if(length(testNAIndices)==0){
	  U <- apply(testGeno_Imp,1,function(x) getU(M.Pdt,x)) 
	}
	
	Test_SortedIndices <- Test_SortedPredictedValues[[2]]
	Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
	U_Sorted <- U[(Test_SortedIndices)]
	
	outputDF <- cbind.data.frame(Sorted_Test_StrainID ,round(Test_SortedPredictedValues[[1]],digits=2),round(U_Sorted,digits=5))
	colnames(outputDF) <- c("LineID",paste("Predicted Value for ",trait,sep=""),"Upper Bound of Reliability")
	
	
   return(outputDF)
	
 }
	
#### 


 getRankedPredictedValuesMT <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModelMT,optTS=NULL){ 
   
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	 			 
	 TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
	 
	
	 rownames(TestGenoTable) <- TestData_Table_Num_Filt[,1]
	 
	 set.seed(125)
	 strainID <- as.character(TrainData_Table_Num_Filt[,1])
	 
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	 
	 # trainSetID <- as.character(as.vector( strainID[sample(c(1:length(strainID)),500)]))
	 # optTS <- cbind(trainSetID,c(1:500))
	 
	 trainSetID <- as.character(as.vector(strainID))
	 optTS <- cbind(trainSetID,c(1:length(trainSetID)))
	
	 if(is.null(optTS)){
	
		 trainPheno <- TrainData_Table_Num_Filt[,trait]
		 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 #Geno <- rownames(TrainData_Table_Num_Filt)
		 Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 }
	 
	 if(!is.null(optTS)){
	     
	     strainID <- as.character(TrainData_Table_Num_Filt[,1])
		 trainSetID <- as.character(as.vector(optTS[,1]))
	     trainIndices <- which(strainID %in% trainSetID)
		 trainPheno <- TrainData_Table_Num_Filt[trainIndices,trait]
		 geno_012 <- apply(TrainData_Table_Num_Filt[trainIndices,c(initCol:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[trainIndices,1])
	 }
	  
## Prepare Data Table for GP 
	 # Data <- cbind.data.frame(Geno,trainPheno)
	 # colnames(Data) <- c("Geno","Pheno")
	 # Geno <- "Geno"
	 # Pheno <- "Pheno"
	

#### Impute trainGeno and train using mixed.solve
	
	if(anyNA(trainGeno)){	 
		trainGeno_Imp <- snpQC(trainGeno,impute=TRUE,remove=FALSE)
	}else{trainGeno_Imp <- trainGeno}
	 
## Kinship Matrix 
    # rownames(trainGeno_Imp) <- Geno
	 #rownames(trainPheno) <- Geno
	 A <- A.mat(trainGeno)	
	 colnames(A) <- Geno
	 rownames(A) <- Geno

	 Y <- trainPheno
	 X <- rbind(trainGeno_Imp,TestGenoTable)
	 #rownames(X) <- c(rownames(trainGeno_Imp),rownames(TestGenoTable))
	 n <- nrow(X)
	 p <- ncol(X)
	 
	 nTst <- nrow(TestGenoTable)
	 Ytst <- matrix(rep(NA,nTst*ncol(Y)),nrow=nTst,ncol=ncol(Y))
	 colnames(Ytst) <- colnames(Y)
	 yNA <- as.matrix(rbind(Y,Ytst))
	 yNA_DF <- as.data.frame(yNA)
	 
	 #yNA_DF$id <- as.factor(rownames(X))
	 yNA_DF$id <- as.factor(paste(c(Geno,rownames(TestGenoTable)),"_",yNA_DF[,1],sep=""))
	 
	 if(GPModelMT == "BRR (BGLR)"){ 
				  
			ETA <- list(list(X=X,model="BRR"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }
	 if(GPModelMT == "RKHS (BGLR)"){ 
	        A.Tot <- A.mat(X)
			ETA <- list(list(K=A.Tot,model="RKHS"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }
	 if(GPModelMT == "Spike-Slab (BGLR)"){ 
				  
			ETA <- list(list(X=X,model="SpikeSlab"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }	
	 if(GPModelMT == "GBLUP (SOMMER)"){ 
	         
			A.Tot <- A.mat(X)
			# rownames(A.Tot) <- rownames(X) 
			# colnames(A.Tot) <- rownames(X)
			rownames(A.Tot) <- yNA_DF$id #c(Geno,rownames(TestGenoTable))
			colnames(A.Tot) <- yNA_DF$id #c(Geno,rownames(TestGenoTable))
			fm3 <- mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
            random=~vs(id,Gu=A.Tot),
            rcov=~units,MVM=TRUE,
            data=yNA_DF,verbose = TRUE)
	 }
	 tst <- c((length(trainSetID)+1):nrow(X))
	 if(GPModelMT != "GBLUP (SOMMER)"){ 
	  Test_PredictedValues <- fm3$ETAHat[tst,]
	 }
	 if(GPModelMT == "GBLUP (SOMMER)"){  
	  UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
	  Test_PredictedValues <- c()
	  for(nT in 1:length(trait)){
	    Test_PredictedValues <- cbind(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
	  }
	 }
	 
	Test_SortedPredictedValues <- sort.int(Test_PredictedValues[,1],decreasing=TRUE,index.return=TRUE) 
	Test_StrainID <- rownames(TestData_Table_Num_Filt)	 
### Upper bound of Reliability

    trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) x+1)  
    cleanData <- cleanREP(trainPheno,trainGeno_Imp2)
    M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M
      
    getU <- function(M.Pdt,v){
	   v.hat <- M.Pdt %*% v
	   U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   return(U)
    }
    U <- apply(TestGenoTable,1,function(x) getU(M.Pdt,x)) 

### Sort Output
	
	Test_SortedIndices <- Test_SortedPredictedValues[[2]]
	Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
	U_Sorted <- U[(Test_SortedIndices)]
	
	Test_SortedPredValuesTab <- apply(Test_PredictedValues[Test_SortedIndices,],2,function(x) round(x,digits=2))
	
### Output DF	
	outputDF <- cbind.data.frame(Sorted_Test_StrainID ,Test_SortedPredValuesTab,round(U_Sorted,digits=5))
	colnames(outputDF) <- c("LineID",colnames(Test_SortedPredValuesTab),"Upper Bound of Reliability")
	
	
   return(outputDF)
	
 }


  
	
  getemCVR <- function(Data_Table_Num_Filt_List,trait,nTraits,k,nIter){
  
     TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	 pheno0 <- as.numeric(as.character(TrainData_Table_Num_Filt[,trait]))
	 
	  
	 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
	 geno <- apply(geno_012,2,function(x) x-1)
	 Geno <- rownames(TrainData_Table_Num_Filt)
  
     if(anyNA(geno)){	 
		geno_Imp0 <- snpQC(geno,impute=TRUE,remove=FALSE)
	 }else{geno_Imp0 <- geno}
	 
	 if(anyNA(pheno0)){ 
	  ph_NA_Indices  <- which(is.na(pheno0))
	  pheno <- pheno0[-ph_NA_Indices] 
	  geno_Imp <- geno_Imp0[-ph_NA_Indices,] 
     } 
     if(!anyNA(pheno0)){ 
	   pheno <- pheno0
	   geno_Imp <- geno_Imp0
     }
	 
     emCVR <- emCV(pheno,geno_Imp,k,nIter)
  
     return(emCVR)
  }
  
 #################################
 
  
 
  getMTCVR <- function(Data_Table_Num_Filt_List,trait,nTraits,k,nIter){
  
     TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	
	 pheno_wNA <- as.numeric(as.character(TrainData_Table_Num_Filt[,trait]))
	 geno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
	 geno_wNA<- apply(geno_012,2,function(x) x-1)
	 Geno_wNA <- rownames(TrainData_Table_Num_Filt)
	 
	 NAIndices <- unlist(apply(pheno_wNA,2,function(x)which(is.na(x))))
  
     pheno <- pheno_wNA[-NAIndices,]
     geno <- geno_wNA[-NAIndices,]
	 Geno <- Geno_wNA[-NAIndices]
	 
     if(anyNA(geno)){	 
		geno_Imp <- snpQC(geno,impute=TRUE,remove=FALSE)
	 }else{geno_Imp <- geno}

     n <- nrow(geno_Imp)
##############################################################
	 
	 # nTr <- round(n*((k-1)/k),digits=0) 
	 # trainIndices <- sample(c(1:n),nTr)
	 # testIndices <- setdiff(c(1:n),trainIndices)
	  
	 # trPheno <- pheno[trainIndices,]
	 # trGeno <- geno_Imp[trainIndices,]
	 # testPheno <- pheno[testIndices,]
	 # testGeno <- geno_Imp[testIndices,]
	 
	 # trGenoNames <- Geno[trainIndices]
	 # tstGenoNames <- Geno[testIndices]
	 
############################################################

  GPModelMT_List <- list("BRR (BGLR)","RKHS (BGLR)","Spike-Slab (BGLR)") #,"GBLUP (SOMMER)")
  PA_Out_GP <- list()
  
  

  for(nGP in 1:length(GPModelMT_List)){ 
    
	GPModelMT <- GPModelMT_List[[nGP]]
	
	PA_List <- list()
	
	for(nrep in 2:nIter){
 
     nK <- floor(n/k)
	 k_List <- list()
	 set.seed(125+nrep) 
	 tot <- c(1:n)
	  for(i in 1:k){
	  	  
	   k_List[[i]] <- sample(tot,nK)
	   tot <- setdiff(tot,k_List[[i]]) 
	    
	  }

	  trIndices_List <- list()
	  tstIndices_List <- list()
      for(i in 1:k){ 
	 
	      trIndices_List[[i]] <- unlist(k_List[-i])
		  tstIndices_List[[i]] <- k_List[[i]]
	  }	  
		 
	PA <- list()
	 
	for(i in 1:k){
	
	 trainIndices <-  trIndices_List[[i]]
	 testIndices <- tstIndices_List[[i]]
	 	
#########################################################
     pheno <- pheno_wNA[-NAIndices,]
	 Y <- pheno
	 nTst <- length(testIndices)
	 
	 pheno[testIndices,] <- matrix(rep(NA,nTst*ncol(pheno)),nrow=nTst,ncol=ncol(pheno))
	 	 
## Kinship Matrix 
     rownames(geno_Imp) <- Geno
	 rownames(pheno) <- Geno
	
	 yNA <- as.matrix(pheno)
	 X <- geno_Imp
	 rownames(X) <- rownames(geno_Imp)
	 n <- nrow(X)
	 p <- ncol(X)
	 
	 yNA_DF <- as.data.frame(yNA)
	 yNA_DF$id <- as.factor(rownames(X))
	 
	 if(GPModelMT == "BRR (BGLR)"){
				  
			ETA <- list(list(X=X,model="BRR"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }
	
	 if(GPModelMT == "RKHS (BGLR)"){ 
	        A.Tot <- A.mat(X)
			ETA <- list(list(K=A.Tot,model="RKHS"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }
	 if(GPModelMT == "Spike-Slab (BGLR)"){ 
				  
			ETA <- list(list(X=X,model="SpikeSlab"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }	
	 if(GPModelMT == "GBLUP (SOMMER)"){ 
	         
			A.Tot <- A.mat(X)
			rownames(A.Tot) <- rownames(X) 
			colnames(A.Tot) <- rownames(X)
			fm3 <- mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
            random=~vs(id,Gu=A.Tot),
            rcov=~units,
            data=yNA_DF,verbose = TRUE)
	 }
	 
	 tst <- testIndices
	
	 if(GPModelMT != "GBLUP (SOMMER)"){ 
	  Test_PredictedValues <- fm3$ETAHat[tst,]
	 }
	 if(GPModelMT == "GBLUP (SOMMER)"){  
	  UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
	  Test_PredictedValues <- c()
	  for(nT in 1:length(trait)){
	    Test_PredictedValues <- cbind(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
	  }
	 }
	 	 
	  PA[[i]] <- diag(cor(Test_PredictedValues,Y[tst,]))
	 
	}
	  
	 PA_List[[nrep]] <- do.call(rbind,lapply(PA,function(x) x))
	  
	}
	PA_Out_GP[[nGP]] <- apply(do.call(rbind,lapply(PA_List,function(x)x)),2,mean)
	
  }
	
	
	PA_Out1 <-  apply(do.call(rbind,lapply(PA_Out_GP,function(x) x)),2,function(x) round(x,digits=2))
	
	
	PA_Out2 <- cbind(unlist(GPModelMT_List),PA_Out1)
	PA_Out <- rbind(c("GPModel",colnames(PA_Out1)),PA_Out2)
	return(PA_Out)
   
 }

####

 getOptimalTS <- function(Data_Table_Num_Filt_List,trait,nTraits,noCandidates,nTrainToSelect,GAParameters){
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	
	# NAIndices <-  which(is.na(TrainData_Table_Num[,trait]))
    # Train_Data_Table_Num_Filt <- NUST_Data_Table_Num[-NAIndices,]
		 
	
	trainPheno <- TrainData_Table_Num_Filt[,trait]
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(2:(ncol(TestData_Table_Num_Filt)))],2,as.numeric)
	
## Complete Genotypes table
	
	geno_012 <- rbind(trainGeno_012,testGeno_012)
    totGeno <-  apply(geno_012,2,function(x) x-1)
	rownames(totGeno) <- c(rownames(TrainData_Table_Num_Filt),rownames(TestData_Table_Num_Filt))
	if(length(trait)==1){
	names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
    }
	if(length(trait)>1){
	rownames(trainPheno) <- rownames(TrainData_Table_Num_Filt)
    }
	
    totCandidates <-  rownames(totGeno)
	
    #Test <- testIds
	Test <- rownames(TestData_Table_Num_Filt)
	testIndices <- which(totCandidates %in% Test)
 
### Set1 with 99 PCs, 100 candidates , 100 test and 50 nTrain 
    set.seed(125)
	
## Reduce or keep candidate genotypic data 
    
    nTrain <- nrow(TrainData_Table_Num_Filt)
	trainSetIndices <- c(1:nTrain)
	# if(noCandidates!= nTrain){
	  # trainSetIndices <-sample(c(1:nTrain),noCandidates)
	# }
	# if(noCandidates== nTrain){
	#}
	
	G <- rbind(totGeno[trainSetIndices,],totGeno[testIndices,])
	rownames(G) <- c(rownames(totGeno)[trainSetIndices],rownames(totGeno)[testIndices])
	Candidates <- rownames(totGeno)[trainSetIndices]
	
	G_Imp <- snpQC(G,remove=FALSE,impute=TRUE)
	rownames(G_Imp) <- rownames(G)
    GenoSVD <- svd(G_Imp,nu=99,nv=99)
    PC99 <- G%*%GenoSVD$v
    rownames(PC99)<-rownames(G_Imp)
    Train_STPGA <- c()
	
    system.time({
	
	  if(length(GAParameters$errorstat)==1){
         Train_STPGA <- GenAlgForSubsetSelection(P=PC99,Candidates,Test,ntoselect=nTrainToSelect,InitPop=GAParameters$InitPop,
         npop=GAParameters$npop, nelite=GAParameters$nelite, mutprob=GAParameters$mutprob, mutintensity = GAParameters$mutintensity,
         niterations=GAParameters$niterations,minitbefstop=GAParameters$minitbefstop, tabu=GAParameters$tabu,
         tabumemsize = GAParameters$tabumemsize,plotiters=GAParameters$plotiters,errorstat=GAParameters$errorstat,
		 lambda=GAParameters$lambda, mc.cores=GAParameters$mc.cores)
		 
	   }
	   if(length(GAParameters$errorstat)>1){
         Train_STPGA <- GenAlgForSubsetSelectionMO(P=PC99,Candidates,Test,ntoselect=nTrainToSelect,InitPop=GAParameters$InitPop,
         npop=GAParameters$npop, nelite=GAParameters$nelite, mutprob=GAParameters$mutprob, mutintensity = GAParameters$mutintensity,
         niterations=GAParameters$niterations,minitbefstop=GAParameters$minitbefstop, tabu=GAParameters$tabu,
         tabumemsize = GAParameters$tabumemsize,plotiters=GAParameters$plotiters,errorstat=GAParameters$errorstat,
		 lambda=GAParameters$lambda, mc.cores=GAParameters$mc.cores)
		 
	   }
	  
    })
	
	Train_STPGA_Indices <-  which(Candidates %in% as.character(Train_STPGA$`Solution with rank 1`))
    	
    return(list(Train_STPGA,Train_STPGA_Indices))
 }
 
 
 getRandomTS <-  function(Data_Table_Num_Filt_List,trait,nTraits,noCandidates,nTrainToSelect){ 
 
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	
	trainPheno <- TrainData_Table_Num_Filt[,trait]
	if(length(trait)==1){
	 names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
    }
	if(length(trait)>1){
	 rownames(trainPheno) <- rownames(TrainData_Table_Num_Filt)
    }
	
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(2:(ncol(TestData_Table_Num_Filt)))],2,as.numeric)
	
## Complete genotypic table of all candidates

    trainGeno <-  apply(trainGeno_012,2,function(x) x-1)
	rownames(trainGeno) <- rownames(TrainData_Table_Num_Filt)
	
	Test <- rownames(TestData_Table_Num_Filt)
	
### Set1 with 99 PCs, 100 candidates , 100 test and 50 nTrain 
   
	
## Reduce or keep Candidate Genotypic data 
    nTrain <- nrow(TrainData_Table_Num_Filt)
	set.seed(125)
	trainSetIndices <- c(1:nTrain)
	
	# if(noCandidates!= nTrain){
	# trainSetIndices <-sample(c(1:nTrain),noCandidates)
	# }
	# if(noCandidates== nTrain){
	#}
	
### Define candidate training set
	
    G_Train <- trainGeno[trainSetIndices,]
	Candidates<- rownames(trainGeno)[trainSetIndices]
		
## Impute G_Train 

	G_Train_Imp <- snpQC(G_Train,remove=FALSE,impute=TRUE)
	rownames(G_Train_Imp) <- rownames(G_Train)
	
	trainRandomIndices <- sample(c(1:nrow(G_Train_Imp)),nTrainToSelect) 
	Train_Random <- rownames(G_Train_Imp)[trainRandomIndices]
		
   return(list(Train_Random,trainRandomIndices))

  }
  
   
 getTSComparisons <- function(Data_Table_Num_Filt_List,Train_STPGA,Train_Random,trait,nTraits,testIds,optTS=NULL){ 
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	
	trainPheno <- TrainData_Table_Num_Filt[,trait]
	names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
   	
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(2:(ncol(TestData_Table_Num_Filt)))],2,as.numeric)
	
	trainGeno <-  apply(trainGeno_012,2,function(x) x-1)
	testGeno <- apply(testGeno_012,2,function(x) x-1)
	rownames(trainGeno) <- rownames(TrainData_Table_Num_Filt)

    ## TotGeno
	
	geno_012 <- rbind(trainGeno_012,testGeno_012)
    totGeno <-  apply(geno_012,2,function(x) x-1)
	rownames(totGeno) <- c(rownames(TrainData_Table_Num_Filt),rownames(TestData_Table_Num_Filt))
	
    Test <- testIds
    Candidates <- as.character(rownames(trainGeno))
 
    ## Train Set for STPGA and Random   
 
	# if(is.null(optTS)){
	  
	   # train_STPGA_Ind <- which(Candidates %in% as.character(Train_STPGA[[1]]$`Solution with rank 1`))
	# } 
	
	if(is.null(optTS)){
	  
	   train_STPGA_Ind <- c(1:length(Candidates))
	} 
    if(!is.null(optTS)){ 
  
       train_STPGA_Ind <- which(Candidates %in% as.character(as.vector(optTS))) 
	}
  
     train_Random_Ind <- which(Candidates %in% as.character(Train_Random[[1]]))
     trainGeno_STPGA <- trainGeno[train_STPGA_Ind,]
	 trainGeno_Random <- trainGeno[train_Random_Ind,]
  
	 
	## Impute missing values 
	
	 trainGeno_STPGA_Imp <- snpQC(trainGeno_STPGA,remove=FALSE,impute=TRUE)
	 trainGeno_Random_Imp <- snpQC(trainGeno_Random,remove=FALSE,impute=TRUE)

     trainPheno_STPGA <- trainPheno[train_STPGA_Ind]
     trainPheno_Random <- trainPheno[train_Random_Ind]
             
     pred_STPGA <- emCV(trainPheno_STPGA,trainGeno_STPGA_Imp,5,5)
	 pred_Random <- emCV(trainPheno_Random,trainGeno_Random_Imp,5,5)
		
	 
	 PA_Table <- rbind(pred_STPGA[c("emRR","emBB","emBL")],pred_Random[c("emRR","emBB","emBL")])
 	 colnames(PA_Table) <- c("emRR","emBB","emBL")
	 rownames(PA_Table) <-  c("STPGA Training Set","Random Training Set")
	 
	 return(PA_Table)
  }
  
   
   
 getTSComparisonsMT <- function(Data_Table_Num_Filt_List,Train_STPGA,Train_Random,trait,nTraits,testIds,optTS=NULL){ 
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
		
	
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(2:(ncol(TestData_Table_Num_Filt)))],2,as.numeric)
	
	trainGeno <-  apply(trainGeno_012,2,function(x) x-1)
	testGeno <- apply(testGeno_012,2,function(x) x-1)
	rownames(trainGeno) <- rownames(TrainData_Table_Num_Filt)

    ## TotGeno
	
	geno_012 <- rbind(trainGeno_012,testGeno_012)
    totGeno <-  apply(geno_012,2,function(x) x-1)
	rownames(totGeno) <- c(rownames(TrainData_Table_Num_Filt),rownames(TestData_Table_Num_Filt))
	
    Test <- testIds
    Candidates <- as.character(rownames(trainGeno))
 
    ## Train Set for STPGA and Random   
 
	# if(is.null(optTS)){
	  
	   # train_STPGA_Ind <- which(Candidates %in% as.character(Train_STPGA[[1]]$`Solution with rank 1`))
	# }
    if(is.null(optTS)){
	  
	   train_STPGA_Ind <- c(1:length(Candidates))
	} 	
    if(!is.null(optTS)){ 
  
       train_STPGA_Ind <- which(Candidates %in% as.character(as.vector(optTS))) 
	}
  
     train_Random_Ind <- which(Candidates %in% as.character(Train_Random[[1]]))
     trainGeno_STPGA <- trainGeno[train_STPGA_Ind,]
	 trainGeno_Random <- trainGeno[train_Random_Ind,]
  
	 
	## Impute missing values 
	
	 trainGeno_STPGA_Imp <- snpQC(trainGeno_STPGA,remove=FALSE,impute=TRUE)
	 trainGeno_Random_Imp <- snpQC(trainGeno_Random,remove=FALSE,impute=TRUE)
	
	### ST Models

  STModels<- FALSE
  if(STModels==TRUE){	
	 
	 PA_TableComb <- c()
	 for(nT in 1:length(trait)){
	 
	  trainPheno <- TrainData_Table_Num_Filt[,trait[nT]]
	  names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
   	
 	  trainPheno_STPGA <- trainPheno[train_STPGA_Ind]
      trainPheno_Random <- trainPheno[train_Random_Ind]
      NAIndices_STPGA <- which(is.na(trainPheno_STPGA))   
      NAIndices_Random<- which(is.na(trainPheno_Random))
	  if(length(NAIndices_STPGA) >= 1){
        pred_STPGA <- emCV(trainPheno_STPGA[-NAIndices_STPGA],trainGeno_STPGA_Imp[-NAIndices_STPGA,],5,5)
	  }
	  if(length(NAIndices_Random) >= 1){
	    pred_Random <- emCV(trainPheno_Random[-NAIndices_Random],trainGeno_Random_Imp[-NAIndices_Random,],5,5)
	  }
	  if(length(NAIndices_STPGA) <1){
        pred_STPGA <- emCV(trainPheno_STPGA,trainGeno_STPGA_Imp,5,5)
	  }
	  if(length(NAIndices_Random) < 1){
	   pred_Random <- emCV(trainPheno_Random,trainGeno_Random_Imp,5,5)
	  }
	  PA_Table <- rbind(rep(trait[nT],3),c("RR","BB","BL"),pred_STPGA[c("emRR","emBB","emBL")],pred_Random[c("emRR","emBB","emBL")])
 	  #colnames(PA_Table) <- c("emRR","emBB","emBL")
	  if(nT==1){
	   PA_Table2 <-  cbind(c("Trait","GPModel","STPGA","Random"),PA_Table)
	  }
	  if(nT>1){
	  
	    PA_Table2 <- PA_Table
	  }
	  
	  PA_TableComb <- cbind(PA_TableComb,PA_Table2)
	 } 
	  
	  PA_Out_Table <- PA_TableComb
    }  

  
	 
###MT Models 

  MTModels <- TRUE

  if(MTModels==TRUE){
	
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 trainPheno <- TrainData_Table_Num_Filt[,trait]
     rownames(trainPheno) <- rownames(TrainData_Table_Num_Filt)
   	
 	 	 
	 NAIndices <- unlist(apply(TrainData_Table_Num_Filt[,trait],2,function(x)which(is.na(x))))
   
    
	 trainPheno_STPGA <- trainPheno[train_STPGA_Ind,]
     trainPheno_Random <- trainPheno[train_Random_Ind,]
	 trainSet_Pheno <- list(trainPheno_STPGA,trainPheno_Random)
	 trainSet_Geno <-  list(trainGeno_STPGA_Imp,trainGeno_Random_Imp)
	 
	 PA_Out_List <- list()
	 
	 nIter <- 2
	 k<- 2
	 
     for(nTS in 1:2){ 	 
	
	  pheno_wNA<- trainSet_Pheno[[nTS]]
      geno_Imp_wNA <- trainSet_Geno[[nTS]]
	
	  
	  pheno <- pheno_wNA[-NAIndices,]
	  geno_Imp <- geno_Imp_wNA[-NAIndices,]
	    
	  Geno <- rownames(geno_Imp)
	  
	  rownames(geno_Imp) <- Geno
	  rownames(pheno) <- Geno
	  
	  
	  Y <- pheno
	  yNA <- as.matrix(pheno)
	  X <- geno_Imp
	  rownames(X) <- rownames(geno_Imp)
	  n <- nrow(X)
	  p <- ncol(X)
	  	 
		 
	  yNA_DF <- as.data.frame(yNA)
	  yNA_DF$id <- as.factor(rownames(X))
	  
	  PA_TableComb <- c()
	
##############################################################
	 
      GPModelMT_List <- list("BRR (BGLR)","RKHS (BGLR)","Spike-Slab (BGLR)") #,"GBLUP (SOMMER)")
      PA_Out_GP <- list()
  
      n <- nrow(trainPheno_STPGA)

      for(nGP in 1:length(GPModelMT_List)){ 
    
		GPModelMT <- GPModelMT_List[[nGP]]
		
		PA_List <- list()
		
		for(nrep in 2:nIter){
	 
		 nK <- floor(n/k)
		 k_List <- list()
		 set.seed(125+nrep) 
		 tot <- c(1:n)
		  for(i in 1:k){
			set.seed(125+nrep+k+5) 
		   k_List[[i]] <- sample(tot,nK)
		   tot <- setdiff(tot,k_List[[i]]) 
			
		  }

		  trIndices_List <- list()
		  tstIndices_List <- list()
		  for(i in 1:k){ 
		 
			  trIndices_List[[i]] <- unlist(k_List[-i])
			  tstIndices_List[[i]] <- k_List[[i]]
		  }	  
			 
		PA <- list()
		 
		for(i in 1:k){
		
		 trainIndices <-  trIndices_List[[i]]
		 testIndices <- tstIndices_List[[i]]
		 
		
	#########################################################
		
		 nTst <- length(testIndices)
		 
		 pheno[testIndices,] <- matrix(rep(NA,nTst*ncol(pheno)),nrow=nTst,ncol=ncol(pheno))
			 
	## Kinship Matrix 
	     
			
		 
		 if(GPModelMT == "BRR (BGLR)"){
					  
				ETA <- list(list(X=X,model="BRR"))
				fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
		 }
		
		 if(GPModelMT == "RKHS (BGLR)"){ 
				A.Tot <- A.mat(X)
				ETA <- list(list(K=A.Tot,model="RKHS"))
				fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
		 }
		 if(GPModelMT == "Spike-Slab (BGLR)"){ 
					  
				ETA <- list(list(X=X,model="SpikeSlab"))
				fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
		 }	
		 if(GPModelMT == "GBLUP (SOMMER)"){ 
				 
				A.Tot <- A.mat(X)
				rownames(A.Tot) <- rownames(X) 
				colnames(A.Tot) <- rownames(X)
				fm3 <- mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
				random=~vs(id,Gu=A.Tot),
				rcov=~units,
				data=yNA_DF,verbose = TRUE)
		 }
		 
		tst <- testIndices
	
		if(GPModelMT != "GBLUP (SOMMER)"){ 
		  Test_PredictedValues <- fm3$ETAHat[tst,]
		}
		if(GPModelMT == "GBLUP (SOMMER)"){  
		  UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
		  Test_PredictedValues <- c()
		  for(nT in 1:length(trait)){
		 
		    Test_PredictedValues <- cbind(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
		  }
		}
	 	 
	    PA[[i]] <- diag(cor(Test_PredictedValues,Y[tst,]))
	 
	}
	  
	   PA_List[[nrep]] <- do.call(rbind,lapply(PA,function(x) x))
	  
	}
	PA_Out_GP[[nGP]] <- apply(do.call(rbind,lapply(PA_List,function(x)x)),2,function(x) mean(x,na.rm=TRUE))
	
  }
  
	PA_Out1 <-  apply(do.call(rbind,lapply(PA_Out_GP,function(x) x)),2,function(x) round(x,digits=2))
	
	
	PA_Out2 <- cbind(unlist(GPModelMT_List),PA_Out1)
	PA_Out <- rbind(c("GPModel",colnames(PA_Out1)),PA_Out2)
    PA_Out_List[[nTS]] <- PA_Out

  }
  
  STPGAOut <- cbind(rep("STPGA",nrow(PA_Out_List[[1]])),PA_Out_List[[1]])
  RandomOut <- cbind(rep("Random",nrow(PA_Out_List[[2]])),PA_Out_List[[2]])
  
  PA_Out_Table1 <- sapply(c(1:nrow(STPGAOut)),function(x) rbind(STPGAOut[,x],RandomOut[,x]))
  PA_Out_Table <- PA_Out_Table1[-1,]
  PA_Out_Table[1,1] <- "TrainSet"
  } 
#########################	 
	 
	return(PA_Out_Table)
 }
  

# pred.tt <- kin.blup(as.data.frame(Data),Geno,Pheno,GAUSS=FALSE,K=A,covariate=NULL,PEV=TRUE,n.core=1,theta.seq=NULL)

  # VarG <- pred.tt$Vg
  # VarE <- pred.tt$Ve
  # A.tt <- A
  
# ## Reliability
 # I.tt <- diag(nrow(A.tt))
 # In.tt <- A.tt + (I.tt*(VarE/VarG))
 # In.tt.inv <- solve(In.tt) 
 
 # totGeno <- rbind(trainGeno_Imp,testGeno)
 
  
  # A.tot <- A.mat(totGeno)
  # trainIndices <- c(1:nrow(trainGeno_Imp))
  # testIndices <- c(nrow(trainGeno_Imp):nrow(totGeno))
  # A.ut <- A.tot[trainIndices,testIndices]
  # A.ut.t <- A.tot[testIndices,trainIndices]
  # R <- diag(A.ut.t%*%In.tt.inv%*%A.ut)

# ### Upper bound of Reliability

    # trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) x+1)  
    # Rep1 <- cleanREP(trainPheno,trainGeno_Imp2)
   
    # M <-  Rep1[[2]]
    # M.Pdt <- t(M)%*% solve(M %*% t(M)) %*% M
    # M.Pdt <- t(MQC)%*% solve(MQC %*% t(MQC)) %*% MQC
    # v.hat <- M.Pdt %*% v
    # U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
   
   # getU <- function(M.Pdt,v){
	   # v.hat <- M.Pdt %*% v
	   # U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   # return(U)
    # }
    # testGeno <- TestGenoTable
    # U <- apply(testGeno,1,function(x) getU(M.Pdt,x)) 
	
	# plot(Test_PredictedValues,U) 
	
	
	
	 # totGeno <- rbind(trainGeno_Imp,testGeno)
	 # rownames(totGeno) <- c(rownames(TrainData_Table_Num_Filt),rownames(Test_Genotypes_Table_Mod_Num_Filt))
	 # Candidates <- rownames(totGeno) 
	 # Test <- rownames(Test_Genotypes_Table_Mod_Num_Filt) 
	 
	 
	 


	
# function to remove repeated genotypes
cleanREPV2 = function(y,gen,fam=NULL,thr=0.95){ 

  if(is.vector(y)) y=matrix(y,ncol=1)
  if(is.null(fam)) fam = rep(1,nrow(y))
    
  ibs = function(gen){
  f1 = function(x,gen) apply(gen,1,function(y,x) mean(y==x,na.rm = T),x=x)
  f2 = apply(gen,1,f1,gen=gen)
  return(f2)}  
  
  GG = function(gen, r = 1) {
    a1 = (gen - 1)
    a1[a1 == -1] = 0
    A1 = (tcrossprod(a1))
    a2 = -(gen - 1)
    a2[a2 == -1] = 0
    A2 = (tcrossprod(a2))
    d = round(exp(-abs(gen - 1)))
    D = tcrossprod(d)
    G = A1 + A2 + D
    G = (G/ncol(gen))^r
    return(G)
  }
  cat("solving identity matrix\n")
  G = GG(gen)
  
  rownames(G) = 1:nrow(G)
  
  lt = G*lower.tri(G) # lower triang
  r = 1* lt>thr # logical matrix: repeatitions
  # starting point of new data
  #rownames(gen) = 1:nrow(gen)
  Ny=y;  Nfam=fam;  Ngen=gen 
  NGen <- gen
  # summary
  cs = colSums(r) # how many times id will be repeated
  while(any(cs>0)){
    i = which(cs>0)[1]
    cat("indiviual",rownames(gen)[i],"had",cs[i],'duplicate(s)\n')
    w = which(r[,i])
    if(ncol(y)>1){y[i,] = colMeans(y[c(i,w),],na.rm=T)
    }else{y[i] = mean(y[c(i,w)],na.rm=T)}
    if(ncol(y)>1){Ny=Ny[-w,]}else{Ny=Ny[-w]}
    Nfam=Nfam[-w]
    Ngen=Ngen[-w,]
	rownames(Ngen) <- rownames(NGen)[-w]
	NGen <- Ngen
    r = r[-w,]
    cs = colSums(r)
  }
  return(list(y=Ny,gen=Ngen,fam=Nfam))
}

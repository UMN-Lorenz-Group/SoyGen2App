	load("TS_Demo_PreComputed_In.RData") 

	source("/home/lorenza/vramasub/GP_Soy/GS_Pipeline_Jan_2022_FnsApp.R") 

	noToReduce <- nrow(NUST_Data_Table_Num_Filt[[1]])
	nTrainToSelect <- 1000
	optTS <- NULL
	 
	if(is.null(optTS)){
	  TS_STPGA <- getOptimalTS(NUST_Data_Table_Num_Filt,trait,nTraits,noToReduce,nTrainToSelect)
	}
	
### 
	if(!is.null(optTS)){
	  TS_STPGA <- as.character(as.vector(optTS))
	}
 
	TS_Random <- getRandomTS(NUST_Data_Table_Num_Filt,trait,nTraits,testIDs,noToReduce,nTrainToSelect) 
	PA_Table <- getTSComparisons(NUST_Data_Table_Num_Filt,TS_STPGA,TS_Random,trait,nTraits,testIDs,optTS)
	 
    save(TS_STPGA,TS_Random,PA_Table,file="TS_Demo_PreComputed_Out.RData") 
	
####################################################################################333
	
	
	load("TS_Demo_PreComputed_In.RData") 

	source("/home/lorenza/vramasub/GP_Soy/GS_Pipeline_Jan_2022_FnsApp.R") 
	
	noToReduce <- nrow(NUST_Data_Table_Num_Filt[[1]])
	nTrainToSelect <- 1500
	optTS <- NULL

  
    Data_Table_Num_Filt_List <- NUST_Data_Table_Num_Filt
  
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	
	trainPheno <- TrainData_Table_Num_Filt[,trait]
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(2:(ncol(TestData_Table_Num_Filt)))],2,as.numeric) 
	
## Complete Genotypes table
	
	geno_012 <- rbind(trainGeno_012,testGeno_012)
    totGeno <-  apply(geno_012,2,function(x) x-1)
	rownames(totGeno) <- c(rownames(TrainData_Table_Num_Filt),rownames(TestData_Table_Num_Filt))
	names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
   
    totCandidates <-  rownames(totGeno)
	
    Test <- rownames(TestData_Table_Num_Filt)
	testIndices <- which(totCandidates %in% Test)
	
	testGeno_012 <- TestData_Table_Num_Filt
 
    geno_012 <- rbind(trainGeno_012,testGeno_012)
### Set1 with 99 PCs, 100 candidates , 100 test and 50 nTrain 

    set.seed(125)
	
## Reduce or keep candidate genotypic data 
    
    nTrain <- nrow(TrainData_Table_Num_Filt)
	
	if(noToReduce!= nTrain){
	  trainSetIndices <-sample(c(1:nTrain),noToReduce)
	}
	
	if(noToReduce== nTrain){
	   trainSetIndices <- c(1:nTrain)
	}
	
	trainGeno_Pre <- totGeno[trainSetIndices,]
	trainPheno_Pre <- trainPheno[trainSetIndices] 
	
	trainGenoPre_Imp <- snpQC(trainGeno_Pre,remove=FALSE,impute=TRUE)
	
	trainGeno_Pre_Imp <- apply(trainGenoPre_Imp,2,function(x) x+1)
	trainClean <- cleanREPV2(trainPheno_Pre,trainGeno_Pre_Imp)
	
		
	trainGeno_Clean <- trainClean[[2]]
	
    #G <-  rbind(trainGeno_Clean,totGeno[testIndices,])
	rownames(G) <- c(rownames(trainGeno_Clean),rownames(totGeno)[testIndices])
	
	
	G_Imp <- snpQC(G,remove=FALSE,impute=TRUE)
	rownames(G_Imp) <- rownames(G) 
	
	
###  

 	
	 trainGeno_Clean <- trainClean[[2]]
	 trainPhenoClean <- trainClean[[1]]
	 Candidates <- rownames(trainGeno_Clean)
     Geno <- Candidates
     A <- A.mat(trainGeno_Clean)	
	 colnames(A) <- Geno
	 rownames(A) <- Geno
	

## Prepare Data Table for GP 
	 Data <- cbind.data.frame(Geno,trainPhenoClean)
	 
	 colnames(Data) <- c("Geno","Pheno")
	 Geno <- "Geno"
	 Pheno <- "Pheno"
  
  
### Reliability
  
    pred.tt <- kin.blup(as.data.frame(Data),Geno,Pheno,GAUSS=FALSE,K=A,covariate=NULL,PEV=TRUE,n.core=1,theta.seq=NULL)

    VarG <- pred.tt$Vg
    VarE <- pred.tt$Ve
    
	R2 <-  1- (pred.tt$PEV /(VarG *diag(A.tt)))  
 
 
    A.tt <- A.mat(trainGeno_Clean)

    I.tt <- diag(nrow(A.tt))
	
    In.tt <- A.tt + (I.tt*(VarE/VarG))
    In.tt.inv <- solve(In.tt) 

    G <-  rbind(trainGeno_Clean,apply(testGeno_012[,-1],2,as.numeric))
   
    totGeno <- G
    A.tot <- A.mat(totGeno)
    trainIndices <- c(1:nrow(trainGeno_Clean))
    testIndices <- c((nrow(trainGeno_Clean)+1):nrow(totGeno))
    A.ut <- A.tot[testIndices,trainIndices]
    A.ut.t <- t(A.ut)
    R <- diag(A.ut%*%In.tt.inv%*%A.ut.t)

    A.uu <- A.tot[testIndices,testIndices]

  

   Geno.tot <- rownames(totGeno)
   Pheno.tot <- c(trainPhenoClean,rep(NA,nrow(testGeno_012)))
   
   Data.tot <- cbind.data.frame(Geno.tot,Pheno.tot)
   colnames(Data.tot) <- c("Geno","Pheno")
   Geno <- "Geno"
   Pheno <- "Pheno"

   pred.tot <- kin.blup(as.data.frame(Data.tot),Geno,Pheno,GAUSS=FALSE,K=A.tot,covariate=NULL,PEV=TRUE,n.core=1,theta.seq=NULL)
   R2.tot <-  1- (pred.tot$PEV /(pred.tot$Vg *diag(A.tot)))  
   summary(R2.tot[testIndices])
   plot(R2.tot[testIndices],pred.tot$pred) 
  
   save(Data.tot,A.tot,R2.tot,pred.tot,file="pred_kin_blup_tot.RData")
   
   
      pred.tot2 <- kin.blup(as.data.frame(Data.tot),Geno,Pheno,GAUSS=FALSE,K=A.tot,covariate=NULL,PEV=TRUE,n.core=1,theta.seq=NULL)


### Upper Bound of Reliability

	testGeno <- apply(testGeno_012[,-1],2,as.numeric)
	
# i) 
	M <-  trainGeno_Clean[1:200,]
   # summary(U)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   # 0.5357  0.5960  0.6156  0.6191  0.6409  0.7046
	
# ii) 
### a) STPGA indices
	
  TS_Rank1 <- TS_STPGA$`Solution with rank 1`
  length(TS_Rank1)
#[1] 1500
  
  STPGA_Rank1Indices <- which(rownames(trainGeno_Clean) %in% TS_Rank1)
  length(STPGA_Rank1Indices)
# [1] 1500

  M <- trainGeno_Clean[STPGA_Rank1Indices,]
	
	
	 (solve(M.train %*% t(M.train) + (diag(nrow(M.train)) *(VarE/VarG))))
	
 # summary(U)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.9484  0.9609  0.9639  0.9642  0.9675  0.9788
 
### b) Random indices sample of 1500 lines 
 
  M <- trainGeno_Clean[sample(c(1:nrow(trainGeno_Clean)),1500),] 
 
 
 # summary(U)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.9493  0.9601  0.9634  0.9637  0.9671  0.9794

 
# iii) 
 
    M <-  trainGeno_Clean[1:500,]
	
 # summary(U)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.7014  0.7375  0.7518  0.7554  0.7730  0.8216
 
 
# iv) 
 
      M <-  trainGeno_Clean
 
 
### Estimate Ui	
	
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M
      
    getU <- function(M.Pdt,v){
	   v.hat <- M.Pdt %*% v
	   U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   return(U)
    }
   
   U <- apply(testGeno,1,function(x) getU(M.Pdt,x)) 
   
   
     M.train <- trainGeno_Clean
     (solve(M.train %*% t(M.train) + (diag(nrow(M.train)) *(VarE/VarG))))
   
	
 # summary(U)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.9752  0.9809  0.9824  0.9826  0.9842  0.9901

	
#############################################################################################
## Parameters: 

 #  nTrainToSelect - 500, 1500 
 #  npop/nelite - 100/5,10 ; 200/10,20
 #  niterations - 100, 200 
   
	
	
	GenoSVD <- svd(G_Imp,nu=99,nv=100)
    PC100 <- G%*%GenoSVD$v
    rownames(PC100)<-rownames(G_Imp)
  	
	STPGA_Select <- function(PC,Candidates,Test,nTrainToSelect){ 
	  
	     
	    Train_STPGA <- GenAlgForSubsetSelection(P=PC,Candidates,Test,ntoselect=nTrainToSelect,InitPop=NULL,
         npop=100, nelite=10, mutprob=.5, mutintensity = 1,
         niterations=100,minitbefstop=50, tabu=TRUE,
         tabumemsize = 1,plotiters=FALSE,
         lambda=1e-6, mc.cores=10) 
		 
		return(Train_STPGA)
	} 
	 
	 
	 library(foreach)
	 library(doParallel)
	 registerDoParallel(10)  
	 
	 
	 nTrainToSelect2 <- 500
	 Train_STPGA_List_500 <- foreach(k=1:10) %dopar%  (STPGA_Select(PC100,Candidates,Test,nTrainToSelect2)) 
	 
	 
	 
	  
	PA_Table_List_500 <- list()
	
	for(k in 1:10){
	
	 TS_Random <- getRandomTS(NUST_Data_Table_Num_Filt,trait,nTraits,testIDs,noToReduce,nTrainToSelect2) 
	
	 PA_Table <- getTSComparisons(NUST_Data_Table_Num_Filt,Train_STPGA_List_500[[k]],TS_Random,trait,nTraits,testIDs,optTS)
	
	 PA_Table_List_500[[k]] <- PA_Table
	}
	
	
	PA_Diff_Tab <- do.call(rbind,lapply(PA_Table_List,function(x) x[1,]-x[2,]))

	 apply(PA_Diff_Tab,2,mean)
	
	
		
	PA_Diff_Tab_500 <- do.call(rbind,lapply(PA_Table_List_500,function(x) x[1,]-x[2,]))

	apply(PA_Diff_Tab_500,2,mean)
	
	
   # Error in if (meansdiff < tolconv) { :
   # missing value where TRUE/FALSE needed
   # Timing stopped at: 4901 81.5 4986

   ## Error fixed after cleanREP (removing duplicates in training geno data
   
######


   system.time({
         Train_STPGA <- GenAlgForSubsetSelection(P=PC100,Candidates,Test,ntoselect=nTrainToSelect,InitPop=NULL,
         npop=100, nelite=10, mutprob=.5, mutintensity = 1,
         niterations=100,minitbefstop=50, tabu=TRUE,
         tabumemsize = 1,plotiters=FALSE,
         lambda=1e-6, mc.cores=10)
    })
    
#############

	GenoSVD <- svd(G_Imp,nu=99,nv=500)
    PC500 <- G%*%GenoSVD$v
    rownames(PC500)<-rownames(G_Imp)

	 system.time({
         Train_STPGA_500 <- GenAlgForSubsetSelection(P=PC500,Candidates,Test,ntoselect=nTrainToSelect,InitPop=NULL,
         npop=200, nelite=10, mutprob=.5, mutintensity = 1,
         niterations=100,minitbefstop=50, tabu=TRUE,
         tabumemsize = 1,plotiters=FALSE,
         lambda=1e-6, mc.cores=10)
    })
	
	
	
	
	
	
	
	
	
##	
	 
#	 500 -200
	 
  # TS	      RR	  BB	  BL
# STPGA	    0.5371	0.4739	0.5884
# Random	0.3777	0.452	0.4731 



####


    nTrainToSelect <- 500
 	GenoSVD <- svd(G_Imp,nu=50,nv=50)
    PC50 <- G%*%GenoSVD$v
    rownames(PC50)<-rownames(G_Imp)

	 system.time({
         Train_STPGA_50 <- GenAlgForSubsetSelection(P=PC50,Candidates,Test,ntoselect=nTrainToSelect,InitPop=NULL,
         npop=100, nelite=10, mutprob=.5, mutintensity = 1,
         niterations=100,minitbefstop=50, tabu=TRUE,
         tabumemsize = 1,plotiters=FALSE,
         lambda=1e-6, mc.cores=10,tolconv=1e-9)
    }) 
	
	
#########################

    nTrainToSelect <- 500
 	GenoSVD <- svd(G_Imp,nu=50,nv=50)
    PC50 <- G%*%GenoSVD$v
    rownames(PC50)<-rownames(G_Imp)

	 system.time({
         Train_STPGA_50 <- GenAlgForSubsetSelection(P=PC50,Candidates,Test,ntoselect=nTrainToSelect,InitPop=NULL,
         npop=300, nelite=30, mutprob=.5, mutintensity = 1,
         niterations=100,minitbefstop=50, tabu=TRUE,
         tabumemsize = 1,plotiters=FALSE,
         lambda=1e-6, mc.cores=10,tolconv=1e-9)
    })
	
	
	
	
###cleanREP 
 
  # y <- trainPheno_Pre
  # gen <- apply(trainGenoPre_Imp,2,function(x)x+1)
  # fam=NULL
  # thr=0.95
	

# trainClean <- list(y=Ny,gen=Ngen,fam=Nfam)
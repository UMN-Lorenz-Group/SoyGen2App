if(!require(shiny)){
	install.packages("shiny")
}
library(shiny) 

### Change source file path to working directory

FN <- paste(getwd(),"/GS_Pipeline_Jan_2022_FnsApp.R",sep="")
source(FN)

PN <-  paste(getwd(),"/GSPipeline.png",sep="")

## Set the file upload limit to X MB
options(shiny.maxRequestSize=500*1024^2)

ui <- fluidPage( 
  #theme = bslib::bs_theme(bootswatch = "solar"),
  tags$style(HTML(".message-text { color: red; }")),
  
  fluidRow(
    column(4),
    column(6, div(
      id = "app-title",
      titlePanel(tags$strong("SOYGEN2 Genomic Selection App")),
      tags$p("An app to perform genomic predictions given genotypic and phenotypic data")
    ))),
  
  tabsetPanel(id="inData",
              
              tabPanel("Home",
                      fluidRow(
                         column(5),column(width=4,tags$h3(tags$strong("SOYGEN2 App")))
                              # column(6,offset=6),actionButton(inputId ="Home_Data", "next")
                      ),
                      fluidRow(
                        column(2),column(width=9,tags$h5("The SOYGEN2 App implements a genomic selection pipeline designed by public soybean breeders in the US MidWest. 
                        The GS pipeline involves steps starting from data management to making selection decisions based on genome estimated breeding values.The steps 
                        in the pipeline are depicted in the schematic below. The current version starts with the upload of filtered and imputed data exported from databases.") 
                      )),
                      tags$br(),
                      tags$br(),
                      fluidRow(
                        column(2),column(width=9,tags$a(tags$img(src="GSPipeline.png",title="Genomic Selection Pipeline",width="750",height="400")))), 
                           
                      fluidRow(
                        column(2),column(width=9,tags$h6(("Contributors: Vishnu Ramasubramanian, Cleiton Wartha, Paolo Vitale, Sushan Ru,and Aaron Lorenz")))),
                      tags$br(),        
                      fluidRow(
                        column(2),column(width=9,tags$h6("Contact: Vishnu Ramasubramanian - vramasub@umn.edu")))
               ),
            ## 
     
  
              
     ## Tab for loading data
     
       tabPanel("Load Data",
                       
                       fluidRow(
                          column(4),column(width=4,tags$h3(tags$strong("Load Data Files"))),   
                          #column(10),column(width=2,actionButton("Data_Home", "prev"),
                          #actionButton("Data_Trait", "next"))
                          ), 
                       tags$br(),
                       fluidRow(
                         column(1),column(width=3,tags$h4(tags$strong("Phenotypic Data"))),
                         column(5,offset=2,tags$h4(tags$strong("Genotypic Data"))) 
                       ),
                       
                      ####  
                     
                     tags$br(),
                      fluidRow(
                         column(1),column(width=5,fileInput("infileBLUEs", "Choose Phenotype File (CSV)", accept = ".csv")),
                         column(6),column(width=5,fileInput("infileVCF", "Choose Genotype File (VCF)", accept = ".vcf"))
                      ),
                      fluidRow(
                        column(1),column(width=5,checkboxInput("header", "Header", TRUE)),
                        column(6),column(width=5,checkboxInput("header", "Header", TRUE))
                      ),
                     fluidRow(
                       column(1),column(width=5,tags$ul(tags$li("A '.csv' file containing phenotypic information. Phenotype file should have the line ID in the first column and phenotypic trait values after that with the trait ID as column name."))),
                       column(width=5,tags$ul(tags$li("A '.vcf' containing the genotypic information of both the candidate lines used to train the GP model and the target lines whose values need to be predicted.")))
                     ),
                     tags$br(),
                     fluidRow(
                           column(1),column(width=5,tags$h4(tags$strong(textOutput("PhenoHeader")))),
                           column(6),column(width=5,tags$h4(tags$strong(textOutput("GenoHeader"))))
                      ),
                      fluidRow(
                          column(1),column(width=5,tags$h6(tableOutput("PhenoTable"))),
                          column(6),column(width=5,tags$h6(tableOutput("GenoTable")))),
                      tags$br()
        ),  
       
      
     ###  
     tabPanel("Filter Genotypic Data",
       
    # mainPanel(
       tags$br(),
       tags$br(),
       fluidRow(column(1),
         column(width=10,"Filter sites (markers) and taxa (lines) in genotype table using rTASSEL (Monier et al. 2021), if you uploaded a raw genotype file or a QC genotype file that requires filtering. For example, it is common to 
                    set minimum number of sites that are not missing to 80% of the lines in the input genotype data and set the MAF threshold to 0.02/0.05. 
                    If you uploaded a QC genotype file, you can choose to skip this step")),
       
       tags$br(),
       tags$br(),
       
       fluidRow(
         column(1),column(width=3,tags$strong(tags$h4("Filter Genotype Table Sites"))),
         #Filt2
         column(5,offset=2,tags$strong(tags$h4("Filter Genotype Table Taxa")))),
       tags$br(),
       fluidRow(
         #Filt1
         column(1),column(width=3,
                          numericInput(inputId="siteMinCnt","Minimum Site Count (TASSEL)",value = 0,min =0, max=0)),
         #Filt2
         column(5,offset=2,numericInput(inputId="minNotMissing","Minimum Not Missing (TASSEL)",value = 0.9,min =0.5, max=1))),
         tags$br(),
    
       tags$br(),
       fluidRow(
         #filt1
         column(1),column(width=3,
                          numericInput(inputId="MAF","Minimum Allele Frequency",value =0.02,min =0, max=0.5)),
         #Filt2
         column(5,offset=2,
           actionButton(inputId="FilterTaxa","Filter Genotype Table Taxa"))),
    
       tags$br(),
       fluidRow(
         #Filt1
         column(1),column(width=3,
                          actionButton(inputId="FilterSites","Filter Genotype Table Sites"))),
         
       tags$br(),
       tags$br(),
       fluidRow(
         column(1),column(width=3,
                          checkboxInput("setGenoFilt1Tas", "Use Filtered Genotypes From Filter Sites For Next Steps", TRUE)),
         column(5,offset=2,
           checkboxInput("setGenoFilt2Tas", "Use Filtered Genotypes From Filter Taxa For Next Steps", TRUE))),  
       tags$br(),
       tags$br(),
       # fluidRow(
       #   column(6),column(width=3,tags$strong(tags$h4("Filter Genotype Table Taxa")))),
       # tags$br(),
       #fluidRow(
         # column(5,offset=2,numericInput(inputId="minNotMissing","Minimum Not Missing (TASSEL)",value = 0.9,min =0.5, max=1))),
       tags$br(),
     # fluidRow(
     #  column(5,offset=0,
     #   actionButton(inputId="FilterTaxa","Filter Genotype Table Taxa"))),
       tags$br(),
       tags$br(),
     # fluidRow(
     #   column(5,offset=0,
     #    checkboxInput("setGenoFilt2Tas", "Use Filtered Genotypes From Filter Taxa For Next Steps", TRUE))),
     # 
       fluidRow(column(1),column(width=4,tags$h4(tags$strong(textOutput("GenoFiltHeader")))),
                column(7),column(width=4,tags$h4(tags$strong(textOutput("GenoFiltHeader2"))))),      
       tags$br(),
       fluidRow(column(1),column(width=4,tags$h6(tableOutput("FilteredGenoTable"))),
                column(7),column(width=4,tags$h6(tableOutput("FilteredGenoTable2")))),
       #fluidRow(column(6),column(width=5,tags$h4(tags$strong(textOutput("GenoFiltHeader2"))))),
       tags$br(),
       #fluidRow(column(6),column(width=5,tags$h6(tableOutput("FilteredGenoTable2")))),
       tags$br(),
       tags$br()
     ),
     
    tabPanel("Impute Genotypic Data",  
              sidebarLayout(
               sidebarPanel(
                tags$strong(tags$h4("Impute Genotype Table ")),
                tags$br(),
                selectInput(inputId="imputeMet","Select Imputation Method",choices=c("LDKNNI","Numeric","AlphaPlantImpute"),multiple=FALSE,selected="LDKNNI"),
                tags$br(),
                
                # 
                conditionalPanel(condition="input.imputeMet == 'LDKNNI'",
                           numericInput(inputId="l","Number of High LD Sites",value = 30,min =30,max=1000),
                           numericInput(inputId="k","Neighboring Samples",value =10,min =10, max=100),
                ),
                conditionalPanel(condition="input.imputeMet == 'Numeric'",
                                 numericInput(inputId="nN","Number of Nearest Neighbors",value = 5,min =2,max=100),
                                 selectInput(inputId="Dist","Distance",choices=c("Euclidean", "Manhattan", "Cosine"),multiple=FALSE,selected="Euclidean"),
                ),
                
                conditionalPanel(condition="input.imputeMet == 'Numeric' || input.imputeMet == 'LDKNNI'",
                actionButton(inputId="Impute","Impute Genotype Scores"),
                tags$br(),
                tags$br(),
                tags$br(),
                
                ),
                tags$br(),
                tags$br(),
                downloadButton("ExportImpGenoDF", "Export Imputed Genotypes"),
                tags$br(),
                tags$br(),
                tags$br()
              ),
              mainPanel(
                
                conditionalPanel(condition="input.imputeMet == 'LDKNNI' || input.imputeMet == 'Numeric'",
                tags$br(),
                tags$br(),
                  fluidRow(
                     column(width=10,tags$h5(
                    "Impute missing genotype scores using rTASSEL (Monier et al. 2021), if you uploaded a raw genotype file or a QC genotype file with missing scores. 
                     Available options include numeric imputation and imputation using LD-K-Nearest Neighbors method are availble. 
                     For LD-KNN imputation set parameters l and K (Money et al. 2015). l corresponds to the number of high LD sites  
                     and k corresponds to the number of neighboring samples that 
                     are used to impute scores."))
                  ),
                tags$br(),
                  
                fluidRow(column(2),column(width=8,tags$h4(tags$strong(textOutput("GenoImpHeader"))))),
                tags$br(),
                fluidRow(column(2),column(width=8,tags$h6(tags$strong(tableOutput("ImputedGenoTable"))))),
                tags$br(),
                tags$br(),
                tags$br(),
                tags$br(),
                checkboxInput("setGenoImpTas", "Use Imputed Genotypes For Next Steps", TRUE), 
                tags$br(),
                tags$br(),
                ),
                conditionalPanel(condition="input.imputeMet == 'AlphaPlantImpute'",
                uiOutput("APIUI"),
                )
              )
            )),
       
## Tab for Trait Selection 
## 
              
            tabPanel("Set Target & Trait",
                   sidebarLayout(
                       sidebarPanel(
                         tags$br(),
                         tags$strong(tags$h4("Export merged training table & target genotypes table after selecting trait")),
                         tags$br(),
                         tags$strong(tags$h5("Download Merged Geno-Pheno Training Table")),
                         tags$br(),
                         downloadButton("ExportMerged_CompTS_DF", "Export Merged Training Data"),
                         tags$br(),
                         tags$br(),
                         tags$strong(tags$h5("Download Filtered Target Genotype Table")),
                         tags$br(),
                         downloadButton("ExportTargetGeno_DF", "Export Target Geno Data"),
                         tags$br(),
                         tags$br()
                       ),
                    
                      mainPanel(
                       tags$br(),
                       fluidRow(
                         column(2),column(width=4,tags$h4(tags$strong("Set Target Population"))),
                         column(6),column(width=4,tags$h4(tags$strong("Select Trait"))),
                       ),
                       tags$br(),
                       fluidRow(
                       column(1),column(width=5,fileInput("infileTargetTable", "Choose Target Data Information File (CSV)", accept = ".csv")),
                       column(5),column(width=5,selectInput(inputId="trait","Choose One or More Traits",choices=NULL,multiple=TRUE))),
                       fluidRow( column(1),column(width=5,checkboxInput("header", "Header", TRUE))),
                       tags$br(),
                       fluidRow(
                         column(1),column(width=5,tags$ul(tags$li(" A '.csv' file containing information on the target data, whose genetic values are to be predicted using the genomic prediction model. 
                       The first column in the target data information file should have the line IDs."))), 
                         column(6),column(width=5,tags$ul(tags$li("Select one trait or multiple traits to train the genomic prediction model.")))),
                       tags$br(), 
                       fluidRow(
                        column(1),column(width=5,div(
                                    tags$h4(tags$strong(textOutput("TargetHeader"))),
                                    tags$h6(tableOutput("TargetTable")))),
                         column(5),column(width=5,div( 
                           tags$h4(tags$strong(textOutput("summaryHeader"))),
                           tags$h6(tableOutput("Summary"))))
                        ),
                       
                       # Define a CSS class for text color
                     
                       # fluidRow(
                       #   column(4),column(width=5,div(
                       #     tags$h4(tags$strong(uiOutput("MssgTarget")))))),
                           #tags$h6(tableOutput("TargetTable"))))
                       tags$br()
                    )
                   )
                 ),
              
### TS optimization tab
              tabPanel("Optimize Training Population",
                       sidebarLayout(
                         sidebarPanel(
                           numericInput(inputId="noCandidates","Candidate Set Size",value = 0,min =0, max=0),
                           tags$br(),
                           numericInput(inputId="noToSelect","Training Set Size",value =0,min =0, max=0),
                           tags$br(),
                           selectInput(inputId="optCriteria","Select Optimization Criteria",choices=c("PEVMEAN2","PEVMAX2","CDMEAN2","CDMAX2","DOpt","AOpt","EOpt"),multiple=TRUE,selected="PEVMEAN2"),
                           tags$br(),
                           actionButton(inputId="Optimize","Optimize Training Sets"),
                           tags$br(),
                           tags$br(),
                           tags$br(),
                           downloadButton("ExportOptTS", "Export Optimal TrainSet"),
                           tags$br(),
                           tags$br(),
                           tags$br(),
                           checkboxInput("setGA", "Use Default Genetic Algorithm Parameters", TRUE),
                           conditionalPanel(condition="input.setGA == false",
                                            numericInput(inputId="noPop","GA Population Size",value = 100,min =1,max=1000),
                                            
                                            numericInput(inputId="noElite","Elite Population Size",value =10,min =1, max=50),
                                            
                                            numericInput(inputId="mutProb","Mutation Probability",value =0.5,min =0.1, max=1),
                                            
                                            numericInput(inputId="mutInt","Mutation Intensity",value =1,min =0, max=1),
                                            
                                            numericInput(inputId="nIterGA","# Iterations",value =100,min =50, max=1000),
                                            
                                            selectInput(inputId="tabu","Tabu Search",choices=c("TRUE","FALSE"),selected="TRUE"),
                                            
                                            numericInput(inputId="tabumemsize","Tabu Memory Size",value =1,min =1, max=10),
                                            
                                            numericInput(inputId="lambda","Shrinkage Parameter (lambda)",value =1e-6,min =1e-9, max=1),
                                            
                                            numericInput(inputId="mcCores","#Cores",value =10,min =1, max=20)
                                            
                          )                          
                          
                         ), mainPanel(
                           fluidRow(
                             column(width=10,tags$h3(tags$strong("Select Optimal Training Population using 'STPGA' "))),
                             #actionButton("TS_Trait", "prev"),
                             #actionButton("TS_CVR", "next")
                           ),
                           fluidRow(
                             
                             column(width=10, tags$p("A training population for genomic prediction is selected by optimizing an objective function using a genetic algorithm.
                              The optimization method uses the 'GenAlgForSubsetSelection' function implemented in the 'STPGA' R package (Akdemir 2017).
                              In this implementation, the parameters for the genetic algorithm are set to default values.") 
                            )),
                           tags$br(),
                           tags$br(),
                           fluidRow(
                             
                             column(width=10, tags$p("The user needs to set the candidate and training population size. Candidate population refers to the set of genotypes with phenotypic data. The default size is set to the number of genotypes in the input genotypic data that are not in the target population set. 
                             Training population size refers to the number of genotypes that are to be selected in the training population. The default value is set to 80% of the candidate set.") 
                             )),
                           tags$br(),
                           tags$br(),
                           
                           fluidRow(
                             column(2),
                             column(6, div(
                               tags$h5(tags$strong(textOutput("tsOptHeader"))),
                               tags$h6(tableOutput("PredAccuracyTable")), 
                              
                             )))
                         )
              )), 
              
              
              ### CV tab
              tabPanel("Cross Validations",
                sidebarLayout(
                    sidebarPanel(
                           numericInput(inputId="k",label = "Enter k for k-fold cross validation",value=2,min=2,max=10),
                           numericInput(inputId="nIter",label = "Enter n for n-iterations of each cycle of cross validation",value=2,min=2,max=10),
                           tags$br(),
                           tags$br(),
                           actionButton("CrossValidationST", "Run Cross Validation for Single Trait Models"),
                           tags$br(),
                           tags$br(),
                           tags$br(),
                           tags$br(),
                           actionButton("CrossValidationMT", "Run Cross Validation for Multi-trait Models"),
                           tags$br()
                         ), 
                        mainPanel(
                         fluidRow(
                           column(1),column(width=7,tags$h3(tags$strong("Cross Validation of GP Models"))), 
                           #actionButton("CVR_TS", "prev"),
                           #actionButton("CVR_GP", "next")
                         
                         ),
                         fluidRow(
                           
                           column(width=10, tags$p("Perform an optional k-fold cross validation with training dataset to identify the model with the best prediction accuracy.
                                           K-fold cross validation is performed using 'emCV' function implemented in the 'bWGR' package (Xavier et al. 2020).
                                           For multivariate models, 'Multitrait' and 'mmer' functions implemented in BGLR and Sommer packages are used for cross validations."))),
                         tags$br(),
                         fluidRow(
                           column(1),column(width=8,tags$h5(tags$strong(textOutput("cvrHeader"))))),
                         tags$br(),
                         
                         fluidRow(
                            column(2),column(width=8,tags$h5(tags$strong(textOutput("Mssg"))))),
                         tags$br(),
                         fluidRow(
                           column(2),column(width=8,tags$h5(tableOutput("emCVRST")))),
                         tags$br(),
                        # fluidRow(
                         #  column(3),column(width=8,tags$h5(tableOutput("emCVRMT")))),
                         tags$br()
                         
                  ) 
              )),
              #   
      ### GP Tab 
             tabPanel("Genomic Prediction",
                      tabsetPanel(id="GPModels",
                                   
                       tabPanel("Single Trait",
                       
                       sidebarLayout(
                         sidebarPanel(
                           selectInput(inputId="TrainSet","Choose Training Set",c("Complete Input Genotype Set","Optimal Train Set from Step 3","Random Train Set from Step 3"),selected="Complete Input Genotype Set"),
                           tags$br(),
                           tags$br(),
                           selectInput(inputId="fixed","Choose Fixed Effect",choices=NULL),
                           tags$br(),
                           tags$br(),
                           selectInput(inputId="GPModelST","Choose Prediction Model for Single Trait",c("rrBLUP (rrBLUP)","rrBLUP (bWGR)","BayesB (bWGR)","BayesLASSO (bWGR)")),
                           tags$br(),
                           tags$br(),
                           actionButton("RunPredictionsST", "Predict Single Trait!"),
                           tags$br(),
                           tags$br(),
                           tags$br(),
                           tags$br(),
                           downloadButton("ExportOut", "Export Output Table")
						   
                         ), mainPanel(
                           fluidRow(
                             column(width=8,tags$h3(tags$strong("Train Genomic Prediction Model"))), 
                            #actionButton("GP_CVR", "prev"),
                            # actionButton("GP_VP", "next")
                           ),
                           fluidRow(
                            
                             column(width=10, tags$p(" Select the statistical method to train the genomic prediction model and predict the values of target lines. rrBLUP method is implemented using the
                                           rrBLUP (Endelman 2011) package. Expectation maximization based RR-BLUP, BayesB and BayesLASSO methods are implemented using the bWGR package (Xavier et al. 2020).
                                           "))),
                           tags$br(),
                           tags$br(),
                           tags$h4(textOutput("RankedLinesHeader")),
                           tableOutput("Ranked_Lines_for_SelectionST"),
                           uiOutput("scatter")
                        )
                       )
                     ),
                     
                     tabPanel("Multi-trait",

                              sidebarLayout(
                                sidebarPanel(
                                  selectInput(inputId="TrainSetMT","Choose Training Set",c("Complete Input Genotype Set","Optimal Train Set from Step 3","Random Train Set from Step 3")),
                                  tags$br(),
                                  tags$br(),
                                  selectInput(inputId="GPModelMT","Choose Prediction Model for Multiple Traits",c("BRR (BGLR)","RKHS (BGLR)","Spike-Slab(BGLR)","Mmer (Sommer)")),
                                  actionButton("RunPredictionsMT", "Predict Multiple Traits!"),
                                  tags$br(),
                                  tags$br(),
                                  downloadButton("ExportOutMT", "Export Output Table")

                                ), mainPanel(
                                  fluidRow(
                                    column(width=8,tags$h3(tags$strong("Train Multi-trait Genomic Prediction Model"))), 
                                    #actionButton("GP_ST", "prev"),
                                    #actionButton("GP_MTME", "next")
                                  ),
                                  fluidRow(
                                    column(width=10, tags$p("Select the statistical method to train the genomic prediction model and predict the values of target lines.Multi-trait predictions are implemented using the BGLR and Sommer packages. The Multitrait function in BGLR implements Bayesian Ridge Regression, RKHS, and Spike-Slab methods (Perez-Rodriguez P and de Los Campos G, 2022).The multi-trait GBLUP is implemented using the 'mmer' function in 'sommer' package (Covarrubias-Pazaran G 2016)."))),
                                  tags$br(),
                                  tags$br(),
                                  #tags$h4(textOutput("RankedLinesHeader")),
                                  tableOutput("Ranked_Lines_for_SelectionMT"),
                                  uiOutput("scatterMT")

                                )
                              )
                     )
                     
                     # tabPanel("Multi-environment",
                     # 
                     #      sidebarLayout(
                                # sidebarPanel(
                                # selectInput(inputId="TrainSetME","Choose Training Set",c("Complete Input Genotype Set","Optimal Train Set from Step 3","Random Train Set from Step 3")),
                                # tags$br(),
                                # tags$br(),
                                # selectInput(inputId="GPModelMT","Choose Prediction Model for Multiple Traits",c("BRR (BGLR)","RKHS (BGLR)","Spike-Slab(BGLR)","Mmer (Sommer)")),
                                # actionButton("RunPredictionsMT", "Predict Multiple Traits!"),
                                # tags$br(),
                                # tags$br(),
                                # downloadButton("ExportOutMT", "Export Output Table")
                                #  
                     #            sidebarPanel(
                     #              selectInput(inputId="Year","Select Year/Years",NULL),
                     #              tags$br(),
                     #              tags$br(),
                     #              selectInput(inputId="Location","Select Location/Locations",NULL),
                     #              tags$br(),
                     #              tags$br(),
                     #              tags$br(),
                     #              selectInput(inputId="fixed","Choose Fixed Effect",choices=NULL),
                     #              tags$br(),
                     #              tags$br(),
                     #              tags$br(),
                     #              selectInput(inputId="MEModel","Choose MultiEnvironment Model ",c("Main Effect","Homogeneous Variance (CS)","Heterogeneous Variance (CS+DG)","Unstructured (US)")),
                     #              actionButton("RunPredictionsME", "Fit Multi-environmental Model!"),
                     #              tags$br(),
                     #              tags$br(),
                     #              tags$br(),
                     #              tags$br()
                     # 
                     #            ), mainPanel(
                     # 
                     #              fluidRow(
                     #                column(width=10, tags$p("Environments are defined as combinations of year and locations in which phenotypic data are collected.
                     #                Multi-environmental models are implemented using the 'mmer' function in 'sommer' package.
                     #                Fit genomic prediction models taking into account only the main effects or the main effects + GxE effects using the 'sommer' package'.
                     #                The user needs to select one or many years and locations and the type of variance-covariance structure for fitting the model"))),
                     #              tags$br(),
                     #              tags$br(),
                     #              tags$h4(textOutput("RankedLinesHeader")),
                     #              tableOutput("Ranked_Lines_for_SelectionSTME"),
                     #              tableOutput("Ranked_Lines_for_SelectionMTME")
                     # 
                     #            )
                     #          )
                     # )
                  )
              )
      )

  )             
                     
         
###

server <- function(input,output,session){
  
  options(shiny.silent.errors=FALSE)
 
  #Geno
  
  # Geno <- reactive({
  #   
  #   genoFile <- input$infileVCF
  #   ext <- tools::file_ext(genoFile$datapath)
  #   req(genoFile)
  #   validate(need(ext == "vcf", "Please upload a vcf file"))
  #   
  #  # withProgress(message = 'Reading Data', value = 0, {
  #  #   NUST_Genotypes_VCF <- read.table(genoFile$datapath)})
  #   
  #   # withProgress(message = 'Converting VCF to Dataframe', value = 0, {
  #   #   gt2d <- VCFtoDF(genoFile$datapath) }) 
  #   
  #   
  #    gt2d
  # })
  
  GenoTas <- reactive({
    genoFile <- input$infileVCF
    ext <- tools::file_ext(genoFile$datapath)
    req(genoFile)
    validate(need(ext == "vcf", "Please upload a vcf file"))
    withProgress(message = 'Reading Data', value = 0, {
    getTasObj(genoFile$datapath)})
  })
  
  Geno <- eventReactive(input$infileVCF,{ 
    withProgress(message = 'Converting VCF to Dataframe', value = 0, {
    gt2d <- getGenoTas_to_DF(GenoTas())})
   gt2d
  
  })
  
  
  
  #markerSet <- reactive(Geno()[,"SNPID"])
  
  #Pheno
  Pheno <- reactive({
    
    phenoFile <- input$infileBLUEs
    
    ext <- tools::file_ext(phenoFile$datapath)
    req(phenoFile)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(phenoFile$datapath, header = input$header)
    
  })
  
  observeEvent(input$infileBLUEs, {
    updateSelectInput(inputId = "trait",choices = colnames(Pheno())[2:ncol(Pheno())])})
  
  observeEvent(input$infileBLUEs, {
    updateSelectInput(inputId = "fixed",choices = c("NULL",colnames(Pheno())[2:ncol(Pheno())]),selected="NULL")})
  
  
#Target
  TargetTab <- reactiveVal(NULL)
  
 
  # eventReactive(input$infileTargetTable,{
  #  
  #   TargetFile <- inTarget()
  #   
  #   ext <- tools::file_ext(TargetFile$datapath)
  #   req(TargetFile)
  #   validate(need(ext == "csv", "Please upload a csv file"))
  #   
  #   TargetTab(read.csv(TargetFile$datapath, header = input$header))
  # 
  #  })
  
  observeEvent(input$infileTargetTable, {
    # Triggered when a file is uploaded
    req(input$infileTargetTable) # Ensure a file is uploaded
    
    TargetFile <- input$infileTargetTable
    
    ext <- tools::file_ext(TargetFile$datapath)
    # Stop execution and show error if the file is not a CSV
    if (ext != "csv") {
      shiny::showNotification("Please upload a CSV file", type = "error")
      TargetTab(NULL) # Reset or clear previous data
      return() # Stop further execution
    }
    
    # If execution reaches here, it means a CSV file is uploaded
    tryCatch({
      # Attempt to read the CSV
      dataTable <- read.csv(TargetFile$datapath, header = input$header)
      TargetTab(dataTable) # Store the read data
    }, error = function(e) {
      # In case of error reading the file, show notification and reset TargetTab
      shiny::showNotification("Error reading the CSV file", type = "error")
      TargetTab(NULL)
    })
  })
#####
  
  phenoHead <- eventReactive(input$infileBLUEs,{ paste("Phenotype Table with ",nrow(Pheno())," lines and ",ncol(Pheno())-1," traits",sep="")})
  output$PhenoHeader <- renderText({phenoHead()})
  output$PhenoTable <- renderTable({as.data.frame((Pheno())[1:5,1:5])})
  
  genoHead <- eventReactive(input$infileVCF,{ paste("Genotype Table with ",ncol(Geno())-5," lines and ",nrow((Geno()))," markers",sep="")})
  output$GenoHeader <- renderText({genoHead()})
  output$GenoTable <- renderTable({as.data.frame((Geno())[1:5,1:5])})
  
  
  
  # TargetHead <- eventReactive(input$infileTargetTable,{
  #    if(!is.null(TargetTab())){
  #       paste("Table with information on ",nrow(TargetTab())," Target lines",sep="")
  #       output$TargetHeader <- renderText({TargetHead()})
  #       
  #       output$TargetTable <-  renderTable({
  #         TargetTableOut <- as.data.frame(TargetTab()[1:5,]) 
  #         colnames(TargetTableOut) <- ""
  #         TargetTableOut
  #       })
  #     }
  #  })
  
  TargetHead <- eventReactive(input$infileTargetTable, {
    if (!is.null(TargetTab())) {
      paste("Information on", nrow(TargetTab()), "Target lines")
    }
  })
  output$TargetHeader <- renderText({
    req(TargetHead()) # Ensure that TargetHead has a value before proceeding
    TargetHead()
  })
  
  #### Target Table 
  
  TargetIDs <- reactive(NULL)
  
  nTraits <- eventReactive(input$infileBLUEs,{(ncol(Pheno())-1)})
  
  TargetIDs <- eventReactive(input$infileTargetTable,{
          as.character(TargetTab()[,1])})
  
  
  TargetHead <- eventReactive(input$infileTargetTable,{ paste("Table with information on ",nrow(TargetTab())," Target lines",sep="")})
  output$TargetHeader <- renderText({TargetHead()})
  
  output$TargetTable <-  renderTable({
    
    TargetTableOut <- as.data.frame(TargetTab()[1:5,]) 
    colnames(TargetTableOut) <- ""
    TargetTableOut
  })
  
  output$TargetTable <- renderTable({
    req(TargetTab()) # Ensure that TargetTab is not NULL
    TargetTableOut <- as.data.frame(TargetTab()[1:5,]) 
    colnames(TargetTableOut) <- "" # Clear column names if needed
    TargetTableOut
  })
 
####  

  
  TargetIDs <- reactiveVal(NULL)
 
  nTraits <- eventReactive(input$infileBLUEs,{(ncol(Pheno())-1)})
  # TargetIDs <- eventReactive(input$infileTargetTable,{
  #         as.character(TargetTab()[,1])})
  # 
  observeEvent(input$infileTargetTable,{
           TargetIDs(as.character(TargetTab()[,1]))})
  
## Trait 
  
  Trait <- reactive(input$trait)
  
  summaryHead <- eventReactive(input$trait,{paste("Summary of distribution of selected phenotypic values")})
  output$summaryHeader <- renderText({summaryHead()})
  nSelTraits <- reactive(length(Trait()))
  
  SummaryTxt <- eventReactive(input$trait,{ 
    do.call(rbind,lapply(Trait(),function(x) round(summary(Pheno()[,x]),digits=2)))
  })
  
  observeEvent(input$trait,{ 
  output$Summary <- renderTable({
   
     trTable <- cbind(unlist(Trait()),nSelTraits(),SummaryTxt())
     print(trTable)
    })
  })
  
  # boxHead <- eventReactive(input$trait,{paste("Boxplot of selected phenotypic values")})
  # output$boxHeader <- renderText({boxHead()})
  # 
  

 

  #%>% bindCache(t_Geno(),Pheno(),TargetIDs(),unlist(Trait())) %>% bindEvent(input$trait,input$infileVCF,input$infileTargetTable,input$infileBLUEs)
 
  
### Filter Data

### Filter 1  
  observeEvent(input$infileVCF, {
    updateNumericInput(inputId = "siteMinCnt",value= round(0.8*ncol(Geno()),digits=0),min = round(0.1*ncol(Geno()),digits=0), max=ncol(Geno()))}) 
  
    
  siteMinCnt <- reactive(input$siteMinCnt)
  MAF <- reactive(input$MAF)
  
  
  # GenoTas <- reactive({
  #   genoFile <- input$infileVCF
  #   ext <- tools::file_ext(genoFile$datapath)
  #   req(genoFile)
  #   validate(need(ext == "vcf", "Please upload a vcf file"))
  #   getTasObj(genoFile$datapath)
  # })
  

  GenoFilt1 <- eventReactive(input$FilterSites,{
    withProgress(message = 'Filtering Sites', value = 0, {
      getFilteredSitesGenoData(GenoTas(),siteMinCnt(),MAF())}) 
  })
  
### Filter 2
  
  minNotMissing <- reactive(input$minNotMissing)
  
  setTasGenoFilt1 <- reactive(input$setGenoFilt1Tas)
 
 
  GenoFilt2 <- eventReactive(input$FilterTaxa,{
    
    if(setTasGenoFilt1()== FALSE){        
      withProgress(message = 'Filtering Taxa', value = 0, {getFilteredTaxaGenoData(GenoTas(),minNotMissing()) })
    }else if(setTasGenoFilt1()== TRUE){
      withProgress(message = 'Filtering Taxa', value = 0, {getFilteredTaxaGenoData(GenoFilt1(),minNotMissing())})
    }
  })
  
   setTasGenoFilt2 <- reactive(input$setGenoFilt2Tas)

## Tas to DF
   # Geno_DF <- reactive(getGenoTas_to_DF(GenoTas(),Geno()))
   # GenoFilt1_DF <- reactive(getGenoTas_to_DF(GenoFilt1(),Geno()))
   # GenoFilt2_DF <- reactive(getGenoTas_to_DF(GenoFilt2(),Geno()))

   # Geno_DF <- reactive(getGenoTas_to_DF(GenoTas()))
   # GenoFilt1_DF <- reactive(getGenoTas_to_DF(GenoFilt1()))
   # GenoFilt2_DF <- reactive(getGenoTas_to_DF(GenoFilt2()))

   Geno_DF <- observeEvent(input$infileVCF, {getGenoTas_to_DF(GenoTas())})
   GenoFilt1_DF <- eventReactive(input$FilterSites,{getGenoTas_to_DF(GenoFilt1())})
   GenoFilt2_DF <- eventReactive(input$FilterTaxa,{getGenoTas_to_DF(GenoFilt2())})
   
   
   # Geno_DF <- reactive(as.data.frame(as.matrix(GenoTas())))
   # GenoFilt1_DF <- reactive(as.data.frame(as.matrix(GenoFilt1())))
   # GenoFilt2_DF <- reactive(as.data.frame(as.matrix(GenoFilt2())))
   # 

####  
  genoFiltHead <- eventReactive(input$FilterSites,{ paste("Filter Sites: Genotype Table with ",ncol(GenoFilt1_DF())-5," lines and ",nrow((GenoFilt1_DF()))," markers",sep="")})
  output$GenoFiltHeader <- renderText({genoFiltHead()})
  output$FilteredGenoTable <- renderTable({as.data.frame((GenoFilt1_DF())[1:5,1:6])})
  
###
  genoFiltHead2 <- eventReactive(input$FilterTaxa,{ paste("Filter Taxa: Genotype Table with ",ncol(GenoFilt2_DF())-5," lines and ",nrow((GenoFilt2_DF()))," markers",sep="")})
  output$GenoFiltHeader2 <- renderText({genoFiltHead2()})
  output$FilteredGenoTable2 <- renderTable({as.data.frame((GenoFilt2_DF())[1:5,1:6])})
  
#####

### Imputation

  FiltGeno <- reactive({ 
      if(setTasGenoFilt1()== FALSE  & setTasGenoFilt2()== FALSE){        
        GenoTas()
      }else if(setTasGenoFilt1()== TRUE  & setTasGenoFilt2()== FALSE){        
        GenoFilt1()
      }else if(setTasGenoFilt2()== TRUE){
        GenoFilt2()
      }
      
  })
   
   l<- reactive(input$l)
   k <- reactive(input$k)
   
   nN <- reactive(input$nN)
   Dist <- reactive(input$Dist)
   impMethod <- reactive(input$imputeMet)
  
   
  GenoImp <-  eventReactive(input$Impute,{
    
    withProgress(message = 'Imputing Genotypic Scores', value = 0, {
     
       print(impMethod())
       if(impMethod()=='Numeric'){ 
       
        df<- getImputedData_Num(FiltGeno(),nN(),Dist())
         
       }else if(impMethod()=='LDKNNI'){ 
        df <-  getImputedData_LDKNNI(FiltGeno(),l(),k())
         
       }
      
       df
    })
     
   },ignoreNULL = TRUE)
   
  
   GenoImp_DF1 <- eventReactive(input$Impute,{
                getGenoTas_to_DF(GenoImp())
    })
   
  
### if you use only reactive, this will throw an error ncol(GenoImp_DF())-5
  
   #genoImpHead <- eventReactive(input$Impute,{paste("Genotype Table with ",ncol(GenoImp_DF())-5," lines and ",nrow((GenoImp_DF()))," markers",sep="")})
   
  # genoImpHead <- reactive({paste("Genotype Table with ",ncol(GenoImp_DF())-5," lines and ",nrow((GenoImp_DF()))," markers",sep="")})
   
   
   genoImpHead <- reactive({
     # Ensure GenoImp_DF is only called when it has updated data
     df <- as.data.frame(GenoImp_DF1()) # This call ensures dependency
     
     if(is.data.frame(df) && ncol(df) > 5) {
   #     # Successful case
       paste("Genotype Table with ", ncol(df) - 5, " lines and ", nrow(df), " markers", sep = "")
     } else {
   #     # Fallback or error message
       "Waiting for data"
     }
    })
   # 
   # UI and server logic to render the header text and imputed genotype table
   output$GenoImpHeader <- renderText({
     genoImpHead()  # This will now wait for GenoImp_DF to complete
   })
   
  
   # output$GenoImpHeader <- renderText({genoImpHead()})
    output$ImputedGenoTable <- renderTable({as.data.frame((GenoImp_DF1())[1:5,1:10])})
   # 
### 
     
   observeEvent(input$imputeMet,{ 
   
     if(impMethod()=="AlphaPlantImpute"){ 
     
     library(reticulate)
     reticulate::virtualenv_create("pyEnv", python = "3.10.0")
     reticulate::virtualenv_install("pyEnv",c("numpy","pandas"))
     reticulate::use_virtualenv("pyEnv", required = TRUE)
     reticulate::py_run_string("import sys")
     
     py_module_available("alphaplantimpute2")
   
     output$APIUI <- renderUI({
     # Check if AlphaPlantImpute is selected
     
       fluidPage(
         
         fluidRow(
           #column(1),
           column(width=4,tags$strong(tags$h3("Build Haplotype Library "))),
           column(7),column(width=4,tags$strong(tags$h3("Impute Genotype Table ")))
         ),
         
         tags$br(),
         
         fluidRow(
           #column(1),
           column(width=4,
                            numericInput(inputId="nHap",label = "Enter number of haplotypes",value=20,min=2,max=100),
           ),
           column(7),column(width=5,
                            tags$strong(("Input founders file for pedigree based imputation")))
         ),
         tags$br(),
         
         fluidRow(
           #column(1)
           column(width=4,
                            numericInput(inputId="nSampRnds",label = "Enter number of sample rounds",value=5,min=2,max=100),
           ),
           column(7),column(width=4, 
                            fileInput("founder_file", "Select Founder Files", multiple = FALSE,accept = ".txt")),
         ),
         tags$br(),
         
         fluidRow( 
           #column(1),
           column(width=3,
                            numericInput(inputId="HDthresh",label = "Enter non-missing threshold for high density genotypes",value=0.9,min=0.5,max=1),
           )),
         tags$br(),
         
         # Build library button
         fluidRow( 
           column(2),column(width=3,
                            actionButton("build_library", "Build Library"),
           ),
           column(7),column(width=3,
                            actionButton("impute_APIdata", "Impute Data")
           ),
         ),   
         tags$br(),
         tags$br(),
         
         tags$head(
           tags$style(HTML("
      #message {
        
        /*max-height: 1200px; Set maximum height */
        overflow-y: scroll; /* Enable vertical scrolling */
        overflow-x: scroll;  /*Hide horizontal scrolling */
        overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
        width: 400px; 
        /*max-width: 100%; */
        padding: 6px 12px;
        height: 300px;
      }
    "))
         ),
         
         tags$head(
           tags$style(HTML("
      #message2 {
        
        overflow-y: scroll; /* Enable vertical scrolling */
        overflow-x: scroll;  /*Hide horizontal scrolling */
        overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
        width: 400px; 
        /*max-width: 100%; */
        padding: 6px 12px;
        height: 300px;
      }
    "))
         ),
         fluidRow(
           
           column(width = 5,
                  verbatimTextOutput("message")
           ),
           column(width = 1),  # Adjust the width or remove if not needed
           column(width = 5,
                  verbatimTextOutput("message2")
           )
         ),
         
   
   
       )
    })
     
  } 
      
  })
      
 
   
  observeEvent(input$imputeMet,{
    if(impMethod()=="AlphaPlantImpute" && (!is.null(FiltGeno()))){
     getGenoData_API(FiltGeno())
    }
  })
   
####  
   
   founder <- reactive({
     
     # Attempt to retrieve file information
     founderFile <- input$founder_file
     
     # Proceed only if a file is uploaded
     if (!is.null(founderFile)) {
       # Extract file extension
       ext <- tools::file_ext(founderFile$datapath)
       
       # Validate file type
       validate(need(ext == "txt", "Please upload a txt file"))
       
       # Return the file path if the file is valid
       return(founderFile$datapath)
     }
     
     # Return NULL if no file is uploaded
     return(NULL)
   })
   
   #### 
   temp_file <- reactiveVal()
   nHap <- reactive(input$nHap)
   nSampRnds <- reactive(input$nSampRnds)
   HDthresh <- reactive(input$HDthresh)
   
   observeEvent(input$build_library,{
     
     # Extract necessary input values
     
     # Path for the temporary file
     temp_file(tempfile())
     
     # Start the process in a separate R process
     rProcess <- callr::r_bg(function(nHap, nSampRnds, HDthresh, temp_file) {
       library(reticulate)
       sys <- import("sys")
       api2 <- import("alphaplantimpute2.alphaplantimpute2")
       
       sys$argv <- c('alphaplantimpute2', '-createlib', '-out', 'lib',
                     '-genotypes', 'current_GenoTable.genotypes',
                     '-n_haplotypes', as.character(nHap),
                     '-n_sample_rounds', as.character(nSampRnds),
                     '-hd_threshold', as.character(HDthresh),
                     '-seed', '42')
       
       
       sink(temp_file)
       api2$main()
       sink()
     }, args = list(nHap(), nSampRnds(), HDthresh(), temp_file()), stdout = temp_file(), stderr = temp_file())
     
     ### Test output 
     
     processAlive <- reactive(rProcess$is_alive())
     # Periodically read the file and update the UI
     output$message <- renderText({
       
       invalidateLater(1000, session) # Update every second
       
       if(file.exists(temp_file())){
         lines <- readLines(temp_file(), warn = FALSE)
         return(paste(lines, collapse = "\n"))
       }else {
         return("Waiting for output...")
       }
       
       if(!processAlive()){
         
         if (file.exists(temp_file())) {
           lines <- readLines(temp_file(), warn = FALSE)
           txt <- paste(lines,collapse="\n")
           return(paste0(txt,"\n","Build Completed"))
         }else {
           return("Waiting for output...")
         }
       }
     })
     
   })
   
##
   
   temp_file2 <- reactiveVal()
   processComplete2 <- reactiveVal(FALSE)
   ouputRead <- reactiveVal(FALSE)
   
   
   # Start the background process and return rProcess
   rProcess2 <- eventReactive(input$impute_APIdata,{
     
     library(reticulate)
     sys <- import("sys")
     api2 <- import("alphaplantimpute2.alphaplantimpute2")
     
     ## Extract necessary input values
     
     founder_file <- founder()
     
     # Path for the temporary file
     
     temp_file2(tempfile())
     
     # Start the process in a separate R process
     callr::r_bg(function(founder_file,temp_file){
       library(reticulate)
       sys <- import("sys")
       api2 <- import("alphaplantimpute2.alphaplantimpute2")
       
       if(is.null(founder_file)){
         sys$argv <- c('alphaplantimpute2', '-impute', '-out', 'imputed_out',
                       '-genotypes', 'current_GenoTable.genotypes',
                       '-libphase','lib.phase'
         )
       }else {
         sys$argv <- c('alphaplantimpute2', '-impute', '-out', 'imputed_out',
                       '-genotypes', 'current_GenoTable.genotypes',
                       '-founders', as.character(founder_file),
                       '-libphase','lib.phase'
         )
       }
       
       sink(temp_file)
       api2$main()
       # print("Imputing")
       # Sys.sleep(2)
       sink()
     }, args = list(founder_file,temp_file2()), stdout = temp_file2(), stderr = temp_file2())
     
   })
   
   
   # Reactive expression to check process status and read temp file
   processStatus2 <- reactive({
     invalidateLater(1000,session) # Update every second
     
     if (!is.null(rProcess2()) && rProcess2()$is_alive() && !ouputRead()){
       if (file.exists(temp_file2())) {
         lines <- readLines(temp_file2(), warn = FALSE)
         return(paste(lines,collapse = "\n"))
       }else {
         return("Waiting for output...")
       }
     }else if(!is.null(rProcess2()) && !rProcess2()$is_alive() && !ouputRead()){
       
       if (file.exists(temp_file2())) {
         lines <- readLines(temp_file2(), warn = FALSE)
         txt <- paste(lines, collapse = "\n")
         processComplete2(TRUE)
         return(paste(txt,"Imputation Completed",collapse="\n"))
         
       }else {
         return("Waiting for output...")
       }
      
     }else if(!is.null(rProcess2()) && !rProcess2()$is_alive() && ouputRead()){
       
       if (file.exists(temp_file2())) {
         lines <- readLines(temp_file2(), warn = FALSE)
         txt <- paste(lines, collapse = "\n") 
         
         return(paste0(txt,"\n","Imputation Completed","\n", "Imputed Geno Table has Dim:",paste(dim(GenoImp_DF()),collapse=" ")))
       }else{
         return("Output not read")
       }
     }else{return("Var Error")}
   })
   
   # Update the UI with process status
   output$message2 <- renderText({
     processStatus2()
   })
   
   
   # Process and return the output once the background process is complete
   observe({
     # Ensure the background process is completed
     
     req(processComplete2())
     ouputRead(TRUE)
   })
   
   ### Read Imputed genotypic data and prepare output 
   
   # #impGeno <- reactive({
   GenoImp_DF2 <- reactive({
     # Ensure outputRead is set to TRUE
     req(ouputRead())
     # ... Your existing code to process and return the output ...
     
     ImpGeno <- getImpGenoData_API()
     ImpGeno
   })

   GenoImp_DF <- reactive({
     if (!is.null(GenoImp_DF1()) && nrow(GenoImp_DF1()) > 0) {
       # GenoImp_DF1 is not NULL and has data
       return(GenoImp_DF1())
     } else if (!is.null(GenoImp_DF2()) && nrow(GenoImp_DF2()) > 0) {
       # GenoImp_DF2 is not NULL and has data, and GenoImp_DF1 was NULL or had no data
       return(GenoImp_DF2())
     } else {
       # Both are NULL or have no data, return NULL or a default value/data frame
       return(NULL)  # Or any default value you deem appropriate
     }
   })
  
   
   # genoImpHead <- reactive(paste("Genotype Table with ",ncol(GenoImp_DF())-5," lines and ",nrow((GenoImp_DF()))," markers",sep=""))
   # GenoImp_DF <- eventReactive(input$Impute,{getGenoTas_to_DF(GenoImp(),Geno())})
   #
   #GenoImp_DF <- reactive(getGenoTas_to_DF(GenoImp(),Geno()))
   #GenoImp_DF <- reactive(getGenoTas_to_DF(GenoImp()))
   #GenoImp_DF <- reactive(as.data.frame(as.matrix(GenoImp())))
   
   
   
### Merge and Process Data
#
   
   
   setTasImpGeno <- reactive(input$setGenoImpTas) 
   
### Check for which genotype matrix to set   
   
   GenoPre <- reactive({
     
     if( setTasGenoFilt1() == FALSE  && setTasGenoFilt2() == FALSE && setTasImpGeno()== FALSE ){        
       Geno_DF()
     }else if(setTasGenoFilt1() == TRUE  && setTasGenoFilt2() == FALSE && setTasImpGeno()== FALSE){        
       GenoFilt1_DF()
     }else if(setTasGenoFilt2() == TRUE && setTasImpGeno()== FALSE){
      GenoFilt2_DF()
     }else if(setTasImpGeno()== TRUE){
       GenoImp_DF()
     }
     
   })
   
   
### Merge geno and pheno data    
  
   
  mergedData <-  eventReactive(input$trait,{
    t_Geno <- reactive(GenoPre())
    withProgress(message = 'Merging Pheno and Geno Data', value = 0, {
     tryCatch({
        # Simulate data merging process
        #browser()
        result <- getMergedData(t_Geno(),Pheno(),TargetIDs())

        if (length(result) != 2) {
          # Throw a custom error if the length condition is not met
          stop("The merged data must have exactly 2 elements.")
        }
        return(result)
      }, error = function(e) {
        # Print error and stop the execution
        print(e)
        stop(e)
      })
    })
  })


###
   
   processedData <- reactive({
     prData <-  getProcessedData(mergedData(),Trait())
     validate(
       need(length(mergedData()) == 2, "Error: The merged data does not have exactly 2 elements.")
     )
     return(prData)
   })
 
###    
     
   observeEvent(input$trait, {
     if(length(processedData())==2){
      updateNumericInput(inputId = "noCandidates",value= nrow(processedData()[[1]]),min =2, max=nrow(processedData()[[1]]))
    }
   })
   
   observeEvent(input$trait,{
     if(length(processedData())==2){
      updateNumericInput(inputId = "noToSelect",value= nrow(processedData()[[1]]),min =2, max=nrow(processedData()[[1]]))
    }
   })
   
   mergedData <-  reactive({
  # mergedData <- eventReactive(input$trait,{
     t_Geno <- reactive(GenoPre())
     withProgress(message = 'Merging Pheno and Geno Data', value = 0, {
       getMergedData(t_Geno(),Pheno(),TargetIDs())
     }) 
   })
   
###   
   processedData <- reactive(getProcessedData(mergedData(),Trait()))
   

   
   # observeEvent(input$trait, {
   #   updateNumericInput(inputId = "noCandidates",value= 250,min =2, max=nrow(processedData()[[1]]))}) 
   # 
   # observeEvent(input$trait,{
   #   updateNumericInput(inputId = "noToSelect",value= 100,min =2, max=nrow(processedData()[[1]]))})
   # 
   

  
   # MssgTargetSet <- eventReactive(input$trait,{
   #   if(nrow(processedData()[[2]])==0){
   #     "The current genotypic file doesn't have target line genotypes. 
   #      Load genotypic data file with Target IDs"
   #   }
   # })
   # 
   # output$MssgTarget <- renderUI({
   #   message <- MssgTargetSet()
   #   HTML(paste0("<div class='message-text'>", message, "</div>"))
   # })
   # 
    

   MssgTargetSet <- eventReactive(input$trait,{
     if(nrow(processedData()[[2]])==0 | is.null(processedData()[[2]])){
       response <- "The current genotypic file doesn't have target line genotypes. 
        Load genotypic data file with Target IDs"
       return(response)
     }
    
   })
   
   output$MssgTarget <- renderUI({
     message <- MssgTargetSet()
     HTML(paste0("<div class='message-text'>", message, "</div>"))
   })
   
   
   
   
###    
     
   observeEvent(input$trait, {
     updateNumericInput(inputId = "noCandidates",value= nrow(processedData()[[1]]),min =2, max=nrow(processedData()[[1]]))}) 
   
   observeEvent(input$trait,{
     updateNumericInput(inputId = "noToSelect",value= 100,min =2, max=nrow(processedData()[[1]]))})
   
   

# TS Optimization
  
  noCandidates <- reactive(input$noCandidates)
  nTrainToSelect <- reactive(input$noToSelect)
  optimCriteria <- reactive(input$optCriteria)
  
  predictionData <- reactive({getPredictionData(processedData(),noCandidates())}) 

  GAParameters <- reactiveValues(npop=100,nelite=10,mutprob=0.5,mutintensity=1,niterations=100,minitbefstop=50,tabu="TRUE",tabumemsize=1,plotiters="FALSE",errorstat="PEVMEAN2",lambda=1e-6,mc.cores=10)
  
  # GAParameters$npop= reactive(input$noPop)
  # GAParameters$nelite=reactive(input$noElite)
  # GAParameters$mutprob=reactive(input$mutProb)
  # GAParameters$mutintensity=reactive(input$mutInt)
  # GAParameters$niterations=reactive(input$nIterGA)
  # GAParameters$tabu= reactive(input$tabu)
  # GAParameters$tabumemsize= reactive(input$tabumemsize)
  # GAParameters$plotiters=reactive(FALSE)
  # GAParameters$errorstat= reactive(input$optCriteria)
  # GAParameters$lambda=reactive(input$lambda)
  # GAParameters$mc.cores= reactive(input$mcCores)
  
  # 
  Train_STPGA <-   eventReactive(input$Optimize,{
    withProgress(message = 'Running Optimizations', value = 0, {
      getOptimalTS(predictionData(),unlist(Trait()),nTraits(),noCandidates(),nTrainToSelect(),isolate(reactiveValuesToList(GAParameters)))})
  })
    
  #%>% bindCache(processedData(),Trait(),nTraits(),noCandidates(),nTrainToSelect())

  Train_Random <-   eventReactive(input$Optimize,{
    withProgress(message = 'Running Random Set', value = 0, { 
      getRandomTS(predictionData(),unlist(Trait()),nTraits(),noCandidates(),nTrainToSelect())})
  })
    
  #%>% bindCache(processedData(),Trait(),nTraits(),noCandidates(),nTrainToSelect())
 
  
 #  TSOptOutputList <- reactive({
 #    if(nSelTraits()==1){
 #      withProgress(message = 'Performing CrossValidations', value = 0, {
 #        getTSComparisons(predictionData(),Train_STPGA(),Train_Random(),unlist(Trait()),nTraits(),optTS())
 #      })
 #    }else if(nSelTraits()>1){
 #      withProgress(message = 'Performing CrossValidations', value = 0, {
 #        getTSComparisonsMT(predictionData(),Train_STPGA(),Train_Random(),unlist(Trait()),nTraits(),optTS())
 #      })
 #    }
 # })
  
 
      
  TSOptOutputList <- eventReactive(input$Optimize,{
      if(nSelTraits()==1){
       withProgress(message = 'Performing CrossValidations', value = 0, {
        getTSComparisons(predictionData(),Train_STPGA(),Train_Random(),unlist(Trait()),nTraits(),TargetIDs())
       })
      }else if(nSelTraits()>1){
       withProgress(message = 'Performing CrossValidations', value = 0, {
        getTSComparisonsMT(predictionData(),Train_STPGA(),Train_Random(),unlist(Trait()),nTraits(),TargetIDs())
       })
      }

    })
  
  tsOptHead <- eventReactive(input$Optimize,{paste("Correlation between observed and predicted values")})
  output$tsOptHeader <- renderText({tsOptHead()}) 
  output$PredAccuracyTable <- renderTable({
  TSTable <- TSOptOutputList()
  
  colnames(TSTable) <- c("Ridge Regression","BayesB","Bayes LASSO")
  rownames(TSTable) <- c("Optimal Set","Random Set")
  TSTable},
  colnames=TRUE,rownames=TRUE,digits=2)
  
  
 
  
## CVR
 
  k <- reactive(input$k)
  nIter <- reactive(input$nIter)
  
### Conflicting statement messages
 
  
  Mssg2CVR <- eventReactive(input$CrossValidationMT,{
    if(nSelTraits()==1){
      "Select more than one trait for cross validation of multi-trait models"
    }
  })
  
  output$Mssg <- renderText({Mssg2CVR()})
  
  
 
  cvrOutputListMT <- eventReactive(input$CrossValidationMT,{ withProgress(message = 'Running CrossValidations', value = 0, {
        getMTCVR(predictionData(),unlist(Trait()),nTraits(),k(),nIter()) })
  })

###
  
  cvrHead <- eventReactive(input$CVEvent,{
  paste("Correlation between observed and predicted values for ",paste(Trait(),collapse=" and "),sep="")})
  
  output$cvrHeader <- renderText({cvrHead()}) 
 
  cvrOutputListST <- eventReactive(input$CrossValidationST,{
  
    if(nSelTraits()==1){
       PATab <- withProgress(message = 'Running CrossValidations', value = 0, {
         getemCVR(predictionData(),unlist(Trait()),nTraits(),k(),nIter())})
       PATable <- round(PATab[c("emRR","emBB","emBL")],digits=2) 
       names(PATable)<- c("Ridge Regression","BayesB","Bayes LASSO")
       PATable <- rbind.data.frame(names(PATable),PATable)
       colnames(PATable) <- c("Ridge Regression","BayesB","Bayes LASSO")
       rownames(PATable) <- c("Model","Accuracy")
       return(PATable[2,])
    }
    if(nSelTraits()>1){
       PATableComb <- c()
       for(nSelT in 1:nSelTraits()){
            i <- reactive(nSelT)
           PATab <- withProgress(message = 'Running CrossValidations', value = 0, {
                    getemCVR(predictionData(),unlist(Trait())[i()],nTraits(),k(),nIter())})
           PATable <- round(PATab[c("emRR","emBB","emBL")],digits=2)
           #PATable2 <- c(unlist(Trait())[i()],PATable)
           PATableComb <- rbind(PATableComb,PATable)
           
       }
       #PATableComb <- rbind(c("Trait","Ridge Regression","BayesB","Bayes LASSO"),PATableComb)
       #colnames(PATableComb) <- rep("",ncol(PATableComb))
       colnames(PATableComb) <- c("Ridge Regression","BayesB","Bayes LASSO")
       rownames(PATableComb) <- unlist(Trait())
       return(PATableComb)
    }
    
  })
     

  output$emCVRST <- renderTable({
   
       cvrOutputListST()
     
   },colnames=TRUE,rownames=TRUE)
  
    
  output$emCVRMT <- renderTable({

    if(nSelTraits()>1){
       PATable <- cvrOutputListMT()
       #colnames(PATable) <- rep("",ncol(PATable))
       print.data.frame(as.data.frame(PATable))
    }

  },colnames=TRUE,rownames=TRUE)
  
  
  observeEvent(input$CrossValidationST, {
    updateTextInput(session,"CVEvent", value = "ST")
  })
  
  observeEvent(input$CrossValidationMT,{
    updateTextInput(session,"CVEvent", value = "MT")
  })
  
#### Optimized Training Sets
  
  TS <- reactive(input$TrainSet) 
  

  optTS <-  reactive({ 
    if(TS() == "Complete Input Genotype Set"){
      optimTS <- NULL
    }else if(TS() == "Optimal Train Set from Step 3"){
      
      optimTS <- Train_STPGA()[[2]]
    }else if(TS() == "Random Train Set from Step 3"){
      
      optimTS <- Train_Random()[[1]]
      
    }
    (optimTS)
  })
  
## GP  Models   /// & Fixed()!= "NULL" ///| Fixed()== "NULL"
  
  GPModelST <- reactive(input$GPModelST)
  Fixed <- reactive(input$fixed)
  
  fixedData_List <- reactive({ 
    
    if(!is.null(Fixed()) & Fixed()!="NULL"){
      getFixedData_List(predictionData(),Trait(),Fixed(),TargetIDs())
    }else if(is.null(Fixed()) | Fixed() == "NULL"){ 
       "NULL" 
    } 
  })
      
  
 
  outputList <- eventReactive(input$RunPredictionsST,{ 
    if(nSelTraits()==1){
      outputDF <-  withProgress(message = 'Running Computations', value = 0, {
        getRankedPredictedValues(predictionData(),nTraits(),unlist(Trait()),GPModelST(),fixedX=Fixed(),fixedData=fixedData_List(),optTS())})
      return(outputDF)
    }else if(nSelTraits()>1){
      outputDFComb <- c()
      for(nSelT in 1:nSelTraits()){
        i <- reactive(nSelT)
        outputDF <-  withProgress(message = 'Running Computations', value = 0, {
        getRankedPredictedValues(predictionData(),nTraits(),unlist(Trait())[i()],GPModelST(),fixedX=Fixed(),fixedData=fixedData_List(),optTS())})
       
       outputDFComb <- cbind(outputDFComb,outputDF[,2])
      } 
      outputDFComb <- cbind(outputDF[,1],outputDFComb,outputDF[,3])
      colnames(outputDFComb) <- c("GermplasmID",unlist(Trait()),"UpperBound of Reliability")
      return(outputDFComb)
    }
    
 })
  
  
### GPModel MT
  
  GPModelMT <- reactive(input$GPModelMT)

## 
  
  
  outputListMT <- eventReactive(input$RunPredictionsMT,{ 
    
    withProgress(message = 'Running Computations', value = 0, {
      getRankedPredictedValuesMT(predictionData(),nTraits(),unlist(Trait()),GPModelMT(),optTS())
    })
  })

### GPModel ME

 GPModelME <- reactive(input$GPModelME)
  
  outputListME <- eventReactive(input$RunPredictionsME,{ 
    
    withProgress(message = 'Running Computations', value = 0, {
      getRankedPredictedValuesME(predictionData(),nTraits(),unlist(Trait()),GPModelMT(),optTS())
    })

 })
  
 
  
## Render Table for ST  
 
  output$Ranked_Lines_for_SelectionST <- renderTable({
     
     
      as.data.frame(outputList()[1:20,])
     
  })

## Render Table for MT
  
  output$Ranked_Lines_for_SelectionMT <- renderTable({ if(nSelTraits()>1){
    
      (outputListMT()[1:20,])
   }
  }) 
  #%>% bindCache(outputListMT(),processedData(),nTraits(),unlist(Trait()),GPModelMT()) %>% bindEvent()
  
  
  output$scatter <- renderUI({ 
    if(nSelTraits()==1){
      plotOutput("plots")
    }
  }) 
  output$scatterMT <- renderUI({  
    if(nSelTraits()>1){
      plotOut_List <- lapply(1:nSelTraits(),function(x) { 
      i <- reactive(x)
      plotOutput(noquote(paste("'","plot",i(),"'",sep="")))})
      do.call(tagList,plotOut_List)
    }
  })
  output$plots <- renderPlot({
      if(nSelTraits()==1){
            plot.default(outputList()[,2],outputList()[,3],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main=paste("Upper bound of Reliability vs Predicted Values for ",Trait(),sep=""),font.lab=2,cex.lab=1.25)
      } 
  })
    
  
    
  observeEvent(input$RunPredictionsST,{
    if(nSelTraits()>1){
      lapply(1:nSelTraits(),function(i){
        loc_i <- reactive(i)
        output[[noquote(paste("'","plot",loc_i(),"'",sep=""))]] <- renderPlot({
          plot.default(outputList()[,(loc_i()+1)],outputList()[,(nSelTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main=paste("Upper Bound of Reliability vs Predicted Values for ",Trait()[loc_i()],sep=""),font.lab=2,cex.lab=1.25)
        })
        
      })
    }
  })  
  
############
  
    observeEvent(input$RunPredictionsMT,{ 
    
        lapply(1:nSelTraits(),function(i){
               loc_i <- reactive(i)
               output[[noquote(paste("'","plot",loc_i(),"'",sep=""))]] <- renderPlot({
                   plot.default(outputListMT()[,(loc_i()+1)],outputListMT()[,(nSelTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main=paste("Upper bound of Reliability vs Predicted Values for ",Trait()[loc_i()],sep=""),font.lab=2,cex.lab=1.25)
                })
         
        })
    
    }) 
       
  
#######
    
  switch_page <- function(i) {
    sel <- reactive(i)
    PageID <- c("Home","Load Data","Set Target & Trait","Optimize Training Population",
                "Cross Validations","Genomic Prediction","Visualize Predictions")
    updateTabsetPanel(inputId = "inData", selected = PageID[sel()])
  }
  
  
########
  output$ExportImpGenoDF <- downloadHandler(
    filename = function() {
      "ImputedGenotypesDF.txt"
    },
    content = function(file) {
      
      write.table(as.data.frame(GenoImp_DF()), file, row.names = FALSE)
    }
  )
  
### Export merged training geno-pheno
  
  output$ExportMerged_CompTS_DF <- downloadHandler(
    filename = function() {
      "Complete_Merged_GenoPheno_TSDF.txt"
    },
    content = function(file) {
       write.table(as.data.frame(processedData()[[1]]), file, row.names = FALSE)
    }
  )
  
### Export target geno set  
  output$ExportTargetGeno_DF <- downloadHandler(
    filename= function(){
      "Filtered_Target_GenoData.txt"
    },
    content= function(file){
      write.table(as.data.frame(processedData()[[2]]),file,row.names=FALSE)
    } 
   )
  
### Export Optimal TS
  
  output$ExportOptTS <- downloadHandler(
    filename = function() {
      "OptimalTS.csv"
    },
    content = function(file) {
      
      write.csv(as.data.frame((Train_STPGA()[[2]])), file, row.names = FALSE)
    }
  )
  

  output$ExportOut <- downloadHandler(
      filename = function() {
        "PredictionsST.csv"
      },
      content = function(file) {
         
        write.csv(as.data.frame(outputList()), file, row.names = FALSE)
      }
  )
  
  output$ExportOutMT <- downloadHandler(
    filename = function() {
      "PredictionsMT.csv"
    },
    content = function(file) {
      
        write.csv(as.data.frame(outputListMT()), file, row.names = FALSE)
    }
  )
  
  
  # session$allowReconnect("force")
  # options(shiny.launch.browser=FALSE)

}

shinyApp(server=server,ui=ui)

 



## 
# observeEvent(input$infileBLUEs, {
#   updateTabsetPanel(session,"inData",
#                     output$PhenoTable <- renderDataTable({as.data.frame(Pheno()[1:5,])})) }
# )
# 



# observeEvent(input$infileBLUEs,{cat(paste("Phenotype Table with ",nrow(Pheno())," lines and ",ncol(Pheno())-1," traits\n",sep=""))})



#Meta 
# MetaTab <- reactive({
# 

#metaFile <- input$infileMetaTable
#   
#   ext <- tools::file_ext(metaFile$datapath)
#   req(metaFile)
#   validate(need(ext == "csv", "Please upload a csv file"))
#   #NUST_Meta_Table <- 
#   read.csv(metaFile$datapath)
#   
# })


#fluidRow(
#column(3),column(width=9,tags$h4(tags$strong(textOutput("boxHeader"))))),
#tags$br(),
#plotOutput("plotBox")

#  
# withProgress(message = 'Processing Data', value = 0, {
#  NUST_Data_Table_Num_Filt <- getProcessedData(Geno(),Pheno(),MetaTab())
#     
#     
#   })

# fileInput("infileMetaTable", "Choose Meta Information File (CSV)", accept = ".csv"),
# checkboxInput("header", "Header", TRUE)

# metaHead <- eventReactive(input$infileMetaTable,{ paste("Table with meta information on ",nrow(MetaTab())," lines",sep="")})
# output$MetaHeader <- renderText({metaHead()})
# output$MetaTable <-  renderTable({as.data.frame(MetaTab()[1:5,1:5])})
# 
# tags$h3(textOutput("MetaHeader")),
# tableOutput("MetaTable"),

# output$Ranked_Lines_for_Selection <- renderDataTable({cbind.data.frame(rownames(outputList())[1:50],outputList()[1:50,])})
# output$scatter <- renderPlot({plot.default(outputList()[,1],outputList()[,2],type="p",xlab="Observed Value",ylab='Predicted Value',main="Predicted vs Observed Values")})


#### Convert VCF to Dataframe
######################################################################
###                         VCFToSimple                            ###
###   a function to generate simple format from vcf                ###
######################################################################
# function to convert vcf into simple data format
# need vcfR library and dplyr library
# input: vcf file path and name



# fluidRow(
#    column(3),column(width=5,tags$li(" A '.csv' file containing information on the target data, whose genetic values are to be predicted using the genomic prediction model. 
#                        The first column in the target data information file should have the line IDs."))),
#      
#  fluidRow(
#    column(3),column(width=5,fileInput("infileTargetTable", "Choose Target Data Information File (CSV)", accept = ".csv"),
#      checkboxInput("header", "Header", TRUE)
#    )),
#   fluidRow(
#        column(3),column(width=5,div(
#          tags$h4(tags$strong(textOutput("TargetHeader"))),
#          tags$h6(tableOutput("TargetTable")))))


#  fluidRow(
#   column(6),column(width=5,tags$h6(tableOutput("PhenoTable")))),
#  tags$br(),




# 
# fluidRow(column(2),column(width=5,checkboxInput("header", "Header", TRUE)),
#   column(6),column(width=5,checkboxInput("header", "Header", TRUE))),
# 
# fluidRow(
#   column(6),column(width=5,div(
#     tags$h4(tags$strong(textOutput("GenoHeader"))),
#     tags$h6(tableOutput("GenoTable"))))),
#  tags$br(),
#  tags$br(),
#  tags$br(),
#  fluidRow(
#  column(6),column(width=5,tags$ul(tags$li("A '.csv' file containing phenotypic information. Phenotype file should have the line ID in the first column and phenotypic trait values after that with the trait ID as column name.")))),
#  tags$br(),
#  fluidRow(
#  column(6),column(width=5,fileInput("infileBLUEs", "Choose Phenotype File (CSV)", accept = ".csv"),
#    checkboxInput("header", "Header", TRUE))),
#      
#  tags$br(),
#  fluidRow(
#   column(6),column(width=5,tags$h4(tags$strong(textOutput("PhenoHeader"))))),
#  fluidRow(
#   column(6),column(width=5,tags$h6(tableOutput("PhenoTable")))),
#  tags$br(),
#  
#  fluidRow(
#    column(3),column(width=5,tags$li(" A '.csv' file containing information on the target data, whose genetic values are to be predicted using the genomic prediction model. 
#                        The first column in the target data information file should have the line IDs."))),
#      
#  fluidRow(
#    column(3),column(width=5,fileInput("infileTargetTable", "Choose Target Data Information File (CSV)", accept = ".csv"),
#      checkboxInput("header", "Header", TRUE)
#    )),
#   fluidRow(
#        column(3),column(width=5,div(
#          tags$h4(tags$strong(textOutput("TargetHeader"))),
#          tags$h6(tableOutput("TargetTable")))))


# if(nSelTraits()>1){
#   PATable <- cvrOutputListMT()
#   colnames(PATable) <- rep("",ncol(PATable))
#   print.data.frame(as.data.frame(PATable))
# } 
# 
#})
# if(nSelTraits()==1){
#   cvrOutputListST <- reactive({
#     withProgress(message = 'Running CrossValidations', value = 0, {
#       getemCVR(processedData(),unlist(Trait()),nTraits(),k(),nIter())})
#     
#   }) 
#   PATable <- rbind(c("RR","BB","BL"),round(cvrOutputListST()[c("emRR","emBB","emBL")],digits=2)) 
#   rownames(PATable)<- c("Prediction Model","Prediction Accuracy")
#   colnames(PATable) <- rep("",ncol(PATable))
#   PATable
#   
#   
# }
# }) # %>% bindEvent(input$CrossValidationST)

#%>% bindCache(processedData(),unlist(Trait()),nTraits(),k(),nIter()) 


#%>% bindCache(cvrOutputListMT(),processedData(),Trait(),nTraits(),k(),nIter(),nSelTraits())

## One plot in renderUI
# if(nSelTraits()==1){
#   plotOutput("plots")
# } 
## Multiple plots in renderUI
# if(nSelTraits()>1){
#   plotOutput("plots1")
#   plotOutput("plots2")
#   # plotOut_List <- lapply(1:nSelTraits(),function(x) {plotname <- paste("plot",x,sep="") 
#   # plotOutput(plotname)})
#   # do.call(tagList,plotOut_List)
# }



## Multiple plots in renderPlot      

# for(i in 1:nSelTraits()){ 
#   local({
#     loc_i <- i
#     plotName <- paste("plot",loc_i,sep="")
#     output[[plotName]] <- renderPlot({
#         plot.default(outputListMT()[,2],outputListMT()[,(nSelTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")
#     })
#   })
# }

#    reactive(
# ### One plot in renderPlot
#      if(nSelTraits()==1){
#         output$plots <- renderPlot({
#          plot.default(outputList()[,2],outputList()[,3],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")
#         })
#      }
#    )



# output$Ranked_Lines_for_Selection <- renderDataTable({as.data.frame(outputListMT()[1:50,])})
# output$scatter <- renderPlot({plot.default(outputListMT()[,2],outputListMT()[,(nTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")})
# 


#   outputListMT <- reactive({
# 
#     withProgress(message = 'Running Computations', value = 0, {
#       getRankedPredictedValuesMT(processedData(),nTraits(),unlist(Trait()),GPModelMT(),optTS=NULL)
#     })
#    })  %>% bindCache(processedData(),Trait(),nTraits(),GPModelMT()) %>% bindEvent(input$RunPredictionsMT)



# output$Ranked_Lines_for_Selection <- renderDataTable({as.data.frame(outputList()[1:50,])})
# output$scatter <- renderPlot({ plot.default(outputList()[,2],outputList()[,3],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")})




#%>% bindCache(TSOptOutputList(),processedData(),unlist(Trait()),nTraits(),noCandidates(),nTrainToSelect()) %>% bindEvent(input$noToSelect,input$Optimize,input$trait,input$infileVCF,input$infileTargetTable,input$infileBLUEs)

# TSTable_Mod <- cbind(c("STPGA","Random"),apply(TSTable,2,function(x) round(as.numeric(x),digits=2)))
#colnames(TSTable_Mod) <- c("TS","RR","BB","BL")
# setGAPar <- reactive(input$setGA)


# output$TSGA <- renderUI({ 
#   if(setGAPar() ==TRUE){} 
#   
#   if(setGAPar() ==FALSE){
#   
#     numericInput(inputId="noPop","GA Population Size",value = 100,min =1,max=1000)
#   
#     numericInput(inputId="noElite","Elite Population Size",value =10,min =1, max=50)
#     numericInput(inputId="mutProb","Mutation Probability",value =0.5,min =0.1, max=1)
#     
#     numericInput(inputId="mutInt","Mutation Intensity",value =1,min =0, max=1)
#    
#     numericInput(inputId="nIter","# Iterations",value =100,min =50, max=1000)
#  
#     selectInput(inputId="tabu","Tabu Search",choices=c("TRUE","FALSE"),selected="TRUE")
#     
#     numericInput(inputId="tabumemsize","Tabu Memory Size",value =1,min =1, max=10)
#    
#     numericInput(inputId="lambda","Shrinkage Parameter (lambda)",value =1e-6,min =1e-9, max=1)
#     
#     numericInput(inputId="mcCores","#Cores",value =10,min =1, max=20)
#    
#     
#   }
#   
# })
# 
# Mssg1CVR <- eventReactive(input$CrossValidationST,{
# if(nSelTraits>1){
#    paste("Performing CrossValidation for the first trait")
#   
#  }
# }) 

# cvrOutputListST <- eventReactive(input$CrossValidationST,{
# 
#       withProgress(message = 'Running CrossValidations', value = 0, {
#       getemCVR(processedData(),unlist(Trait()),nTraits(),k(),nIter())})
# 
# })
# 

# %>% bindEvent(input$CrossValidationMT)

#nSelTraits <- eventReactive(input$CrossValidationST,{(1)})


# observeEvent(input$Home_Data,{ switch_page(2)})
# observeEvent(input$Data_Home, switch_page(1))
# 
# observeEvent(input$Data_Trait, switch_page(3))
# observeEvent(input$Trait_Data, switch_page(2))
# 
# observeEvent(input$Trait_TS, switch_page(4))
# observeEvent(input$TS_Trait, switch_page(3))
# 
# observeEvent(input$TS_CVR, switch_page(5))
# observeEvent(input$CVR_TS, switch_page(4))
# 
# observeEvent(input$CVR_GP, switch_page(6))
# observeEvent(input$GP_CVR,switch_page(5))
# 
# observeEvent(input$GP_VP, switch_page(7))
# observeEvent(input$VP_GP, switch_page(6))

# output$ExportOut <- downloadHandler(
#   filename = function() {
#     "Predictions.csv"
#   },
#   content = function(file) {
#     if(nSelTraits()==1){
#       write.csv(as.data.frame(outputList()), file, row.names = FALSE)
#     }else if(nSelTraits()>1){
#       write.csv(as.data.frame(outputListMT()), file, row.names = FALSE)
#     }
#   }
# )

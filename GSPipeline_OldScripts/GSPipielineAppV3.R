
library(shiny) 


## Set the file upload limit to 40 MB

options(shiny.maxRequestSize=40*1024^2)
source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/GS_Pipeline_Jan_2022_FnsApp.R")

ui <- fluidPage( 
  fluidRow(
    column(4),
    column(6, div(
      id = "app-title",
      titlePanel(tags$strong("SOYGEN2 Genomic Selection App")),
      tags$p("An app to perform genomic predictions given genotypic and phenotypic data")
    ))),
  
  tabsetPanel(id="inData",
              
              ## Tab for loading data
              
              tabPanel("Load Data",
                       sidebarLayout(
                         
                         sidebarPanel(
                           fileInput("infileBLUEs", "Choose Phenotype File (CSV)", accept = ".csv"),
                           checkboxInput("header", "Header", TRUE),
                           
                           fileInput("infileVCF", "Choose Genotype File (VCF)", accept = ".vcf"),
                           checkboxInput("header", "Header", TRUE),
                           
                           fileInput("infileTargetTable", "Choose Target Data Information File (CSV)", accept = ".csv"),
                           checkboxInput("header", "Header", TRUE)
                           
                         ),
                         
                         mainPanel(
                           fluidRow(
                             column(3),column(width=4,tags$h3(tags$strong("Load Data Files")))),
                           fluidRow(
                             column(width=11,tags$h4("Load data files containing phenotypic data, genotypic data and information on target lines."),
                                    tags$br(),
                                    tags$ul( tags$li("A '.csv' file containing phenotypic information. Phenotype file should have the line ID in the first column and phenotypic trait values after that with the trait ID as column name."),
                                             tags$br(),
                                             tags$li("A '.vcf' containing the genotypic information of both the candidate lines used to train the GP model and the target lines whose values need to be predicted."), 
                                             tags$br(),
                                             tags$li(" A '.csv' file containing information on the target data, whose genetic values are to be predicted using the genomic prediction model. 
                     The first column in the target data information file should have the line IDs.")
                                    ))),
                           tags$br(),
                           
                           fluidRow(
                             column(2),column(width=9,tags$h4(tags$strong(textOutput("PhenoHeader"))))),
                           fluidRow(
                             column(2),column(width=9,tags$h6(tableOutput("PhenoTable")))),
                           tags$br(),
                           fluidRow(
                             column(2),column(width=9,div(
                               tags$h4(tags$strong(textOutput("GenoHeader"))),
                               tags$h6(tableOutput("GenoTable"))))),
                           tags$br(),
                           fluidRow(
                             column(2),column(width=9,div(
                               tags$h4(tags$strong(textOutput("TargetHeader"))),
                               tags$h6(tableOutput("TargetTable")))))
                         ))
              ),
              
              ## Tab for Trait Selection 
              
              tabPanel("Choose Trait",
                       
                       sidebarLayout(
                         sidebarPanel(
                           
                           selectInput(inputId="trait","Choose Trait",choices=NULL)
                           
                         ),
                         mainPanel(
                           fluidRow(
                             column(3),column(width=4,tags$h3(tags$strong("Select Trait")))),
                           fluidRow(
                             column(1), column(width=11,tags$h5("Select one trait to train the genomic prediction model."))),
                           tags$br(),
                           tags$br(),
                           tags$br(),
                           tags$br(),
                           fluidRow(
                             column(1),column(width=11,tags$h4(tags$strong(textOutput("summaryHeader"))))),
                           tags$br(),
                           fluidRow(
                             column(2),column(width=10,tableOutput("Summary"))),
                           tags$br(),
                           tags$br(),
                           tags$br()
                           
                         ))
              ),
              
              #   ## TS optimization tab
              tabPanel("Optimize Training Population",
                       sidebarLayout(
                         sidebarPanel(
                           
                           numericInput(inputId="noToSelect","Training Set Size to Select",value =200,min =100, max=200),
                           tags$br(),
                           actionButton(inputId="Optimize","Optimize Training Sets"),
                           tags$br(),
                           tags$br(),
                           
                           #fileInput("infileOptTS", "Choose Optimized Training Set  File (csv)", accept = ".csv"),
                           #checkboxInput("header", "Header", TRUE)
                           
                           
                         ), mainPanel(
                           fluidRow(
                             column(1),column(width=11,tags$h3(tags$strong("Select Optimal Training Population using 'STPGA' R Pkg")))),
                           fluidRow(
                             column(1),
                             column(width=10, tags$p("A training population for genomic prediction is selected by optimizing an objective function using a genetic algorithm.
                              The optimization method uses the 'GenAlgForSubsetSelection' function implemented in the 'STPGA' R package (Akdemir 2017).The user needs to enter the number of lines to select for the training population.
                              The parameters for the genetic algorithm are set to default values"),
                                    tags$br(),
                                    tags$br()
                             )),
                           
                           fluidRow(
                             column(2),
                             column(6, div(
                               tags$h5(tags$strong(textOutput("tsOptHeader"))),
                               tags$h6(tableOutput("PredAccuracyTable")), 
                               tags$br()
                             )))
                         ))
              ), 
              
              #   
              ### GP Tab 
              tabPanel("Genomic Prediction",
                       sidebarLayout(
                         sidebarPanel(
                           numericInput(inputId="k",label = "Enter k for k-fold cross validation",value=5,min=2,max=10),
                           numericInput(inputId="nIter",label = "Enter n for n-iterations of each cycle of cross validation",value=5,min=2,max=10),
                           actionButton("CrossValidation", "Run Cross Validation (Optional)"),
                           
                           tags$br(),  
                           tags$br(),
                           tags$br(),  
                           
                           selectInput(inputId="GPModel","Choose Prediction Model",c("rrBLUP (rrBLUP)","rrBLUP (bWGR)","BayesB (bWGR)","BayesLASSO (bWGR)")),
                           actionButton("RunPredictions", "Predict!"),
                           tags$br()
                           
                         ), mainPanel(
                           fluidRow(
                             column(3),column(width=8,tags$h3(tags$strong("Train Genomic Prediction Model")))
                           ),
                           fluidRow(
                             column(1),
                             column(width=10, tags$p("Perform an optional k-fold cross validation with training dataset to identify the model with the best prediction accuracy.
                                           K-fold cross validation is performed using 'emCV' function implemented in the 'bWGR' package (Xavier et al. 2020).
                                           You can also skip the cross validation and select the method to train the genomic prediction model and predict the values of target lines. rrBLUP method is implemented using the
                                           rrBLUP (Endelman 2011) package. Expectation maximization based RR-BLUP, BayesB and BayesLASSO methods are implemented using the bWGR package (Xavier et al. 2020)."))),
                           tags$br(),
                           fluidRow(
                             column(2),column(width=8,tags$h5(tags$strong(textOutput("cvrHeader"))))) ,
                           tags$br(),
                           fluidRow(
                             column(3),column(width=8,tags$h5(tableOutput("emCVR"))))
                         )
                       )
              ),
              #   
              #   ## Tab for Results
              tabPanel("Visualize Predictions",
                       sidebarLayout(
                         sidebarPanel(
                           
                           downloadButton("ExportOut", "Export Output Table"),
                           
                         ),
                         mainPanel(
                           
                           fluidRow(
                             column(3),column(width=8,tags$h3(tags$strong("Visualize and Explore Predictions")))
                           ),
                           tags$br(),
                           tags$h4(textOutput("RankedLinesHeader")),
                           dataTableOutput("Ranked_Lines_for_Selection"),
                           #tags$h4("Predicted vs Observed Values"),
                           plotOutput("scatter") 
                           
                         )   
                       ) 
              )
  )
)


#))

server <- function(input,output){
  
  #Geno
  
  Geno <- reactive({
    
    genoFile <- input$infileVCF
    
    
    ext <- tools::file_ext(genoFile$datapath)
    req(genoFile)
    validate(need(ext == "vcf", "Please upload a vcf file"))
    
    withProgress(message = 'Reading Data', value = 0, {
      NUST_Genotypes_VCF <- read.table(genoFile$datapath)})
    
    withProgress(message = 'Converting VCF to Dataframe', value = 0, {
      gt2d <- VCFtoDF(genoFile$datapath) }) 
    
    gt2d
  })
  
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
  
  
  #Target
  TargetTab <- reactive({
    TargetFile <- input$infileTargetTable
    
    ext <- tools::file_ext(TargetFile$datapath)
    req(TargetFile)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(TargetFile$datapath, header = input$header)
    
  })
  
  
  
  phenoHead <- eventReactive(input$infileBLUEs,{ paste("Phenotype Table with ",nrow(Pheno())," lines and ",ncol(Pheno())-1," traits",sep="")})
  output$PhenoHeader <- renderText({phenoHead()})
  output$PhenoTable <- renderTable({as.data.frame((Pheno())[1:5,1:5])})
  
  genoHead <- eventReactive(input$infileVCF,{ paste("Genotype Table with ",nrow(t(Geno()))," lines and ",ncol(t(Geno()))," markers",sep="")})
  output$GenoHeader <- renderText({genoHead()})
  output$GenoTable <- renderTable({as.data.frame(t(Geno())[1:5,1:5])})
  
  
  TargetHead <- eventReactive(input$infileTargetTable,{ paste("Table with information on ",nrow(TargetTab())," Target lines",sep="")})
  output$TargetHeader <- renderText({TargetHead()})
  output$TargetTable <-  renderTable({as.data.frame(TargetTab()[1:5,c("GermplasmId","Parentage","Generation","Reason","Source")])})
  
  # 
  # 
  
  nTraits <- eventReactive(input$infileBLUEs,{(ncol(Pheno())-1)})
  TargetIDs <- eventReactive(input$infileTargetTable,{as.character(TargetTab()[,1])})
  
  
  
  ## Trait 
  
  
  Trait <- reactive({input$trait})
  
  summaryHead <- eventReactive(input$trait,{paste("Summary of distribution of selected phenotypic values")})
  output$summaryHeader <- renderText({summaryHead()})
  
  SummaryTxt <-  eventReactive(input$trait,{summary(Pheno()[,Trait()])})
  output$Summary <- renderTable({trTable <- rbind(as.vector(names(SummaryTxt())),as.vector(round(SummaryTxt(),digits=2)))
  colnames(trTable)<- rep("",length(as.vector(SummaryTxt())))
  trTable
  })
  
  boxHead <- eventReactive(input$trait,{paste("Boxplot of selected phenotypic values")})
  output$boxHeader <- renderText({boxHead()})
  
  phenoValues <- eventReactive(input$trait,{Pheno()[,Trait()]})
  output$plotBox <- renderPlot({boxplot(phenoValues(),xlab=Trait())})
  
  
  processedData <-  eventReactive(input$trait,{
    withProgress(message = 'Merging Pheno and Geno Data', value = 0, {
      getProcessedData(Geno(),Pheno(),TargetIDs(),Trait())
    }) 
  }) 
  # }
  
  ## TS Optimization
  
  
  nTrainToSelect <- reactive(input$noToSelect)
  # inOptTS <- reactive(input$infileOptTS)
  # 
  # 
  #  optTS <- reactive({
  #   tsFile <- input$infileOptTS
  #   
  #   ext <- tools::file_ext(tsFile$datapath)
  #   req(tsFile)
  #   validate(need(ext == "csv", "Please upload a csv file"))
  #   read.csv(tsFile$datapath)
  #   
  # })
  #  optTS <- reactive({ 
  #  if(inOptTS() =="No file selected"){
  #     NULL
  #  }
  #  })
  # 
  # 
  noToReduce <- reactive(500)
  TSOptOutputList <- eventReactive(input$Optimize,{
    withProgress(message = 'Running Optimizations', value = 0, {
      Train_STPGA <- reactive(getOptimalTS(processedData(),Trait(),nTraits(),noToReduce(),nTrainToSelect()))
      
      Train_Random <- reactive(getRandomTS(processedData(),Trait(),nTraits(),noToReduce(),nTrainToSelect()))
      getTSComparisons(processedData(),Train_STPGA(),Train_Random(),Trait(),nTraits(),NULL)
      
    })
  })
  
  tsOptHead <- eventReactive(input$Optimize,{paste("Correlation between observed and predicted values")})
  output$tsOptHeader <- renderText({tsOptHead()}) 
  output$PredAccuracyTable <- renderTable({
    TSTable <- TSOptOutputList()
    TSTable_Mod <- cbind(c("STPGA","Random"),TSTable)
    colnames(TSTable_Mod) <- c("TS","RR","BB","BL")
    TSTable_Mod})
  
  ## CVR  
  
  k <- reactive(input$k)
  nIter <- reactive(input$nIter)
  
  cvrOutputList <- eventReactive(input$CrossValidation,{  
    
    withProgress(message = 'Running CrossValidations', value = 0, {
      getemCVR(processedData(),Trait(),nTraits(),k(),nIter()) 
    })
    
  })
  
  cvrHead <- eventReactive(input$CrossValidation,{paste("Correlation between observed and predicted values for ",Trait(),sep="")})
  output$cvrHeader <- renderText({cvrHead()})  
  output$emCVR <- renderTable({
    PATable <- rbind(c("RR","BB","BL"),round(cvrOutputList()[c("emRR","emBB","emBL")],digits=2)) 
    rownames(PATable)<- c("Prediction Model","Prediction Accuracy")
    colnames(PATable) <- rep("",3)
    PATable
  })
  
  
  ## GP  
  GPModel <- reactive(input$GPModel)
  outputList <- eventReactive(input$RunPredictions,{ 
    
    withProgress(message = 'Running Computations', value = 0, {
      getRankedPredictedValues(processedData(),nTraits(),Trait(),GPModel(),optTS=NULL)
    })
    
  })
  
  output$Ranked_Lines_for_Selection <- renderDataTable({as.data.frame(outputList()[1:50,])})
  output$scatter <- renderPlot({plot.default(outputList()[,2],outputList()[,3],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")})
  
  output$ExportOut <- downloadHandler(
      filename = function() {
        "Predictions.csv"
      },
      content = function(file) {
        write.csv(as.data.frame(outputList()), file, row.names = FALSE)
      }
    )

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

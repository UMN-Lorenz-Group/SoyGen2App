library(shiny) 


## Set the file upload limit to 40 MB

options(shiny.maxRequestSize=40*1024^2)
source("C:/Users/ivanv/Desktop/UMN_Projects/IMS_Soy/GS_Pipeline_rrBLUP_Jan_2022_FnsApp.R")

ui <- fluidPage( 
  fluidRow(
    column(5),
    column(7, div(
      id = "app-title",
      titlePanel(tags$strong("Soygen2 Genomic Selection App")),
      tags$p("An app to perform genomic predictions given genotypic and phenotypic data")
    ))),
  
 
  
  sidebarLayout(
    sidebarPanel(
      
      fileInput("infileBLUEs", "Choose Phenotype File (CSV)", accept = ".csv"),
      checkboxInput("header", "Header", TRUE),
      
      fileInput("infileVCF", "Choose Genotype File (VCF)", accept = ".vcf"),
      checkboxInput("header", "Header", TRUE),
      
      fileInput("infileMetaTable", "Choose Meta Information File (CSV)", accept = ".csv"),
      checkboxInput("header", "Header", TRUE),
      
      selectInput(inputId="Trait","Choose Trait",character(0)),
    
      selectInput(inputId="GPModel","Choose Prediction Model",c("rrBLUP","BayesB","BayesLASSO")),
      
      actionButton("RunCVR", "Run Optional Cross Validations"),
      
      actionButton("RunPredictions", "Predict!")
      
      # fileInput("infileLinkMap", "Choose Linkage Map File (CSV)", accept = ".csv"),
      # checkboxInput("header", "Header", TRUE)
      
    ),
    mainPanel(
      
      
       # tags$h4("Phenotype Table"),
      # textOutput("PhenoHeader"),
       dataTableOutput("PhenoTable"),
       #tags$h4("Genotype Table"),
       #textOutput("Phenotype Table"),
       dataTableOutput("GenoTable"),
      
        #tags$h4("MetaInformation Table"),
       #textOutput("Phenotype Table"),
       dataTableOutput("MetaTable"),
     
      #tags$h4("Predicted vs Observed Values"),
      plotOutput("scatter"),
      #tags$h4("Ranked Lines for Selection"),
      dataTableOutput("Ranked_Lines_for_Selection")
      ##tableOutput("contents"))
    )
  )
)

server <- function(input,output){
     
  #Geno
  
  Geno <- reactive({
   
    genoFile <- input$infileVCF
    
   
    ext <- tools::file_ext(genoFile$datapath)
    req(genoFile)
    validate(need(ext == "vcf", "Please upload a vcf file"))
    NUST_Genotypes_VCF <- read.table(genoFile$datapath) 
    fluidRow(
      column(7,
             withProgress(message = 'Preparing Data', value = 0, {
               gt2d <- VCFtoDF(genoFile$datapath) }) ))
    
    gt2d
  })
  
  #Pheno
  Pheno <- reactive({
    
    phenoFile <- input$infileBLUEs
  
    ext <- tools::file_ext(phenoFile$datapath)
    req(phenoFile)
    validate(need(ext == "csv", "Please upload a csv file"))
    #NUST_BLUEs <- 
    read.csv(phenoFile$datapath, header = input$header)
    
  })
  
  observeEvent(input$infileBLUEs, {
    updateSelectInput(inputId = "Trait", choices = colnames(Pheno())[2:ncol(Pheno())])})
  
  #Meta 
  MetaTab <- reactive({
     metaFile <- input$infileMetaTable
    
    ext <- tools::file_ext(metaFile$datapath)
    req(metaFile)
    validate(need(ext == "csv", "Please upload a csv file"))
    #NUST_Meta_Table <- 
    read.csv(metaFile$datapath)
   
  })
    # inputTables$Pheno <- NUST_BLUEs
    # inputTables$Geno <- t(gt2d)
    # inputTables$Meta <-  NUST_Meta_Table
    # inputTables$VCFName <- genoFile$datapath
    # 
  #output$PhenoHeader <- eventReactive(input$infileBLUEs,{renderText({paste("Phenotype Table with ",nrow(Pheno())," lines and ",ncol(Pheno())-1," traits",sep="")})})
   output$PhenoTable <- renderDataTable({as.data.frame(Pheno()[1:5,])}) 
   output$GenoTable <- renderDataTable({as.data.frame(t(Geno())[1:5,1:5])})
   output$MetaTable <-  renderDataTable({as.data.frame(MetaTab()[1:5,1:5])})
  
    # output$GenoTable <- eventReactive(input$infileVCF,{renderDataTable({as.data.frame(t(Geno())[1:10,1:10])}) })
    # output$PhenoTable <- eventReactive(input$infileBLUEs,{renderDataTable({as.data.frame(Pheno()[1:10,1:10])}) })
    # output$MetaTable <- eventReactive(input$infileMetaTable,{renderDataTable({as.data.frame(MetaTab()[1:10,1:10])})})
    # 
   outputList <- eventReactive(input$RunPredictions,{ 
      
      gt2d <- Geno()
      # print(Geno())
      # print(Pheno())
      # print(Meta())
      #Msssg <- paste("Processing Data with ",ncol(gt2d)," genotypes and ",nrow(gt2d)," markers",sep = "")
     
     
      
      withProgress(message = 'Processing Data', value = 0, {
        NUST_Data_Table_Num_Filt <- getProcessedData(Geno(),Pheno(),MetaTab())
        
        
      })
      
      withProgress(message = 'Running Computations', value = 0, {
        getRankedPredictedValues(NUST_Data_Table_Num_Filt) 
      })
    
   })
  
  output$Ranked_Lines_for_Selection <- renderDataTable({cbind.data.frame(rownames(outputList())[1:50],outputList()[1:50,])})
  output$scatter <- renderPlot({plot.default(outputList()[,1],outputList()[,2],type="p",xlab="Observed Value",ylab='Predicted Value',main="Predicted vs Observed Values")})
}

shinyApp(server=server,ui=ui)



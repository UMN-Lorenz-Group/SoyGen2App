#### Improve alignment of file input elements in the app 
## ChatGPT 
if(!require(shiny)){
  install.packages("shiny")
}
library(shiny) 

### Better ways to organize the output tables in UI elements 
### Change source file path to working directory

FN <- paste(getwd(),"/GS_Pipeline_Jan_2022_FnsApp.R",sep="")
source(FN)

PN <-  paste(getwd(),"/GSPipeline.png",sep="")

## Set the file upload limit to X MB
options(shiny.maxRequestSize=500*1024^2)

ui <- fluidPage( 
  #theme = bslib::bs_theme(bootswatch = "solar"),
  
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
           column(4, offset = 4,
                  div(
                    h3(strong("Load Data Files")),
                    br(),
                    fluidRow(
                      column(6,
                             div(
                               fileInput("infileBLUEs", "Choose Phenotype File (CSV)", accept = ".csv"),
                               checkboxInput("headerPheno", "Header", TRUE),
                               tags$ul(
                                 tags$li("A '.csv' file containing phenotypic information."),
                                 tags$li("Phenotype file should have the line ID in the first column and phenotypic trait values after that with the trait ID as column name.")
                               ),
                               h4(strong(textOutput("PhenoHeader"))),
                               h6(tableOutput("PhenoTable"))
                             )
                      ),
                      column(6,
                             div(
                               fileInput("infileVCF", "Choose Genotype File (VCF)", accept = ".vcf"),
                               checkboxInput("headerGeno", "Header", TRUE),
                               tags$ul(
                                 tags$li("A '.vcf' containing the genotypic information of both the candidate lines used to train the GP model and the target lines whose values need to be predicted.")
                               ),
                               h4(strong(textOutput("GenoHeader"))),
                               h6(tableOutput("GenoTable"))
                             )
                      )
                    )
                  )
           )
         )
)

)) 



server <- function(input,output){
  
  #Geno
  
  Geno <- reactive({
    
    genoFile <- input$infileVCF
    ext <- tools::file_ext(genoFile$datapath)
    req(genoFile)
    validate(need(ext == "vcf", "Please upload a vcf file"))
    
    # withProgress(message = 'Reading Data', value = 0, {
    #   NUST_Genotypes_VCF <- read.table(genoFile$datapath)})
    
    withProgress(message = 'Converting VCF to Dataframe', value = 0, {
      gt2d <- VCFtoDF(genoFile$datapath) }) 
    
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
  
  phenoHead <- eventReactive(input$infileBLUEs,{ paste("Phenotype Table with ",nrow(Pheno())," lines and ",ncol(Pheno())-1," traits",sep="")})
  output$PhenoHeader <- renderText({phenoHead()})
  output$PhenoTable <- renderTable({as.data.frame((Pheno())[1:5,1:5])})
  
  genoHead <- eventReactive(input$infileVCF,{ paste("Genotype Table with ",ncol(Geno())-5," lines and ",nrow((Geno()))," markers",sep="")})
  output$GenoHeader <- renderText({genoHead()})
  output$GenoTable <- renderTable({as.data.frame((Geno())[1:5,1:5])})
  
}

shinyApp(server=server,ui=ui)


#####



library(shiny)
library(DT)
library(dplyr)
library(pROC)
library(caret)
library(fastAdaboost)
library(kernlab)
library(randomForest)

app = function(model, algorithm, valDT) {

# valDT=read.table("/media/thong/sda/RCID/CoHoai/Tuberculosis/ShinyR/model/example.csv", header=T, sep=",")
rownames(valDT)=valDT$Sample
dt=valDT[,grep("cg", names(valDT))]

download.file(paste0("https://github.com/nthongle317/TB_pred/raw/main/model/", model, ".txt"), destfile=paste0(model, ".txt"))
CpG=read.table(paste0(model, ".txt"), header=T, sep="\t")$CpGs

REQ=data.frame(Required=CpG)
UI=data.frame(RequiredInput=names(dt), UserInput=names(dt))

check=merge(REQ, UI, by.x="Required", by.y="RequiredInput", all.x=T)
n=sum(is.na(check$UserInput))

if (n>0) {
check[is.na(check$UserInput), "UserInput"] = "Missing value in this CpG, check your file"
print(check)
} else if (n==0) {

select=paste0(model, "_", algorithm)

download.file(paste0("https://github.com/nthongle317/TB_pred/raw/main/model/", select, ".rds"), destfile=paste0(select, ".rds"))
rf.m1=readRDS(paste0(select, ".rds"))

pre.m1s <- predict(rf.m1, dt, type="prob")

pre.m1s$`TB possibility` <- round(pre.m1s$HC,2)
pre.m1s$`Prediction` <- ifelse(pre.m1s$HC>0.5, "TB", "Normal")

outDF=cbind(pre.m1s, dt)
outDF$Sample=rownames(outDF)

outDF=outDF[,c("Sample", "TB possibility", "Prediction", CpG)]
rownames(outDF)=NULL
print(outDF)

}

}

ui <- fluidPage(

   navbarPage("TB prediction",
           tabPanel("About", 
            mainPanel(
              h1("About the study"),
              h3("A machine learning approach utilizing DNA methylation as a classifier for Pulmonary Tuberculosis screening"),
              h6("Nhat-Thong Le, Hien Thi-Thu Do, Doan-Minh-Trung Duong, Hong-Ngoc-Doan Tran, Thuc-Quyen Huynh, Khon Huynh, Phuong-Thao Nguyen, Minh-Thong Le & Thi-Thu-Hoai Nguyen"), 
              br(),
              h4(strong("A simplified tools for TB infectious prediction")),
              p("Tuberculosis (TB) caused by Mycobacterium tuberculosis (Mtb) is one of the most lethal infectious diseases with estimates of approximately 1.6 million human deaths in 2022. Recently, non-sputum DNA methylation biomarkers emerged as a promising approach for rapid detection of TB infection. However, comprehensive work to explore potential of DNA methylation in TB prediction has been less common. Here, we aim to introduce a novel set of DNA methylation biomarkers based on characterizing blood DNA methylation. We conducted a pooled analysis of several datasets including 290 samples and used feature selection and five machine learning algorithms to identify the best CpG sites combination that could predict Mtb infectious. MUVR using random forest core algorithm (12-CpG model) combined with Random forest classifier showed outperformance to previous published methods with perfect performance in the discovery cohort and reached the sensitivity of 90%, the specificity of 82% and AUC of 0.91 (95% CI: 0.85-0.97) in the validation cohort. Further differential analysis of Bacillus Calmette-Guerin (BCG) vaccination groups and HIV - Mtb coinfection showed clear differences between BCG and non-BCG groups, as well as between HIV - Mtb coinfection, Mtb infection and healthy samples; which shown high potential to overcome traditional methods. Collectively, DNA methylation was a promising method for early detection of Tuberculosis and potentially a clinical tool for TB diagnostic biomarkers. Further external validation studies are needed to confirm the impact of our tool in daily practice."),
              br(),
              p("The detail report of our study can be found here: DOI/@@@"),
              br(),
              h1("About RCID"),
              h5(strong("Vision:")),
              p("In the short term, Research Center for Infectious Disease will become one of the leading research centers on infectious diseases in Vietnam, contributing to developing high-quality human resources for the locality and the whole country."),
              p("In the long term, the Center aims to develop into a Center of excellent research on Infectious Diseases in Southeast Asia, capable of creating and providing technological solutions and services in response to the epidemic diseases that may develop in our society."),
              br(),
              h5(strong("Mission:")),
              p("Improve the quality of basic research in parallel with applied research and technology transfer in order to meet the needs of society in responding to infectious diseases."),
              p("Enhancing the position and affirming the pioneering role of VNU-HCM in general and VNU in particular in research and community service."),
              p("Improve the quality of basic research in parallel with applied research and technology transfer in order to meet the needs of society in responding to infectious diseases."),
              p("Enhancing the position and affirming the pioneering role of VNU-HCM in general and VNU in particular in research and community service."),
              br(),
              h5(strong("Contact:")),
              p("O2.101A, International University – Vietnam National University HCM City, Quarter 6, Linh Trung Ward, Thu Duc City, Ho Chi Minh City, Vietnam."),
              p("Website: https://rcid.hcmiu.edu.vn/"),
              p("Tel: (+84)283 3724 4270 – (ext:3396)"),
              p("Email: rcid@hcmiu.edu.vn")

   )
            ),

           tabPanel("Guideline",
            mainPanel(
              h1("User guideline"),
              p("We provide a Shiny app tool equipped with a simple and user-friendly interface for screening Tuberculosis infectious condition in individuals. Only best classification model of Random forest using 12-CpG has been included in app, other classification model will be available upon request."),
              br(),
              h3("Data requirement"),
              p("Classification model is generated by using methylation beta value of 12 CpG sites including (cg00782174, cg01737507, cg06497752, cg10140638, cg10453758, cg11166252, cg14095850, cg15705999, cg18644543, cg19616230, cg21184174, cg23181133). Therefore, CSV file must contain beta value of 12 CpG sites and Sample."),
              br(),
              p(strong("Example input is provided here:")),
              DT::dataTableOutput("Example"),
              br(),
              h3("Prediction"),
              p("1.Import data"),
              p("Click on the “Browse” button and select your file."),
              p("2.Start prediction"),
              p("Please wait for a while. The waiting time will depend on the size of your data."),
              p("3.Observe the results"),
              p("Output table provide 2 additional columns including TB possibility and Prediction which responsible for possibility of TB infection (range from 0 to 1) and predicted TB condition (≥0.5 means TB infection, and <0.5 means Normal sample(s))")

   )
       ),

           tabPanel("Study wolkflow",
            mainPanel(
              h1("Worlkflow"),
              p("We built and followed an in-house pipeline to process the pooled DNA methylation profiling data of 9 datasets including 151 TB patients (TBs); 139 healthy controls (HCs) and 167 HIV, HIV-TB coinfection. Pooled data was then splitted into discovery and validation for independent indentifying candidate biomarkers and peformance evaluation"),
              img(src="https://github.com/nthongle317/TB_pred/blob/main/www/figure1.png?raw=true", height="110%", width="110%", align="left"),
              h6(strong("Flow chart of study procedure.")),
              h6("Abbreviations: HC: Healthy controls, TB: tuberculosis, HIV: Human immunodeficiency virus, BMA: Bayesian model averaging, LASSO: Least Absolute Shrinkage and Selection.")
   )
       ),           

           tabPanel("Prediction",
             sidebarLayout(
                sidebarPanel(
                  fileInput("file1", "Choose CSV File", 
                    accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
                    ),
                  checkboxInput("header", "Header", TRUE)
                  ),
               
               mainPanel(
                 DT::dataTableOutput("tablePred")
             )
           )
       )
   )
)
# Define server logic ----
server <- function(input, output) {
  output$tablePred=DT::renderDataTable ({
    x <- "MUVR_rf"
    y <- "rf"
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    valDT=read.csv(file$datapath, header = input$header)
    my_data=app(x, y, valDT)
    DT::datatable(my_data, 
      extensions = "Buttons",
      options = list(
        paging = TRUE,
        rownames = FALSE,
        scrollX=TRUE,
        searching = TRUE,
        ordering = TRUE,
        dom = 'l<"sep">Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        pageLength=10,
        lengthMenu=c(10,20,50,100)
        )
      )
    })

  output$Example <- DT::renderDataTable ({
    DT::datatable(
      read.csv(url("https://raw.githubusercontent.com/nthongle317/TB_pred/main/config/example.csv")),
      extensions = "Buttons",
      options = list(
        paging = TRUE,
        rownames = FALSE,
        scrollX=TRUE,
        searching = TRUE,
        ordering = TRUE,
        dom = 'l<"sep">Bfrtip',
        buttons = c('copy', 'csv'),
        pageLength=10,
        lengthMenu=c(10,20,50,100)
        )
      ) 

    })
}

# Run the app ----
shinyApp(ui = ui, server = server)  
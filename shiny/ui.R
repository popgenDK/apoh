
library(shiny)
library(markdown, lib.loc="/home/genis/R/x86_64-pc-linux-gnu-library/4.1")
library(xfun, lib.loc="/home/genis/R/x86_64-pc-linux-gnu-library/4.1")


shinyUI(fluidPage(

    # Application title
    titlePanel("apoh - Admixture Pedigrees Of Hybrids"),

    # Sidebar with a slider for different input 
    sidebarLayout(
        sidebarPanel(
            fileInput("anccoefFile", "Upload your input file with model estimates (NGSremix .anccoef output)"),

            fileInput("indfile", "Upload file with IDs for the individuals (optional)"),

            uiOutput("chooseInd"),
            
            numericInput("ncases", "Number of top compatible pedigrees to show", value=1, min=1, max=12),
        
            actionButton("run", label="Run"),

            helpText("Download the output tables:"),            
#            downloadButton("downloadPedigreePlots", "Download pdf with all possible pedigrees."),
            downloadButton("downloadIndicesTable", "Download table with summary indices for all samples."),
            downloadButton("downloadBootSupportTable", "Download table with bootstrap support for the pedigrees for each sample.")

#            downloadButton("downloadPairedAncPlots", "Download pdf with all paired ancestry given all possible pedigrees.")
            
          ),

        mainPanel(
            tabsetPanel(
                tabPanel("Usage", includeMarkdown("README.md")),
                tabPanel("Parental admixture plot", plotOutput("parentalAdmixPlot")),
                tabPanel("Summary indices table", dataTableOutput("indicesTable")),
                tabPanel("Individual ordered paired ancestries plots", plotOutput("orderedPairAncPlot")),
                tabPanel("Individual pedigrees", plotOutput("pedigreePlot", width="auto", height="auto")),
                tabPanel("Individual unordered paired ancestries plots", plotOutput("unorderedPairAncPlot")),
                tabPanel("Pedigrees bootstrap support table", dataTableOutput("bootSupport"))
            
        ))
    )
))

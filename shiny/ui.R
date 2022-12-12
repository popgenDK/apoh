
library(shiny)
library(markdown, lib.loc="/newHome/genis/R/x86_64-pc-linux-gnu-library/4.1")
library(xfun, lib.loc="/newHome/genis/R/x86_64-pc-linux-gnu-library/4.1")


shinyUI(fluidPage(

    # Application title
    titlePanel("Inference of recent admixture pedigrees"),

    # Sidebar with a slider for different input 
    sidebarLayout(
        sidebarPanel(
            fileInput("parentalQFile", "Upload your input file with parental admixture proportions"),
            fileInput("pairAncFile", "Upload your input file with paired ancestry proportions (optional)"),
            
            fileInput("indfile", "Upload file with IDs for the individuals (optional)"),

            uiOutput("chooseInd"),
            
            numericInput("ncases", "Number of top compatible pedigrees to show", value=1, min=1, max=12),
        
            actionButton("run", label="Run"),

            helpText("Download the output:"),            
            downloadButton("downloadPedigreePlots", "Download pdf with all possible pedigrees."),
            downloadButton("downloadSummaryTable", "Download table with summary of all near pedigrees."),
            downloadButton("downloadPairedAncPlots", "Download pdf with all paired ancestry given all possible pedigrees.")
            
          ),

        mainPanel(
            tabsetPanel(
                tabPanel("Usage", includeMarkdown("README.md")),
                tabPanel("Ordered paired ancestries plots", plotOutput("orderedPairAncPlot")),
                tabPanel("Pedigrees", plotOutput("pedigreePlot", width="auto", height="auto")),
                tabPanel("Summary Table", dataTableOutput("summaryTable")),
                tabPanel("Unordered paired ancestries plots", plotOutput("unorderedPairAncPlot")),
            
        ))
    )
))

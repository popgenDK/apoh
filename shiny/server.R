
library(shiny)
source("https://raw.githubusercontent.com/popgenDK/apoh/main/apohFuns.R")
source("/home/genis/github/apoh/apohFuns.R")

error <- function(txt)
    validate(  need(FALSE, txt)  ) 

shinyServer(function(input, output) {

    # functions to read paired ancestries from input files
    allPairedAnc <- reactive({
        
        f <- input$anccoefFile
        indf <- input$indfile

        x <- read_ancestries(f$datapath)[[2]]
        
        #inds <- 1:nrow(x)
        if(!is.null(indf)){
            inds <- scan(indf$datapath, what="fd")
            rownames(x) <- inds
            }
        
        return(x)
    })

    
    allParentalQ  <- reactive({
        # function to read parental ancestry from input files
        f <- input$anccoefFile
        indf <- input$indfile

        x <- read_ancestries(f$datapath)[[1]]

        if(!is.null(indf)){
            inds <- scan(indf$datapath, what="fd")
            names(x) <- inds
        }
        
        return(x)
    })

    ## get optimal widht for parental admixture plot given number of samples
    parentalPlotwidth <- reactive({

        ninds <- length(allParentalQ())

        return(min(ninds * 100,1600))
    })

    ## make plot of parental admixture proportions for all samples
    output$parentalAdmixPlot <- renderPlot({

        parentalQ <- allParentalQ()

        indlabels <- TRUE
        if(length(parentalQ) > 20) indlabels <- FALSE
        
        par(oma = c(3,0.5,0,0))
        plotParentalAdmixture(parentalQ, indlabels= indlabels)
        
    }, width=parentalPlotwidth, height=800)


    ## make table with summary indices, etc. for all samples
    output$indicesTable <- renderDataTable({
        
        nshow <- input$ncases
        indf <- input$indfile
        #nshow <- min(nshow, length(pedigrees()))

        parentalQ <- allParentalQ()
        pairedAnc <- allPairedAnc()
        
        if(!is.null(indf)){ inds <- scan(indf$datapath, what="fd")} else {inds <- NULL}

        doManySamplesIndexesTable(parentalQs = allParentalQ(), pairedAncss = allPairedAnc(), ids = inds, npeds=nshow)
    })

    
    # tab to choose wich indivdual to do
    output$chooseInd <- renderUI({
        
        indf <- input$indfile

        inds <- NULL
        if(!is.null(input$parentalQFile)){
            n <- length(allParentalQ())
            inds <- 1:n
            }
        if(!is.null(indf)) inds <- scan(indf$datapath, what="fd")
 

        selectInput('ind', 'Select individual',  c(Choose="", inds), selectize=TRUE)

    })
    
    
    # get pedigrees compatible with parental Q
    # pedigrees <- eventReactive(input$run, {
    pedigrees <- reactive({
        
        validate(
            need(!is.null(input$anccoefFile), "Upload file with model estimates!"),
            need(input$ind != "", "Select individual")
        )


        q <- allParentalQ()[[input$ind]]

        # this is a list of lists of K lenght integer vectors
        ped <- getAllSortedPedigrees(q)

        return(ped)
    })

    heightPlot1 <- reactive({

        nshow <- min(input$ncases, length(pedigrees())) + 1

        ndim <- getdim(nshow)
        
        return(ndim[1] * 400)

    })

    
    widthPlot1 <- reactive({

        nshow <- min(input$ncases, length(pedigrees())) + 1

        ndim <- getdim(nshow)
    
        return(ndim[2] * 400)

    })


    
    
    output$pedigreePlot <- renderPlot({

        nshow <- input$ncases
        nshow <- min(nshow, length(pedigrees()))

        par(mfrow=getdim(nshow + 1))

        parentalQ <- allParentalQ()[[input$ind]]
        
        plotIndependentPedigree(parentalQ, title=paste( "Independent ancestries\npedigree sample ",input$ind))
        
        for(i in 1:nshow) plotPedigree(pedigrees()[[i]], title=paste("Pedigree", i, "\nsample ", input$ind))
    }, width=widthPlot1, height=heightPlot1)
    

    
    
    widthPlot2 <- reactive({

        nshow <- min(input$ncases, length(pedigrees()))

    
        return(300 + 100 * nshow)

    })


    output$unorderedPairAncPlot <- renderPlot({

        validate(
            #need(!is.null(input$parentalQFile), "Upload file with parental admixture proportions"),
            need(!is.null(input$anccoefFile), "Upload file with model estimates!"),
            need(input$ind != "", "Select individual")
        )

        nshow <- input$ncases
        nshow <- min(nshow, length(pedigrees()))
        parentalQ <- allParentalQ()[[input$ind]]
        pairedAnc <- allPairedAnc()[input$ind,]
        #par(mfrow=getdim(nshow))

        plotEstimates(pairedAnc, parentalQ, pedigrees()[1:nshow],
                      main_title = paste("Unordered paired ancestry proportions\nsample", input$ind))
    }, width = widthPlot2, height=400)


    
    widthPlot3 <- reactive({

        nshow <- min(input$ncases, length(pedigrees()))

        
        return(400 + 200 * nshow)

    })


    
    output$orderedPairAncPlot <- renderPlot({

        
        validate(
            need(!is.null(input$anccoefFile), "Upload file with model estimates!"),
            need(input$ind != "", "Select individual")
        )

        nshow <- input$ncases
        nshow <- min(nshow, length(pedigrees()))

        parentalQ <- allParentalQ()[[input$ind]]
        #par(mfrow=getdim(nshow))

        plotOrderedEstimates(parentalQ, pedigrees()[1:nshow],
                       main_title = paste("Ordered paired ancestry proportions\nsample", input$ind))
    }, width = widthPlot3, height=400)



    output$summaryTable <- renderDataTable({
       
        validate(
            need(!is.null(input$anccoefFile), "Upload file with model estimates!"),
            need(input$ind != "", "Select individual")
        )

        parentalQ <- allParentalQ()[[input$ind]]

        d <- makePedigreesSummaryTable(parentalQ, pedigrees())

        d
    }) 




    output$downloadPedigreePlots <- downloadHandler(
        filename = function(){"admixturePedigrees.pdf"},
        content = function(fname){

            pdf(fname)
            for(i in 1:length(pedigrees()))  plotPedigree(pedigrees()[[i]],
                                                           title=paste("Pedigree", i))
            dev.off()

        }
    )
            

            
    output$downloadPairedAncPlots <- downloadHandler(
        filename = function(){"pairedAncestries.pdf"},
        content = function(fname){
        parentalQ <- allParentalQ()[[input$ind]]
        pairedAnc <- allPairedAnc()[input$ind,]

            
            pdf(fname)
            for(i in 1:length(pedigrees()))  plotEstimates(pairedAnc,
                                                            parentalQ,
                                                            pedigrees()[[i]],
                                                            main_title = paste("Paired ancestry proportions\n pedigree", i))
            dev.off()

        }
    )


    
  output$downloadSummaryTable <- downloadHandler(
       filename = function(){"allAdmixPedigreesSummary.tsv"}, 
      content = function(fname){

          parentalQ <- allParentalQ()[[input$ind]]
          pairedAnc <- allPairedAnc()[input$ind,]

          
           d <- makePedigreesSummaryTable(pairedAnc, parentalQ, pedigrees(), is.ord=TRUE)
           write.table(d, fname, sep="\t", col.names=T, row.names=F,quote=F)
    }
  )

            

    


})




library(shiny)
source("/home/genis/github/recAdmix/recAdmixFuns.R") # CHANGE TO GITHUB VERISON ONCE REP IS PUBLIC


error <- function(txt)
    validate(  need(FALSE, txt)  ) 

shinyServer(function(input, output) {

    # functions to read in the input files
    allPairedAnc <- reactive({
        
        f <- input$pairAncFile
        indf <- input$indfile
        
        x <- read_pairedancs(f$datapath)
        
        #inds <- 1:nrow(x)
        if(!is.null(indf)){
            inds <- scan(indf$datapath, what="fd")
            rownames(x) <- inds
            }
        
        return(x)
    })

    
    allParentalQ  <- reactive({
        
        f <- input$parentalQFile
        indf <- input$indfile

        x <- read_parentalancs(f$datapath)

        if(!is.null(indf)){
            inds <- scan(indf$datapath, what="fd")
            names(x) <- inds
        }
        
        return(x)
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
    #pedigrees <- eventReactive(input$run, {
    pedigrees <- reactive({
        
        validate(
            need(!is.null(input$parentalQFile), "Upload file with parental admixture proportions"),
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
        
        for(i in 1:nshow) plotPedigree3(pedigrees()[[i]], title=paste("Pedigree", i, "\nsample ", input$ind))
    }, width=widthPlot1, height=heightPlot1)
    

    
    
    widthPlot2 <- reactive({

        nshow <- min(input$ncases, length(pedigrees()))

    
        return(300 + 100 * nshow)

    })


    output$unorderedPairAncPlot <- renderPlot({

        validate(
            #need(!is.null(input$parentalQFile), "Upload file with parental admixture proportions"),
            need(!is.null(input$pairAncFile), "Upload file with paired ancestry proportions"),
            need(input$ind != "", "Select individual")
        )

        nshow <- input$ncases
        nshow <- min(nshow, length(pedigrees()))
        parentalQ <- allParentalQ()[[input$ind]]
        pairedAnc <- allPairedAnc()[input$ind,]
        #par(mfrow=getdim(nshow))

        plotEstimates3(pairedAnc, parentalQ, pedigrees()[1:nshow],
                       main_title = paste("Unordered paired ancestry proportions\nsample", input$ind))
    }, width = widthPlot2, height=400)


    
    widthPlot3 <- reactive({

        nshow <- min(input$ncases, length(pedigrees()))

        
        return(400 + 200 * nshow)

    })


    
    output$orderedPairAncPlot <- renderPlot({

        
        validate(
            need(!is.null(input$parentalQFile), "Upload file with parental admixture proportions"),
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
            need(!is.null(input$parentalQFile), "Upload file with parental admixture proportions"),
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
            for(i in 1:length(pedigrees()))  plotPedigree3(pedigrees()[[i]],
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
            for(i in 1:length(pedigrees()))  plotEstimates2(pairedAnc,
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




library(shiny)
source("https://raw.githubusercontent.com/popgenDK/apoh/main/apohFuns.R")
#source("/home/genis/github/apoh/apohFuns.R")

error <- function(txt)
    validate(  need(FALSE, txt)  ) 

shinyServer(function(input, output) {

    # functions to read paired ancestries from input files
    allPairedAnc <- reactive({
        
        f <- input$anccoefFile
        indf <- input$indfile

        if(has_boot(f$datapath)){
            x <- read_ancestries_boot(f$datapath)[[1]][[2]]
        } else {
            x <- read_ancestries(f$datapath)[[2]]
        }
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

        if(has_boot(f$datapath)){
            x <- read_ancestries_boot(f$datapath)[[1]][[1]]
        } else {
            x <- read_ancestries(f$datapath)[[1]]
        }

        if(!is.null(indf)){
            inds <- scan(indf$datapath, what="fd")
            names(x) <- inds
        }
        
        return(x)
    })

    parentalQBoots <- reactive({

        f <- input$anccoefFile

        validate(
            need(has_boot(f$datapath), "Input file has no bootstrap replicates!")
        )

        x <- read_ancestries_boot(f$datapath)[[2]]
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
        
        par(oma = c(3,2,0,0))
        plotParentalAdmixture(parentalQ, indlabels= indlabels)
        
    }, width=1000, height=600)


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

    
    ## DO TABLE WITH BOOTSTRAP SUPPORT PER PEDIGREE
    output$bootSupport <- renderDataTable({

        f <- input$anccoefFile
        nshow <- input$ncases
        indf <- input$indfile

        validate(
            need(has_boot(f$datapath), "Input file does not contain bootstrap replicates!")
        )
        #nshow <- min(nshow, length(pedigrees()))

        maxpeds <- 8
        parentalQ <- allParentalQ()
        boots <- parentalQBoots()
        pedigrees <- lapply(parentalQ, function(x) getAllSortedPedigrees(x, maxped=maxpeds))

        
        if(!is.null(indf)){ inds <- scan(indf$datapath, what="fd")} else {inds <- names(parentalQ)}

            bootsupport <- bootSupportPedigree(boots, pedigrees, parentalQ, maxped=maxpeds, showped=nshow)
            data.frame(ID=inds, bootsupport)
    })


    ###############################
    ### DOWNLOAD OUTPUT TABLES ####
    ###############################
    
  output$downloadIndicesTable <- downloadHandler(
      filename = function(){"summaryIndicesTable.tsv"}, 
      content = function(fname){

          nshow <- input$ncases
          indf <- input$indfile

          parentalQ <- allParentalQ()
          pairedAnc <- allPairedAnc()

                  
        if(!is.null(indf)){ inds <- scan(indf$datapath, what="fd")} else {inds <- NULL}

          d <- doManySamplesIndexesTable(parentalQs = allParentalQ(), pairedAncss = allPairedAnc(), ids = inds, npeds=nshow)

          write.table(d, fname, sep="\t", col.names=T, row.names=F,quote=F)
    }
  )


    
    output$downloadBootSupportTable <- downloadHandler(
       filename = function(){"bootSupportTable.tsv"}, 
      content = function(fname){

          f <- input$anccoefFile
          nshow <- input$ncases
          indf <- input$indfile

          validate(
            need(has_boot(f$datapath), "Input file does not contain bootstrap replicates!")
        )
        #nshow <- min(nshow, length(pedigrees()))
          
          maxpeds <- 8
          parentalQ <- allParentalQ()
          boots <- parentalQBoots()
          pedigrees <- lapply(parentalQ, function(x) getAllSortedPedigrees(x, maxped=maxpeds))

        
          if(!is.null(indf)){ inds <- scan(indf$datapath, what="fd")} else {inds <- names(parentalQ)}

          bootsupport <- bootSupportPedigree(boots, pedigrees, parentalQ, maxped=maxpeds, showped=nshow)
          d <- data.frame(ID=inds, bootsupport)

          write.table(d, fname, sep="\t", col.names=T, row.names=F,quote=F)
    }
  )


    


})



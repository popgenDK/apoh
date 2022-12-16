
# this is default palette
#tol_bright palette From Paul Tol: https://personal.sron.nl/~pault/
colorpal <- c('#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#4477AA','#BBBBBB',
          "#332288", "#CC3311", "#EE7733", "#EE3377", "#000000") # this line I added to have colors for higher K

printHelp <- function(){

    cat("Admixture Pedigrees Of Hybrids (apoh).\nToolset to explore recent admixture pedigrees from paired ancestry proportions. \n\n
Arguments:\n\n
\t\t-i --infile:\t path to .anccoeff file with inferred paired ancestries (as outputted by NGSremix)\n
\t\t-o --outdir:\t directory to save output (must not exist).\n
\t\t--ids:\t file with id for individual (optional).\n
\t\t--npedigrees:\t number of top compatible recent admixture pedigrees to show (default 2).\n
\t\t--forcedir:\t use outdir even if it already exists (default FALSE).\n
\t\t--colpal:\t color palette as comma separated list of hexadecimal color codes (e.g. '#EE6677,#228833,#CCBB44...') \n
\t\t-h:\t print help and exit.\n
")

}


readArgs <- function(args){

    pars <- list(infile=NULL,
                 outdir=NULL,
                 ids = NULL,
                 npedigrees = 2,
                 forcedir = FALSE,
                 colpal = NULL
                 )

    for(i in seq(1, length(args), 2)
        ){

        if(args[i] == "--infile" | args[i] == "-i"){
            pars$infile <- args[i+1]
        } else if(args[i] == "--outdir" | args[i] == "-o"){
            pars$outdir <- args[i+1] 
        } else if(args[i] == "--ids"){
            pars$ids <- args[i+1]   
        } else if(args[i] == "--npedigrees"){
            pars$npedigrees <- as.integer(args[i+1])
        } else if(args[i] == "--forcedir"){
            pars$forcedir <- as.logical(args[i+1])
        } else if(args[i] == "--colpal"){
            pars$colpal <- args[i+1]
        } else if(args[i] == "-h" | args[i] == "--help"){
            printHelp()
            stop("Printed help and exited due to -h flag, not really an error.")
        } else {
            printHelp()
            stop("Unkonwn argument ", args[i], ", see above accepted arguments.\n")
        }
    }

    
    if(is.null(pars$infile)){
        printHelp()
        stop("Missing path to input .anccoeff file (--infile is mandatory).")
    } else if(is.null(pars$outdir)){
        stop("Missing path to output directory (--outdir is mandatory).")
    } else if(file.exists(pars$outdir) & pars$forcedir == FALSE){
        stop("Output directory specified by --outdir already exists. Please remove it, use a new directory or run with --forcedir TRUE.")
    }

    return(pars)

}



# some short helper funcitons
combi <- function(K){
  gtools::combinations(K,2,repeats.allowed=T)
}

permi <- function(K){
  gtools::permutations(K,2,repeats.allowed=T)
}

num <- function(x)
  (sqrt(8*x+1)-1)/2



#### to read paired and parental ancestries from ngsremix

read_ancestries <- function(f){

    d <- read.table(f, h=T, stringsAsFactors=F)

    # extract paired ancestries
    pairAnc <- as.numeric(unlist(strsplit(d$paired_est, split=",")))
    pairAnc <- matrix(pairAnc, ncol=length(pairAnc) / nrow(d), byrow=T)
    rownames(pairAnc) <- d$ind
    
    # extract parental admixture
    parentalAnc <- lapply(strsplit(d$parental_est, split=","), as.numeric)
    k <- length(parentalAnc[[1]]) / 2
    parentalAnc <- lapply(parentalAnc, function(x) list(x[1:k], x[(k+1):(k*2)]))
    names(parentalAnc) <- d$ind
    
    return(list(parentalAnc, pairAnc))

}



getAncestorsQ <- function(q, g=3){
    # takes admixture proporitons, discretizes it so they fit to
    # unadmixed individuals g generations back (default 3),
    # returns list of vector of K length with close combinations of ancestros for each K
    
    n <- 2^g
    nanc_raw <- q * n
    nanc <- round(nanc_raw)

    if(sum(nanc) < n){
        i <- which.min(round(nanc_raw) - nanc_raw)
        nanc[i] <- nanc[i] + 1
    }else if(sum(nanc) < n){
        i <- which.max(round(nanc_raw) - nanc_raw)
        nanc[i] <- nanc[i] - 1
    }

    return(nanc)

}



getAncestorSet <- function(q, g=3){
    # get close cases of n unadmixed individual for
    # each of the K ancestral populations
    # g generations ago (default 3)
    # q is K length vector of admixtue proportions
    

    n <- 2^g
    k <- length(q)

    nanc_raw <- q * n
    nanc_base<- floor(nanc_raw)

    nanc_add <- n - sum(nanc_base)
    if(nanc_add > 0){
        w <- gtools::combinations(k, nanc_add)
    
        nanc <- t(apply(w, 1, function(x){ y <- nanc_base; y[x] <- y[x] + 1;return(y) }))
    } else {
        nanc <- matrix(nanc_base, nrow=1)
    }
    return(nanc)
    
}



getAllPedigrees <- function(parentalQ){
    # get all possible pedigrees compatible with
    # parental admixture proporitons
    
    
    nanc_all <- lapply(parentalQ, getAncestorSet)

    combs <- expand.grid(1:nrow(nanc_all[[1]]), 1:nrow(nanc_all[[2]]))
    
    pedigrees <- lapply(1:nrow(combs), function(x) list(nanc_all[[1]][combs[x,1],], nanc_all[[2]][combs[x,2],]))

    return(pedigrees)

}



pAFromPQ <- function(parentalQ){
    # estiamte paired Ancestry from parental Q
    # there must be a most elegant way to do it but I was tired when I wrote it
   
    q1 <- parentalQ[[1]]
    q2 <- parentalQ[[2]]

    k <- length(q1)

    pairedAnc <- rep(0, k*(k-1)/2+k)
    names(pairedAnc) <- apply(combi(k),1,paste0, collapse="")

    for(i in 1:k) for(j in 1:k) pairedAnc[paste0(sort(c(i,j)), collapse="")] <-  pairedAnc[paste0(sort(c(i,j)), collapse="")] + q1[i] * q2[j]

    pairedAnc
    
}



pAFromIQ <- function(q){
    # calculate expected paired ancestry from individual
    # admixture proporitons assuming independence

    k <- length(q)

    pairedAnc <- rep(0, k*(k-1)/2+k)
    names(pairedAnc) <- apply(combi(k),1,paste0, collapse="")

    for(i in 1:k) for(j in 1:k) pairedAnc[paste0(sort(c(i,j)), collapse="")] <- pairedAnc[paste0(sort(c(i,j)), collapse="")] + q[i] * q[j]

    pairedAnc

    
}


parentalQFromPed <- function(pedigree){
    # calculate expected parental admixture pedigree proporitons
    # given a pedigree, defined as a set of unadmixed ancestors for each parent
    # in practice list of two N vectors, where N is number of ancestors
    # 8 usually when we go 5 generations up

    parentalQ <- lapply(pedigree, function(x) x/sum(x))

    parentalQ
}



pAFromPed <- function(pedigree){
    # calculate expected paired ancestry proporitons
    # given a pedigree, defined as a set of unadmixed ancestors for each parent
    # in practice list of two N vectors, where N is number of ancestors
    # 8 usually when we go 5 generations up

    parentalQ <- parentalQFromPed(pedigree)
    pairedAnc <- pAFromPQ(parentalQ)

    pairedAnc

}



orderedPAFromPed <- function(pedigree){
    # takes pedigree
    # returns ordered paired ancestry proporitons (K * K) length

    k <- length(pedigree[[1]])
    per <- permi(k)
    parentalQ <- parentalQFromPed(pedigree)
    qq  <- parentalQ[[1]][per[,1]] * parentalQ[[2]][per[,2]]

    qq
}



orderedPAFromParentalQ <- function(parentalQ){
    # takes parental admixutre proportions
    # returns ordered paired ancestry proporitons (K * K) length

    k <- length(parentalQ[[1]])
    per <- permi(k)
    qq  <- parentalQ[[1]][per[,1]] * parentalQ[[2]][per[,2]]

    qq
}



orderedPAIndependent <- function(parentalQ){
    # takes parental admixture
    # return ordered paired ancestry assuming independence
    # ie assuming both parents have equal ancestry

    q <- do.call("+", parentalQ) / 2
    k <- length(q)

    a <- t(permi(k))
    qq <- q[a[1,]] * q[a[2,]]

    qq

}




pedigreeToCharacter <- function(pedigree){
    ## converts a pedigree object (list of k lenght integer vectors)
    # to character in format "n11,n12,...n1k|n21,n22,...n2k"

    char <- paste(lapply(pedigree, paste, collapse=","), collapse="|")

    char

}



 parentalQToCharacter <- function(parentalQ, decimals=3){
    # converts a parental admixture object (list of k lenght numeric vectors)
    # to character in format "n11,n12,...n1k|n21,n22,...n2k"

    parentalQ <- lapply(parentalQ, round, digits=decimals)

    char <- paste(lapply(parentalQ, paste, collapse=","), collapse="|")

    char
}



pairedAncToCharacter <- function(pairedAnc, decimals=3){
    # converts paired ancestry numeric vector to character comma delimited
    # rounding to decimals (default 3)

    pairedAnc <- lapply(pairedAnc, round, digits=decimals)

    char <- paste(pairedAnc, collapse=",")

    char
    
}



entropy <- function(x){
    # calculate entroyp of discrete probability distribution

    h <- - sum(x * log2(x))
    return(h)
}



js <- function(p, q){
    # general function to calculate jensen-shannon divergence
    # between two discrete probability distributions
    # p and q must be equal lenght vectors

    js <- sqrt(entropy((p + q) / 2) - ((entropy(p) + entropy(q))/2))

    return(js)
}



jsRecentOldParentalQ <- function(parentalQ){
    # given estimate of parental admixture proportions, returns js divergence
    # between the paired ancestries expected given parentalQ, and the paired ancestry
    # expected if independent paired ancestry (older admixture)

    p <- orderedPAIndependent(parentalQ) + 1e-8
    p <- p / sum(p)
    q <- orderedPAFromParentalQ(parentalQ) + 1e-8
    q <- q / sum(q)
    
    jss <- js(p, q)
    return(jss)
}



jsParentalQPed <- function(parentalQ, pedigree){
    # give js divergence between expected paired ancestries given certain
    # estimate of parentalQ, and between the expected given a certain
    # recent admixture pedigree

    p <- orderedPAFromParentalQ(parentalQ) + 1e-8
    p <- p / sum(p)

    q <- orderedPAFromPed(pedigree) + 1e-8
    q <- q / sum(q)
    
    jss <- js(p, q)
    return(jss)

}



jsBestFit <- function(parentalQ){

    pedigrees <- getAllPedigrees(parentalQ)
    dists <- sapply(pedigrees, jsParentalQPed, parentalQ = parentalQ)

    dists <- c(dists, jsRecentOldParentalQ(parentalQ))

    return(min(abs(dists)))
}



jsRatio <- function(parentalQ){

    js1 <- jsBestFit(parentalQ)
    js2 <- jsRecentOldParentalQ(parentalQ)

    idx <- 1 - (js1/js2)
    return(idx)
}



minHetAnc <- function(parentalQ){

    panc <- orderedPAFromParentalQ(parentalQ)
    k <- length(parentalQ[[1]])
    a <- permi(k)

    hets <- panc[a[,1] != a[,2]]
    return(min(hets[hets>0.01]))
}


tFromParentalQ <- function(parentalQ){

    p <- minHetAnc(parentalQ)
    t <- -log2(p) + 1

    return(t)
}


inconsistencyIndex <- function(parentalQ, pairedAnc){
    # use difference between estiamted unordered paired ancestries and
    # paired ancestries expected from estiamted parental ancestries.
    # difference measured as jensen-shannon divergence
    
    
    pairedAnc2 <- pAFromPQ(parentalQ)

    # add and norm so zero values don't mess it up
    pairedAnc <- pairedAnc + 1e-8
    pairedAnc <- pairedAnc/sum(pairedAnc)

    pairedAnc2 <- pairedAnc2 + 1e-8
    pairedAnc2 <- pairedAnc2/sum(pairedAnc2)
    
    idx <- js(pairedAnc, pairedAnc2)

    return(idx)
}


rankPedigreesParentalQdist <- function(parentalQ, pedigrees){


    dists <- sapply(pedigrees, jsParentalQPed, parentalQ=parentalQ)
    return(order(dists, decreasing=F))
}


getAllSortedPedigrees <- function(parentalQ){
    # given parental Q, get all compatible pedigrees by order of
    # more to less distant to paretnalQ

    pedigrees <- getAllPedigrees(parentalQ)

    ord <- rankPedigreesParentalQdist(parentalQ, pedigrees)

    pedigrees <- pedigrees[ord]
    return(pedigrees)

}


doManySamplesIndexesTable <- function(parentalQs, pairedAncss = NULL, npeds=2, ids=NULL){
    # fucntion to summarise results for many sample using
    # 2 or 3 informative indexes

    if(is.null(ids)) ids <- names(parentalQs)
    names(parentalQs) <- ids
    rownames(pairedAncss) <- ids
    
    pedigrees <- lapply(parentalQs, function(x) getAllSortedPedigrees(x)[1:npeds])#[1:2]
    pedSupport <- t(sapply(names(pedigrees), function(x) sapply(pedigrees[[x]], jsParentalQPed, parentalQ=parentalQs[[x]])))

    df <- data.frame("SampleID" = ids)

    if(!is.null(pairedAncss)) df["InconsistencyIndex"] <- format(sapply(ids, function(x) inconsistencyIndex(parentalQs[[x]], pairedAncss[x,])), digits=2,nsmall=2,scientific=F)
    df["Distance to independent pedigree"] <- format(sapply(parentalQs, jsRecentOldParentalQ), digits=2,nsmall=2,scientific=F)
    for(i in 1:npeds) df[paste("Distance to pedigree",i)] <- format(pedSupport[,i], digits=2,nsmall=2,scientific=F)
    
    df["Admixture index"] <-  sapply(1:length(parentalQs), function(x) ifelse(df[x,"Distance to independent pedigree"] > df[x,"Distance to pedigree 1"], tFromParentalQ(parentalQs[[x]]), ""))
    
    return(df)
}



getdim <- function(nshow){
    # function to get number of rows and col
    # for multipanle plot (ie argument for par mfrow)
    # based on number of cases we want to show
    # i.e. number of plots to make

    
    nrow <- ceiling(sqrt(nshow))
    ncol <- ceiling(nshow/nrow)
    
    return(c(nrow,ncol))
}



plotOrderedPairedAnc <- function(ordPairedAnc, x = 0.25, y = 0, title = "Paired ancestry proportions",
                          write_vals = FALSE, new_plot=FALSE,
                          cex.title=1.5, cex.labs = 1,
                          colpal=colorpal){
    
    qq <- ordPairedAnc
    K<- sqrt(length(qq))
    com <- permi(K)

    if(new_plot)
        plot.new()
  
    l1 <- x - 0.2
    r1 <- x - 0.05
    l2 <- x + 0.05
    r2 <- x + 0.2
      
    rect(l1, y + c(0,cumsum(qq[1:(length(qq)-1)])),r1, cumsum(qq),
         col=colpal[com[,1]],border=NA)
    
    rect(l2, y + c(0,cumsum(qq[1:(length(qq)-1)])),r2,cumsum(qq),
         col=colpal[com[,2]],border=NA)
    
    text(x, y + 1.05, title ,xpd=T,font=2,cex=cex.title)
    
    if(write_vals)
        text(x+0.25,cumsum(qq)-qq/2,round(qq*100,1), cex=cex.labs) 
    
}



plotOrderedEstimates <- function(parentalQ, pedigrees,
                           main_title = "Paired ancestry proportions\n(ordered)",
                           cex.main=1.5,
                           cex.lab=1,
                           colpal=colorpal, showsupport = FALSE){

    on.exit(par(mar=c(5,4,4,2) +0.1))

    #par(mar=c(5.1, 4.1, 6.1, 4.1))
 
    k <- length(parentalQ[[1]])

    qq1 <- orderedPAFromParentalQ(parentalQ)
    qq2 <- orderedPAIndependent(parentalQ)
    
    xlim <- c(0, 1 + 0.5 * length(pedigrees))
    
    plot(1, type="n", xlab="", ylab="", xlim=xlim, ylim=c(0, 1), xaxt="n", yaxt="n", bty="n")

    xlim2 <- xlim + c(0.1, -0.1)
    x <- seq(xlim2[1], xlim2[2], length.out= 2 + length(pedigrees))# + length(pedigrees))

    text(x = mean(xlim), y = 1.1, labels=main_title, cex=cex.main, xpd=NA)

    par(xpd=NA)
    plotOrderedPairedAnc(qq1, x=x[1], title="", colpal=colpal)
    plotOrderedPairedAnc(qq2, x=x[2], title="", colpal=colpal)
    
    for(i in 1:length(pedigrees)) plotOrderedPairedAnc(orderedPAFromPed(pedigrees[[i]]), x=x[i+2], title="", colpal=colpal)

    axis(side=2, line=1, cex.axis=1.2)
    text(x=xlim[1] - 0.35 - 0.02 * length(pedigrees), y=0.5, srt=90, labels="Paired ancestry proportions", cex=1.4, xpd=NA)
    
    xlabs <- c("Model estimates",
               "Expected under\nindependent pedigree",
               paste("Expected under\npedigree", 1:length(pedigrees)))

    text(x=x, y=-0.085, labels=xlabs, cex=cex.lab)

    if(showsupport){
        
        pedsupport <- sapply(pedigrees, jsParentalQPed, parentalQ = parentalQ)
        indepsupport <- jsRecentOldParentalQ(parentalQ)

        xlabs2 <- paste("Distance to\nestimate = ", round(c(indepsupport,pedsupport),2))
        text(x=x[-1], y=1.1, labels=xlabs2, cex=cex.lab)
    }

}



plotPairedAnc <- function(pairedAnc, x = 0.25, y = 0, title = "Paired ancestry proportions",
                          write_vals = FALSE, new_plot=FALSE,
                          cex.title=1.5, cex.labs = 1,
                          colpal=colorpal){
    # plot a single case of paired ancestry proportions
    qq <- pairedAnc
    K<- num(length(qq))
    com <- combi(K)

    if(new_plot)
        plot.new()
  
    l1 <- x - 0.1
    r1 <- x
    l2 <- x
    r2 <- x + 0.1
      
    rect(l1, y + c(0,cumsum(qq[1:(length(qq)-1)])),r1, cumsum(qq),
         col=colpal[com[,1]],border=NA)
    
    rect(l2, y + c(0,cumsum(qq[1:(length(qq)-1)])),r2,cumsum(qq),
         col=colpal[com[,2]],border=NA)
    
    text(x, y + 1.05, title ,xpd=T,font=2,cex=cex.title)
    
    if(write_vals)
        text(x+0.25,cumsum(qq)-qq/2,round(qq*100,1), cex=cex.labs) 

}
 


plotEstimates <- function(pairedAnc, parentalQ, pedigrees,
                           main_title = "Paired ancestry proportions",
                           cex.main=1.5,
                          cex.lab=1,
                          colpal=colorpal,
                          showConsistency=FALSE){
    ## wrappter of plotPairedancestries that plots 3 + n cases:
    # 1. paired ancestry assuming independence (i.e. no recent admixture)
    # 2. pared ancestry props inferred from the data
    # 3. paired ancestry props expected given the inferred parental admixture proprotions
    # 4. paired ancestry props expected given the proposed admixture pedigree

    on.exit(par(mar=c(5,4,4,2) +0.1))

    #par(mar=c(7.1, 6.1, 4.1, 4.1))

    k <- length(parentalQ[[1]])

    qq1 <- pairedAnc
    qq2 <- pAFromPQ(parentalQ)

    a <- combi(k)
    q <- tapply(c(qq1,qq1),a,sum)/2

    qq3 <- pAFromIQ(q)
    
    if(length(qq1) != length(qq2)) stop("Paired ancestry and parental Q imply different K.")

    
    xlim <- c(0, 1 + 0.5 * length(pedigrees))
    
    plot(1, type="n", xlab="", ylab="", xlim=xlim, ylim=c(0, 1), xaxt="n", yaxt="n", bty="n")

    xlim2 <- xlim + c(0.1, -0.1)
    x <- seq(xlim2[1], xlim2[2], length.out= 3 + length(pedigrees))# + length(pedigrees))

    text(x = mean(xlim), y = 1.1, labels=main_title, cex=cex.main, xpd=NA)

    par(xpd=NA)
    plotPairedAnc(qq1, x=x[1], title="", colpal=colpal)
    plotPairedAnc(qq2, x=x[2], title="", colpal=colpal)
    plotPairedAnc(qq3, x=x[3], title="", colpal=colpal)
    
    for(i in 1:length(pedigrees)) plotPairedAnc(pAFromPed(pedigrees[[i]]), x=x[i+3], title="", , colpal=colpal)

    axis(side=2, line=0.5, cex.axis=1.2)
    text(x=xlim[1] - 0.25 - 0.02 * length(pedigrees), y=0.5, srt=90, labels="Paired ancestry proportions", cex=1.4)
    
    xlabs <- c("Paired\nancestry\nmodel\nestimates", "Parental\nadmixture\nmodel\nestimates",
               "Expected\nunder\nindependent\npedigree",
               paste("Expected\nunder\npedigree", 1:length(pedigrees)))

    text(x=x, y=-0.2, labels=xlabs, cex=cex.lab)

    if(showConsistency){

        idx <- inconsistencyIndex(pairedAnc = pairedAnc, parentalQ = parentalQ)
        text(x=mean(x[1:2]), y=1.1, labels=paste("Inconsistency index =", round(idx,2)), xpd=NA, cex=cex.lab)

    }
}



nextXpoints <- function(xpoints)
    tapply(xpoints, rep(1:(length(xpoints)/2), each=2), mean)



nextQanc <- function(qanc)
    t(sapply(1:(nrow(qanc)/2), function(x) qanc[x*2-1,]/2 + qanc[x*2,]/2))



plotPedigree <- function(pedigree, title = "", cex.title=1.5,
                         colpal=colorpal){
    # plot a single admixture pedigree, taking as input the number of unadmixed indivduals
    # for each ancestral population in each branch
    # pedigree: list of two vectors of length K that sum to 8

    nanc1 <- pedigree[[1]]
    nanc2 <- pedigree[[2]]

    on.exit(par(mar=c(5,4,4,2) +0.1))

    par(mar=c(0, 2.1, 4.1, 2.1))
    plot.new()

    text(x=0.5, y=1.12, labels=title, xpd=NA, cex=cex.title)
    
    if(sum(nanc1) != sum(nanc2)) stop("Different number of ancestors for each parent")
    if(length(nanc1) != length(nanc2)) stop("Different number of K for each parent")
    
    n <- sum(nanc1)
    g <- log2(n * 2) + 1
    k <- length(nanc1)

    qanc1 <- matrix(unlist(sapply(1:k, function(x) rep(diag(k)[x,], nanc1[x]*2))), ncol=k, nrow=n*2, byrow=T)
    qanc2 <- matrix(unlist(sapply(1:k, function(x) rep(diag(k)[x,], nanc2[x]*2))), ncol=k, nrow=n*2, byrow=T)
    qanc <- rbind(qanc1, qanc2)

    xpoints <- seq(0, 1, length.out= n*2)
    ypoints <- seq(1, 0, length.out = g + 1)
    ydist <- 1/(g*2)
    
    # first generation
    rect(xleft = rep(rep(xpoints, each=2) + c(-1/(n * 6), 1/(n*32)), each=k),
         xright = rep(rep(xpoints, each=2) + c(-1/(n * 32), 1/(n*6)), each=k),
         ybottom = ypoints[1] - ydist/2 +
             rbind(rep(0,nrow(qanc)), apply(t(qanc[,1:(k-1)]), 2, cumsum)) * ydist,
         ytop=ypoints[1]  - ydist/2 + apply(t(qanc), 2, cumsum) * ydist,
         col=colpal[1:k], xpd=NA)

    if(g != 2) {
        # recurisvely plot remaining generations up to parents
        for(i in 2:(g-1)){
            
            xpoints0 <- xpoints
            xpoints <- nextXpoints(xpoints)
            qanc <- nextQanc(qanc)
            
            # draw lines to join parents to next gen kids
            segments(x0=xpoints0, x1= rep(xpoints, each=2), y0=ypoints[i-1]-ydist/2, y1=ypoints[i] + ydist/2)
            
            # plot generation
            rect(xleft = rep(rep(xpoints, each=2) + c(-1/(n * 6), 1/(n*32)), each=k),
                 xright = rep(rep(xpoints, each=2) + c(-1/(n * 32), 1/(n*6)), each=k),
                 ybottom = ypoints[i] - ydist/2 + rbind(rep(0,nrow(qanc)), apply(t(qanc[,1:(k-1)]), 2, cumsum)) * ydist,
                 ytop=ypoints[i]  - ydist/2 + apply(t(qanc), 2, cumsum) * ydist,
                 col=colpal[1:k], xpd=NA)
        }
    }
    
    # plot individual paired ancestry proportions given pedigree    
    xpoints0 <- xpoints
    xpoints <- nextXpoints(xpoints)

    segments(x0=xpoints0, x1= rep(xpoints, each=2), y0=ypoints[g-1]-ydist/2, y1=ypoints[g] + ydist/2)
  
    l1 <- xpoints - 1/(n * 6)
    r1 <- xpoints - 1/(n * 32)
    l2 <- xpoints + 1/(n*32)
    r2 <- xpoints + 1/(n*6)

    y <- ypoints[g] -  ydist/2

    per <- permi(k)
    parentalQ <- parentalQFromPed(pedigree)
    qq  <- parentalQ[[1]][per[,1]] * parentalQ[[2]][per[,2]]
    
    rect(l1, y + c(0,cumsum(qq[1:(length(qq)-1)])) * ydist,r1, y + cumsum(qq) * ydist,
         col=colpal[per[,1]], border=NA)
    rect(l1, y, r1, y+ydist, col=NA)
    
    rect(l2, y + c(0,cumsum(qq[1:(length(qq)-1)])) * ydist,r2, y + cumsum(qq) * ydist,
         col=colpal[per[,2]], border=NA)
    rect(l2, y, r2, y+ydist, col=NA)
    
}



plotIndependentPedigree <- function(parentalQ, n=8, title = "", cex.title=1.5,
                        colpal=colorpal){
    ### plots pedigree assuming ancient admxiutre where all
    # ancestors have equal admixture proporitons
    # so paired anceistreis are independent given admixutre porporitons

    on.exit(par(mar=c(5,4,4,2) +0.1))

    par(mar=c(0, 2.1, 4.1, 2.1))
    plot.new()

    text(x=0.5, y=1.12, labels=title, xpd=NA, cex=cex.title)

    
    g <- log2(n * 2) + 1
    q <- do.call("+", parentalQ) / 2
    k <- length(q)

    a <- t(permi(k))
    qq <- q[a[1,]] * q[a[2,]]

    qanc <- matrix(rep(q, times=n*2), nrow=n*2, byrow=T)

    xpoints <- seq(0, 1, length.out= n*2)
    ypoints <- seq(1, 0, length.out = g + 1)
    ydist <- 1/(g*2)

    # plot first generation
    rect(xleft = rep(rep(xpoints, each=2) + c(-1/(n * 6), 1/(n*32)), each=k),
         xright = rep(rep(xpoints, each=2) + c(-1/(n * 32), 1/(n*6)), each=k),
         ybottom = ypoints[1] - ydist/2 + rbind(rep(0,nrow(qanc)), apply(t(qanc[,1:(k-1)]), 2, cumsum)) * ydist,
         ytop=ypoints[1]  - ydist/2 + apply(t(qanc), 2, cumsum) * ydist,
         col=colpal[1:k], xpd=NA)

       # recurisvely plot remaining generations up to parents
    for(i in 2:(g-1)){
        
        xpoints0 <- xpoints
        xpoints <- nextXpoints(xpoints)
        qanc <- nextQanc(qanc)
        
        # draw lines to join parents to next gen kids
        segments(x0=xpoints0, x1= rep(xpoints, each=2), y0=ypoints[i-1]-ydist/2, y1=ypoints[i] + ydist/2)
        
        # plot generation
        rect(xleft = rep(rep(xpoints, each=2) + c(-1/(n * 6), 1/(n*32)), each=k),
         xright = rep(rep(xpoints, each=2) + c(-1/(n * 32), 1/(n*6)), each=k),
         ybottom = ypoints[i] - ydist/2 + rbind(rep(0,nrow(qanc)), apply(t(qanc[,1:(k-1)]), 2, cumsum)) * ydist,
         ytop=ypoints[i]  - ydist/2 + apply(t(qanc), 2, cumsum) * ydist,
         col=colpal[1:k], xpd=NA)
    }

        # plot individual paired ancestry proportions given pedigree    
    xpoints0 <- xpoints
    xpoints <- nextXpoints(xpoints)

    segments(x0=xpoints0, x1= rep(xpoints, each=2), y0=ypoints[g-1]-ydist/2, y1=ypoints[g] + ydist/2)
  
    l1 <- xpoints - 1/(n * 6)
    r1 <- xpoints - 1/(n * 32)
    l2 <- xpoints + 1/(n*32)
    r2 <- xpoints + 1/(n*6)

    y <- ypoints[g] -  ydist/2
    
    rect(l1, y + c(0,cumsum(qq[1:(length(qq)-1)])) * ydist,r1, y + cumsum(qq) * ydist,
         col=colpal[t(a)[,1]], border=NA)
    rect(l1, y, r1, y+ydist, col=NA)
    
    rect(l2, y + c(0,cumsum(qq[1:(length(qq)-1)])) * ydist,r2, y + cumsum(qq) * ydist,
         col=colpal[t(a)[,2]], border=NA)
    rect(l2, y, r2, y+ydist, col=NA)


}




plotParentalAdmixture <- function(parentalQ,
                                  inds=NULL,
                                  indlabels=TRUE,
                                  main.title="Parental admixture proportions"){
    ### do plot of parental admixture for many samples

    if(is.null(inds)) inds <- names(parentalQ)

    parentalQ <- do.call("rbind", lapply(parentalQ, unlist))
    
    k <- ncol(parentalQ) / 2
    barplot(t(parentalQ), col=colorpal[1:k], space=0,
            xlab="", yaxt="n", xaxt="n", border=NA,
            main="")
    title(main=main.title, line=1, cex.main=1.5, xpd=NA)

    abline(v=1:nrow(parentalQ), col="white", lwd=1, xpd=FALSE)

    axis(side=2, at=c(0,0.5,0.98), labels=c(0,0.5,"1 0"), line=-0, cex.axis=1.2)
    axis(side=2, at=c(1.02,1.5,2), labels=c("",0.5,1), line=-0, cex.axis=1.2)
    text(label="Parent 1", x=-1 * nrow(parentalQ) / 10, y=0.5, srt=90, xpd=NA, cex=1.5)
    text(label="Parent 2", x=-1 * nrow(parentalQ) / 10, y=1.5, srt=90, xpd=NA, cex=1.5)

    if(indlabels) text((1:length(inds)) - 0.5, y=-0.1, labels=inds, xpd=NA, cex=1.3,srt=90)

}


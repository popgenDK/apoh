
# this is default palette
#tol_bright palette From Paul Tol: https://personal.sron.nl/~pault/
colorpal <- c('#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#4477AA','#BBBBBB',
          "#332288", "#CC3311", "#EE7733", "#EE3377", "#000000") # this line I added to have colors for higher K

printHelp <- function(){

    cat("Admixture Pedigrees Of Hybrids (APOH).\nToolset to explore recent admixture pedigrees from paired ancestry proportions. \n\n
Arguments:\n\n
\t\t-i --inancestries:\t path to file with inferred paired ancestries (as outputted by NGSremix)\n
\t\t-o --outprefix:\t prefix to write output files to.\n
\t\t--ind:\t file with id for individual.\n
\t\t-h:\t print help and exit.\n
")

}


readArgs <- function(args){

    pars <- list(inancestries,
                 outprefix=NULL,
                 ind = NULL
                 )

    for(i in seq(1, length(args), 2)){

        if(args[i] == "--"){
            pars$inancestries <- args[i+1]
        } else if(args[i] == "--outprefix"){
            pars$outprefix <- args[i+1]
        } else if(args[i] == "--ind"){
            pars$ind <- args[i+1]
        } else if(args[i] == "-h"){
            printHelp()
            stop("Printed help and exited due to -h flag, not really an error.")
        } else {
            printHelp()
            stop("Unkonwn argument ", args[i], ", see above accepted arguments.\n")
        }
    }

    
    if(is.null(pars$parentalQ)){
        printHelp()
        stop("Missing path to file with parental admixture proportions (--parentalQ).")
    } else if(is.null(pars$pairedAnc)){
        stop("Missing path to file with paired ancestry proportions (--pairedAnc).")
    } else if(is.null(pars$outprefix)){
        stop("Missing prefix to save output to (--outprefix).")
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

read_parentalancs <- function(f){
    # takes path to paired ancestries file
    # reads and returns parental admxiture proportions as a list

    l <- read_ancestries(f)
    return(l[[1]])

}




read_pairedancs <- function(f){
    # takes path to paired ancestries file
    # reads and returns paired ancestrie proportions as a list

    l <- read_ancestries(f)
    return(l[[1]])

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



distParentalQPed <- function(parentalQ, pedigree){
    # euclidean distance between given parental admixture proportions and
    # parental admixture proprotions expected given a pedigree

    
    eParentalQ <- parentalQFromPed(pedigree)

    d <- dist(rbind(unlist(parentalQ), unlist(eParentalQ)), method="euclidean")

    d

}



distPairedAncPed <- function(pairedAnc, pedigree){
    # euclidean distance between given paired ancestry and
    # paired ancestry expected given a pedigree

    ePairedAnc <- pAFromPed(pedigree)

    d <- dist(rbind(pairedAnc, ePairedAnc), method="euclidean")
    
    d
    
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



rankPedigreesParentalQdist <- function(parentalQ, pedigrees){
    # return order of pedigrees from less to more distance
    # of expected parental Q to the observed ones

    dists <- sapply(pedigrees, distParentalQPed, parentalQ=parentalQ)
    
    ord <- order(dists)

    ord

}



kl <- function(p, q){
    # general function to calculate kullback-leiber divergence
    # between two discrete probability distributions
    # p and q must be equal lenght vectors
    
    a <- 0
    for(i in 1:length(p)) a = a + p[i] * log2(p[i] / q[i])
    return(a)
}


js <- function(p, q){
    # general function to calculate jensen-shannon divergence
    # between two discrete probability distributions
    # p and q must be equal lenght vectors
    
    m <- (p + q) / 2
    jss <- 0.5 * kl(p, m) + 0.5 * kl(q,m)
    return(jss)
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


consistencyIndex <- function(parentalQ, pairedAnc){
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


getAllSortedPedigrees <- function(parentalQ){
    # given parental Q, get all compatible pedigrees by order of
    # more to less distant to paretnalQ

    pedigrees <- getAllPedigrees(parentalQ)

    ord <- rankPedigreesParentalQdist(parentalQ, pedigrees)

    pedigrees <- pedigrees[ord]
    return(pedigrees)

}



makePedigreesSummaryTable <- function(parentalQ, pedigrees=NULL){#function(pairedAnc, parentalQ, pedigrees=NULL){

    if(is.null(pedigrees)) pedigrees <- getAllPedigrees(parentalQ)

    #distsPairedAnc <- sapply(pedigrees, distPairedAncPed, pairedAnc=pairedAnc)
    distsParentalQ <- sapply(pedigrees, distParentalQPed, parentalQ=parentalQ)

    
    charPeds <- sapply(pedigrees, pedigreeToCharacter)
    charParQ <- sapply(lapply(pedigrees, parentalQFromPed), parentalQToCharacter)
    #charPairAnc <- sapply(lapply(pedigrees, pAFromPed), pairedAncToCharacter)
    ids <- paste("pedigree", 1:length(pedigrees), sep="")

    # add independent pedigree
    q <- do.call("+", parentalQ) / 2
    indParentalQ <- list(q,q)
    
    distsParentalQ <- c(distsParentalQ, dist(rbind(unlist(parentalQ), unlist(indParentalQ)), method="euclidean"))
    charPeds <- c(charPeds, "Not recent")
    charParQ <- c(charParQ, parentalQToCharacter(indParentalQ))
    ids <- c(ids, "Independent pedigree")

    # ordered paired ancestries from parental 
    jsDivergence <- unlist(sapply(pedigrees, jsParentalQPed, parentalQ=parentalQ))
    jsDivergence <- c(jsDivergence, jsRecentOldParentalQ(parentalQ))
    paParental <- orderedPAFromPed(parentalQ)

    # order
    ord <- order(distsParentalQ)
    #ord <- order(jss)
    
    pedSummary <- data.frame(`Pedigree ID`=ids[ord], Pedigree=charPeds[ord],
                               `Expected parental Q` = charParQ[ord],
                               `Observed parental Q` = parentalQToCharacter(parentalQ),
                             `Distance expected to observed parental Q` = distsParentalQ[ord],
                             `JS divergence expected to estimated` = jsDivergence[ord])
                               #`Expected paired ancestry` = charPairAnc[ord],
                               #`Observed paired ancestry` = pairedAncToCharacter(pairedAnc),
                               #`Distance expected to observed paired ancestry` = distsPairedAnc[ord])

    pedSummary

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
                           colpal=colorpal){


    on.exit(par(mar=c(5,4,4,2) +0.1))

    par(mar=c(5.1, 4.1, 6.1, 4.1))
 
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

    axis(side=2, line=2.2)
    text(x=xlim[1] - 0.5 - 0.02 * length(pedigrees), y=0.5, srt=90, labels="Paired ancestry proportions", cex=cex.lab)
    
    xlabs <- c("Model estimates",
               "Expected under\nindependent ancestries",
               paste("Expected under\npedigree", 1:length(pedigrees)))

    text(x=x, y=-0.1, labels=xlabs, cex=cex.lab)



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



plotEstimates <- function(pairedAnc, parentalQ, cex.titles=1){
    ## wrappter of plotPairedancestries that plots 3 cases:
    # 1. paired ancestry assuming independence (i.e. no recent admixture)
    # 2. pared ancestry props inferred from the data
    # 3. paired ancestry props expected given the inferred parental admixture proprotions
    k <- length(parentalQ[[1]])

    qq1 <- pairedAnc
    qq2 <- pAFromPQ(parentalQ)

    a <- combi(k)
    q <- tapply(c(qq1,qq1),a,sum)/2

    qq3 <- pAFromIQ(q)
    
    if(length(qq1) != length(qq2)) stop("Paired ancestry and parental Q imply different K.")

    plot(1, type="n", xlab="", ylab="", xlim=c(-0.5, 2.5), ylim=c(0, 1), xaxt="n", yaxt="n", bty="n")

    par(xpd=NA)
    plotPairedAnc(qq3, x=0, title="Paired ancestry given\nindividual admixture proportion", cex.title=cex.titles)
    plotPairedAnc(qq1, x=1, title="Inferred paired ancestry", cex.title=cex.titles)
    plotPairedAnc(qq2, x=2, title="Paired ancestry given inferred\nparental admixture proportions", cex.title=cex.titles)

}


plotEstimates2 <- function(pairedAnc, parentalQ, pedigree,
                           main_title = "Paired ancestry proportions",
                           cex.main=1.5,
                           cex.subtitles=0.8){
    ## wrappter of plotPairedancestries that plots 4 cases:
    # 1. paired ancestry assuming independence (i.e. no recent admixture)
    # 2. pared ancestry props inferred from the data
    # 3. paired ancestry props expected given the inferred parental admixture proprotions
    # 4. paired ancestry props expected given the proposed admixture pedigree

    on.exit(par(mar=c(5,4,4,2) +0.1))

    par(mar=c(5.1, 4.1, 6.1, 4.1))
 
    k <- length(parentalQ[[1]])

    qq1 <- pairedAnc
    qq2 <- pAFromPQ(parentalQ)

    a <- combi(k)
    q <- tapply(c(qq1,qq1),a,sum)/2

    qq3 <- pAFromIQ(q)
    qq4 <- pAFromPed(pedigree)
    
    if(length(qq1) != length(qq2)) stop("Paired ancestry and parental Q imply different K.")

    xlim <- c(-0.5, 2.5)
    
    plot(1, type="n", xlab="", ylab="", xlim=xlim, ylim=c(0, 1), xaxt="n", yaxt="n", bty="n")

    x <- seq(xlim[1], xlim[2], length.out=4)

    text(x = mean(xlim), y = 1.2, labels=main_title, cex=cex.main, xpd=NA)
    
    par(xpd=NA)
    plotPairedAnc(qq3, x=x[1], title="Given individual\nadmixture proportions", cex.title=cex.subtitles)
    plotPairedAnc(qq1, x=x[2], title="Model estimates", cex.title=cex.subtitles)
    plotPairedAnc(qq2, x=x[3], title="Given parental\nadmixture proportions", cex.title=cex.subtitles)
    plotPairedAnc(qq4, x=x[4], title="Given admixture\npedigree", cex.title=cex.subtitles)


}



plotEstimates3 <- function(pairedAnc, parentalQ, pedigrees,
                           main_title = "Paired ancestry proportions",
                           cex.main=1.5,
                           cex.lab=1){
    ## wrappter of plotPairedancestries that plots 3 + n cases:
    # 1. paired ancestry assuming independence (i.e. no recent admixture)
    # 2. pared ancestry props inferred from the data
    # 3. paired ancestry props expected given the inferred parental admixture proprotions
    # 4. paired ancestry props expected given the proposed admixture pedigree

    on.exit(par(mar=c(5,4,4,2) +0.1))

    par(mar=c(7.1, 6.1, 4.1, 4.1))

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
    plotPairedAnc(qq1, x=x[1], title="")
    plotPairedAnc(qq2, x=x[2], title="")
    plotPairedAnc(qq3, x=x[3], title="")
    
    for(i in 1:length(pedigrees)) plotPairedAnc(pAFromPed(pedigrees[[i]]), x=x[i+3], title="")

    axis(side=2, line=0.5)
    text(x=xlim[1] - 0.35 - 0.02 * length(pedigrees), y=0.5, srt=90, labels="Paired ancestry proportions", cex=cex.lab)
    
    xlabs <- c("Paired\nancestry", "Parental\nadmixture",
               "Independent\npedigree",
               paste("Pedigree", 1:length(pedigrees)))

    text(x=x, y=-0.1, labels=xlabs, cex=cex.lab)
    
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

    # recurisvely plot remaining generations
    for(i in 2:g){
        
        xpoints0 <- xpoints
        xpoints <- nextXpoints(xpoints)
        qanc <- nextQanc(qanc)
        
        # jdraw oin parents to next gen kids by arrows
        segments(x0=xpoints0, x1= rep(xpoints, each=2), y0=ypoints[i-1]-ydist/2, y1=ypoints[i] + ydist/2)
        
        # plot generation 2  
        rect(xleft = rep(rep(xpoints, each=2) + c(-1/(n * 6), 1/(n*32)), each=k),
         xright = rep(rep(xpoints, each=2) + c(-1/(n * 32), 1/(n*6)), each=k),
         ybottom = ypoints[i] - ydist/2 + rbind(rep(0,nrow(qanc)), apply(t(qanc[,1:(k-1)]), 2, cumsum)) * ydist,
         ytop=ypoints[i]  - ydist/2 + apply(t(qanc), 2, cumsum) * ydist,
         col=colpal[1:k], xpd=NA)

    }
}




plotPedigree2 <- function(pedigree, title = "", cex.title=1.5,
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
    
    qq <- pAFromPed(pedigree)

    l1 <- xpoints - 1/(n * 6)
    r1 <- xpoints - 1/(n * 32)
    l2 <- xpoints + 1/(n*32)
    r2 <- xpoints + 1/(n*6)

    y <- ypoints[g] -  ydist/2

    com <- combi(k)

    rect(l1, y + c(0,cumsum(qq[1:(length(qq)-1)])) * ydist,r1, y + cumsum(qq) * ydist,
         col=colpal[com[,1]], border=NA)
    rect(l1, y, r1, y+ydist, col=NA)
    
    rect(l2, y + c(0,cumsum(qq[1:(length(qq)-1)])) * ydist,r2, y + cumsum(qq) * ydist,
         col=colpal[com[,2]], border=NA)
    rect(l2, y, r2, y+ydist, col=NA)
    
}



plotPedigree3 <- function(pedigree, title = "", cex.title=1.5,
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

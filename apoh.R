

whereDir <- function(){
    # function to get directory where scripts are, so apohFuns.R can be sourced when run from any folderfrom outside. Assumes apoh.R and apohFuns.R are in the same folder
    # made by krishang
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file"
    match <- grep(needle, cmdArgs)
    tf <- unlist(strsplit(cmdArgs[match], "="))[2]
    d <- dirname(tf)
    return(d)
}


d <- whereDir()
source(paste(d, "apohFuns.R", sep="/"))

args <- commandArgs(trailingOnly=T)
pars <- readArgs(args)
f <- pars$infile
outdir <- pars$outdir
npeds <- pars$npedigrees
if(!is.null(pars$colpal)) colorpal <- unlist(strsplit(pars$colpal, ","))

cat("Starting to run apoh, will use ", f, " as input\n\n\n")

# read paired ancestry proporitons
if(!has_boot(f)){
    res <- read_ancestries(f)
} else {
    res0 <- read_ancestries_boot(f)
    res <- res0[[1]]
    boots <- res0[[2]]
}
    
parentalAnc <- res[[1]]
pairAnc <- res[[2]]

ids <- names(parentalAnc)
if(!is.null(pars$ids)) ids <- scan(pars$ids, what="fld")


# create output directory
dir.create(outdir)


# make plot of parentaladmxiture proportions
outpng1 <- paste0(outdir, "/parentalAdmixture.png")
bitmap(outpng1, w=8, h=4, res=300)
plotParentalAdmixture(parentalAnc, inds=ids)
dev.off()


# make tabe of summary indices for all samples
outtsv <- paste0(outdir, "/allSummaryIndices.tsv")
d <- doManySamplesIndexesTable(parentalQs = parentalAnc,
                               pairedAncss = pairAnc,
                               npeds = npeds,
                               ids=ids)
write.table(x=d,file=outtsv, row.names=F, col.names=T, quote=F, sep="\t")


# for each individual, make directory and save different stuff
for(i in 1:length(ids)){

    id <- ids[i]
    outdir2 <- paste0(outdir, "/", id)
    dir.create(outdir2)

    pedigrees <- getAllSortedPedigrees(parentalAnc[[i]])[1:npeds]
    # make ordered paied ancestry plot
    outpng <- paste0(outdir2, "/", id, "_orderedPairAnc.png")
    bitmap(outpng, w= 4 + 2 * npeds, h=4, res=300)
    plotOrderedEstimates(parentalAnc[[i]], pedigrees,
                         main_title = paste("Ordered paired ancestry proportions sample", id),
                         colpal=colorpal,
                         showsupport=FALSE)
    dev.off()
    
    # make unordered paired ancestry proportions plot
    outpng <- paste0(outdir2, "/", id, "_unorderedPairAnc.png")
    bitmap(outpng, w= 4.5 + 1.5 * npeds, h=4, res=300)
    plotEstimates(pairAnc[i,], parentalAnc[[i]], pedigrees,
                         main_title = paste("Unordered paired ancestry proportions sample", id),
                         colpal=colorpal,
                         showConsistency=FALSE)
    dev.off()
    
    # make pedigree plots
    for(j in 1:npeds){
        outpng <- paste0(outdir2, "/", "pedigree",j,"_sample_",id,".png")
        bitmap(outpng, w=4,h=4, res=300)
        plotPedigree(pedigree=pedigrees[[j]],
                     title=paste("Pedigree",j,"sample",i),
                     colpal=colorpal)
        dev.off()

    }
}


## if has bootstrap, save table with boostrap support for pedigrees
if(has_boot(f)){

    outboot <- paste0(outdir, "/pedigreeBootSupport.tsv")
    maxpeds <- 8

    pedigrees <- lapply(parentalAnc, function(x) getAllSortedPedigrees(x, maxped=maxpeds))

    bootsupport <- bootSupportPedigree(boots, pedigrees, maxped=maxpeds)
    bootsupport <- data.frame(ID=ids, bootsupport)
    write.table(x=bootsupport,file=outboot, row.names=F, col.names=T, quote=F, sep="\t")
    cat("Wrote bootstrap results in ", outboot, "\n")

}


cat("Finished running apoh, all ouput can be found in ", outdir,"\n")

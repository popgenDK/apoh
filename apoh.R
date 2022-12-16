

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

f <- pars$inancestries
outdir <- pars$outdir
dir.create(outdir)

res <- read_ancestries(f)
parentalAnc <- res[[1]]
pairAnc <- res[[2]]

ids <- names(parentalAnc)
if(!is.null(pars$ids)) ids <- scan(pars$ids, what="fld")

# make plot of parentaladmxiture proportions
outpng1 <- paste0(outdir, "/parantalAdmixture.png")
bitmap(outpng1, w=8, h=6, res=300)
plotParentalAdmixture(parentalAnc, inds=ids)
dev.off()


# make tabe of summary indices for all samples


# for each individual, make directory and save different stuff
for(i in 1:length(ids)){

    id <- ids[i]
    outdir2 <- paste0(outdir, "/", id)
    dir.create(outdir2)

    # make parental ancestires plot

    # make paired ancestry proportions plot

    # 

}



# prefix path for ouptu files
outprefix <- pars$outprefix

# get set of compatible pedigrees
pedigrees <- getAllSortedPedigrees(parentalQ)

# make summary table
summaryTable <- makePedigreesSummaryTable(pairedAnc, parentalQ, pedigrees, is.ord=TRUE)

# write output table
outtable <- paste0(outprefix, "_summaryAdmixPedigrees.tsv")
#cat(outtable)
write.table(summaryTable, outtable, col.names=T, row.names=F, quote=F, sep="\t")

plotEstimates3(

# make plot of all pedigrees
outpdf1 <- paste0(outprefix, "_allAdmixPedigrees.pdf")
pdf(outpdf1)
for(i in 1:length(pedigrees)) plotPedigree3(pedigrees[[i]], title=paste("Pedigree", i, ind))
dev.off()


# make plot of all paired ancestries for pedigrees
outpdf2 <- paste0(outprefix, "_allPairedAncObsExp.pdf")
pdf(outpdf2)
for(i in 1:length(pedigrees)) plotEstimates2(pairedAnc, parentalQ, pedigrees[[i]], main_title=paste("Paired ancestries\nind",ind,"pedigree", i))
dev.off()



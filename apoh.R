

whereDir <- function(){
    # function to get directory where scripts are, so recAdmixFuns.R can be sourced when run from any folderfrom outside. Assumes recAdmix.R and recAdmixFuns.R are in the same folder
    # made by krishang
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file"
    match <- grep(needle, cmdArgs)
    tf <- unlist(strsplit(cmdArgs[match], "="))[2]
    d <- dirname(tf)
    return(d)
}



d <- whereDir()
source(paste(d, "recAdmixFuns.R", sep="/"))


args <- commandArgs(trailingOnly=T)
pars <- readArgs(args)

# expected format: parental Q 2 rows K columns tab or space delimited file
# paired ancestry 1 row K*(K-1)/2+K length tab or space delimited
parentalQ <-  as.list(as.data.frame(t(as.matrix(read.table(pars$parentalQ)))))
pairedAnc <- scan(pars$pairedAnc, what=.3)

ind <- pars$ind

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



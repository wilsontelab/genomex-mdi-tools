#--------------------------------------------------------------------
# utilities related to parsing parts of the genome, coverage, bins, etc.
#--------------------------------------------------------------------

# genome span coverage
readCountsToRpkm <- function(
    count, 
    sizeBp = NULL, sizeKb = NULL,         # one must be provided
    totalReads = NULL, totalReadsM = NULL # one must be provided
){
    if(is.null(sizeKb)) sizeKb <- sizeBp / 1e3
    if(is.null(totalReadsM)) totalReadsM <- totalReads / 1e6
    count / sizeKb / totalReadsM
}

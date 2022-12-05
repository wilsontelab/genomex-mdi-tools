#=====================================================================================
# create a set of pseudo-randomized fragment endpoints that samples every genome
# exactly twice, once as a fragment start, once as a fragment end
# output the paired reads for each fragment
#=====================================================================================

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# parse environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
source(file.path(rUtilDir, 'utilities.R'))
checkEnvVars(list(
    string = c(
        'DATA_NAME',
        'ACTION_DIR',
        'MODULES_DIR',
        'GENOME',
        'GENOME_CHROMS',
        'GENOME_FASTA',
        'CHROM_FASTA_DIR'     
    ),
    integer = c(
        "READ_LENGTH"
    )
))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
# options(warn = 2)
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'genome'), c('chroms', 'faidx'))
#=====================================================================================

#=====================================================================================
# initialize a set of distributions of genomic fragment insert sizes
#-------------------------------------------------------------------------------------
minFragStart <- 0
maxFragStart <- 1 * env$READ_LENGTH - 1
minFragEnd   <- 2 * env$READ_LENGTH - 2 # prevents read overruns in too-small molecules
maxFragEnd   <- 3 * env$READ_LENGTH - 3
fragStarts   <- minFragStart:maxFragStart # 0-referenced
fragEnds     <- minFragEnd:maxFragEnd     
read1Starts  <- fragStarts + 1 # 1-referenced
read1Ends    <- read1Starts + env$READ_LENGTH - 1
nDist <- 1000
dists <- lapply(1:nDist, function(i) {
    read2Ends <- sample(fragEnds) + 1
    list(
        read2Starts = read2Ends - env$READ_LENGTH + 1,
        read2Ends   = read2Ends
    )
})
qual <- paste0(rep("F", env$READ_LENGTH), collapse = "")
#=====================================================================================

#=====================================================================================
# pull and print all matching read pair mimics
#-------------------------------------------------------------------------------------
setCanonicalChroms()
loadFaidx()
for(chrom in canonicalChroms){
    chromLength <- getChromLength(chrom)
    blockStarts <- seq(1, chromLength - maxFragEnd, env$READ_LENGTH) # jump around the chromosome to help unsure mappability in bwa stream
    blockStarts <- blockStarts[sample(length(blockStarts))]

    blockStarts <- blockStarts[1]

    for(blockStart in blockStarts){
        blockSeq <- toupper(getRefSeq(chrom, blockStart + c(minFragStart, maxFragEnd)))
        if(grepl('[^N]', blockSeq)){ # skip blocks with nothing but N bases
            dist <- dists[[sample.int(nDist, 1)]]
            ids <- paste0(
                "@", 
                chrom, 
                "_", 
                blockStart + read1Starts - 1, 
                "_", 
                blockStart + dist$read2Ends - 1
            )
            read1 <- paste(
                ids,
                substring(blockSeq, read1Starts, read1Ends),
                "+",
                qual,
                sep = "\n"
            )
            read2 <- paste(
                ids,
                substring(blockSeq, dist$read2Starts, dist$read2Ends),
                "+",
                qual,
                sep = "\n"
            )
            writeLines(as.vector(rbind(read1, read2))) # output interleaved read pairs
        }
    }
}
#=====================================================================================

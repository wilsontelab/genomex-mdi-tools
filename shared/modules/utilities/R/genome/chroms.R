#----------------------------------------------------------
# genome manipulations
#----------------------------------------------------------
canonicalChroms <- character()
chromIndex      <- list() # 1-reference chromosome indices, e.g., chr3 => 3 ...
revChromIndex   <- list() # ... and vice versa

# parsers for restricting work to properly ordered canonical chromosomes
setCanonicalChroms <- function(){ 
    GENOME_CHROMS <- strsplit(env$GENOME_CHROMS, "\\s+")[[1]]
    addChrom <- function(chr){
        chrom <- paste0("chr", chr)
        if(chrom %in% GENOME_CHROMS) canonicalChroms <<- c(canonicalChroms, chrom)
    } 
    sapply(c(1:90, 'X'), addChrom)
    if(is.null(env$SUPPRESS_CHR_Y)) addChrom('Y') # chrY included unless specifically excluded
    if(!is.null(env$USE_CHR_M))     addChrom('M') # chrM excluded unless specifically included
    maxI <- length(canonicalChroms)
    for(i in 1:maxI){
        chrom <- canonicalChroms[i]
        chromIndex[[chrom]] <<- i # special handling of unmapped reads
        revChromIndex[[i]] <<- chrom
    }
    chromIndex[['*']] <<- 99 # special handling of unmapped reads
    revChromIndex[[99]] <<- '*'
}

# construct a map of all concatenated chromosomes from a genome fai file
# must be called after setCanonicalChroms
# requires library(bit64)
loadChromSizes <- function(windowSize = 1000){ 
    message("loading chrom sizes")
    index_file <- paste(env$GENOME_FASTA, 'fai', sep = ".")
    x <- fread(
        index_file, 
        header = FALSE, 
        sep = "\t", 
        stringsAsFactors = FALSE,
        col.names  = c('chrom',     'nChromBases', 'offset',  'lineSeqLen', 'lineLen'), 
        colClasses = c('character', 'integer64',   'numeric', 'integer',    'integer')
    ) 
    x[, nChromWindows := as.integer64(floor((nChromBases - 1) / windowSize) + 1)]
    getRunningCum <- function(n){
        x <- cumsum(n)
        c(as.integer64(0), x[1:(length(x) - 1)])        
    }
    x[, ":="(
        chromIndex = unlist(chromIndex[chrom]),
        nBasesBefore   = getRunningCum(nChromBases),
        nWindowsBefore = getRunningCum(nChromWindows)
    )]
    x[, ":="(
        nBasesThrough   = nBasesBefore   + nChromBases,
        nWindowsThrough = nWindowsBefore + nChromWindows
    )]
    x <- x[, .SD, .SDcols = c(
        "chromIndex","chrom",
        "nChromBases",  "nBasesBefore",  "nBasesThrough",
        "nChromWindows","nWindowsBefore","nWindowsThrough"
    )]
    setkey(x, chromIndex)
    x
}
expandChromWindows <- function(chromSizes){
    message("expanding chrom windows for faster lookup")
    chromSizes[, .(
        chrom = chrom,
        windowIndex = 1:as.integer(nChromWindows) - 1
    ), keyby = .(chromIndex)]
}

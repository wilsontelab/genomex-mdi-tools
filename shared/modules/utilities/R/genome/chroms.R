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

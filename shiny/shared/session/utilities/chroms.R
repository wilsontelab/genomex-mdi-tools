# nodes are codified into an integer64 for streamlined comparison
# this function expands integer64 nodes out to chrom/strand/pos, returned as a data.table
parseSignedNodes <- function(chromSizes, nodes, side = NULL, canonical = FALSE) {
    genomeIs <- abs(nodes) # 1-referenced, per initialize_windows.pl
    chromIs <- Vectorize(function(i) which(chromSizes$nBasesThrough >= i)[1])(genomeIs) # sapply does not work with integer64!
    dt <- chromSizes[chromIs][, .(
        chrom = chrom, 
        chromIndex = chromIndex,
        refPos = as.integer(genomeIs - nBasesBefore),
        strand = ifelse(nodes > 0, "+", "-")
    )]    
    if(canonical) setnames(dt, c("cChrom","cChromIndex","cRefPos","cStrand"))
    if(!is.null(side)) setnames(dt, paste0(names(dt), side))
    dt
}

# collape chrom/strand/pos to integer64 nodes, returned as a vector of integer64
getSignedNode <- function(chromSizes, chromIndex, coordinate, strand, add = 0) {
    chromSizes[chromIndex, nBasesBefore] + coordinate + add * ifelse(strand == "+", 1, -1)
}

# use integer64 nodes to determine if a pair of nodes is on the canonical genome strand
isCanonicalStrand <- function(nodePair){
    o <- order(abs(nodePair))
    nodePair[o[1]] > 0 # always returns the same ordered nodes independent of strand (it is ~unimportant what that order is)
}
isCanonicalStrand_vectorized <- Vectorize(function(node1, node2){
    isCanonicalStrand(c(node1, node2))
})

# chromosome lists and sizes
listSourceChromosomes <- function(genome, metadata = NULL, force = FALSE){
    if(genome$source == "UCSC") listUcscChromosomes(  genome = genome$genome, force = force)
                           else listCustomChromosomes(list(genome = genome, metadata = metadata))
}
listCanonicalChromosomes <- function(genome, force = FALSE){
    chromosomes <- listSourceChromosomes(genome, force = force)
    req(chromosomes)
    chroms <- chromosomes$chromosome
    arabic <- c(1:100, "X", "Y", "M")
    roman <- c(as.character(as.roman(1:100)), "Y", "M")
    isUpper  <- any(startsWith("CHR", chroms))
    isArabic <- any(grepl("^CHR\\d", toupper(chroms)))
    allowed <- if(!isUpper &&  isArabic) paste0("chr", arabic)
          else if( isUpper &&  isArabic) paste0("CHR", arabic)
          else if(!isUpper && !isArabic) paste0("chr", roman)
          else if( isUpper && !isArabic) paste0("CHR", roman)
          else chroms
    standardized <- allowed[allowed %in% chroms]
    required <- sort(chroms[!grepl('_', chroms)])
    if(length(required) == 0) return(chroms) # weird genomes where every chrom fails, just return them all as is
    unique(c(required[!(required %in% standardized)], standardized))
}
getChromosomeSizes <- function(genome, metadata = NULL){
    chroms <- listCanonicalChromosomes(genome)
    sizes <- listSourceChromosomes(genome)[, .(chromosome, size)]
    sizes <- sapply(chroms, function(chrom) sizes[chromosome == chrom, size])
    ends <- cumsum(as.numeric(sizes))
    starts <- c(1, ends - 1)[1:length(ends)]
    stripes <- rep(c("odd", "even"), 60)[1:length(chroms)]
    gieStain <- if(is.null(metadata) || length(metadata) == 0){
        genomes <- genome$genome
        stripes
    } else if({
        x <- getCustomCompositeType(list(genome = genome, metadata = metadata))
        isTruthy(x) && x == "UCSC"
    }) {
        delimiter <- getCustomCompositeDelimiter(metadata)
        genomes <- sapply(chroms, function(chrom) strsplit(chrom, delimiter)[[1]][2])
        uniqueGenomes <- unique(genomes)
        c(
            rep(c("odd1", "even1"), 60)[1:sum(genomes == uniqueGenomes[1])],
            rep(c("odd2", "even2"), 60)[1:sum(genomes == uniqueGenomes[2])]
        )
    } else {
        genomes <- genome$genome
        stripes
    }
    data.table(
        genome = genomes,
        chromStart = starts,
        chromEnd = ends,
        size = sizes,
        name = chroms,
        gieStain = gieStain
    )   
}
getGenomeSize <- function(genome){
    req(genome)
    chroms <- listSourceChromosomes(genome)
    canonical <- listCanonicalChromosomes(genome)
    chroms[chromosome %in% canonical, sum(as.integer64(size))]
}
isProperChromosome <- function(chrom_){
    startsWith(toupper(chrom_), "CHR")
}

# for plots that aggregate chromosome bins, mask the last bins that may have partial data
maskLastChromBins <- function(x, chromCol = "chrom"){
    I <- x[, .I[length(.I)], by = chromCol][[2]]
    x[I, y := NA]
}

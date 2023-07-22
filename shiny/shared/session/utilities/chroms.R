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

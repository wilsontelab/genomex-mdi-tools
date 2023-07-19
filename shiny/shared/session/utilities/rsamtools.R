#----------------------------------------------------------------------
# support functions for reading browser track data from TABIX-indexed files
#----------------------------------------------------------------------

# read a tbi index file into an Rsamtools tabix object
# cache it since the index load is relatively slow
getCachedTabix <- function(bgzFile, cacheDir = NULL, create = FALSE, force = FALSE, ttl = CONSTANTS$ttl$day){
    startSpinner(session, message = "loading tabix")
    fileName <- basename(bgzFile)
    if(is.null(cacheDir)) cacheDir <- dirname(bgzFile)
    rdsFile <- file.path(cacheDir, paste(fileName, "rds", sep = "."))
    if(!file.exists(rdsFile) || create){
        startSpinner(session, message = "loading tabix...")
        unlink(rdsFile)
        saveRDS(Rsamtools::TabixFile(bgzFile), file = rdsFile)
    }
    rdsFile <- loadPersistentFile(
        file = rdsFile,
        force = create || force,
        ttl = ttl
    )
    persistentCache[[rdsFile]]$data
}

# query a bgz or other tabixed file by the current window coordinates
# pass in an object from getCachedTabix
# always use standardized column names when applicable: chrom, start, end, name, score, strand
getTabixRangeData <- function(tabix, coord, col.names = NULL, colClasses = NULL){
    req(coord$chromosome, coord$chromosome != "all")
    gRange <- GenomicRanges::GRanges(
        seqnames = coord$chromosome, 
        ranges = IRanges::IRanges(start = as.integer(coord$start), end = as.integer(coord$end))
    )
    fread(
        text = Rsamtools::scanTabix(tabix, param = gRange)[[1]],
        col.names = col.names,
        colClasses = colClasses,
        header = FALSE
    )
}

# for files compressed as runs of bins with the same score, expand the runs out to individual bins
# pass in an object from getTabixRangeData
expandTabixBinRuns <- function(dt, binSize, stranded = TRUE, na.value = 0, 
                               useBinCenter = TRUE, minusStrandNeg = TRUE){
    if(nrow(dt) == 0) return(data.table(strand = character(), x = integer(), y = double()))
    bins <- if(stranded){
        dt <- dt[, .(
            start_ = as.integer(seq(start, end, binSize)),
            score = score
        ), by = .(strand, start)]
        allBins <- as.data.table(dt[, expand.grid(
            start_ = seq(min(start_, na.rm = TRUE), max(start_, na.rm = TRUE), binSize), 
            strand = c("+", "-"), 
            stringsAsFactors = FALSE
        )])
        merge(allBins, dt, by = c("strand", "start_"), all.x = TRUE)
    } else {
        dt <- dt[, .(
            start_ = as.integer(seq(start, end, binSize)),
            score = score,
            strand = "."
        ), by = .(start)]
        allBins <- dt[, .(
            start_ = seq(min(start_, na.rm = TRUE), max(start_, na.rm = TRUE), binSize)
        )]
        merge(allBins, dt, by = c("start_"), all.x = TRUE)
    }
    bins[is.na(score), score := na.value]
    if(useBinCenter) bins[, start_ := as.integer(start_ + binSize / 2)]
    if(stranded && minusStrandNeg) bins[strand == "-", score := -score]
    bins[, .(strand = strand, x = start_, y = score)][order(x)]
}

# if too many bins, reduce the number of plotted XY points
# pass in an object from expandTabixBinRuns, or any similar format with no missing bins
aggregateTabixBins <- function(bins, track, coord, scalar, aggFn = mean){
    if(nrow(bins) == 0) return(data.table(strand = character(), point = integer(), x = integer(), y = double()))
    bins[, point := as.integer((1:.N - 1) / scalar), by = .(strand)]
    bins[, 
        .(
            x = median(as.double(x)), 
            y = aggFn(as.double(y), na.rm = FALSE)
        ), 
        by = .(strand, point)
    ]    
}

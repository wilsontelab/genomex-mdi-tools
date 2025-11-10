#----------------------------------------------------------------------
# support functions for reading browser track data from TABIX-indexed files
#----------------------------------------------------------------------

# if needed, bgzip a file to file.bgz; does nothing if
#   file is a bgz file, i.e., ends with .bgz
#   file.bgz already exists
# return the bgz file path for getCachedTabix
bgzipForTabix <- function(file){
    if(endsWith(file, "bgz")) return(file)
    bgzFile <- paste0(file, ".bgz")
    if(!file.exists(bgzFile)) Rsamtools::bgzip(file)
    bgzFile
} 

# read (and if needed, create) a tbi index file into an Rsamtools tabix object
# cache it since the index load is relatively slow
# expects at least BED3 format if we are to create the tabix index
# otherwise, getCachedTabix and getTabixRangeData enforce no requirements on column content
getCachedTabix <- function(bgzFile, cacheDir = NULL, create = FALSE, index = FALSE, force = FALSE, ttl = CONSTANTS$ttl$day){
    req(file.exists(bgzFile))
    tryCatch({
        startSpinner(session, message = "loading tabix")
        fileName <- basename(bgzFile)
        if(is.null(cacheDir)) cacheDir <- dirname(bgzFile)
        indexFile <- paste(bgzFile, "tbi", sep = ".")
        rdsFile <- file.path(cacheDir, paste(fileName, "rds", sep = "."))
        indexExists <- file.exists(indexFile)
        rdsExists   <- file.exists(rdsFile)
        indexNewerThanRds <- indexExists && rdsExists && (file.info(indexFile)$mtime > file.info(rdsFile)$mtime)
        bgzNewerThanIndex <- indexExists && file.info(bgzFile)$mtime > file.info(indexFile)$mtime
        newIndexNeeded <- !indexExists || index  || bgzNewerThanIndex
        newRdsNeeded   <- !rdsExists   || create || indexNewerThanRds || newIndexNeeded
        if(newIndexNeeded || newRdsNeeded){
            if(newIndexNeeded){
                unlink(indexFile)
                startSpinner(session, message = "indexing bgz...")
                Rsamtools::indexTabix(
                    bgzFile, 
                    seq = 1,
                    start = 2,
                    end = 3
                )
            }
            startSpinner(session, message = "loading tabix...")
            unlink(rdsFile)
            saveRDS(Rsamtools::TabixFile(bgzFile), file = rdsFile)
        }
        rdsFile <- loadPersistentFile(
            file = rdsFile,
            force = newRdsNeeded || force,
            ttl = ttl
        )
        persistentCache[[rdsFile]]$data
    }, error = function(e){
        print(e)
        stopSpinner(session)
        req(FALSE)
    })
}

# query and filter a bgz or other tabixed file against the current window coordinates
# pass in an object from getCachedTabix
# always use standardized column names when provided: chrom, start, end, name, score, strand
# otherwise, caller is responsible for parsing columns of the returned data.table
getTabixRangeData <- function(tabix, coord, col.names = NULL, colClasses = NULL, skipChromCheck = FALSE){
    req(coord$chromosome, skipChromCheck || startsWith(toupper(coord$chromosome), "CHR"))
    gRange <- GenomicRanges::GRanges(
        seqnames = coord$chromosome, 
        ranges = IRanges::IRanges(start = as.integer(coord$start), end = as.integer(coord$end))
    )

    lines <- Rsamtools::scanTabix(tabix, param = gRange)[[1]]
    if(length(lines) == 1) lines <- paste0(lines, "\n")
    if(is.null(col.names)) fread( # fread does not seem to offer any acceptable NA/NULL value for col.names
        text = lines,
        colClasses = colClasses,
        header = FALSE,
        sep = "\t"
    ) else fread(
        text = lines,
        col.names = col.names,
        colClasses = colClasses,
        header = FALSE,
        sep = "\t"
    )
}

# when the column content of getTabixRangeData is not enforced a priori by an app or track,
# use getTabixRangeData() %>% parseTabixBedFormat() to coerce data.table to one of:
#   chrom, start, end (BED3)
#   chrom, start, end, name  (BED4 v1)
#   chrom, start, end, score (BED4 v2)
#   chrom, start, end, name, score (BED5)
#   chrom, start, end, name, score, strand (BED6)
# where the input data.table must already be one of the above formats (plus any extra columns, which are dropped)
# i.e., we simply recognize, name and trim the BED column format based on column data types, we do not create it by guessing at column identities
parseTabixBedFormat <- function(dt){
    if(nrow(dt) == 0) return(data.table(
        chrom   = character(), 
        start   = integer(), 
        end     = integer(), 
        name    = character(), 
        score   = double(), 
        strand  = character()
    )) 
    ncol <- ncol(dt)  
    req(ncol >= 3)
    colNames <- c("chrom", "start", "end")
    checkCol4 <- function(dt){
        if(ncol(dt) < 4){
            dt
        } else if(typeof(dt[[4]]) == "character"){
            colNames <<- c(colNames, "name")   
            dt     
        } else {
            colNames <<- c(colNames, "score")   
            dt[, 1:4]         
        }
    }
    checkCol5 <- function(dt){
        if(ncol(dt) < 5){
            dt
        } else if(typeof(dt[[5]]) != "character"){
            colNames <<- c(colNames, "score")   
            dt     
        } else {   
            dt[, 1:4]         
        }
    }
    checkCol6 <- function(dt){
        if(ncol(dt) < 6){
            dt
        } else if(typeof(dt[[6]]) == "character"){
            colNames <<- c(colNames, "strand")   
            dt[, 1:6]          
        } else {   
            dt[, 1:5]         
        }
    }
    dt <- checkCol4(dt) %>% checkCol5() %>% checkCol6()
    setnames(dt, colNames)
    dt
}

# for files compressed as runs of bins with the same score, expand the runs out to individual bins
# pass in an object from getTabixRangeData
expandTabixBinRuns <- function(dt, binSize, stranded = TRUE, na.value = 0, minusStrandNeg = TRUE){
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
    if(stranded && minusStrandNeg) bins[strand == "-", score := -score]
    bins[, .(strand = strand, x = start_, y = score)][order(x)]
}

# if too many bins, reduce the number of plotted XY points
# pass in an object from expandTabixBinRuns, or any similar format with no missing bins
aggregateTabixBins <- function(bins, track, coord, plotBinSize, aggFn = mean, 
                               asFractionOfMax = FALSE, maxY = NULL, limitToOne = TRUE){
    if(nrow(bins) == 0) return(data.table(strand = character(), x = integer(), y = double()))
    bins[, x := floor(x / plotBinSize) * plotBinSize + 1] # thus, x is the leftmost coordinate of the plot bin
    hasZ <- "z" %in% names(bins)
    bins <- bins[, 
        .(
            y = aggFn(as.double(y), na.rm = FALSE),
            z = if(hasZ) aggFn(as.double(z), na.rm = FALSE) else NA_real_
        ), 
        by = .(strand, x)
    ] 
    if(asFractionOfMax) {
        if(is.null(maxY) || maxY == 0) maxY <- bins[, max(y, na.rm = FALSE)]
        bins[, y := y / maxY]
        if(limitToOne) bins[, y := pmin(1, y)]
    }
    bins
}

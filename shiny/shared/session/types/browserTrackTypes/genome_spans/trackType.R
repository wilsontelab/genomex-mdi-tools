#----------------------------------------------------------------------
# genome_spans is a track type for creating many plot types from standarized BED format span files
# files do not need to actually be named ".bed", but they will functionally be BED if they can be read by this track type
# expects one or more input bed[.bgz] files, i.e., with at least columns (in this order):
#   chrom
#   start (0-indexed)
#   end   (1-indexed)
# plus variable additional columns name, score and strand, as per parseTabixBedFormat(),
# silent errors are thrown if a column required for a plot type is not found in the data.table
#----------------------------------------------------------------------
# genome_spans is not a complete track class, it is used by calling function from your track class
#----------------------------------------------------------------------

# support for span packing
getPackedSpans <- function(track, coord, spansFamily, itemsList, itemData){
    Pack_Padding_Fraction <- getTrackSetting(track, spansFamily, "Pack_Padding_Fraction", 0.025)
    padding <- coord$width * Pack_Padding_Fraction
    ends <- 0
    itemData[, y := NA_integer_]
    for(i in 1:nrow(itemData)){
        for(j in 1:length(ends)) if(itemData[i, x1] >= ends[j] + padding) itemData[i, y := j]
        if(is.na(itemData[i, y])) itemData[i, y := length(ends) + 1]
        ends[itemData[i, y]] <- itemData[i, x2]
    }
    itemsList$ymin <- 0.5
    itemsList$ymax <- length(ends) + 0.5
    list(
        itemsList = itemsList,
        itemData = itemData
    )
}
getUnpackedSpans <- function(itemsList, itemData){
    itemData[, y := 1:.N]
    itemsList$ymin <- 0.5
    itemsList$ymax <- nrow(itemData) + 0.5  
    list(
        itemsList = itemsList,
        itemData = itemData
    )
}

# track build function
build.genome_spans_track <- function(track, reference, coord, layout, dataFn, trackBuffer = NULL,
                                     spansFamily = "Spans", scoresFamily = "Scores", yAxisFamily = "Y_Axis"){

    # collect all individual bed tracks
    itemsList <- getItemsData(track, reference, coord, dataFn, parseXY = FALSE)
    if(!itemsList$hasData) return(trackInfo(track, coord, layout, "no spans in region"))
    itemNames <- sub(".bed.bgz", "", sapply(names(itemsList$d), basename))
    nItems <- length(itemNames)

    # parse the data into XY coordinates base on requested plot type
    Plot_Spans_As <- trimws(getTrackSetting(track, spansFamily, "Plot_Spans_As", "scores"))
    Stranded <- getTrackSetting(track, spansFamily, "Stranded", TRUE)
    nullStrand <- if(Stranded) "+" else "."

    # score plot types
    if(Plot_Spans_As == "heat_map"){
        Heat_Map_Bins <- getTrackSetting(track, spansFamily, "Heat_Map_Bins", 100) 
        Heat_Map_Max  <- getTrackSetting(track, spansFamily, "Heat_Map_Max", 0) 
        plotBinSize <- round(coord$width / Heat_Map_Bins, 0)
        itemsList$d <- lapply(itemsList$d, function(dt){
            if(!("score" %in% names(dt))) return(NA) else {
                if(!("strand" %in% names(dt))) dt[, strand := nullStrand]
                dt[, start := start + 1]
                bins <- expandTabixBinRuns(dt, binSize = 1, stranded = Stranded, na.value = 0, minusStrandNeg = FALSE) %>%
                        aggregateTabixBins(track, coord, plotBinSize, aggFn = mean, 
                                           asFractionOfMax = TRUE, maxY = Heat_Map_Max, limitToOne = TRUE)
                bins[, .(
                    strand = strand,
                    x1 = x,
                    x2 = x + plotBinSize - 1,
                    alpha = y
                )]
            }
        })
        missingScores <- is.na(itemsList$d)
        if(any(missingScores)) {
            return(trackInfo(track, coord, layout, paste("no score column:", paste(itemNames[missingScores], collapse = ", ")), isError = TRUE))
        }
        if(!is.null(trackBuffer)) trackBuffer[[track$id]] <- itemsList$d
        buildHeatMapTrackImage(
            track, coord, layout,
            itemsList, itemNames,
            stranded = Stranded, ylab = NULL,
            dataFamily = spansFamily
        )  

    } else if(Plot_Spans_As == "scores"){
        Score_Position <- getTrackSetting(track, spansFamily, "Score_Position", "center") 
        itemsList$d <- lapply(itemsList$d, function(dt){
            if(!("score" %in% names(dt))) return(NA) else {
                x <- data.table(
                    strand = if("strand" %in% names(dt)) dt$strand else nullStrand,                
                    x = switch(
                        Score_Position,
                        start   = dt$start + 1,
                        center  = dt$start + 1 + (dt$end - dt$start - 1) / 2,
                        end     = dt$end
                    ),
                    y = dt$score
                )    
                itemsList$ymin <<- min(itemsList$ymin, x$y, na.rm = TRUE)
                itemsList$ymax <<- max(itemsList$ymax, x$y, na.rm = TRUE)    
                x 
            }
        })
        missingScores <- is.na(itemsList$d)
        if(any(missingScores)) {
            return(trackInfo(track, coord, layout, paste("no score column:", paste(itemNames[missingScores], collapse = ", ")), isError = TRUE))
        }
        itemData <- mergeXYTrackItems(itemsList, itemNames)
        if(!is.null(trackBuffer)) trackBuffer[[track$id]] <- itemData
        buildXYTrackImage(
            track, coord, layout,
            itemsList, itemNames, itemData,
            stranded = Stranded, allowNeg = TRUE, ylab = NULL,
            dataFamily = scoresFamily, yAxisFamily = yAxisFamily
        )            

    # span plot types (including when the Y-axis is the score value)
    } else {
        # Multi_Sample <- trimws(getTrackSetting(track, spansFamily, "Multi_Sample", "admixed"))
        itemsList$hasScore <- logical()
        itemsList$d <- lapply(1:nItems, function(i){
            dt <- itemsList$d[[i]]
            hasScore <- "score" %in% names(dt)
            x <- data.table(
                source = itemNames[i],
                strand = if("strand" %in% names(dt)) dt$strand else nullStrand,                
                x1 = dt$start + 1,
                x2 = dt$end,
                y = if(!hasScore) NA else dt$score
            )    
            itemsList$ymin <<- min(itemsList$ymin, x$y, na.rm = TRUE)
            itemsList$ymax <<- max(itemsList$ymax, x$y, na.rm = TRUE)   
            itemsList$hasScore <<- c(itemsList$hasScore, hasScore)
            x 
        })
        missingScores <- !itemsList$hasScore
        if(any(missingScores) && Plot_Spans_As == "scored_spans") {
            return(trackInfo(track, coord, layout, paste("no score column:", paste(itemNames[missingScores], collapse = ", ")), isError = TRUE))
        }
        itemData <- do.call(rbind, itemsList$d)[order(x1, x2)]
        yaxt <- "s"
        if(Plot_Spans_As == "packed_spans"){
            x <- getPackedSpans(track, coord, spansFamily, itemsList, itemData)
            itemsList <- x$itemsList
            itemData <- x$itemData
            yaxt <- "n"
        } else if(Plot_Spans_As == "unpacked_spans"){
            x <- getUnpackedSpans(itemsList, itemData)
            itemsList <- x$itemsList
            itemData <- x$itemData
            yaxt <- "n"
        } # default from above is scored_spans
        if(!is.null(trackBuffer)) trackBuffer[[track$id]] <- itemData
        buildSpanTrackImage (
            track, coord, layout,
            itemsList, itemNames, itemData,
            stranded = Stranded, allowNeg = TRUE, ylab = NULL, ylim = c(itemsList$ymin, itemsList$ymax), yaxt = yaxt,
            dataFamily = spansFamily, yAxisFamily = yAxisFamily, hLines = Plot_Spans_As == "scored_spans"
        )
    }
}

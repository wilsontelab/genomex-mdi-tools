#----------------------------------------------------------------------
# binned_XY is a track type for reproducibly plotting XY scatter plots with regularly sequenced X-axis values
# expects:
#   one or more track$settings$items, yielding one (set of) overlaid traces per item
#   dataFn(track, reference, coord, itemName, item) must return binned data for one item in format of expandTabixBinRuns, aggregateTabixBins, etc.
#       i.e., as a data.table(strand = character(), x = seq(min, max, binSize), y = numeric())
#   contiguous bin runs
# supports:
#   stranded(+/-) or unstranded(.) data types
#   dynamic rebinning based on user settings
#   various plot representations (points, lines, etc.)
#   various Y axis limits types
#   aggregation of traces by difference, mean, etc.
#----------------------------------------------------------------------
# binned_XY is not a complete track class, it is used by calling function from your track class
#----------------------------------------------------------------------
# TODO: custom color palettes
#----------------------------------------------------------------------

# track build function
build.binned_XY_track <- function(track, reference, coord, layout, dataFn, 
                                  stranded = TRUE, allowNeg = FALSE, ylab = NULL,
                                  center = FALSE, binSize = NULL,
                                  dataFamily = "Data", yAxisFamily = "Y_Axis",
                                  highlightsFn = NULL, highlightsStyle = "backgroundShading", 
                                  defaultItems = NULL){

    # collect all individual items
    itemsList <- getItemsData(track, reference, coord, dataFn, defaultItems = defaultItems) 
    if(!itemsList$hasData) return(trackInfo(track, coord, layout, "no usable data to plot"))
    itemNames <- names(itemsList$d)
    nItems <- length(itemNames)
    originalItemNames <- itemNames 
    itemData <- mergeXYTrackItems(itemsList, itemNames)

    # if requested, aggregate multiple items together or plot as a change relative to the first item
    Aggregate <- getTrackSetting(track, dataFamily, "Aggregate", "none")
    if(nItems > 1 && Aggregate != "none"){
        if(Aggregate == "difference"){
            for(i in 2:nItems) itemData[[i + 2]] <- itemData[[i + 2]] - itemData[[3]]
            itemData[[3]] <- NULL
            itemNames <- paste(itemNames[2:nItems], itemNames[1], sep = " - ")
            setnames(itemData, c("strand", "x", itemNames))
            itemsList$ymin <- itemData[, min(.SD, na.rm = TRUE), .SDcols = itemNames]
            itemsList$ymax <- itemData[, max(.SD, na.rm = TRUE), .SDcols = itemNames]
            nItems <- nItems - 1
        } else {
            itemData[[3]] <- apply(itemData[, .SD, .SDcols = itemNames], 1, get(Aggregate), na.rm = TRUE)
            itemData <- itemData[, .SD, .SDcols = c("strand", "x", itemNames[1])]
            itemNames <- Aggregate   
            setnames(itemData, c("strand", "x", itemNames))
            itemsList$ymin <- min(itemData[[3]], na.rm = TRUE) 
            itemsList$ymax <- max(itemData[[3]], na.rm = TRUE) 
            nItems <- 1
        }        
    }

    # if requested, center the plot point in their bins
    if(center) itemData[, x := x + binSize / 2]

    # collect any highlight regions
    highlights <- if(is.null(highlightsFn)) NULL 
                  else getItemsData(track, reference, coord, highlightsFn, 
                                    parseXY = FALSE, optional = TRUE) 

    # make and return the plot
    buildXYTrackImage(
        track, reference, coord, layout,
        itemsList, itemNames, itemData,
        stranded = stranded, allowNeg = allowNeg, ylab = ylab,
        dataFamily = dataFamily, yAxisFamily = yAxisFamily, 
        Aggregate = Aggregate, originalItemNames = originalItemNames,
        highlights = highlights, highlightsStyle = highlightsStyle
    )
}

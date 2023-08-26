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
                                  center = FALSE, binSize = NULL){

    # collect all individual items
    itemsList <- getItemsData(track, reference, coord, dataFn) 
    if(!itemsList$hasData) return(trackInfo(track, coord, layout, "no usable data to plot"))
    itemNames <- names(itemsList$d)
    nItems <- length(itemNames)
    originalItemNames <- itemNames
    originialNItems <- nItems  

    # merge samples up to a single table to get all possible x values
    itemData <- Reduce(function(dt1, dt2) merge(dt1, dt2, by = c("strand","x"), all=TRUE), itemsList$d)
    setnames(itemData, c("strand", "x", itemNames)) 

    # if requested, aggregate multiple items together or plot as a change relative to the first item
    Aggregate <- getTrackSetting(track, "Data", "Aggregate", "none")
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

    # set the plot frame
    padding <- padding(track, layout)
    height <- height(track, 2) + padding$total # or set a known, fixed height in inches
    ylab <- ylab(track,
        if(is.null(ylab)) ""
        else if(is.function(ylab)) ylab()
        else ylab
    )

    # set the dynamic Y-axis
    Min_Y         <- trimws(getTrackSetting(track, "Y_Axis", "Min_Y", ""))
    Max_Y         <- trimws(getTrackSetting(track, "Y_Axis", "Max_Y", ""))
    Force_To_Zero <- getTrackSetting(track, "Y_Axis", "Force_To_Zero", TRUE)
    Symmetric     <- getTrackSetting(track, "Y_Axis", "Symmetric", TRUE)
    if(is.infinite(itemsList$ymin)) itemsList$ymin <- 0
    if(is.infinite(itemsList$ymax)) itemsList$ymax <- 0.01
    ylim <- if(stranded || allowNeg) {
        extremeValue <- max(abs(c(
            itemsList$ymin, 
            itemsList$ymax
        )))
        if(Symmetric) c(
            if(Min_Y == "") -extremeValue else as.numeric(Min_Y), 
            if(Max_Y == "")  extremeValue else as.numeric(Max_Y)
        ) else if(Force_To_Zero) c(
            if(Min_Y == "") min(0, itemsList$ymin) else as.numeric(Min_Y), 
            if(Max_Y == "") max(0, itemsList$ymax) else as.numeric(Max_Y)
        ) else c(
            if(Min_Y == "") itemsList$ymin else as.numeric(Min_Y), 
            if(Max_Y == "") itemsList$ymax else as.numeric(Max_Y)
        )                
    } else {
        if(Force_To_Zero) c(
            if(Min_Y == "") 0              else as.numeric(Min_Y), 
            if(Max_Y == "") itemsList$ymax else as.numeric(Max_Y)
        ) else c(
            if(Min_Y == "") itemsList$ymin else as.numeric(Min_Y), 
            if(Max_Y == "") itemsList$ymax else as.numeric(Max_Y)
        )  
    }

    # if requested, center the plot point in their bins
    if(center) itemData[, x := x + binSize / 2]

    # make the plot
    mai <- NULL    
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = ylab, # yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        # add horizontal rules if requested  
        zeroLine(track)  
        hLines(track, ylim)   

        # randomize point order to avoid overplotting
        I <- 1:nItems
        palette <- CONSTANTS$palettes[[getTrackSetting(track, "Data", "Color_Palette", "plotly")]]        
        if(getTrackSetting(track, "Data", "Plot_As", "points") == "points"){
            dd <- do.call(rbind, lapply(I, function(i) itemData[, 
                .(strand = strand, x = x, y = .SD[[itemNames[i]]], color = palette[[i]])
            ]))[order(sample(.N))]
            if(stranded) for(strand_ in c("+", "-")){
                ddd <- dd[strand == strand_]
                plotXY(track, ddd, color = ddd$color, family = "Data")
            } else {
                plotXY(track, dd,  color = dd$color,  family = "Data")
            }

        # other track types are inherently overplotted    
        # order in items list determines plot order         
        } else {
            for(i in I){ 
                dd <- itemData[, .(strand = strand, x = x, y = .SD[[itemNames[i]]])]
                if(stranded) for(strand_ in c("+", "-")){
                    plotXY(track, dd[strand == strand_], color = palette[[i]], family = "Data")
                } else {
                    plotXY(track, dd, color = palette[[i]], family = "Data")
                }
            }            
        }

        # add a legend
        if(Aggregate == "none"){
            legend <- itemNames
            colors <- unlist(palette[I])
        } else if(Aggregate == "difference"){
            legend <- originalItemNames
            legend[1] <- paste("x - ", legend[1])
            colors <- c(NA, unlist(palette[I]))
        } else {
            legend <- c(Aggregate, originalItemNames)
            colors <- c(palette[[1]], rep(NA, originialNItems))
        }
        trackLegend(track, coord, ylim, legend = legend, pch = 19, cex = 1, col = colors)
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

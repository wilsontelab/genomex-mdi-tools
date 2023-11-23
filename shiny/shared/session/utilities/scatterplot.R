# generic support function for XY browser track types
# expects:
#   itemsList as returned by getItemsData()
#   itemData as returned by mergeXYTrackItems()
#   associated itemNames
#   all parsing of XY values has been done by caller
# supports:
#   stranded(+/-) or unstranded(.) data types
#   various plot representations (points, lines, etc.)
#   various Y axis limits types

# merge XY data samples up to a single table to get all possible x values
# thus parses itemsList to itemData, with one named column per y data series
mergeXYTrackItems <- function(itemsList, itemNames){
    itemData <- Reduce(function(dt1, dt2) merge(dt1, dt2, by = c("strand","x"), all=TRUE), itemsList$d)
    setnames(itemData, c("strand", "x", itemNames)) 
    itemData
}

# support functions for buildXYTrackImage and other tracks builders with similar needs
getDynamicYLim <- function(track, yAxisFamily, itemsList, stranded, allowNeg){
    Min_Y         <- trimws(getTrackSetting(track, yAxisFamily, "Min_Y", ""))
    Max_Y         <- trimws(getTrackSetting(track, yAxisFamily, "Max_Y", ""))
    Force_To_Zero <- getTrackSetting(track, yAxisFamily, "Force_To_Zero", TRUE)
    Symmetric     <- getTrackSetting(track, yAxisFamily, "Symmetric", TRUE)
    if(is.infinite(itemsList$ymin)) itemsList$ymin <- 0
    if(is.infinite(itemsList$ymax)) itemsList$ymax <- 0.01
    if(stranded || allowNeg) {
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
}

# parse standardized input to create consistently formatted XY plots of various types
buildXYTrackImage <- function(track, reference, coord, layout,
                               itemsList, itemNames, itemData,
                               stranded = TRUE, allowNeg = FALSE, ylab = NULL,
                               dataFamily = "Data", yAxisFamily = "Y_Axis", 
                               Aggregate = "none", originalItemNames = NULL,
                               highlights = NULL, highlightsStyle = "backgroundShading"){
    nItems <- length(itemNames)

    # set the plot frame
    padding <- padding(track, layout)
    height <- height(track, 2) + padding$total # or set a known, fixed height in inches
    ylab <- ylab(track,
        if(is.null(ylab)) {
            if(nItems == 1) itemNames else ""
        }
        else if(is.function(ylab)) ylab()
        else ylab
    )

    # set the dynamic Y-axis
    ylim <- getDynamicYLim(track, yAxisFamily, itemsList, stranded, allowNeg)

    # make the plot
    mai <- NULL    
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = ylab, # yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        # add rules if requested/required
        zeroLine(track)  
        hLines(track, ylim)   
        chromLines(track, reference, coord)

        # add any region highlights underneath the plot points
        if(!is.null(highlights)) for(sample in names(highlights$d)){
            d <- highlights$d[[sample]]
            if(nrow(d) == 0) next
            rect(d$x1, ylim[1], d$x2, ylim[2], 
                 col = if("color" %in% names(d)) d$color else "grey50", border = NA)
        }

# List of 4
#  $ d      :List of 2
#   ..$ NA12878:Classes ‘data.table’ and 'data.frame':    0 obs. of  3 variables:
#   .. ..$ x1   : int(0)
#   .. ..$ x2   : int(0)
#   .. ..$ color: logi(0)
#   .. ..- attr(*, ".internal.selfref")=<externalptr>
#   ..$ HCT116 :Classes ‘data.table’ and 'data.frame':    2 obs. of  3 variables: 
#   .. ..$ x1   : int [1:2] 69810000 70510000
#   .. ..$ x2   : int [1:2] 70410000 81910000
#   .. ..$ color: chr [1:2] "blue" NA
#   .. ..- attr(*, ".internal.selfref")=<externalptr>
#  $ ymin   : logi NA
#  $ ymax   : logi NA
#  $ hasData: logi TRUE


        # randomize point order to avoid overplotting
        I <- 1:nItems
        palette <- CONSTANTS$palettes[[getTrackSetting(track, dataFamily, "Color_Palette", "plotly")]]        
        if(getTrackSetting(track, dataFamily, "Plot_As", "points") == "points"){
            dd <- do.call(rbind, lapply(I, function(i) itemData[, 
                .(strand = strand, x = x, y = .SD[[itemNames[i]]], color = palette[[i]])
            ]))[order(sample(.N))]
            if(stranded) for(strand_ in c("+", "-")){
                ddd <- dd[strand == strand_]
                plotXY(track, ddd, color = ddd$color, family = dataFamily)
            } else {
                plotXY(track, dd,  color = dd$color,  family = dataFamily)
            }

        # other track types are inherently overplotted    
        # order in items list determines plot order         
        } else {
            for(i in I){ 
                dd <- itemData[, .(strand = strand, x = x, y = .SD[[itemNames[i]]])]
                if(stranded) for(strand_ in c("+", "-")){
                    plotXY(track, dd[strand == strand_], color = palette[[i]], family = dataFamily)
                } else {
                    plotXY(track, dd, color = palette[[i]], family = dataFamily)
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
            colors <- c(palette[[1]], rep(NA, length(originalItemNames)))
        }
        trackLegend(track, coord, ylim, legend = legend, pch = 19, cex = 1, col = colors)
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}

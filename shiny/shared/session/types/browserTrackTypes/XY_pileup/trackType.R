#----------------------------------------------------------------------
# XY_pileup is a track type for reproducibly plotting pileup tracks as stacked barplots
# the prototypical case is read_pileup but any data type with categories per span can be used
# expects:
#   a list of one or more pileup data sets returned by dataFn(), each plotted as one horizontal subTrack, in format:
#       list(list(ylab = character, pileup = data.table(start, end, category1, category2, ...)))
#   allValues parameter as a character vector of all allowed categories, even if missing from table
#       allValues dictates the order of the bar plot stacking
#   colorPalette, a named list with a plot color for every possible category in allValues
# where:
#   spans (start:end) MUST be contiguous throughout the browser window to maintain vertical register
#   spans do NOT have to be equal width; bar widths are adjusted to reflect span widths
#   start and end are both 1-referenced (not BED style)
#   category values may or may not add up to the same value in each span, as determined by the caller
# note:
#   unlike most other browser tracks, XY_pileup, via its call to barplot, uses an arbitrary X axis
#   that always fills the available horizontal space (hence the need for contiguous spans)
#----------------------------------------------------------------------
# XY_pileup is not a complete track class, it is used by calling function from your track class
#----------------------------------------------------------------------

# track build function
build.XY_pileup_track <- function(track, reference, coord, layout, 
                                  dataFn, allValues, colorPalette){

    # check for a valid plot window
    Max_Width <- getTrackSetting(track, "Pileup", "Max_Width", 1000)
    if(!isTruthy(coord$width <= Max_Width)) return(trackInfo(track, coord, layout, "window too wide to plot pileup"))

    # parse the input data
    pileups <- dataFn()
    if(!allAreTruthy(pileups, length(pileups) > 0)) 
        return(trackInfo(track, coord, layout, "no usable data to plot"))
    pileups <- lapply(pileups, function(pileup){
        if(nrow(pileup$pileup) == 0) return(pileup)
        pileup$pileup <- pileup$pileup[end >= coord$start & start <= coord$end] # flush the ends out, masking the first and last as unreliable
        pileup$pileup[1, ":="(start = coord$start)]
        pileup$pileup[nrow(pileup$pileup), ":="(end = coord$end)]
        usedValues <- allValues[allValues %in% names(pileup$pileup)]
        c(pileup, list(
            usedValues = usedValues,
            maxY = pileup$pileup[, max(rowSums(.SD, na.rm = TRUE)), .SDcols = usedValues]
        ))      
    })

    # set the layout
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    ylim <- NULL # set below 

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        layout(matrix(1:length(pileups), ncol = 1, byrow = TRUE))
        par(xpd = TRUE)        
        for(pileup in pileups){
            if(nrow(pileup$pileup) == 0){
                ylim <<- c(0, 1)
                plot(0, 0, type = "n", bty = "n",
                    xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
                    ylim = ylim,  ylab = "", yaxt = "n",
                    xaxs = "i", yaxs = "i")
                trackNoData(coord, ylim, paste(pileup$ylab, "has no pileup data in window"))
            } else {
                ylim <<- c(0, max(sapply(pileups, function(x) x$maxY)))
                colors <- unlist(colorPalette[pileup$usedValues])
                barplot(
                    height = t(as.matrix(pileup$pileup[, .SD, .SDcols = pileup$usedValues])), 
                    width = pileup$pileup[, end - start + 1],
                    col = colors,
                    space = 0, border = NA, 
                    xlab = "", xaxt = "n",
                    ylim = ylim, ylab = pileup$ylab, 
                    bty = "n",
                    xaxs = "i", yaxs = "i",
                    legend.text = TRUE, args.legend = list(
                        bty = "n", 
                        x = coord$width * 1.1,
                        y = ylim[2],
                        border = rev(colors)
                    )
                )                
            }
        }
        par(xpd = FALSE)
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}

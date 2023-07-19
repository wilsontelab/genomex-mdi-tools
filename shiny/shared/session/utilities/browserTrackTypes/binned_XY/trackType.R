#----------------------------------------------------------------------
# binned_XY is a track type for reproducibly plotting XY scatter plots with regularly sequenced X-axis values
# expects:
#   one or more track$settings$items, yielding one (set of) overlaid traces per item
#   dataFn(track, reference, coord, itemName, item) must return binned data in format of expandTabixBinRuns, aggregateTabixBins, etc.
#       i.e., as data.table(strand = character(), x = seq(min, max, binSize), y = numeric())
# supports:
#   stranded(+/-) or unstranded(.) data types
#   dynamic rebinning based on user settings
#   various plot representations (points, lines, etc.)
#----------------------------------------------------------------------
# binned_XY is not a complete track class, it is used by calling function from your track class
#----------------------------------------------------------------------
# TODO: custom color palettes
#----------------------------------------------------------------------

# track build function
build.binned_XY_track <- function(track, reference, coord, layout, 
                                  dataFn, stranded = TRUE, allowNeg = TRUE, ylab = ""){
    items <- getItemsData(track, reference, coord, dataFn, stranded = stranded)
    padding <- padding(track, layout)
    height <- height(track, 2) + padding$total # or set a known, fixed height in inches
    ylim <- ylim(track, if(stranded) c(-items$ymax, items$ymax) else c(items$ymin, items$ymax))
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

        # plot traces on top of each other by item
        # order in items list determines plot order 
        I <- 1:length(items$d)
        palette <- CONSTANTS$plotlyColors
        for(i in I){
            if(stranded) for(strand_ in c("+", "-")){
                plotXY(track, items$d[[i]][strand == strand_], palette[[i]])
            } else {
                plotXY(track, items$d[[i]], palette[[i]])
            }
        }

        # add a legend
        trackLegend(track, coord, ylim, 
                    legend = names(items$d), pch = 19, cex = 1, col = unlist(palette[I]))
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

# parse standardized input to create consistently heatmap tracks
# expect elements of itemsList$d to be data.table(strand, x1 = coord, x2 = coord, alpha = alpha)
buildHeatMapTrackImage <- function(track, coord, layout,
                                   itemsList, itemNames,
                                   stranded = TRUE, ylab = NULL, 
                                   dataFamily = "Data"){
    nItems <- length(itemNames)
    Heat_Map_Exponent  <- getTrackSetting(track, dataFamily, "Heat_Map_Exponent", 1) 

    # set the plot frame
    padding <- padding(track, layout)
    height <- height(track, 2) + padding$total # or set a known, fixed height in inches

    # set the dynamic Y-axis
    ylim <- c(0.5, nItems * (stranded + 1) + 0.5)

    # make the plot
    mai <- NULL    
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 
        I <- 1:nItems
        colorPalette <- getTrackSetting(track, dataFamily, "Spans_Color_Palette", "plotly")
        isGreyscale <- colorPalette == "greyscale"
        palettes <- CONSTANTS$palettes[[colorPalette]]        
        for(i in I){ 
            color <- if(isGreyscale) "#000000" else palettes[[i]]
            if(stranded) for(strand_ in c("+", "-")){
                dd <- itemsList$d[[i]][strand == strand_]
                dd[, y := ylim[2] + 0.5 - i * ((strand_ == "-") + 1) ] # plot 1st sample on top, + on top of - strands
                plotHeatMap(track, dd, color = color, exponent = Heat_Map_Exponent, family = dataFamily)
            } else {
                dd <- itemsList$d[[i]]
                dd[, y := ylim[2] + 0.5 - i ]
                plotHeatMap(track, dd, color = color, exponent = Heat_Map_Exponent, family = dataFamily)
            }
        }            

        # add a legend
        if(!isGreyscale){
            legend <- itemNames
            colors <- unlist(palettes[I])
            trackLegend(track, coord, ylim, legend = legend, pch = 19, cex = 1, col = colors)
        }
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}

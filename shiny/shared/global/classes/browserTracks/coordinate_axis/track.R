#----------------------------------------------------------------------
# coordinate_axis trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_coordinate_axisTrack <- function(trackId) {
    list(  
        click = TRUE,
        hover = FALSE,
        brush = FALSE,
        items = FALSE
    )
}

# build method for the S3 class
build.coordinate_axisTrack <- function(track, reference, coord, layout){
    req(objectHasData(reference$genome))
    padding <- padding(track, layout)
    height <- 3.3 / layout$linesPerInch + padding$total
    unit <- parseUnit(scaleUnit(track), coord$end)
    xlab <- paste0(reference$genome$genome, " ", coord$chromosome, " (", unit$unit, ")")   
    axis <- getBrowserTrackSetting(track, "Plot_Options", "Axis_Orientation", default = "bottom")        
    lwd <- lwd(track)
    xlim <- coord$range / unit$multiplier
    ylim <- c(0, 1)
    bty <- switch(axis, top = "7", bottom = "l", "n")
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mar <- if(axis == "top") list(top = 3.1, bottom = 0.1)
               else if(axis == "bottom") list(top = 0.1, bottom = 3.1)
               else list(top = 1.25, bottom = 1.25)
        mai <<- setMdiTrackMai(layout, padding, mar = mar)
        plot(0, 0, type = "n", bty = bty,
            xlim = xlim, xlab = "", xaxt = "n",
            ylim = ylim, ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") 
        if(axis == "top"){
            axis(3, lwd = lwd)
            mtext(xlab, side = 3, line = 2, cex = 1.1)
            box(typ = "o", col = "red")
        } else if(axis == "bottom") {
            axis(1, lwd = lwd)
            mtext(xlab, side = 1, line = 2, cex = 1.1)
        } else {
            at <- axis(3, labels = FALSE, tick = FALSE)
            lines(xlim, rep(0.5, 2), lwd = lwd)        
            segments(at, 0, at, 1, lwd = lwd)
        }
    })
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click or track$hover is TRUE, above
click.coordinate_axisTrack <- function(track, click){
    app$browser$center(click$coord$x)
}

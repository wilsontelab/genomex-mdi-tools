#----------------------------------------------------------------------
# support functions that help browser tracks report in-track text feedback messages
#   trackInfo()   = a non-failure message in place of the normal track, in black font
#   trackError()  = a failure message in place of the normal track, in red font
#   trackNoData() = add a non-failure message into the middle of the normal but empty track, in black font
#----------------------------------------------------------------------
# usage from within build.track functions:
#   if(!allAreTruthy(...)) return( trackInfo( track, coord, layout, "my info message") )
#   if(!allAreTruthy(...)) return( trackError(track, coord, layout, "my error message") )
#   plot(...)
#   if(trackHasNoData) trackNoData(coord, ylim, "my noData message") else <plot the track data>
#   <return track as per normal>
#----------------------------------------------------------------------
trackInfo <- function(track, coord, layout, message, isError = FALSE, cex = 1.1){
    padding <- padding(track, layout)
    height <- 2.2 / layout$linesPerInch + padding$total
    ylim <- c(0, 1)
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 2.1, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n",
            ylim = ylim,  ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") 
        mtext(
            paste0(track$type, ": ", message), 
            side = 3, 
            line = 0.5, 
            outer = FALSE, 
            at = NA,
            adj = NA, 
            padj = NA, 
            cex = cex,
            col = if(isError) "red3" else "grey10"
        )
    })
    list(ylim = ylim, mai = mai, image = image)
}
trackError <- function(track, coord, layout, message, cex = 1.1){
    trackInfo(track, coord, layout, message, isError = TRUE, cex = cex)
}
trackNoData <- function(coord, ylim, message, y = NULL, isError = FALSE, cex = 1.1){
    if(is.null(y)) y <- ylim[1] + diff(ylim) / 2
    text(
        coord$start + coord$width / 2,
        y,
        message, 
        cex = cex,
        col = if(isError) "red3" else "grey10"
    )
}

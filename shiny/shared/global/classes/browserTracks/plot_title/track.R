#----------------------------------------------------------------------
# plot_title trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_plot_titleTrack <- function(trackId) {
    list(
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = FALSE
    )
}

# build method for the S3 class
build.plot_titleTrack <- function(track, reference, coord, layout){
    padding <- padding(track, layout)
    height <- 2.2 / layout$linesPerInch + padding$total
    ylim <- c(0, 1)
    mai <- NULL
    title <- getBrowserTrackSetting(track, "Track", "Plot_Title", "auto")
    if(title == "auto") title <- paste(reference$genome$description, coord$region)
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 2.1, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n",
            ylim = ylim,  ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") 
        mtext(
            title, 
            side = 3, 
            line = 0.5, 
            outer = FALSE, 
            at = NA,
            adj = NA, 
            padj = NA, 
            cex = 1.2
        )
        if(getBrowserTrackSetting(track, "Track", "Show_Line", FALSE))
            lines(c(coord$start, coord$end), c(0, 0), col = "black", lwd = 2)
    })
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

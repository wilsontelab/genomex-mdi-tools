#----------------------------------------------------------------------
# scale_bar trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_scale_barTrack <- function(trackId) {
    list(
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = FALSE
    )
}

# build method for the S3 class
build.scale_barTrack <- function(track, reference, coord, layout){
    padding <- padding(track, layout)
    height <- 1.5 / layout$linesPerInch + padding$total
    ylim <- c(0, 1)
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0.1, bottom = 0.1))
        plot(0, 0, typ = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n",
            ylim = ylim, ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i")
        strw <- strwidth(" 888 Mb  ") * 1.1 / 2         
        exponent <- 9
        while(10^exponent > coord$width - strw) exponent <- exponent - 1
        barWidth <- 10^exponent
        unit <- parseUnit(scaleUnit(track), barWidth)
        lab <- paste(" ", barWidth / unit$multiplier, " ", unit$unit, "  ", sep = "")
        strw <- strwidth(lab) * 1.1
        rect(coord$start, 0.6, coord$start + barWidth, 0.4, col = "black")

        # margin text label
        tmp <- par('xpd')
        par(xpd = TRUE, cex = 1.15)
        text(coord$start - strw / 2, 0.5, lab)
        par(xpd = tmp, cex = 1)
    })
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

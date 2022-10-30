#----------------------------------------------------------------------
# generic support functions for drawing browserTracks
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# MDI (as opposed to UCSC) track image assembly
#----------------------------------------------------------------------
mdiTrackImage <- function(layout, height, plotFn, ...){
    # NB: do not use magick::image_graph; it is buggy, especially on Linux
    pngFile <- tempfile("trackBrowser", fileext = ".png")
    png(
        filename = pngFile,
        width = layout$width, 
        height = height * layout$dpi, 
        units = "px", 
        pointsize = layout$pointsize,
        bg = "white",  
        res = layout$dpi, 
        type = "cairo"
    )
    plotFn(...)
    dev.off()
    image <- magick::image_read(pngFile)
    unlink(pngFile)
    image
}
setMdiTrackMai <- function(layout, padding, mai = NULL, mar = NULL){ # layout set by browser, mai or mar set by track
    if(is.null(mai)) mai <- list( # track allowed to provide either mai or mar
        top    = mar$top    / layout$linesPerInch,
        bottom = mar$bottom / layout$linesPerInch
    )
    mai <- c(
        mai$bottom + padding$bottom,
        layout$mai$left,
        mai$top + padding$top,
        layout$mai$right
    )
    par(mai = mai)
    mai
}

#----------------------------------------------------------------------
# convert browser inputs to graphic options
# helper methods for the browserInput S3 class
#----------------------------------------------------------------------
coordinates.browserInput <- function(input){
    # TODO: implement conversion of jumpTo region or gene to coordinates
    start <- trimws(input$start)
    end   <- trimws(input$end)
    if(length(start) == 0 || start == "") start <- 1
    if(length(end)   == 0 || end   == "") end   <- as.integer(start) + 10000 
    start <- as.integer(start)
    end   <- as.integer(end)
    getBrowserCoord(input$chromosome, start, end)
}

#----------------------------------------------------------------------
# override browser input values and settings for tracks where adjustsWidth == TRUE
#----------------------------------------------------------------------
getBrowserCoord <- function(chrom, start, end){
    list(
        chromosome = chrom,
        start  = start, 
        end    = end,
        width  = end - start + 1,
        range  = c(start, end),
        region = paste0(chrom, ":", start, "-", end)
    )
}
adjustLayoutWidth <- function(layout, plotWidth){
    layout$plotWidth <- plotWidth / layout$dpi
    layout$browserWidth <- layout$plotWidth + layout$mai$left + layout$mai$right
    layout$width <- layout$browserWidth * layout$dpi
    layout
}

#----------------------------------------------------------------------
# convert track settings to graphic options
# typically, give priority to user setting with fallback to track default
# helper methods for the browserTrack S3 class
#----------------------------------------------------------------------
getBrowserTrackSetting <- function(track, optionFamily, option, default = NULL){
    user <- track$settings[[optionFamily]]()[[option]]$value # track developer must ensure option is offered
    if(is.null(user)) return(default)
    if(typeof(user) == "character"){
        user <- trimws(user)
        if(user == "" || user == "auto") return(default)
    }
    user
}

# Track_Options options family
padding.browserTrack <- function(track, layout){
    top    <- getInches(track$settings$get("Track_Options", "Top_Padding"),    layout$lengthUnit)
    bottom <- getInches(track$settings$get("Track_Options", "Bottom_Padding"), layout$lengthUnit)
    list(
        top = top,
        bottom = bottom,
        total = top + bottom
    )
}
height.browserTrack <- function(track, default){
    getBrowserTrackSetting(track, "Track_Options", "Height", default)
}
ylab.browserTrack <- function(track, default = ""){
    getBrowserTrackSetting(track, "Track_Options", "Y_Axis_Label", default)
}
ylim.browserTrack <- function(track, y){
    user <- getBrowserTrackSetting(track, "Track_Options", "Y_Limits")
    if(is.null(user)) return(paddedRange(y))
    user <- gsub('\\s', '', user)
    user <- gsub('to', '-', user)
    user <- gsub('::', ':', user)
    user <- as.numeric(strsplit(user, '[-:,]')[[1]])
    if(is.na(user[2])) user[2] <- -user[1]
    sort(user)
}
bty.browserTrack <- function(track, default){
    user <- getBrowserTrackSetting(track, "Track_Options", "Bounding_Box")
    if(is.null(user)) return(default)
    if(user) "o" else "n"
}

# Plot_Options options family
scaleUnit.browserTrack <- function(track, default = "auto"){
    getBrowserTrackSetting(track, "Plot_Options", "Scale_Unit", default)
}
typ.browserTrack <- function(track, default){
    user <- getBrowserTrackSetting(track, "Plot_Options", "Plot_Type")
    if(is.null(user)) return(default)
    switch(
        user,
        Points = "p",
        Lines = "l",
        Both = "b",
        Histogram = "h"
        # TODO: implement special handling of area, histogram, 
    )
}
pch.browserTrack <- function(track, default){
    user <- getBrowserTrackSetting(track, "Plot_Options", "Point_Symbol")
    if(is.null(user)) return(default)
    switch(
        user,
        "Open Circles" = 1,
        "Filled Circles" = 19,
        "Open Squares" = 0,
        "Filled Squares" = 15
    )
}
lwd.browserTrack <- function(track, default){
    user <- getBrowserTrackSetting(track, "Plot_Options", "Line_Width")
    if(is.null(user)) return(default)
    user
}
cex.browserTrack <- function(track, default){
    user <- getBrowserTrackSetting(track, "Plot_Options", "Point_Size")
    if(is.null(user)) return(default)
    user
}
col.browserTrack <- function(track, default){
    user <- getBrowserTrackSetting(track, "Plot_Options", "Color")
    if(is.null(user)) return(default)
    CONSTANTS$plotlyColor[[user]]
}

#----------------------------------------------------------------------
# generic track plotting functions
# these are not S3 methods but are named similarly for clarity
#----------------------------------------------------------------------

# plot XY data tracks
plotXY.browserTrack <- function(
    track,
    input,
    x, 
    y, 
    xaxt = "n",  # most tracks do not show their X axis
    xlab = "",   # or X axis label, but track can override this behavior
    ylab = NULL,
    typ  = NULL,
    bty  = NULL,
    pch  = NULL,
    lwd  = NULL,
    cex  = NULL,
    col  = NULL,
    ... # additional arguments passed to plot()
){
    if(is.null(ylab)) ylab <- ylab(track, "")
    if(is.null(typ))  typ  <-  typ(track, "p")
    if(is.null(bty))  bty  <-  bty(track, "n")
    if(is.null(pch))  pch  <-  pch(track, 16)
    if(is.null(lwd))  lwd  <-  lwd(track, 1)
    if(is.null(cex))  cex  <-  cex(track, 1)
    if(is.null(col))  col  <-  col(track, "black")
    plot(
        x = x,
        y = y,
        xaxt = xaxt,
        xlab = xlab,
        ylab = ylab,
        typ = typ,
        bty = bty,
        pch = pch,
        lwd = lwd,
        cex = cex,
        col = col,
        xaxs = "i",
        yaxs = "i",
        ...
    )
}

#----------------------------------------------------------------------
# generic support functions for drawing browserTracks
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# MDI (as opposed to UCSC) track image assembly
#----------------------------------------------------------------------
mdiTrackImage <- function(layout, height, plotFn, ...){
    image <- magick::image_graph( 
        width = layout$width,
        height = height * layout$dpi, # provided by track in inches
        bg = "white", 
        pointsize = layout$pointsize, 
        res = layout$dpi, 
        clip = FALSE, 
        antialias = TRUE 
    )
    plotFn(...)
    dev.off()
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
    list(
        chromosome = input$chromosome,
        start  = start, 
        end    = end,
        width  = end - start + 1,
        range  = c(start, end),
        region = paste0(input$chromosome, ":", start, "-", end)
    )
}

#----------------------------------------------------------------------
# convert track settings to graphic options
# typically, give priority to user setting with fallback to track default
# helper methods for the trackSettings S3 class
#----------------------------------------------------------------------
getBrowserTrackSetting <- function(settings, optionFamily, option, default = NULL){
    user <- settings[[optionFamily]]()[[option]]$value # track developer must ensure option is offered
    if(is.null(user)) return(default)
    if(typeof(user) == "character"){
        user <- trimws(user)
        if(user == "" || user == "auto") return(default)
    }
    user
}

# Track_Options options family
padding.trackSettings <- function(settings, layout){
    top    <- getInches(settings$get("Track_Options", "Top_Padding"),    layout$lengthUnit)
    bottom <- getInches(settings$get("Track_Options", "Bottom_Padding"), layout$lengthUnit)
    list(
        top = top,
        bottom = bottom,
        total = top + bottom
    )
}
height.trackSettings <- function(settings, default){
    getBrowserTrackSetting(settings, "Track_Options", "Height", default)
}
ylab.trackSettings <- function(settings, default = ""){
    getBrowserTrackSetting(settings, "Track_Options", "Y_Axis_Label", default)
}
ylim.trackSettings <- function(settings, y){
    user <- getBrowserTrackSetting(settings, "Track_Options", "Y_Limits")
    if(is.null(user)) return(paddedRange(y))
    user <- gsub('\\s', '', user)
    user <- gsub('to', '-', user)
    user <- gsub('::', ':', user)
    user <- as.numeric(strsplit(user, '[-:,]')[[1]])
    if(is.na(user[2])) user[2] <- -user[1]
    sort(user)
}
bty.trackSettings <- function(settings, default){
    user <- getBrowserTrackSetting(settings, "Track_Options", "Bounding_Box")
    if(is.null(user)) return(default)
    if(user) "o" else "n"
}

# Plot_Options options family
scaleUnit.trackSettings <- function(settings, default = "auto"){
    getBrowserTrackSetting(settings, "Plot_Options", "Scale_Unit", default)
}
typ.trackSettings <- function(settings, default){
    user <- getBrowserTrackSetting(settings, "Plot_Options", "Plot_Type")
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
pch.trackSettings <- function(settings, default){
    user <- getBrowserTrackSetting(settings, "Plot_Options", "Point_Symbol")
    if(is.null(user)) return(default)
    switch(
        user,
        "Open Circles" = 1,
        "Filled Circles" = 19,
        "Open Squares" = 0,
        "Filled Squares" = 15
    )
}
lwd.trackSettings <- function(settings, default){
    user <- getBrowserTrackSetting(settings, "Plot_Options", "Line_Width")
    if(is.null(user)) return(default)
    user
}
cex.trackSettings <- function(settings, default){
    user <- getBrowserTrackSetting(settings, "Plot_Options", "Point_Size")
    if(is.null(user)) return(default)
    user
}
col.trackSettings <- function(settings, default){
    user <- getBrowserTrackSetting(settings, "Plot_Options", "Color")
    if(is.null(user)) return(default)
    CONSTANTS$plotlyColor[[user]]
}

#----------------------------------------------------------------------
# generic track plotting functions
# these are not S3 methods but are named similarly for clarity
#----------------------------------------------------------------------

# plot XY data tracks
plotXY.trackSettings <- function(
    settings,
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
    if(is.null(ylab)) ylab <- ylab(settings, "")
    if(is.null(typ)) typ <- typ(settings, "p")
    if(is.null(bty)) bty <- bty(settings, "n")
    if(is.null(pch)) pch <- pch(settings, 16)
    if(is.null(lwd)) lwd <- lwd(settings, 1)
    if(is.null(cex)) cex <- cex(settings, 1)
    if(is.null(col)) col <- col(settings, "black")
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

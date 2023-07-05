#----------------------------------------------------------------------
# generic support functions for drawing browserTracks
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# MDI (as opposed to UCSC) track image assembly
#----------------------------------------------------------------------
mdiTrackImage <- function(layout, height, plotFn, message = NULL, ...){ # for track that generate plots
    # NB: do not use magick::image_graph; it is buggy, especially on Linux
    startSpinner(session, message = message)
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
    stopSpinner(session)
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
pngToMdiTrackImage <- function( # for tracks that generate images, not plots
    pngFile, 
    layout, 
    verticalPadding = 0L, # in pixels
    hasLeft  = FALSE, # set to TRUE if pngFile already has a region for left or right margin labels
    hasRight = FALSE
){
    image <- magick::image_read(pngFile)
    info  <- magick::image_info(image)
    heightOut <- info$height + 2 * verticalPadding
    missingLeftPixels <- if(hasLeft) 0 else layout$mai$left * layout$dpi
    image <- magick::image_extent( # add space for the right legend region and top and bottom padding
        image, 
        geometry = paste0(layout$width - missingLeftPixels, "x", heightOut), 
        gravity = "West", 
        color = "white"
    )
    if(!hasLeft) image <- magick::image_extent( # add space for the right legend region and top and bottom padding
        image, 
        geometry = paste0(layout$width, "x", heightOut), 
        gravity = "East", 
        color = "white"
    )
    image
}

#----------------------------------------------------------------------
# automated track naming
#----------------------------------------------------------------------
getTrackDisplayName <- function(track){
    trackName <- track$settings$get("Track_Options", "Track_Name")
    if(
        is.null(trackName) || 
        trackName == "" || 
        trackName == "auto" || 
        trackName == track$type
    ){
        if(
            isTruthy(track$items) && 
            !is.null(track$settings$items)
        ){
            items <- track$settings$items() # list of lists
            if(length(items) == 1) {
                item <- items[[1]]
                itemFields <- names(item)
                allowedFields <- itemFields %in% c("name","Name","Sample_ID","sample","Sample","track","Track") # could add more recognized name defaults here
                if(any(allowedFields)) {
                    trackName <- item[[itemFields[which(allowedFields)[1]]]]
                } else {
                    trackName <- track$type 
                }
            } else {
                trackName <- "multi-sample"
            }
        } else {
            trackName <- track$type 
        }
    }
    trackName  
}

#----------------------------------------------------------------------
# convert browser inputs to graphic options
# helper methods for the browserInput S3 class
#----------------------------------------------------------------------
coordinates.browserInput <- function(input){
    start <- trimws(input$start)
    end   <- trimws(input$end)
    req(start, end)
    # if(length(start) == 0 || start == "") start <- 1
    # if(length(end)   == 0 || end   == "") end   <- as.integer64(start) + 10000 
    start <- as.integer64(start)
    end   <- as.integer64(end)
    req(start < end, end > 1)
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
        range  = as.numeric(c(start, end)), # using as numeric allows the range to be used in plot(xlim)
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
    ylab <- getBrowserTrackSetting(track, "Track_Options", "Y_Axis_Label", default)
    if(ylab == "auto") default else if(ylab == "none") "" else ylab
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
# other track plotting functions
#----------------------------------------------------------------------
trackLegend.browserTrack <- function(track, coord, ylim, bty = "n", ...){
    par(xpd = TRUE)
    legend(coord$end + coord$width * 0.01, ylim[2], bty = bty, ...)
    par(xpd = FALSE)
}

#----------------------------------------------------------------------
# trackNav builder support functions
#----------------------------------------------------------------------
trackNavObservers__ <- list()
initTrackNav.browserTrack <- function(track, session, inputName, actionFn = NULL) { # set actionFn for input, but not a table
    req(getBrowserTrackSetting(track, "Track_Options", "Show_Navigation", "hide") != "hide")
    navName <- paste(track$type, track$id, inputName, sep = "_")
    if(!is.null(trackNavObservers__[[navName]])) trackNavObservers__[[navName]]$destroy()
    if(!is.null(actionFn)) trackNavObservers__[[navName]] <<- observeEvent(session$input[[navName]], { 
        actionFn(session$input[[navName]]) 
    })
    navName
}
trackNavInput.browserTrack <- function(track, session, navName, shinyInputFn, 
                          value = NULL, selected = NULL, ...){
    default <- if(is.null(value)) selected else value
    x <- isolate( session$input[[navName]] )
    if(is.null(x)) x <- default
    tags$div(
        class = "trackBrowserInput",
        if(is.null(value)) shinyInputFn(
            session$ns(navName), 
            selected = x,
            ...
        ) else shinyInputFn(
            session$ns(navName), 
            value = x,
            ...
        )
    )
}
trackNavTable.browserTrack <- function(track, session, browserId, navName, 
                                       actionFn, options = list(), ...){
    # x <- isolate( session$input[[navName]] ) # get rows already selected?
    if(is.null(options$lengthMenu)) options$lengthMenu <- c(5,10,15,20,50,100)
    if(is.null(options$pageLength)) options$pageLength <- options$lengthMenu[1]
    bufferedTableServer(
        navName,
        browserId,
        session$input,
        selectionFn = actionFn,
        filterable = TRUE,
        settings = track$settings,
        options = options,
        ...
    )
    bufferedTableUI(
        session$ns(navName),
        title = getTrackDisplayName(track), 
        downloadable = TRUE,
        settings = TRUE,
        width = 12,
        selection = "single",
        style = "display: block;",
        collapsible = TRUE
    )      
}
trackNavCanNavigate.browserTrack <- function(track){
    getBrowserTrackSetting(track, "Track_Options", "Show_Navigation", "hide") %in% c("navigate", "navigate_and_expand")
}
trackNavCanExpand.browserTrack <- function(track){
    getBrowserTrackSetting(track, "Track_Options", "Show_Navigation", "hide") %in% c("expand", "navigate_and_expand")
}
handleTrackNavTableClick <- function(track, chrom, start, end, expandFn = NULL){
    navigate <- trackNavCanNavigate(track)
    expand   <- trackNavCanExpand(track)
    if(navigate && expand){
        app$browser$jumpToCoordinates(chrom, start, end, then = expandFn)
    } else if(navigate){
        app$browser$jumpToCoordinates(chrom, start, end)
    } else if(expand){
        expandFn()
    }
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

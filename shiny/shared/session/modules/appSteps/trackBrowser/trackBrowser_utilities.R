#----------------------------------------------------------------------
# generic support functions for drawing browserTracks
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# MDI (as opposed to UCSC) track image assembly
#----------------------------------------------------------------------
mdiTrackImage <- function(layout, height, plotFn, message = NULL, ...){ # for track that generate plots
    # NB: do not use magick::image_graph; it is buggy, especially on Linux
    if(!is.null(message)) startSpinner(session, message = message)
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
    image <- tryCatch({
        plotFn(...) 
        dev.off()
        image <- magick::image_read(pngFile)
        unlink(pngFile)
        image
    }, error = function(e){
        message(paste("mdiTrackImage: plotFn error:"))
        print(e)
        dev.off()
        req(FALSE)
    })
    # stopSpinner(session)
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
    trackName <- track$settings$get("Track", "Track_Name")
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
                allowedFields <- itemFields %in% c(
                    "name","Name","
                    Sample_ID","sample","Sample",
                    "track","Track",
                    "library","Library"
                ) # could add more recognized name defaults here
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
    family <- track$settings[[optionFamily]]
    if(is.null(family)) return(default)
    value <- family()[[option]]$value
    if(is.null(value)) return(default)
    if(typeof(value) == "character"){
        value <- trimws(value)
        if(value == "" || value == "auto") return(default)
    }
    value
}
getTrackSetting <- getBrowserTrackSetting

# Track options family
padding.browserTrack <- function(track, layout){
    top    <- getInches(track$settings$get("Track", "Top_Padding"),    layout$lengthUnit)
    bottom <- getInches(track$settings$get("Track", "Bottom_Padding"), layout$lengthUnit)
    list(
        top = top,
        bottom = bottom,
        total = top + bottom
    )
}
height.browserTrack <- function(track, default){
    getBrowserTrackSetting(track, "Track", "Height", default)
}
ylab.browserTrack <- function(track, default = ""){
    ylab <- getBrowserTrackSetting(track, "Track", "Y_Axis_Label", default)
    if(ylab == "auto") default else if(ylab == "none") "" else ylab
}
ylim.browserTrack <- function(track, y, family = "Track", setting = "Y_Limit"){
    user <- trimws(getBrowserTrackSetting(track, family, setting, ""))
    if(user == "") return(paddedRange(y))
    user <- gsub('\\s', '', user)
    user <- gsub('to', '-', user)
    user <- gsub('::', ':', user)
    tryCatch({
        user <- as.numeric(strsplit(user, '[-:,]')[[1]])
        if(is.na(user[2])) user[2] <- -user[1]
        sort(user)          
    }, error = function(e){
        paddedRange(y)
    })
}
bty.browserTrack <- function(track, default){
    user <- getBrowserTrackSetting(track, "Track", "Bounding_Box")
    if(is.null(user)) return(default)
    if(user) "o" else "n"
}

# Plot_Options options family
scaleUnit.browserTrack <- function(track, default = "auto", family = "Plot_Options"){
    getBrowserTrackSetting(track, family, "Scale_Unit", default)
}
typ.browserTrack <- function(track, default, family = "Plot_Options"){
    user <- getBrowserTrackSetting(track, family, "Plot_Type")
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
pch.browserTrack <- function(track, default = 19, family = "Plot_Options"){
    user <- getBrowserTrackSetting(track, family, "Point_Symbol")
    if(is.null(user)) return(default)
    switch(
        user,
        "open_circles" = 1,
        "filled_circles" = 19,
        "open_squares" = 0,
        "filled_squares" = 15
    )
}
lwd.browserTrack <- function(track, default = 1.5, family = "Plot_Options"){
    user <- getBrowserTrackSetting(track, family, "Line_Width")
    if(is.null(user)) return(default)
    user
}
lty.browserTrack <- function(track, default = 1, family = "Plot_Options"){
    user <- getBrowserTrackSetting(track, family, "Line_Type")
    if(is.null(user)) return(default)
    user
}
cex.browserTrack <- function(track, default = 1, family = "Plot_Options"){
    user <- getBrowserTrackSetting(track, family, "Point_Size")
    if(is.null(user)) return(default)
    user
}
col.browserTrack <- function(track, default = CONSTANTS$plotlyColors$blue, family = "Plot_Options"){
    user <- getBrowserTrackSetting(track, family, "Color")
    if(is.null(user)) return(default)
    CONSTANTS$plotlyColor[[user]]
}

#----------------------------------------------------------------------
# plot frame functions
#----------------------------------------------------------------------
zeroLine.browserTrack <- function(track, color = CONSTANTS$plotlyColors$black, lwd = 1){
    abline(h = 0, col = color, lwd = lwd)
}
hLines.browserTrack <- function(track, ylim, color = CONSTANTS$plotlyColors$grey, lwd = 0.25){
    doLines <- getTrackSetting(track, "Track", "Horizontal_Lines", TRUE)
    if(!doLines) return()
    unit <- 10 ** floor(log10(max(abs(ylim))))
    getY <- function(){
        y <- seq(-unit * 10, unit * 10, unit)
        y[y >= ylim[1] & y <= ylim[2]]        
    }
    y <- getY()
    if(length(y) <= 3) {
        unit <- unit / 10
        y <- getY()
    }
    if(length(y) > 10) {
        unit <- unit * 2
        y <- getY()
    }
    abline(h = y, col = color, lwd = lwd)
}
trackLegend.browserTrack <- function(track, coord, ylim, bty = "n", ...){
    par(xpd = TRUE)
    legend(coord$end + coord$width * 0.01, ylim[1] + diff(ylim) / 2, bty = bty, yjust = 0.5, ...)
    par(xpd = FALSE)
}
sideLabelPlacement <- function(layout, coord){
    isLeft <- layout$nRegions == 1 || layout$regionI == 1 || layout$arrangement == "stacked"
    list(
        isLeft = isLeft,
        side = if(isLeft) 2 else 4,
        coord = if(isLeft) coord$start else coord$end,
        marginWidth = if(isLeft) layout$mai$left else layout$mai$right
    )
}

#----------------------------------------------------------------------
# data retrieval and plotting functions
#----------------------------------------------------------------------
getBgzFilesFromItems <- function(track, reference, type){ # return a set of typically large bgz files found on a server but not packaged
    items <- track$settings$items()
    selectedSources <- getSourcesFromTrackSamples(items)
    c(sapply(names(selectedSources), function(sourceId){
        genome <- getSourcePackageOption(sourceId, "genome", "genome") 
        if(genome != reference$genome$genome) return(character())  
        outputDir <- getSourcePackageOption(sourceId, "output", "output-dir")
        dataName <- getSourcePackageOption(sourceId, "output", "data-name")
        selectedSamples <- selectedSources[[sourceId]]$Sample_ID
        bgzFiles <- paste(selectedSamples, genome, type, "bgz", sep = ".")  
        bgzFiles <- c(
            file.path(outputDir, selectedSamples, bgzFiles), # hits when the data package corresponds directly to the pileup file location
            file.path(outputDir, dataName, selectedSamples, bgzFiles) # hits when the data package aggregated a set of sample found in sub-folders
        )
        bgzFiles[file.exists(bgzFiles)]
    }))
}
getItemsData.browserTrack <- function(track, reference, coord, dataFn, parseXY = TRUE){
    items <- track$settings$items()
    req(items, length(items) > 0)
    itemKeys <- names(items) # derived from keyColumn
    ymin <- NA
    ymax <- NA
    d <- lapply(itemKeys, function(itemKey){
        dd <- dataFn(track, reference, coord, itemKey, items[[itemKey]]) # typically, expect dataFn to filter data to the browser window
        if(parseXY){
            y <- if("y1" %in% names(dd)) dd[, c(y1, y2)] else dd[, y]
            ymin <<- min(ymin, y, na.rm = TRUE)
            ymax <<- max(ymax, y, na.rm = TRUE)
        }
        dd
    })
    names(d) <- itemKeys
    list(
        d = d, # at this stage, different items are still separate in a named list
        ymin = ymin,
        ymax = ymax,
        hasData = any(sapply(itemKeys, function(itemKey) {
            nrow(d[[itemKey]] > 0) && 
            (!parseXY || any(!is.na(d[[itemKey]]$x)))
        }))
    )
}
plotXY.browserTrack <- function(track, d, color = NULL, family = "Plot_Options", ...){
    if(is.null(color)) color <- col(track, family = family) 
    if(nrow(d) == 0) return()
    switch(
        getTrackSetting(track, family, "Plot_As", "lines"),
        points = points(
            d$x, 
            d$y, 
            pch = pch(track, family = family), 
            cex = cex(track, family = family),
            col = color, 
            ...
        ),
        both =  points(
            d$x, 
            d$y, 
            typ = "b", 
            pch = pch(track, family = family), 
            cex = cex(track, family = family),
            lwd = lwd(track, family = family), 
            lty = lty(track, family = family),
            col = color, 
            ...
        ),
        histogram = points(
            d$x, 
            d$y, 
            typ = "h", 
            cex = cex(track, family = family),
            col = color, 
            ...
        ),
        area = polygon(
            c(d$x[1], d$x, d$x[length(d$x)]), 
            c(0, d$y, 0), 
            col = color, 
            border = "grey20",
            ...
        ),
        lines(
            d$x, 
            d$y, 
            lwd = lwd(track, family = family), 
            lty = lty(track, family = family),
            col = color, 
            ...
        )
    )
}
plotSpans.browserTrack <- function(track, d, color = NULL, family = "Plot_Options", ...){
    if(nrow(d) == 0) return()
    if(is.null(color)) color <- col(track, family = family) 
    Span_Line_Width <- getTrackSetting(track, family, "Span_Line_Width", 2)
    for(i in 1:nrow(d)){
        span <- d[i]
        lines(c(span$x1 - 0.5, span$x2 + 0.5), rep(span$y, 2), col = color, lwd = Span_Line_Width)
    }
}
plotHeatMap.browserTrack <- function(track, d, color = NULL, exponent = 1, family = "Plot_Options", ...){
    if(is.null(color)) color <- col(track, family = family) 
    if(nrow(d) == 0) return()
    rect(
        d$x1 - 0.5, # expects integer base values in x1 and x2
        d$y  - 0.5, # expects integer line numbers in y
        d$x2 + 0.5, 
        d$y  + 0.5, 
        col = addAlphaToColor(color, 1 - (1 - d$alpha)**exponent), # expects alpha as a fraction of 1
        border = NA
    )
} 

#----------------------------------------------------------------------
# trackNav builder support functions
#----------------------------------------------------------------------
trackNavObservers__ <- list()
initTrackNav.browserTrack <- function(track, session, inputName, actionFn = NULL) { # set actionFn for input, but not a table
    req(getBrowserTrackSetting(track, "Track", "Show_Navigation", "hide") != "hide")
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
    getBrowserTrackSetting(track, "Track", "Show_Navigation", "hide") %in% c("navigate", "navigate_and_expand")
}
trackNavCanExpand.browserTrack <- function(track){
    getBrowserTrackSetting(track, "Track", "Show_Navigation", "hide") %in% c("expand", "navigate_and_expand")
}
handleTrackNavTableClick <- function(regionI, track, chrom, start, end, expandFn = NULL){
    navigate <- trackNavCanNavigate(track)
    expand   <- trackNavCanExpand(track)
    if(navigate && expand){
        app$browser$jumpToCoordinates(regionI, chrom, start, end, then = expandFn)
    } else if(navigate){
        app$browser$jumpToCoordinates(regionI, chrom, start, end)
    } else if(expand){
        expandFn()
    }
}

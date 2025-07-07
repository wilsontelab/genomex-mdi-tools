# trackBrowser server module for assembling the composite browser output image
# there may be one or multiple region images rendered by a single browser instance
trackBrowserImageServer <- function(id, browser, regionI) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------
default <- list(
    arrangement = "stacked",
    labelWidth = 0.75,
    legendWidth = 1.25,
    browserWidth = 7.5
)
observers <- list()

#----------------------------------------------------------------------
# access and parse external data
#----------------------------------------------------------------------
# determine our plot dimensions
browserOptions <- browser$settings$Browser_Options
arrangement <- reactive({
    arrangement <- browserOptions()$Region_Arrangement$value
    if(isTruthy(arrangement)) arrangement else default$arrangement
})
observers$arrangement <- observeEvent(arrangement(), {
    arrangement <- arrangement()
    display <- if(arrangement == "stacked") "block" else "inline-block"
    runjs(paste0('$(".browserTrackImageWrapper").css("display", "', display, '")'))
})
labelWidth <- reactive({
    x <- browserOptions()$Label_Width$value
    if(!isTruthy(x)) x <- default$labelWidth
    if(arrangement() == "stacked" || regionI == 1) x else 0 # show label on left-most plot only
})
legendWidth <- reactive({
    x <- browserOptions()$Legend_Width$value
    if(!isTruthy(x)) x <- default$legendWidth
    if(arrangement() == "stacked" || regionI == browser$nRegions()) x else 0 # show legend on right-most plot only
})
browserWidth <- reactive({
    x <- browserOptions()$Browser_Width$value
    if(!isTruthy(x)) x <- default$browserWidth
    if(arrangement() == "stacked") return(x)
    labelWidth <- browserOptions()$Label_Width$value
    if(!isTruthy(x)) x <- default$labelWidth
    legendWidth <- browserOptions()$Legend_Width$value
    if(!isTruthy(x)) x <- default$legendWidth
    regionWidth <- (x - labelWidth - legendWidth) / browser$nRegions() # disperse all but legends equally between regions
    regionWidth + labelWidth() + legendWidth()
})
#----------------------------------------------------------------------
# genome region to plot
coordinates <- browser$coordinates[[regionI]]
reference <- reactive({ # parse the target genome region
    reference <- browser$reference$reference()
    req(reference$genome)
    reference
})
coord <- reactive({ # parse the plot window
    coord <- coordinates(coordinates$input)
    req(coord)
    req(coord$width > 0)
    coord 
})

#----------------------------------------------------------------------
# plot methods shared between browser and expansion tracks
#----------------------------------------------------------------------
trackCache <- list( # cache for track images to prevent slow replotting of unchanged tracks
    browser   = list(), # names = trackIds, values = list(trackHash, contents)
    expansion = list()
)
forceTrackRefresh <- reactiveValues() # force a track to refresh
buildAllTracks <- function(trackIds, fnName, type, layout, externalCoord = NULL){
    reference <- reference()
    coord <- if(is.null(externalCoord)) coord() else externalCoord  
    sharedHash <- digest(list(reference, coord, layout))
    nImages <- 0 # may be less than or equal to nTracks, depending on build successes
    tracks <- browser$tracks$tracks()
    builds <- lapply(trackIds, function(trackId) {
        tryCatch({
            track <- tracks[[trackId]]$track
            trackHash <- digest(list(
                track$settings$all_(), 
                if(is.null(track$settings$items)) NA else track$settings$items(),
                if(type == "expansion") expandingTrack() else NA,
                sharedHash,
                forceTrackRefresh[[trackId]]
            ))
            if(is.null(trackCache[[type]][[trackId]]) || 
               trackCache[[type]][[trackId]]$trackHash != trackHash) {
                 trackCache[[type]][[trackId]] <<- list(
                    trackHash = trackHash,
                    contents = track[[fnName]](reference, coord, layout)
                 )
            }
            nImages <<- nImages + 1
            trackCache[[type]][[trackId]]$contents
        }, error = function(e) {
            if(serverEnv$IS_DEVELOPER){
                message(paste("buildAllTracks: track error:", tracks[[trackId]]$track$type, trackId))
                print(e)                
            }
            NULL
        })
    })
    if(nImages == 0) {
        stopSpinner(session)
        return(NULL)
    }
    builds
}
assembleCompositeImage <- function(builds, type, pngFile){
    fileName <- paste0(app$NAME, "-", id, "-", type, "-image.png")       
    if(is.null(pngFile)) pngFile <- file.path(sessionDirectory, fileName)
    composite <- magick::image_append(
        magick::image_join(lapply(builds[!sapply(builds, is.null)], function(build) build$image)), 
        stack = TRUE
    )
    magick::image_write(composite, path = pngFile, format = "png")
    stopSpinner(session)
    pngFile
}
adjustLayoutForIP <- function(layout, builds){ # adjust layout for mdiInteractivePlotServer
    layout$heights <- sapply(builds, function(build) if(is.null(build$image)) 0 else magick::image_info(build$image)$height)
    layout$height <- sum(layout$heights)
    layout$browserHeight <- layout$height / layout$dpi
    layout$ylim <- lapply(builds, function(build) build$ylim)
    layout$mai  <- lapply(builds, function(build) build$mai)  
    layout    
}

#----------------------------------------------------------------------
# render the composite browser plot image using base graphics
#----------------------------------------------------------------------
browserLayout <- list()
createBrowserPlot <- function(pngFile = NULL, externalCoord = NULL){ # called to generate plot for both screen and image file download
    isPrint <- !is.null(pngFile) && is.null(externalCoord)

    # collect the track list
    trackIds <- browser$tracks$orderedTrackIds()
    nTracks <- length(trackIds)
    req(nTracks > 0)
    tracks <- browser$tracks$tracks()

    # inherit plot targets from the browser inputs
    reference <- reference()
    coord <- if(is.null(externalCoord)) coord() else externalCoord

    # parse the plot layout based on plots alone
    browserOptions <- browserOptions()
    lengthUnit <- browserOptions$Length_Unit$value
    dpi <- if(isPrint) browser$printDpi else browser$screenDpi
    linesPerInch <- if(isPrint) browser$linesPerInchPrint() else browser$linesPerInchScreen()   
    browserWidth <- getInches(browserWidth(), lengthUnit)
    labelWidth   <- getInches(labelWidth(),   lengthUnit)
    legendWidth  <- getInches(legendWidth(),  lengthUnit)
    plotWidth <- browserWidth - labelWidth - legendWidth

    req(plotWidth > 0)
    startSpinner(session)    
    layout <- list(
        isPrint = isPrint,
        printMultiplier = if(isPrint) browser$printDpi / browser$screenDpi else 1,
        dpi = dpi,
        linesPerInch = linesPerInch,
        browserWidth = browserWidth,
        width = browserWidth * dpi,
        plotWidth = plotWidth,
        mai = list(
            left  = labelWidth, 
            right = legendWidth,
            lengthUnit = lengthUnit 
        ),
        pointsize = as.integer(browserOptions$Font_Size$value),
        arrangement = arrangement(),
        nRegions = browser$nRegions(),
        regionI = regionI
    )

    # adjust the label width for UCSC track behavior
    tracks <- browser$tracks$tracks()
    trackTypes <- sapply(tracks, function(x) x$type)
    if("ucsc_tracks" %in% trackTypes) layout <- adjustLayoutForUcsc(layout) 

    # override browser coordinates and width if exactly 1 track has adjustsWidth = TRUE
    adjustingTrackIds <- unlist(sapply(trackIds, function(trackId) {
        x <- tracks[[trackId]]$track$adjustsWidth
        if(is.null(x) || !x) character() else trackId
    }))
    req(length(adjustingTrackIds <= 1))
    if(length(adjustingTrackIds) == 1){
        tryCatch({
            x <- tracks[[adjustingTrackIds]]$track$adjustWidth(reference, coord, layout)
            # coord <- x$coord
            layout <- x$layout
        }, error = function(e) {
            print(e)
            stopSpinner(session)
            req(FALSE)
        })
    }

    # build all tracks using reactives
    builds <- buildAllTracks(trackIds, "buildTrack", "browser", layout, externalCoord)
    if(is.null(builds)) return(NULL)

    # assemble and return the composite browser image
    list(
        pngFile = assembleCompositeImage(builds, "browser", pngFile),
        preBuildLayout = layout,  
        layout = adjustLayoutForIP(layout, builds),      
        parseLayout = yPixelToTrack
    ) 
}
browserPlot <- mdiInteractivePlotServer(
    "image", 
    contents = reactive({
        req(!browser$isInitializing())
        contents <- createBrowserPlot()
        browserLayout <<- contents[c("preBuildLayout", "layout")]
        isolate({ browserIsDone( browserIsDone() + 1 ) })
        contents
    })
)
browserIsDone <- reactiveVal(0)

# ----------------------------------------------------------------------
# render the expansion plot image using base graphics
# ----------------------------------------------------------------------
expandingTrack <- reactiveVal(NULL)
createExpansionPlot <- function(trackId, pngFile = NULL){ # called to generate plot for both screen and image file download

    # build all expansion track images using reactives
    builds <- buildAllTracks(trackId, "buildExpansion", "expansion", browserLayout$preBuildLayout)
    if(is.null(builds)) return(NULL)

    # assemble and return the composite expansion image
    list(
        pngFile = assembleCompositeImage(builds, "expansion", pngFile),       
        layout = adjustLayoutForIP(browserLayout$layout, builds),  
        parseLayout = yPixelToTrack
    )
}
expansionPlot <- mdiInteractivePlotServer(
    "expansionImage", 
    contents = reactive({
        expandingTrack <- expandingTrack()
        if(is.list(expandingTrack)) {
            shinyjs::show(selector = ".expansionImageWrapper")
            createExpansionPlot(expandingTrack$trackId)
        } else {
            shinyjs::hide(selector = ".expansionImageWrapper")
            NULL 
        }
    })
)

#----------------------------------------------------------------------
# handle user interactions with the browser plot
#----------------------------------------------------------------------
# convert a Y pixel to the plot it landed in for use by mdiInteractivePlotServer
interactingTrack <- NULL
yPixelToTrack <- function(x, y){
    req(browserLayout)
    req(y)
    layout <- browserLayout$layout
    cumHeights <- cumsum(layout$heights)
    i <- which(cumHeights > y)[1]
    y <- if(i == 1) y else y - cumHeights[i - 1] # distance from top of track
    interactingTrack <<- browser$tracks$tracks()[[ browser$tracks$trackOrder()[i, trackId] ]]$track
    coord <- coord()
    list(
        x = x, # no transformation required, tracks are never side by side
        y = y,
        layout = list(
            width  = layout$width,
            height = layout$heights[i],
            dpi    = layout$dpi,
            xlim   = as.integer64(c(coord$start, coord$end)),
            ylim   = layout$ylim[[i]],            
            mai    = layout$mai[[i]]
        )        
    )
}
#----------------------------------------------------------------------
# transmit the click and hover actions to the track
doTrackClick <- function(click){
    req(interactingTrack$click)
    click(interactingTrack, click, regionI)
}
observers$click <- observeEvent(browserPlot$click(), {
    click <- browserPlot$click()
    if(click$keys$alt){ # on all tracks, Alt-click opens the track settings (leaves no-key, Ctrl and Shift for track use)
        interactingTrack$settings$open()
    } else {
        doTrackClick(click)
    }
}, ignoreInit = FALSE)
observers$hover <- observeEvent(browserPlot$hover(), {
    req(interactingTrack$hover)
    hover(interactingTrack, browserPlot$hover(), regionI)
}, ignoreInit = FALSE)
#----------------------------------------------------------------------
# transmit brush action to in-window zoom by default
observers$brush <- observeEvent(browserPlot$brush(), {
    brush <- browserPlot$brush() 
    if(isTruthy(interactingTrack$forceBrush)) {
        brush(interactingTrack, brush, regionI)
    } else if(all(!as.logical(brush$keys))){ # on all tracks, a no-key brush creates a window coordinate zoom
        x <- sort(c(brush$coord$x1, brush$coord$x2))
        if(x[1] < x[2]){ # make sure this isn't an accidental brush that was supposed to be a click
            coordinates$jumpToCoordinates(NA, x[1], x[2])
        } else { # handle zero-width brushes as clicks
            brush$coord$x <- mean(x)
            brush$coord$y <- mean(brush$coord$y1, brush$coord$y2)
            doTrackClick(brush)
        }  
    } else if(isTruthy(interactingTrack$brush)) { # tracks optionally handles key+brush events
        brush(interactingTrack, brush, regionI)
    }
}, ignoreInit = FALSE)

#----------------------------------------------------------------------
# initialization
#----------------------------------------------------------------------
initialize <- function(jobId, loadData, loadSequence = NULL){
    # do nothing at present, here for future considerations
    if(!is.null(loadSequence)) doNextLoadSequenceItem(loadData, loadSequence)
}

# module return value
list(
    expandingTrack = expandingTrack,
    browserIsDone = browserIsDone,
    createBrowserPlot = createBrowserPlot,
    initialize = initialize,
    destroy = function(){
        for(observer in observers) observer$destroy()
        observers <<- list()
    },
    forceTrackRefresh = function(trackId){
        if(is.null(forceTrackRefresh[[trackId]])) forceTrackRefresh[[trackId]] <- 1
        else forceTrackRefresh[[trackId]] <- forceTrackRefresh[[trackId]] + 1
    }
)
#----------------------------------------------------------------------
})}

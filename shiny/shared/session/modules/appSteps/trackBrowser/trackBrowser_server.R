#----------------------------------------------------------------------
# server components for the trackBrowser appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
trackBrowserServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# set initialization defaults
#----------------------------------------------------------------------
defaultGenome <- "hg38"
defaultAnnotation <- "wgEncodeGencodeBasicV41"
defaultTrackTypes <- c(
    "plot_title",
    "chromosome",
    "scale_bar",
    "coordinate_axis"
)

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
suite  <- 'genomex-mdi-tools'
module <- 'trackBrowser'
class(input) <- unique(append("browserInput", class(input)))
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    download = downloadHandler(
        filename = paste0(app$NAME, "-", id, ".png"), 
        content = function(pngFile) createBrowserPlot(pngFile)
    ),
    settings = id, # for step-level settings
    size = "m",
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)
browserIsInitialized <- reactiveVal(FALSE)
isLoadingDefaultGenome <- FALSE
confirmBrowserInit <- function(...) {
    hide("initMessage")
    browserIsInitialized(TRUE) 
    genome <- genome()
    if(is.null(genome)    || genome == "")      genomeInput(defaultGenome)
    annotation <- annotation()
    if(is.null(annotation) || annotation == "") {
        isLoadingDefaultGenome <<- TRUE
        annotationInput(defaultAnnotation)
    }
    if(length(names(tracks)) == 0) isolate({ # set default header tracks if not loading a bookmark
        for(i in seq_along(c(defaultTrackTypes, options$defaultTrackTypes))){
            trackId <- as.character(i)
            cssId <- paste("track", trackId, sep = "_")
            tracks[[trackId]] <- initTrack(cssId, trackId, defaultTrackTypes[i])
        }
    })
}
setTimeout(confirmBrowserInit, delay = 1000)
addTrackPrompt <- "-- add a new track --"
duplicateTrackPrompt <- "-- duplicate a track --"
duplicateTrackPromptId <- "0"

#----------------------------------------------------------------------
# fill cascading genome > chromosomes
#----------------------------------------------------------------------
genomes <- reactiveVal( listUcscGenomes() )
observeEvent(input$getUcscGenomes, {
    genomes( listUcscGenomes(force = TRUE) )
})
genomeI <- reactiveVal(NULL)
genomesTable <- bufferedTableServer(
    "genomes",
    id,
    input,
    tableData = genomes,
    select = genomeI,
    options = list( searchDelay = 0 )
)
genome <- reactive({
    genomeInput <- tryCatch({ get("genomeInput") }, error = function(e) NULL)
    if(!is.null(genomeInput)) genomeInput() else NULL
})
genomeInput <- popupInputServer(
    "genome", 
    "Set Working Genome", 
    callback = function(...){
        rowI <- genomesTable$selectionObserver()
        if(is.na(rowI)) {
            genomeI(NULL)
            NA
        } else {
            genomeI(rowI)
            genomes()[rowI]$genome # the popup's return value
        }
    },
    tags$p(tags$strong({
        x <- genome()
        if(is.null(x) || is.na(x)) "" else paste("current selection = ", genome())
    })),
    bufferedTableUI(session$ns("genomes")),
    actionLink(session$ns("getUcscGenomes"), "Reload from UCSC")
)
#----------------------------------------------------------------------
annotations <- reactiveVal( NULL )
observeEvent(input$getUcscAnnotations, {
    genome <- genome()
    req(genome)
    genomes( listUcscAnnotations(genome = genome, force = TRUE) )
})
annotationI <- reactiveVal(NULL)
annotationsTable <- bufferedTableServer(
    "annotations",
    id,
    input,
    tableData = annotations,
    select = annotationI,
    options = list( searchDelay = 0 )
)
annotation <- reactive({
    annotationInput <- tryCatch({ get("annotationInput") }, error = function(e) NULL)
    if(!is.null(annotationInput)) annotationInput() else NULL
})
annotationInput <- popupInputServer(
    "annotation", 
    "Set Working Annotation", 
    callback = function(...){
        rowI <- annotationsTable$selectionObserver()
        if(is.na(rowI)) {
            annotationI(NULL)
            NA
        } else {
            annotationI(rowI)
            annotations()[rowI]$track # the popup's return value
        }
    },
    tags$p(tags$strong({
        x <- annotation()
        if(is.null(x) || is.na(x)) "" else paste("current selection = ", annotation())
    })),
    bufferedTableUI(session$ns("annotations")),
    actionLink(session$ns("getUcscAnnotations"), "Reload from UCSC"),
    active = reactive({ 
        genome <- genome()
        !is.null(genome) && !is.na(genome)
    })    
)
#----------------------------------------------------------------------
chromosomes <- reactiveVal(NULL)
observeEvent(genome(), {
    genome <- genome()
    req(genome)
    annotations(listUcscAnnotations(genome))
    if(loadingBookmark){
        loadingBookmark <<- FALSE
    } else {
        if(isLoadingDefaultGenome) isLoadingDefaultGenome <- FALSE
        else annotationInput(NULL)            
        chromosomes(c(listCanonicalChromosomes(genome), "all"))
        freezeReactiveValue(input, "chromosome")
        updateSelectInput(session, "chromosome", choices = chromosomes(), selected = NULL)        
    }
})
chromosomeSize <- reactive({
    genome <- genome()
    req(genome)
    chrom <- input$chromosome
    req(chrom)
    if(chrom == "all") getGenomeSize(genome)
    else listUcscChromosomes(genome)[chromosome == chrom, size]
})

#----------------------------------------------------------------------
# manage browser tracks
#----------------------------------------------------------------------
labelTrackTypes <- c(
    "plot_title",
    "chromosome",
    "scale_bar",
    "coordinate_axis"
)
trackTypes <- list()       # key = trackType, value = settings template file path
tracks <- reactiveValues() # key = trackId,   value = browserTrackServer()
nullTrackOrder <- data.table(trackId = character(), order = integer())
trackOrder <- reactiveVal(nullTrackOrder)

# assemble the track types available to this app
addTrackType <- function(trackType, tracksFolder){
    trackTypes[[trackType]] <<- file.path(tracksFolder, trackType, "settings.yml")
}
addTrackTypes <- function(dir, classPath){
    tracksFolder <- file.path(dir, classPath)
    trackTypes <- list.dirs(tracksFolder, full.names = FALSE, recursive = FALSE)
    sapply(trackTypes, addTrackType, tracksFolder)
}
classPath <- "shiny/shared/global/classes/browserTracks"
genomexDirs <- parseExternalSuiteDirs(suite)
addTrackTypes(genomexDirs$suiteDir, classPath)
addTrackTypes(gitStatusData$suite$dir, classPath)
classPath <- "classes/browserTracks"
addTrackTypes(app$DIRECTORY, classPath)

# initialize available tracks
initTrackTypes <- observe({ 
    names <- names(trackTypes)
    sortedTrackTypes <- c(labelTrackTypes, sort(names[!(names %in% labelTrackTypes)]))
    updateSelectInput(
        session, 
        "addTrack", 
        choices = c(
            addTrackPrompt,
            sortedTrackTypes
        ),
        selected = addTrackPrompt
    )
    initTrackTypes$destroy()
})

# handle track addition from select input or bookmark
initTrack <- function(cssId, trackId, trackType){
    track <- browserTrackServer(
        cssId = cssId,
        trackId = trackId,
        trackType = trackType,
        settingsFile = trackTypes[[trackType]],
        browserInput = input,
        genome = genome,
        annotation = annotation,
        # size = NULL,
        # cacheKey = NULL, # a reactive/reactiveVal that returns an id for the current settings state
        # fade = FALSE,
        title = paste("Track parameters (", trackType, ")"),
        # immediate = FALSE, # if TRUE, setting changes are transmitted in real time
        # resettable = TRUE  # if TRUE, a Reset All Setting link will be provided    
    )

    # append the new track at the bottom of the track list
    insertUI(
        paste0("#", session$ns("trackList"), " .rank-list"),
        where = "beforeEnd",
        multiple = FALSE,
        immediate = TRUE,
        browserTrackUI(session$ns(cssId), track)
    )
    trackOrder <- trackOrder()
    trackOrder <- rbind(
        trackOrder, 
        data.table(trackId = trackId, order = nrow(trackOrder) + 1)
    )
    trackOrder(trackOrder)
    track
}
observeEvent(input$addTrack, {

    # parse the track request
    trackType <- input$addTrack
    req(trackType)
    req(trackType != addTrackPrompt)
    updateSelectInput(session, "addTrack", selected = addTrackPrompt) # reset the prompt

    # create the new track
    trackId <- as.character(max(0, as.integer(names(tracks))) + 1)
    cssId <- paste("track", trackId, sep = "_")
    tracks[[trackId]] <- initTrack(cssId, trackId, trackType)
}, ignoreInit = TRUE)

# handle track addition from duplication of an existing track
getTrackNames <- function(trackIds){
    sapply(trackIds, function(trackId){
        track <- tracks[[trackId]]
        trackType <- track$type
        trackName <- track$track$settings$get("Track_Options", "Track_Name")
        if(!is.null(trackName) && trackName != "auto" && trackName != "") trackName
        else trackType
    })
}
output$duplicateTrack <- renderUI({
    trackIds <- plotTrackIds()
    req(trackIds)
    names(trackIds) <- getTrackNames(trackIds)
    promptId <- duplicateTrackPromptId
    names(promptId) <- duplicateTrackPrompt
    selectInput(session$ns("duplicateTrackSelect"), NULL, choices = c(promptId, trackIds))
})
observeEvent(input$duplicateTrackSelect, {

    # parse the track request
    dupTrackId <- input$duplicateTrackSelect
    req(dupTrackId != duplicateTrackPromptId)
    updateSelectInput(session, "duplicateTrackSelect", selected = duplicateTrackPromptId) # reset the prompt

    # create the new track
    dupTrack <- tracks[[dupTrackId]]
    trackId <- as.character(max(0, as.integer(names(tracks))) + 1)
    cssId <- paste("track", trackId, sep = "_")
    tracks[[trackId]] <- initTrack(cssId, trackId, dupTrack$type)
    tracks[[trackId]]$track$settings$replace(dupTrack$track$settings$all_())
}, ignoreInit = TRUE)

# handle track reordering and deletion
isRankListInit <- FALSE
observeEvent({
    input$trackRankList
    input$deleteRankList
}, {

    # parse the request
    currentTrackIds <- trackOrder()[, trackId]    
    newTrackIds <- sapply(strsplit(input$trackRankList, '\\s+'), function(x) x[length(x)])

    # declare the new track order
    nTracks <- length(newTrackIds)
    if(isRankListInit) {
        trackOrder(if(nTracks > 0){
            data.table(trackId = newTrackIds, order = 1:nTracks)
        } else nullTrackOrder)

        # delete tracks as needed
        for(trackId in currentTrackIds) 
            if(!(trackId %in% newTrackIds)) tracks[[trackId]] <- NULL
        removeUI(".trackDeleteTarget .browserTrack")
    } else isRankListInit <<- TRUE
})

#----------------------------------------------------------------------
# browser-level tracks metadata
#----------------------------------------------------------------------
browserLayout <- list()
plotTrackIds <- reactive({ # the current track ids, in plotting order
    trackOrder <- trackOrder()
    if(nrow(trackOrder) > 0) trackOrder[order(order), trackId] else character()
})
coordinateWidth <- reactive({
    as.integer(input$end) - as.integer(input$start) + 1
})

#----------------------------------------------------------------------
# render the composite plot image using base graphics
#----------------------------------------------------------------------
screenDpi <- 96 # TODO: expose as settings?
printDpi  <- 300
getLinesPerInch <- function(dpi){ # conversion between lines and inches based on font size
    fileName <- paste0(app$NAME, "-", id, "-getLinesPerInch.png")     
    pngFile <- file.path(sessionDirectory, fileName)
    png(
        pngFile,
        width = 5, # dimensions not important to result
        height = 5,
        units = "in",
        pointsize = settings$Browser_Options()$Font_Size$value,
        res = dpi,
        type = "cairo"
    )
    linesPerInch <- par("mar")[1] / par("mai")[1]
    dev.off()
    unlink(pngFile)     
    linesPerInch   
}
linesPerInchScreen <- reactive( getLinesPerInch(screenDpi) )
linesPerInchPrint  <- reactive( getLinesPerInch(printDpi) )

createBrowserPlot <- function(pngFile = NULL){ # called to generate plot for both screen and image file download
    isPrint <- !is.null(pngFile)

    # collect the track list
    trackIds <- plotTrackIds()
    nTracks <- length(trackIds)    
    req(nTracks > 0)

    # parse the target genome region
    reference <- list(
        genome = genome(),
        annotation = annotation()
    )
    req(reference$genome)
    coord <- coordinates(input)
    req(coord)
    req(coord$width > 0)

    # parse the plot layout based on plots alone
    browserSettings <- settings$Browser_Options()
    lengthUnit <- browserSettings$Length_Unit$value
    dpi <- if(isPrint) printDpi else screenDpi
    linesPerInch <- if(isPrint) linesPerInchPrint() else linesPerInchScreen()
    browserWidth <- getInches(browserSettings$Browser_Width$value, lengthUnit)
    labelWidth   <- getInches(browserSettings$Label_Width$value,   lengthUnit)
    legendWidth  <- getInches(browserSettings$Legend_Width$value,  lengthUnit)
    plotWidth <- browserWidth - labelWidth - legendWidth
    req(plotWidth > 0)
    startSpinner(session)    
    layout <- list(
        isPrint = isPrint,
        printMultiplier = if(isPrint) printDpi / screenDpi else 1,
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
        pointsize = as.integer(browserSettings$Font_Size$value)
    )

    # adjust the label width for UCSC track behavior
    # TODO: only do this if there are UCSC tracks?
    layout <- adjustLayoutForUcsc(layout) 

    # override browser coordinates and width if exactly 1 track has adjustsWidth = TRUE
    adjustingTrackIds <- unlist(sapply(trackIds, function(trackId) {
        x <- tracks[[trackId]]$track$adjustsWidth
        if(is.null(x) || !x) character() else trackId
    }))
    req(length(adjustingTrackIds <= 1))
    if(length(adjustingTrackIds) == 1){
        tryCatch({
            x <- tracks[[adjustingTrackIds]]$track$adjustWidth(reference, coord, layout)
            coord <- x$coord
            layout <- x$layout
        }, error = function(e) {
            print(e)
            stopSpinner(session)
            req(FALSE)
        })
    }

    # build all tracks using reactives
    builds <- lapply(trackIds, function(trackId) {
        tryCatch({
            tracks[[trackId]]$track$build(reference, coord, layout)
        }, error = function(e) {
            print(e)
            NULL
        })
    })

    # purge null tracks that depend on information not yet available
    builds <- builds[!sapply(builds, is.null)]
    if(length(builds) == 0) {
        stopSpinner(session)
        return(NULL)
    }

    # adjust layout for mdiInteractivePlotServer
    layout$heights <- sapply(builds, function(build) magick::image_info(build$image)$height)
    layout$height <- sum(layout$heights)
    layout$browserHeight <- layout$height / dpi
    layout$ylim <- lapply(builds, function(build) build$ylim)
    layout$mai  <- lapply(builds, function(build) build$mai)

    # assemble to composite browser image
    fileName <- paste0(app$NAME, "-", id, "-image.png")       
    pngFile <- if(isPrint) pngFile else file.path(sessionDirectory, fileName) 
    composite <- magick::image_append(
        magick::image_join(lapply(builds, function(build) build$image)), 
        stack = TRUE
    )
    magick::image_write(composite, path = pngFile, format = "png")

    # finish up
    stopSpinner(session)
    list(
        pngFile = pngFile,
        layout = layout,
        parseLayout = yPixelToTrack
    ) 
}
browser <- mdiInteractivePlotServer(
    "image", 
    contents = reactive({
        req(browserIsInitialized())
        contents <- createBrowserPlot()
        browserLayout <<- contents$layout
        contents
    })
)

#----------------------------------------------------------------------
# handle user interactions with plot
#----------------------------------------------------------------------

# convert a Y pixel to the plot it landed in for use by mdiInteractivePlotServer
interactingTrack <- NULL
yPixelToTrack <- function(x, y){
    req(browserLayout)
    req(y)
    cumHeights <- cumsum(browserLayout$heights)
    i <- which(cumHeights > y)[1]
    y <- if(i == 1) y else y - cumHeights[i - 1] # distance from top of track
    interactingTrack <<- tracks[[ trackOrder()[i, trackId] ]]$track
    list(
        x = x, # no transformation required, tracks are never side by side
        y = y,
        layout = list(
            width  = browserLayout$width,
            height = browserLayout$heights[i],
            dpi    = browserLayout$dpi,
            xlim   = as.integer(c(input$start, input$end)),
            ylim   = browserLayout$ylim[[i]],            
            mai    = browserLayout$mai[[i]]
        )        
    )
}

# transmit the click and hover actions to the track
observeEvent(browser$click(), {
    req(interactingTrack$click)
    d <- browser$click()
    # TODO: react to ctrl click via d$keys$ctrl, open track settings
    click(interactingTrack, d$coord$x, d$coord$y)
})
observeEvent(browser$hover(), {
    req(interactingTrack$hover)
    d <- browser$hover()
    hover(interactingTrack, d$coord$x, d$coord$y)
})

# transmit brush action to in-window zoom by default
observeEvent(browser$brush(), {
    d <- browser$brush()$coord   
    brush <- interactingTrack$brush
    if(is.null(brush) || !brush){
        jumpToCoordinates(input$chromosome, d$x1, d$x2)
    } else {
        brush(
            interactingTrack, 
            x1 = min(d$x1, d$x2), 
            y1 = min(d$y1, d$y2), 
            x2 = max(d$x1, d$x2), 
            y2 = max(d$y1, d$y2)
        )
    }
})

#----------------------------------------------------------------------
# browser navigation support functions
#----------------------------------------------------------------------
jumpToCoordinates <- function(chromosome, start, end, strict = FALSE){ # arguments are strict coordinates
    start <- as.integer(start)
    end   <- as.integer(end)
    if(start > end){
        tmp <- start
        start <- end
        end <- tmp
    }
    if(!input$strict && !strict){
        padding <- (end - start + 1) * 0.05
        start <- as.integer(start - padding)
        end   <- as.integer(end   + padding)
    }
    chromosomeSize <- chromosomeSize()
    req(chromosomeSize)
    if(start < 1) start <- 1
    if(end > chromosomeSize) end <- chromosomeSize
    updateSelectInput(session, "chromosome", selected = chromosome)
    updateTextInput(session, "start", value = as.character(start))
    updateTextInput(session, "end",   value = as.character(end))
}
doZoom <- function(exp){
    start  <- as.integer(input$start)
    end    <- as.integer(input$end)
    factor <- as.integer(input$zoomFactor)
    width  <- (end - start + 1)
    center <- start + width / 2
    newHalfWidth <- (width * factor ** exp) / 2
    jumpToCoordinates(
        input$chromosome, 
        center - newHalfWidth, 
        center + newHalfWidth, 
        strict = TRUE
    )
}
observeEvent(input$zoomOut, { doZoom( 1) }, ignoreInit = TRUE)
observeEvent(input$zoomIn,  { doZoom(-1) }, ignoreInit = TRUE)
doMove <- function(factor, direction){
    start  <- as.integer(input$start)
    end    <- as.integer(input$end)
    width  <- (end - start + 1)
    increment <- width * factor * direction
    jumpToCoordinates(
        input$chromosome, 
        start + increment, 
        end + increment, 
        strict = TRUE
    )
}
observeEvent(input$moveLeft,   { doMove(1,    -1) }, ignoreInit = TRUE)
observeEvent(input$nudgeLeft,  { doMove(0.05, -1) }, ignoreInit = TRUE)
observeEvent(input$nudgeRight, { doMove(0.05,  1) }, ignoreInit = TRUE)
observeEvent(input$moveRight,  { doMove(1,     1) }, ignoreInit = TRUE)
observeEvent(input$all,        { 
    jumpToCoordinates(input$chromosome, 1, 1e9, strict = TRUE) }, 
    ignoreInit = TRUE
)
center <- function(x){
    coord <- coordinates(input)
    halfWidth <- coord$width / 2
    jumpToCoordinates(
        coord$chromosome, 
        x - halfWidth, 
        x + halfWidth, 
        strict = TRUE
    )
}

#----------------------------------------------------------------------
# additional within-track navigation actions, e.g., scrolling through a stack
#----------------------------------------------------------------------
output$trackNavs <- renderUI({

    # process track list
    trackIds <- plotTrackIds()
    nTracks <- length(trackIds)    
    req(nTracks > 0)
    trackNames <- getTrackNames(trackIds)
    names(trackNames) <- trackIds

    # extract any required navigation input rows
    nNavs <- 0
    navs <- lapply(trackIds, function(trackId) {
        track <- tracks[[trackId]]$track
        hasNav <- !is.null(track$navigation) && track$navigation
        if(!hasNav) return(NULL)
        nNavs <<- nNavs + 1
        list(
            ui = navigation(track, session),
            name = trackNames[[trackId]]
        )
    })

    # if needed, populate the trackNav rows
    if(nNavs == 0) NULL else lapply(navs, function(nav){
        if(is.null(nav)) return(NULL)
        tags$div(
            style = "margin-bottom: 8px;",
            tags$p(
                style = "display: inline-block; margin: 30px 10px 0px 10px",
                tags$strong(nav$name)
            ),
            nav$ui
        )
    })
})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
loadingBookmark <- FALSE
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    genomeInput(bm$outcomes$genome)
    annotationInput(bm$outcomes$annotation)
    chromosomes(bm$outcomes$chromosomes)
    updateSelectInput(session,   'chromosome',  choices = bm$outcomes$chromosomes, selected = bm$input$chromosome)
    updateTextInput(session,     'start',       value    = bm$input$start)
    updateTextInput(session,     'end',         value    = bm$input$end)
    updateTextInput(session,     'zoomFactor',  value    = bm$input$zoomFactor)
    updateCheckboxInput(session, 'strict',      value    = bm$input$strict)
    trackIds <- bm$outcomes$trackOrder[order(order), trackId]
    isolate({
        lapply(trackIds, function(trackId){
            track <- bm$outcomes$tracks[[trackId]]        
            tracks[[trackId]] <- initTrack(track$cssId, trackId, track$type)
            tracks[[trackId]]$track$settings$replace(track$settings)
            if(!is.null(track$items)) tracks[[trackId]]$track$settings$items(track$items)
        })
    })
    loadingBookmark <<- TRUE
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = list(
        genome = genomeInput,
        annotation = annotationInput,        
        chromosomes = chromosomes,
        trackOrder = trackOrder,
        tracks = reactive({
            trackIds <- trackOrder()[, trackId]
            x <- lapply(trackIds, function(trackId){
                track <- tracks[[trackId]]
                list(
                    cssId = track$cssId,  
                    type = track$type,
                    settings = track$track$settings$all_(),
                    items = if(is.null(track$track$settings$items)) NULL 
                            else track$track$settings$items()
                )
            })
            names(x) <- trackIds
            x
        })
    ),
    jumpToCoordinates = jumpToCoordinates,
    center = center,
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

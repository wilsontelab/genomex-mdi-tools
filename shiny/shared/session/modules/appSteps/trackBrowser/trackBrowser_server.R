#----------------------------------------------------------------------
# server components for the trackBrowser appStep module
#----------------------------------------------------------------------
library(bit64)

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
trackBrowserServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module options and settings
#----------------------------------------------------------------------
suite  <- 'genomex-mdi-tools'
module <- 'trackBrowser'
class(input) <- unique(append("browserInput", class(input)))
appStepDir <- getAppStepDir(module)
browser <- list( # the single object passed to modules for access to appStep and sub-module-level objects
    isInitializing = reactiveVal(TRUE),
    id = id,
    input = input,
    session = session
)
browser$options <- setDefaultOptions(options, stepModuleInfo[[module]])
browser$settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    download = downloadHandler(
        filename = paste0(app$NAME, "-", id, ".png"), 
        content = constructDownloadComposite
    ),
    settings = id, # for step-level settings
    size = "m",
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# initialize the browser sub-modules
#----------------------------------------------------------------------
browser$reference   <- trackBrowserReferenceServer("reference", browser)
browser$coordinates <- list()
browser$tracks      <- trackBrowserTracksServer("tracks", browser)
browser$images      <- list()

#----------------------------------------------------------------------
# support simultaneous plotting of one or multiple genome regions
# each region has its own set of coordinate inputs and associated output image
#----------------------------------------------------------------------
browser$nRegions <- reactive({ 
    browser$settings$get("Browser_Options", "Number_of_Regions", 1) 
})
addRegionCoordinates <- function(regionI, initialize){
    elementId <- paste0("coordinates", regionI)
    insertUI(".trackBrowserCoordinateInputs", "beforeEnd", immediate = TRUE, ui = trackBrowserCoordinatesUI(session$ns(elementId), regionI))
    browser$coordinates[[regionI]] <<- trackBrowserCoordinatesServer(elementId, browser, regionI)
    initFn <- browser$coordinates[[regionI]]$initialize
    if(initialize) initFn(NA, list()) else initFn
}
removeRegionCoordinates <- function(regionI){
    elementId <- paste0("coordinates", regionI)
    removeUI(paste0("#", session$ns(elementId)), immediate = TRUE)
    browser$coordinates[[regionI]]$destroy()
    browser$coordinates[[regionI]] <<- NULL
}
addRegionImage <- function(regionI, initialize){
    elementId <- paste0("image", regionI)
    insertUI(".trackBrowserImages", "beforeEnd", immediate = TRUE, ui = trackBrowserImageUI(session$ns(elementId), regionI))
    browser$images[[regionI]] <<- trackBrowserImageServer(elementId, browser, regionI)
    initFn <- browser$images[[regionI]]$initialize
    if(initialize) initFn(NA, list()) else initFn
}
removeRegionImage <- function(regionI){
    elementId <- paste0("image", regionI)
    removeUI(paste0("#", session$ns(elementId)), immediate = TRUE)
    browser$images[[regionI]]$destroy()
    browser$images[[regionI]] <<- NULL
}
observeEvent(browser$nRegions(), {
    req(!browser$isInitializing())
    nRegions <- browser$nRegions()
    nCurrentRegions <- length(browser$coordinates)
    req(nRegions, nCurrentRegions > 0)
    if(nRegions < 1) nRegions <- 1
    if(nRegions < nCurrentRegions) {
        for(regionI in (nRegions + 1):nCurrentRegions){
            removeRegionCoordinates(regionI)
            removeRegionImage(regionI)
        }
    } else if(nRegions > nCurrentRegions) {
        for(regionI in (nCurrentRegions + 1):nRegions){
            addRegionCoordinates(regionI, TRUE)
            addRegionImage(regionI, TRUE)
        }
    }
}, ignoreInit = TRUE)

#----------------------------------------------------------------------
# image scaling support
# calculation only needs to be made once for all regions whenever font size changes
#----------------------------------------------------------------------
browser$screenDpi <- 96 # TODO: expose as settings?
browser$printDpi  <- 300
getLinesPerInch <- function(dpi){ # conversion between lines and inches based on font size
    fileName <- paste0(app$NAME, "-", id, "-getLinesPerInch.png") 
    pngFile <- file.path(sessionDirectory, fileName)
    png(
        pngFile,
        width = 5, # dimensions not important to result
        height = 5,
        units = "in",
        pointsize = browser$settings$Browser_Options()$Font_Size$value,
        res = dpi,
        type = "cairo"
    )
    linesPerInch <- par("mar")[1] / par("mai")[1]
    dev.off()
    unlink(pngFile)   
    linesPerInch   
}
browser$linesPerInchScreen <- reactive( getLinesPerInch(browser$screenDpi) )
browser$linesPerInchPrint  <- reactive( getLinesPerInch(browser$printDpi) )

#----------------------------------------------------------------------
# additional track-level navigation actions, e.g., tabulating a feature list
# these start hidden even if offered by a track; user must enable them, typically for just one track
#----------------------------------------------------------------------
output$trackNavs <- renderUI({

    # process track list
    trackIds <- browser$tracks$orderedTrackIds()
    nTracks <- length(trackIds)    
    req(nTracks > 0)
    trackNames <- browser$tracks$getTrackNames(trackIds)
    names(trackNames) <- trackIds

    # extract any required navigation input rows
    nNavs <- 0
    tracks <- browser$tracks$tracks()
    navs <- lapply(trackIds, function(trackId) {
        track <- tracks[[trackId]]$track
        hasNav <- isTruthy(track$navigation)
        if(!hasNav) return(NULL)
        ui <- tryCatch({ navigation(track, session, id, browser) }, error = function(e) NULL)
        if(is.null(ui)) return(NULL)
        nNavs <<- nNavs + 1
        list(
            ui = ui,
            name = trackNames[[trackId]]
        )
    })

    # if needed, populate the trackNav rows
    if(nNavs == 0) NULL else lapply(navs, function(nav){
        if(is.null(nav)) return(NULL)
        class <- nav$ui[[1]]$attribs$class
        if(!is.null(class) && class == "trackBrowserInput") tags$div(
            style = "margin-bottom: 8px;",
            tags$p(
                style = "display: inline-block; margin: 30px 10px 5px 10px",
                tags$strong(nav$name)
            ),
            nav$ui
        ) else  nav$ui
    })
})

# ----------------------------------------------------------------------
# create a single object description table - any click that sets objectTableData() replaces the table contents
# ----------------------------------------------------------------------
expandedObjectData <- reactiveVal(NULL)
objectTableData <- reactiveVal(NULL)
objectTable <- bufferedTableServer(
    "objectTable",
    id,
    input,
    tableData = objectTableData,
    selection = 'none',
    options = list(
        searching = FALSE,
        paging = FALSE,
        info = FALSE
    )
)
observeEvent(objectTableData(), {
    toggle(selector = ".objectTableWrapper", condition = isTruthy(objectTableData()))
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------
# create a single expansion table - any click that sets expansionTableData() replaces the table contents
# ----------------------------------------------------------------------
expansionTableData <- reactiveVal(NULL)
expandingTrack <- reactiveVal(NULL)
expansionTable <- bufferedTableServer(
    "expansionTable",
    id,
    input,
    tableData = expansionTableData,
    selection = 'single',
    selectionFn = function(selectedRow){
        expandingTrack <- expandingTrack()   
        req(selectedRow, expandingTrack)
        track <- browser$tracks$tracks()[[expandingTrack$trackId]]
        req(track, track$track$expand2)
        ui <- expand2(track$track, browser, expansionTableData()[selectedRow])
        req(ui)
        expansionUI(ui)
    },
    options = list()
)
observeEvent(expansionTableData(), {
    toggle(selector = ".expansionTableWrapper", condition = isTruthy(expansionTableData()))
}, ignoreNULL = FALSE)
clearObjectExpansions <- function(){
    hide(selector = ".browserExpansionWrapper")
    expandingTrack(NULL)
    for(regionI in 1:browser$nRegions()) browser$images[[regionI]]$expandingTrack(NULL)
    objectTableData(NULL)
    expansionTableData(NULL)
    expansionUI(NULL)
}

# ----------------------------------------------------------------------
# futher enable tracks to add add arbitrary bottom content in response to expand2 or other actions
# ----------------------------------------------------------------------
expansionUI <- reactiveVal(NULL)
output$expansionUI <- renderUI({
    ui <- expansionUI()
    req(ui)
    ui
})

#----------------------------------------------------------------------
# construct a composite print quality image for download
#----------------------------------------------------------------------
constructDownloadComposite <- function(pngFile) {
    startSpinner(session, message = "rendering download")
    nRegions <- browser$nRegions()
    png <- if(nRegions == 1) browser$images[[1]]$createBrowserPlot(pngFile) else {
        arrangement <- browser$settings$Browser_Options()$Region_Arrangement$value
        blankImage <-  magick::image_blank(30, 30, color = "white")
        getImage <- function(regionI){
            browser$images[[regionI]]$createBrowserPlot(pngFile)
            magick::image_read(pngFile)
        }
        images <- getImage(1)
        for(regionI in 2:nRegions) images <- c(images, blankImage, getImage(regionI))
        composite <- magick::image_append(images, stack = arrangement == "stacked")
        magick::image_write(composite, path = pngFile, format = "png")
    }
    stopSpinner(session)
    png
}

#----------------------------------------------------------------------
# define bookmarking actions and initialize the browser load sequence
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks, fail = FALSE)  
    loadingBookmark <- isTruthy(bm)
    if(loadingBookmark) browser$settings$replace(bm$settings)    
    loadData <- if(loadingBookmark) bm else list()
    nRegions <- loadData$settings$Browser_Options$Number_of_Regions$value
    if(is.null(nRegions)) nRegions <- 1
    loadSequence <- c(
        list(
            browser$reference$initializeGenome,
            browser$reference$initializeAnnotation
        ),
        lapply(1:nRegions, addRegionCoordinates, FALSE),
        list(
            browser$tracks$initialize
        ),
        lapply(1:nRegions, addRegionImage, FALSE),
        function(...) reportProgress("browser is initialized")
    )
    doNextLoadSequenceItem(loadData, loadSequence)
    browser$isInitializing(FALSE)
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    loadSourceFile = function(...) NULL,  # enable the trackBrowser for use as the one and only app step
    session = browser$session,
    input = browser$input,
    settings = browser$settings$all_,
    # TODO: add trackOutcomes functionality to allow tracks to send information to downstream app steps
    outcomes = list(
        analysisSetName = reactive({ "trackBrowser" }),
        genome      = browser$reference$genomeInput,
        annotation  = browser$reference$annotationInput,  
        coordinates = reactive({ lapply(browser$coordinates, function(x) reactiveValuesToList(x$input)) }),
        trackOrder  = browser$tracks$trackOrder,
        tracks      = browser$tracks$bookmarkTracks
    ),
    browserIsDoneReactive = function(regionI) reactive({ browser$images[[regionI]]$browserIsDone() }),
    jumpToCoordinates = function(regionI, ...) browser$coordinates[[regionI]]$jumpToCoordinates(...),
    center = function(regionI, ...) browser$coordinates[[regionI]]$center(...),
    expandingTrack = function(regionI, trackData){ # to describe the source of the expansion image/table as expandingTrack(list(trackId = trackId, object = xxx))
        expandingTrack(trackData)
        browser$images[[regionI]]$expandingTrack(trackData)
    },
    expandedObjectData = expandedObjectData,   # a potentially larger version of the expanded object than objectTableData
    objectTableData = objectTableData,         # to populate the object description table (e.g., a gene)
    expansionTableData = expansionTableData,   # to populate the object expansion table   (e.g., a gene's transripts)
    expansionUI = expansionUI,                 # arbitrary expansion UI content passed to renderUI
    clearObjectExpansions = clearObjectExpansions,
    externalTrackSuites = browser$tracks$externalTrackSuites,
    forceTrackTypeRefresh = function(trackType, regionI = 1){
        for(trackId in browser$tracks$getTrackIdsByType(trackType)){
            browser$images[[regionI]]$forceTrackRefresh(trackId)
        }
    },
    forceTrackRefresh = function(trackId, regionI = 1) browser$images[[regionI]]$forceTrackRefresh(trackId),
    getTrackSettings = function(trackId) browser$tracks$tracks()[[trackId]]$track$settings,
    isReady = reactive({ getStepReadiness(options$source) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

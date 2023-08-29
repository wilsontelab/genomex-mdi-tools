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
browser <- list(
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
        content = function(pngFile) createBrowserPlot(pngFile) # TODO: build composite image
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
# support simultaneous plotting of 1 or more genome regions
#----------------------------------------------------------------------
browser$nRegions <- reactive({ browser$settings$get("Browser_Options", "Number_of_Regions", 1) })
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

# browserIsInitialized <- reactiveVal(FALSE)
# isLoadingDefaultGenome <- FALSE
# confirmBrowserInit <- function(...) {
#     hide("initMessage")
#     browserIsInitialized(TRUE) 
#     genome <- genome()
#     if(!objectHasData(genome)) {
#         genomes <- genomes()
#         req(objectHasData(genomes))
#         genomeInput(genomes[genome == defaultGenome])
#     }
#     annotation <- annotation()
#     if(!objectHasData(annotation)) {
#         isLoadingDefaultGenome <<- TRUE
#         annotations <- annotations()
#         req(objectHasData(annotations))
#         annotationInput(annotations[track == defaultAnnotation])
#     }
#     if(length(names(tracks)) == 0) isolate({ # set default header tracks if not loading a bookmark
#         for(i in seq_along(c(defaultTrackTypes, options$defaultTrackTypes))){
#             trackId <- getTrackId()
#             cssId <- paste("track", trackId, sep = "_")
#             tracks[[trackId]] <<- initTrack(cssId, trackId, defaultTrackTypes[i])
#             createTrackSettingsObserver(trackId)
#         }
#     })
# }
# setTimeout(confirmBrowserInit, delay = 1000)

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

# # ----------------------------------------------------------------------
# # create a single object description table - any click that sets objectTableData() replaces the table contents
# # ----------------------------------------------------------------------
# objectTableData <- reactiveVal(NULL)
# objectTable <- bufferedTableServer(
#     "objectTable",
#     id,
#     input,
#     tableData = objectTableData,
#     selection = 'none',
#     options = list(
#         searching = FALSE,
#         paging = FALSE,
#         info = FALSE
#     )
# )
# observeEvent(objectTableData(), {
#     toggle(selector = ".objectTableWrapper", condition = isTruthy(objectTableData()))
# }, ignoreNULL = FALSE)

# # ----------------------------------------------------------------------
# # create a single expansion table - any click that sets expansionTableData() replaces the table contents
# # ----------------------------------------------------------------------
# expansionTableData <- reactiveVal(NULL)
# expansionTable <- bufferedTableServer(
#     "expansionTable",
#     id,
#     input,
#     tableData = expansionTableData,
#     selection = 'single',
#     selectionFn = function(selectedRow){
#         expandingTrack <- expandingTrack()     
#         req(selectedRow, expandingTrack)
#         track <- tracks[[expandingTrack$trackId]]
#         req(track, track$track$expand2)
#         expand2(track$track, reference(), coord(), expansionTableData()[selectedRow])
#     },
#     options = list()
# )
# observeEvent(expansionTableData(), {
#     toggle(selector = ".expansionTableWrapper", condition = isTruthy(expansionTableData()))
# }, ignoreNULL = FALSE)
# clearObjectExpansions <- function(){
#     hide(selector = ".browserExpansionWrapper")
#     expandingTrack(NULL)
#     objectTableData(NULL)
#     expansionTableData(NULL)
#     expansionUI(NULL)
# }

# # ----------------------------------------------------------------------
# # futher enable tracks to add add arbitrary bottom content in response to expand[2] actions
# # ----------------------------------------------------------------------
# expansionUI <- reactiveVal(NULL)
# output$expansionUI <- renderUI({
#     ui <- expansionUI()
#     req(ui)
#     ui
# })

#----------------------------------------------------------------------
# define bookmarking actions and initialize the browser
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks, fail = FALSE)  
    loadingBookmark <- isTruthy(bm)
    if(loadingBookmark) browser$settings$replace(bm$settings)    
    loadData <- if(loadingBookmark) bm else list()
    nRegions <- if(loadingBookmark && !is.null(bm$settings$Browser_Options$Number_of_Regions)) bm$settings$Browser_Options$Number_of_Regions$value else 1
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
    initializeNextTrackBrowserElement(loadData, loadSequence)
    browser$isInitializing(FALSE)
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    loadSourceFile = function(...) NULL,  # enable the trackBrowser for use as the one and only app step
    input = browser$input,
    settings = browser$settings$all_,
    # TODO: add trackOutcomes functionality to allow tracks to send information to downstream app steps
    outcomes = list(
        analysisSetName = reactive({ "trackBrowser" }),
        genome      = browser$reference$genomeInput,
        annotation  = browser$reference$annotationInput,  
        coordinates = reactive({ lapply(browser$coordinates, function(x) reactiveValuesToList(x$input)) }),
        trackOrder  = browser$tracks$trackOrder,
        tracks      = browser$tracks$bookmarkTracks,
        NULL
    ),
    # jumpToCoordinates = jumpToCoordinates,
    # center = center,
    # expandingTrack = expandingTrack,          # to describe the source of the expansion image/table as expandingTrack(trackId = trackId, object = xxx)
    # objectTableData = objectTableData,        # to populate the object description table (e.g., a gene)
    # expansionTableData = expansionTableData,  # to populate the object expansion table   (e.g., a gene's transripts)
    # expansionUI = expansionUI,                # arbitrary expansion UI content passed to renderUI
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

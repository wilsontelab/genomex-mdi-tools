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
# set initialization defaults
#----------------------------------------------------------------------
# defaultTrackTypes <- c(
#     "plot_title",
#     "chromosome",
#     "scale_bar",
#     "coordinate_axis",
#     "genes"
# )
# addTrackPrompt <- "-- add a new track --"
# duplicateTrackPrompt <- "-- duplicate a track --"
# duplicateTrackPromptId <- "0"

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
        content = function(pngFile) createBrowserPlot(pngFile) # TODO: build composite image
    ),
    settings = id, # for step-level settings
    size = "m",
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# initialize the browser in an ordered sequence of events
# the browser loads empty, then fills either from bookmark or defaults
#----------------------------------------------------------------------
reference   <- trackBrowserReferenceServer("reference",         id, options, input, settings)
coordinates <- list()
# browserIsInitialized <- reactiveVal(FALSE)
# isLoadingDefaultGenome <- FALSE
# getTrackId <- function() gsub("( |:|-)", "_", paste(as.character(Sys.time()), sample.int(1e8, 1)))
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

# currentChromosomeSize <- reactive({ # size of the currently active chromosome (not the one we may be switching to...)
#     reference$getChromosomeSize(input$chromosome)
# })

#----------------------------------------------------------------------
# support simultaneous plotting of 1 or more genome regions
#----------------------------------------------------------------------
nRegions <- reactive({
    settings$get("Browser_Options", "Number_of_Regions", 1)
})
addRegionCoordinates <- function(regionI, initialize){
    coordinatesId <- paste0("coordinates", regionI)
    insertUI(".trackBrowserCoordinateInputs", "beforeEnd", immediate = TRUE, ui = trackBrowserCoordinatesUI(session$ns(coordinatesId), regionI))
    coordinates[[regionI]] <<- trackBrowserCoordinatesServer(coordinatesId, id, options, input, settings, reference, regionI)
    initFn <- coordinates[[regionI]]$initialize
    if(initialize) initFn(NA, list()) else initFn
}
removeRegionCoordinates <- function(regionI){
    coordinatesId <- paste0("coordinates", regionI)
    removeUI(paste0("#", session$ns(coordinatesId)), immediate = TRUE)
    coordinates[[regionI]] <<- NULL
}
observeEvent(nRegions(), {
    nRegions <- nRegions()
    nCurrentRegions <- length(coordinates)
    req(nRegions, nCurrentRegions > 0)
    if(nRegions < 1) nRegions <- 1
    if(nRegions < nCurrentRegions) {
        for(regionI in (nRegions + 1):nCurrentRegions){
            removeRegionCoordinates(regionI)
        }
    } else if(nRegions > nCurrentRegions) {
        for(regionI in (nCurrentRegions + 1):nRegions){
            addRegionCoordinates(regionI, TRUE)
        }
    }
}, ignoreInit = TRUE)

# #----------------------------------------------------------------------
# # manage browser tracks
# #----------------------------------------------------------------------
# trackTypes <- list() # key = trackType, value = settings template file path
# tracks <- list()     # key = trackId,   value = browserTrackServer()
# nullTrackOrder <- data.table(trackId = character(), order = integer())
# trackOrder <- reactiveVal(nullTrackOrder)

# # assemble the track types available to this app
# parentAppTrackTypes <- character()
# addTrackType <- function(trackType, tracksFolder){
#     path <- file.path(tracksFolder, trackType, "settings.yml")
#     if(file.exists(path)) trackTypes[[trackType]] <<- path
# }
# addTrackTypes <- function(dir, classPath, isParentApp = FALSE){
#     tracksFolder <- file.path(dir, classPath)
#     trackTypes <- list.dirs(tracksFolder, full.names = FALSE, recursive = FALSE)
#     if(isParentApp) parentAppTrackTypes <<- sort(trackTypes)
#     sapply(trackTypes, addTrackType, tracksFolder)
# }
# classPath <- "classes/browserTracks"
# globalClassPath <- file.path("shiny/shared/global", classPath) 
# genomexDirs <- parseExternalSuiteDirs(suite)
# addTrackTypes(genomexDirs$suiteDir, globalClassPath) # apps always offer global tracks from genomex-mdi-tools
# addTrackTypes(app$DIRECTORY, classPath, isParentApp = TRUE) # apps always offer tracks they define themselves
# if(!is.null(options$tracks)) for(trackType in options$tracks){ # apps can additionally offer global tracks declared in the app's config.yml that are...
#     if(grepl("//", trackType)){ # ... defined in an external suite, if that suite is actually installed ...
#         x <- strsplit(trackType, "//")[[1]]
#         trackSuiteDirs <- parseExternalSuiteDirs(x[1])
#         if(
#             isTruthy(gitStatusData$dependencies[[x[1]]]$loaded) &&
#             !is.null(trackSuiteDirs)
#         )  addTrackType(x[2], file.path(trackSuiteDirs$suiteDir, globalClassPath)) 
#     } else { # ... or in the parent suite of the app
#         addTrackType(trackType, file.path(gitStatusData$suite$dir, globalClassPath))
#     }
# }

# # initialize available tracks
# initTrackTypes <- observe({ 
#     names <- names(trackTypes)
#     firstTrackTypes <- c(defaultTrackTypes, parentAppTrackTypes)
#     sortedTrackTypes <- c(firstTrackTypes, sort(names[!(names %in% firstTrackTypes)]))
#     updateSelectInput(
#         session, 
#         "addTrack", 
#         choices = c(
#             addTrackPrompt,
#             sortedTrackTypes
#         ),
#         selected = addTrackPrompt
#     )
#     initTrackTypes$destroy()
# })

# # handle track addition from select input or bookmark
# initTrack <- function(cssId, trackId, trackType){
#     track <- browserTrackServer(
#         cssId = cssId,
#         trackId = trackId,
#         trackType = trackType,
#         settingsFile = trackTypes[[trackType]], # includes any presets defined by the trackType
#         presets = options$presets[[trackType]], # add any presets defined by the calling app
#         browserInput = input,
#         genome = genome,
#         annotation = annotation,
#         # size = NULL,
#         # cacheKey = NULL, # a reactive/reactiveVal that returns an id for the current settings state
#         # fade = FALSE,
#         title = paste("Track parameters (", trackType, ")"),
#         # immediate = FALSE, # if TRUE, setting changes are transmitted in real time
#         # resettable = TRUE  # if TRUE, a Reset All Setting link will be provided    
#     )

#     # append the new track at the bottom of the track list
#     insertUI(
#         paste0("#", session$ns("trackList"), " .rank-list"),
#         where = "beforeEnd",
#         multiple = FALSE,
#         immediate = TRUE,
#         browserTrackUI(session$ns(cssId), track)
#     )
#     trackOrder <- trackOrder()
#     trackOrder <- rbind(
#         trackOrder, 
#         data.table(trackId = trackId, order = nrow(trackOrder) + 1)
#     )
#     trackOrder(trackOrder)
#     track
# }
# observeEvent(input$addTrack, {
#     # parse the track request
#     trackType <- input$addTrack
#     req(trackType)
#     req(trackType != addTrackPrompt)
#     updateSelectInput(session, "addTrack", selected = addTrackPrompt) # reset the prompt

#     # create the new track
#     trackId <- getTrackId()
#     cssId <- paste("track", trackId, sep = "_")
#     tracks[[trackId]] <<- initTrack(cssId, trackId, trackType)
#     createTrackSettingsObserver(trackId)
# }, ignoreInit = TRUE)

# # handle track addition from duplication of an existing track
# getTrackNames <- function(trackIds){
#     sapply(trackIds, function(trackId) getTrackDisplayName(tracks[[trackId]]$track))
# }
# output$duplicateTrack <- renderUI({
#     trackIds <- plotTrackIds()
#     req(trackIds)
#     names(trackIds) <- paste0(
#         getTrackNames(trackIds), 
#         " (",
#         sapply(trackIds, function(x) tracks[[x]]$type),
#         ")"
#     )
#     promptId <- duplicateTrackPromptId
#     names(promptId) <- duplicateTrackPrompt
#     selectInput(session$ns("duplicateTrackSelect"), NULL, choices = c(promptId, trackIds))
# })
# observeEvent(input$duplicateTrackSelect, {
#     # parse the track request
#     dupTrackId <- input$duplicateTrackSelect
#     req(dupTrackId != duplicateTrackPromptId)
#     updateSelectInput(session, "duplicateTrackSelect", selected = duplicateTrackPromptId) # reset the prompt

#     # create the new track
#     dupTrack <- tracks[[dupTrackId]]
#     trackId <- getTrackId()
#     cssId <- paste("track", trackId, sep = "_")
#     tracks[[trackId]] <<- initTrack(cssId, trackId, dupTrack$type)
#     tracks[[trackId]]$track$settings$replace(dupTrack$track$settings$all_())    
#     createTrackSettingsObserver(trackId)
# }, ignoreInit = TRUE)

# # handle track reordering and deletion
# isRankListInit <- FALSE
# observeEvent({
#     input$trackRankList
#     input$deleteRankList
# }, {

#     # parse the request
#     currentTrackIds <- trackOrder()[, trackId] 
#     newTrackIds <- if(isTruthy(input$trackRankList)) sapply(strsplit(input$trackRankList, '\\s+'), function(x) x[length(x)])
#                    else currentTrackIds # 082623, unclear why input$trackRankList suddenly starting hitting here with NULL value

#     # declare the new track order
#     nTracks <- length(newTrackIds)
#     if(isRankListInit) {
#         trackOrder(if(nTracks > 0){
#             data.table(trackId = newTrackIds, order = 1:nTracks)
#         } else nullTrackOrder)

#         # delete tracks as needed
#         for(trackId in currentTrackIds) 
#             if(!(trackId %in% newTrackIds)) {
#                 tracks[[trackId]] <<- NULL
#                 trackSettingsObservers[[trackId]] <<- NULL
#                 if(!is.null(trackSettingsUndoId) && trackSettingsUndoId == trackId) trackSettingsUndoId <<- NULL
#             }
#         removeUI(".trackDeleteTarget .browserTrack")
#     } else isRankListInit <<- TRUE
# }, ignoreInit = TRUE)

# # undo the last track settings change, intended for disaster recover, not a complete history tracking
# trackSettingsObservers <- list()
# trackSettingsUndoId <- NULL
# createTrackSettingsObserver <- function(trackId){
#     trackSettingsObservers[[trackId]] <<- observeEvent(tracks[[trackId]]$track$settings$all_(), {
#         trackSettingsUndoId <<- trackId
#         clearObjectExpansions()
#     })
# }
# observeEvent(input$undoTrackSettings, {
#     req(trackSettingsUndoId, tracks[[trackSettingsUndoId]])
#     tracks[[trackSettingsUndoId]]$track$settings$undo()
# })

# #----------------------------------------------------------------------
# # browser-level tracks metadata
# #----------------------------------------------------------------------
# plotTrackIds <- reactive({ # the current track ids, in plotting order
#     trackOrder <- trackOrder()
#     if(nrow(trackOrder) > 0) trackOrder[order(order), trackId] else character()
# })
# coordinateWidth <- reactive({
#     as.integer64(input$end) - as.integer64(input$start) + 1
# })

# #----------------------------------------------------------------------
# # image scaling support
# #----------------------------------------------------------------------
# screenDpi <- 96 # TODO: expose as settings?
# printDpi  <- 300
# getLinesPerInch <- function(dpi){ # conversion between lines and inches based on font size
#     fileName <- paste0(app$NAME, "-", id, "-getLinesPerInch.png")     
#     pngFile <- file.path(sessionDirectory, fileName)
#     png(
#         pngFile,
#         width = 5, # dimensions not important to result
#         height = 5,
#         units = "in",
#         pointsize = settings$Browser_Options()$Font_Size$value,
#         res = dpi,
#         type = "cairo"
#     )
#     linesPerInch <- par("mar")[1] / par("mai")[1]
#     dev.off()
#     unlink(pngFile)     
#     linesPerInch   
# }
# linesPerInchScreen <- reactive( getLinesPerInch(screenDpi) )
# linesPerInchPrint  <- reactive( getLinesPerInch(printDpi) )

# #----------------------------------------------------------------------
# # plot methods shared between browser and expansion tracks
# #----------------------------------------------------------------------
# reference <- reactive({ # parse the target genome region
#     reference <- list(
#         genome = genome(),
#         annotation = annotation()
#     )
#     req(reference$genome)
#     reference
# })
# coord <- reactive({ # parse the plot window
#     coord <- coordinates(input)
#     req(coord)
#     req(coord$width > 0)
#     coord 
# })
# trackCache <- list( # cache for track images to prevent slow replotting of unchanged tracks
#     browser   = list(), # names = trackIds, values = list(trackHash, contents)
#     expansion = list()
# )
# buildAllTracks <- function(trackIds, fnName, type, reference, coord, layout){
#     sharedHash <- digest(list(reference, coord, layout))
#     nImages <- 0 # may be less than or equal to nTracks, depending on build successes

#     builds <- lapply(trackIds, function(trackId) {
#         tryCatch({
#             track <- tracks[[trackId]]$track
#             trackHash <- digest(list(
#                 track$settings$all_(), 
#                 if(is.null(track$settings$items)) NA else track$settings$items(),
#                 if(type == "expansion") expandingTrack() else NA,
#                 sharedHash
#             ))
#             if(is.null(trackCache[[type]][[trackId]]) || 
#                trackCache[[type]][[trackId]]$trackHash != trackHash) {
#                  trackCache[[type]][[trackId]] <<- list(
#                     trackHash = trackHash,
#                     contents = track[[fnName]](reference, coord, layout)
#                  )
#             }
#             nImages <<- nImages + 1
#             trackCache[[type]][[trackId]]$contents
#         }, error = function(e) {
#             message(paste("buildAllTracks: track error:", tracks[[trackId]]$track$type, trackId))
#             print(e)
#             NULL
#         })
#     })
#     if(nImages == 0) {
#         stopSpinner(session)
#         return(NULL)
#     }
#     builds
# }
# assembleCompositeImage <- function(builds, type, pngFile){
#     fileName <- paste0(app$NAME, "-", id, "-", type, "-image.png")       
#     if(is.null(pngFile)) pngFile <- file.path(sessionDirectory, fileName)
#     composite <- magick::image_append(
#         magick::image_join(lapply(builds[!sapply(builds, is.null)], function(build) build$image)), 
#         stack = TRUE
#     )
#     magick::image_write(composite, path = pngFile, format = "png")
#     stopSpinner(session)
#     pngFile
# }
# adjustLayoutForIP <- function(layout, builds){ # adjust layout for mdiInteractivePlotServer
#     layout$heights <- sapply(builds, function(build) if(is.null(build$image)) 0 else magick::image_info(build$image)$height)
#     layout$height <- sum(layout$heights)
#     layout$browserHeight <- layout$height / layout$dpi
#     layout$ylim <- lapply(builds, function(build) build$ylim)
#     layout$mai  <- lapply(builds, function(build) build$mai)  
#     layout    
# }

# #----------------------------------------------------------------------
# # render the composite browser plot image using base graphics
# #----------------------------------------------------------------------
# browserLayout <- list()
# createBrowserPlot <- function(pngFile = NULL){ # called to generate plot for both screen and image file download
#     isPrint <- !is.null(pngFile)

#     # collect the track list
#     trackIds <- plotTrackIds()
#     nTracks <- length(trackIds)    
#     req(nTracks > 0)

#     # inherit plot targets from the browser inputs
#     reference <- reference()
#     coord <- coord()

#     # parse the plot layout based on plots alone
#     browserSettings <- settings$Browser_Options()
#     lengthUnit <- browserSettings$Length_Unit$value
#     dpi <- if(isPrint) printDpi else screenDpi
#     linesPerInch <- if(isPrint) linesPerInchPrint() else linesPerInchScreen()
#     browserWidth <- getInches(browserSettings$Browser_Width$value, lengthUnit)
#     labelWidth   <- getInches(browserSettings$Label_Width$value,   lengthUnit)
#     legendWidth  <- getInches(browserSettings$Legend_Width$value,  lengthUnit)
#     plotWidth <- browserWidth - labelWidth - legendWidth
#     req(plotWidth > 0)
#     startSpinner(session)    
#     layout <- list(
#         isPrint = isPrint,
#         printMultiplier = if(isPrint) printDpi / screenDpi else 1,
#         dpi = dpi,
#         linesPerInch = linesPerInch,
#         browserWidth = browserWidth,
#         width = browserWidth * dpi,
#         plotWidth = plotWidth,
#         mai = list(
#             left  = labelWidth, 
#             right = legendWidth,
#             lengthUnit = lengthUnit 
#         ),
#         pointsize = as.integer(browserSettings$Font_Size$value)
#     )

#     # adjust the label width for UCSC track behavior
#     # TODO: only do this if there are UCSC tracks?
#     layout <- adjustLayoutForUcsc(layout) 

#     # override browser coordinates and width if exactly 1 track has adjustsWidth = TRUE
#     adjustingTrackIds <- unlist(sapply(trackIds, function(trackId) {
#         x <- tracks[[trackId]]$track$adjustsWidth
#         if(is.null(x) || !x) character() else trackId
#     }))
#     req(length(adjustingTrackIds <= 1))
#     if(length(adjustingTrackIds) == 1){
#         tryCatch({
#             x <- tracks[[adjustingTrackIds]]$track$adjustWidth(reference, coord, layout)
#             coord <- x$coord
#             layout <- x$layout
#         }, error = function(e) {
#             print(e)
#             stopSpinner(session)
#             req(FALSE)
#         })
#     }

#     # build all tracks using reactives
#     builds <- buildAllTracks(trackIds, "buildTrack", "browser", reference, coord, layout)
#     if(is.null(builds)) return(NULL)

#     # assemble and return the composite browser image
#     list(
#         pngFile = assembleCompositeImage(builds, "browser", pngFile),
#         preBuildLayout = layout,  
#         layout = adjustLayoutForIP(layout, builds),      
#         parseLayout = yPixelToTrack
#     ) 
# }
# browser <- mdiInteractivePlotServer(
#     "image", 
#     contents = reactive({
#         req(browserIsInitialized())
#         contents <- createBrowserPlot()
#         browserLayout <<- contents[c("preBuildLayout", "layout")]
#         isolate({ browserIsDone( browserIsDone() + 1 ) })
#         contents
#     })
# )
# browserIsDone <- reactiveVal(0)

# # ----------------------------------------------------------------------
# # render the expansion plot image using base graphics
# # ----------------------------------------------------------------------
# expandingTrack <- reactiveVal(NULL)
# createExpansionPlot <- function(trackId, pngFile = NULL){ # called to generate plot for both screen and image file download

#     # build all expansion track images using reactives
#     builds <- buildAllTracks(trackId, "buildExpansion", "expansion", 
#                              reference(), coord(), browserLayout$preBuildLayout)
#     if(is.null(builds)) return(NULL)

#     # assemble and return the composite expansion image
#     list(
#         pngFile = assembleCompositeImage(builds, "expansion", pngFile),       
#         layout = adjustLayoutForIP(browserLayout$layout, builds),  
#         parseLayout = yPixelToTrack
#     )
# }
# expansionImage <- mdiInteractivePlotServer(
#     "expansionImage", 
#     contents = reactive({
#         expandingTrack <- expandingTrack()
#         if(is.list(expandingTrack)) {
#             show(selector = ".expansionImageWrapper")
#             createExpansionPlot(expandingTrack$trackId)            
#         } else {
#             hide(selector = ".expansionImageWrapper")
#             NULL 
#         }
#     })
# )

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

# #----------------------------------------------------------------------
# # handle user interactions with the browser plot
# #----------------------------------------------------------------------

# # convert a Y pixel to the plot it landed in for use by mdiInteractivePlotServer
# interactingTrack <- NULL
# yPixelToTrack <- function(x, y){
#     req(browserLayout)
#     req(y)
#     layout <- browserLayout$layout
#     cumHeights <- cumsum(layout$heights)
#     i <- which(cumHeights > y)[1]
#     y <- if(i == 1) y else y - cumHeights[i - 1] # distance from top of track
#     interactingTrack <<- tracks[[ trackOrder()[i, trackId] ]]$track
#     list(
#         x = x, # no transformation required, tracks are never side by side
#         y = y,
#         layout = list(
#             width  = layout$width,
#             height = layout$heights[i],
#             dpi    = layout$dpi,
#             xlim   = as.integer64(c(input$start, input$end)),
#             ylim   = layout$ylim[[i]],            
#             mai    = layout$mai[[i]]
#         )        
#     )
# }

# # transmit the click and hover actions to the track
# doTrackClick <- function(click){
#     req(interactingTrack$click)     
#     click(interactingTrack, click)    
# }
# observeEvent(browser$click(), {
#     click <- browser$click()
#     if(click$keys$alt){ # on all tracks, Alt-click opens the track settings (leaves no-key, Ctrl and Shift for track use)
#         interactingTrack$settings$open()
#     } else {
#         doTrackClick(click)
#     }
# })
# observeEvent(browser$hover(), {
#     req(interactingTrack$hover)
#     hover(interactingTrack, browser$hover())
# })

# # transmit brush action to in-window zoom by default
# observeEvent(browser$brush(), {
#     brush <- browser$brush() 
#     if(isTruthy(interactingTrack$forceBrush)) {
#         brush(interactingTrack, brush)
#     } else if(all(!as.logical(brush$keys))){ # on all tracks, a no-key brush creates a window coordinate zoom
#         x <- sort(c(brush$coord$x1, brush$coord$x2))
#         if(x[1] < x[2]){ # make sure this isn't an accidental brush that was supposed to be a click
#             jumpToCoordinates(input$chromosome, x[1], x[2])
#         } else { # handle zero-width brushes as clicks
#             brush$coord$x <- mean(x)
#             brush$coord$y <- mean(brush$coord$y1, brush$coord$y2)
#             doTrackClick(brush)
#         }  
#     } else if(isTruthy(interactingTrack$brush)) { # tracks optionally handles key+brush events
#         brush(interactingTrack, brush)
#     }
# })




# #----------------------------------------------------------------------
# # additional within-track navigation actions, e.g., scrolling through/tabulating a feature list
# #----------------------------------------------------------------------
# output$trackNavs <- renderUI({

#     # process track list
#     trackIds <- plotTrackIds()
#     nTracks <- length(trackIds)    
#     req(nTracks > 0)
#     trackNames <- getTrackNames(trackIds)
#     names(trackNames) <- trackIds

#     # extract any required navigation input rows
#     nNavs <- 0
#     navs <- lapply(trackIds, function(trackId) {
#         track <- tracks[[trackId]]$track
#         hasNav <- isTruthy(track$navigation)
#         if(!hasNav) return(NULL)
#         ui <- tryCatch({ navigation(track, session, id, reference, coord) }, error = function(e) NULL)
#         if(is.null(ui)) return(NULL)
#         nNavs <<- nNavs + 1
#         list(
#             ui = ui,
#             name = trackNames[[trackId]]
#         )
#     })

#     # if needed, populate the trackNav rows
#     if(nNavs == 0) NULL else lapply(navs, function(nav){
#         if(is.null(nav)) return(NULL)
#         class <- nav$ui[[1]]$attribs$class
#         if(!is.null(class) && class == "trackBrowserInput") tags$div(
#             style = "margin-bottom: 8px;",
#             tags$p(
#                 style = "display: inline-block; margin: 30px 10px 5px 10px",
#                 tags$strong(nav$name)
#             ),
#             nav$ui
#         ) else  nav$ui
#     })
# })


#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
# loadingBookmark <- FALSE
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks, fail = FALSE)   
    if(isTruthy(bm)) settings$replace(bm$settings)    
    loadData <- if(isTruthy(bm)) bm else list()
    nRegions <- if(isTruthy(bm) && !is.null(bm$settings$Browser_Options$Number_of_Regions)) bm$settings$Browser_Options$Number_of_Regions$value else 1
    loadSequence <- c(
        list(
            reference$initializeGenome,
            reference$initializeAnnotation
        ),
        lapply(1:nRegions, addRegionCoordinates, FALSE),
        list(
            function(...) dprint("OK TO HERE")
        )
    )

    initializeNextTrackBrowserElement(loadData, loadSequence)

    # trackIds <- bm$outcomes$trackOrder[order(order), trackId]
    # isolate({
    #     lapply(trackIds, function(trackId){
    #         track <- bm$outcomes$tracks[[trackId]]        
    #         tracks[[trackId]] <<- initTrack(track$cssId, trackId, track$type)
    #         tracks[[trackId]]$track$settings$replace(track$settings)
    #         createTrackSettingsObserver(trackId)
    #         if(!is.null(track$items)) tracks[[trackId]]$track$settings$items(track$items)
    #     })
    # })
    # loadingBookmark <<- TRUE
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    loadSourceFile = function(...) NULL,  # enable the trackBrowser for use as the one and only app step
    input = input,
    settings = settings$all_,
    # TODO: add trackOutcomes functionality to allow tracks to send information to downstream app steps
    outcomes = list(
        analysisSetName = reactive({ "trackBrowser" }),
        genome = reference$genomeInput,
        annotation = reference$annotationInput,    
    #     trackOrder = trackOrder,
    #     tracks = reactive({
    #         trackIds <- trackOrder()[, trackId]
    #         x <- lapply(trackIds, function(trackId){
    #             track <- tracks[[trackId]]
    #             list(
    #                 cssId = track$cssId,  
    #                 type = track$type,
    #                 settings = track$track$settings$all_(),
    #                 items = if(is.null(track$track$settings$items)) NULL 
    #                         else track$track$settings$items()
    #             )
    #         })
    #         names(x) <- trackIds
    #         x
    #     }),
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

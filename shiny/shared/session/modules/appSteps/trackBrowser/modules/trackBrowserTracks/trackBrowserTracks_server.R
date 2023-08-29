# trackBrowser server module for selecting and ordering browser tracks
# there is always a single set of tracks selected at a time and applied to all output images
trackBrowserTracksServer <- function(id, browser) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------
defaultTrackTypes <- c(
    "plot_title",
    "chromosome",
    "scale_bar",
    "coordinate_axis",
    "genes"
)
addTrackPrompt <- "-- add a new track --"
duplicateTrackPrompt <- "-- duplicate a track --"
duplicateTrackPromptId <- "0"
suite  <- 'genomex-mdi-tools'
module <- 'trackBrowser'

#----------------------------------------------------------------------
# manage browser track types availabe and in use
#----------------------------------------------------------------------
trackTypes <- list() # key = trackType, value = settings template file path
tracks <- reactiveVal(list())     # key = trackId,   value = browserTrackServer()
nullTrackOrder <- data.table(trackId = character(), order = integer())
trackOrder <- reactiveVal(nullTrackOrder)
orderedTrackIds <- reactive({ # the current track ids, in plotting order
    trackOrder <- trackOrder()
    if(nrow(trackOrder) > 0) trackOrder[order(order), trackId] else character()
})
#----------------------------------------------------------------------
# assemble the track types available to this app
parentAppTrackTypes <- character()
addTrackType <- function(trackType, tracksFolder){
    path <- file.path(tracksFolder, trackType, "settings.yml")
    if(file.exists(path)) trackTypes[[trackType]] <<- path
}
addTrackTypes <- function(dir, classPath, isParentApp = FALSE){
    tracksFolder <- file.path(dir, classPath)
    trackTypes <- list.dirs(tracksFolder, full.names = FALSE, recursive = FALSE)
    if(isParentApp) parentAppTrackTypes <<- sort(trackTypes)
    sapply(trackTypes, addTrackType, tracksFolder)
}
classPath <- "classes/browserTracks"
globalClassPath <- file.path("shiny/shared/global", classPath) 
genomexDirs <- parseExternalSuiteDirs(suite)
addTrackTypes(genomexDirs$suiteDir, globalClassPath) # apps always offer global tracks from genomex-mdi-tools
addTrackTypes(app$DIRECTORY, classPath, isParentApp = TRUE) # apps always offer tracks they define themselves
if(!is.null(browser$options$tracks)) for(trackType in browser$options$tracks){ # apps can additionally offer global tracks declared in the app's config.yml that are...
    if(grepl("//", trackType)){ # ... defined in an external suite, if that suite is actually installed ...
        x <- strsplit(trackType, "//")[[1]]
        trackSuiteDirs <- parseExternalSuiteDirs(x[1])
        if(
            isTruthy(gitStatusData$dependencies[[x[1]]]$loaded) &&
            !is.null(trackSuiteDirs)
        )  addTrackType(x[2], file.path(trackSuiteDirs$suiteDir, globalClassPath)) 
    } else { # ... or in the parent suite of the app
        addTrackType(trackType, file.path(gitStatusData$suite$dir, globalClassPath))
    }
}
#----------------------------------------------------------------------
# initialize available track types
initTrackTypes <- function(){
    names <- names(trackTypes)
    firstTrackTypes <- c(defaultTrackTypes, parentAppTrackTypes)
    sortedTrackTypes <- c(firstTrackTypes, sort(names[!(names %in% firstTrackTypes)]))
    updateSelectInput(
        session, 
        "addTrack", 
        choices = c(
            addTrackPrompt,
            sortedTrackTypes
        ),
        selected = addTrackPrompt
    ) 
}

#----------------------------------------------------------------------
# add, delete and reorder tracks
#----------------------------------------------------------------------
# track identifiers
getTrackId <- function() gsub("( |:|-)", "_", paste(as.character(Sys.time()), sample.int(1e8, 1)))
getTrackNames <- function(trackIds){
    tracks <- tracks()
    sapply(trackIds, function(trackId) getTrackDisplayName(tracks[[trackId]]$track))
}
#----------------------------------------------------------------------
# handle track addition from select input or bookmark
addTrack <- function(trackType, trackId = NULL){
    ns <- if(is.null(trackId)) session$ns else browser$session$ns
    if(is.null(trackId)) trackId <- getTrackId()
    cssId <- paste("track", trackId, sep = "_")
    track <- browserTrackServer(
        cssId = cssId,
        trackId = trackId,
        trackType = trackType,
        settingsFile = trackTypes[[trackType]], # includes any presets defined by the trackType
        presets = browser$options$presets[[trackType]], # add any presets defined by the calling app
        browserInput = input, 
        genome = browser$reference$genome,
        annotation = browser$reference$annotation,
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
        browserTrackUI(ns(cssId), track) # see above; unclear why different ns is required when addTrack is called from init vs. user action

    )
    trackOrder <- trackOrder()
    trackOrder <- rbind(
        trackOrder, 
        data.table(trackId = trackId, order = nrow(trackOrder) + 1)
    )
    trackOrder(trackOrder)
    tracks_ <- tracks()
    tracks_[[trackId]] <- track
    tracks(tracks_)
    trackId
}
observeEvent(input$addTrack, {
    trackType <- input$addTrack
    req(trackType)
    req(trackType != addTrackPrompt)
    updateSelectInput(session, "addTrack", selected = addTrackPrompt) # reset the prompt
    trackId <- addTrack(trackType)
    createTrackSettingsObserver(trackId)
}, ignoreInit = TRUE)
#----------------------------------------------------------------------
# handle track addition from duplication of an existing track
output$duplicateTrack <- renderUI({
    trackIds <- orderedTrackIds()
    req(trackIds)
    tracks <- tracks()
    names(trackIds) <- paste0(
        getTrackNames(trackIds), 
        " (",
        sapply(trackIds, function(x) tracks[[x]]$type),
        ")"
    )
    promptId <- duplicateTrackPromptId
    names(promptId) <- duplicateTrackPrompt
    selectInput(session$ns("duplicateTrackSelect"), NULL, choices = c(promptId, trackIds))
})
observeEvent(input$duplicateTrackSelect, {
    dupTrackId <- input$duplicateTrackSelect
    req(dupTrackId != duplicateTrackPromptId)
    updateSelectInput(session, "duplicateTrackSelect", selected = duplicateTrackPromptId) # reset the prompt
    dupTrack <- tracks()[[dupTrackId]]
    trackId <- addTrack(dupTrack$type)
    tracks()[[trackId]]$track$settings$replace(dupTrack$track$settings$all_())
    createTrackSettingsObserver(trackId)
}, ignoreInit = TRUE)
#----------------------------------------------------------------------
# handle track reordering and deletion
isRankListInit <- FALSE
observeEvent({
    input$trackRankList
    input$deleteRankList
}, {
    if(isRankListInit) {

        # declare the new track order
        currentTrackIds <- trackOrder()[, trackId] 
        newTrackIds <- sapply(strsplit(input$trackRankList, '\\s+'), function(x) x[length(x)])
        nTracks <- length(newTrackIds)
        trackOrder(if(nTracks > 0){
            data.table(trackId = newTrackIds, order = 1:nTracks)
        } else nullTrackOrder)

        # delete tracks as needed
        tracks_ <- tracks()
        for(trackId in currentTrackIds) {
            if(!(trackId %in% newTrackIds)) {
                tracks_[[trackId]] <- NULL
                trackSettingsObservers[[trackId]] <<- NULL
                if(!is.null(trackSettingsUndoId) && trackSettingsUndoId == trackId) trackSettingsUndoId <<- NULL
            }
        }
        tracks(tracks_)
        removeUI(".trackDeleteTarget .browserTrack")
    } else isRankListInit <<- TRUE
}, ignoreInit = TRUE)

#----------------------------------------------------------------------
# handle track settings
#----------------------------------------------------------------------
# undo the last track settings change, intended for disaster recover, not a complete history tracking
trackSettingsObservers <- list()
trackSettingsUndoId <- NULL
createTrackSettingsObserver <- function(trackId){
    track <- tracks()[[trackId]]
    trackSettingsObservers[[trackId]] <<- observeEvent(track$track$settings$all_(), {
        trackSettingsUndoId <<- trackId
        # clearObjectExpansions()
    })
}
observeEvent(input$undoTrackSettings, {
    tracks <- tracks()
    req(trackSettingsUndoId, tracks[[trackSettingsUndoId]])
    tracks[[trackSettingsUndoId]]$track$settings$undo()
})

#----------------------------------------------------------------------
# initialization
#----------------------------------------------------------------------
initialize <- function(jobId, loadData, loadSequence){
    initTrackTypes()
    if(is.null(loadData$outcomes$trackOrder)){
        for(trackType in defaultTrackTypes) addTrack(trackType)
    } else {
        trackIds <- loadData$outcomes$trackOrder[order(order), trackId]
        lapply(trackIds, function(trackId){
            track <- loadData$outcomes$tracks[[trackId]]         
            addTrack(track$type, trackId)
            tracks()[[trackId]]$track$settings$replace(track$settings)
            createTrackSettingsObserver(trackId)
            if(!is.null(track$items)) tracks()[[trackId]]$track$settings$items(track$items)
        })
    }
    initializeNextTrackBrowserElement(loadData, loadSequence)
}

#----------------------------------------------------------------------
# module return value
list(
    tracks = tracks,
    orderedTrackIds = orderedTrackIds,
    trackOrder = trackOrder,
    bookmarkTracks = reactive({
        trackIds <- trackOrder()[, trackId]
        tracks <- tracks()
        x <- lapply(trackIds, function(trackId){
            track <- tracks[[trackId]]
            list(
                # cssId = track$cssId,  
                type = track$type,
                settings = track$track$settings$all_(),
                items = if(is.null(track$track$settings$items)) NULL 
                        else track$track$settings$items()
            )
        })
        names(x) <- trackIds
        x
    }),
    initialize = initialize    
)
#----------------------------------------------------------------------
})}

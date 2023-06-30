#----------------------------------------------------------------------
# server components for the browserTrack widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
browserTrackServer <- function(
    cssId, 
    trackId,    
    trackType, # track class appends "Track" to trackType
    settingsFile,
    browserInput,
    genome, # genome and annotation are zero or single-row data.tables
    annotation,
    size = "m",
    ... # additional arguments passed to settingsServer
) { 
    moduleServer(cssId, function(input, output, session) {    
#----------------------------------------------------------------------

# parse requested track settings
suite  <- 'genomex-mdi-tools'
module <- "browserTrack"
widgetDir <- getWidgetDir(module, suite = suite)
library <- read_yaml(file.path(widgetDir, "settings.yml")) 
template <- list(Track_Options = library$allTracks$Track_Options) # so it is always first
request <- if(is.null(settingsFile)) list() else read_yaml(settingsFile)
if(!is.null(request$include)) for(include in request$include){
    include <- library[[include]]
    if(is.null(include)) next
    for(family in names(include)){
        if(is.null(template[[family]])) template[[family]] <- list()
        for(option in names(include[[family]])) template[[family]][[option]] <- include[[family]][[option]]
    }
}
loadTrackSettings <- function(settings){
    for(family in names(settings)){
        if(is.null(template[[family]])) template[[family]] <<- list()
        for(option in names(settings[[family]])) 
            template[[family]][[option]] <<- settings[[family]][[option]]
    } 
}
if(!is.null(request$shared)) for(sharedFile in request$shared){
    dir <- dirname(settingsFile)
    file <- file.path(dir, sharedFile)
    loadTrackSettings(read_yaml(file))
}
if(!is.null(request$settings)) loadTrackSettings(request$settings)
template$Track_Options$Track_Name <- list( # override and prepend universal Track_Name option
    type = "textInput",
    value = trackType
)

# initialize the track settings
trackClass <- paste0(trackType, "Track")
settings <- settingsServer(
    "settings", 
    cssId, 
    templates = list(template), 
    size = size,
    s3Class = "trackSettings",
    ...
)
class(browserInput) <- append("browserInput", class(browserInput))
constructor <- paste0("new_", trackClass)
track <- get(constructor)(trackId)
class(track) <- unique(append(c("browserTrack", trackClass), class(track)))

# initialize the track list items, e.g., samples or UCSC tracks
if(!is.null(track$items) && track$items) {
    settings$items <- reactiveVal(NULL) # to be filled by items(settings) method
    observeEvent(input$items, {
        reference <- list(
            genome = genome(),
            annotation = annotation()
        )
        tryCatch({ items(track, session, input, reference) }, error = function(e) print(e))
    })
}

# extend the track with required properties, methods, and classes
track$id   <- trackId
track$type <- trackType
track$settings <- settings
track$browser <- browserInput
track$adjustWidth <- function(reference, coord, layout) adjustWidth(track, reference, coord, layout)
track$buildTrack <- function(reference, coord, layout) build(track, reference, coord, layout)
track$buildExpansion <- function(reference, coord, layout) expand(track, reference, coord, layout)

# update the track display name based on settings changes
observeEvent(settings$Track_Options(), {
    html("name", settings$Track_Options()$Track_Name$value)
})

# activate the track level list (e.g., samples) and settings wrappers for click events
selfDestruct <- observe({
    session$sendCustomMessage(
        type = 'updateTrackDialogs',
        message = trackId
    )
    selfDestruct$destroy()  
})

# return value
list(
    id = trackId,
    cssId = cssId,
    type = trackType,
    track = track
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

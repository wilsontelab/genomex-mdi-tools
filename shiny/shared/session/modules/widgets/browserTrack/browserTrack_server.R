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
    presets = NULL, # passed from an app's config.yml to supplement the presets specified by the track itself
    trackSuiteName = NULL, # the tool suite that hosts the track definitions
    ... # additional arguments passed to settingsServer
) { 
    moduleServer(cssId, function(input, output, session) {    
#----------------------------------------------------------------------

# parse requested track settings
suite  <- 'genomex-mdi-tools'
module <- "browserTrack"
browserTrackTypesDir <- "types/browserTrackTypes"
appBrowserTrackTypesDir <- file.path("apps", app$NAME, browserTrackTypesDir)
sharedBrowserTrackTypesDir <- file.path("shared/session", browserTrackTypesDir)
widgetDir <- getWidgetDir(module, suite = suite)
library <- read_yaml(file.path(widgetDir, "settings.yml")) 
template <- list(Track = library$allTracks$Track) # so it is always first
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
    file <- file.path(gitStatusData$suite$dir, "shiny", sharedBrowserTrackTypesDir, sharedFile) # first check this suite's shared trackType settings
    if(!file.exists(file)){ # next, check for files relative to the settings.yml file being processed
        dir <- dirname(settingsFile)
        file <- file.path(dir, sharedFile)
    }
    if(!file.exists(file)) for(suiteDir in app$browser$externalTrackSuites) { # finally check for shared settings files in external track suites
        file <- file.path(suiteDir, "shiny", sharedBrowserTrackTypesDir, sharedFile)
        if(file.exists(file)) break
    }
    if(file.exists(file)) loadTrackSettings(read_yaml(file))
}
if(!is.null(request$trackType)) for(trackType_ in request$trackType){ # allows multiple track types for a given track
    # look for the requested trackType settings ...
    ttFileName <- file.path(trackType_, "settings.yml")    
    trackTypeRequest <- loadExternalYml(gitStatusData$suite$name, file.path(appBrowserTrackTypesDir, ttFileName)) # ... in app ...
    if(is.null(trackTypeRequest)) trackTypeRequest <- loadExternalYml(gitStatusData$suite$name, file.path(sharedBrowserTrackTypesDir, ttFileName)) # ... in the parent suite ...
    if(is.null(trackTypeRequest)) trackTypeRequest <- loadExternalYml("genomex-mdi-tools", file.path(sharedBrowserTrackTypesDir, ttFileName)) # ... in genomex-mdi-tools ...
    if(is.null(trackTypeRequest) && !is.null(trackSuiteName)) trackTypeRequest <- loadExternalYml(trackSuiteName, file.path(sharedBrowserTrackTypesDir, ttFileName)) # ... in the track's host suite ...
    if(!is.null(trackTypeRequest)) {
        if(!is.null(trackTypeRequest$include)) for(include in trackTypeRequest$include){
            include <- library[[include]]
            if(is.null(include)) next
            for(family in names(include)){
                if(is.null(template[[family]])) template[[family]] <- list()
                for(option in names(include[[family]])) template[[family]][[option]] <- include[[family]][[option]]
            }
        }
        if(!is.null(trackTypeRequest$settings)) loadTrackSettings(trackTypeRequest$settings)
    }
}
if(!is.null(request$override)) for(family in names(request$override)){ # allow track to override the defaults of their parent trackTypes
    if(is.null(template[[family]])) next
    for(option in names(request$override[[family]])) {
        if(is.null(template[[family]][[option]])) next
        template[[family]][[option]]$value <- request$override[[family]][[option]]
    }
}
if(!is.null(request$settings)) loadTrackSettings(request$settings)
template$Track$Track_Name <- list( # override and prepend universal Track_Name option
    type = "textInput",
    value = "auto"
)

# initialize the track settings
trackClass <- paste0(trackType, "Track")
class(browserInput) <- append("browserInput", class(browserInput))
constructor <- paste0("new_", trackClass)
track <- get(constructor)(trackId)
class(track) <- unique(append(c("browserTrack", trackClass), class(track)))
if(isTruthy(track$navigation)){
    template$Track$Show_Navigation <- list(
        type = "selectInput",
        choices = c(
            "hide",
            "table_only",
            "navigate",
            if(isTruthy(track$expand)) "expand" else NULL,
            if(isTruthy(track$expand)) "navigate_and_expand" else NULL            
        ),
        value = if(is.logical(track$navigation)) "hide" else track$navigation
    ) 
}
presets <- if(is.null(request$presets)) presets       # presets were defined by app in config.yml
           else if(is.null(presets)) request$presets  # request$presets were defined by track in settings.yml
           else c(presets, request$presets)           # one or both could be concatenated in final track
settings <- settingsServer(
    "settings", 
    cssId, 
    templates = list(template), 
    size = size,
    s3Class = "trackSettings",
    presets = presets,
    ...
)

# initialize the track list items, e.g., samples or UCSC tracks
trackHasItems <- !is.null(track$items) && track$items
if(trackHasItems) {
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
track$adjustWidth    <- function(reference, coord, layout) adjustWidth(track, reference, coord, layout)
track$buildTrack     <- function(reference, coord, layout) build(track, reference, coord, layout)
track$buildExpansion <- function(reference, coord, layout) expand(track, reference, coord, layout)

# update the track display name based on settings and items changes
observeEvent({
    settings$Track()
    if(trackHasItems && !is.null(settings$items)) settings$items()
}, {
    html("name", getTrackDisplayName(track))
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

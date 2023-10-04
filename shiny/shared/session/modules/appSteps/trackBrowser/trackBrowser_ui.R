#----------------------------------------------------------------------
# UI components for the trackBrowser appStep module
#----------------------------------------------------------------------
trackBrowserUI <- function(id, options) {
    module <- 'trackBrowser'
    appStepDir <- getAppStepDir(module)
    ns <- NS(id)
    options <- setDefaultOptions(options, stepModuleInfo$trackBrowser)
    standardSequentialTabItem(

        # page header
        options$longLabel,
        options$leaderText,
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        download = TRUE,
        settings = TRUE,

        # track browser styles
        tags$style(slurpFile(file.path(appStepDir, "trackBrowser.css"))),
        tags$script(slurpFile(file.path(appStepDir, "trackBrowser.js"))),

        # top row of nav inputs
        tags$div(
            style = "white-space: nowrap;",
            class = "trackBrowserSelectizeWrapper",

            # singular genome reference + annotation
            trackBrowserReferenceUI(ns("reference")),

            # one or more sets of region inputs filled by trackBrowserServer
            tags$div(
                style = "display: inline-block;",
                tags$div(
                    class = "trackBrowserCoordinateInputs"
                ),

                # singular additional track-level navigation options; if multiple regions, track must handle dispersal
                tags$div(
                    id = ns("trackNavWrapper"),
                    class = "trackNavWrapper",
                    uiOutput(ns("trackNavs"))
                )
            )
        ),

        # second row of tracks and image outputs
        tags$div(
            style = "white-space: nowrap;",
            class = "trackBrowserSelectizeWrapper",

            # track configuration
            trackBrowserTracksUI(ns("tracks")),

            # the browser output area
            tags$div(
                class = "trackBrowserImages",
                style = "display: inline-block; white-space: nowrap;"
            )
        ),

        # tables for displaying data describing either...
        # ... an object itself ...
        fluidRow(
            class = "objectTableWrapper browserExpansionWrapper",
            style = "display: none;",
            bufferedTableUI(
                id = ns("objectTable"), 
                title = NULL, 
                downloadable = TRUE, 
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE
            )
        ),

        # ... or expansion details about that object
        # row selection in expansionTable might populate expansionUI
        fluidRow(
            class = "expansionTableWrapper browserExpansionWrapper",
            style = "display: none;",
            bufferedTableUI(
                id = ns("expansionTable"), 
                title = NULL, 
                downloadable = TRUE, 
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE
            )  
        ),

        # a place for arbitrary, track-defined UI content based on expand2 or other click actions
        fluidRow(
            uiOutput(ns("expansionUI"))
        )
    )
}

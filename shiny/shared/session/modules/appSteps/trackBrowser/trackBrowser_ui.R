#----------------------------------------------------------------------
# UI components for the trackBrowser appStep module
#----------------------------------------------------------------------

# module ui function
trackBrowserUI <- function(id, options) {

    # initialize namespace
    module <- 'trackBrowser'
    appStepDir <- getAppStepDir(module)
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$trackBrowser)

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        options$leaderText,

        # page header links, uncomment as needed
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

        #----------------------------------------------------------------------
        # top row of browser-level options and navigation
        #----------------------------------------------------------------------
        tags$div(
            style = "white-space: nowrap;",
            class = "trackBrowserSelectizeWrapper",
            trackBrowserReferenceUI(ns("reference")),
            tags$div(
                class = "trackBrowserCoordinateInputs",
                style = "display: inline-block;" # one or more sets of span inputs filled by trackBrowserServer
            ),
            NULL
        ),
        tags$div(
            style = "white-space: nowrap;", #; margin-top: 4px
            class = "trackBrowserSelectizeWrapper",
            trackBrowserTracksUI(ns("tracks")),
        #     #----------------------------------------------------------------------
        #     # the browser output area, with additional track nav options
        #     #----------------------------------------------------------------------
        #     tags$div(
        #         style = "display: inline-block; white-space: normal;",

        #         # additional within-track navigation options, e.g., scrolling through a stack
        #         tags$div(
        #             id = ns("trackNavWrapper"),
        #             class = "trackNavWrapper",
        #             uiOutput(ns("trackNavs"))
        #         ),

        #         # initialization message
        #         tags$p(
        #             id = ns("initMessage"),
        #             style = "padding: 15px;",
        #             tags$strong("Please wait 2 seconds for the browser to initialize.")
        #         ),

        #         # the browser output image
        #         mdiInteractivePlotUI(id = ns("image")),

        #         # image for a track (one at a time) to illustrate expanded details, e.g., on a feature  object click
        #         # aligned with browser tracks to allow expansion to have the same X axis (or not...)
        #         div( 
        #             class = "expansionImageWrapper browserExpansionWrapper",
        #             style = "display: none;",
        #             mdiInteractivePlotUI(id = ns("expansionImage")) 
        #         ),
        #     ),
        #     NULL
        ),

        # # tables for display data describing either an object itself and/or its expansion details
        # # put these below to provide maximum width for complex tables
        # fluidRow(
        #     class = "objectTableWrapper browserExpansionWrapper",
        #     style = "display: none;",
        #     bufferedTableUI(
        #         id = ns("objectTable"), 
        #         title = NULL, 
        #         downloadable = TRUE, 
        #         width = 12,
        #         collapsible = TRUE,
        #         collapsed = FALSE
        #     )
        # ),
        # fluidRow(
        #     class = "expansionTableWrapper browserExpansionWrapper",
        #     style = "display: none;",
        #     bufferedTableUI(
        #         id = ns("expansionTable"), 
        #         title = NULL, 
        #         downloadable = TRUE, 
        #         width = 12,
        #         collapsible = TRUE,
        #         collapsed = FALSE
        #     )  
        # ),

        # # a further place for arbitrary, track-defined UI content based on expand[2] actions
        # fluidRow(
        #     uiOutput(ns("expansionUI"))
        # ),
        NULL
    )
}

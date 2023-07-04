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
    trackRankListId <- "trackRankListGroup"

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
            tags$div(
                class = "trackBrowserInput ucscInput genomeInput",
                popupInputUI(ns('genome'), "Genome")
            ),
            tags$div(
                class = "trackBrowserInput ucscInput annotationInput",
                popupInputUI(ns('annotation'), "Annotation")
            ),
            tags$div(
                class = "trackBrowserInput",
                selectInput(ns('chromosome'), "Chromosome", choices = c(), selectize = FALSE),
            ),
            tags$div(
                class = "trackBrowserInput coordinateInput",
                textInput(ns('start'), "Start", value = 1),
            ),
            tags$div(
                class = "trackBrowserInput coordinateInput",
                textInput(ns('end'), "End", value = 10000),
            ),
            tags$div(
                class = "trackBrowserInput",
                style = "margin: 32px 5px 0px 4px",
                actionLink(ns('all'), "all"),
            ),
            tags$div(
                class = "trackBrowserInput",
                actionButton(ns('zoomOut'), "-"),
            ),
            tags$div(
                class = "trackBrowserInput",
                actionButton(ns('zoomIn'), "+"),
            ),
            tags$div(
                class = "trackBrowserInput",
                textInput(ns('zoomFactor'), "Zoom", value = 10, width = "41px"),
            ),
            tags$div(
                class = "trackBrowserInput",
                actionButton(ns('moveLeft'), "<<"),
            ),
            tags$div(
                class = "trackBrowserInput",
                actionButton(ns('nudgeLeft'), "<"),
            ),
            tags$div(
                class = "trackBrowserInput",
                actionButton(ns('nudgeRight'), ">"),
            ),
            tags$div(
                class = "trackBrowserInput",
                actionButton(ns('moveRight'), ">>"),
            ),
            tags$div(
                class = "trackBrowserInput backButton",
                actionLink(ns('back'), label = "", icon = icon("arrow-left"))
            ),
            tags$div(
                class = "trackBrowserInput",
                textInput(ns('jumpTo'), "Jump To", value = ""),
            ),
            tags$div(
                class = "trackBrowserInput annotationSearchInput",
                popupInputUI(ns('annotationSearch'), label = "", value = NULL, icon = icon("search"), buttonFn = actionLink)
            ),
        ),
        tags$div(
            style = "white-space: nowrap; margin-top: 4px;",

            #----------------------------------------------------------------------
            # the vertical, sortable list of tracks
            #----------------------------------------------------------------------
            tags$div(
                style = "width: 243px;",
                class = "browserContentPanel addTrack",
                tags$div(
                    class = "browserTrack trackDeleteTarget",                    
                    rank_list( 
                        text = "drop track to delete", 
                        labels = NULL, 
                        input_id = ns("deleteRankList"), 
                        css_id = NULL, 
                        options = sortable_options(
                            group = ns(trackRankListId),
                            multiDrag = FALSE 
                        ), 
                        class = "default-sortable"
                    )   
                ),
                tags$div(
                    id = ns("trackList"),
                    class = "browserTrackList",
                    rank_list( 
                        text = "", 
                        labels = NULL, 
                        input_id = ns("trackRankList"), 
                        css_id = NULL, 
                        options = sortable_options(
                            group = ns(trackRankListId),
                            multiDrag = FALSE 
                        ), 
                        class = "default-sortable"
                    )                    
                ),
                uiOutput(ns("duplicateTrack")),
                selectInput(ns("addTrack"), NULL, choices = c()),
                actionLink(ns("undoTrackSettings"), "undo track settings")            
            ),

            #----------------------------------------------------------------------
            # the browser output area, with additional track nav options
            #----------------------------------------------------------------------
            tags$div(
                style = "display: inline-block; white-space: normal;",

                # additional within-track navigation options, e.g., scrolling through a stack
                tags$div(
                    id = ns("trackNavWrapper"),
                    class = "trackNavWrapper",
                    uiOutput(ns("trackNavs"))
                ),

                # initialization message
                tags$p(
                    id = ns("initMessage"),
                    style = "padding: 15px;",
                    tags$strong("Please wait 2 seconds for the browser to initialize.")
                ),

                # the browser output image
                mdiInteractivePlotUI(id = ns("image")),

                # image for a track (one at a time) to illustrate expanded details, e.g., on a feature  object click
                # aligned with browser tracks to allow expansion to have the same X axis (or not...)
                div( 
                    class = "expansionImageWrapper browserExpansionWrapper",
                    style = "display: none;",
                    mdiInteractivePlotUI(id = ns("expansionImage")) 
                ),
            ),
            NULL
        ),

        # tables for display data describing either an object itself and/or its expansion details
        # put these below to provide maximum width for complex tables
        fluidRow(
            class = "objectTableWrapper browserExpansionWrapper",
            style = "display: none;",
            bufferedTableUI(
                id = ns("objectTable"), 
                title = NULL, 
                downloadable = TRUE, 
                width = 12
            )
        ),
        fluidRow(
            class = "expansionTableWrapper browserExpansionWrapper",
            style = "display: none;",
            bufferedTableUI(
                id = ns("expansionTable"), 
                title = NULL, 
                downloadable = TRUE, 
                width = 12
            )  
        )      
    )
}

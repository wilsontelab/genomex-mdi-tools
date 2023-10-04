# trackBrowser UI module for setting plot coordinates in genome via various inputs
# there may be one or multiple distinct sets of coordinates, i.e., regions, plotted by a single browser instance
trackBrowserCoordinatesUI <- function(id, regionI) {
    ns <- NS(id) 
    getLabel <- function(label) if(regionI == 1) label else NULL
    buttonClass <- if(regionI == 1) "trackBrowserInput" else "trackBrowserInput2"
    linkTopMargin <- if(regionI == 1) "32px" else "7px"
    popupInputTopMargin <- if(regionI == 1) "12px" else "0px"
    tags$div(
        id = id,
        tags$div(
            class = "trackBrowserInput",
            selectInput(ns('chromosome'), getLabel("Chromosome"), choices = c()),
        ),
        tags$div(
            class = "trackBrowserInput coordinateInput",
            textInput(ns('start'), getLabel("Start"), value = 1),
        ),
        tags$div(
            class = "trackBrowserInput coordinateInput",
            textInput(ns('end'), getLabel("End"), value = 10000),
        ),
        tags$div(
            class = "trackBrowserInput",
            style = paste("margin:", linkTopMargin, "5px 0px 4px"),
            actionLink(ns('all'), "all"),
        ),
        tags$div(
            class = buttonClass,
            actionButton(ns('zoomOut'), "-"),
        ),
        tags$div(
            class = buttonClass,
            actionButton(ns('zoomIn'), "+"),
        ),
        tags$div(
            class = buttonClass,
            actionButton(ns('moveLeft'), "<<"),
        ),
        tags$div(
            class = buttonClass,
            actionButton(ns('nudgeLeft'), "<"),
        ),
        tags$div(
            class = buttonClass,
            actionButton(ns('nudgeRight'), ">"),
        ),
        tags$div(
            class = buttonClass,
            actionButton(ns('moveRight'), ">>"),
        ),
        tags$div(
            class = "trackBrowserInput backButton",
            style = paste("margin:", linkTopMargin, "5px 0px 5px"),
            actionLink(ns('back'), label = "", icon = icon("arrow-left"))
        ),
        tags$div(
            class = "trackBrowserInput",
            textInput(ns('jumpTo'), getLabel("Jump To"), value = ""),
        ),
        tags$div(
            class = paste0("trackBrowserInput annotationSearchInput", min(regionI, 2)),
            popupInputUI(ns('annotationSearch'), label = NULL, value = NULL, icon = icon("search"), buttonFn = actionLink)
        )
    )
}

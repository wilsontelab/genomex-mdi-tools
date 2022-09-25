#----------------------------------------------------------------------
# UI components for the browserTrack widget module
#----------------------------------------------------------------------

# module ui function
browserTrackUI <- function(id, track) {
    ns <- NS(id)
    tags$div(
        class = "rank-list-item browserTrack",  
        style = "padding: 0px !important;", # must do here since Sortable uses !important in CSS
        tags$div(
            class = "browserTrackLabel",
            tags$p(tags$strong(
                id = ns("name"),
                track$type
            )),
            tags$p(
                track$type,
                tags$span(
                    style = "color: transparent; font-size: 0.1em;",
                    track$id
                )
            )
        ),
        tags$div(
            class = "browserTrackSettings",
            settingsUI(ns("settings"))
        ),
        if(!is.null(track$track$items) && track$track$items) tags$div(
            class = "browserTrackItems",
            actionLink(
                ns('items'), 
                '', 
                icon('list', verify_fa = FALSE)
            )
        ) else ""  
    )
}

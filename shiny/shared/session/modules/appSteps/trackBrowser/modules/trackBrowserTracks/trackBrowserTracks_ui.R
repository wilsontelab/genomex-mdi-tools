# trackBrowser UI module for selecting and ordering browser tracks
# there is always a single set of tracks selected at a time and applied to all output images
trackBrowserTracksUI <- function(id, ...) {
    ns <- NS(id) 
    trackRankListId <- "trackRankListGroup"
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
    )
}

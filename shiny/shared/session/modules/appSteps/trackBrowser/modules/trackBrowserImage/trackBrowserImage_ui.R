# trackBrowser UI module for assembling the composite browser output image
# there may be one or multiple region images rendered by a single browser instance
trackBrowserImageUI <- function(id, regionI) {
    ns <- NS(id) 
    tags$div(
        id = id,
        class = "browserTrackImageWrapper",
        style = "margin-right: 5px; vertical-align: top;",
        
        # the main browser image, all tracks in one genome region
        mdiInteractivePlotUI(id = ns("image")),

        # image for a track (one at a time) to illustrate expanded details for the region, e.g., on a feature object click
        # aligned with browser tracks to allow expansion to have the same X axis (or not...)
        div( 
            class = "expansionImageWrapper browserExpansionWrapper",
            style = "display: none;",
            mdiInteractivePlotUI(id = ns("expansionImage")) 
        )
    )
}

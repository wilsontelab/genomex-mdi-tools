# trackBrowser UI module for loading a reference genome and assoicated annotation
# there is always a single reference genome selected at a time
trackBrowserReferenceUI <- function(id, ...) {
    ns <- NS(id) 
    tags$div(
        style = "display: inline-block; vertical-align: top;",
        tags$div(
            class = "trackBrowserInput ucscInput genomeInput",
            popupInputUI(ns('genome'), "Genome")
        ),
        tags$div(
            class = "trackBrowserInput ucscInput annotationInput",
            popupInputUI(ns('annotation'), "Annotation")
        ),
        NULL
    )
}

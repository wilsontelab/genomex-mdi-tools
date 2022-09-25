#----------------------------------------------------------------------
# ucsc_tracks trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_ucsc_tracksTrack <- function() {
    list(
        click = FALSE,
        hover = FALSE,
        items = TRUE
    )
}

# build method for the S3 class
build.ucsc_tracksTrack <- function(settings, input, reference, coord, layout){
    items <- settings$items()
    ucscTracks <- if(!is.null(items)) lapply(items, function(x) x$display) else list()
    image <- ucscTrackImage(reference, coord, layout, ucscTracks)
    list(image = image)
}

# method for the list icon = track multi-select
items.ucsc_tracksTrack <- function(settings, session, input, reference, track){
    showTrackItemsDialog(
        settings,
        session,
        title = "Select UCSC Tracks",
        itemTypePlural = "Tracks",
        tableData = function() listUcscTracks(reference$genome),
        keyColumn = "track",
        extraColumns = c("shortLabel"),
        options = list(
            display = list(
                type = "selectInput",
                args = list(
                    choices = c("dense", "pack", "full", "squish", "hide"),
                    selected = "dense",
                    width = "90px"                  
                )
            )
        ),
        size = "xl"
    )
}

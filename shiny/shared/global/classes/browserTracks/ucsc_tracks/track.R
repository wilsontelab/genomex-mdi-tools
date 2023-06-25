#----------------------------------------------------------------------
# ucsc_tracks trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_ucsc_tracksTrack <- function(trackId) {
    list(
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE
    )
}

# build method for the S3 class
build.ucsc_tracksTrack <- function(track, reference, coord, layout){
    req(objectHasData(reference$genome))
    items <- track$settings$items()
    ucscTracks <- if(!is.null(items)) lapply(items, function(x) x$display) else list()
    ruler <- getBrowserTrackSetting(track, "Track_Options", "Ruler", default = FALSE)
    image <- ucscTrackImage(reference$genome$genome, coord, layout, ucscTracks, ruler = ruler)
    list(image = image)
}

# method for the list icon = track multi-select
items.ucsc_tracksTrack <- function(track, session, input, reference){
    req(objectHasData(reference$genome))
    showTrackItemsDialog(
        track$settings,
        session,
        title = "Select UCSC Tracks",
        itemTypePlural = "Tracks",
        tableData = function() listUcscTracks(reference$genome$genome),
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

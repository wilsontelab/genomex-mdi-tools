#----------------------------------------------------------------------
# genome_spans trackBrowser track (i.e., a browserTrack)
# genome_spans is a wholly generic wrapper around the trackType of the same name
# it allows a user to visualize any spans files without app-specific developer support
# developers may want to create a similar custom wrapper around the trackType to support object filtering, expansion, etc.
#----------------------------------------------------------------------

# constructor for the S3 class; REQUIRED
new_genome_spansTrack <- function(trackId) {
    list(
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        navigation = TRUE, # whether the track offers a custom, additional row of within-track navigation inputs
        expand = FALSE,
        expand2 = FALSE
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.genome_spansTrack <- function(...) showTrackFilesDialog(...)

# build method for the S3 class; REQUIRED
build.genome_spansTrack <- function(...){
    build.genome_spans_track(..., dataFn = genomex_loadWindowGenomeSpans)
}

# # method for the S3 class to populate one or more trackNav inputs above the browser output
# navigation.genome_spansTrack <- function(...){
#     genomex_genomeSpansNavTable(..., genomex_loadAllGenomeSpans)
# }

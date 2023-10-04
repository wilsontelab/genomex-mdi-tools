#----------------------------------------------------------------------
# genome_spans trackBrowser track (i.e., a browserTrack)
# genome_spans is a wholly generic wrapper around the trackType of the same name
# it allows a user to visualize any spans files without app-specific developer support
# developers may want to create a similar wrapper around the trackType to support customized object filtering, expansion, etc.
#----------------------------------------------------------------------
genome_spansTrackBuffer <- reactiveValues()

# constructor for the S3 class; REQUIRED
new_genome_spansTrack <- function(trackId) {
    list(
        click = TRUE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        navigation = TRUE, 
        expand = FALSE, # TODO: enable expand on plotted spans
        expand2 = FALSE
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.genome_spansTrack <- function(...) showTrackFilesDialog(..., extensions = c("bed", "bed.bgz"))

# build method for the S3 class; REQUIRED
build.genome_spansTrack <- function(...){
    build.genome_spans_track(..., dataFn = genomex_loadWindowGenomeSpans, trackBuffer = genome_spansTrackBuffer)
}

# plot interaction methods for the S3 class
click.genome_spansTrack <- function(track, click, regionI){
    Plot_Spans_As <- getTrackSetting(track, "Spans", "Plot_Spans_As", "scores")
    d <- genome_spansTrackBuffer[[track$id]]
    req(nrow(d) > 0)
    d <- if(Plot_Spans_As != "scores"){
        d <- d[x1 <= click$coord$x & x2 >= click$coord$x]
        req(nrow(d) > 0)
        yDist <- d[, abs(click$coord$y - y)]
        d[yDist == min(yDist)][1]
    }
    if(click$key$ctrl) { 
        app$browser$objectTableData(d) # TODO: instead, could retrieve the original bed file row
    } else {
        app$browser$jumpToCoordinates(regionI, NA, d$x1 - 1, d$x2 + 1)        
    }
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.genome_spansTrack <- function(...){
    genomex_genomeSpansNavTable(..., genomex_loadAllGenomeSpans)
}

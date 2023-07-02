#----------------------------------------------------------------------
# __MODULE_NAME__ trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new___MODULE_NAME__Track <- function(trackId) {
    list(
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = FALSE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        genome = FALSE # whether the track supports a whole genome view
    )
}

# build method for the S3 class; REQUIRED
build.__MODULE_NAME__Track <- function(track, reference, coord, layout){

    # use `req()` to determine if track is ready for plotting
    # req(xyz)

    # use generic methods and any other custom code to determine the track's (dynamic) Y span
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    ylim <- c(0, 1)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "My Track Label", yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        # plotting actions go here, e.g. points(x, y)
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.__MODULE_NAME__Track <- function(track, click){
    # custom actions, use str(click) to explore
}
hover.__MODULE_NAME__Track <- function(track, hover){
    # custom actions, use str(hover) to explore
}
brush.__MODULE_NAME__Track <- function(track, brush){
    # custom actions, use str(brush) to explore
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.__MODULE_NAME__Track <- function(track, session, input, reference){
    showTrackItemsDialog(
        track$settings,
        session,
        title = "Select XYZ",
        itemTypePlural = "Items",
        tableData = function() data.frame(xxx = 1, yyy = 2),
        keyColumn = "xxx",
        extraColumns = c("yyy"),
        options = list(
            XXX = list(
                type = "selectInput", # or textInput, etc.
                args = list(
                    choices = c("aaa", "bbb"),
                    selected = "aaa",
                    width = "50px"                  
                )
            )
        ),
        size = "l"
    )
}

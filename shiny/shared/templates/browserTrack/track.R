#----------------------------------------------------------------------
# __MODULE_NAME__ trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class; REQUIRED
new___MODULE_NAME__Track <- function(trackId) {
    list(
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = FALSE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        expand = FALSE,
        expand2 = FALSE
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
# regionI indicates the region plot the user interacted with, the value must be passed to app$browser$jumpToCoordinates, etc.
click.__MODULE_NAME__Track <- function(track, click, regionI){
    # custom actions, use str(click) to explore
}
hover.__MODULE_NAME__Track <- function(track, hover, regionI){
    # custom actions, use str(hover) to explore
}
brush.__MODULE_NAME__Track <- function(track, brush, regionI){
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

# expand method for the S3 class
# one expansion image can be shown per region, with same width as the main plots
# regionI must be passed to app$browser$expandingTrack
expand.__MODULE_NAME__Track <- function(track, reference, coord, layout, regionI){
    # build a track image the same way as for the main track build
    # typically, need to pass data from "build" to "expand" via a variable scoped to the app
    # expansion tracks may follow the browser genome coordinate, or have an entirely different layout
}

# expand2 method for the S3 class
# a single output at browser bottom, replaced on call to any expand2 function
expand2.__MODULE_NAME__Track <- function(track, browser, selectedRowData){
    # return a tagList of arbitrary UI content in response to an expansion table click
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
# only one navigation set is shown per track, your navigation should decide how to handle multiple regions
navigation.__MODULE_NAME__Track <- function(track, session, id, browser){

    # initialize the trackNavs, including input observers to handle user actions
    # initTrackNav will fail silenty if setting Track/Show_Navigation is not set or =="hide"
    navName1 <- initTrackNav(track, session, "inputName1", function(inputValue1){
        # do work as needed based on the input value, e.g., make a call to app$browser$jumpToCoordinates()
    })
    navName2 <- initTrackNav(track, session, "tableName2") # table reactive functions are provided below
    # etc.

    # as needed, create a reactive with the data to enumerate in the trackNavs
    trackNavData <- reactive({
        data.table( # use track settings or other information to populate dynamically
            chrom = c("chr1","chr2","chr3"),
            start = 1e8,
            end = 2e8
        )
    })

    # return the associated navigation UI elements
    tagList(
        trackNavInput(
            track, 
            session, 
            navName1, # the name as provided by initTrackNav
            radioButtons, # any valid Shiny UI input function
            label = "Input Name 1", # all further arguments must be named and are passed to the UI function
            choices = c("foo", "bar"),            
            selected = "bar", # the default value
            inline = TRUE,
            width = "150px"
            # add other argument to pass to the Shiny UI input function
        ),
        trackNavTable(
            track, 
            session, 
            browser$id,
            navName2, # the name as provided by initTrackNav
            tableData = trackNavData, # populate a table based on track settings, etc.
            actionFn = function(selectedRow){
                req(selectedRow)
                d <- trackNavData()[selectedRow]
                # do other work as needed based on the input value (e.g., open a modal, navigate)
                # you should honor the value of trackNavCanNavigate(track) and trackNavCanExpand(track)
                # the easiest way is to use handleTrackNavTableClick(regionI, track, chrom, start, end, expandFn)
            }
            # add other argument to pass to bufferedTableServer, but omit "selection" and "width"
        )
        # etc.
    )
}

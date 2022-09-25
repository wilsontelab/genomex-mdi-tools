#----------------------------------------------------------------------
# __TRACK_TYPE__ trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_testTrack <- function() {
    list(      
        click = TRUE,
        hover = FALSE,
        items = FALSE
    )
}

# build method for the S3 class
build.testTrack <- function(settings, input){

    d <- settings$Data_Options()

    coord <- coordinates(input)
    x <- coord$start:coord$end
    y <- x * 0.1 + rnorm(coord$width, sd = d$Noise_Factor$value)
    ylim <- ylim(settings, y)

    list(
        mar = list(top = 0, bottom = 0),
        args = list(
            settings = settings,
            input = input,
            x = x,
            y = y,
            ylim = ylim,
            ylab = ylab(settings, "Test Label")        
        ),
        ylim = ylim, # last values used for click and hover handling
        height = height(settings, 1)
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click or track$hover is TRUE, above
click.testTrack <- function(track, x, y){

    dmsg()
    dprint(track$id)
    dprint(paste("x = ", x))
    dprint(paste("y = ", y))

}
hover.testTrack <- function(track, x, y){

    dmsg()
    dprint(track$id)
    dprint(paste("x = ", x))
    dprint(paste("y = ", y))

}

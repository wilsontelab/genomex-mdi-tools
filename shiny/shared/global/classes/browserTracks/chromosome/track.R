#----------------------------------------------------------------------
# chromosome_simple trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_chromosomeTrack <- function() {
    list(
        click = TRUE,
        hover = FALSE,
        items = FALSE
    )
}

# build method for the S3 class
chromosomeTrackData <- list()
chromBandColors <- list(
    acen    = 'darkgoldenrod4',
    stalk   = 'darkgoldenrod4',
    gneg    = 'grey95',
    gpos25  = 'grey75',
    gpos33  = 'grey66',
    gpos50  = 'grey50',
    gpos66  = 'grey33',
    gpos75  = 'grey25',
    gpos100 = 'grey5',
    gvar    = 'grey5'
) 
build.chromosomeTrack <- function(settings, input, reference, coord, layout){
    req(reference$genome)
    req(input$chromosome)
    track <- settings$get("Plot_Options", "Feature_Track")
    features <- if(track == "none") NULL else getChromosomeFeatures(reference$genome, track)
    chroms <- listUcscChromosomes(reference$genome)
    chromSize <- chroms[chromosome == input$chromosome, size]
    padding <- padding(settings, layout)
    height <- height(settings, 0.25) + padding$total
    unit <- parseUnit(scaleUnit(settings), chromSize)
    ylim <- c(0, 1)
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n",
            ylim = ylim,  ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") 
        getX <- function(pos){ pos / chromSize * coord$width + coord$start }

        # whole chromosome
        y1 <- 0.25
        y2 <- 0.75
        rect(coord$start, y1, coord$end, y2, col="grey", border=NA)

        # major features
        if(!is.null(features)){
            features <- features[chrom == input$chromosome]
            if(nrow(features) > 0) rect(
                getX(features$chromStart), 
                y1, 
                getX(features$chromEnd), 
                y2, 
                col = if(track == "gap") "black" 
                      else sapply(features$gieStain, function(x) chromBandColors[[x]]), 
                border = NA
            )
        }

        # window location rectangle
        x1 <- getX(coord$start)
        x2 <- getX(coord$end)
        y1 <- 0.1
        y2 <- 0.9
        rect(x1, y1, x2, y2, col = NA, border = "red", lwd = 1.25)

        # chromosome text label
        tmp <- par('xpd')
        par(xpd = TRUE, cex = 1.15)
        size <- paste0(round(chromSize / unit$multiplier, 0), " ", unit$unit)
        lab <- paste0(coord$chromosome, ", ", size)
        text(coord$start, 0.5, lab, pos = 2)
        par(xpd = tmp, cex = 1)
    })
    chromosomeTrackData <<- list(
        chromSize = chromSize,
        coord = coord
    )
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click or track$hover is TRUE, above
click.chromosomeTrack <- function(track, x, y){
    d <- chromosomeTrackData
    req(d$coord)
    app$browser$center((x - d$coord$start) / d$coord$width * d$chromSize)  
}

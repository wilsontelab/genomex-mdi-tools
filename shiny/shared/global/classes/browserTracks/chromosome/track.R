#----------------------------------------------------------------------
# chromosome_simple trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_chromosomeTrack <- function(trackId) {
    list(
        click = TRUE,
        hover = FALSE,
        brush = TRUE,
        forceBrush = TRUE, # special setting for genome-scale tracks to override standard no-key brush
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
    gvar    = 'grey5',
    odd     = 'grey40',
    even    = 'grey80',
    odd1    = "#4444dd",
    even1   = "#bbbbee",
    odd2    = "#dd4444",
    even2   = "#eebbbb"
) 
build.chromosomeTrack <- function(track, reference, coord, layout){
    req(objectHasData(reference$genome))
    genome <- reference$genome$genome
    req(coord$chromosome)
    featuresTrack <- getBrowserTrackSetting(track, "Plot_Options", "Feature_Track", "cytoBand") # i.e., a UCSC track
    compositeGenomes <- listCompositeGenomes(reference)
    isAll <- coord$chromosome == "all"
    isCompositeGenome <- coord$chromosome %in% compositeGenomes
    features <- if(isAll) getChromosomeSizes(reference$genome, reference$metadata)
                else if(isCompositeGenome) getChromosomeSizes(reference$genome, reference$metadata)[genome == coord$chromosome]
                else if(featuresTrack == "none") NULL 
                else getUcscChromosomeFeatures(reference, coord$chromosome, featuresTrack)
    if(isAll || isCompositeGenome){
        chromSize <- features[, sum(chromEnd - chromStart + 1)]
        chromStart_ <- features[, min(chromStart)] # since 2nd composite genome doesn't start at 1
    } else {
        chroms <- listSourceChromosomes(reference$genome, metadata = reference$metadata)
        chromSize <- chroms[chromosome == coord$chromosome, size]
        chromStart_ <- 1
    }
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total
    unit <- parseUnit(scaleUnit(track), chromSize)
    ylim <- c(0, 1)
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n",
            ylim = ylim,  ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") 
        getX <- function(pos){
            (pos - chromStart_ + 1)  / chromSize * as.numeric(coord$width) + as.numeric(coord$start)
        }

        # whole chromosome (or genome)
        y1 <- 0.25
        y2 <- 0.75
        rect(coord$start, y1, coord$end, y2, col="grey", border=NA)

        # major features
        if(!is.null(features)){
            if(nrow(features) > 0) rect(
                getX(features$chromStart), 
                y1, 
                getX(features$chromEnd), 
                y2, 
                col = if(featuresTrack == "gap") "black" 
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
        placement <- sideLabelPlacement(layout, coord)
        size <- paste0(round(chromSize / unit$multiplier, 0), " ", unit$unit)
        lab <- paste0(coord$chromosome, ", ", size)
        labelWidth <- strwidth(lab, font = layout$pointsize, units = "in")
        if(labelWidth > placement$marginWidth * 0.9) lab <- coord$chromosome
        text(placement$coord, 0.5, lab, pos = placement$side)
        par(xpd = tmp, cex = 1)
    })
    chromosomeTrackData <<- list(
        chromStart = chromStart_,
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
click.chromosomeTrack <- function(track, click, regionI){
    d <- chromosomeTrackData
    req(d$coord)
    app$browser$center(regionI, d$chromStart + (click$coord$x - d$coord$start) / d$coord$width * d$chromSize)  
}
brush.chromosomeTrack <- function(track, brush, regionI){
    d <- chromosomeTrackData
    req(d$coord)
    getX <- function(x) d$chromStart + (x - d$coord$start) / d$coord$width * d$chromSize
    app$browser$jumpToCoordinates(
        regionI,
        d$coord$chromosome, 
        getX(brush$coord$x1), 
        getX(brush$coord$x2), 
        strict = TRUE
    )  
}

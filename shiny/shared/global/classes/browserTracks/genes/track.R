#----------------------------------------------------------------------
# genes trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_genesTrack <- function(trackId) {
    list(
        click = TRUE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = FALSE,
        navigation = FALSE # whether the track offers a custom, additional row of within-track navigation inputs
    )
}

# build method for the S3 class; REQUIRED
genesTrackBuffer <- list()
build.genesTrack <- function(track, reference, coord, layout){
    req(objectHasData(reference$annotation))
    padding <- padding(track, layout)
    lwd <- getBrowserTrackSetting(track, "Plot_Options", "Line_Weight", 2)
    Max_Genes_BP <- getBrowserTrackSetting(track, "Plot_Options", "Max_Genes_BP", 50000000)
    req(coord$width <= Max_Genes_BP)
    Max_Exons_BP <- getBrowserTrackSetting(track, "Plot_Options", "Max_Exons_BP", 1000000)
    Min_Gene_Size_Label <- getBrowserTrackSetting(track, "Plot_Options", "Min_Gene_Size_Label", 25000)
    Force_Gene_Labels <- getBrowserTrackSetting(track, "Plot_Options", "Force_Gene_Labels", "")
    Force_Gene_Labels <- strsplit(trimws(Force_Gene_Labels), '(,|\\s+)')[[1]]
    palette <- switch(
        getBrowserTrackSetting(track, "Plot_Options", "Palette", "green/red"),
        "blue/orange" = list("+" = "blue", "-" = "orange3"),
        list("+" = "green4", "-" = "red3")
    )
    height <- getBrowserTrackSetting(track, "Track", "Height", 0.8)
    height <- height + padding$total # or set a known, fixed height in inches
    ylim <- c(0, 4)
    genes <- getRegionGenes(reference$genome, reference$annotation, coord, force = FALSE) %>%
             setUcscFeatureEndpoints(reference$annotation)
    genesTrackBuffer[[track$id]] <<- genes

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = ylab(track, ""), yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        # process genes
        req(nrow(genes) > 0)
        geneY   <- genes[, ifelse(strand == "-", 1, 2)]
        geneCol <- genes[, unlist(palette[strand])]
        rect(
            genes[, start], 
            geneY, 
            genes[, end], 
            geneY + 1, 
            col = NA, 
            border = geneCol, 
            lty = 1, 
            lwd = lwd
        )

        # only add details if view is sufficiently narrow
        if(coord$width < Max_Exons_BP){

            # add exons
            exons <- expandUcscExons(genes, reference$annotation)
            exonY   <- exons[, ifelse(strand == "-", 1, 2)]
            exonCol <- exons[, unlist(palette[strand])]
            rect(
                exons[, start], 
                exonY, 
                exons[, end], 
                exonY + 1, 
                col = exonCol, 
                border = NA
            )
        }

        # add gene name labels
        I <- genes[, {
            size <- end - start
            I <- size >= Min_Gene_Size_Label
            if(sum(I) <= 5) I
            else {
                sizeThreshold <- rev(sort(size[I]))[5]
                size >= sizeThreshold | name2 %in% Force_Gene_Labels
            }
        }]       
        labels <- genes[I, name2]
        labelHalfWidths <- strwidth(labels) / 2 * 1.2
        x <- genes[I, start + (end - start) / 2]
        y <- genes[I, ifelse(strand == "-", 0.5, 3.5)]
        if(sum(I) > 0) text(
            pmin(coord$end - labelHalfWidths, pmax(coord$start + labelHalfWidths, x)), 
            y = y, 
            labels = labels, 
            cex = 1.15, 
            col = geneCol[I]
        )
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
click.genesTrack <- function(track, click){
    strand_ <- if(click$coord$y < 2) "-" else "+"
    genes <- genesTrackBuffer[[track$id]][strand == strand_ & 
                                          start <= click$coord$x & 
                                          end >= click$coord$x]
    req(nrow(genes) > 0)
    app$browser$jumpToCoordinates(genes[1, chrom], genes[, min(start)], genes[, max(end)])
}
hover.genesTrack <- function(track, hover){
    # custom actions
}

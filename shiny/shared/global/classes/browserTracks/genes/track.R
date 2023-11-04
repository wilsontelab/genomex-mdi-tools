#----------------------------------------------------------------------
# genes trackBrowser track (i.e., a browserTrack)
# provides click navigation and Ctrl-click gene card expansion
#----------------------------------------------------------------------

# constructor for the S3 class
new_genesTrack <- function(trackId) {
    list(
        click = TRUE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = FALSE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        expand = FALSE,
        expand2 = FALSE
    )
}

# build method for the S3 class; REQUIRED
genesTrackBuffer <- list()
build.genesTrack <- function(track, reference, coord, layout){
    reference <- fillCompositeAnnotation(reference, coord$chromosome) # does nothing for a standard genome
    if(!is.null(reference$chromosome)) coord$chromosome <- reference$chromosome
    req(objectHasData(reference$annotation), coord$chromosome != "all")
    genes <- getRegionGenes(reference, coord, force = FALSE) %>%
             setUcscFeatureEndpoints(reference)
    genesTrackBuffer[[track$id]] <<- list(reference = reference, genes = genes)

    padding <- padding(track, layout)
    lwd <- getBrowserTrackSetting(track, "Plot_Options", "Line_Weight", 2)
    Max_Genes_BP <- getBrowserTrackSetting(track, "Plot_Options", "Max_Genes_BP", 50000000)
    if(!isTruthy(coord$width <= Max_Genes_BP)) return(trackInfo(track, coord, layout, "window too wide to plot genes"))
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

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = ylab(track, ""), yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        # Watson/Crick strand lines
        abline(h = 1:2 + 0.5, col = "grey50", lwd = 0.5)
        placement <- sideLabelPlacement(layout, coord) 
        tmp <- par('xpd')
        par(xpd = TRUE, cex = 1.25)
        text(placement$coord, 2.5, "+", pos = placement$side)
        text(placement$coord, 1.5, "\U2012", pos = placement$side)
        par(xpd = tmp, cex = 1)

        # process genes
        if(nrow(genes) == 0) trackNoData(coord, ylim, "no genes in window", y = 0.5) else {
            geneY   <- genes[, ifelse(strand == "-", 1, 2)]
            geneCol <- genes[, unlist(palette[strand])]
            rect(
                genes[, start], 
                geneY, 
                genes[, end], 
                geneY + 1, 
                col = "white", 
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
        }
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.genesTrack <- function(track, click, regionI){

    # get the target gene
    strand_ <- if(click$coord$y < 2) "-" else "+"
    x <- genesTrackBuffer[[track$id]]
    genes <- x$genes[
        strand == strand_ & 
        start <= click$coord$x & 
        end >= click$coord$x
    ][1]
    req(nrow(genes) > 0)

    # use NCBI gene search as a highly generalized, species agnostic way of retrieving gene information
    # (with the additional benefit that the site allows being embedded in an iFrame, not all do)
    if(click$key$ctrl) { 
        genome <- x$reference$genome
        if("name2" %in% names(genes)) app$browser$expansionUI(tags$iframe(
            src = paste0(
                "https://www.ncbi.nlm.nih.gov/gene?term=(", 
                if(is.null(genome$scientificName)) genome$organism else genome$scientificName, 
                "[Organism]) AND ", genes[1, name2], "[Gene Name]"
            ),
            style = "width: 100%; height: 1000px; margin-top: 5px;"
        ))  

    # execute simple click navigation to gene limits
    } else {
        app$browser$jumpToCoordinates(regionI, NA, genes[, min(start)], genes[, max(end)])
    }
}

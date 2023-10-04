# load and plot genome spans for genome_spans trackType, and potentially others

# load windows spans, for plotting
genomex_loadWindowGenomeSpans <- function(track, reference, coord, bedFile, ...){
    bgzipForTabix(bedFile) %>% 
    getCachedTabix(create = FALSE, index = FALSE, force = FALSE) %>% 
    getTabixRangeData(coord) %>% 
    parseTabixBedFormat()
}

# load and cache all spans, for nav table
genomex_cacheAllGenomeSpans <- function(bgzFile, permanent = TRUE, from = "disk", create = "asNeeded"){
    sessionCache$get(
        'spans', 
        key = digest(bgzFile), 
        permanent = permanent, 
        from = from, 
        create = create, 
        createFn = function(...) {
            startSpinner(session, message = "loading spans")
            cmd <- paste("zcat ", bgzFile)
            spans <- fread(cmd = cmd)
            stopSpinner(session)
            spans
        }
    )$value     
}
genomex_loadAllGenomeSpans <- function(bedFile){
    bgzipForTabix(bedFile) %>% 
    genomex_cacheAllGenomeSpans() %>% 
    parseTabixBedFormat()
}

# generic BED track nav function showing all available columns up to, but not past, column 6 (strand)
genomex_genomeSpansNavTable <- function(track, session, id, browser, loadFn, ...){
    navTableName <- initTrackNav(track, session, "navTable") # table reactive functions are provided below  
    outCols <- c("source","chrom","start","end","name","score","strand")
    trackNavData <- reactive({
        spanFiles <- names(track$settings$items())
        spans <- do.call(rbind, lapply(spanFiles, function(spanFile){
            source_ <- gsub(".bed", "", gsub(".bed.bgz", "", basename(spanFile)))
            spans <- loadFn(spanFile)
            spans[, source := source_]
            for(col in outCols) if(!(col %in% names(spans))) spans[[col]] <- NA
            spans[, .SD, .SDcols = outCols]
        }))
        req(nrow(spans) <= 10000)
        spans
    })
    handleRowClick <- function(selectedRow){
        req(selectedRow)
        span <- trackNavData()[selectedRow]
        if(span$end - span$start < 10){
            span$start <- span$start - 10
            span$end   <- span$end   + 10
        }
        handleTrackNavTableClick(NULL, track, span$chrom, span$start + 1, span$end)
    }
    tagList(
        trackNavTable(
            track, 
            session, 
            browser$id,
            navTableName, # the name as provided by initTrackNav
            tableData = trackNavData, # populate a table based on track settings, etc.
            actionFn = handleRowClick
        )
    )  
}

# parse standardized input to create consistently formatted XY plots of various types
getPackedSpans 
buildSpanTrackImage <- function(track, coord, layout,
                               itemsList, itemNames, itemData,
                               stranded = TRUE, allowNeg = FALSE, ylab = NULL, ylim = NULL, yaxt = "s",
                               dataFamily = "Data", yAxisFamily = "Y_Axis", hLines = FALSE){
    nItems <- length(itemNames)

    # set the plot frame
    padding <- padding(track, layout)
    height <- height(track, 2) + padding$total # or set a known, fixed height in inches
    ylab <- ylab(track,
        if(is.null(ylab)) {
            if(nItems == 1) itemNames else ""
        }
        else if(is.function(ylab)) ylab()
        else ylab
    )

    # set the dynamic Y-axis
    if(is.null(ylim)) ylim <- getDynamicYLim(track, yAxisFamily, itemsList, stranded, allowNeg)
    yPadding <- diff(ylim) * 0.05
    ylim <- c(ylim[1] - yPadding, ylim[2] + yPadding)

    # make the plot
    mai <- NULL    
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = ylab, yaxt = yaxt,
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        # add horizontal rules if requested  
        zeroLine(track)  
        if(hLines) hLines(track, ylim)   

        # plot spans as lines
        I <- 1:nItems
        palette <- CONSTANTS$palettes[[getTrackSetting(track, dataFamily, "Spans_Color_Palette", "plotly")]]   
        for(i in I){ 
            dd <- itemData[source == itemNames[i]]
            if(nrow(dd) == 0) next
            dstr(dd)
            if(stranded) for(strand_ in c("+", "-")){
                ddd <- dd[strand == strand_]
                if(nrow(ddd) == 0) next
                if(strand_ == "-" && ddd[1, y] > 0) ddd[, y := -y]
                plotSpans(track, ddd, color = palette[[i]], family = dataFamily)
            } else {
                plotSpans(track, dd, color = palette[[i]], family = dataFamily)
            }
        }            

        # add a legend
        trackLegend(track, coord, ylim, legend = itemNames, pch = 19, cex = 1, col = unlist(palette[I]))
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}

# trackBrowser server module for setting plot coordinates in genome via various inputs
# there may be one or multiple distinct sets of coordinates plotted by a single browser instance
trackBrowserCoordinatesServer <- function(
    id, 
    browserId, browserOptions, browserInput, browserSettings,
    reference, regionI
) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# browser navigation history, to support the back button
#----------------------------------------------------------------------
coordinateHistory <- list()
pushCoordinateHistory <- function(coord){ # technically, these act as shift and unshift...
    coord$strict <- TRUE # any object padding was applied when the window was loaded
    coord$history <- FALSE # i.e., don't re-record this view
    coordinateHistory <<- c(list(coord), coordinateHistory)
    if(length(coordinateHistory) > 100) coordinateHistory <<- coordinateHistory[1:100]
}
popCoordinateHistory <- function(){
    req(length(coordinateHistory) > 1) # current window is in slot 1
    coordinateHistory <<- coordinateHistory[2:length(coordinateHistory)]    
    coordinateHistory[[1]]
}
observeEvent(input$back,  { 
    do.call(jumpToCoordinates, popCoordinateHistory())
}, ignoreInit = TRUE)  

#----------------------------------------------------------------------
# browser navigation support functions
#----------------------------------------------------------------------
isStrict <- function(){
    x <- browserSettings$Browser_Options()$Strict_Coordinates$value
    if(is.null(x)) FALSE else x
}
jumpToCoordinates <- function(chromosome, start, end, strict = FALSE, history = TRUE, then = NULL){ # arguments are strict coordinates
    req(chromosome, start, end)
    start <- as.integer64(start)
    end   <- as.integer64(end)
    if(start > end){
        tmp <- start
        start <- end
        end <- tmp
    }
    if(!isStrict() && !strict){
        padding <- (end - start + 1) * 0.05
        start <- as.integer64(start - padding)
        end   <- as.integer64(end   + padding)
    }
    genome <- reference$genome()
    chromosomeSize <- getTargetChromosomeSize(chromosome) 
    if(start < 1) start <- 1
    if(end > chromosomeSize) end <- chromosomeSize
    clearObjectExpansions()
    if(history) pushCoordinateHistory(list(chromosome = chromosome, start = start, end = end))
    updateSelectInput(session, "chromosome", selected = chromosome)
    updateTextInput(session, "start", value = as.character(start))
    updateTextInput(session, "end",   value = as.character(end))
    if(!is.null(then)) thenObserver <- observeEvent(browserIsDone(), {
        setTimeout(then, delay = 100)
        thenObserver$destroy()
    }, ignoreInit = TRUE)
}
doZoom <- function(exp){
    start  <- as.integer64(input$start)
    end    <- as.integer64(input$end)
    factor <- browserSettings$get("Browser_Options", "Zoom_Factor", 10)
    width  <- (end - start + 1)
    center <- start + width / 2
    newHalfWidth <- (width * factor ** exp) / 2
    jumpToCoordinates(
        input$chromosome, 
        center - newHalfWidth, 
        center + newHalfWidth, 
        strict = TRUE
    )
}
observeEvent(input$zoomOut, { doZoom( 1) }, ignoreInit = TRUE)
observeEvent(input$zoomIn,  { doZoom(-1) }, ignoreInit = TRUE)
doMove <- function(factor, direction){
    start  <- as.integer64(input$start)
    end    <- as.integer64(input$end)
    width  <- (end - start + 1)
    increment <- width * factor * direction
    jumpToCoordinates(
        input$chromosome, 
        start + increment, 
        end + increment, 
        strict = TRUE
    )
}
observeEvent(input$moveLeft,   { doMove(1,    -1) }, ignoreInit = TRUE)
observeEvent(input$nudgeLeft,  { doMove(0.05, -1) }, ignoreInit = TRUE)
observeEvent(input$nudgeRight, { doMove(0.05,  1) }, ignoreInit = TRUE)
observeEvent(input$moveRight,  { doMove(1,     1) }, ignoreInit = TRUE)
observeEvent(input$all, { 
    jumpToCoordinates(input$chromosome, 1, currentChromosomeSize(), strict = TRUE) 
}, ignoreInit = TRUE)
center <- function(x){
    coord <- coordinates(input)
    halfWidth <- coord$width / 2
    jumpToCoordinates(
        coord$chromosome, 
        x - halfWidth, 
        x + halfWidth, 
        strict = TRUE
    )
}

#----------------------------------------------------------------------
# jumpTo coordinates
#----------------------------------------------------------------------
checkJumpChrom <- function(chrom_){  
    chrom <- chromosomeSizes()[name == chrom_, .(name, size)]
    req(nrow(chrom) == 1)
    chrom
}
checkJumpCenter <- function(chrom, center_){
    req(!grepl('\\D', center_))
    center <- as.integer64(center_)
    req(center >= 1, center <= chrom$size)
    center
}
checkJumpStart <- function(start_){
    req(!grepl('\\D', start_))
    start <- as.integer64(start_)
    req(start >= 1)
    start
}
checkJumpEnd <- function(chrom, start, end_){
    req(!grepl('\\D', end_))
    end <- as.integer64(end_)
    req(end > start, end <= chrom$size)
    end
}
checkJumpGene <- function(gene){
    genome <- reference$genome()
    annotation <- reference$annotation()
    req(gene, objectHasData(genome), objectHasData(annotation))
    gene <- getGene(genome, annotation, gene) %>% setUcscFeatureEndpoints(annotation)
    req(nrow(gene) == 1)
    list(
        chromosome = gene$chrom, 
        start      = gene$start, 
        end        = gene$end,
        strict = FALSE
    )
}
executeJumpTo <- function(action = NULL){
    req(action)
    updateTextInput(session, "jumpTo", value = "")
    updateTextInput(session, "end", value = "")
    do.call(jumpToCoordinates, action)
}
observeEvent(input$jumpTo,  { 
    req(input$jumpTo)    
    jumpTo <- trimws(input$jumpTo)
    req(jumpTo)
    parts <- strsplit(jumpTo, '(:|,|-|\\s+)')[[1]]
    action <- tryCatch({
        switch(
            length(parts),
            { # a single word, assumed to be an annotation feature name, i.e., a gene
                checkJumpGene(parts[1])
            }, { # two parts, assumed to be chrom + center, at strict current width
                chrom <- checkJumpChrom(parts[1])
                center <- checkJumpCenter(chrom, parts[2])
                coord <- coordinates(input)
                halfWidth <- coord$width / 2
                list( 
                    chromosome = chrom$name, 
                    start = center - halfWidth, 
                    end   = center + halfWidth,
                    strict = TRUE
                )
            }, { # three parts, assumbed to be a region, i.e., chrom + start + end
                chrom <- checkJumpChrom(parts[1])
                start <- checkJumpStart(parts[2])
                end <- checkJumpEnd(chrom, start, parts[3])
                list(
                    chromosome = chrom$name, 
                    start = start, 
                    end   = end # subjected to Browser_Options$Strict_Coordinates
                )
            } 
        )
    }, error = function(e) NULL)
    executeJumpTo(action)
}, ignoreInit = TRUE)

#----------------------------------------------------------------------
# gene search popup navigation
#----------------------------------------------------------------------
genes <- reactive({
    genome <- reference$genome()
    annotation <- reference$annotation()
    req(nrow(genome) > 0, nrow(annotation) > 0)
    getGenomeGenes(genome, annotation) %>% setUcscFeatureEndpoints(annotation) 
})
geneI <- reactiveVal(NULL)
genesTable <- bufferedTableServer(
    "genes",
    id,
    input,
    tableData = reactive({ genes()[, .(name2, chrom, strand, start, end)] }),
    select = geneI,
    options = list( searchDelay = 0 )
)
annotationSearchInput <- popupInputServer(
    "annotationSearch", 
    "Navigate to a Gene", 
    callback = function(...){
        rowI <- genesTable$selectionObserver()
        if(is.na(rowI)) {
            geneI(NULL)
            NA
        } else {
            geneI(rowI)
            genes()[rowI]$name2 # the popup's return value, a character gene name
        }
    },
    updateLabel = FALSE,
    bufferedTableUI(session$ns("genes"))
    # ,
    # actionLink(session$ns("getUcscGenomes"), "Reload from UCSC")
)
observeEvent(annotationSearchInput(),  { 
    gene <- annotationSearchInput()
    action <- checkJumpGene(gene)
    executeJumpTo(action)
}, ignoreInit = TRUE)  

#----------------------------------------------------------------------
# initialization
#----------------------------------------------------------------------
initialize <- function(jobId, loadData, loadSequence = NULL){
    chromosomes <- reference$chromosomes()
    loadData$chromosome  <- if(is.null(loadData$input$chromosome)) chromosomes[1] else loadData$input$chromosome
    loadData$start       <- if(is.null(loadData$input$start)) 1   else loadData$input$start
    loadData$end         <- if(is.null(loadData$input$end)) 10000 else loadData$input$end
    updateSelectInput(session, 'chromosome', choices = chromosomes, selected = loadData$chromosome)
    updateTextInput(session,   'start',      value    = loadData$start)
    updateTextInput(session,   'end',        value    = loadData$end)
    pushCoordinateHistory(list(chromosome = loadData$chromosome, start = loadData$start, end = loadData$end))
    if(!is.null(loadSequence)) initializeNextTrackBrowserElement(loadData, loadSequence)
}

# module return value
list(
    initialize = initialize,
    jumpToCoordinates = jumpToCoordinates,
    center = center,
    NULL
)
#----------------------------------------------------------------------
})}

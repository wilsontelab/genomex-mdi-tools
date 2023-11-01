# trackBrowser server module for loading a reference genome and assoicated annotation
# there is always a single reference genome selected at a time
trackBrowserReferenceServer <- function(id, browser) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------
defaults <- list(
    source = "UCSC",
    genome = "hg38",
    annotation = "wgEncodeGencodeBasicV41"
)

#----------------------------------------------------------------------
# fill cascading genome > annotations, chromosomes
#----------------------------------------------------------------------
# genome
ucscGenomes   <- reactiveVal( listUcscGenomes() )
customGenomes <- reactiveVal( listCustomGenomes() )
genomeMetadata <- reactiveVal(list())
observeEvent(input$reloadGenomes, {
    ucscGenomes(  listUcscGenomes(  force = TRUE))
    customGenomes(listCustomGenomes(force = TRUE))
})
ucscGenomeI   <- reactiveVal(NULL)
customGenomeI <- reactiveVal(NULL)
ucscGenomesTable <- bufferedTableServer(
    "ucscGenomes",
    id,
    input,
    tableData = ucscGenomes, # don't change this, needs to retain source
    select = ucscGenomeI,
    options = list( searchDelay = 0 )
)
customGenomesTable <- bufferedTableServer(
    "customGenomes",
    id,
    input,
    tableData = customGenomes,
    select = customGenomeI,
    options = list( searchDelay = 0 )
)
genome <- reactive({
    genomeInput <- tryCatch({ get("genomeInput") }, error = function(e) data.table())
    if(!is.null(genomeInput)) genomeInput() else data.table()
})
genomeInput <- popupInputServer(
    "genome", 
    "Set Working Genome", 
    callback = function(...){
        source <- input$genomeInputTabPanel
        if(source == "UCSC"){
            rowI <- ucscGenomesTable$selectionObserver()
            if(is.na(rowI)) {
                ucscGenomeI(NULL)
                customGenomeI(NULL)
                genomeMetadata(list())
                data.table()
            } else {
                ucscGenomeI(rowI)
                customGenomeI(NULL)
                genomeMetadata(list())
                ucscGenomes()[rowI] # the popup's return value, one genome row
            }
        } else {
            rowI <- customGenomesTable$selectionObserver()
            if(is.na(rowI)) {
                ucscGenomeI(NULL)
                customGenomeI(NULL)
                genomeMetadata(list())
                data.table()
            } else {
                ucscGenomeI(NULL)
                customGenomeI(rowI)
                genome <- customGenomes()[rowI]
                genomeMetadata(loadCustomGenomeMetadata(genome$genome))
                genome # the row carries a source column telling whether it is UCSC or Custom
            }
        }
    },
    labelCol = "genome",
    tags$p(tags$strong({
        x <- genome()
        if(!objectHasData(x)) "" else paste0("current selection = ", x$genome, " (", x$source, ")")
    })),
    tabsetPanel(
        id = session$ns("genomeInputTabPanel"),
        tabPanel("UCSC",   bufferedTableUI(session$ns("ucscGenomes"))),
        tabPanel("Custom", bufferedTableUI(session$ns("customGenomes")))
    ),
    actionLink(session$ns("reloadGenomes"), "Reload Genomes")
)
#----------------------------------------------------------------------
# annotation
annotations <- reactiveVal( data.table() )
setAnnotationsFromGenome <- function(genome, force = FALSE){
    req(genome, nrow(genome) > 0)
    annotations( 
        if(genome$source == "UCSC") listUcscAnnotations(  genome = genome$genome, force = force)
                               else listCustomAnnotations(genome = genome$genome, force = force)
    )
}
observeEvent(input$reloadAnnotations, {
    setAnnotationsFromGenome(genome(), force = TRUE)
})
annotationI <- reactiveVal(NULL)
annotationsTable <- bufferedTableServer(
    "annotations",
    id,
    input,
    tableData = annotations,
    select = annotationI,
    options = list( searchDelay = 0 )
)
annotation <- reactive({
    annotationInput <- tryCatch({ get("annotationInput") }, error = function(e) data.table())
    if(!is.null(annotationInput)) annotationInput() else data.table()
})
annotationInput <- popupInputServer(
    "annotation", 
    "Set Working Annotation", 
    callback = function(...){
        rowI <- annotationsTable$selectionObserver()
        if(is.na(rowI)) {
            annotationI(NULL)
            data.table()
        } else {
            annotationI(rowI)
            annotations()[rowI] # the popup's return value, one annotation row
        }
    },
    labelCol = "track",
    tags$p(tags$strong({
        x <- annotation()
        if(!objectHasData(x)) "" else paste("current selection = ", x$track)
    })),
    bufferedTableUI(session$ns("annotations")),
    actionLink(session$ns("reloadAnnotations"), "Reload Annotations"),
    active = reactive({ 
        genome <- genome()
        isTruthy(genome) && nrow(genome) > 0 
    })
)
#----------------------------------------------------------------------
# chromosomes
chromosomes <- reactiveVal(NULL)
observeEvent({
    genome()
    genomeMetadata()
}, {
    annotationInput(data.table())
    genome <- genome()
    req(genome, nrow(genome) > 0)
    setAnnotationsFromGenome(genome)
    chromosomes(c(
        listCanonicalChromosomes(genome), 
        listCompositeGenomes(reference = list(genome = genome, metadata = genomeMetadata())),
        "all"
    ))
})

#----------------------------------------------------------------------
# additional metadata
#----------------------------------------------------------------------
genomeSize <- reactive({
    genome <- genome()
    req(genome, nrow(genome) > 0)
    getGenomeSize(genome)
})
chromosomeSizes <- reactive({ # size of the currently active chromosome (not the one we may be switching to...)
    genome <- genome()
    req(genome, nrow(genome) > 0)
    getChromosomeSizes(genome, genomeMetadata())
})
getChromosomeStart <- function(chrom_){
    req(chrom_)
    if(isProperChromosome(chrom_) || chrom_ == "all") 1
    else chromosomeSizes()[genome == chrom_, min(chromStart)]
}
getChromosomeEnd <- function(chrom_){
    req(chrom_)
    if(isProperChromosome(chrom_) || chrom_ == "all") getChromosomeSize(chrom_)
    else chromosomeSizes()[genome == chrom_, max(chromEnd)]
}
getChromosomeSize <- function(chrom_){
    req(chrom_)
    if(isProperChromosome(chrom_)) chromosomeSizes()[name == chrom_, size] 
    else if(chrom_ == "all") genomeSize()
    else chromosomeSizes()[genome == chrom_, sum(size)]
}

#----------------------------------------------------------------------
# initialization
#----------------------------------------------------------------------
initializeGenome <- function(jobId, loadData, loadSequence){
    if(is.null(loadData$outcomes$genome)) loadData$outcomes$genome <- ucscGenomes()[genome == defaults$genome]
    setAnnotationsFromGenome(loadData$outcomes$genome)
    genomeInput(loadData$outcomes$genome)
    genomeMetadata(
        if(loadData$outcomes$genome$source == "UCSC") list()
        else loadCustomGenomeMetadata(loadData$outcomes$genome$genome)
    )
    doNextLoadSequenceItem(loadData, loadSequence)
}
initializeAnnotation <- function(jobId, loadData, loadSequence){
    if(is.null(loadData$outcomes$annotation)) loadData$outcomes$annotation <- annotations()[track == defaults$annotation]
    annotationInput(loadData$outcomes$annotation)
    doNextLoadSequenceItem(loadData, loadSequence)
}
#----------------------------------------------------------------------
# module return value
list(
    genomeMetadata = genomeMetadata,
    genomeInput = genomeInput,
    genome = genome,
    annotationInput = annotationInput,
    annotation = annotation,
    chromosomes = chromosomes,
    genomeSize = genomeSize,
    chromosomeSizes = chromosomeSizes,
    getChromosomeStart = getChromosomeStart,
    getChromosomeEnd = getChromosomeEnd,
    getChromosomeSize = getChromosomeSize,
    initializeGenome = initializeGenome,
    initializeAnnotation = initializeAnnotation
)
#----------------------------------------------------------------------
})}

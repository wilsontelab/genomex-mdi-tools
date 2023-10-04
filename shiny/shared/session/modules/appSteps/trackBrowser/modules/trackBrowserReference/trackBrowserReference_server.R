# trackBrowser server module for loading a reference genome and assoicated annotation
# there is always a single reference genome selected at a time
trackBrowserReferenceServer <- function(id, browser) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------
defaults <- list(
    genome = "hg38",
    annotation = "wgEncodeGencodeBasicV41"
)

#----------------------------------------------------------------------
# fill cascading genome > annotations, chromosomes
#----------------------------------------------------------------------
# genome
genomes <- reactiveVal( listUcscGenomes() )
observeEvent(input$getUcscGenomes, {
    genomes( listUcscGenomes(force = TRUE) )
})
genomeI <- reactiveVal(NULL)
genomesTable <- bufferedTableServer(
    "genomes",
    id,
    input,
    tableData = genomes,
    select = genomeI,
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
        rowI <- genomesTable$selectionObserver()
        if(is.na(rowI)) {
            genomeI(NULL)
            data.table()
        } else {
            genomeI(rowI)
            genomes()[rowI] # the popup's return value, one genome row
        }
    },
    labelCol = "genome",
    tags$p(tags$strong({
        x <- genome()
        if(!objectHasData(x)) "" else paste("current selection = ", x$genome)
    })),
    bufferedTableUI(session$ns("genomes")),
    actionLink(session$ns("getUcscGenomes"), "Reload from UCSC")
)
#----------------------------------------------------------------------
# annotation
annotations <- reactiveVal( data.table() )
observeEvent(input$getUcscAnnotations, {
    genome <- genome()
    req(genome)
    genomes( listUcscAnnotations(genome = genome$genome, force = TRUE) )
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
    actionLink(session$ns("getUcscAnnotations"), "Reload from UCSC"),
    active = reactive({ nrow(genome()) > 0 })    
)
#----------------------------------------------------------------------
# chromosomes
chromosomes <- reactiveVal(NULL)
observeEvent(genome(), {
    genome <- genome()
    req(nrow(genome()) > 0)
    annotations(listUcscAnnotations(genome$genome))
    annotationInput(data.table())            
    chromosomes(c(listCanonicalChromosomes(genome$genome), "all"))     
})

#----------------------------------------------------------------------
# additional metadata
#----------------------------------------------------------------------
genomeSize <- reactive({
    genome <- genome()
    req(nrow(genome()) > 0)
    getGenomeSize(genome$genome)  
})
chromosomeSizes <- reactive({ # size of the currently active chromosome (not the one we may be switching to...)
    genome <- genome()
    req(nrow(genome()) > 0)
    getChromosomeSizes(genome$genome)
})
getChromosomeSize <- function(chrom_){
    req(chrom_)
    if(chrom_ == "all") genomeSize()
    else chromosomeSizes()[name == chrom_, size]    
}

#----------------------------------------------------------------------
# initialization
#----------------------------------------------------------------------
initializeGenome <- function(jobId, loadData, loadSequence){
    if(is.null(loadData$outcomes$genome)) loadData$outcomes$genome <- genomes()[genome == defaults$genome]
    annotations(listUcscAnnotations(loadData$outcomes$genome$genome))
    genomeInput(loadData$outcomes$genome)
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
    genomeInput = genomeInput,
    genome = genome,
    annotationInput = annotationInput,
    annotation = annotation,
    chromosomes = chromosomes,
    genomeSize = genomeSize,
    chromosomeSizes = chromosomeSizes,
    getChromosomeSize = getChromosomeSize,
    initializeGenome = initializeGenome,
    initializeAnnotation = initializeAnnotation
)
#----------------------------------------------------------------------
})}

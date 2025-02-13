#----------------------------------------------------------------------
# pahts and support functions for reading genome data from the user-defined custom genomes
#----------------------------------------------------------------------
customGenomesDirs <- c(
    file.path(serverEnv$RESOURCES_DIR, "genomes"),
    file.path(serverEnv$RESOURCES_DIR, "genomes/custom")
)
nullCustomGenome <- data.table(
    source = character(), 
    genome = character(), 
    organism = character(), 
    scientificName = character(), 
    description = character()
)
nullCustomAnnotation <- data.table(
    track = character(), 
    type = character(), 
    group = character(), 
    shortLabel = character(), 
    longLabel = character()
)
customGenomeCols <- names(nullCustomGenome)
getCustomGenomeFile <- function(genome, extension){
    filename <- paste(genome, extension, sep = ".")
    file <- file.path(customGenomesDirs[1], genome, filename)
    if(!file.exists(file)) file <- file.path(customGenomesDirs[2], genome, filename)
    file
}

#----------------------------------------------------------------------
# retrieve the metadata about a specific custome genome
#----------------------------------------------------------------------
loadCustomGenomeMetadata <- function(genome){
    metadataFile <- getCustomGenomeFile(genome, "yml")
    if(!file.exists(metadataFile)) return(NULL)
    read_yaml(metadataFile)
}

#----------------------------------------------------------------------
# custom genome metadata
#----------------------------------------------------------------------
isCompositeGenome <- function(reference){
    reference$genome$source == "Custom" &&
    isTruthy(reference$metadata$composite) && 
    isTruthy(reference$metadata$compositeType)
}
isCompositeGenome2 <- function(metadata){
    isTruthy(metadata$composite) && 
    isTruthy(metadata$compositeType)
}
getCustomCompositeType <- function(reference){
    if(isCompositeGenome(reference)) reference$metadata$compositeType 
    else NA
}
getCustomCompositeDelimiter <- function(metadata){
    x <- metadata$compositeDelimiter
    if(is.null(x)) "_" else x
}

#----------------------------------------------------------------------
# retrieve genome-level metadata from user-defined custom genomes
#----------------------------------------------------------------------
listCustomGenomes <- function(...){
    candidateGenomes <- list.dirs(
        path = customGenomesDirs[1], 
        full.names = FALSE, 
        recursive = FALSE
    )
    if(dir.exists(customGenomesDirs[2])) candidateGenomes <- c(candidateGenomes, list.dirs(
        path = customGenomesDirs[2], 
        full.names = FALSE, 
        recursive = FALSE
    ))
    if(length(candidateGenomes) == 0) return(nullCustomGenome)
    do.call(rbind, lapply(candidateGenomes, function(genome){
        fai <- getCustomGenomeFile(genome, "fa.fai")
        if(!file.exists(fai)) return(NULL)
        metadata <- loadCustomGenomeMetadata(genome)
        if(is.null(metadata)) return(NULL)
        metadata$source <- "Custom"
        if(is.null(metadata$genome)) metadata$genome <- genome
        for(col in customGenomeCols) if(is.null(metadata[[col]])) metadata[[col]] <- NA
        as.data.table(metadata[customGenomeCols])
    }))
}
listCustomAnnotations <- function(genome, ...){ # this is used to populate the annotation popupInput
    nullCustomAnnotation # at present, only support is for composite UCSC genomes
}                        # in future, could allow custom genomes to define their own annotations in their source folder
listCompositeGenomes <- function(reference){ # the set of individual genomes that comprise the composite ...
    isComposite <- reference$genome$source == "Custom" &&
                   isTruthy(reference$metadata$composite) && 
                   isTruthy(reference$metadata$compositeSources)
    if(isComposite) reference$metadata$compositeSources else list() # ... here as a named list of genome-specific annotations
}
listCompositeGenomeNames <- function(reference) names(listCompositeGenomes(reference)) # ... here as just a vector of genome names
listCustomChromosomes <- function(reference){
    fai <- getCustomGenomeFile(reference$genome$genome, "fa.fai") # read custom chromosomes from the requisite fai index
    fai <- fread(fai)[, 1:2]
    setnames(fai, c("chromosome","size"))
    fai
}
fillCompositeAnnotation_ <- function(reference, compositeSource, genome, chromosome){
    if(compositeSource$source == "UCSC"){
        reference$genome <- getUscsGenome(genome)
        reference$annotation <- as.data.table(compositeSource)
        reference$metadata <- list()
        reference$chromosome <- chromosome
    } # TODO: implement user-built custom annotations
    reference
}
fillCompositeAnnotation <- function(reference, chromosome){ # when user select one custom genome as chromosome, retrieve its dedicated annotation
    compositeSources <- listCompositeGenomes(reference)
    if(objectHasData(compositeSources)) {
        compositeSource <- compositeSources[[chromosome]]
        if(objectHasData(compositeSource)) {
            reference <- fillCompositeAnnotation_(reference, compositeSource, genome = chromosome, chromosome = "all")
        } else {
            delimiter <- getCustomCompositeDelimiter(reference$metadata)
            chrom <- strsplit(chromosome, delimiter)[[1]]
            if(isTruthy(chrom[2])){
                compositeSource <- compositeSources[[chrom[2]]]
                if(objectHasData(compositeSource)) {
                    reference <- fillCompositeAnnotation_(reference, compositeSource, genome = chrom[2], chromosome = chrom[1])
                }
            }
        }
    }
    reference
}

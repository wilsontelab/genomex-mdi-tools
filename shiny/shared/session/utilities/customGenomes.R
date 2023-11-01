#----------------------------------------------------------------------
# pahts and support functions for reading genome data from the user-defined custom genomes
#----------------------------------------------------------------------
customGenomesDir <- file.path(serverEnv$RESOURCES_DIR, "genomes/custom")
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
    file.path(customGenomesDir, genome, filename)
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
getCustomCompositeType <- function(reference){
    isComposite <- reference$genome$source == "Custom" &&
                   isTruthy(reference$metadata$composite) && 
                   isTruthy(reference$metadata$compositeType)
    if(isComposite) reference$metadata$compositeType
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
    if(!dir.exists(customGenomesDir)) return(nullCustomGenome)
    candidateGenomes <- list.dirs(
        path = customGenomesDir, 
        full.names = FALSE, 
        recursive = FALSE
    )
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
listCustomAnnotations <- function(genome, ...){
    nullCustomAnnotation # at present, only support is for composite UCSC genomes
}                        # in future, could allow custom genomes to define their own annotations in their source folder
listCompositeGenomes <- function(reference){
    isComposite <- reference$genome$source == "Custom" &&
                   isTruthy(reference$metadata$composite) && 
                   isTruthy(reference$metadata$compositeSources)
    if(isComposite) names(reference$metadata$compositeSources) else character()
}
listCustomChromosomes <- function(reference){
    fai <- getCustomGenomeFile(reference$genome$genome, "fa.fai")
    fai <- fread(fai)[, 1:2]
    setnames(fai, c("chromosome","size"))
    fai
}

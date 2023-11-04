# generic function for either UCSC or custom genome retrieval
getChromGenes <- function(reference, chromosome = "all", force = FALSE){
    fileName <- paste("genes", chromosome, sep = ".")
    file <- getUcscRdsFile(file.path("genomes", reference$genome$genome, "annotations", reference$annotation$track), fileName)
    file <- loadPersistentFile(
        file = file,
        force = force,
        unlink = force,
        ttl = CONSTANTS$ttl$year,
        create = function(file){
            transcripts <- if(!is.null(reference$annotation$source) && reference$annotation$source == "Custom"){
                NULL # not implemented yet; when it is, custom annotations must conform the UCSC annotation table formats
            } else {
                getUcscChromTranscripts(reference$genome$genome, reference$annotation, chromosome = chromosome, force = force)
            }
            startSpinner(session, message = paste("tabulating genes"))
            saveRDS(aggregateUcscTranscriptsToGenes(reference$annotation, transcripts), file = file)
            stopSpinner(session)
        }
    )
    persistentCache[[file]]$data 
}
getGenomeGenes <- function(reference, force = FALSE){
    getChromGenes(reference, chromosome = "all", force = force)
}
getRegionGenes <- function(reference, coord, force = FALSE){
    genes <- getChromGenes(reference, chromosome = coord$chromosome, force = force)
    if(isBigGenePred(reference$annotation)){
        genes[chromStart <= coord$end & coord$start <= chromEnd]
    } else {
        genes[txStart    <= coord$end & coord$start <= txEnd]
    }
}
getGene <- function(reference, gene, force = FALSE){
    genes <- getGenomeGenes(reference, force = force)
    genes[toupper(name2) == toupper(gene)]
}

#----------------------------------------------------------------------
# support functions for incorporating UCSC browser tracks into MDI trackBrowsers
# or otherwise retrieving genomic data from UCSC resources
#----------------------------------------------------------------------
ucscTimeout <- 30 # seconds
ucscGet <- function(url, ...){
    tryCatch({
        httr::GET(url, ..., httr::timeout(ucscTimeout))
    }, error = function(e){
        message()
        message("UCSC RETRIEVAL FAILURE")
        message(url)
        print(e)
        stopSpinner(session)
        req(FALSE)
    })
}

#----------------------------------------------------------------------
# rendering UCSC track images
# 'genome' arguments take a single character genome name
#----------------------------------------------------------------------
ucscRenderTracksUrl <- "http://genome.ucsc.edu/cgi-bin/hgRenderTracks"
ucscPixelsPerChar <- list(
    "8"   = 6,  # determined empirically
    "10"  = 6,
    "12"  = 7,
    "14"  = 8,
    "18"  = 10, # 18 never used for screen or print given Font_Size offerings
    "24"  = 14,
    "34"  = 18  # used for Font_Size 8,10 and 12 since UCSC max font size is 34
)
ucscLabelExtraPixels <- 16 # fixed extra space+lines added in UCSC label region not included in hgt.labelWidth
ucscFontSizes <- c(6, 8, 10, 12, 14, 18, 24, 34)
adjustLayoutForUcsc <- function(layout){

    # increment the UCSC font size - similar point size numbers as R are too small
    # account for scaling of UCSC request to account for print resolution
    layout$uscsPointsize <- {
        i <- which(ucscFontSizes >= layout$pointsize * layout$printMultiplier + 2)
        if(length(i) == 0) 34 else ucscFontSizes[i[1]]
    }

    # resolve the final widths
    pixelsPerChar <- ucscPixelsPerChar[[as.character(layout$uscsPointsize)]]
    pixels <- layout$mai$left * layout$dpi
    layout$uscsLabelWidthChars <- max(3, ceiling(pixels / pixelsPerChar))
    pixels <- layout$uscsLabelWidthChars * pixelsPerChar + ucscLabelExtraPixels
    layout$mai$left <- pixels / layout$dpi
    layout$mai$right <- max(1 / layout$dpi, layout$mai$right) # UCSC browser adds 1px right margin
    layout$plotWidth <- layout$browserWidth - layout$mai$left - layout$mai$right
    layout
}
ucscTrackImage <- function(genome, coord, layout, tracks = list(), ruler = FALSE){
    req(tracks, length(tracks) > 0)    
    startSpinner(session, message = "getting UCSC tracks")

    # assemble the get url query
    query <- list(
        db = genome,
        position = coord$region,
        guidelines = "off",
        hgt.labelWidth = layout$uscsLabelWidthChars,
        pix = (layout$mai$left + layout$plotWidth) * layout$dpi + 1,
        textSize = layout$uscsPointsize,
        hideTracks = 1
    )

    # allow user to determine whether to use MDI or UCSC coordinate rulers
    if(!isTruthy(ruler)) query$ruler = "hide" 

    # enforce vertical track ordering (does this always work? might be buggy on UCSC end)
    trackOrder <- as.list(as.character(1:length(tracks)))
    names(trackOrder) <- paste0(names(tracks), "_imgOrd")

    # one MDI track might stack multiple UCSC tracks
    ucsc <- ucscGet(ucscRenderTracksUrl, query = c( 
        query, 
        tracks, # name = track, value = full|dense|pack|hide 
        trackOrder   
    ))

    # TODO: use a better null error image if UCSC fails
    if(is.null(ucsc) || !isTruthy(ucsc$status_code) || ucsc$status_code != 200) {
        print(ucsc$url)
        ucsc <- ucscGet(ucscRenderTracksUrl)
    }
    image <- magick::image_read(ucsc$content)
    info  <- magick::image_info(image)
    image <- magick::image_extent( # add space for the right legend region and top and bottom padding
        image, 
        geometry = paste0(layout$width, "x", info$height + 2 * 10), 
        gravity = "West", 
        color = "white"
    )
    stopSpinner(session)
    image
}
# hubUrl=<url>
# <trackName>.heightPer=<number> - sets the height of the a bigWig track in pixels

#----------------------------------------------------------------------
# URLs and support functions for reading genome data from the UCSC API
#----------------------------------------------------------------------
ucscResourceDir <- file.path(serverEnv$RESOURCES_DIR, "UCSC")
ucscApiPrefix  <- "https://api.genome.ucsc.edu/"
ucscListPrefix <- paste0(ucscApiPrefix, "list/")
ucscGetPrefix  <- paste0(ucscApiPrefix, "getData/")
createUcscDirectory <- function(relPath){
    dir <- if(is.null(relPath)) ucscResourceDir 
           else file.path(ucscResourceDir, relPath)
    if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    dir
}
getUcscRdsFile <- function(relPath, fileName){
    dir <- createUcscDirectory(relPath)
    file.path(dir, paste(fileName, "rds", sep = "."))    
}

#----------------------------------------------------------------------
# retrieve genome-level metadata
# 'genome' arguments take a single character genome name
#----------------------------------------------------------------------
listUcscGenomes <- function(force = FALSE){
    file <- getUcscRdsFile("genomes", "genomes")
    file <- loadPersistentFile(
        file = file,
        force = force,
        unlink = force,
        ttl = CONSTANTS$ttl$month,
        create = function(file){
            startSpinner(session, message = "getting UCSC genomes")
            target <- "ucscGenomes"
            url <- paste0(ucscListPrefix, target)
            ucsc <- ucscGet(url, httr::accept_json()) 
            ucsc <- httr::content(ucsc, type = "application/json")[[target]]
            dt <- data.table(do.call(rbind, ucsc))
            dt[, genome := names(ucsc)]
            cols <- c("organism", "genome", "scientificName", "description", "orderKey")
            dt <- dt[, .SD, .SDcols = cols][, lapply(.SD, unlist, recursive = FALSE)] 
            dt[, source := "UCSC"]
            saveRDS(
                dt[order(organism, orderKey), .(source, genome, organism, scientificName, description)], 
                file = file
            )
            stopSpinner(session)
        }
    )
    persistentCache[[file]]$data
}
getUscsGenome <- function(genome_){
    genomes <- listUcscGenomes()
    genomes[genome == genome_]
}
listUcscTracks <- function(genome, force = FALSE){
    file <- getUcscRdsFile(file.path("genomes", genome), "tracks")
    file <- loadPersistentFile(
        file = file,
        force = force,
        unlink = force,
        ttl = CONSTANTS$ttl$month,
        create = function(file){
            startSpinner(session, message = paste("getting", genome, "tracks"))
            target <- "tracks"
            url <- paste0(ucscListPrefix, target)
            ucsc <- ucscGet(url, httr::accept_json(), 
                              query = list(genome = genome, trackLeavesOnly = 1))
            ucsc <- httr::content(ucsc, type = "application/json")[[genome]]
            cols <- c("type", "group", "shortLabel", "longLabel")   
            dt <- as.data.table(t(sapply(ucsc, function(track){
                sapply(cols, function(col) if(is.null(track[[col]])) NA else track[[col]] )
            })))
            dt <- dt[, track := names(ucsc)]
            saveRDS(
                dt[, .SD, .SDcols = c("track", cols)], 
                file = file
            )   
            stopSpinner(session)
        }
    )
    persistentCache[[file]]$data
}
listUcscAnnotations <- function(genome, force = FALSE){
    listUcscTracks(genome, force)[type == "genePred" | type == "bigGenePred"]
}
listUcscChromosomes <- function(genome, force = FALSE){
    file <- getUcscRdsFile(file.path("genomes", genome), "chromosomes")
    file <- loadPersistentFile(
        file = file,
        force = force,
        unlink = force,
        ttl = CONSTANTS$ttl$month,
        create = function(file){
            startSpinner(session, message = paste("getting", genome, "chromosomes"))
            target <- "chromosomes"
            url <- paste0(ucscListPrefix, target)
            ucsc <- ucscGet(url, httr::accept_json(), query = list(genome = genome)) 
            ucsc <- httr::content(ucsc, type = "application/json")[[target]]
            saveRDS(
                data.table(chromosome = names(ucsc), size = unlist(ucsc)), 
                file = file
            )   
            stopSpinner(session)           
        }
    )
    persistentCache[[file]]$data
}
#----------------------------------------------------------------------
# generic track get function; caller decides whether to work by chromosome or whole genome
# 'genome' and 'track' arguments take a single character entity name
#----------------------------------------------------------------------
getUcscTrackTable <- function(genome, track, chromosome = NULL,
                              col.names = NULL, colClasses = NULL,   
                              force = FALSE, ttl = "year", 
                              maxItemsOutput = -1){ # default always tries to get (hopefully) all data
    fileName <- if(is.null(chromosome)) track else paste(track, chromosome, sep = ".")
    file <- getUcscRdsFile(file.path("genomes", genome, "tracks", track), fileName)
    file <- loadPersistentFile(
        file = file,
        force = force,
        unlink = force,
        ttl = CONSTANTS$ttl[[ttl]],
        create = function(file){
            startSpinner(session, message = paste("downloading", genome, "track", track))
            url <- paste0(ucscGetPrefix, "track")
            ucsc <- ucscGet(url, httr::accept_json(), query = list(
                genome = genome, 
                track  = track,
                chrom  = chromosome,
                maxItemsOutput = maxItemsOutput
            ))
            if(!(ucsc$status_code %in% c(200, 206))) { # code 206 mean partial content
                print(paste("UCSC failed:", ucsc$status_code, genome, track, chromosome))
                stopSpinner(session)
                return(data.table()) # failure returns an empty data.table, but does not save it
            }

            startSpinner(session, message = paste("parsing", genome, "track", track))
            ucsc <- httr::content(ucsc, type = "application/json")[[track]]

            startSpinner(session, message = paste("tabulating", genome, "track", track))
            parseFeatures <- function(rows){
                if(length(rows) == 0) return("")
                paste0(
                    paste(
                        sapply(rows, function(row){
                            paste(
                                sapply(col.names, function(col) if(is.null(row[[col]])) "NA" else row[[col]]),
                                collapse = "\t"
                            )
                        }),
                        collapse = "\n"             
                    ), 
                    "\n"
                ) 
            }
            dt <- fread(
                text = paste(
                    if(!is.null(names(ucsc))){ # handle the more typical format with features itemized by chromosome
                        sapply(names(ucsc), function(chrom) parseFeatures(ucsc[[chrom]]) )
                    } else { # handle the bigGenePred format, with a single list of all features
                        parseFeatures(ucsc)
                    },
                    collapse = ""
                ), 
                colClasses = colClasses,
                col.names = col.names
            )

            startSpinner(session, message = paste("saving", genome, "track", track))
            saveRDS(dt, file = file)   
            stopSpinner(session)           
        }
    )
    persistentCache[[file]]$data 
}

#----------------------------------------------------------------------
# retrieve low-resolution chromosomeone metadata, for chromosome track
#----------------------------------------------------------------------
getUcscChromosomeFeatures <- function(reference, chromosome, track){
    genome <- if(reference$genome$source == "UCSC"){
        reference$genome$genome
    } else if({
        x <- getCustomCompositeType(reference)
        isTruthy(x) && x == "UCSC"
    }) {
        delimiter <- getCustomCompositeDelimiter(reference$metadata)
        x <- strsplit(chromosome, delimiter)[[1]] # expects composite chromosome names in form ChromDelimiterGenome, e.g., chr3_hg38
        chromosome <- x[1]
        x[2]
    } else {
        return(NULL)
    }
    availableTracks <- listUcscTracks(genome)
    req(availableTracks)
    if(!(track %in% availableTracks$track)) return(NULL)
    tryCatch({
        dt <- getUcscTrackTable(
            genome, 
            track,
            col.names = switch(
                track,
                cytoBandIdeo = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
                cytoBand     = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
                gap          = c("chrom", "chromStart", "chromEnd", "type", "size")
            ),
            colClasses = switch(
                track,
                cytoBandIdeo = c("character", "integer", "integer", "character", "character"),
                cytoBand     = c("character", "integer", "integer", "character", "character"),
                gap          = c("character", "integer", "integer", "character", "integer")
            )
        )
        dt[chrom == chromosome]
    }, error = function(e) NULL)
}

#----------------------------------------------------------------------
# standardized subsets of UCSC columns used by genomeBrowser
# post-processing functions below, or the consuming code, must distinguish and handle the genePred formats
#----------------------------------------------------------------------
ucscGenePredTypes <- list(
    genePred = c(
        "name",  # unique name for the transcript isoform or other feature
        "name2", # gene name shared across isoforms
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
        "exonFrames"           
    ),
    bigGenePred = c(
        "name", 
        "name2", 
        "chrom",
        "strand",
        "chromStart", # synonym for "txStart"
        "chromEnd",
        "blockCount", # synonum for "exonCount"
        "blockSizes", # alternative encoding of exon positions
        "chromStarts",
        "exonFrames"   
    )
)
isBigGenePred <- function(annotation) annotation$type == "bigGenePred"

#----------------------------------------------------------------------
# retrieve annotation-level metadata
# 'genome' and 'annotation' arguments take a single-row data.table describing one entity
# 'chromosome' and 'gene' take a single character entity name
#----------------------------------------------------------------------
# annotation metadata for a genome, chromosome, region or gene
getUcscChromTranscripts <- function(genome, annotation, chromosome = "all", force = FALSE){
    fileName <- paste("transcripts", chromosome, sep = ".")
    file <- getUcscRdsFile(file.path("genomes", genome, "annotations", annotation$track), fileName)
    file <- loadPersistentFile(
        file = file,
        force = force,
        unlink = force,
        ttl = CONSTANTS$ttl$year,
        create = function(file){
            startSpinner(session, message = paste("tabulating transcripts"))
            dt <- getUcscTrackTable(genome, annotation$track, col.names = ucscGenePredTypes[[annotation$type]], force = FALSE)
            saveRDS(if(chromosome == "all") dt else dt[chrom == chromosome], file = file)   
            stopSpinner(session)           
        }
    )
    persistentCache[[file]]$data 
}
getUcscGenomeTranscripts <- function(genome, annotation, force = FALSE){
    getUcscChromTranscripts(genome, annotation, chromosome = "all", force = force)
}
getUcscRegionTranscripts <- function(genome, annotation, coord, force = FALSE){
    transcripts <- getUcscChromTranscripts(genome, annotation, chromosome = coord$chromosome, force = force)
    if(isBigGenePred(annotation)){
        transcripts[chromStart <= coord$end & coord$start <= chromEnd]
    } else {
        transcripts[txStart    <= coord$end & coord$start <= txEnd]
    }
}
getUcscGeneTranscripts <- function(genome, annotation, gene, coord = NULL){
    transcripts <- getUcscGenomeTranscripts(genome, annotation, force = force)
    transcripts[toupper(name2) == toupper(gene)]
}
aggregateUcscTranscriptsToGenes <- function(annotation, transcripts){
    sizeLimit <- 5e6 # size limit correct for naming problems in annotation where name2 applies to widely disparte chromosome regions
    x <- if(isBigGenePred(annotation)) transcripts[, {
        blocks <- unique(data.table(
            sizes  = unlist(strsplit(blockSizes,  ",")),
            starts = unlist(strsplit(chromStarts, ","))            
        ))
        .(
            names       = paste(sort(name), collapse = ","),
            chrom       = chrom[1],
            strand      = strand[1],
            chromStart  = min(chromStart), # the widest annotated span of the gene
            chromEnd    = max(chromEnd),
            blockCount  = max(blockCount),
            blockSizes  = list(as.integer(blocks$sizes)),
            chromStarts = list(as.integer(blocks$starts))
        ) 
    }, by = .(name2)][chromStart - chromEnd < sizeLimit] else transcripts[, {
        exons <- unique(data.table(
            starts = unlist(strsplit(exonStarts, ",")),
            ends   = unlist(strsplit(exonEnds,  ","))            
        )) 
        .(
            names       = paste(sort(name), collapse = ","),
            chrom       = chrom[1],
            strand      = strand[1],
            txStart     = min(txStart), # the widest annotated span of the gene
            txEnd       = max(txEnd),
            exonCount   = max(exonCount),
            exonStarts  = list(as.integer(exons$starts)),
            exonEnds    = list(as.integer(exons$ends))
        )  
    }, by = .(name2)][txEnd - txStart < sizeLimit]
    x[order(name2)]
}

#----------------------------------------------------------------------
# perform post-processing on genes and transcripts
#----------------------------------------------------------------------
setUcscFeatureEndpoints <- function(features, reference){ # convert endpoints to standardized names
    if(isBigGenePred(reference$annotation)){
        features[, ":="(
            start = chromStart,
            end   = chromEnd
        )]
    } else {
        features[, ":="(
            start = txStart,
            end   = txEnd
        )]
    }
    features
}
expandUcscExons <- function(genes, annotation){
    if(nrow(genes) == 0) return(data.table())
    if(isBigGenePred(annotation)){
        genes[, {
            exonStarts <- chromStart + unlist(chromStarts)
            .(
                chrom  = chrom,
                strand = strand, 
                start  = exonStarts, # thus, yields one row per exon
                end    = exonStarts + unlist(blockSizes) - 1
            )
        }, by = .(name2)]
    } else {
        genes[, .(
            chrom  = chrom,
            strand = strand, 
            start  = unlist(exonStarts), # thus, yields one row per exon
            end    = unlist(exonEnds)
        ), by = .(name2)]
    }
}

# list schema from specified track in UCSC database genome -
# api.genome.ucsc.edu/list/schema?genome=hg38;track=wgEncodeGencodeBasicV41
# { "downloadTime": "2022:09:13T15:49:23Z", "downloadTimeStamp": 1663084163, "genome": "hg38", "track": "wgEncodeGencodeBasicV41", "dataTime": "2022-07-12T10:52:01", "dataTimeStamp": 1657648321, "columnTypes": [ { "name": "bin", "sqlType": "smallint(5) unsigned", "jsonType": "number", "description": "Indexing field to speed chromosome range queries"} , { "name": "name", "sqlType": "varchar(255)", "jsonType": "string", "description": "Name of gene (usually transcript_id from GTF)"} , { "name": "chrom", "sqlType": "varchar(255)", "jsonType": "string", "description": "Reference sequence chromosome or scaffold"} , { "name": "strand", "sqlType": "char(1)", "jsonType": "string", "description": "+ or - for strand"} , { "name": "txStart", "sqlType": "int(10) unsigned", "jsonType": "number", "description": "Transcription start position (or end position for minus strand item)"} , { "name": "txEnd", "sqlType": "int(10) unsigned", "jsonType": "number", "description": "Transcription end position (or start position for minus strand item)"} , { "name": "cdsStart", "sqlType": "int(10) unsigned", "jsonType": "number", "description": "Coding region start (or end position for minus strand item)"} , { "name": "cdsEnd", "sqlType": "int(10) unsigned", "jsonType": "number", "description": "Coding region end (or start position for minus strand item)"} , { "name": "exonCount", "sqlType": "int(10) unsigned", "jsonType": "number", "description": "Number of exons"} , { "name": "exonStarts", "sqlType": "longblob", "jsonType": "string", "description": "Exon start positions (or end positions for minus strand item)"} , { "name": "exonEnds", "sqlType": "longblob", "jsonType": "string", "description": "Exon end positions (or start positions for minus strand item)"} , { "name": "score", "sqlType": "int(11)", "jsonType": "number", "description": "score"} , { "name": "name2", "sqlType": "varchar(255)", "jsonType": "string", "description": "Alternate name (e.g. gene_id from GTF)"} , { "name": "cdsStartStat", "sqlType": "enum('none','unk','incmpl','cmpl')", "jsonType": "string", "description": "Status of CDS start annotation (none, unknown, incomplete, or complete)"} , { "name": "cdsEndStat", "sqlType": "enum('none','unk','incmpl','cmpl')", "jsonType": "string", "description": "Status of CDS end annotation (none, unknown, incomplete, or complete)"} , { "name": "exonFrames", "sqlType": "longblob", "jsonType": "string", "description": "Exon frame {0,1,2}, or -1 if no frame for exon"} ] , "shortLabel": "Basic", "type": "genePred", "longLabel": "Basic Gene Annotation Set from GENCODE Version 41 (Ensembl 107)", "itemCount": 108676, "parent": "wgEncodeGencodeV41ViewGenes", "parentParent": "wgEncodeGencodeV41", "shortLabel": "Basic", "trackHandler": "wgEncodeGencode", "priority": "1", "subGroups": "view=aGenes name=Basic", "parent": "wgEncodeGencodeV41ViewGenes on", "longLabel": "Basic Gene Annotation Set from GENCODE Version 41 (Ensembl 107)", "type": "genePred"} 

# Get DNA sequence from specified chromosome in UCSC database genome -
# api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrM
# Get DNA sequence from specified chromosome and start,end coordinates in UCSC database genome -
# api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrM;start=4321;end=5678

# Get track data for specified track in UCSC database genome -
# api.genome.ucsc.edu/getData/track?genome=hg38;track=wgEncodeGencodeBasicV41;maxItemsOutput=100
# { "downloadTime": "2022:09:13T15:52:23Z", "downloadTimeStamp": 1663084343, "genome": "hg38", "dataTime": "2022-07-12T10:52:01", "dataTimeStamp": 1657648321, "trackType": "genePred", "track": "wgEncodeGencodeBasicV41", "wgEncodeGencodeBasicV41": { "chr1": [ { "bin": 0, "name": "ENST00000684719.1", "chrom": "chr1", "strand": "-", "txStart": 67092164, "txEnd": 67134970, "cdsStart": 67093004, "cdsEnd": 67127240, "exonCount": 8, "exonStarts": "67092164,67095234,67096251,67115351,67125751,67127165,67131141,67134929,", "exonEnds": "67093604,67095421,67096321,67115464,67125909,67127257,67131227,67134970,", "score": 0, "name2": "C1orf141", "cdsStartStat": "cmpl", "cdsEndStat": "cmpl", "exonFrames": "0,2,1,2,0,0,-1,-1,"} , { "bin": 0, "name": "ENST00000371007.6", "chrom": "chr1", "strand": "-", "txStart": 67092164, "txEnd": 67231852, "cdsStart": 67093004, "cdsEnd": 67127240, "exonCount": 8, "exonStarts": "67092164,67095234,67096251,67115351,67125751,67127165,67131141,67231845,", "exonEnds": "67093604,67095421,67096321,67115464,67125909,67127257,67131227,67231852,", "score": 0, "name2": "C1orf141", "cdsStartStat": "cmpl", "cdsEndStat": "cmpl", "exonFrames": "0,2,1,2,0,0,-1,-1,"} , 

# Get track data for specified track and chromosome in UCSC database genome -
# api.genome.ucsc.edu/getData/track?genome=hg38;track=wgEncodeGencodeBasicV41;chrom=chr21
# { "downloadTime": "2022:09:13T15:53:27Z", "downloadTimeStamp": 1663084407, "genome": "hg38", "dataTime": "2022-07-12T10:52:01", "dataTimeStamp": 1657648321, "trackType": "genePred", "track": "wgEncodeGencodeBasicV41", "chrom": "chr21", "start": 0, "end": 46709983, "wgEncodeGencodeBasicV41": [ { "bin": 623, "name": "ENST00000624081.1", "chrom": "chr21", "strand": "+", "txStart": 5011798, "txEnd": 5017145, "cdsStart": 5011798, "cdsEnd": 5011798, "exonCount": 4, "exonStarts": "5011798,5012547,5014385,5016934,", "exonEnds": "5011874,5012687,5014471,5017145,", "score": 0, "name2": "ENSG00000279493", "cdsStartStat": "none", "cdsEndStat": "none", "exonFrames": "-1,-1,-1,-1,"} ,

# Wiggle track data for specified track, chromosome with start and end limits in an assembly hub genome -
# api.genome.ucsc.edu/getData/track?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/hub.txt;genome=CAST_EiJ;track=gc5Base;chrom=chr1;start=4321;end=5678
# Wiggle track data for specified track in a UCSC database genome -
# api.genome.ucsc.edu/getData/track?genome=galGal6;track=gc5BaseBw;maxItemsOutput=100
# bigBed data from a UCSC database, chrom and start,end limits -
# api.genome.ucsc.edu/getData/track?genome=galGal6;track=ncbiRefSeqOther;chrom=chr1;start=750000;end=55700000

# list public hubs - api.genome.ucsc.edu/list/publicHubs
# list genomes from specified hub - api.genome.ucsc.edu/list/hubGenomes?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/hub.txt
# list tracks from specified hub and genome - api.genome.ucsc.edu/list/tracks?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/hub.txt;genome=CAST_EiJ
# list chromosomes from specified track in UCSC database genome - api.genome.ucsc.edu/list/chromosomes?genome=hg38;track=gold
# list chromosomes from assembly hub genome -
# api.genome.ucsc.edu/list/chromosomes?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/hub.txt;genome=CAST_EiJ
# list chromosomes from specified track in assembly hub genome -
# api.genome.ucsc.edu/list/chromosomes?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/hub.txt;genome=CAST_EiJ;track=assembly

# Get track data for specified track, chromosome and start,end coordinates in UCSC database genome -
# api.genome.ucsc.edu/getData/track?genome=hg38;track=gold;chrom=chr1;start=47000;end=48000




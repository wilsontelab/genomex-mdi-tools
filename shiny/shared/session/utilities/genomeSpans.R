# # load genome spans for genome_spans trackType
# genomex_loadAllGenomeSpans <- function(bgzFile, ...){
#     getCachedTabix(bgzFile, create = TRUE, force = TRUE)
# }

genomex_loadWindowGenomeSpans <- function(track, reference, coord, bgzFile, ...){

    dmsg("genomex_loadWindowGenomeSpans")
    dmsg(bgzFile)

    x <- getCachedTabix(bgzFile, create = TRUE, force = TRUE)
    #  %>% 
    # getTabixRangeData(coord)  
    
    data.table(chrom = character(), start= integer(), end = integer())

}

# # construct a genome spans name table (BED6, regardless of input)
# genomex_genomeSpansNavTable <- function(track, session, browserId, reference, coord, 
#                                         loadFn){
#     navTableName <- initTrackNav(track, session, "navTable") # table reactive functions are provided below  
#     trackNavDataUnformatted <- reactive({
#         selectedFilePaths <- names(track$settings$items())


#         itemsList <- getItemsData(track, reference, coord, dataFn, parseXY = FALSE)
#         if(!itemsList$hasData) return(trackInfo(track, coord, layout, "no usable data to plot"))


#         spans <- svx_getTrackJunctions(track, selectedSources, loadFn, chromOnly = FALSE)

#         req(nrow(spans) <= 10000)
#         spans
#     })
#     trackNavDataFormatted <- reactive({
#         outCols <- c("chrom","start","end","name","score","strand")
#         spans <- trackNavDataUnformatted()
#         for(cols %in% outCols) if(!(col %in% names(spans))) spans[[col]] <- NA
#         spans[, .SD, .SDcols = outCols]
#     })
#     handleRowClick <- function(selectedRow){
#         req(selectedRow)
#         span <- trackNavDataUnformatted()[selectedRow]
#         handleTrackNavTableClick(track, span$chrom, span$start + 1, span$end)
#     }
#     tagList(
#         trackNavTable(
#             track, 
#             session, 
#             browserId,
#             navTableName, # the name as provided by initTrackNav
#             tableData = trackNavDataFormatted, # populate a table based on track settings, etc.
#             actionFn = handleRowClick
#         )
#     )  
# }

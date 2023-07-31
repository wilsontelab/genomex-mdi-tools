#----------------------------------------------------------------------
# read_pileup trackBrowser track (i.e., a browserTrack)
# expects that input bgz files, and associated tabix bgx.tbi indices, are in format:
#   chromShort(1,2,3,...,'X','Y','M'), startPos(integer), runLength(integer), readBases(##M##A##+##-)...
# where readBases == 11M2A3+4- at a position means:
#   11 reads Matched reference at all bases in the run
#   2 reads have an A base (similar for C,G,T)
#   3 reads were followed by an insertion
#   4 reads were deleted at the base
# the order of operations in readBases is not important, it is simply a base-specific coverage tally
# the logic of the format is to be as compact as possible in tabixed bgz format
# appropriate read pileup files are created the genomex-mdi-tools/pileup pipeline module
#----------------------------------------------------------------------

# constructor for the S3 class; REQUIRED
new_read_pileupTrack <- function(trackId) {
    list(
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        expand = FALSE,
        expand2 = FALSE
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.read_pileupTrack <- function(...) showTrackSamplesDialog(...)

# build method for the S3 class; REQUIRED
build.read_pileupTrack <- function(track, reference, coord, layout){
    build.XY_pileup_track(
        track, reference, coord, layout,
        allValues = c("A","C","G","T","-","+","M"), # places M on top of alt bases
        colorPalette = CONSTANTS$baseColors,
        dataFn = function(){

            # parse the input bgz files
            Pileup_File <- getTrackSetting(track, "Pileup", "Pileup_File", NULL)
            items <- track$settings$items()
            isFile <- !is.null(Pileup_File) && file.exists(Pileup_File)
            isItems <- length(items) > 0
            bgzFiles <- c(
                if(isFile) Pileup_File else character(),
                if(isItems) getBgzFilesFromItems(track, reference, "pileup.txt") else character()
            )
            if(length(bgzFiles) == 0) return(FALSE)

            # adjust chromosome name
            coordShort <- coord
            coordShort$chromosome <- sub("chr", "", coordShort$chromosome)             

            # use tabix to extract data, then process into suitable categories table
            lapply(bgzFiles, function(bgzFile){
                sample <- strsplit(basename(bgzFile), "\\.")[[1]][1]
                pileup <- getCachedTabix(bgzFile) %>% # start is 1-referenced
                          getTabixRangeData(coordShort, col.names = c("chrom","start","length","bases"), 
                                                        colClasses = c("character","integer","integer","character"))
                pileup <- if(is.null(pileup) || nrow(pileup) == 0) data.table(
                    start = integer(),
                    end = integer(),
                    M = integer()
                ) else {
                    pileup[, end := start + length - 1]
                    dcast(
                        pileup[, {
                            counts <- as.numeric(unlist(regmatches(bases, gregexpr('\\d+', bases))))
                            values <-            unlist(regmatches(bases, gregexpr('\\D',  bases)))
                            data.table(base = values, count = counts)
                        }, keyby = .(start, end)],
                        start + end ~ base, 
                        fun.aggregate = sum,
                        value.var = "count",
                        fill = 0
                    )                    
                }
                list( # as required by build.XY_pileup_track
                    pileup = pileup,
                    ylab = sample
                )
            })
        }
    )
}

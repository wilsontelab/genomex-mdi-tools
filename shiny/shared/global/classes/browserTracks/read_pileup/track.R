#----------------------------------------------------------------------
# read_pileup trackBrowser track (i.e., a browserTrack)
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
    coordShort <- coord
    coordShort$chromosome <- sub("chr", "", coordShort$chromosome)    
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
            req(isFile || isItems)
            bgzFiles <- c(
                if(isFile) Pileup_File else character(),
                if(isItems) getBgzFilesFromItems(track, reference, "pileup.txt") else character()
            )
            req(length(bgzFiles) > 0)

            # used tabix to read data, then process into suitable categories table
            lapply(bgzFiles, function(bgzFile){
                pileup <- getCachedTabix(bgzFile) %>% # start is 1-referenced
                getTabixRangeData(coordShort, col.names = c("chrom","start","length","bases"), colClasses = c("character","integer","integer","character"))
                req(pileup, nrow(pileup) > 0)
                pileup[, end := start + length - 1]
                pileup <- dcast(
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
                list( # as required by build.XY_pileup_track
                    pileup = pileup,
                    ylab = strsplit(basename(bgzFile), "\\.")[[1]][1]
                )
            })
        }
    )
}

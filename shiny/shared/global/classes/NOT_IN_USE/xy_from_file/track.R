#----------------------------------------------------------------------
# __TRACK_TYPE__ trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_xy_from_fileTrack <- function() {
    list(      
        click = TRUE,
        hover = FALSE,
        items = FALSE
    )
}

# build method for the S3 class




build.xy_from_fileTrack <- function(settings, input){
    coord <- coordinates(input)
    # d <- settings$Data_Options()
    # req(d$Data_File$value)

    # dmsg(serverEnv$UPLOADS_DIR)
    # dmsg(d$Data_File$value)
    # dataFile <- file.path(serverEnv$UPLOADS_DIR, d$Data_File$value)
    # extension <- tools::file_ext(dataFile)


    # data <- if(extension == "bgz"){

    #     Rsamtools::scanTabix(

    #     )
    #     dmsg("running tabix")
    #     bedr::tabix(
    #         coord$chromosome,
    #         dataFile,
    #         params = NULL,
    #         tmpDir = NULL,
    #         deleteTmpDir = TRUE,
    #         outputDir = NULL,
    #         outputFile = NULL,
    #         check.zero.based = FALSE,
    #         check.chr = TRUE,
    #         check.valid = TRUE,
    #         check.sort = TRUE,
    #         check.merge = TRUE,
    #         verbose = TRUE
    #     )        
    # }

    # dstr(data)

# params
# A string that includes all the extra parameters and arguments to the bedtools commmand. For example if you wanted to do a left outer join you would specificy method as intersect and use params = c("-loj -header"). If you leave input and method as defaults then this is this string represents the full command.




    # inputFile <- if(extension == "bgz"){
    #     R.utils::gunzip(dataFile, ext = "bgz", remove = FALSE)
    # } else dataFile

    # featureFile <- loadPersistentFile(
    #     file = inputFile,
    #     ttl = CONSTANTS$ttl$week,


    #     force = TRUE,

    #     header = FALSE,
    #     col.names = col.names
    # )
    # if(extension == "bgz") unlink(inputFile)

    # dstr(persistentCache[[featureFile]]$data)

    

#     R.utils::gunzip()



# 


#     List of 3
#  $ type  : chr "fileInput"
#  $ accept: chr [1:2] ".bed" ".bgz"
#  $ value : chr "hg38.ideogram.bed.bgz"

    x <- 1:10
    ylim <- range(x)

    list(
        mar = list(top = 0, bottom = 0),
        args = list(
            settings = settings,
            input = input,
            x = x,
            y = x,
            ylim = ylim,
            ylab = ylab(settings, "Test Label")        
        ),
        ylim = ylim, # last values used for click and hover handling
        height = height(settings, 1)
    )



    # coord <- coordinates(input)
    # x <- coord$start:coord$end
    # y <- x * 0.1 + rnorm(coord$width, sd = d$Noise_Factor$value)
    # ylim <- ylim(settings, y)

    # list(
    #     mar = list(top = 0, bottom = 0),
    #     args = list(
    #         settings = settings,
    #         input = input,
    #         x = x,
    #         y = y,
    #         ylim = ylim,
    #         ylab = ylab(settings, "Test Label")        
    #     ),
    #     ylim = ylim, # last values used for click and hover handling
    #     height = height(settings, 1)
    # )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click or track$hover is TRUE, above
click.xy_from_fileTrack <- function(track, x, y){

    dmsg()
    dprint(track$id)
    dprint(paste("x = ", x))
    dprint(paste("y = ", y))

}
hover.xy_from_fileTrack <- function(track, x, y){

    dmsg()
    dprint(track$id)
    dprint(paste("x = ", x))
    dprint(paste("y = ", y))

}

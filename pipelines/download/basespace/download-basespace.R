
# documentation
#   https://developer.basespace.illumina.com/docs/content/documentation/faq/developer-faq
#   https://developer.basespace.illumina.com/docs/content/documentation/rest-api/api-reference
#   https://bioconductor.org/packages/release/bioc/vignettes/BaseSpaceR/inst/doc/BaseSpaceR.pdf
#   https://bioconductor.org/packages/release/bioc/manuals/BaseSpaceR/man/BaseSpaceR.pdf

# load dependencies
library(BaseSpaceR)

# get all environment variables as list
env <- as.list(Sys.getenv())
forceDownload <- as.logical(as.integer(env$FORCE_DOWNLOAD))
listFilesOnly <- as.logical(as.integer(env$LIST_FILES))
maxListLength <- as.integer(env$MAX_LIST_LENGTH)

# log in to BaseSpace with the user's private access token
appAuth <- AppAuth(access_token = env$ACCESS_TOKEN)

# function to retrieve table of Ids and Names
getIdNames <- function(x) {
    x <- as.data.frame(t( sapply(seq_along(x), function(i) c(x[i]$Id, x[i]$Name) ) ))
    names(x) <- c('Id', 'Name')
    x <- x[order(x$Name), ]
    print(x)
}

# check if the user knows the project Id; if not, show them a list
if(is.na(env$PROJECT_ID) || env$PROJECT_ID == "NA"){
    getIdNames(listProjects(appAuth, Limit = maxListLength))
} else {
    
    # check if the user knows the sample name(s); if not, show them a list
    samples <- listSamples(appAuth, projectId = env$PROJECT_ID, Limit = maxListLength)
    if(is.na(env$SAMPLE_NAMES) || env$SAMPLE_NAMES == "NA"){
        message(paste('All Samples for Project', env$PROJECT_ID))
        getIdNames(samples) 
    } else {
        totalFileSize <- 0

        # process one sample name at a time
        sampleNames <- strsplit(env$SAMPLE_NAMES, ',')[[1]]
        for(sampleName in sampleNames){
            outDir <- file.path(env$TASK_DIR, sampleName)
            dir.create(outDir, showWarnings = FALSE)
            
            # collect all Samples for the requested sample name (can be more than one depending on lanes)
            for(i in seq_along(samples)){
                if(samples[i]$Name == sampleName){
                    
                    # act on all Files for the Sample
                    files <- listFiles( samples[i], Limit = 1000 )
                    for(j in seq_along(files)){
                        totalFileSize <- totalFileSize + files[j]$Size
                        localFile <- file.path(outDir, files[j]$Name)
                        
                        # don't re-download existing files unless requested
                        if(listFilesOnly) {
                            print(paste(files[j]$Name, '=', round(files[j]$Size / 1e9, 1), 'GB'))
                        } else if (!file.exists(localFile) || forceDownload){
                            
                            # abort if any file fails
                            tryCatch(
                                getFiles(appAuth, id = files[j]$Id, destDir = outDir, verbose = TRUE), 
                                warning = function(w) stop(paste( "Warning:", conditionMessage(w))), 
                                error   = function(e) stop(paste( "Error:",   conditionMessage(e)))
                            );
                        } else {
                            message(paste(localFile, 'already exists'))
                        }
                    }          
                }
            }   
        }
        
        # if listing, report the total size that will be downloaded, i.e., local space required
        if(listFilesOnly) print(paste('TOTAL', '=', round(totalFileSize / 1e12, 3), 'TB'))
    }
}

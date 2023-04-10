# use a previously constructed fastq index to give random access to a read's SEQ and QUAL

# ba59c9b1-7fbd-43c3-b7b0-bc662bf9770a    PAK57679_pass_909c8bfb_7f82fa6d_0       0
# ccc98ece-eeab-4b5d-bed7-7e86552cc6ec    PAK57679_pass_909c8bfb_7f82fa6d_0       1137

# load the index
loadFastqIndex <- function(env){
    indexFile <- paste0(env$DATA_FILE_PREFIX, ".align.fastq_index.txt.gz")
    index <- fread(indexFile)
    setnames(index, c("qName", "filePrefix", "offset"))
    index[, length := c(diff(offset), 1e7)]
    setkey(index, qName)
    index
}

# retrieve a single read
workingFastqFile <- ""
workingFastqCon <- NULL
getIndexedRead <- function(index, qName_){
    index <- index[qName_]
    fileName <- paste0(index$filePrefix, ".fastq.gz")

    tmpFile <- paste0(env$TMP_DIR_WRK, "/", fileName)
    tmpFile <- "/scratch/wilsonte_root/wilsonte0/wilsonte/7754-SA_P2_SOLO/pilot/HCT_untargeted/DEBUG.fastq";

    if(!(workingFastqFile == tmpFile)){ # can generally expect callers to process reads in order, so keep handles open
        if(!is.null(workingFastqCon)) close(workingFastqCon)
        fastqFile <- paste0(env$INPUT_DIR, "/", fileName)
        # if(!file.exists(tmpFile)) system2("zcat", c(fastqFile, ">", tmpFile))
        workingFastqCon <<- file(tmpFile, "r")
        workingFastqFile <<- tmpFile
    }
    readLines(workingFastqCon, n = 4)
}

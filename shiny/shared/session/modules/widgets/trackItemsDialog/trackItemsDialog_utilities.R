#----------------------------------------------------------------------
# launch a dialog for populating a single browser track's items from a source table
#----------------------------------------------------------------------
showTrackItemsDialog <- function(
    settings, # as passed to the track type's `items` method
    session, 
    title, # title for the dialog
    itemTypePlural = "Items", # for labeling the 'Selected' and Available' boxes
    isFiles = FALSE, # if TRUE, offers a file loader instead of a table of selections
    extensions = character(), # allowed file types when isFiles == TRUE
    tableData, # a reactive, or function with no arguments, that returns a data.frame or data.table
    keyColumn, # the unique identifier column in `tableData()`
    extraColumns = list(), # additional `tableData()` columns to display in the selection panel
    options = list(), # settable values that control the display of each indiviual item
    size = "l" # the requested size of the dialog
){
    id <- "trackItemsDialog"
    nsId <- session$ns(id)
    dialog <- trackItemsDialogServer(
        id,
        isFiles, 
        extensions,
        tableData,
        keyColumn,
        extraColumns,
        options,
        selected = settings$items()
    )
    showUserDialog(
        title = title,
        trackItemsDialogUI(
            nsId,
            itemTypePlural,
            isFiles, 
            keyColumn,
            extraColumns,
            options
        ),
        size = size, 
        type = 'okCancel', 
        easyClose = FALSE,
        fade = FALSE,
        observers = dialog$observers,
        callback = function(...) settings$items( dialog$selected() )
    )
}

# common form of a track items dialog that shows and returns Samples as Sample_ID and Project
showTrackSamplesDialog <- function(track, session, ...){
    showTrackItemsDialog(
        track$settings,
        session,
        title = "Select Samples",
        itemTypePlural = "Samples",
        tableData = reactive({
            uploadName <- appStepNamesByType$upload
            x <- as.data.table(app[[uploadName]]$outcomes$samples())
            req(x)
            x[, .(Sample_ID, Project)]
        }),
        keyColumn = "Sample_ID",
        extraColumns = c("Project"),
        size = "l" # xl
    )
}

# common form of a track items dialog that shows and returns entire data package sources (not their sub-samples)
showTrackSourcesDialog <- function(track, session, ...){
    showTrackItemsDialog(
        track$settings,
        session,
        title = "Select Data Sources",
        itemTypePlural = "Sources",
        tableData = reactive({
            uploadName <- appStepNamesByType$upload
            sourceIds <- names(app[[uploadName]]$outcomes$sources())
            data.table(
                Source_ID = sourceIds,
                Data_Package = sapply(sourceIds, getSourceFilePackageName)
            )
        }),
        keyColumn = "Source_ID",
        extraColumns = c("Data_Package"),
        size = "l"
    )
}

# common form of a track items dialog that allows users to select one or more system files
showTrackFilesDialog <- function(track, session, extensions = character(), ...){
    ext <- if(length(extensions > 0)) paste0("(", paste(extensions, collapse = ", "), ")") else ""
    showTrackItemsDialog(
        track$settings,
        session,
        title = paste("Select Files", ext),
        itemTypePlural = "Files",
        isFiles = TRUE,
        extensions = extensions,
        tableData = NA,
        keyColumn = "fileName",
        extraColumns = "parentDir",
        size = "l"
    )
}

# parse a track's selected samples into named list where sourceId -> data.table of samples
getSourcesFromTrackSamples <- function(selectedSamples){ # selectedSamples is a list of sample lists, each with Sample_ID and Project
    sources <- list()
    uploadName <- appStepNamesByType$upload
    uploadSamples <- as.data.table(app[[uploadName]]$outcomes$samples())
    for(selectedSample in selectedSamples){   
        sourceId <- uploadSamples[Sample_ID == selectedSample$Sample_ID & Project == selectedSample$Project, Source_ID]
        sources[[sourceId]] <- rbind(sources[[sourceId]], as.data.table(selectedSample))
    }    
    sources
}

# assign a unique color to each unique selected sample
getColorsBySelectedSample <- function(selectedTargets, isMultiSample = TRUE, dt = NULL){ # selected sources as returned by getSourcesFromTrackSamples
    allSamples <- if(isMultiSample) unique(unlist(lapply(selectedTargets, function(x) x$Sample_ID)))
                  else unique(dt[, unlist(strsplit(samples, ","))])
    allSamples <- allSamples[allSamples != ""]
    sampleCols <- as.list(1:length(allSamples))
    names(sampleCols) <- paste0(",", allSamples, ",") 
    sampleCols   
}
dt_colorBySelectedSample <- function(dt, selectedTargets, isMultiSample = TRUE){ # dt expected to have samples columns, and usually nSamples, columns
    if(nrow(dt) == 0){
        dt$color <- character()
        return(dt)
    }
    sampleCols <- getColorsBySelectedSample(selectedTargets, isMultiSample, dt)
    if("nSamples" %in% names(dt)){
        dt[, color := sapply(1:.N, function(i){
            if(nSamples[i] > 1) "black" else CONSTANTS$plotlyColors[[sampleCols[[samples[i]]]]]
        })]      
    } else {
        dt[, color := unlist(CONSTANTS$plotlyColors[unlist(sampleCols[samples])])]   
    }
    dt
}

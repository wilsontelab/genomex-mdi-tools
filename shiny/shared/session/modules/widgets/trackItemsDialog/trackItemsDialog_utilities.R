#----------------------------------------------------------------------
# launch a dialog for populating a single browser track's items from a source table
#----------------------------------------------------------------------
showTrackItemsDialog <- function(
    settings, # as passed to the track type's `items` method
    session, 
    title, # title for the dialog
    itemTypePlural = "Items", # for labeling the 'Selected' and Available' boxes
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

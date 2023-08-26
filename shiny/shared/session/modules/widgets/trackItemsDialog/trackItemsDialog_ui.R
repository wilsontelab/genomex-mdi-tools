#----------------------------------------------------------------------
# UI components for the trackItemsDialog widget module 
# for populating a single browser track's items from a source table
#----------------------------------------------------------------------

# module ui function
trackItemsDialogUI <- function(
    id, 
    itemTypePlural, 
    isFiles,
    keyColumn,
    extraColumns,
    options
) {
    ns <- NS(id)
    columnNames <- c("Action", keyColumn, extraColumns, names(options))
    selectedItemsId <- "selectedItems"
    selectedItemsCssId <- ns(selectedItemsId)
    tagList(
        fluidRow(
            box(
                width = 12,
                title = paste("Selected", itemTypePlural),
                tags$table(
                    class = "trackItemsTable",
                    tags$thead(
                        tags$tr(
                            lapply(columnNames, function(col){
                                tags$th( tags$strong(col) ) 
                            })   
                        )     
                    ),
                    tags$tbody(
                        id = selectedItemsCssId,
                        class = selectedItemsId
                    )
                ),
                sortable_js( 
                    css_id = selectedItemsCssId, 
                    options = sortable_options( 
                        onSort = sortable_js_capture_input(input_id = selectedItemsCssId) 
                    )
                )
            )
        ),
        fluidRow(
            if(isFiles) box(
                width = 12,
                serverSourceFilesButtonUI(
                    ns("shinyFilesButton"),
                    multiple = TRUE,
                    buttonType = "primary",
                    style = "width: 250px;"
                )
            ) else bufferedTableUI(
                ns("availableItems"), 
                title = paste("Available", itemTypePlural),
                width = 12
            )
        )
    )
}

#----------------------------------------------------------------------
# UI components for the trackItemsDialog widget module 
# for populating a single browser track's items from a source table
#----------------------------------------------------------------------

# module ui function
trackItemsDialogUI <- function(
    id, 
    itemTypePlural, 
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
            bufferedTableUI(
                ns("availableItems"), 
                title = paste("Available", itemTypePlural),
                width = 12
            )
        )
    )
}

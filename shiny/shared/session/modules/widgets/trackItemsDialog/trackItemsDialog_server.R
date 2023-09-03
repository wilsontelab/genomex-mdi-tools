#----------------------------------------------------------------------
# server components for the trackItemsDialog widget module
# for populating a single browser track's items from a source table
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
trackItemsDialogServer <- function(
    id,
    isFiles,
    tableData,
    keyColumn,
    extraColumns,
    options,
    selected = list()
) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
selectItemLinkId   <- "selectItem"
removeItemLinkId   <- "removeItem"
changeItemOptionId <- "changeItemOption"
observers <- list()
staticColumns <- c(keyColumn, extraColumns)
selected <- reactiveVal(selected) # name = item key, value = option settings
defaults <- lapply(options, function(x){
    if(is.null(x$args$selected)) x$args$value else x$args$selected
})

#----------------------------------------------------------------------
# a single file input button for searching for a system file
#----------------------------------------------------------------------
<<<<<<< HEAD
# serverSourceFilesButtonServer(
#     "shinyFilesButton", 
#     input, 
#     session, 
#     rw = "read", 
#     filetypes = "bgz",
#     loadFn = function(files) {
#         req(files, nrow(files) > 0)
#         files <- data.table(files)
#         for(i in nrow(files)){
#             path <- files[i, datapath]
#             selected <- selected()
#             if(!is.null(selected[[path]])) return()
#             selected[[path]] <- c(files[i, name], basename(dirname(path)))
#             names(selected[[path]]) <- staticColumns
#             selected(selected)    
#             insertSelectedRow(path)
#             updateRowActions()
#         }
#     }
# )
=======
serverFileSelector <- serverFileSelectorServer(
    "serverFileSelector",
    extensions = "mdi"
)
observers$serverFileSelector <- observeEvent(serverFileSelector$selectedFile(), {
    path <- serverFileSelector$selectedFile()
    req(path)
    selected <- selected()
    if(!is.null(selected[[path]])) return()
    selected[[path]] <- c(basename(path), basename(dirname(path)))
    names(selected[[path]]) <- staticColumns
    selected(selected)    
    insertSelectedRow(path)
    updateRowActions()
})
>>>>>>> main

#----------------------------------------------------------------------
# the table of available source items
#----------------------------------------------------------------------
selectHandlerId <- session$ns(selectItemLinkId)
sourceData <- reactive({ # add row selection links to the caller's table
    if(isFiles) return(NULL)
    tableData <- tableData()
    req(tableData)
    cbind(
        Action = tableActionLinks(selectHandlerId, nrow(tableData), "Select", useMdiSharedHandler = TRUE),
        as.data.table(tableData)
    )
})
sourceTable <- if(isFiles) list(observers = list()) else bufferedTableServer(
    "availableItems",
    id,
    input,
    tableData = sourceData,
    selection = 'none',
    options = list(
        searchDelay = 0,
        lengthMenu = c(5, 10, 20),
        pageLength = 5 # TODO: is 10 better? 
    )
)

#----------------------------------------------------------------------
# initialize cached item selections
#----------------------------------------------------------------------
insertSelectedRow <- function(key){
    row <- selected()[[key]]
    insertUI(
        # paste0("#", session$ns("selectedItems")), # a version of Shiny stopped properly handling dom id's for some elements in modals 
        ".trackItemsTable tbody", 
        where = "beforeEnd",
        ui = tags$tr(
            "data-key" = key,
            tags$td(actionLink(session$ns(removeItemLinkId), "Remove")),
            lapply(staticColumns, function(col) tags$td(row[[col]])),
            lapply(names(options), function(option){
                args <- c(
                    list(
                        inputId = session$ns(option),
                        label = NULL
                    ),
                    options[[option]]$args
                )
                if(is.null(args$selected)) args$value <- row[[option]]
                                   else args$selected <- row[[option]]
                tags$td( 
                    class = changeItemOptionId,
                    do.call(options[[option]]$type, args)
                )
            }) 
        ),
        immediate = TRUE
    )    
}
updateRowActions <- function(){
    session$sendCustomMessage("updateRemoveItemLinks", session$ns(removeItemLinkId))
    session$sendCustomMessage("updateItemOptionChanges", list(
        id = session$ns(changeItemOptionId),
        options = options
    ))    
}
initSelected <- observe({
    for(key in names(selected())) insertSelectedRow(key)
    updateRowActions()    
    initSelected$destroy()
})

#----------------------------------------------------------------------
# update item selections
#----------------------------------------------------------------------

# add a selected item
# observers$selectItemLinkId <- observeEvent(input[[selectItemLinkId]], {
#     i <- getTableActionLinkRow(input, selectItemLinkId)
#     row <- tableData()[i]
#     key <- row[[keyColumn]]
#     selected <- selected()
#     if(!is.null(selected[[key]])) return()
#     selected[[key]] <- c(
#         lapply(staticColumns, function(col) row[[col]]),
#         defaults
#     )
#     names(selected[[key]]) <- c(staticColumns, names(defaults))
#     selected(selected)    
#     insertSelectedRow(key)
#     updateRowActions()
# }, ignoreInit = TRUE)
selectItemHandler <- function(x){
    i <- getTableActionLinkRow2(x)
    row <- tableData()[i]
    key <- row[[keyColumn]]
    selected <- selected()
    if(!is.null(selected[[key]])) return()
    selected[[key]] <- c(
        lapply(staticColumns, function(col) row[[col]]),
        defaults
    )
    names(selected[[key]]) <- c(staticColumns, names(defaults))
    selected(selected)    
    insertSelectedRow(key)
    updateRowActions()
}
addMdiSharedEventHandler(selectHandlerId, selectItemHandler)

# remove a selected item
observers$removeItemLinkId <- observeEvent(input[[removeItemLinkId]], {
    selected <- selected()
    selected[[input[[removeItemLinkId]]]] <- NULL
    selected(selected)
}, ignoreInit = TRUE)

# reorder the selected items
observers$selectedItems <- observeEvent(input$selectedItems, {
    keys <- sapply(input$selectedItems, function(x) strsplit(x, "\t")[[1]][2])
    selected <- selected()
    selected(selected[keys])
}, ignoreInit = TRUE)

#----------------------------------------------------------------------
# handle item option changes
#----------------------------------------------------------------------
observers$changeItemOptionId <- observeEvent(input[[changeItemOptionId]], {
    selected <- selected()
    d <- input[[changeItemOptionId]]
    selected[[d$key]][[d$option]] <- d$value
    selected(selected)
}, ignoreInit = TRUE)

#----------------------------------------------------------------------
# return value
#----------------------------------------------------------------------
list(
    selected = selected,
    observers = c(
        sourceTable$observers,
        observers
    )
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

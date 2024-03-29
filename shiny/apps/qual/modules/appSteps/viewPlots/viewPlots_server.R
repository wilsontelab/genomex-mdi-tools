#----------------------------------------------------------------------
# server components for the viewPlots appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
viewPlotsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'viewPlots'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    # settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# add server code sections as needed
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
) 
mdiInteractivePlotServer(
    "qualPlots",
    hover = FALSE,    
    click = FALSE,
    brush = FALSE,
    delay = 500,  
    contents = reactive({
        sourceId <- sourceId()      
        req(sourceId)
        list(
            pngFile = getSourceFilePath(sourceId, "qualPlots"),
            layout = list(
                width = 900,
                height = 900,
                pointsize = 8, # defaults to 8
                dpi = 96
            )
        )   
    })
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

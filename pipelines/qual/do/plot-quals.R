
# make distribution plots to understand the range of read QUAL scores

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(data.table)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'PLOTS_DIR',
        'PLOT_PREFIX',
        'INTERIM_FILE'        
    ),
    integer = c(
        'N_READS',
        'N_TERMINAL_BASES'
    )
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) ########################
#=====================================================================================

#=====================================================================================
# load the data
message("loading and parsing molecules")
d <- fread(env$INTERIM_FILE)
dd <- d[sample(.N, min(.N, 100000))]

# report on read lengths
message()
message("Read Lengths")
print(data.frame(
    q5  = quantile(d$length, 0.05),
    q50 = quantile(d$length, 0.5),
    q95 = quantile(d$length, 0.95)
))
message()

# setup composite plots
dir.create(env$PLOTS_DIR, showWarnings = FALSE)
png(
    filename = paste0(env$PLOT_PREFIX, ".qual-plots.png"),
    width = 6, 
    height = 6, 
    units = "in", 
    pointsize = 8,
    res = 600, 
    type = "cairo"
)
layout(matrix(1:4, ncol = 2, byrow = TRUE))
lim <- range(0, 40)
lab <- list(
    read  = "Whole Read",
    first = paste("First", env$N_TERMINAL_BASES, "Bases"),
    last  = paste("Last",  env$N_TERMINAL_BASES, "Bases")
)
col <- list(
    read  = "blue",
    first = "green3",
    last  = "red3" 
)

# histograms of read, first, last
message("plotting distributions")
plot(NA, NA, typ="n", xlim = lim, ylim = c(0, 0.5), xlab = "Average Base Phred Score", ylab = "Density")
for(span in names(lab)) lines(density(d[[span]]), col = col[[span]])
legend("topleft", legend = unlist(lab), col = unlist(col), lwd = 2) 

# density plot of first vs. last QUAL
message("plotting last vs. first")
plot(
    jitter(d$first, amount = 0.5), jitter(d$last, amount = 0.5), pch = ".", col = rgb(0,0,0,0.05),
    xlim = lim, ylim = lim, 
    xlab = lab$first, ylab = lab$last
)

# histograms of last / first
message("plotting ratio")
plot(density(d$last / d$first), xlim = c(0, 1), xlab = paste(lab$last, "/", lab$first), ylab = "Density", main = "")
#=====================================================================================

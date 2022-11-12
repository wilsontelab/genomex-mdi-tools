parseUnit <- function(unit, value){
    if(unit == "auto"){
        if(value >= 1e9){
            unit <- "Gb"
        } else if(value >= 1e6){
            unit <- "Mb"
        } else if (value >= 1e3){
            unit <- "kb"
        }
    } 
    if(unit == "Gb") {
        multiplier <- 1e9
    } else if(unit == "Mb") {
        multiplier <- 1e6
    } else if(unit == "kb"){
        multiplier <- 1e3
    } else {
        multiplier <- 1
        unit <- "bp"
    }
    list(unit = unit, multiplier = multiplier)
}

#--------------------------------------------------------------
# define global constants
#--------------------------------------------------------------
# you must EXTEND not REPLACE the framework's global CONSTANTS object!
# e.g., CONSTANTS$name <- value
#--------------------------------------------------------------
CONSTANTS$baseColors <- list( # generally follow IGV base color conventions
    M = rgb(0.75, 0.75, 0.75), # any base match = light grey
    A = rgb(0,    0.8,    0), # green  
    C = rgb(0,    0,    1), # blue
    G = rgb(0.82, 0.43, 0),
    T = rgb(1,    0,    0), # red
    N = rgb(0.9, 0.9, 0.9), # N, treated as M
    X = NA,                 # X, treated as missing, no color
    '-' = rgb(0.1, 0.1, 0.1),     # deleted/missing = black
    '+' = rgb(0.75,   0,    0.75), # insertion = purple
    D = rgb(0.1, 0.1, 0.1),     # deleted/missing = black
    I = rgb(0.75,   0,    0.75), # insertion = purple
    ' ' = rgb(1,   1,   1)  # white-space
)

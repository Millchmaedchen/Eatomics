## The code is developed by Milena Kraus and Mariet Stephen. 

# This is a Shiny web application (ui + server) for quantitative proteomics data analysis.
# You can run the application by clicking
# the 'Run App' button above.

# Install or load neccessary R packages 
install.packages("pacman") # for easy install and update #+#
library(pacman) #+#

p_load("rstudioapi")

# Shiny libraries
p_load("shiny") #+#
p_load("shinythemes")
p_load('shinycssloaders')
p_load("shinyalert")
p_load('shinyFiles')
p_load("shinyWidgets")
p_load("shinycssloaders")
p_load("shinyBS")

# Load data
#p_load("readxl")
#p_load("readr")
#p_load(openxlsx)

# Visualization
p_load("pheatmap")
p_load("RColorBrewer")
p_load("plotly")
p_load("ggplot2")
p_load("gridExtra")
p_load("ggthemes")
#install_github("vqv/ggbiplot")
#library(ggbiplot)
#p_load("ggiraph")
p_load("autoplotly")
#p_load("EnhancedVolcano")
#p_load(kableExtra)
p_load("ggrepel")
p_load("gtools")


# Tidy data
p_load("tidyverse")
#p_load("dplyr")
#p_load("data.table")
#p_load(DT)
p_load("janitor")
#p_load(broom)

# Analysis logic
p_load("sva")
p_load("imputeLCMD")
p_load("modelr")
p_load("limma")

# Report
p_load(markdown)

## no clue
#p_load("stringi")
#p_load("rgl")
#p_load("crosstalk")
#p_load("markdown")
#p_load(xtable)
#p_load("rlang")
#p_load(nlme)


# Needed for server setup
#p_load("RJDBC")
#p_load("devtools")

# Load non-reactive helper functions 
source('R/helpers.R')

# Load dependency on ssGSEA algorithm 
source('R/ssGSEA_PSEA.R')


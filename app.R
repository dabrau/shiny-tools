library(shiny)
setwd("~/Projects/shiny_tools")
source("./ui.R")
source("./server.R")

# See above for the definitions of ui and server
shinyApp(ui = ui, server = server)
library(shiny)
setwd("~/Projects/shiny_tools")
source("./server.R")
source("./ui.R")

# See above for the definitions of ui and server
shinyApp(ui = ui, server = server)
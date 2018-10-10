library(shiny)
setwd("~/Projects/shiny_tools")
source("./server.R")
source("./ui.R")

# options(repos = BiocInstaller::biocinstallRepos())
shinyApp(ui = ui, server = server)
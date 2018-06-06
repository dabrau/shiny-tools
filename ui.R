library(plotly)

# Define UI for app that draws a histogram ----
ui <- fluidPage(# App title ----
                titlePanel("HEATMAP"),
                
                # Sidebar layout with input and output definitions ----
                sidebarLayout(
                  # Sidebar panel for inputs ----
                  sidebarPanel(
                    a(href="https://www.google.com", "Example Matrix"),
                    fileInput(
                      "matrix",
                      "Upload Matrix as Tab delimited .txt File",
                      multiple = FALSE,
                      accept = c("text/tsv",
                                 "text/tab-separated-values,text/plain",
                                 ".txt")
                    ),
                    
                    br(),
                    
                    tabsetPanel(
                      tabPanel(
                        "Options",
                        h4("Scaling / Normalization"),
                        checkboxInput("tmm", "TMM normalize", value = TRUE),
                        checkboxInput("log", "log values", value = FALSE),
                        checkboxInput("z_score", "Z - score rows", value = TRUE),
                        
                        br(),
                        
                        h4("Gene Name Mapping"),
                        a(href="https://www.google.com", "Ensembl to HUGO"),
                        checkboxInput("hugo_names", "HUGO names", value = TRUE),
                        
                        br(),
                        
                        h4("Clustering"),
                        selectInput(
                          "dendrogram",
                          h5("Dendrogram"),
                          choices = list("none", "row", "column", "both"),
                          selected = "both"
                        ),
                        
                        selectInput(
                          "dist_method",
                          h5("Distance Method"),
                          choices = list(
                            "euclidean",
                            "maximum",
                            "manhattan",
                            "canberra",
                            "binary",
                            "minkowski"
                          ),
                          selected = "euclidean"
                        ),
                        
                        selectInput(
                          "hclust_method",
                          h5("Clustering Method"),
                          choices = list(
                            "ward.D",
                            "ward.D2",
                            "single",
                            "complete",
                            "average",
                            "mcquitty",
                            "median",
                            "centroid"
                          ),
                          selected = "ward.D2"
                        )
                      ),
                      tabPanel(
                        "Columns",
                        a(href="https://www.google.com", "Example Column List"),
                        fileInput(
                          "column_list",
                          "Upload Column list as .txt File",
                          multiple = FALSE,
                          accept = c("text/tsv",
                                     "text/tab-separated-values,text/plain",
                                     ".txt")
                        ),
                        
                        br(),
                        
                        a(href="https://www.google.com", "Example Column Color Labels"),
                        fileInput(
                          "column_labels",
                          "Upload Column Color labels as Tab delimited .txt File",
                          multiple = FALSE,
                          accept = c("text/tsv",
                                     "text/tab-separated-values,text/plain",
                                     ".txt")
                        )
                      ),
                      tabPanel(
                        "Rows",
                        a(href="https://www.google.com", "Example Row List"),
                        fileInput(
                          "row_list",
                          "Upload Row list as .txt File",
                          multiple = FALSE,
                          accept = c("text/tsv",
                                     "text/tab-separated-values,text/plain",
                                     ".txt")
                        )
                      )
                    )
                  ),
                  
                  # Main panel for displaying outputs ----
                  mainPanel(# Output: Histogram ----
                            plotlyOutput(outputId = "heatmap"))
                ))
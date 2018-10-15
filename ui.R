library(plotly)
library(colourpicker)

ui <- fluidPage(
                tags$head(
                  tags$style(HTML("

                  .matrix-col {
                    height: 40px;
                    display: flex;
                    margin-top: 5px;
                    border-bottom: 1px ridge #eee;
                  }

                  .matrix-col > div:first-child {
                    flex-grow: 1;
                  }
                  
                  .matrix-col > div:nth-child(2) {
                    max-width: 90px;
                  }

                  .matrix-col > .checkbox {
                    margin-top: 26px;
                  }

                  .matrix-col > div:nth-child(2) > .colourpicker {
                    position: relative;
                    bottom: 20px;
                  }
              
                  "))
                ),
  
                # App title ----
                titlePanel("HEATMAP"),
                
                # Sidebar layout with input and output definitions ----
                sidebarLayout(
                  # Sidebar panel for inputs ----
                  sidebarPanel(
                    downloadLink("example_matrix", "download matrix example"),
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
                        checkboxInput("tmm", "TMM normalize", value = FALSE),
                        checkboxInput("log", "log values", value = FALSE),
                        checkboxInput("z_score", "Z - score rows", value = FALSE),
                        
                        br(),
                        
                        h4("Gene Name Mapping"),
                        downloadLink("hugo_mapping", "download Ensembl to HUGO reference"),
                        checkboxInput("hugo_names", "HUGO names", value = FALSE),
                        
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
                        # downloadLink("all_columns_list", "download list of all columns"),
                        # fileInput(
                        #   "column_list",
                        #   "Upload Column list as .txt File",
                        #   multiple = FALSE,
                        #   accept = c("text/tsv",
                        #              "text/tab-separated-values,text/plain",
                        #              ".txt")
                        # ),
                        uiOutput("col_checkbox")
                      ),
                      tabPanel(
                        "Rows",
                        # downloadLink("all_rows_list", "download list of all rows"),
                        # fileInput(
                        #   "row_list",
                        #   "Upload Row list as .txt File",
                        #   multiple = FALSE,
                        #   accept = c("text/tsv",
                        #              "text/tab-separated-values,text/plain",
                        #              ".txt")
                        # ),
                        uiOutput("row_checkbox")
                      )
                    )
                  ),
                  
                  # Main panel for displaying outputs ----
                  mainPanel(# Output: Heatmap
                    downloadButton("download_hm_matrix", label = "Download Heatmap Matrix"),
                    plotlyOutput(outputId = "heatmap"))
                ))
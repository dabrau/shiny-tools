library(plotly)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  
  # App title ----
  titlePanel("HEATMAP"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    
    
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      fileInput("matrix", "Upload Tab delimited .txt File",
                multiple = FALSE,
                accept = c("text/tsv",
                           "text/tab-separated-values,text/plain",
                           ".txt")),
      
      br(),
      
      tabsetPanel(
        tabPanel("Options",
          h4("Scaling / Normalization"),
          checkboxInput("tmm", "TMM normalize", value = TRUE),
          checkboxInput("log", "log counts", value = TRUE),
          checkboxInput("z_score", "Z - score rows", value = TRUE),
      
          br(),
      
          h4("Gene Name Mapping"),
          checkboxInput("hugo_names", "HUGO names", value = TRUE),
      
          br(),
      
          h4("Clustering"),
          selectInput("dendrogram", h5("Dendrogram"), 
                  choices = list("none", "row", "column", "both"), selected = "both"),
      
          selectInput("dist_method", h5("Distance Method"), 
                  choices = list("euclidean", "maximum", "manhattan", "canberra",
                                 "binary", "minkowski"), selected = "euclidean"),
      
          selectInput("hclust_method", h5("Clustering Method"), 
                  choices = list("ward.D", "ward.D2", "single", "complete",
                                 "average", "mcquitty", "median", "centroid"), selected = "ward.D2")
      ),
      tabPanel("Columns", h1("summary")),
      tabPanel("Rows", h1("table"))
      )
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotlyOutput(outputId = "heatmap")
      
    )
  )
)
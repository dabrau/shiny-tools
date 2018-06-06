library(tidyverse)
library(edgeR)
library(RColorBrewer)
library(plotly)
library(shinyHeatmaply)

scale_gene_matrix <- function(matrix) {
  t(
    scale(
      t(log(matrix[rowSums(matrix) > 0, ] + .0001))
    )
  )
}

tmm_norm <- function(counts_matrix) {
  dge <- DGEList(counts = counts_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  logCPM <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0.25)
  logCPM
}

tsv_to_matrix <- function(path) {
  df <- read_tsv(path)
  mat <- as.matrix(df[, 2:ncol(df)])
  rownames(mat) <- df[[1]]
  mat
}

#stub_data
# counts_df <- read_tsv("./stub_data.txt")
# counts_matrix <- as.matrix(counts_df[1:20, 2:ncol(counts_df)])
# rownames(counts_matrix) <- counts_df$gene_id[1:20]
# heatmap_data <- counts_matrix %>% tmm_norm %>% scale_gene_matrix

server <- function(input, output) {
  heatmap_data <- reactive({
    if (!is.null(input$matrix)) {
      hm_matrix <- tsv_to_matrix(input$matrix$datapath)
      hm_matrix <- hm_matrix[rowSums(hm_matrix) > 0, ]
      
      
      if (ncol(hm_matrix) > 20) {
        hm_matrix <- hm_matrix[, 1:20]
      }
      
      if (input$tmm == TRUE) {
        hm_matrix <- tmm_norm(hm_matrix)
      }
      
      if(nrow(hm_matrix) > 300) {
       hm_matrix <- hm_matrix[1:300, ]
      }

      if (input$log == TRUE) {
        hm_matrix <- log(hm_matrix + .0001)
      }
      
      if (input$z_score == TRUE) {
        hm_matrix <- t(scale(t(hm_matrix)))
      }
      
      hm_matrix
    } else matrix(c(c(1,1), c(1,1)), nrow = 2, ncol = 2)
  })
  
  output$heatmap <- renderPlotly({
    plot <- heatmaply(
      heatmap_data(),
      colors = colorRampPalette(c("green", "black", "red"))(80),
      margins = c(200, 200),
      dendrogram = input$dendrogram,
      dist_method = input$dist_method,
      hclust_method = input$hclust_method,
      grid_gap = 1,
      branches_lwd = 0.3
    ) %>% layout (height = 800)
  })
}
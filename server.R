library(tidyverse)
library(edgeR)
library(RColorBrewer)
library(plotly)
library(shinyHeatmaply)

scale_gene_matrix <- function(matrix) {
  t(scale(t(log(matrix[rowSums(matrix) > 0,] + .0001))))
}

tmm_norm <- function(counts_matrix) {
  dge <- DGEList(counts = counts_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  logCPM <-
    cpm(
      dge,
      normalized.lib.sizes = TRUE,
      log = TRUE,
      prior.count = 0.25
    )
  logCPM
}

tsv_to_matrix <- function(path) {
  df <- read_tsv(path)
  mat <- as.matrix(df[, 2:ncol(df)])
  rownames(mat) <- df[[1]]
  mat
}

hugo_mapping <- read_tsv("./ensembl_hugo_mapping.tsv")

ensembl_to_hugo_rownames <- function(hm_matrix) {
  ensembl_ids <- rownames(hm_matrix)
  hugo <- hugo_mapping %>% filter(ensembl %in% ensembl_ids)
  renamed <- hm_matrix[hugo$ensembl,]
  rownames(renamed) <- hugo$hugo
  renamed
}


server <- function(input, output) {
  heatmap_data <- reactive({
    if (!is.null(input$matrix)) {
      # read in uploaded matrix
      hm_matrix <- tsv_to_matrix(input$matrix$datapath)
      
      if (is.null(input$column_list)) {
        # limit number of columns if column list not provided
        if (ncol(hm_matrix) > 20) {
          hm_matrix <- hm_matrix[, 1:20]
        }
      } else {
        cols <- read_tsv(input$column_list$datapath, col_names = FALSE)[[1]]
        cols <- cols[cols %in% colnames(hm_matrix)]
        hm_matrix <- hm_matrix[, cols]
      }
      
      # remove rows with all 0, assumption all values are positive maybe should check explicitly
      hm_matrix <- hm_matrix[rowSums(hm_matrix) > 0,]
      
      if (input$tmm == TRUE) {
        hm_matrix <- tmm_norm(hm_matrix)
      }
      
      if (input$hugo_names) {
        hm_matrix <- ensembl_to_hugo_rownames(hm_matrix)
      }
      
      if (is.null(input$row_list)) {
        # limit number of rows if row list not provided
        if (nrow(hm_matrix) > 300) {
          hm_matrix <- hm_matrix[1:300,]
        }
      } else {
        rows <- read_tsv(input$row_list$datapath, col_names = FALSE)[[1]]
        rows <- rows[rows %in% rownames(hm_matrix)]
        hm_matrix <- hm_matrix[rows,]
      }
      
      if (input$log == TRUE) {
        hm_matrix <- log(hm_matrix + .0001)
      }
      
      if (input$z_score == TRUE) {
        hm_matrix <- t(scale(t(hm_matrix)))
      }
      
      hm_matrix
    } else {
      # placeholder
      matrix(c(c(1, 1), c(1, 1)), nrow = 2, ncol = 2)
    }
  })
  
  output$heatmap <- renderPlotly({
    # prevent reordering when dendrogram is 'none', 'row', or 'column'
    dendrogram_params <- list("Rowv" = TRUE, "Colv" = TRUE)
    if (input$dendrogram == "none") {
      dendrogram_params$Rowv = FALSE
      dendrogram_params$Colv = FALSE
    } else if (input$dendrogram == 'row') {
      dendrogram_params$Colv = FALSE
    } else if (input$dendrogram == 'column') {
      dendrogram_params$Rowv = FALSE
    }
    
    hm <- heatmaply(
      heatmap_data(),
      colors = colorRampPalette(c("green", "black", "red"))(80),
      margins = c(200, 200),
      Rowv = dendrogram_params$Rowv,
      Colv = dendrogram_params$Colv,
      dendrogram = input$dendrogram,
      dist_method = input$dist_method,
      hclust_method = input$hclust_method,
      grid_gap = 1,
      branches_lwd = 0.3
    ) %>% layout(height = 800)
  })
}
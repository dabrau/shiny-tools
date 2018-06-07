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

hugo_mapping <- read_tsv("./ensembl-hugo-mapping.txt")

ensembl_to_hugo_rownames <- function(hm_matrix) {
  ensembl_ids <- rownames(hm_matrix)
  hugo <- hugo_mapping %>% filter(ensembl %in% ensembl_ids)
  renamed <- hm_matrix[hugo$ensembl,]
  rownames(renamed) <- hugo$hugo
  renamed
}

# example input files for download
matrix_example <- read_tsv("./example-inputs/example-counts-matrix.txt")
column_list_example <- read_tsv("./example-inputs/example-column-list.txt", col_names = FALSE)
column_labels_example <- read_tsv("./example-inputs/example-column-labels.txt")
row_list_example <- read_tsv("./example-inputs/example-row-list.txt", col_names = FALSE)

tsv_dl_handler <- function(df, filename, col_names) {
  downloadHandler(
    filename = filename,
    content = function(con) {
      write_tsv(df, con, col_names = col_names)
    },
    contentType = "text/tab-separated-values"
  )
}

server <- function(input, output) {
  output$example_matrix <- tsv_dl_handler(matrix_example, "example-counts-matrix.txt", TRUE)
  output$example_column_list <- tsv_dl_handler(column_list_example, "example-column-list.txt", FALSE)
  output$example_column_labels <- tsv_dl_handler(column_labels_example, "example-column-labels.txt", TRUE)
  output$example_row_list <- tsv_dl_handler(row_list_example, "example-row-list.txt", FALSE)
  
  output$hugo_mapping <- tsv_dl_handler(hugo_mapping, "ensembl-hugo-mapping.txt", TRUE)
  
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
      
      # remove rows with all 0
      rows <- apply(hm_matrix, 1, function(row) sum(row != 0) != 0)
      hm_matrix <- hm_matrix[rows, ]
      
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
      placeholder <- matrix(c(c(1, 1), c(1, 1)), nrow = 2, ncol = 2)
      colnames(placeholder) <- c("item1", "item2")
      placeholder
    }
  })
  
  column_color_labels <- reactive({
    if (!is.null(input$column_labels)) {
      df <- read_tsv(input$column_labels$datapath)
      
      colors_map <- list()
      for (label_group in 2:ncol(df)) {
        color <- df[[label_group]]
        names(color) <- df[[1]]
        colors_map[[colnames(df)[[label_group]]]] <- color
      }
      
      colors_map
    } else list()
    
  })
  
  output$heatmap <- renderPlotly({
    hm_params <- list(
      x = heatmap_data(),
      colors = colorRampPalette(c("green", "black", "red"))(80),
      margins = c(200, 200),
      Rowv = TRUE,
      Colv = TRUE,
      dendrogram = input$dendrogram,
      dist_method = input$dist_method,
      hclust_method = input$hclust_method,
      grid_gap = 1,
      branches_lwd = 0.3
    )
    
    # prevent reordering when dendrogram is 'none', 'row', or 'column'
    if (input$dendrogram == "none") {
      hm_params$Rowv = FALSE
      hm_params$Colv = FALSE
    } else if (input$dendrogram == 'row') {
      hm_params$Colv = FALSE
    } else if (input$dendrogram == 'column') {
      hm_params$Rowv = FALSE
    }
    
    if (!is.null(input$column_labels)) {
      # generate dataframe for column labels
      col_labels <- map(column_color_labels(), function(label_map) {
        map_chr(colnames(heatmap_data()), function(col) {
          if (col %in% names(label_map)) {
            if (is.na(label_map[[col]])) {
              "unknown"
            } else label_map[[col]] 
          } else "unknown"
        }) %>% as.factor
      }) %>% as.data.frame
      
      hm_params$col_side_colors = col_labels
    }
    
    do.call(heatmaply, hm_params) %>% layout(height = 800)
  })
}
library(tidyverse)
library(edgeR)
library(RColorBrewer)
library(plotly)
library(shinyHeatmaply)
library(scales)

options(shiny.maxRequestSize=30*1024^2) 

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
  
  #remove empty columns
  cols_to_keep <- map_lgl(colnames(df), function(col_name) {
    sum(is.na(df[[col_name]])) != nrow(df)
  })
  
  #remove empty rows
  rows_to_keep <- map_lgl(rownames(df), function(row_name) {
    sum(is.na(df[row_name, ])) != ncol(df)
  })
  
  df <- df[rows_to_keep, cols_to_keep]
  
  mat <- as.matrix(df[, 2:ncol(df)])
  rownames(mat) <- df[[1]]
  mat
}

hugo_mapping <- read_tsv("./ensembl-hugo-mapping.txt")

ensembl_to_hugo_rownames <- function(hm_matrix) {
  ensembl_ids <- rownames(hm_matrix)
  labels <- map_chr(ensembl_ids, function(id) {
    name <- hugo_mapping %>% filter(ensembl == id) %>% dplyr::select(hugo) %>% unlist
    if (length(name) == 0) { id } else { name[1] }
  })
  
  rownames(hm_matrix) <- labels
  hm_matrix
}

# example input files for download
matrix_example <- read_tsv("./example-inputs/example-counts-matrix.txt")

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
  hm_matrix <- reactive({
    
    if (!is.null(input$matrix)) {
      mat <- tsv_to_matrix(input$matrix$datapath)
      # remove rows with all 0
      rows <- apply(mat, 1, function(row) sum(row != 0) != 0)
      mat <- mat[rows, ]
    } else {
      # placeholder
      mat <- matrix(c(c(1, 1), c(1, 1)), nrow = 2, ncol = 2)
      colnames(mat) <- c("item1", "item2")
    }
    
    mat
  })
  
  
  default_color <- "#63B8FF"
  output$col_checkbox <- renderUI({
    choice <-  unique(colnames(hm_matrix()))
    cols <- lapply(choice, FUN = function(col_name) {
      tags$div(class = "matrix-col",
               checkboxInput(paste("column", col_name, sep = "_"), col_name, TRUE),
               colourInput(paste("color", col_name, sep = "_"), "", value = default_color,
                           allowedCols = c(default_color, "#1E90FF", "#FFFFFF"),
                           palette = "limited")
      )
    })
    
    do.call(tags$div, cols)
  })
  
  output$row_checkbox <- renderUI({
    mat <- hm_matrix()
    if (input$hugo_names) {
      mat <- ensembl_to_hugo_rownames(hm_matrix())
    }
    choice <-  unique(rownames(mat))
    checkboxGroupInput("row_checkbox","", choices = choice, selected = choice)
  })
  
  columns <- reactive({
    col_inputs <- names(input)[str_detect(names(input), "column_")]
    col_names <- str_sub(col_inputs, start = 8)
    values <- map_lgl(col_inputs, function(inp) input[[inp]])
    col_names[values]
  })

  heatmap_data <- reactive({
    hm_mat <- hm_matrix()
    #row and column selectors
    if (!is.null(input$row_checkbox)) {
      hm_mat <- hm_mat[rownames(hm_mat) %in% input$row_checkbox, ]
    }
    
    if (length(columns()) != 0) {
      hm_mat <- hm_mat[, colnames(hm_mat) %in% columns()]
    }
    
    if (input$tmm == TRUE) {
      hm_matrix <- tmm_norm(hm_mat)
    }
    
    if (input$hugo_names) {
      hm_mat <- ensembl_to_hugo_rownames(hm_mat)
    }
    
    if (input$log == TRUE) {
      hm_mat <- log(hm_mat + .0001)
    }
    
    if (input$z_score == TRUE) {
      hm_mat <- t(scale(t(hm_mat)))
    }
    
    hm_mat
  })
  
  col_colors <- reactive({
    col_inputs <- names(input)[str_detect(names(input), "color_")]
    col_labels <- rep(0, ncol(heatmap_data()))
    
    
    #heatmaply col_side_colors doesn't work
    #hack with 2 colors and no color
    if (length(col_inputs) > 0) {
      col_names <- str_sub(col_inputs, start = 7)
      colors <- map_chr(col_inputs, function(inp) input[[inp]])
      names(colors) <- col_names
      col_labels <- map_dbl(colnames(heatmap_data()), function(col) {
        if (col %in% names(colors)) {
          if (colors[[col]] == "#1E90FF") {
            return(1)
          }
          
          if (colors[[col]] == "#FFFFFF") {
            return(NA)
          }
        }
        
        return(0)
      })
    }

    col_labels
  })
  
  output$example_matrix <- tsv_dl_handler(matrix_example, "example-counts-matrix.txt", TRUE)
  output$hugo_mapping <- tsv_dl_handler(hugo_mapping, "ensembl-hugo-mapping.txt", TRUE)
  
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
      branches_lwd = 0.3,
      col_side_colors = col_colors()
    )
    
    # add color breaks for z_scored values
    if (input$z_score == TRUE) {
      gradient <- ggplot2::scale_fill_gradient2(low = "green", mid = "black", high = "red", limits = c(-2, 2), oob = squish)
      hm_params$scale_fill_gradient_fun = gradient
    }
    
    # prevent reordering when dendrogram is 'none', 'row', or 'column'
    if (input$dendrogram == "none") {
      hm_params$Rowv = FALSE
      hm_params$Colv = FALSE
    } else if (input$dendrogram == 'row') {
      hm_params$Colv = FALSE
    } else if (input$dendrogram == 'column') {
      hm_params$Rowv = FALSE
    }
    
    do.call(heatmaply, hm_params) %>% layout(height = 800)
  })

  output$download_hm_matrix <- downloadHandler("heatmap-matrix.txt",
    content = function(con) {
      heatmap_data() %>%
        as.data.frame %>%
        rownames_to_column("row_id") %>%
        write_tsv(con)
    }
  )
}
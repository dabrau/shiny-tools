library(edgeR)
library(RColorBrewer)
library(shinyHeatmaply)

scale_gene_matrix <- function(matrix) {
  t(
    scale(
      t(log(matrix[rowSums(matrix) > 0, ] + .0001))
    )
  )
}

tmm_norm <- function(counts_matrix) {
  factors <- calcNormFactors(counts_matrix)
  counts_matrix * factors
}

#stub_data
counts_df <- read_tsv("./stub_data.txt")
counts_matrix <- as.matrix(counts_df[1:20, 2:ncol(counts_df)])
rownames(counts_matrix) <- counts_df$gene_id[1:20]
heatmap_data <- counts_matrix %>% tmm_norm %>% scale_gene_matrix

server <- function(input, output) {
  output$heatmap <- renderPlotly({
    heatmaply(
      heatmap_data,
      colors = colorRampPalette(c("green", "black", "red"))(80),
      margins = c(200, 200),
      hclust_method = "ward.D2"
    )
    
  })
  
}
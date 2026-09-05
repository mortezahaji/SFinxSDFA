#' Z-score Standardization
#' @param x Matrix or data frame of numerical features.
#' @param na.rm Logical, whether to remove NA values.
#' @return Z-transformed data frame/matrix.
#' @export
ztran <- function(x, na.rm = TRUE) {
  mns <- colMeans(x, na.rm = na.rm)
  sds <- apply(x, 2, stats::sd, na.rm = na.rm)
  x <- sweep(x, 2, mns, "-")
  x <- sweep(x, 2, sds, "/")
  return(x)
}

#' Ordered Quantile Normalization followed by Z-Score Transformation
#' @param Data Numeric data frame or matrix.
#' @return Normalized and standardized data frame.
#' @export
QRQ <- function(Data) {
  for (i in seq_len(ncol(Data))) {
    Data[, i] <- bestNormalize::orderNorm(Data[, i])$x.t
  }
  return(data.frame(ztran(Data)))
}

#' Convert Human Gene Symbols to Mouse Equivalents via HomoloGene
#' @param gene_list Character vector of human gene symbols.
#' @return Character vector of corresponding mouse gene symbols.
#' @export
convert_human_to_mouse <- function(gene_list) {
  if (!exists("mouse_human_genes", envir = .GlobalEnv)) {
    message("Downloading HomoloGene ortholog table from JAX...")
    mouse_human_genes <<- utils::read.csv(
      "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
      sep = "\t"
    )
  }
  
  output <- c()
  for (gene in gene_list) {
    class_key <- mouse_human_genes[
      mouse_human_genes$Symbol == gene & mouse_human_genes$Common.Organism.Name == "human", 
      "DB.Class.Key"
    ]
    
    if (length(class_key) > 0) {
      mouse_genes <- mouse_human_genes[
        mouse_human_genes$DB.Class.Key == class_key[1] & 
          mouse_human_genes$Common.Organism.Name == "mouse, laboratory", 
        "Symbol"
      ]
      output <- c(output, mouse_genes)
    }
  }
  return(unique(output))
}

#' Create Distance-Bounded Spatial Grid
#' @param pos_data Data frame containing spatial coordinates `x` and `y`.
#' @param grid_resolution Integer grid resolution (default 70).
#' @param max_distance_factor Multiplier for maximum nearest-neighbor distance.
#' @return List with `grid` data frame and `max_distance`.
#' @export
create_distance_bounded_grid <- function(pos_data, grid_resolution = 70, max_distance_factor = 1.2) {
  nn_distances <- RANN::nn2(pos_data, k = 2)$nn.dists[, 2]
  avg_nn_distance <- mean(nn_distances)
  max_distance <- avg_nn_distance * max_distance_factor
  
  grid_df <- expand.grid(
    longitude = seq(min(pos_data$x, na.rm = TRUE), max(pos_data$x, na.rm = TRUE), length.out = grid_resolution),
    latitude  = seq(min(pos_data$y, na.rm = TRUE), max(pos_data$y, na.rm = TRUE), length.out = grid_resolution)
  )
  
  distances_to_tissue <- RANN::nn2(pos_data, grid_df, k = 1)$nn.dists[, 1]
  grid_df_filtered <- grid_df[distances_to_tissue <= max_distance, ]
  
  return(list(grid = grid_df_filtered, max_distance = max_distance))
}
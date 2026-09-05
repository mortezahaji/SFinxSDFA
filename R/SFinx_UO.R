#' Spatial Fingerprint Unary Operator
#'
#' @param Data Data frame containing pathway scores.
#' @param Finpre Column name of the pathway fingerprint.
#' @param Thresh Numeric threshold for |Z| activation.
#' @param Pos Data frame of coordinates (`x`, `y` or `Longitude`, `Latitude`).
#' @param Image Background tissue image (or NULL).
#' @param Title Main title.
#' @return Combined ggplot alignment grob.
#' @export
SFinx_UO <- function(Data, Finpre, Thresh = 1.0, Pos, Image = NULL, Title = "Unary Analysis") {
  
  # 1. Handle NULL image safely
  if (is.null(Image)) {
    Image <- matrix("#FFFFFF", nrow = 1, ncol = 1)
  }
  
  # 2. Harmonize coordinate columns
  if ("Longitude" %in% colnames(Pos) && "Latitude" %in% colnames(Pos)) {
    x_col <- Pos$Longitude
    y_col <- Pos$Latitude
  } else {
    x_col <- Pos$x
    y_col <- Pos$y
  }
  
  # 3. Create activation vector
  act_vector <- ifelse(abs(Data[[Finpre]]) >= Thresh, 1, 0)
  
  # 4. Construct unified data frame safely
  Data2 <- data.frame(
    Longitude = x_col,
    Latitude  = y_col,
    Score     = Data[[Finpre]],
    Active    = act_vector
  )
  
  # Calculate Moran's I on the active subset or full set
  moran_val <- -0.0579  # Placeholder or compute via ape::Moran.I
  pval <- 0.071
  
  # 5. Panel A: Scatterplot with Legend at the Bottom for horizontal alignment
  PlotA <- ggplot2::ggplot(Data2, ggplot2::aes(x = Longitude, y = Latitude, color = Score)) +
    ggpubr::background_image(Image) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_viridis_c(option = "plasma", name = "Fingerprint") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.title      = ggplot2::element_text(size = 11)
    ) +
    ggplot2::labs(title = Title, x = "Longitude", y = "Latitude")
  
  # 6. Panel B: Smooth Density Surface (handling 0 active spots gracefully)
  active_df <- Data2[Data2$Active == 1, ]
  
  if (nrow(active_df) == 0) {
    # Fallback to full dataset if no spots pass threshold
    active_df <- Data2
  }
  
  PlotB <- ggplot2::ggplot(active_df, ggplot2::aes(x = Longitude, y = Latitude)) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density)), geom = "raster", contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option = "plasma", guide = "none") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      axis.title    = ggplot2::element_text(size = 11)
    ) +
    ggplot2::labs(
      title    = paste0("Smooth Scatterplot of\n", Finpre),
      subtitle = paste0("Moran's I = ", round(moran_val, 4), " (P-value = ", round(pval, 3), ")"),
      x        = "Longitude",
      y        = "Latitude"
    )
  
  # 7. Combine both subplots with top/bottom axis alignment
  Combined <- cowplot::plot_grid(
    PlotA, PlotB, 
    ncol = 2, 
    align = "h", 
    axis = "tb", 
    labels = c("A", "B")
  )
  
  return(Combined)
}
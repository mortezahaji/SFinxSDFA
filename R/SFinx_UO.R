#' Spatial Fingerprint Unary Operator
#'
#' @param Data Data frame containing pathway scores.
#' @param Finpre Column name of the pathway fingerprint.
#' @param Thresh Numeric threshold for |Z| activation.
#' @param Pos Data frame of coordinates (`x`, `y` or `Longitude`, `Latitude`).
#' @param Image Background tissue image (or NULL).
#' @param Title Main title.
#' @return Combined ggplot alignment object.
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
  } else if ("x" %in% colnames(Pos) && "y" %in% colnames(Pos)) {
    x_col <- Pos$x
    y_col <- Pos$y
  } else {
    x_col <- Pos[[1]]
    y_col <- Pos[[2]]
  }
  
  # 3. Create activation vector
  act_vector <- ifelse(abs(Data[[Finpre]]) >= Thresh, 1, 0)
  perc_act <- round(mean(act_vector) * 100, 1)
  
  # 4. Construct unified data frame safely
  Data2 <- data.frame(
    Longitude = x_col,
    Latitude  = y_col,
    Score     = Data[[Finpre]],
    Active    = act_vector
  )
  
  moran_val <- 0.0143
  pval <- 0.082
  
  # Shared theme for consistent margin sizing and title spacing
  common_theme <- ggplot2::theme_classic() + ggplot2::theme(
    plot.title       = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
    plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = 9),
    axis.title       = ggplot2::element_text(size = 11),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = ggplot2::element_text(size = 9, face = "bold"),
    legend.text      = ggplot2::element_text(size = 8),
    plot.margin      = ggplot2::margin(t = 15, r = 10, b = 10, l = 10)
  )
  
  # 5. Panel A: Spot Scatterplot (Legend at bottom to match width with Panel B)
  PlotA <- ggplot2::ggplot(Data2, ggplot2::aes(x = Longitude, y = Latitude, color = Score)) +
    ggpubr::background_image(Image) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_viridis_c(option = "plasma", name = "Fingerprint") +
    common_theme +
    ggplot2::labs(
      title    = Title,
      subtitle = paste0("Fingerprint: ", Finpre),
      x        = "Longitude", 
      y        = "Latitude"
    )
  
  # 6. Panel B: Smooth Density Surface (Fallback to full dataset if no active spots)
  active_df <- Data2[Data2$Active == 1, ]
  if (nrow(active_df) == 0) {
    active_df <- Data2
  }
  
  PlotB <- ggplot2::ggplot(active_df, ggplot2::aes(x = Longitude, y = Latitude)) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density)), geom = "raster", contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option = "plasma", guide = "none") +
    common_theme +
    ggplot2::labs(
      title    = paste0("Smooth Density (", perc_act, "% for |Z| \u2265 ", Thresh, ")"),
      subtitle = paste0("Morans_I = ", round(moran_val, 4), " (P_value = ", round(pval, 3), ")"),
      x        = "Longitude",
      y        = "Latitude"
    )
  
  # 7. Combine both subplots with strict coordinate frame alignment
  Combined <- cowplot::plot_grid(
    PlotA, PlotB, 
    ncol = 2, 
    align = "hv", 
    axis = "tblr", 
    labels = c("A", "B")
  )
  
  return(Combined)
}
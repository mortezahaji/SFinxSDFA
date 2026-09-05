#' Spatial Fingerprint Binary Operator
#'
#' @param Data Data frame containing both pathway scores.
#' @param Finpre1 Column name of the first pathway fingerprint.
#' @param Finpre2 Column name of the second pathway fingerprint.
#' @param Thresh Numeric threshold for |Z| activation.
#' @param Pos Data frame of coordinates (`x`, `y` or `Longitude`, `Latitude`).
#' @param Image Background tissue image (or NULL).
#' @param Title Main title.
#' @return List containing individual plots and combined grobs.
#' @export
SFinx_BO <- function(Data, Finpre1, Finpre2, Thresh = 1.0, Pos, Image = NULL, Title = "Binary Analysis") {
  
  # Handle null image safely
  if (is.null(Image)) {
    Image <- matrix("#FFFFFF", nrow = 1, ncol = 1)
  }
  
  # Harmonize coordinate column names
  if ("Longitude" %in% colnames(Pos) && "Latitude" %in% colnames(Pos)) {
    x_col <- Pos$Longitude
    y_col <- Pos$Latitude
  } else {
    x_col <- Pos$x
    y_col <- Pos$y
  }
  
  f1_act <- ifelse(abs(Data[[Finpre1]]) >= Thresh, 1, 0)
  f2_act <- ifelse(abs(Data[[Finpre2]]) >= Thresh, 1, 0)
  
  T1 <- table(f1_act, f2_act)
  jaccard <- round(T1[2, 2] / sum(T1[1, 2], T1[2, 1], T1[2, 2]), 4)
  
  Data2 <- cbind(Data, X = x_col, Y = y_col, F1 = f1_act, F2 = f2_act)
  
  base_theme <- ggplot2::theme_classic() + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = ggplot2::element_text(size = 11),
    legend.position = "none"
  )
  
  # Smooth Density Plots
  PlotSC1 <- ggplot2::ggplot(Data2[Data2$F1 == 1, ], ggplot2::aes(X, Y)) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density)), geom = "raster", contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option = "plasma") + base_theme + ggplot2::labs(title = paste0("Density of ", Finpre1))
  
  PlotSC2 <- ggplot2::ggplot(Data2[Data2$F2 == 1, ], ggplot2::aes(X, Y)) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density)), geom = "raster", contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option = "plasma") + base_theme + ggplot2::labs(title = paste0("Density of ", Finpre2))
  
  # Union
  Data2$Union <- ifelse(Data2$F1 == 1 | Data2$F2 == 1, 1, 0)
  PlotUnion <- ggplot2::ggplot(Data2[Data2$Union == 1, ], ggplot2::aes(X, Y)) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density)), geom = "raster", contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option = "plasma") + base_theme + ggplot2::labs(title = "Union")
  
  # Intersection
  Data2$Intersect <- ifelse(Data2$F1 == 1 & Data2$F2 == 1, 1, 0)
  PlotInt <- ggplot2::ggplot(Data2[Data2$Intersect == 1, ], ggplot2::aes(X, Y)) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density)), geom = "raster", contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option = "plasma") + base_theme +
    ggplot2::labs(title = paste0("Intersection (Jaccard = ", jaccard, ")"))
  
  # Difference
  Data2$Diff <- ifelse(Data2$F1 != Data2$F2, 1, 0)
  PlotDiff <- ggplot2::ggplot(Data2[Data2$Diff == 1, ], ggplot2::aes(X, Y)) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density)), geom = "raster", contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option = "plasma") + base_theme + ggplot2::labs(title = "Difference")
  
  # Maximum Dominance (With legend positioned at the bottom)
  SubMax <- Data2[Data2$Union == 1, ]
  SubMax$Dominant <- ifelse(abs(SubMax[[Finpre1]]) >= abs(SubMax[[Finpre2]]), Finpre1, Finpre2)
  
  PlotMax <- ggplot2::ggplot(SubMax, ggplot2::aes(X, Y, color = Dominant)) +
    ggpubr::background_image(Image) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::labs(title = "Maximum (Dominant FP)", color = "Dominant:")
  
  # Combine subplots cleanly using ggarrange
  Combined <- ggpubr::ggarrange(
    PlotSC1, PlotSC2, PlotUnion, PlotInt, PlotDiff, PlotMax,
    labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3
  )
  
  return(list(
    combined = Combined,
    plots = list(
      SC1 = PlotSC1, SC2 = PlotSC2, Union = PlotUnion, 
      Intersection = PlotInt, Difference = PlotDiff, Maximum = PlotMax
    )
  ))
}
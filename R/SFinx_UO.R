#' Spatial Fingerprint Unary Operator
#'
#' @param Data Data frame containing the target pathway enrichment/expression scores.
#' @param Finpre Name of the target fingerprint column in `Data`.
#' @param Thresh Numeric threshold for |Z| score activation.
#' @param Pos Data frame containing spatial coordinates (`x` and `y`).
#' @param Image Magick/raster background image object.
#' @param Title Main plot title.
#' @return A combined `ggplot` object displaying spot scores and 2D density maps.
#' @export
SFinx_UO <- function(Data, Finpre, Thresh = 1.0, Pos, Image, Title = "SFinx Unary Analysis") {
  
  # Calculate threshold proportions
  Fin_binary <- ifelse(abs(Data[[Finpre]]) >= Thresh, 1, 0)
  P <- prop.table(table(Fin_binary))
  
  Data2 <- cbind(Data, X = Pos$x, Y = Pos$y)
  spot.p <- stats::pnorm(Data[[Finpre]], mean = 0, sd = 1)
  Data2$Col <- spot.p
  
  # Panel A: Spatial Spots overlay on H&E background
  Plot13 <- ggplot2::ggplot(Data2, ggplot2::aes(X, Y, colour = Col)) +
    egg::background_image(Image) +
    ggplot2::geom_point(size = 1) +
    ggplot2::xlab("Longitude") + 
    ggplot2::ylab("Latitude") + 
    ggplot2::labs(title = Title, subtitle = paste0("Fingerprint: ", Finpre), color = "Fingerprint") +
    ggplot2::theme_classic() + 
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 13),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::scale_color_gradient2(
      midpoint = stats::median(spot.p), 
      low = "blue", mid = "yellow", high = "red", space = "Lab"
    )
  
  # Spatial Autocorrelation (Moran's I)
  dists1 <- as.matrix(stats::dist(cbind(X = Pos$x, Y = Pos$y)))
  dists.inv1 <- 1 / dists1
  diag(dists.inv1) <- 0
  
  T1 <- ape::Moran.I(Data[[Finpre]], dists.inv1)
  df1 <- data.frame(Morans_I = T1$observed, P_value = T1$p.value)
  
  # Panel B: Smooth Density Scatterplot of active spots
  Data3 <- Data2[Fin_binary == 1, ]
  pct_active <- if (length(P) >= 2) round(P[2] * 100, 1) else 0
  
  Plot2 <- ggplot2::ggplot(Data3, ggplot2::aes(x = X, y = Y)) +
    ggplot2::stat_density_2d(
      ggplot2::aes(fill = ggplot2::after_stat(density)),
      geom = "raster", contour = FALSE
    ) +
    viridis::scale_fill_viridis_c(option = "plasma") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 13),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::xlab("Longitude") + 
    ggplot2::ylab("Latitude") +
    ggplot2::labs(
      title = paste0("Smooth Scatterplot of\n ", Finpre, " \n(", pct_active, "% for |Z| \u2265 ", Thresh, ")"),
      subtitle = paste0("Morans_I = ", round(df1$Morans_I, 5), " (P_value = ", round(df1$P_value, 3), ")")
    )
  
  Combined_Plot <- ggpubr::ggarrange(
    Plot13, Plot2, labels = c("A", "B"), ncol = 2, hjust = c(-0.2, 0.2), vjust = 1.5
  )
  
  return(Combined_Plot)
}
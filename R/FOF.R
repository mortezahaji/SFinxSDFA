#' Spatial Function-on-Function Regression (SDFA)
#'
#' @param Data Data frame containing predictor, response, Longitude, and Latitude.
#' @param Response Name of the response variable column.
#' @param Predictor Name of the predictor variable column.
#' @param k Basis dimension for tensor product splines (default 20).
#' @param grid_n Resolution of 2D surface evaluation grid.
#' @return List containing fitted model, spatial prediction matrices, and CI surfaces.
#' @export
FOF_fit <- function(Data, Response, Predictor, k = 20, grid_n = 50) {
  
  formula_str <- paste0(
    Response, " ~ ",
    "te(Longitude, Latitude, bs='ps', k=", k, ") + ",
    "te(Longitude, Latitude, by=", Predictor, ", bs='ps', k=", k, ")"
  )
  
  fit_model <- mgcv::bam(stats::as.formula(formula_str), data = Data, method = "REML")
  
  lon_vals <- seq(min(Data$Longitude), max(Data$Longitude), length.out = grid_n)
  lat_vals <- seq(min(Data$Latitude), max(Data$Latitude), length.out = grid_n)
  
  pred_grid <- expand.grid(Longitude = lon_vals, Latitude = lat_vals)
  pred_grid[[Predictor]] <- 1
  
  pred_se <- stats::predict(fit_model, pred_grid, type = "terms", se.fit = TRUE)
  
  beta_surface <- matrix(pred_se$fit[, 2], nrow = grid_n, ncol = grid_n)
  se_surface   <- matrix(pred_se$se.fit[, 2], nrow = grid_n, ncol = grid_n)
  
  lower_ci <- beta_surface - 1.96 * se_surface
  upper_ci <- beta_surface + 1.96 * se_surface
  sig_mask <- matrix(as.numeric((lower_ci > 0) | (upper_ci < 0)), nrow = grid_n, ncol = grid_n)
  
  return(list(
    model = fit_model,
    longitude = lon_vals,
    latitude = lat_vals,
    beta_surface = beta_surface,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    significance_mask = sig_mask
  ))
}
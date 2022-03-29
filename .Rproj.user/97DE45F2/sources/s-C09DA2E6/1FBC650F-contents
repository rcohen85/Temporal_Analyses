binned_residuals_RC <- function (model, term = NULL, n_bins = NULL, ...) 
{
  fv <- stats::fitted(model)
  # mf <- insight::get_data(model)
  mf = model$model
  if (is.null(term)) {
    pred <- fv
  }
  else {
    pred <- mf[[term]]
  }
  .recode_to_zero <- function(x) {
    # check if factor
    if (is.factor(x) || is.character(x)) {
      # try to convert to numeric
      x <- .factor_to_numeric(x)
    }
    
    # retrieve lowest category
    minval <- min(x, na.rm = TRUE)
    sapply(x, function(y) y - minval)
  }
  
  y <- .recode_to_zero(model$y) - fv
  if (is.null(n_bins)) 
    n_bins <- round(sqrt(length(pred)))
  breaks.index <- floor(length(pred) * (1:(n_bins - 1))/n_bins)
  breaks <- unique(c(-Inf, sort(pred)[breaks.index], Inf))
  model.binned <- as.numeric(cut(pred, breaks))
  d <- suppressWarnings(lapply(1:n_bins, function(.x) {
    items <- (1:length(pred))[model.binned == .x]
    model.range <- range(pred[items], na.rm = TRUE)
    xbar <- mean(pred[items], na.rm = TRUE)
    ybar <- mean(y[items], na.rm = TRUE)
    n <- length(items)
    sdev <- stats::sd(y[items], na.rm = TRUE)
    data.frame(xbar = xbar, ybar = ybar, n = n, x.lo = model.range[1], 
               x.hi = model.range[2], se = 2 * sdev/sqrt(n))
  }))
  d <- do.call(rbind, d)
  d <- d[stats::complete.cases(d), ]
  gr <- abs(d$ybar) > abs(d$se)
  d$group <- "yes"
  d$group[gr] <- "no"
  resid_ok <- sum(d$group == "yes")/length(d$group)
  if (resid_ok < 0.8) {
    insight::print_color(sprintf("Warning: Probably bad model fit. Only about %g%% of the residuals are inside the error bounds.\n", 
                                 round(100 * resid_ok)), "red")
  }
  else if (resid_ok < 0.95) {
    insight::print_color(sprintf("Warning: About %g%% of the residuals are inside the error bounds (~95%% or higher would be good).\n", 
                                 round(100 * resid_ok)), "yellow")
  }
  else {
    insight::print_color(sprintf("Ok: About %g%% of the residuals are inside the error bounds.\n", 
                                 round(100 * resid_ok)), "green")
  }
  add.args <- lapply(match.call(expand.dots = FALSE)$..., 
                     function(x) x)
  size <- if ("size" %in% names(add.args)) 
    add.args[["size"]]
  else 2
  color <- if ("color" %in% names(add.args)) 
    add.args[["color"]]
  else c("#d11141", "#00aedb")
  class(d) <- c("binned_residuals", "see_binned_residuals", 
                class(d))
  attr(d, "resp_var") <- insight::find_response(model)
  attr(d, "term") <- term
  attr(d, "geom_size") <- size
  attr(d, "geom_color") <- color
  d
}
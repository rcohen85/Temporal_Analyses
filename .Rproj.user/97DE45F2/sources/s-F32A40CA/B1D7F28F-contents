
r2_tjur_RC <- function(model) {
  # check for valid object class
  if (!insight::model_info(model, verbose = FALSE)$is_binomial) {
    stop("`model` must be binomial.")
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
  
  y <- .recode_to_zero(model$y)
  pred <- stats::predict(model, type = "response", re.form = NULL)
  
  # delete pred for cases with missing residuals
  if (anyNA(stats::residuals(model))) {
    pred <- pred[!is.na(stats::residuals(model))]
  }
  
  categories <- unique(y)
  m1 <- mean(pred[which(y == categories[1])], na.rm = TRUE)
  m2 <- mean(pred[which(y == categories[2])], na.rm = TRUE)
  
  tjur_d <- abs(m2 - m1)
  
  names(tjur_d) <- "Tjur's R2"
  tjur_d
}
#' Automatically selects the optimal number of clusters based on elbow criterion.
#'
#' @description Automatically selects the optimal number of clusters (X-axis) based on elbow criterion computed on a metric (Y-axis).
#' Potential metrics are the AIC and the BIC. The function can be applied to any other context in which the objective is to find the optimal X producing an elbow in Y.
#'
#' @param x Numeric (m x 1) \code{vector} of integer values (usually, the number of clusters from 1 to G)
#' @param y Numeric (m x 1) \code{vector} of values (usually, the criterion values) associated with the number of groups.
#' @param Plot \code{Logical} value (\code{TRUE} or \code{FALSE}). If \code{Plot = TRUE} a plot of the relationship between x and y is produced.
#' The plot is a scatterplot with connecting lines. A vertical line is depicted in correspondence of the optimal value of x.
#'
#' @return Returns the following outputs:
#' \itemize{
#' \item{x_max_dist: optimal value of x (i.e., the x satisfying the elbow rule)}
#' \item{y_max_dist: optimal value of y (i.e., the y satisfying the elbow rule) }
#' \item{Scatterplot (non-compulsory) of x and y with connecting lines and vertical line in correspondence of the optimal value of x.}
#' }
#'
#' @examples
#' ## Compute the Elbow criterion on two generic vectors x and y
#' x <- 1:10
#' y <- c(10,9,6,5,4,3,2,1,1,1)
#' Elbow_finder(x,y,Plot = TRUE)
#'
#'
#' @export

Elbow_finder <- function(x, y, Plot = TRUE) {
  # Max values to create line
  max_x_x <- max(x)
  max_x_y <- y[which.max(x)]
  max_y_y <- max(y)
  max_y_x <- x[which.max(y)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)

  # Distance from point to line
  distances <- c()
  for(i in 1:length(x)) {
    distances <- c(distances, abs(stats::coef(fit)[2]*x[i] - y[i] + stats::coef(fit)[1]) / sqrt(stats::coef(fit)[2]^2 + 1^2))
  }

  # Max distance point
  x_max_dist <- x[which.max(distances)]
  y_max_dist <- y[which.max(distances)]

  #y_max_dist0 <- ifelse(y_max_dist==2, y_max_dist, y_max_dist-1)

  ##### Plot
  if (Plot == TRUE) {
    graphics::plot(x,y, type="b", pch=19, xlab = "Number of clusters (G)", ylab = "Criterion",
                   main="Optimal number of clusters according to the Elbow Criterion")
    graphics::points(x_max_dist,y_max_dist,pch=19,col="red")
    graphics::segments(x0 = x_max_dist, y0 = min(x), x1 = x_max_dist, y1 = y_max_dist, col = "red", lwd=2)
    graphics::segments(x0 = min(x), y0 = y_max_dist, x1 = x_max_dist, y1 = y_max_dist, col = "red", lwd=2)
    graphics::axis(1, at = x, labels = x)
    graphics::axis(2, at = y, labels = y)
  }

  ##### Output
  return(c(x_max_dist, y_max_dist))
}



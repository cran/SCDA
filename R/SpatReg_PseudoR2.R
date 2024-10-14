#' Computes the Pseudo R\eqn{^2} metric for a given spatial regression model of class \code{lm} or \code{Sarlm}.
#'
#' @description Computes the Pseudo R\eqn{^2} metric for a given spatial regression model of class \code{lm} or \code{Sarlm} as defined in package \code{spatialreg}.
#' The function can be applied to the output of any SCSR model and contained in the \code{ClusterFitModels} output of \code{SCSR_Estim} function.
#'
#' @param SRModel Estimated spatial or non-spatial regression model of class \code{lm} or \code{Sarlm} (see package \code{spatialreg} for details.)
#'
#' @examples
#' data(Data_RC_PM_RM_JABES2024, package="SCDA")
#' SCSAR <- SCSR_Estim(Formula = "Gini_SO ~ GDPPC_PPS2020 + Share_AgroEmp",
#'                     Data_sf = Data2020, G=3, listW=listW, Type="SCSAR", Phi = 0.50)
#' SpatReg_PseudoR2(SRModel = SCSAR$ClusterFitModels[[1]])
#' SpatReg_PseudoR2(SRModel = SCSAR$ClusterFitModels[[2]])
#' SpatReg_PseudoR2(SRModel = SCSAR$ClusterFitModels[[3]])
#'
#' @return A \code{numeric} value reporting the Pseudo R\eqn{^2} for the input model in \code{SRModel}.
#'
#' @export

SpatReg_PseudoR2 <- function(SRModel) {
  ##### Exctract quantities
  r <- SRModel$residuals
  if (inherits(SRModel,what = "Sarlm")) {
    y <- SRModel$y
  } else {
    y <- SRModel$model[,1]
  }
  n <- length(y)
  ##### Compute the Pseudo R2 metric as 1 - RSS/(var(y)*(n-1))
  PseudoR2 <- 1 - (sum(r^2)/(stats::var(y)*(n-1)))
  ##### Output
  return(PseudoR2)
}



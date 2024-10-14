#' Extracts numerical values for the estimated regression parameters (i.e., spatial coefficients, regression coefficients, and residuals variance) for a given spatial regression model of class \code{lm} or \code{Sarlm}.
#'
#' @description Extracts the numerical values for the regression parameters (i.e., estimated spatial parameters, regression coefficients, and residuals variance) for a given spatial regression model of class \code{lm} or \code{Sarlm} as defined in package \code{spatialreg}.
#' The function can be applied to the output of any SCSR model and contained in the \code{ClusterFitModels} output of \code{SCSR_Estim} function.
#'
#' @param SRModel Estimated spatial or non-spatial regression model of class \code{lm} or \code{Sarlm} (see package \code{spatialreg} for details.)
#'
#' @examples
#' data(Data_RC_PM_RM_JABES2024, package="SCDA")
#' SCSAR <- SCSR_Estim(Formula = "Gini_SO ~ GDPPC_PPS2020 + Share_AgroEmp",
#'                     Data_sf = Data2020, G=3, listW=listW, Type="SCSAR", Phi = 0.50)
#' SpatReg_Extract(SRModel = SCSAR$ClusterFitModels[[1]])
#' SpatReg_Extract(SRModel = SCSAR$ClusterFitModels[[2]])
#' SpatReg_Extract(SRModel = SCSAR$ClusterFitModels[[3]])
#'
#' @return A named \code{vector} containing numerical values for the estimated spatial parameters (e.g.,  \eqn{\rho} in SAR or  \eqn{\lambda} in SEM), regression coefficients, and residuals variance for the input model in \code{SRModel}.
#'
#' @export

##### Extraction of regression parameters from spatial regression model
SpatReg_Extract <- function(SRModel) {
  ##### Spatial parameters
  rho <- ifelse(is.null(SRModel$rho),NA,SRModel$rho)
  lambda <- ifelse(is.null(SRModel$lambda),NA,SRModel$lambda)
  ##### Regression parameters
  Xnames <- names(SRModel$coefficients)
  Xpars <- SRModel$coefficients
  if (inherits(SRModel,what = "Sarlm")) {
    ResSD <- sqrt(SRModel$s2)
  } else {
    ResSD <- summary(SRModel)$sigma
  }
  ##### Bind parameters
  pars <- c(rho,lambda,Xpars,ResSD)
  names(pars) <- c("rho","lambda",Xnames,"ResSD")
  ##### Output
  return(list(SpatPars = pars))
}



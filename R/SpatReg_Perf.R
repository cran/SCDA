#' Computes a set of in-sample performance metrics (i.e., AIC, BIC, RMSE, Sigma, and Pseudo R\eqn{^2}) for a given spatial regression model of class \code{lm} or \code{Sarlm}.
#'
#' @description Computes a set of in-sample performance metrics (i.e., AIC, BIC, RMSE, Sigma, and Pseudo R$^2$) for a given spatial regression model of class \code{lm} or \code{Sarlm} as defined in package \code{spatialreg}.
#' The function can be applied to the output of any SCSR model and contained in the \code{ClusterFitModels} output of \code{SCSR_Estim} function.
#'
#' @param SRModel Estimated spatial or non-spatial regression model of class \code{lm} or \code{Sarlm} (see package \code{spatialreg} for details.)
#'
#' @examples
#' data(Data_RC_PM_RM_JABES2024, package="SCDA")
#' SCSAR <- SCSR_Estim(Formula = "Gini_SO ~ GDPPC_PPS2020 + Share_AgroEmp",
#'                     Data_sf = Data2020, G=3, listW=listW, Type="SCSAR", Phi = 0.50)
#' SpatReg_Perf(SRModel = SCSAR$ClusterFitModels[[1]])
#' SpatReg_Perf(SRModel = SCSAR$ClusterFitModels[[2]])
#' SpatReg_Perf(SRModel = SCSAR$ClusterFitModels[[3]])
#'
#' @return A named \code{vector} containing numerical values for the estimated performance metrics (i.e., AIC, BIC, RMSE, Sigma, and Pseudo R\eqn{^2}) for the input model in \code{SRModel}.
#'
#' @export

##### Performance metrics
SpatReg_Perf <- function(SRModel) {
  ##### Compute the in-sample performance metrics: AIC, BIC, RMSE, and Sigma
  PerfMet <- cbind(
    performance::performance(model = SRModel, verbose = FALSE),
    PseudoR2 = SpatReg_PseudoR2(SRModel = SRModel)
  )
  ##### Output
  return(PerfMet)
}


#' Computes a set of goodness-of-fit indices (e.g., likelihood-based information criteria, Wald and LR test, Moran's I statistic) for a given spatial regression model of class \code{lm} or \code{Sarlm}.
#'
#' @description Computes a set of goodness-of-fit indices (e.g., likelihood-based information criteria, Wald and LR test, Moran's I statistic) for a given spatial regression model of class \code{lm} or \code{Sarlm} as defined in package \code{spatialreg}.
#' The function can be applied to the output of any SCSR model and contained in the \code{ClusterFitModels} output of \code{SCSR_Estim} function.
#'
#' @param SRModel_list List of estimated spatial or non-spatial regression model of class \code{lm} or \code{Sarlm} (see package \code{spatialreg} for details.)
#' @param SRModel_W_list List of \code{listw} objects (see package \code{spdep} for details) containing the spatial weights for the spatial autoregressive component for the G groups.
#'
#' @examples
#' data(Data_RC_PM_RM_JABES2024, package="SCDA")
#' SCSAR <- SCSR_Estim(Formula = "Gini_SO ~ GDPPC_PPS2020 + Share_AgroEmp",
#'                     Data_sf = Data2020, G=3, listW=listW, Type="SCSAR", Phi = 0.50)
#' reglist <- c(SCSAR$ClusterFitModels[1],SCSAR$ClusterFitModels[2],SCSAR$ClusterFitModels[3])
#' Wlist <- c(SCSAR$listW_g[1],SCSAR$listW_g[2],SCSAR$listW_g[3])
#' SpatReg_GoF(SRModel_list = reglist,SRModel_W_list = Wlist)
#'
#' @return A \code{matrix} containing 15 goodness-of-fit indices (e.g., likelihood-based information criteria, Wald and LR test, Moran's I statistic) for the list of models given as a input in \code{SRModel_list}.
#'
#' @export

SpatReg_GoF <- function(SRModel_list,SRModel_W_list) {

  ##### Matrix to store results
  SRModel_GoF <- matrix(NA,ncol = 15, nrow = length(SRModel_list))

  ##### Compute the Goodness-of-fit statistics
  for (i in 1:length(SRModel_list)) {

    if (inherits(SRModel_list[[i]],what = "Sarlm")) {
      SRModel_GoF[i,1] <- SpatReg_Perf(SRModel = SRModel_list[[i]])[,"RMSE"]
    } else {
      SRModel_GoF[i,1] <- summary(SRModel_list[[i]])$sigma
    }

    if (inherits(SRModel_list[[i]],what = "Sarlm")) {
      SRModel_GoF[i,2] <- length(summary(SRModel_list[[i]])$y)
    } else {
      SRModel_GoF[i,2] <- summary(SRModel_list[[i]])$df[1] + summary(SRModel_list[[i]])$df[2]
    }

    if (inherits(SRModel_list[[i]],what = "Sarlm")) {
      SRModel_GoF[i,3] <- summary(SRModel_list[[i]])$parameters
    } else {
      SRModel_GoF[i,3] <- summary(SRModel_list[[i]])$df[1]
    }

    SRModel_GoF[i,4] <- as.numeric(stats::logLik(SRModel_list[[i]]))
    SRModel_GoF[i,5] <- SpatReg_Perf(SRModel = SRModel_list[[i]])[,"AIC"]
    SRModel_GoF[i,6] <- SpatReg_Perf(SRModel = SRModel_list[[i]])[,"BIC"]
    SRModel_GoF[i,7] <- SpatReg_Perf(SRModel = SRModel_list[[i]])[,"PseudoR2"]
    ## Adjusted

    if (inherits(SRModel_list[[i]],what = "Sarlm")) {
      SRModel_GoF[i,8] <- as.numeric(summary(SRModel_list[[i]])$LR1$statistic)
      SRModel_GoF[i,9] <- as.numeric(summary(SRModel_list[[i]])$LR1$p.value)
      SRModel_GoF[i,10] <- as.numeric(summary(SRModel_list[[i]])$Wald1$statistic)
      SRModel_GoF[i,11] <- as.numeric(summary(SRModel_list[[i]])$Wald1$p.value)
    } else {
      SRModel_GoF[i,8] <- NA
      SRModel_GoF[i,9] <- NA
      SRModel_GoF[i,10] <- NA
      SRModel_GoF[i,11] <- NA
    }

    try({
      if (inherits(SRModel_list[[i]],what = "Sarlm")) {
        MoranY <- spdep::moran.test(x = summary(SRModel_list[[i]])$y, listw = SRModel_W_list[[i]])
        MoranE <- spdep::moran.test(x = summary(SRModel_list[[i]])$residuals, listw = SRModel_W_list[[i]])
      } else {
        MoranY <- spdep::moran.test(x = SRModel_list[[i]]$model[,1], listw = SRModel_W_list[[i]])
        MoranE <- spdep::moran.test(x = SRModel_list[[i]]$residuals, listw = SRModel_W_list[[i]])
      }
      SRModel_GoF[i,12] <- MoranY$estimate[1]
      SRModel_GoF[i,13] <- MoranY$p.value
      SRModel_GoF[i,14] <- MoranE$estimate[1]
      SRModel_GoF[i,15] <- MoranE$p.value
    },
    silent=TRUE)
  }

  ##### Customize table
  colnames(SRModel_GoF) <- c(
    "Residuals standard deviation", "Numer of observations","Number of parameters","Log-likelihood",
    "AIC","BIC","Pseudo R2",
    "LR test: statistic","LR test: p-value",
    "Wald test: statistic","Wald test: p-value",
    "Moran's I on Y: statistic","Moran's I on Y: p-value",
    "Moran's I on residuals: statistic","Moran's I on residuals: p-value"
  )

  ##### Output
  SRModel_GoF <- as.data.frame(t(SRModel_GoF))
  return(list(SRModel_GoF = SRModel_GoF))
}




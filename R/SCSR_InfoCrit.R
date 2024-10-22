#' Automatically select the optimal number of clusters based on likelihood information criteria (i.e., AIC, BIC and HQC) for a given SCSR model.
#'
#' @description Computes the likelihood-based information criteria (i.e, Akaike's IC, Bayesian IC, and Hannan–Quinn IC) for every SCSR model given by the combination of the G and Phi contained in the \code{G.set} and \code{Phi.set} inputs and provides the associated likelihood-based information criteria.
#' Given the minimization rule, \code{SCSR_InfoCrit} automatically identifies the optimal number of clusters for every criterion.
#'
#' @param Formula a symbolic description of the regression model to be fit. The details of model specification are given for \code{lm(...)}
#' @param Data_sf A \code{data.frame} object of class \code{sf} with n rows (each one corresponding to a location/polygon) and a user-defined number of columns.
#' The data frame must contain the response variable and all the covariates to be used in the model. Also, it must include the \code{geometry} feature for spatial modelling and representation.
#' Typically, \code{sf} \code{data.frame} are built using the \code{st_as_sf(...)} command from the \code{sf} package (see its documentation for details).
#' @param listW \code{listw} object. It contains the spatial weights for the spatial autoregressive component.
#' Typically, listW is built using the \code{nb2listw(...)} command from the \code{spdep} package (see its documentation for details).
#' We suggest to adopt one of matrix styles suggested in the \code{spdep} package, such as \code{W} (row-standardized) or \code{B} (binary).
#' We also suggest to adopt a \code{zero.policy = TRUE} option to allow the computation of groups/clusters with isolated units. In this regard, we recall that if \code{zero.policy = FALSE} and \code{Type = "SCSAR"} causes \code{SCSR_Estim(...)} to terminate with an error.
#' See package \code{spatialreg} for details on the \code{zero.policy} input.
#' @param G.set Integer vector. Sequence of clusters to be considered. Default is \code{G = c(2,3,4)}.
#' @param Phi.set Non-negative (>=0) real-valued vector. Sequence of spatial penalty parameter. Default is \code{Phi = c(0.50,1)}.
#' @param Type Character. Declares which model specification has to be estimated. Admitted strings are:
#' \itemize{
#'  \item{\code{"SCLM"} for linear regression model without spatial effects (LM);}
#'  \item{\code{"SCSAR"} for spatial autoregressive (SAR) model;}
#'  \item{\code{"SCSEM"} for linear regression model with spatial autoregressive error term or spatial Durbin model (SEM);}
#'  \item{\code{"SCSLX"} for linear regression model with spatially-lagged response variable and covariates (SLX);}
#' }
#' @param CenterVars \code{Logical} value (\code{TRUE} or \code{FALSE}) stating whether the response variable and the covariates have to be centered around the mean in the iterative algorithm to update memberships and group-wise parameters.
#' Centering is only use in the iterative procedure, while final estimates provided to the user are computed original (i.e., non-centered) variables.
#' @param ScaleVars \code{Logical} value (\code{TRUE} or \code{FALSE}) stating whether the response variable and the covariates have to be scaled with respect to their standard deviation in the iterative algorithm to update memberships and group-wise parameters.
#' Scaling is only used in the iterative procedure, while final estimates provided to the user are computed original (i.e., non-scaled) variables.
#' @param Maxitr Integer value. Maximum number of iterations for the iterative algorithm. Convergence criterion is fixed to  \eqn{\varepsilon} = 10^(-5).
#' @param RelTol Tolerance for the relative improvement in the log-likelihood (exit criterion) from iteration k to k+1. Default is \eqn{\varepsilon_{Rel}} = 10^-6
#' @param AbsTol Tolerance for the absolute improvement in the log-likelihood (exit criterion) from iteration k to k+1. Default is \eqn{\varepsilon_{Abs}} = 10^-5
#' @param Seed Integer value. Define the random number generator (RNG) state for random number generation in R.
#' Deafult is \code{seed = 123456789}.
#' @param Verbose \code{Logical} value (\code{TRUE} or \code{FALSE}). Toggle warnings and messages. If \code{verbose = TRUE} (default) the function
#' prints on the screen some messages describing the progress of the tasks. If \code{verbose = FALSE} any message about the progression is suppressed.
#'
#' @author Paolo Maranzano <>
#' @author Raffaele Mattera <>
#'
#' @return A list object containing the following outputs:
#' \itemize{
#' \item{IC: a \code{data.frame} object containing one row for each combination of the supplied vectors G.set and Phi.set and 5 columns (G,Phi,AIC,BIC,HQC).}
#' \item{OptimPars: a \code{data.frame} object with 3 rows (criteria) and 2 columns (Parameters) with the optimal combination of G and Phi for every criterion.}
#' }
#'
#'
#' @details Given the vectors G.set = c(2,3,4) and Phi.set = c(0.50,1), the function 'SCSR_InfoCrit' will compute 3x2=6 models, each at a given combination of G and Phi.
#' For computional details on the spatially-clustered models, we kindly refer to
#' Cerqueti, R., Maranzano, P. & Mattera, R. "Spatially-clustered spatial autoregressive models with application to agricultural market concentration in Europe". arXiv preprints (<doi:10.48550/arXiv.2407.15874>)
#'
#'
#' @examples
#' \donttest{
#' data(Data_RC_PM_RM_JABES2024, package="SCDA")
#' SCSAR_IC <- SCSR_InfoCrit(Formula = "Gini_SO ~ GDPPC_PPS2020 + Share_AgroEmp",
#'                           Data_sf = Data2020, listW=listW, Type="SCSAR",
#'                           Maxitr = 100, Phi.set = c(0.50,1), G.set=c(2,3))
#' }
#'
#' @export

SCSR_InfoCrit <- function(Formula, Data_sf, listW, Phi.set=c(0.50,1), G.set=c(2,3,4),
                          Type=c("SCLM","SCSAR","SCSEM","SCSLX"),
                          CenterVars = TRUE, ScaleVars = TRUE,
                          Maxitr=200,RelTol = 10^-6, AbsTol = 10^-5,
                          Verbose = TRUE, Seed = 123456789){

  ########## Intro message
  if (Verbose == TRUE) {
    TypeName <- switch(EXPR = Type,
                       "SCLM" = "linear regression model without spatial effects (LM)",
                       "SCSAR" = "Spatial autoregressive model (SAR)",
                       "SCSEM" = "linear regression model with spatial autoregressive error term (SEM)",
                       "SCSLX" = "linear regression model with spatially-lagged response variable and covariates (SLX)"
    )
    message(paste0("Computing information criteria for the spatially-clustered spatial autoregressive model started at ",round(Sys.time()),"\n"))
    message(paste0("Selected model: ",TypeName,"\n"))
  }

  ########## Setup
  Comb <- expand.grid(G = G.set,Phi = Phi.set)

  ##### Computing information criteria for each combination
  for(l in 1:dim(Comb)[1]) {
    try({
      message(paste0("Computing combination ",l," of ",dim(Comb)[1],": g=",Comb$G[l]," and phi=",Comb$Phi[l]," started at ",round(Sys.time()),"\n"))
      fit <- SCSR_Estim(Formula = Formula, Data_sf = Data_sf, listW = listW,
                        G=Comb$G[l], Phi=Comb$Phi[l],
                        CenterVars = CenterVars, ScaleVars = ScaleVars,
                        Type=Type, Maxitr=Maxitr,RelTol = RelTol, AbsTol = AbsTol,
                        Verbose = Verbose, Seed = Seed)
      ## Number of regression parameters (covariates + intercept + spatial parameter rho + residual variance sigma squared)
      pp <- summary(fit$ClusterFitModels[[1]])$parameters
      # pp <- dim(fit$Beta)[1]
      ## Total number of estimated parameters (covariates + group-wise variances)*groups
      k <- Comb$G[l]*pp
      ## Total number of observations/units
      n <- dim(fit$sBeta)[1]
      # "Knee Point Detection in BIC for Detecting the Number of Clusters" (Qinpei Zhao, Ville Hautamaki & Pasi Fränti)
      ## Clustering BIC, AIC, and HQC using MLE
      Comb$BIC[l] <- -2*fit$MLE + k*log(n)
      Comb$AIC[l] <- -2*fit$MLE + k*2
      Comb$HQC[l] <- -2*fit$MLE + k*2*log(log(n))
      # ## Clustering BIC, AIC, and HQC using PMLE
      # Comb$BIC_PMLE[l] <- -2*fit$PMLE + k*log(n)
      # Comb$AIC_PMLE[l] <- -2*fit$PMLE + k*2
      # Comb$HQC_PMLE[l] <- -2*fit$PMLE + k*2*log(log(n))
    },
    silent=TRUE)
  }

  ##### Selection
  OptimPars <- matrix(data = NA, nrow = 3, ncol = 2)
  OptimPars[1,] <- c(Comb$G[which.min(Comb$BIC)],Comb$Phi[which.min(Comb$BIC)])
  OptimPars[2,] <- c(Comb$G[which.min(Comb$AIC)],Comb$Phi[which.min(Comb$AIC)])
  OptimPars[3,] <- c(Comb$G[which.min(Comb$HQC)],Comb$Phi[which.min(Comb$HQC)])
  # OptimPars[4,] <- c(Comb$G[which.min(Comb$BIC)],Comb$Phi[which.min(Comb$BIC)])
  # OptimPars[5,] <- c(Comb$G[which.min(Comb$AIC)],Comb$Phi[which.min(Comb$AIC)])
  # OptimPars[6,] <- c(Comb$G[which.min(Comb$HQC)],Comb$Phi[which.min(Comb$HQC)])
  rownames(OptimPars) <- c("BIC","AIC","HQC")
  # "BIC - PMLE","AIC - PMLE","HQC - PMLE"
  colnames(OptimPars) <- c("G_star","Phi_star")
  OptimPars <- data.frame(OptimPars)

  ##### Exit message
  if (Verbose == TRUE) {
    message(paste0("Optimal hyperparameters (minimum BIC): g* = ",OptimPars[1,1]," and phi* = ",OptimPars[1,2],"\n"))
    message(paste0("Optimal hyperparameters (minimum AIC): g* = ",OptimPars[2,1]," and phi* = ",OptimPars[2,2],"\n"))
    message(paste0("Optimal hyperparameters (minimum HQC): g* = ",OptimPars[3,1]," and phi* = ",OptimPars[3,2],"\n"))
    # cat(paste0("Optimal hyperparameters (minimum BIC with PMLE): g* = ",OptimPars[4,1]," and phi* = ",OptimPars[4,2],"\n"))
    # cat(paste0("Optimal hyperparameters (minimum AIC with PMLE): g* = ",OptimPars[5,1]," and phi* = ",OptimPars[5,2],"\n"))
    # cat(paste0("Optimal hyperparameters (minimum HQC with PMLE): g* = ",OptimPars[6,1]," and phi* = ",OptimPars[6,2],"\n"))
    message(paste0("Computing information criteria for the spatially-clustered spatial autoregressive model ended at ",round(Sys.time()),"\n"))
  }

  ##### Results
  Result <- list(IC = Comb,OptimPars = OptimPars)
  return(Result)
}

#' Estimate spatially-clustered spatial regression models
#'
#'
#' @description Estimates spatially-clustered spatial regression (SCSR) models, such as the spatially-clustered linear regression model (SCLM),
#' the spatially-clustered spatial autoregressive model (SCSAR), the spatially-clustered spatial durbin model (SCSEM),
#' and the spatially-clustered linear regression model with spatially-lagged exogenous covariates and response variable (SCSLX).
#' Estimation is performed via cluster-wise maximum likelihood as presented in <https://arxiv.org/abs/2407.15874>.
#'
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
#' @param G Integer value. Number of clusters to be considered. When 'G=1', the pooled regression (no clusterwise) is estimated. Default is 'G = 2'.
#' @param Phi Non-negative (>=0) real value. Spatial penalty parameter. Default is 'Phi = 1'.
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
#'
#' @return A list object containing the following outputs:
#' \itemize{
#' \item{ClusterFitModels: G-dimensional list containing the estimated clustered regression models of class \code{lm} or \code{Sarlm}}
#' \item{Beta: (G x p) matrix of cluster-wise or pooled regression coefficients}
#' \item{Sig: G-dimensional vector of cluster-wise standard deviations}
#' \item{VCov: (p x p x G) array of cluster-wise variance-covariance matrices of coefficients}
#' \item{W_g: G-dimensional list containing for the g-th cluster with cardinality n_g a (n_g x n_g) spatial weighting matrix}
#' \item{listW_g: G-dimensional list containing for the g-th cluster the weights list}
#' \item{Group: (n x 1) vector of group assignment}
#' \item{sBeta: (n x p) matrix of location-wise regression coefficients}
#' \item{sSig: (n x 1) vector of location-wise standard deviations}
#' \item{MLE: Estimated maximum log-likelihood}
#' \item{Iter: The number of iteration needed to satisfy the convergence criterion and end up the clustering iterative loop}
#' }
#'
#'
#' @details The package \code{SCSR} computes the spatially-clustered spatial regression models based on the \code{spatialreg} package (see <https://cran.r-project.org/web/packages/spatialreg/index.html>).
#' SCSAR model is estimated using the function \code{lagsarlm}; SCSEM model is estimated using the function \code{errorsarlm}; SCSLX model is estimated using the function \code{lmSLX}.
#' SCLM model is estimated using the \code{lm} function from package \code{stats}.
#' Thus, estimated SCSAR, SCSEM and SCSLX models belong to class \code{Sarlm}, while estimated SCLM belongs to class \code{lm}.
#' We kindly refer to the package \code{spatialreg} for any detail regarding computational aspects (e.g., optimization).
#' Also, we refer to the package \code{spdep} for computational details on the spatial weighting matrix via \code{listw2mat(...)}, \code{nb2listw(...)} and \code{nb2mat(...)} from the \code{spdep} package.
#' For computional details on the spatially-clustered models, we kindly refer to Cerqueti, R., Maranzano, P. & Mattera, R. "Spatially-clustered spatial autoregressive models with application to agricultural market concentration in Europe". arXiv preprints (<doi:10.48550/arXiv.2407.15874>)
#'
#'
#' @examples
#' data(Data_RC_PM_RM_JABES2024, package="SCDA")
#' SCSAR <- SCSR_Estim(Formula = "Gini_SO ~ GDPPC_PPS2020 + Share_AgroEmp",
#'                     Data_sf = Data2020, G=3, listW=listW, Type="SCSAR", Phi = 0.50)
#' SCLM <- SCSR_Estim(Formula = "Gini_SO ~ GDPPC_PPS2020 + Share_AgroEmp",
#'                    Data_sf = Data2020, G=3, listW=listW, Type="SCLM", Phi = 0.50)
#'
#'
#' @export

SCSR_Estim <- function(Formula, Data_sf, listW, G = 2, Phi = 1,
                       Type=c("SCLM","SCSAR","SCSEM","SCSLX"),
                       CenterVars = FALSE, ScaleVars = FALSE,
                       Maxitr = 100, RelTol = 10^-6, AbsTol = 10^-5,
                       Verbose = TRUE, Seed = 123456789){

  ############################
  ########## Checks ##########
  ############################

  # See: https://r-pkgs.org/dependencies-in-practice.html#sec-dependencies-in-suggests-r-code

  ##### Check if package 'spatialreg' is installed
  rlang::check_installed("spatialreg",
                         reason = "Package \"spatialreg\" must be installed to estimate SCSR models.")

  ##### Check if package 'spdep' is installed
  rlang::check_installed("spdep",
                         reason = "Package \"spdep\" must be installed to estimate SCSR models.")

  ##### Define %notin%
  '%notin%' <- Negate('%in%')

  ##### Checks on the inputs
  if (Type %notin% c("SCLM","SCSAR","SCSEM","SCSLX")) {
    stop("Wrong setup for 'Type'. Use one of 'SCLM','SCSAR','SCSEM' or 'SCSLX'.",
         call. = FALSE)
  }
  if (Phi < 0) {
    stop("Wrong setup for 'Phi', which must be non-negative (i.e., Phi >= 0).",
         call. = FALSE)
  }
  if (Verbose %notin% c(FALSE,TRUE)) {
    stop("Wrong setup for 'Verbose'. Use TRUE or FALSE.", call. = FALSE)
  }
  if (!inherits(listW,what = "listw")) {
    stop("Wrong setup for 'listW'. Use an object of class 'listw'.", call. = FALSE)
  }
  Formula <- stats::as.formula(Formula)
  if (!inherits(Formula,what = "formula")) {
    stop("Wrong setup for 'Formula'. Use an object of class 'formula'.", call. = FALSE)
  }





  ###########################
  ########## Setup ##########
  ###########################

  if (Verbose == TRUE) {
    TypeName <- switch(EXPR = Type,
                       "SCLM" = "linear regression model without spatial effects (LM)",
                       "SCSAR" = "Spatial autoregressive model (SAR)",
                       "SCSEM" = "linear regression model with spatial autoregressive error term (SEM)",
                       "SCSLX" = "linear regression model with spatially-lagged response variable and covariates (SLX)"
    )
    message(paste0("Spatially-clustered spatial regression estimation started at ",round(Sys.time()),"\n"))
    message(paste0("Selected model: ",TypeName," with g=",G," groups and spatial penalty parameter phi=",Phi,"\n"))
  }

  ### Spatial information
  W <- spdep::listw2mat(listW)
  Wstyle <- listW$style
  Zpolicy <- attr(listW, "zero.policy")
  Sp_object <- methods::as(object = Data_sf, Class = "Spatial")
  Sp_coords <- sp::coordinates(Sp_object)

  ### Seed for replication
  set.seed(Seed)

  ### Standardization
  YX_names <- colnames(stats::model.frame(formula = Formula, data = Data_sf))
  Data_sf_trs <- dplyr::mutate(
    .data = Data_sf,
    dplyr::across(YX_names, ~ as.numeric(scale(x = .x, center = CenterVars, scale = ScaleVars)))
  )

  ### Extract necessary data in matrix form and corresponding names
  YX <- stats::model.frame(formula = Formula, data = Data_sf_trs)
  Y <- YX[,1]
  if (sub(pattern = "~",replacement = "",x = Formula)[3] == "1") {
    Xnames <- c("Intercept")
    # Augment data with constant column for predictions
    X <- XX <- as.matrix(rep(1,length(Y)))
  } else {
    Xnames <- c("Intercept",colnames(YX)[-1])
    X <- as.matrix(YX[,-1])
    # Augment data with constant column for predictions
    XX <- cbind(Intercept = rep(1,length(Y)),X)
  }

  ### Standardization
  Y_orig <- Y
  Y <- scale(x = Y, center = CenterVars, scale = ScaleVars)
  if (sub(pattern = "~",replacement = "",x = Formula)[3] == "1") {
    X_orig <- X
  } else {
    X_orig <- X
    X <- apply(X = X, MARGIN = 2, FUN = function(x) scale(x,center = CenterVars, scale = ScaleVars))
  }

  ### Augment data with constant column for predictions
  if (sub(pattern = "~",replacement = "",x = Formula)[3] == "1") {
    XX <- as.matrix(rep(1,length(Y)))
  } else {
    XX <- cbind(Intercept = rep(1,length(Y)),X)
  }

  ### W as a sparse matrix
  W <- methods::as(W, "sparseMatrix")
  ### Binary W matrix (for penalty term)
  if (TRUE) {
    Wbin <- (W != 0) + 0
    Wbin <- methods::as(Wbin, "sparseMatrix")
  } else {
    Wbin <- W
  }





  ################################
  ########## Dimensions ##########
  ################################

  ### n: Number of observations
  n <- length(Y)

  ### p: Number of regression coefficients (excluding residual variance)
  if (Type=="SCLM"){
    # covariates + intercept
    p <- ncol(XX)
    XXnames <- c(Xnames)
  }
  if (Type=="SCSAR"){
    # covariates  + intercept + rho
    p <- ncol(XX) + 1
    XXnames <- c("rho",Xnames)
  }
  if (Type=="SCSEM"){
    # covariates + intercept + lambda
    p <- ncol(XX) + 1
    XXnames <- c("rho",Xnames)
  }
  if (Type=="SCSLX"){
    # if (Wstand==FALSE){
    #   # covariates + intercept + lagged covariates + lambda + rho
    #   p <- 2*(dim(X)[2]+1)
    # } else {
    #   p <- (2*dim(X))+1
    # }
    p <- (2*dim(X))+1
  }





  ####################################################
  ########## Iterative clustering algorithm ##########
  ####################################################

  ##### Initial random partitioning with K-means algorithm on coordinates
  M <- 20
  WSS <- c()
  CL <- list()
  for(k in 1:M){
    CL[[k]] <- stats::kmeans(Sp_coords, G)
    WSS[k] <- CL[[k]]$tot.withinss
  }
  Ind <- CL[[which.min(WSS)]]$cluster

  ##### Objects to store
  Pen <- rep(0, G)
  Beta <- matrix(0, p, G)
  colnames(Beta) <- paste0("G=",1:G)
  rownames(Beta) <- XXnames
  Sig <- rep(1, G)
  VCov <- array(data = 0, dim = c(p,p,G))
  VCov <- vector(mode = "list", length = G)
  fit <- list()

  ##### Debug
  # k <- 2
  # Ind_check <- Ind
  # Beta_check <- Beta[,1]
  # val_check <- LL_check <- val3_check <- 0
  # Pen_check <- Ind
  # Q_check <- Ind

  ##### Iterative algorithm
  val <- 0
  for(k in 1:Maxitr){

    cval <- val

    ### Penalty term
    Ind.mat <- matrix(0, n, G)
    for(g in 1:G){
      Ind.mat[Ind==g, g] <- 1
    }
    Ind.mat <- methods::as(Ind.mat, "sparseMatrix")
    Pen <- Wbin%*%Ind.mat

    ### Model parameters (clusterwise and pooled case)
    for(g in 1:G){
      ### Control for minimum number of observations in each cluster
      if(length(Ind[Ind==g]) > p+1){
        if(Type=="SCLM"){
          suppressWarnings(
            listWg <- spdep::subset.listw(x = listW, subset = (1:length(Ind) %in% which(Ind == g)))
          )
          fit[[g]] <- stats::lm(formula = Formula, data = Data_sf_trs[which(Ind == g),])
          Beta[,g] <- as.vector( stats::coef(fit[[g]]) )
          VCov[[g]] <- stats::vcov(fit[[g]])
          resid <- Y - as.vector(XX%*%Beta[,g])
          Sig[g] <- sqrt(mean(resid[Ind==g]^2))
          Sig[g] <- max(Sig[g], 0.1)
        }
        if (Type=="SCSAR"){
          suppressWarnings(
            listWg <- spdep::subset.listw(x = listW, subset = (1:length(Ind) %in% which(Ind == g)))
          )
          fit[[g]] <- spatialreg::lagsarlm(formula = Formula, data = Data_sf_trs[which(Ind == g),],
                                           listw=listWg, zero.policy=Zpolicy)
          Beta[,g] <- as.vector(stats::coef(fit[[g]]))
          VCov[[g]] <- stats::vcov(fit[[g]])
          Sig[g] <- sqrt(fit[[g]]$s2)
        }
        if (Type=="SCSEM"){
          suppressWarnings(
            listWg <- spdep::subset.listw(x = listW, subset = (1:length(Ind) %in% which(Ind == g)))
          )
          fit[[g]] <- spatialreg::errorsarlm(formula = Formula, data = Data_sf_trs[which(Ind == g),],
                                             listw=listWg, zero.policy=Zpolicy)
          Beta[,g] <- as.vector(stats::coef(fit[[g]]))
          VCov[[g]] <- stats::vcov(fit[[g]])
          Sig[g] <- sqrt(fit[[g]]$s2)
        }
        if (Type=="SCSLX"){
          # if (Wstand==FALSE){
          #   WI <- spdep::lag.listw(spdep::mat2listw(W), XX[,1], zero.policy = zero.policy)
          #   WX <- spatialreg::create_WX(x = XX[,-1], listw = spdep::mat2listw(W), zero.policy=zero.policy)
          #   XXX <- cbind(XX,WI,WX)
          # } else {
          #   WX <- spatialreg::create_WX(x = XX[,-1], listw = spdep::mat2listw(W), zero.policy=zero.policy)
          #   XXX <- cbind(XX,WX)
          # }
          suppressWarnings(
            listWg <- spdep::subset.listw(x = listW, subset = (1:length(Ind) %in% which(Ind == g)))
          )
          fit[[g]] <- spatialreg::lmSLX(formula = Formula, data = Data_sf_trs[which(Ind == g),],
                                        listw=listWg, zero.policy=Zpolicy)
          Beta[,g] <- as.vector(stats::coef(fit[[g]]))
          VCov[[g]] <- stats::vcov(fit[[g]])
          WX <- spatialreg::create_WX(x = XX[,-1], listw = listW, zero.policy=Zpolicy)
          XXX <- cbind(XX,WX)
          resid <-  Y - as.vector(XXX%*%Beta[,g])
          Sig[g] <- max(sqrt(mean(resid[Ind==g]^2)),0.10)
        }
        if (k == 1) {
          mval <- sum(unlist(lapply(fit,function(x){as.numeric(x$LL)})))
        }
      }
      # else print error/warning "N_k < P, estimation not possible!"
    }

    ##### Step 2: Update group membership
    # New assignment is performed by maximizing the element-wise sum of the likelihood at the old groupings and the penalty term
    if(Type=="SCLM"){
      Mu <- XX%*%Beta      # (n,G)-matrix
      ESig <- t(matrix(rep(Sig,n), G, n))    # (n,G)-matrix
      LL <- stats::dnorm(Y, Mu, ESig, log=T)
      Q <- LL + Phi*Pen
    }
    if (Type=="SCSAR"){
      LL <- matrix(NA, n, G)
      for (g in 1:G) {
        LL[,g] <- as.vector(sarlogLik_i(Y=Y,X=XX,W=W,beta=Beta[-1,g],rho=Beta[1,g],sigma2=Sig[g]^2))
      }
      Q <- LL + Phi*Pen
    }
    if (Type=="SCSEM"){
      LL <- matrix(NA, n, G)
      for (g in 1:G) {
        LL[,g] <- as.vector(semlogLik_i(Y=Y,X=XX,W=W,beta=Beta[-1,g],lamb=Beta[1,g],sigma2=Sig[g]^2))
      }
      Q <- LL + Phi*Pen
    }
    if (Type=="SCSLX"){
      Mu <- XXX%*%Beta     # (n,G)-matrix
      ESig <- t(matrix(rep(Sig,n), G, n))    # (n,G)-matrix
      LL <- stats::dnorm(Y, Mu, ESig, log=T)
      Q <- LL + Phi*Pen
    }

    Ind <- apply(Q, 1, which.max)

    ##### Value of objective function: sum of the G maximum penalized log-likelihood
    # if (Type == "SCLM") {
    #   val <- sum(unlist(lapply(fit,function(x){as.numeric(stats::logLik(fit[[x]]))})))
    # } else {
    #   val <- sum(unlist(lapply(fit,function(x){as.numeric(x$LL)})))
    # }
    LL_g <- c()
    Q_g <- c()
    for(g in 1:G){
      LL_g[g] <- sum(LL[which(Ind == g), g])
      Q_g[g] <- sum(Q[Ind == g, g])
    }
    val <- sum(Q_g)


    ##### Check the exit conditions
    # Max log-likehood
    if (k > 1) {
      mval <- max(mval,val)
    }
    # Relative improvement
    RelImp <- abs((cval-mval)/mval)
    # Absolute improvement
    AbsImp <- abs(mval-cval)
    # Message
    if (Verbose == TRUE) {
      message(paste0("* Iteration ",k,": log-likelihood improvements -- Absolute:",
                     round(AbsImp,6)," ; Relative: ",round(RelImp,4),"\n"))
    }
    # AbsImp < AbsTol
    if( RelImp < RelTol | AbsImp < AbsTol ) {
      break
    }

    # ##### Debug
    # Ind_check <- cbind(Ind,Ind_check)
    # Beta_check <- cbind(Beta[,1],Beta_check)
    # val_check <- cbind(val,val_check)
    # LL_check <- cbind(LL[,1],LL_check)
    # Pen_check <- cbind(Pen[,1],Pen_check)
    # Q_check <- cbind(Q[,1],Q_check)
    # cat(as.matrix(c(k,RelImp,AbsImp,cval,val,mval,cval-val)))
    # cat("\n")
  }
  cat("\n")


  ##### Check convergence
  if (k < Maxitr) {
    Convergence <- "Convergence reached"
  }
  if (k == Maxitr) {
    Convergence <- "Maximum number of iterations reached"
  }



  ########## Final regression results (re-estimation)
  W_g <- listW_g <- vector(mode = "list", length = G)
  if (Type=="SCLM"){
    for (g in 1:G) {
      suppressWarnings(
        listW_g[[g]] <- spdep::subset.listw(x = listW, subset = (1:length(Ind) %in% which(Ind == g)))
      )
      W_g[[g]] <- spdep::listw2mat(listW_g[[g]])
      fit[[g]] <- stats::lm(formula = Formula, data = Data_sf[which(Ind == g),])
      Beta[,g] <- as.vector(stats::coef(fit[[g]]))
      VCov[[g]] <- stats::vcov(fit[[g]])
      resid <- Y - as.vector(XX%*%Beta[,g])
      Sig[g] <- sqrt(mean(resid[Ind==g]^2))
      Sig[g] <- max(Sig[g], 0.1)
    }
  }
  if (Type=="SCSAR"){
    for (g in 1:G) {
      suppressWarnings(
        listW_g[[g]] <- spdep::subset.listw(x = listW, subset = (1:length(Ind) %in% which(Ind == g)))
      )
      W_g[[g]] <- spdep::listw2mat(listW_g[[g]])
      fit[[g]] <- spatialreg::lagsarlm(formula = Formula, data = Data_sf[which(Ind == g),],
                                       listw=listW_g[[g]], zero.policy=Zpolicy)
      Beta[,g] <- as.vector(stats::coef(fit[[g]]))
      VCov[[g]] <- stats::vcov(fit[[g]])
      Sig[g] <- max(sqrt(fit[[g]]$s2), 0.1)
    }
  }
  if (Type=="SCSEM"){
    for (g in 1:G) {
      suppressWarnings(
        listW_g[[g]] <- spdep::subset.listw(x = listW, subset = (1:length(Ind) %in% which(Ind == g)))
      )
      W_g[[g]] <- spdep::listw2mat(listW_g[[g]])
      fit[[g]] <- spatialreg::errorsarlm(formula = Formula, data = Data_sf[which(Ind == g),],
                                         listw=listW_g[[g]], zero.policy=Zpolicy)
      Beta[,g] <- as.vector(stats::coef(fit[[g]]))
      VCov[[g]] <- stats::vcov(fit[[g]])
      Sig[g] <- max(sqrt(fit[[g]]$s2), 0.1)
    }
  }
  if (Type=="SCSLX"){
    for (g in 1:G) {
      suppressWarnings(
        listW_g[[g]] <- spdep::subset.listw(x = listW, subset = (1:length(Ind) %in% which(Ind == g)))
      )
      W_g[[g]] <- spdep::listw2mat(listW_g[[g]])
      fit[[g]] <- spatialreg::lmSLX(formula = Formula, data = Data_sf[which(Ind == g),],
                                    listw=listW_g[[g]], zero.policy=Zpolicy)
      Beta[,g] <- as.vector(stats::coef(fit[[g]]))
      VCov[[g]] <- stats::vcov(fit[[g]])
      WX <- spatialreg::create_WX(x = XX[,-1], listw = listW, zero.policy=Zpolicy)
      XXX <- cbind(XX,WX)
      resid <-  Y - as.vector(XXX%*%Beta[,g])
      Sig[g] <- max(sqrt(mean(resid[Ind==g]^2)),0.10)
    }
  }


  ##### Final value of objective function: sum of the G maximum log-likelihood
  MLE <- sum(unlist(lapply(fit, logLik)))
  MLE2 <- sum(LL_g)
  PMLE <- sum(Q_g)
  # if (Type == "SCLM") {
  #   ML <- sum(unlist(lapply(fit,function(x){as.numeric(stats::logLik(fit[[x]]))})))
  # } else {
  #   MLE <- sum(unlist(lapply(fit,function(x) {as.numeric(logLik(fit[[x]]))})))
  # }


  ##### Spatially-varying parameters
  if (sub(pattern = "~",replacement = "",x = Formula)[3] == "1" & Type == "SCLM") {
    # Location-wise regression coefficients
    sBeta <- as.matrix(Beta[,Ind])
    colnames(sBeta) <- "Intercept"
    # Location-wise error variance
    sSig <- Sig[Ind]
  } else {
    # Location-wise regression coefficients
    sBeta <- t(Beta[,Ind])
    # Location-wise error variance
    sSig <- Sig[Ind]
  }


  ##### Change names to coefficients and var-cov matrix
  for (g in 1:G) {
    names(fit[[g]]$coefficients) <- Xnames
  }


  ##### Exit message
  if (Verbose == TRUE) {
    message(paste0("Spatially-clustered spatial regression estimation ended at ",round(Sys.time()),"\n"))
  }


  ########## Results
  result <- list(ClusterFitModels = fit, Beta = Beta, Sig = Sig, VCov = VCov,
                 Group = Ind, sBeta = sBeta, sSig = sSig,
                 MLE = MLE, MLE2 = MLE2, PMLE = PMLE, Iter = k,
                 W_g = W_g, listW_g = listW_g, Convergence = Convergence)
  return(result)

}

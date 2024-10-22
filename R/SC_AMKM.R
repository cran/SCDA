#' Spatial Clustering for sf data
#'
#'
#' @description Perform spatial clustering using K-means and AMKM (Adjacent Matrix K-Means) algorithms on sf data.
#'
#'
#' @param Data_sf A \code{data.frame} object of class \code{sf} with n rows (each one corresponding to a location) and a user-defined number of columns.
#' It must include the \code{geometry} feature for spatial modelling and representation.
#' Typically, \code{sf} \code{data.frame} are built using the \code{st_as_sf(...)} command from the \code{sf} package (see its documentation for details).
#' @param IndexCol \code{Integer} value. Number of the dataset ID column. If there isn't an ID column \code{IndexCol=0}.
#' @param Method \code{Character}. Must be one of: 'AMKM' or 'K-means'. If \code{method='AMKM'}, the Adjacent Matrix K-Means clustering is performed. If \code{method='K-means'}, K-means clustering is performed.
#' @param Distance \code{Character}. The distance measure to be used to compute the dissimilarity matrix. This must be one of: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski". By default, \code{distance='euclidean'}.
#' @param MinNc \code{Integer} value. Minimal number of clusters, between 1 and (number of objects - 1). Default is \code{MinNc=2}.
#' @param MaxNc \code{Integer} value. Maximal number of clusters, between 2 and (number of objects - 1), greater or equal to MinNc. Default is \code{MaxNc=10}.
#' @param Metric \code{Character}. The validation index to be calculated for the selection of the optimal clustering partition. This should be one of : "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw", "all" (all indices except GAP, Gamma, Gplus and Tau), "alllong" (all indices with Gap, Gamma, Gplus and Tau included).
#' Default is \code{Metric='silhouette'}.
#' @param RidDim \code{Character}.The dimensionality reduction method. This should be one of : 'pca' or 'laplacian'. if 'RidDim='pca'', a principal component analysis is performed. if 'RidDim='laplacian'' the laplacian matrix dimensionality reduction method is performed . Default is \code{RidDim='pca'}.
#' @param CenterVars \code{Logical} value (\code{TRUE} or \code{FALSE}) stating whether the features have to be centered around the mean. Default is \code{TRUE}.
#' @param ScaleVars \code{Logical} value (\code{TRUE} or \code{FALSE}) stating whether the features have to be scaled with respect to their standard deviation. Default is \code{TRUE}.
#' @param MakePlot \code{Logical} value (\code{TRUE} or \code{FALSE}) stating whether the plot must be displayed. Default is \code{TRUE}.
#' @param ExplainedVariance \code{numeric}. cumulate percentage of the variance explained by the eigenvalues of the dimesionality reduction method. Must be between 0 and 1. Default is \code{ExplainedVariance=0.9}.
#' @param KeepCoord \code{Logical} value (\code{TRUE} or \code{FALSE}) stating whether the coordinate must be taken into account in K-means algorithm. Available only when 'method='K-means''. Default is \code{TRUE}.
#' @param CRS \code{Integer} value. Coordinate reference system. something suitable as input to st_crs.command from the \code{sf} package (see its documentation for details). Default is \code{CRS=4326}
#' @param Verbose \code{Logical} value (\code{TRUE} or \code{FALSE}). Toggle warnings and messages. If \code{verbose = TRUE} (default) the function
#' prints on the screen some messages describing the progress of the tasks. If \code{verbose = FALSE} any message about the progression is suppressed. Default is \code{TRUE}.
#' @param Seed \code{Integer} value. Define the random number generator (RNG) state for random number generation in R.
#' Deafult is \code{seed = 123456789}.
#'
#' @return A list object containing the following outputs:
#' \itemize{
#' \item{df: n row dataframe with the following columns : ID, Longitude, Latitude and Cluster (the optimal partition)}
#' \item{plot: Display cluster partition in a map.}
#' }
#'
#' @details AMKM calculations is done decomposing the input dataset in two subset. The first one contains the features while the second one contains the coordinates (longitude and latitude).
#' A dissimilarity matrix is calculated on both subset using the parameter distance for the feature and the Great Circle distance for coordinates.
#' Then an adjacent matrix (n x n) is computed on every dissimilarity matrix using gaussian kernel.
#' To reduce the dimensionality of the adjacent matrix a dimentionality reduction method is necessary (see RidDim param. for more)
#' K-means is applied with no modification at its original algorithm.
#'
#'@author Camilla Lionetti <lionetticamilla511@gmail.com>, Francesco Caccia <francesco.caccia2000@gmail.com>
#'
#' @examples
#' library(sp)
#' library(sf)
#' data("meuse")
#' dati<-meuse
#' dati<-subset(dati,select=sapply(dati,is.numeric))
#' dati<-st_as_sf(dati, coords = c("x", "y"),crs =28992)
#' SC <- SC_AMKM(Data_sf=dati,IndexCol=0, Method="AMKM",MinNc = 5,MaxNc = 5 ,CRS=28992)
#'
#' @export

SC_AMKM <- function(Data_sf, IndexCol, Method, Distance="euclidean",MinNc=2, MaxNc=10, Metric="silhouette", RidDim="pca", CenterVars=T, ScaleVars=T, MakePlot=T, ExplainedVariance=0.9, KeepCoord=T, Seed=123456789, Verbose=T, CRS=4326){



  ############################
  ########## Checks ##########
  ############################

  # See: https://r-pkgs.org/dependencies-in-practice.html#sec-dependencies-in-suggests-r-code


  rlang::check_installed("sf",
                         reason = "Package \"sf\" must be installed to handle the dataset.")

  rlang::check_installed("NbClust",
                         reason = "Package \"NbClust\" must be installed to perform clustering.")

  if (MakePlot==T){
    rlang::check_installed("ggplot2",
                           reason = "Package \"ggplot2\" must be installed to make the plot.")
    rlang::check_installed("ggspatial",
                           reason = "Package \"ggspatial\" must be installed to make the plot.")

  }

  ##### Define %notin%
  '%notin%' <- Negate('%in%')

  ##### Checks on the inputs

  if (Verbose %notin% c(FALSE,TRUE)) {
    stop("Wrong setup for 'Verbose'. Use TRUE or FALSE.", call. = FALSE)
  }

  if (Method %notin% c("AMKM","K-means")) {
    stop("Wrong setup for 'Method'. Use AMKM or K-means.", call. = FALSE)
  }

  if (RidDim %notin% c("pca","laplacian")) {
    stop("Wrong setup for 'RidDim'. Use pca or laplacian.", call. = FALSE)
  }

  if (MakePlot %notin% c(FALSE,TRUE)) {
    stop("Wrong setup for 'MakePlot'. Use TRUE or FALSE.", call. = FALSE)
  }

  if (CenterVars %notin% c(FALSE,TRUE)) {
    stop("Wrong setup for 'CenterVars'. Use TRUE or FALSE.", call. = FALSE)
  }

  if (ScaleVars %notin% c(FALSE,TRUE)) {
    stop("Wrong setup for 'ScaleVars'. Use TRUE or FALSE.", call. = FALSE)
  }

  if ((ExplainedVariance <0) |(ExplainedVariance >=1)) {
    stop("Wrong setup for 'ExplainedVariance'. Must be between 0 and 1", call. = FALSE)
  }

  if (KeepCoord %notin% c(FALSE,TRUE)) {
    stop("Wrong setup for 'KeepCoord'. Use TRUE or FALSE.", call. = FALSE)
  }

  #creation of two subset:
  set.seed(Seed)
  data <- stats::na.omit(Data_sf) #remove NA's

  if (Method=="AMKM"){
    coord<-as.data.frame(st_coordinates(data))
    coord = sf::st_as_sf(coord, coords = c("X", "Y"),crs = CRS)
    data<-sf::st_drop_geometry(data)

    if (IndexCol>0){
      index<-data[,IndexCol]
      data<-data[,-IndexCol]
      if ((CenterVars==T) & (ScaleVars==T)){
        data<-scale(data, center=T, scale=T)
      }
      if ((CenterVars==F) & (ScaleVars==T)){
        data<-scale(data, center=F, scale=T)
      }
      if ((CenterVars==T) & (ScaleVars==F)){
        data<-scale(data, center=T, scale=F)
      }


      #AM for coord:
      D<-as.matrix(sf::st_distance(coord))
      attr(D,"class") <- NULL
      A <- exp(-D/(2*max(D)))#gaussian kernel

      #AM for data:
      D2<-as.matrix(stats::dist(data, method = Distance))
      attr(D2,"class") <- NULL
      A2 <- exp(-D2/(2*max(D2)))#gaussian kernel
      if (RidDim=="pca"){

        #pca coord
        if (Verbose==T){
          message(paste0(("PCA in  progress...")))
        }

        pca_D<-stats::prcomp(A, center = TRUE, scale. = TRUE)
        var_explained <- pca_D$sdev^2 / sum(pca_D$sdev^2)
        cum_var_explained <- cumsum(var_explained)
        num_components <- which(cum_var_explained >= ExplainedVariance)[1]
        spcomp <- scale(as.matrix(pca_D$x[,1:num_components]), center=T, scale=T)

        #pca data
        pca_D2<-stats::prcomp(A2, center = TRUE, scale. = TRUE)
        var_explained <- pca_D2$sdev^2 / sum(pca_D2$sdev^2)
        cum_var_explained <- cumsum(var_explained)
        num_components <- which(cum_var_explained >= ExplainedVariance)[1]
        spcomp2 <- scale(as.matrix(pca_D2$x[,1:num_components]), center=T, scale=T)
      }

      if (RidDim=="laplacian"){
        if (Verbose==T){
          message(paste0(("Laplacian dim reduction in progress...")))
        }

        B<-diag(rowSums(A))
        B_inv_sqrt<- diag(1 / sqrt(diag(B)))
        L <- B_inv_sqrt%*% (B - A) %*% B_inv_sqrt
        eigenvector<-as.data.frame(eigen(L)$vectors)
        eigenvalues<-(eigen(L)$values)

        var_explained <- eigenvalues/sum(eigenvalues)
        cum_var_explained <- cumsum(var_explained)
        num_components <- which(cum_var_explained >= ExplainedVariance)[1]
        spcomp <- scale(as.matrix(eigenvector[,1:num_components]), center=T, scale=T)

        B2<-diag(rowSums(A2))
        B2_inv_sqrt<- diag(1 / sqrt(diag(B2)))
        L2 <- B2_inv_sqrt%*% (B2 - A2) %*% B2_inv_sqrt
        eigenvector2<-as.data.frame(eigen(L2)$vectors)
        eigenvalues2<-(eigen(L2)$values)

        var_explained <- eigenvalues2/sum(eigenvalues2)
        cum_var_explained <- cumsum(var_explained)
        num_components <- which(cum_var_explained >= ExplainedVariance)[1]
        spcomp2 <- scale(as.matrix(eigenvector[,1:num_components]), center=T, scale=T)
      }

      #AMKM
      if (Verbose==T){
        message(paste0(("AMKM in progress...")))
      }

      X<-as.matrix(cbind(spcomp, spcomp2))
      D_full <- as.matrix(stats::dist(X, method = Distance))
      km_sp<-NbClust::NbClust(data=X, method="kmeans", distance=NULL, diss=D_full, min.nc = MinNc, max.nc= MaxNc, index=Metric)
      out<-as.data.frame(cbind(index,st_coordinates(coord), km_sp$Best.partition))
      colnames(out)<-c("ID","longitude", "latitude", "cluster")
      out$cluster<-as.factor(out$cluster)
      #plot:
      if(MakePlot==T){
        if(Verbose==T){
          message(paste0("Making plot..."))
        }

        g0=ggplot2::ggplot()
        g1=ggspatial::annotation_map_tile(type="osm",zoomin = 0)
        g2=ggplot2::geom_sf(data = sf::st_as_sf(out[,-1], coords = c("longitude", "latitude"),crs = CRS), mapping=ggplot2::aes(color = .data$cluster), size = 2) # aggiunge i punti
        g3=ggspatial::annotation_north_arrow(which_north = "true", location="tl")
        g4=ggspatial::annotation_scale(location="br") # punti cardinali e scala
        output_plot<-g0+g1+g2+g3+g4+
          ggplot2::geom_sf()+
          ggplot2::labs(title = "AMKM")+
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4))
        l<-list(df=out, plot=output_plot)
        return(l)
      }
      return(df=out)
    } else{
      if ((CenterVars==T) & (ScaleVars==T)){
        data<-scale(data, center=T, scale=T)
      }
      if ((CenterVars==F) & (ScaleVars==T)){
        data<-scale(data, center=F, scale=T)
      }
      if ((CenterVars==T) & (ScaleVars==F)){
        data<-scale(data, center=T, scale=F)
      }

      #AM for coord:
      D<-as.matrix(sf::st_distance(coord))
      attr(D,"class") <- NULL
      A <- exp(-D/(2*max(D)))#gaussian kernel

      #AM for data:
      D2<-as.matrix(stats::dist(data, method = Distance))
      attr(D2,"class") <- NULL
      A2 <- exp(-D2/(2*max(D2)))#gaussian kernel
      if (RidDim=="pca"){
        if(Verbose==T){
          message(paste0("PCA in progress..."))
        }

        #pca coord
        pca_D<-stats::prcomp(A, center = TRUE, scale. = TRUE)
        var_explained <- pca_D$sdev^2 / sum(pca_D$sdev^2)
        cum_var_explained <- cumsum(var_explained)
        num_components <- which(cum_var_explained >= ExplainedVariance)[1]
        spcomp <- scale(as.matrix(pca_D$x[,1:num_components]), center=T, scale=T)

        #pca data
        pca_D2<-stats::prcomp(A2, center = TRUE, scale. = TRUE)
        var_explained <- pca_D2$sdev^2 / sum(pca_D2$sdev^2)
        cum_var_explained <- cumsum(var_explained)
        num_components <- which(cum_var_explained >= ExplainedVariance)[1]
        spcomp2 <- scale(as.matrix(pca_D2$x[,1:num_components]), center=T, scale=T)
      }

      if (RidDim=="laplacian"){
        if(Verbose==T){
          message(paste0("Laplacian dim rdeuction in progress..."))
        }

        B<-diag(rowSums(A))
        B_inv_sqrt<- diag(1 / sqrt(diag(B)))
        L <- B_inv_sqrt%*% (B - A) %*% B_inv_sqrt
        eigenvector<-as.data.frame(eigen(L)$vectors)
        eigenvalues<-(eigen(L)$values)

        var_explained <- eigenvalues/sum(eigenvalues)
        cum_var_explained <- cumsum(var_explained)
        num_components <- which(cum_var_explained >= ExplainedVariance)[1]
        spcomp <- scale(as.matrix(eigenvector[,1:num_components]), center=T, scale=T)

        B2<-diag(rowSums(A2))
        B2_inv_sqrt<- diag(1 / sqrt(diag(B2)))
        L2 <- B2_inv_sqrt%*% (B2 - A2) %*% B2_inv_sqrt
        eigenvector2<-as.data.frame(eigen(L2)$vectors)
        eigenvalues2<-(eigen(L2)$values)

        var_explained <- eigenvalues2/sum(eigenvalues2)
        cum_var_explained <- cumsum(var_explained)
        num_components <- which(cum_var_explained >= ExplainedVariance)[1]
        spcomp2 <- scale(as.matrix(eigenvector[,1:num_components]), center=T, scale=T)
      }
      #AMKM###
      if(Verbose==T){
        message(paste0("AMKM in progress..."))
      }

      X<-as.matrix(cbind(spcomp, spcomp2))
      D_full <- as.matrix(stats::dist(X, method = Distance))
      km_sp<-NbClust::NbClust(data=X, method="kmeans", distance=NULL, diss=D_full, min.nc = MinNc, max.nc= MaxNc, index=Metric)
      out<-as.data.frame(cbind(st_coordinates(coord), km_sp$Best.partition))
      colnames(out)<-c("longitude", "latitude", "cluster")
      out$cluster<-as.factor(out$cluster)
      #plot:
      if(MakePlot==T){
        if(Verbose==T){
          message(paste0(("Making plot...")))
        }

        g0=ggplot2::ggplot()
        g1=ggspatial::annotation_map_tile(type="osm",zoomin = 0)
        g2=ggplot2::geom_sf(data = st_as_sf(out, coords = c("longitude", "latitude"),crs = CRS), mapping=ggplot2::aes(color = .data$cluster), size = 2) # aggiunge i punti
        g3=ggspatial::annotation_north_arrow(which_north = "true", location="tl")
        g4=ggspatial::annotation_scale(location="br") # punti cardinali e scala
        output_plot<-g0+g1+g2+g3+g4+
          ggplot2::geom_sf()+
          ggplot2::labs(title = "AMKM")+
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4))
        l<-list(df=out, plot=output_plot)
        return(l)
      }
      return(df=out)

    }
  }

  #METHOD= K_MEANS#####

  if(Method=="K-means"){
    coord<-as.data.frame(sf::st_coordinates(data))
    data<-sf::st_drop_geometry(data)


    ifelse(KeepCoord==T, data<-cbind(data, coord$X, coord$Y), data<-data)



    if (IndexCol>0){
      index<-data[,IndexCol]
      data<-data[,-IndexCol]
      if ((CenterVars==T) & (ScaleVars==T)){
        data<-scale(data, center=T, scale=T)
      }
      if ((CenterVars==F) & (ScaleVars==T)){
        data<-scale(data, center=F, scale=T)
      }
      if ((CenterVars==T) & (ScaleVars==F)){
        data<-scale(data, center=T, scale=F)
      }

      #K-means
      if(Verbose==T){
        message(paste0(("K-Means in progress...")))
      }

      km_sp<-NbClust::NbClust(data=data, method="kmeans", distance=Distance, min.nc = MinNc, max.nc= MaxNc, index=Metric)
      out<-as.data.frame(cbind(index,coord, km_sp$Best.partition))
      colnames(out) <- c("ID","longitude", "latitude", "cluster")
      out$cluster <- as.factor(out$cluster)
      #plot:
      if(MakePlot==T){
        if(Verbose==T){
          message(paste0("Making plot..."))
        }

        g0=ggplot2::ggplot()
        g1=ggspatial::annotation_map_tile(type="osm",zoomin = 0)
        g2=ggplot2::geom_sf(data = st_as_sf(out[,-1], coords = c("longitude", "latitude"),crs = CRS), mapping=ggplot2::aes(color = .data$cluster), size = 2) # aggiunge i punti
        g3=ggspatial::annotation_north_arrow(which_north = "true", location="tl")
        g4=ggspatial::annotation_scale(location="br") # punti cardinali e scala
        output_plot<-g0+g1+g2+g3+g4+
          ggplot2::geom_sf()+
          ggplot2::labs(title = "KM")+
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4))
        l<-list(df=out, plot=output_plot)
        return(l)
      }
      return(df=out)
    } else{
      if ((CenterVars==T) & (ScaleVars==T)){
        data<-scale(data, center=T, scale=T)
      }
      if ((CenterVars==F) & (ScaleVars==T)){
        data<-scale(data, center=F, scale=T)
      }
      if ((CenterVars==T) & (ScaleVars==F)){
        data<-scale(data, center=T, scale=F)
      }

      #K-Means###
      if(Verbose==T){
        message(paste0("K-Means in progress..."))
      }

      km_sp<-NbClust::NbClust(data=data, method="kmeans", distance=Distance, min.nc = MinNc, max.nc= MaxNc, index=Metric)
      out<-as.data.frame(cbind(coord, km_sp$Best.partition))
      colnames(out)<-c("longitude", "latitude", "cluster")
      out$cluster<-as.factor(out$cluster)
      #plot:
      if(MakePlot==T){
        if(Verbose==T){
          message(paste0("Making plot..."))
        }

        g0=ggplot2::ggplot()
        g1=ggspatial::annotation_map_tile(type="osm",zoomin = 0)
        g2=ggplot2::geom_sf(data = st_as_sf(out, coords = c("longitude", "latitude"),crs = CRS), mapping=ggplot2::aes(color = .data$cluster), size = 2) # aggiunge i punti
        g3=ggspatial::annotation_north_arrow(which_north = "true", location="tl")
        g4=ggspatial::annotation_scale(location="br") # punti cardinali e scala
        output_plot<-g0+g1+g2+g3+g4+
          ggplot2::geom_sf()+
          ggplot2::labs(title = "KM")+
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4))
        l<-list(df=out, plot=output_plot)
        return(l)
      }
      return(df=out)

    }
  }


}






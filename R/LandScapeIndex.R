require(dplyr)
require(ggplot2)
require(data.table)
require(rgdal)
require(rgeos)

getData <- function(){
  Sh <- readOGR('DATA/Bauges_Landscape/forestStandsSimulatedPolygons.shp')
  Da <- read.csv('DATA/Bauges_Landscape/forestStands_diameterDistribution.txt', sep=';')

  ListID <- unique(Sh@data$STAND_ID)
  DF <- data.frame(ID=rep(NA, length(ListID)), XCenter=NA, YCenter=NA, Area=NA)
  j <- 0
  for (i in ListID){
     j <- j + 1
     tt <- Sh[Sh@data$STAND_ID==i, ]
     Area <- tt@polygons[[1]]@area
     J <- as.vector(gCentroid(tt)@coords)
     DF[j,] <- c(i, J[1], J[2], Area)
  }

  DF <- left_join(DF, Da, by=c('ID'='STAND_ID'))
}

#' Compute Shannon Index
#'
#' This function compute the Shannon index for a gridded landscape 
#' @param DFQuad, data.frame of gridded landscape from Quad function
#' @return alpha, beta and gamma diversity based on Shannon Index
#' @export
ComputeSh <- function(DFQuad){
  plog <- function(x){
    x[is.na(x)] <- 0
    x <- x/sum(x)
    -sum(x*log(x), na.rm=TRUE)
  }
  iClass  <- which(substr(names(DFQuad),1,7)=='N_Class')
  SH <- apply(DFQuad[, ..iClass], 1, plog)
  DFtot <- DFQuad[, ..iClass] * matrix(rep(DFQuad$Area, length(iClass)), nrow=dim(DFQuad)[1])
  DFtot <- apply(DFtot, 2, sum, na.rm=TRUE) / sum(DFtot, na.rm=TRUE)
  ### ?????????????????????
  alpha <- sum(SH*DFQuad$Area) / sum(DFQuad$Area)
  ### ???????????????????
  gamma <- -sum(DFtot * log(DFtot), na.rm=TRUE) 
  beta <- gamma/alpha
  return(c(alpha, beta, gamma))
}

#' Compute Shannon Index for a given scale
#'
#' This function compute the Shannon index for a landscape and a resolution
#' @param DF, data.frame of landscape data
#' @param Res, numeric resolution (in km) of the grid
#' @return alpha, beta and gamma diversity based on Shannon Index
#' @export
ComputeShScale <- function(DF, Res=1){
  SH <- data.frame(N=1:9, Res=Res, Ngrid=NA, alpha=NA, beta=NA, gamma=NA)
  for (N in 1:9){
    DFi <- Quad(DF, Res=Res, N=N)
    SH[N, ] <- c(N, Res, dim(DFi)[1], ComputeSh(DFi))
  }
  return(SH)
}

#' Compute a gridded landscape
#'
#' This function compute the gridded landscape
#' @param DF, data.frame of landscape data
#' @param Res, numeric resolution (in km) of the grid
#' @param N, integer from 1 to 9; here to include different grid configurations
#' @return data.frame of gridded landscape
#' @export
Quad <- function(DF, Res=1, N=9){
    xrange <- range(DF$XCenter) * 1e-3
    yrange <- range(DF$YCenter) * 1e-3
    dgrid <- expand.grid(dx=seq(Res/10, Res/2, length.out=3), dy=seq(Res/10, Res/2, length.out=3))
    dgrid <- dgrid[N, ]
    xgrid <- seq((xrange[1]-dgrid$dx), (xrange[2]+Res/2), by=Res)
    ygrid <- seq((yrange[1]-dgrid$dy), (yrange[2]+Res/2), by=Res)
#    xgrid <- seq(floor(xrange[1]-Res/2), ceiling(xrange[2]+Res/2), by=Res)
#    ygrid <- seq(floor(yrange[1]-Res/2), ceiling(yrange[2]+Res/2), by=Res)
    DF <- mutate(DF, xind=findInterval(XCenter * 1e-3, xgrid),
            yind=findInterval(YCenter * 1e-3, ygrid), Iind=paste(xind, yind, sep='/'))
    DF <- as.data.table(DF)
    pArea <- DF[, pArea:=Area/sum(Area), by=list(xind,yind)][, pArea]
    iClass  <- which(substr(names(DF),1,7)=='N_Class')
    DFClass <- DF[, ..iClass]
    DFClass <- DFClass * matrix(rep(pArea, dim(DFClass)[2]), nrow=dim(DF)[1])
    DFClass <- mutate(DFClass, Iind=paste(DF$xind, DF$yind, sep='/'))
    DF2 <- DFClass[, lapply(.SD, sum, na.rm=TRUE), by=Iind] 
    DF3 <- DF[,list(
	    N=.N,
	    DG_1=sum(DG_1*Area*NHA_1, na.rm=TRUE)/sum(Area*NHA_1, na.rm=TRUE),
	    DG_2=sum(DG_2*Area*NHA_2, na.rm=TRUE)/sum(Area*NHA_2, na.rm=TRUE),
	    DG_tot=sum(Dgtot*Area*(NHA_1+NHA_2), na.rm=TRUE)/sum(Area*(NHA_1+NHA_2), na.rm=TRUE),
	    Gsp1=sum(Gsp1*Area, na.rm=TRUE)/sum(Area, na.rm=TRUE),
	    Gsp2=sum(Gsp2*Area, na.rm=TRUE)/sum(Area, na.rm=TRUE),
	    G=sum(G*Area, na.rm=TRUE)/sum(Area, na.rm=TRUE),
	    NHA_1=sum(Area*NHA_1, na.rm=TRUE)/sum(Area, na.rm=TRUE),
	    NHA_2=sum(Area*NHA_2, na.rm=TRUE)/sum(Area, na.rm=TRUE),
	    Area=sum(Area, na.rm=TRUE),
	    XCenter=mean(XCenter), YCenter=mean(YCenter),
	    xgrid=xgrid[mean(xind)]+Res/2, ygrid=ygrid[mean(yind)]+Res/2), by=Iind]
    DF <- merge(DF2, DF3, by="Iind")
    return(DF)
}


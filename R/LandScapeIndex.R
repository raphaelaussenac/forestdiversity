require(dplyr)
require(ggplot2)
require(data.table)
require(rgdal)
require(raster)
require(rgeos)
require(tidyr)

#' Retrieve Landscape dataset
#'
#' This function gather landscape data into on data.table 
#' @return DF, data.table
#' @export
getdataBaugesInit <- function(){
  LandscapeRaw <- read.asciigrid('DATA/Bauges_Landscape/cellID.asc')
  if (exists('DATA/Bauges_Landscape/trees75.Rds')){
    DF <- readRDS('DATA/Bauges_Landscape/trees75.Rds')
  }else{
    DF <- read.csv('DATA/Bauges_Landscape/trees75.csv')
    DF <- data.table::as.data.table(DF)
  }
  DF <- DF[(!is.na(dbh))]
  DF <- dplyr::mutate(DF, Area=25*25*1e-4)
  Landscape <- data.table::as.data.table(coordinates(LandscapeRaw))
  Landscape <- dplyr::mutate(Landscape, cellID=LandscapeRaw@data[[1]])
  Landscape <- Landscape[(cellID %in% DFstand$cellID)]
  names(Landscape)[(1:2)] <- c('XCenter', 'YCenter')
  DF <- merge(DF, Landscape, by='cellID')
  return(DF)
}

#' Compute class Landscape
#'
#' This function takes a data.table of landscape and return in a class form 
#' @param DF, data.table, the landscape with dbh for each tree
#' @return DFclass, data.table, the equivalent of DF in class form
#' @export
LandscapeDBHtoClass <- function(DF){
  IntClassDBH <- seq(5, 255, by=10)
  DF <- dplyr::mutate(DF, Cat=findInterval(dbh, IntClassDBH), Area=25*25*1e-4)
  DFstand <- DF[, .(NHA=sum(n/Area), BA=sum((n*pi*(dbh/200)^2)/Area),
      Dg=sqrt(sum(n*dbh^2)/sum(n)), Area=mean(Area), XCenter=XCenter[1], YCenter=YCenter[1]), by="cellID"]
  DFclass <- DF[, .(n=sum(n)), by=list(cellID, Cat)]
  DFclass <- pivot_wider(DFclass, values_from=n, names_from=Cat,
         names_prefix='Nclass_', values_fill=list(n=0))
  DFclass <- merge(DFstand, DFclass, by='cellID')
  return(DFclass)
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
  iClass  <- which(substr(names(DFQuad),1,7)=='Nclass_')
  SH <- apply(DFQuad[, ..iClass], 1, plog)
  DFtot <- DFQuad[, ..iClass] * matrix(rep(DFQuad$Area, length(iClass)), nrow=dim(DFQuad)[1])
  DFtot <- apply(DFtot, 2, sum, na.rm=TRUE) / sum(DFtot, na.rm=TRUE)
  alpha <- exp(sum(DFQuad$Area*log(SH)*is.finite(log(SH)), na.rm=TRUE)/sum(DFQuad$Area)) 
  gamma <- -sum(DFtot * log(DFtot), na.rm=TRUE) 
  beta <- gamma/alpha
  CVG <- sum(DFQuad$Area * (DFQuad$G - mean(DFQuad$G))^2 / sum(DFQuad$Area))
  CVDg <- sum(DFQuad$Area * (DFQuad$Dg - mean(DFQuad$Dg))^2 / sum(DFQuad$Area))
  return(c(alpha, beta, gamma, CVG, CVDg))
}

#' Compute Shannon Index for a given scale
#'
#' This function compute the Shannon index for a landscape and a resolution
#' @param DF, data.frame of landscape data
#' @param Res, numeric resolution (in km) of the grid
#' @return alpha, beta and gamma diversity based on Shannon Index
#' @export
ComputeShScale <- function(DF, Res=1){
  SH <- data.frame(N=1:9, Res=Res, Ngrid=NA, alpha=NA, beta=NA, gamma=NA, CVG=NA,CVDg=NA)
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
    DF <- dplyr::mutate(DF, xind=findInterval(XCenter * 1e-3, xgrid),
            yind=findInterval(YCenter * 1e-3, ygrid),
	    Iind=paste(xind, yind, sep='/'))
    pArea <- DF[, pArea:=Area/sum(Area), by=list(xind,yind)][, pArea]
    iClass  <- which(substr(names(DF),1,7)=='Nclass_')
    DFClass <- DF[, ..iClass]
    DFClass <- DFClass * matrix(rep(pArea, dim(DFClass)[2]), nrow=dim(DF)[1])
    DFClass <- mutate(DFClass, Iind=paste(DF$xind, DF$yind, sep='/'))
    DF2 <- DFClass[, lapply(.SD, sum, na.rm=TRUE), by=Iind] 
    DF3 <- DF[,list(
	    Dg=sum(Dg*Area*NHA, na.rm=TRUE)/sum(Area*NHA, na.rm=TRUE),
	    BA=sum(BA*Area, na.rm=TRUE)/sum(Area, na.rm=TRUE),
	    Area=sum(Area, na.rm=TRUE),
	    XCenter=mean(XCenter), YCenter=mean(YCenter),
	    xgrid=xgrid[mean(xind)]+Res/2, ygrid=ygrid[mean(yind)]+Res/2), by=Iind]
    DF <- merge(DF2, DF3, by="Iind")
    return(DF)
}

Plot_Grid <- function(DF, Res=1, N=9){
	require(viridis)
  DFQuad <- Quad(DF, Res=Res, N=9)
  plog <- function(x){
        x[is.na(x)] <- 0
        x <- x/sum(x)
       -sum(x*log(x), na.rm=TRUE)
  }
  iClass  <- which(substr(names(DFQuad),1,7)=='Nclass_')
  SH <- apply(DFQuad[, ..iClass], 1, plog)
  DFQuad <- mutate(DFQuad, Sh=SH)
  pl <- ggplot(DFQuad, aes(xmin=xgrid-Res/2, xmax=xgrid+Res/2,
	ymin=ygrid-Res/2, ymax=ygrid+Res/2, fill=Sh)) +
        scale_fill_viridis(begin=0, end=2) + geom_rect() +
	theme(text=element_text(size=24)) +
       	xlab('') + ylab('')
  png(file=paste0('Grid',Res,'.png'), width=1800, height=1200,res=120);print(pl);dev.off()
}

Plot_Shannon_Scale <- function(DF, Res=c(0.05, 0.1, 1, 2, 5)){
    SHscale <- do.call(rbind, lapply(Res, ComputeShScale, DF=DF))
    SHm <- group_by(SHscale, Res) %>% summarise(alphaM=mean(alpha),
	alphaS=sd(alpha), betaM=mean(beta), betaS=sd(beta), 
	CVGM=mean(CVG), CVDgM=mean(CVDg),
	CVGS=sd(CVG), CVDgS=sd(CVDg),
	gamma=mean(gamma)) %>% ungroup()
    p1 <- ggplot(SHm, aes(x=Res)) + geom_pointrange(aes(y=alphaM,
	ymin=alphaM-2*alphaS, ymax=alphaM+2*alphaS)) + 
        geom_line(aes(y=alphaM)) + theme(text=element_text(size=24)) +
        xlab('Grain') + ylab('Alpha')
    p2 <- ggplot(SHm, aes(x=Res)) + geom_pointrange(aes(y=betaM,
	ymin=betaM-2*betaS, ymax=betaM+2*betaS), col='red') +
        geom_line(aes(y=betaM), col='red') + theme(text=element_text(size=24)) +
        xlab('Grain') + ylab('Beta')
    pl <- multiplot(p1, p2, cols=2)
    p3 <- ggplot(SHm, aes(x=Res)) + geom_pointrange(aes(y=CVGM,
	ymin=CVGM-2*CVGS, ymax=CVGM+2*CVGS), col='red') +
        geom_line(aes(y=CVGM), col='red') + theme(text=element_text(size=24)) +
        xlab('Grain') + ylab('Basal area variance')
    p4 <- ggplot(SHm, aes(x=Res)) + geom_pointrange(aes(y=CVDgM,
	ymin=CVDgM-2*CVDgS, ymax=CVDgM+2*CVDgS), col='red') +
        geom_line(aes(y=CVDgM), col='red') + theme(text=element_text(size=24)) +
        xlab('Grain') + ylab('Quadratic diameter variance')
    return(pl)
}






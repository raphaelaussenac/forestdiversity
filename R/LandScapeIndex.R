#' Retrieve Landscape dataset
#'
#' This function gather landscape data into on data.table 
#' @return DF, data.table
#' @export
getdataBaugesInit <- function(){
  LandscapeRaw <- sp::read.asciigrid('DATA/Bauges_Landscape/cellID.asc')
  if (file.exists('DATA/Bauges_Landscape/trees75.Rds')){
    DF <- readRDS('DATA/Bauges_Landscape/trees75.Rds')
  }else{
    DF <- read.csv('DATA/Bauges_Landscape/trees75.csv')
    DF <- data.table::as.data.table(DF)
    DF <- DF[(!is.na(dbh))]
    DF <- saveRDS(file='DATA/Bauges_Landscape/trees75.Rds', DF)
  }
  DF <- dplyr::mutate(DF, Area=25*25*1e-4)
  Landscape <- data.table::as.data.table(sp::coordinates(LandscapeRaw))
  Landscape <- dplyr::mutate(Landscape, cellID=LandscapeRaw@data[[1]])
  names(Landscape)[(1:2)] <- c('XCenter', 'YCenter')
  DF <- merge(DF, Landscape, by='cellID', all.x=TRUE)
  return(DF)
}

#' Compute class Landscape
#'
#' This function takes a data.table of landscape and return in a class form 
#' @param DF, data.table, the landscape with dbh for each tree
#' @param Nvar, string, the variable used to build class, whehter 'species' or 'D_cm'
#' @param ClassInter num, Class size (same units as the variable chosen)
#' @param ClassIni num, min dbh of the first class
#' @return DFclass, data.table, the equivalent of DF in class form
#' @export
LandscapeDBHtoClass <- function(DF, Nvar='D_cm', ClassInter=10, ClassIni=7.5){
    if (Nvar=='D_cm'){
        ClassInterval <- seq(ClassIni, max(DF$dbh) + 2* ClassInter, by=ClassInter)
        ClassName <- ClassInterval + mean(diff(ClassInterval)) / 2
        DF <- dplyr::mutate(DF, Class=ClassName[findInterval(dbh, ClassInterval)],
   	    Area= 25*25*1e-4)
    }else if(Nvar=='species'){
        DF <- dplyr::mutate(DF, Class=sp, Area=25*25*1e-4)
    }
    DFstand <- DF[, .(NHA=sum(n/Area), BA=sum((n*pi*(dbh/200)^2)/Area),
        Dg=sqrt(sum(n*dbh^2)/sum(n)), Area=mean(Area), XCenter=XCenter[1], YCenter=YCenter[1]), by="cellID"]
    DFclass <- DF[, .(N=sum(n), BA=sum(n*pi*(dbh/200)^2/Area)), by=list(cellID, Class)]
    DFNclass <- tidyr::pivot_wider(dplyr::select(DFclass, cellID, Class, N),
	 values_from=N, names_from=Class, names_prefix='Nclass_', values_fill=list(N=0))
    DFBAclass <- tidyr::pivot_wider(dplyr::select(DFclass, cellID, Class, BA),
	 values_from=BA, names_from=Class, names_prefix='BAclass_', values_fill=list(BA=0))
    DFclass <- merge(DFNclass, DFBAclass, by='cellID')
    DFclass <- merge(DFstand, DFclass, by='cellID')
    return(DFclass)
}

#' Compute Shannon Index
#'
#' This function compute the Shannon index for a gridded landscape 
#' @param DFgrid, Grid object from GridLandcsape function
#' @return diversity indices
#' @export
ComputeSh <- function(DFgrid){
  plog <- function(x){
        x[is.na(x)] <- 0
        x <- x/sum(x)
       -sum(x*log(x), na.rm=TRUE)
  }
  Nclass <- function(x){
        x[is.na(x)] <- 0
        sum(x>0)
  }
## Density proportion
  iClass  <- which(substr(names(DFgrid),1,7)=='Nclass_')
  SH <- apply(DFgrid[, ..iClass], 1, plog)
  DFtot <- DFgrid[, ..iClass] * matrix(rep(DFgrid$Area, length(iClass)), nrow=dim(DFgrid)[1])
  DFtot <- apply(DFtot, 2, sum, na.rm=TRUE) / sum(DFtot, na.rm=TRUE)
  alpha <- exp(sum(DFgrid$Area*log(SH)*is.finite(log(SH)), na.rm=TRUE)/sum(DFgrid$Area)) 
  gamma <- -sum(DFtot * log(DFtot), na.rm=TRUE) 
  beta <- gamma/alpha
## BA proportion
  iClass  <- which(substr(names(DFgrid),1,8)=='BAclass_')
  SH <- apply(DFgrid[, ..iClass], 1, plog)
  DFtot <- DFgrid[, ..iClass] * matrix(rep(DFgrid$Area, length(iClass)), nrow=dim(DFgrid)[1])
  DFtot <- apply(DFtot, 2, sum, na.rm=TRUE) / sum(DFtot, na.rm=TRUE)
  alphaBA <- exp(sum(DFgrid$Area*log(SH)*is.finite(log(SH)), na.rm=TRUE)/sum(DFgrid$Area)) 
  gammaBA <- -sum(DFtot * log(DFtot), na.rm=TRUE) 
  betaBA <- gammaBA/alphaBA
## Nb of class
  alphaNcl <- sum(DFgrid$Area * (apply(DFgrid[, ..iClass], 1, Nclass))) / sum(DFgrid$Area)
  gammaNcl <- Nclass(apply(DFgrid[, ..iClass], 2, sum))
  betaNcl <- gammaNcl / alphaNcl
  CVBA <- sum(DFgrid$Area * (DFgrid$BA - mean(DFgrid$BA))^2 / sum(DFgrid$Area))
  CVDg <- sum(DFgrid$Area * (DFgrid$Dg - mean(DFgrid$Dg))^2 / sum(DFgrid$Area))
  return(data.frame(alphaNcl=alphaNcl, betaNcl=betaNcl, gammaNcl=gammaNcl, alpha=alpha, beta=beta, gamma=gamma, 
      alphaBA=alphaBA, betaBA=betaBA, gammaBA=gammaBA, CVBA=CVBA, CVDg=CVDg, Res=DFgrid$Res[1]))
}

#' Compute Shannon Index for a given scale
#'
#' This function compute the Shannon index for a landscape and a resolution
#' @param DF, data.frame of landscape data
#' @param Res, numeric resolution (in km) of the grid
#' @return Diversity indices
#' @export
ComputeShScale <- function(DF, Res=1){
  SH <- NULL
    for (N in 1:9){
    DFi <- GridLandscape(DF, Res=Res, N=N)
    SH <- rbind(SH, ComputeSh(DFi))
  }
  return(SH)
}

#' Compute a gridded landscape
#'
#' This function compute the gridded landscape
#' @param DFclass, data.frame of landscape data
#' @param Res, numeric resolution (in km) of the grid
#' @param N, integer from 1 to 9; here to include different grid configurations
#' @return data.frame of gridded landscape
#' @export
GridLandscape <- function(DFclass, Res=1, N=9){
    xrange <- range(DFclass$XCenter) * 1e-3
    yrange <- range(DFclass$YCenter) * 1e-3
    dgrid <- expand.grid(dx=seq(Res/10, Res/2, length.out=3), dy=seq(Res/10, Res/2, length.out=3))
    dgrid <- dgrid[N, ]
    xgrid <- seq((xrange[1]-dgrid$dx), (xrange[2]+Res/2), by=Res)
    ygrid <- seq((yrange[1]-dgrid$dy), (yrange[2]+Res/2), by=Res)
    DFclass <- dplyr::mutate(DFclass, xind=findInterval(XCenter * 1e-3, xgrid),
            yind=findInterval(YCenter * 1e-3, ygrid),
	    Iind=paste(xind, yind, sep='/'))
    DFclass <- data.table::as.data.table(DFclass)
    pArea <- DFclass[, pA:=Area/sum(Area), by=list(xind,yind)]
    pArea <- pArea[, pA]
    iClass  <- which(substr(names(DFclass),1,7)=='Nclass_')
    DFiNclass <- DFclass[, ..iClass]
    DFiNclass <- DFiNclass * matrix(rep(pArea, dim(DFiNclass)[2]), nrow=dim(DFclass)[1])
    DFiNclass <- dplyr::mutate(DFiNclass, Iind=paste(DFclass$xind, DFclass$yind, sep='/'))
    DF2 <- DFiNclass[, lapply(.SD, sum, na.rm=TRUE), by=Iind] 
    iClass  <- which(substr(names(DFclass),1,8)=='BAclass_')
    DFiBAclass <- DFclass[, ..iClass]
    DFiBAclass <- DFiBAclass * matrix(rep(pArea, dim(DFiBAclass)[2]), nrow=dim(DFclass)[1])
    DFiBAclass <- dplyr::mutate(DFiBAclass, Iind=paste(DFclass$xind, DFclass$yind, sep='/'))
    DF2b <- DFiBAclass[, lapply(.SD, sum, na.rm=TRUE), by=Iind] 
    DF2 <- merge(DF2, DF2b ,by='Iind')
    DF3 <- DFclass[, list(
	    Dg=sum(Dg*Area*NHA, na.rm=TRUE)/sum(Area*NHA, na.rm=TRUE),
	    BA=sum(BA*Area, na.rm=TRUE)/sum(Area, na.rm=TRUE),
	    Ncells=.N,
	    Area=sum(Area, na.rm=TRUE),
	    XCenter=mean(XCenter), YCenter=mean(YCenter),
	    xgrid=xgrid[mean(xind)]+Res/2, ygrid=ygrid[mean(yind)]+Res/2), by=Iind]
    DF <- merge(DF2, DF3, by="Iind")
    DF <- dplyr::mutate(DF, Res=Res, N=N)
    class(DF) <- append('Grid', class(DF))
    return(DF)
}

plot.Grid <- function(DFgrid){
    plog <- function(x){
        x[is.na(x)] <- 0
        x <- x/sum(x)
       -sum(x*log(x), na.rm=TRUE)
    }
    iClass  <- which(substr(names(DFgrid),1,7)=='Nclass_')
    SH <- apply(DFgrid[, ..iClass], 1, plog)
    DFgrid <- dplyr::mutate(DFgrid, Sh=SH)
    Res <- DFgrid$Res[1]
    pl <- ggplot2::ggplot(DFgrid, ggplot2::aes(xmin=xgrid-Res/2, xmax=xgrid+Res/2,
	ymin=ygrid-Res/2, ymax=ygrid+Res/2, fill=Sh)) +
        viridis::scale_fill_viridis() + #begin=0, end=2) +
       	ggplot2::geom_rect() +
	ggplot2::theme(text=ggplot2::element_text(size=24)) +
       	ggplot2::xlab('') + ggplot2::ylab('')
#    png(file=paste0('Grid',Res,'.png'), width=1800, height=1200,res=120)
    print(pl)#;dev.off()
}

plotAllGrid <- function(DFclass, Res=1){
    plog <- function(x){
        x[is.na(x)] <- 0
        x <- x/sum(x)
       -sum(x*log(x), na.rm=TRUE)
    }
    DFt <- NULL
    for (i in 1:9){
        DFgrid <-  GridLandscape(DFclass, Res=Res, N=i)
        iClass  <- which(substr(names(DFgrid),1,7)=='Nclass_')
        SH <- apply(DFgrid[, ..iClass], 1, plog)
        DFgrid <- dplyr::mutate(DFgrid, Sh=SH)
	DFt <- rbind(DFt, DFgrid)
    }
    DFt <- dplyr::filter(DFt, Area > (100 * Res^2 *.5))
    pl <- ggplot2::ggplot(DFt, ggplot2::aes(xmin=xgrid-Res/2, xmax=xgrid+Res/2,
	ymin=ygrid-Res/2, ymax=ygrid+Res/2, fill=Sh)) +
        viridis::scale_fill_viridis() + #begin=0, end=2) +
       	ggplot2::geom_rect() + ggplot2::facet_wrap(~N) +
	ggplot2::theme(text=ggplot2::element_text(size=24)) +
       	ggplot2::xlab('') + ggplot2::ylab('')
    print(pl)
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






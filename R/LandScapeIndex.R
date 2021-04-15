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
	    NHA=sum(NHA*Area, na.rm=TRUE)/sum(Area, na.rm=TRUE),
	    Ncells=.N,
	    Area=sum(Area, na.rm=TRUE),
	    XCenter=mean(XCenter), YCenter=mean(YCenter),
	    xgrid=xgrid[mean(xind)]+Res/2, ygrid=ygrid[mean(yind)]+Res/2), by=Iind]
    ### Other varaibles (biodiversity)
    iOther <- which(!(substr(names(DFclass),1,7) %in% substr(c(names(DF2), names(DF3), 'cellID', 'xind', 'yind','pA'),1,7)))
    DF4 <- DFclass[, lapply(.SD, meanW, W=Area), by=Iind, .SDcols=iOther]
    DF <- merge(DF2, DF3, by="Iind")
    DF <- merge(DF, DF4, by="Iind")
    DF <- dplyr::mutate(DF, Res=Res, N=N)
    class(DF) <- append('Grid', class(DF))
    return(DF)
}

#' Plot a gridded landscape
#'
#' This function compute the gridded landscape
#' @param DFclass, data.frame of gridded landscape data
#' @return plot of Basal area
#' @export
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
	ymin=ygrid-Res/2, ymax=ygrid+Res/2, fill=BA)) +
        viridis::scale_fill_viridis() + #begin=0, end=2) +
       	ggplot2::geom_rect() +
	ggplot2::theme(text=ggplot2::element_text(size=24)) +
       	ggplot2::xlab('') + ggplot2::ylab('')
#    png(file=paste0('Grid',Res,'.png'), width=1800, height=1200,res=120)
    print(pl)#;dev.off()
}
#' Compute heterogeneity Index for a gridded landscape
#'
#' This function compute the heterogeneity indices for a gridded landscape 
#' @param DFgrid, Grid object from GridLandcsape function
#' @return diversity indices
#' @export
ComputeHeterogeneity <- function(DFgrid){
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
  return(data.frame(alphaNcl=alphaNcl, betaNcl=betaNcl, gammaNcl=gammaNcl, alpha=alpha, beta=beta, gamma=gamma, 
      alphaBA=alphaBA, betaBA=betaBA, gammaBA=gammaBA, CVBA=CVBA, Res=DFgrid$Res[1]))
}

#' Compute biodiversity Index for a gridded landscape
#'
#' This function compute the heterogeneity indices for a gridded landscape 
#' @param DFgrid, Grid object from GridLandcsape function
#' @return diversity indices
#' @export
ComputeBiodiversity <- function(DFGrid){
    iC <- which(!(substr(names(DFGrid),1,7) %in% c('Nclass_', 'BAclass', 'Iind',
	 'Ncells', 'XCenter', 'YCenter', 'xgrid', 'ygrid', 'Res', 'N', 'Area', 'BA', 'Dg', 'NHA')))
    DFg <- DFGrid[, ..iC]
    DFMean <- apply(DFg, 2, meanW, W=DFGrid$Area)
    DF10 <- apply(DFg, 2, reldist::wtd.quantile, na.rm=TRUE, q=0.1, weight=DFGrid$Area)
    DF90 <- apply(DFg, 2, reldist::wtd.quantile, na.rm=TRUE, q=0.9, weight=DFGrid$Area)
    names(DF10) <- paste0(names(DF10), 'P10')
    names(DF90) <- paste0(names(DF90), 'P90')
    OUTThresh <- NULL
    if ('VDW' %in% names(DFGrid)){
        VDW20 <- sum(DFGrid$Area[DFGrid$VDW>=20]) * 1e-2
        VDW60 <- sum(DFGrid$Area[DFGrid$VDW>=60]) * 1e-2
	temp <- c(VDW20, VDW60);names(temp) <- c('VDW20', 'VDW60')
        OUTThresh <- c(OUTThresh, temp)
    }
    if ('LSDTN' %in% names(DFGrid)){
        LSDTN1 <- sum(DFGrid$Area[DFGrid$LSDTN>=1]) / sum(DFGrid$Area)
        LSDTN3 <- sum(DFGrid$Area[DFGrid$LSDTN>=3]) / sum(DFGrid$Area)
	temp <- c(LSDTN1, LSDTN3);names(temp) <- c('LSDTN1', 'LSDTN3')
        OUTThresh <- c(OUTThresh, temp)
    }
    if ('LLDTN' %in% names(DFGrid)){
        LLDTN2 <- sum(DFGrid$Area[DFGrid$VLLTN>=2]) / sum(DFGrid$Area)
        LLDTN6 <- sum(DFGrid$Area[DFGrid$VLLTN>=6]) / sum(DFGrid$Area)
	temp <- c(LLDTN2, LLDTN6);names(temp) <- c('LLDTN2', 'LLDTN6')
        OUTThresh <- c(OUTThresh, temp)
    }
    if ('Cover' %in% names(DFGrid)){
	Cover50 <- sum(DFGrid$Area[DFGrid$Cover>=50]) / sum(DFGrid$Area)
        temp <- Cover50; names(temp) <- 'Cover50'
	OUTThresh <- c(OUTThresh, temp)
    }
    if (is.null(OUTThresh)){
        OUT <- data.frame(t(DFMean), t(DF10), t(DF90))
    }else{
        OUT <- data.frame(t(DFMean), t(DF10), t(DF90),  t(OUTThresh))
    }
    return(OUT)
}

#' Compute heterogeneity and biodiversity metrics for a distribution landscape and a given resolution
#'
#' This function compute the heterogeneity and biodiversity indices for a landscape and a resolution
#' @param DF, data.frame of landscape data
#' @param Res, numeric resolution (in km) of the grid
#' @param Out, str, if equal Mean return mean over grid configuration, otherwise return for each configuration
#' @return Diversity indices
#' @export
ComputeHeterogeneityScale <- function(DF, Res=1, Out='Mean'){
    SH <- NULL
        for (N in 1:9){
        DFi <- GridLandscape(DF, Res=Res, N=N)
        SH <- rbind(SH, cbind(ComputeHeterogeneity(DFi), ComputeBiodiversity(DFi)))
    }
    if (Out=='Mean'){
        return(apply(SH, 2, mean))
    }else{
        return(dplyr::mutate(SH, N=1:9))
    }
}

#' Compute Heterogeneity and biodiversitymetrics for a different scales
#'
#' This function compute the indices index for a landscape and a resolution
#' @param DF, data.frame of distribution landscape data
#' @param Res, numeric vector, resolution (in km) of the grid
#' @return Plot
#' @export
ComputeHeterogeneityMultiScale <- function(DF, Res=c(0.05, 0.1, 1, 2, 5), PLOT=FALSE){
    HetMultiscale <- do.call(rbind, lapply(Res, ComputeHeterogeneityScale, DF=DF, Out='All'))
    HetMultiscale <- data.table::as.data.table(HetMultiscale)
    HetMultiscaleMean <- HetMultiscale[, lapply(.SD, mean, na.rm=TRUE), by=Res]
    names(HetMultiscaleMean) <- paste0('Mean', names(HetMultiscaleMean))
    HetMultiscaleSD <- HetMultiscale[, lapply(.SD, sd, na.rm=TRUE), by=Res]
    names(HetMultiscaleSD) <- paste0('SD', names(HetMultiscaleSD))
    OUT <- cbind(HetMultiscaleMean, HetMultiscaleSD)
    p1 <- ggplot2::ggplot(OUT, ggplot2::aes(x=Res)) + ggplot2::geom_pointrange(ggplot2::aes(y=Meanalpha,
	ymin=Meanalpha-2*SDalpha, ymax=Meanalpha+2*SDalpha)) + 
        ggplot2::geom_line(ggplot2::aes(y=Meanalpha)) + ggplot2::theme(text=ggplot2::element_text(size=24)) +
        ggplot2::xlab('Grain') + ggplot2::ylab('Alpha')
    p2 <- ggplot2::ggplot(OUT, ggplot2::aes(x=Res)) + ggplot2::geom_pointrange(ggplot2::aes(y=Meanbeta,
	ymin=Meanbeta-2*SDbeta, ymax=Meanbeta+2*SDbeta), col='red') +
        ggplot2::geom_line(ggplot2::aes(y=Meanbeta), col='red') + ggplot2::theme(text=ggplot2::element_text(size=24)) +
        ggplot2::xlab('Grain') + ggplot2::ylab('Beta')
    p3 <- ggplot2::ggplot(OUT, ggplot2::aes(x=Res)) + ggplot2::geom_pointrange(ggplot2::aes(y=MeanCVBA,
	ymin=MeanCVBA-2*SDCVBA, ymax=MeanCVBA+2*SDCVBA), col='red') +
        ggplot2::geom_line(ggplot2::aes(y=MeanCVBA), col='red') + ggplot2::theme(text=ggplot2::element_text(size=24)) +
        ggplot2::xlab('Grain') + ggplot2::ylab('Basal area variance')
    p4 <- ggplot2::ggplot(OUT, ggplot2::aes(x=Res)) + ggplot2::geom_pointrange(ggplot2::aes(y=MeanCoverP90,
	ymin=MeanCoverP90-2*SDCoverP90, ymax=MeanCoverP90+2*SDCoverP90), col='red') +
        ggplot2::geom_line(ggplot2::aes(y=MeanCoverP90), col='red') + ggplot2::theme(text=ggplot2::element_text(size=24)) +
        ggplot2::xlab('Grain') + ggplot2::ylab('9th Decile of canopy cover')

    pl <- multiplot(p1, p3, p2, p4, cols=2)
    if (PLOT==TRUE){print(pl)}
    return(OUT)
}

####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

meanW <- function(X, W){
    iNA <- which(!is.na(X) & !is.na(W))
    W <- X[iNA]
    X <- X[iNA]
    return(sum(X*W)/sum(W))
}


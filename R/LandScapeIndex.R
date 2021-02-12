require(dplyr)
require(ggplot2)
require(data.table)
require(rgdal)
require(raster)
require(rgeos)
require(tidyr)

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
    DF <- mutate(DF, xind=findInterval(XCenter * 1e-3, xgrid),
            yind=findInterval(YCenter * 1e-3, ygrid),
	    Iind=paste(xind, yind, sep='/'))
    DF <- as.data.table(DF)
    pArea <- DF[, pArea:=Area/sum(Area), by=list(xind,yind)][, pArea]
    iClass  <- which(substr(names(DF),1,7)=='Nclass_')
    DFClass <- DF[, ..iClass]
    DFClass <- DFClass * matrix(rep(pArea, dim(DFClass)[2]), nrow=dim(DF)[1])
    DFClass <- mutate(DFClass, Iind=paste(DF$xind, DF$yind, sep='/'))
    DF2 <- DFClass[, lapply(.SD, sum, na.rm=TRUE), by=Iind] 
    DF3 <- DF[,list(
	    Dg=sum(Dg*Area*NHA, na.rm=TRUE)/sum(Area*NHA, na.rm=TRUE),
	    G=sum(G*Area, na.rm=TRUE)/sum(Area, na.rm=TRUE),
	    Area=sum(Area, na.rm=TRUE),
	    XCenter=mean(XCenter), YCenter=mean(YCenter),
	    xgrid=xgrid[mean(xind)]+Res/2, ygrid=ygrid[mean(yind)]+Res/2), by=Iind]
    DF <- merge(DF2, DF3, by="Iind")
    return(DF)
}

###########Temporary functions
##############INI
getdatadbh <- function(){
  LandscapeRaw <- read.asciigrid('DATA/Bauges_Landscape/BaugesInit/cellID.asc')
  if (exists('DATA/Bauges_Landscape/BaugesInit/trees75.Rds')){
      DF <- readRDS('DATA/Bauges_Landscape/BaugesInit/trees75.Rds') %>%
	  filter(!is.na(dbh))
  }else{
      DF <- read.csv('DATA/Bauges_Landscape/BaugesInit/trees75.csv') %>%
	  filter(!is.na(dbh))
  }
  DF <- mutate(DF, Area=25*25*1e-4)
  Landscape <- as.data.frame(coordinates(LandscapeRaw)) %>%
	  mutate(cellID=LandscapeRaw@data[[1]]) %>%
	  filter((cellID %in% unique(DF$cellID))) %>%
	  mutate(XCenter=s1, YCenter=s2) %>% dplyr::select(-s1, -s2)
  DF <- left_join(DF, Landscape, by='cellID')
  return(DF)
}

getdataClass <- function(){
  LandscapeRaw <- read.asciigrid('DATA/Bauges_Landscape/BaugesInit/cellID.asc')
  if (exists('DATA/Bauges_Landscape/BaugesInit/trees75.Rds')){
      DF <- readRDS('DATA/Bauges_Landscape/BaugesInit/trees75.Rds') %>%
	  filter(!is.na(dbh))
  }else{
      DF <- read.csv('DATA/Bauges_Landscape/BaugesInit/trees75.csv') %>%
	  filter(!is.na(dbh))
  }
  IntClassDBH <- seq(7.5, 150, by=10)
  DF <- mutate(DF, Cat=findInterval(dbh, IntClassDBH)) %>%
	  mutate(Area=25*25*1e-4)
  DFstand <- group_by(DF, cellID) %>% summarise(NHA=sum(n / Area),
	G=sum((n*pi*(dbh/200)^2)/Area), Dg=sqrt(sum(n*dbh^2)/sum(n)),
       	Area=mean(Area)) %>% ungroup()
  DF <- group_by(DF, cellID, Cat) %>% summarise(n=sum(n)) %>% ungroup()
  DF <- pivot_wider(DF, values_from=n, names_from=Cat,
         names_prefix='Nclass_', values_fill=list(n=0))
  Landscape <- as.data.frame(coordinates(LandscapeRaw)) %>%
	  mutate(cellID=LandscapeRaw@data[[1]]) %>%
	  filter((cellID %in% unique(DF$cellID))) %>%
	  mutate(XCenter=s1, YCenter=s2) %>% dplyr::select(-s1, -s2) %>%
	  left_join(DFstand, by='cellID')
  DF <- left_join(DF, Landscape, by='cellID')
  return(DF)
}

###########################
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
###############################
################################
Analyse_Gini <- function(DF){
    Sp1 <- 'Abies alba'
    Sp2 <- 'Picea abies'
    Sp3 <- 'Fagus sylvatica'
    ListSp <- c(paste0('Mono', c(Sp1, Sp2, Sp3)),
		paste0('Bi', c(Sp1, Sp2, Sp3), c(Sp3, Sp1, Sp2)),
		paste0('Bi', c(Sp3, Sp1, Sp2), c(Sp1, Sp2, Sp3)),
		paste0('Tri', Sp1, Sp2, Sp3), paste0('Tri', Sp1, Sp3, Sp2),
		paste0('Tri', Sp2, Sp1, Sp3), paste0('Tri', Sp2, Sp3, Sp1),
		paste0('Tri', Sp3, Sp1, Sp2), paste0('Tri', Sp3, Sp2, Sp1))
    A <- readRDS('outGini3.Rds') %>% filter(Area==1)
#### Find MonoSp
    A <- mutate(A, Comp=NA)
    A$Comp[A$PrSp1>0.75] <- paste0('Mono', A$Sp1[A$PrSp1>0.75])
#### Find Bisp
    ind <- which(A$PrSp1<=0.75 & A$PrSp1>0.25 & A$PrSp2>=0.25 & (A$PrSp1+A$PrSp2)>0.75)
    A$Comp[ind] <- paste0('Bi', A$Sp1[ind], A$Sp2[ind])
### Find Trisp
    ind <- which(A$PrSp1>=0.25 & A$PrSp2>=0.25 & A$PrSp3>=0.25 & (A$PrSp1+A$PrSp2)<=0.75 & (A$PrSp1+A$PrSp2+A$PrSp3)>0.75)
    A$Comp[ind] <- paste0('Tri', A$Sp1[ind], A$Sp2[ind], A$Sp3[ind])
####
    AFilt <- filter(A, Comp %in% ListSp)
    AFilt$Comp[substr(AFilt$Comp, 1, 3)=='Tri'] <- 'Tri'
    AFilt$Comp[AFilt$Comp==paste0('Bi', Sp2, Sp1)] <- paste0('Bi', Sp1, Sp2)
    AFilt$Comp[AFilt$Comp==paste0('Bi', Sp3, Sp1)] <- paste0('Bi', Sp1, Sp3)
    AFilt$Comp[AFilt$Comp==paste0('Bi', Sp3, Sp2)] <- paste0('Bi', Sp2, Sp3)
### Clusters
    AFilt <- mutate(AFilt, cl=NA)
    B <- AFilt
    B[,3:5]=apply(B[,3:5],2,scale)
    ListComp <- names(table(B$Comp))
    for (K in ListComp){
         Bk <- filter(B, Comp==K)
         Kcl <- kmeans(Bk[, c(3,4,5)],centers=10);print(Kcl$between/Kcl$totss)
         AFilt$cl[AFilt$Comp==K] <- Kcl$cluster
    }
    CompoCl <- group_by(AFilt, Comp, cl) %>% summarise(Dgm=median(Dg), 
	NHAm=median(NHA), Ginim=median(Gini),# BAm=median(BA),
	#BA1=quantile(BA,.1), BA9=quantile(BA,.9),
       	Dg1=quantile(Dg,.1), Dg9=quantile(Dg,.9),
	Gini1=quantile(Gini,.1), Gini9=quantile(Gini,.9),
       	NHA1=quantile(NHA, .1), NHA9=quantile(NHA, .9), Ncells=length(unique(cellID))) %>% ungroup()
    pl <- ggplot(CompoCl,aes(x=Dgm, y=Ginim, size=NHAm)) + geom_point() + facet_wrap(~Comp, scales='free') +
	    scale_size_binned(range=c(1,10))
### Find cellID examples
    CompoCl <- mutate(CompoCl, cellEx=NA)
    for (k in 1:dim(CompoCl)[1]){
      Examp <- CompoCl[k, ]
      tt <- filter(AFilt, Comp==Examp$Comp, cl==Examp$cl) %>% mutate(dGini=Gini-Examp$Ginim, 
	dNHA=NHA-Examp$NHAm, dDg=Dg-Examp$Dgm)
      tt[, 14:16] <- apply(tt[, 14:16], 2, scale)
      tt <- mutate(tt, dAll=sqrt(dGini^2+dNHA^2+dDg^2)) %>% arrange(dAll)
      CompoCl$cellEx[k] <- tt$cellID[1]
    }
### Extract cellID
    LandscapeRaw <- read.asciigrid('DATA/Bauges_Landscape/BaugesInit/cellID.asc')
    MatLandscape <- raster::as.matrix(raster(LandscapeRaw))
    DF <- as.data.table(DF)
    DFAll <- NULL
    for (k in 1:dim(CompoCl)[1]){
	   Examp <- CompoCl[k, ]
	   Focalij <- which(MatLandscape==Examp$cellEx, arr.ind=TRUE)
	   Neighboor <- as.vector(MatLandscape[Focalij[1]:(Focalij[1]+3), 
		   Focalij[2]:(Focalij[2]+3)])
	   DFt <- DF[cellID %in% Neighboor, ]	  
	   Dg <- sqrt(sum(DFt$dbh^2*DFt$n)/sum(DFt$n))
	   Are <- length(unique(DFt$cellID)) * 0.25^2
	   Ba <- pi * sum(DFt$n*(DFt$dbh/200)^2) / Are
	   Gin <- GiniPop(Size=DFt$dbh, BA=DFt$n*DFt$dbh^2, weight=DFt$n)
	   NHa <- sum(DFt$n) / Are
	   DFt <- mutate(DFt, Dg=Dg, BA=Ba, Gini=Gin, NHA=NHa)
	   DFt <- mutate(DFt, Compo=Examp$Comp, Cluster=Examp$cl, Ncells=Examp$Ncells, Area=Are)
	   DFAll <- rbind(DFAll, DFt)
   }
   write.table(file='InitBaugesExamples.csv', DFAll, sep='\t')
   png(file='ClusterInitBauges.png', width=2000,height=2000,res=140)
   print(pl);dev.off()
}

ComputeGini3 <- function(DF){
	require(tictoc)
	require(parallel)
    require(data.table)
    DF <- as.data.table(DF)
    LandscapeRaw <- read.asciigrid('DATA/Bauges_Landscape/BaugesInit/cellID.asc')
    DFstand <- unique(dplyr::select(DF, cellID, Area, XCenter, YCenter))
    MatLandscape <- raster::as.matrix(raster(LandscapeRaw))
    NRow <- nrow(MatLandscape) - 5
    getOutRow <- function(MatLandscape, DF, DFstand, i){
        NCol <- ncol(MatLandscape) - 5
        k <- 1
        JJ <- data.frame(cellID=rep(NA, (NCol-5+1)), Area=NA, NHA=NA, Gini=NA,
            Dg=NA, Sp1= NA, PrSp1=NA, Sp2=NA, PrSp2=NA, Sp3=NA, PrSp3=NA)
        for (j in 5:NCol){
            FocalID <- MatLandscape[i, j]
            if (FocalID %in% DFstand$cellID){
                NeighboorID <- as.vector(MatLandscape[(i+3):i,(j+3):j])
                PopI <- DF[cellID %in% NeighboorID, ]
  	        PopI <- PopI[!is.na(dbh), ]
 	        if (dim(PopI)[1]>0){
                    GiniI <- GiniPop(Size=PopI$dbh, BA=PopI$dbh^2*PopI$n, weight=PopI$n)
   	            Psp <- PopI[, .(pBA=sum(n*dbh^2)), by=sp] %>% mutate(pBA=pBA/sum(pBA))
                    Psp <- arrange(Psp, desc(pBA))
 	            OUT <- data.frame(cellID=MatLandscape[i, j], Area=length(unique(PopI$cellID)) * 0.25^2) %>%
                        mutate(NHA=sum(PopI$n)/Area, Gini=GiniI, Dg=sqrt(sum(PopI$n*PopI$dbh^2)/sum(PopI$n))) %>%
  	                mutate(Sp1=as.character(Psp$sp[1]), PrSp1=Psp$pBA[1],
			Sp2=as.character(Psp$sp[2]), PrSp2=Psp$pBA[2],
		       	Sp3=as.character(Psp$sp[3]), PrSp3=Psp$pBA[3])
	                JJ[k, ] <- OUT
	        }
	        k <- k+1
	    }  
        }
        return(filter(JJ, !is.na(cellID)))
    }
    tic()
    cl <- makeForkCluster(6)
    JJ <- do.call(rbind, parLapply(cl, 5:NRow, getOutRow,
        MatLandscape=MatLandscape, DF=DF, DFstand=DFstand))
    stopCluster(cl)
    toc()
    saveRDS(file='outGini3.Rds', JJ)
}

dec_Matrix <- function(Mat){
    Nrow <- nrow(Mat)
    Ncol <- ncol(Mat)
### 1 grid point away
   Mat10 <- rbind(matrix(rep(rep(NA, Ncol), 1), nrow=1), Mat[(1:(Nrow-1)), 1:Ncol])
   MatN10 <- rbind(Mat[(2:(Nrow)), ], matrix(rep(rep(NA, Ncol), 1), nrow=1))
   Mat01 <- cbind(matrix(rep(rep(NA, Nrow), 1), ncol=1), Mat[, (1:(Ncol-1))]) 
   MatN01 <- cbind(Mat[, (2:Ncol)], matrix(rep(rep(NA, Nrow), 1), ncol=1))
   Mat11 <- cbind(matrix(rep(rep(NA, Nrow), 1), ncol=1), Mat10[, (1:(Ncol-1))]) 
   MatW11 <- cbind(Mat10[, (1:(Ncol-1))], matrix(rep(rep(NA, Nrow), 1), ncol=1)) 
   MatN11 <- cbind(MatN10[, (2:Ncol)], matrix(rep(rep(NA, Nrow), 1), ncol=1))
   MatWN11 <- cbind(matrix(rep(rep(NA, Nrow), 1), ncol=1), MatN10[, (1:(Ncol-1))]) 

### 2 grid points away
   Mat20 <- rbind(matrix(rep(rep(NA, Ncol), 2), nrow=2), Mat[(1:(Nrow-2)), 1:Ncol])
   Mat02 <- cbind(matrix(rep(rep(NA, Nrow), 2), ncol=2), Mat[, (1:(Ncol-2))])
   MatN20 <- rbind(Mat[(3:(Nrow)), ], matrix(rep(rep(NA, Ncol), 2), nrow=2))
   MatN02 <- cbind(Mat[, (3:Ncol)], matrix(rep(rep(NA, Nrow), 2), ncol=2))
   Mat22 <- cbind(matrix(rep(rep(NA, Nrow), 2), ncol=2), Mat20[, (1:(Ncol-2))]) 
   Mat21 <- cbind(matrix(rep(rep(NA, Nrow), 1), ncol=1), Mat20[, (1:(Ncol-1))]) 
   Mat12 <- cbind(matrix(rep(rep(NA, Nrow), 2), ncol=2), Mat10[, (1:(Ncol-2))]) 
   MatN21 <- cbind(MatN20[, (2:Ncol)], matrix(rep(rep(NA, Nrow), 1), ncol=1))
   MatN12 <- cbind(MatN10[, (3:Ncol)], matrix(rep(rep(NA, Nrow), 2), ncol=2))
   MatN22 <- cbind(MatN20[, (3:Ncol)], matrix(rep(rep(NA, Nrow), 2), ncol=2))
   MatW12 <- cbind(Mat10[, (3:Ncol)], matrix(rep(rep(NA, Nrow), 2), ncol=2)) 
   MatW21 <- cbind(Mat20[, (2:Ncol)], matrix(rep(rep(NA, Nrow), 1), ncol=1)) 
   MatW22 <- cbind(Mat20[, (3:Ncol)], matrix(rep(rep(NA, Nrow), 2), ncol=2)) 
   MatWN21 <- cbind(matrix(rep(rep(NA, Nrow), 2), ncol=2), MatN10[, (1:(Ncol-2))]) 
   MatWN12 <- cbind(matrix(rep(rep(NA, Nrow), 1), ncol=1), MatN20[, (1:(Ncol-1))]) 
   MatWN22 <- cbind(matrix(rep(rep(NA, Nrow), 2), ncol=2), MatN20[, (1:(Ncol-2))]) 
   
   ListMat <- c('Mat', 'Mat10', 'MatN10', 'Mat01', 'MatN01',
	'Mat11', 'MatW11', 'MatN11', 'MatWN11', 'Mat20', 'Mat02', 'MatN20', 'MatN02',
	'Mat22', 'Mat21', 'Mat12', 'MatN21', 'MatN12', 'MatN22', 'MatW12', 'MatW21',
	'MatW22', 'MatWN21', 'MatWN12', 'MatWN22')
   MatAll <- NULL
   for (k in 1:length(ListMat)){
       MatAll[[k]] <- eval(parse(text=ListMat[k]))
   }
   return(MatAll)
}
















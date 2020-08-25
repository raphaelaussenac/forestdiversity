#' Create a dummy virtual plot
#'
#' @return A data.frame of a virtual plot
#' @export
BuildDummy <- function(){
    Inter <- 10
    o <- read.csv('DATA/all_profound.csv')
    o <- filter(o, src=='4c',year==1967)
    DF <- mutate(o, X=runif(dim(o)[1],-10,10), Y=runif(dim(o)[1],-10, 10), ClassSize=(1+floor((D_cm-Inter/2)/Inter))*Inter)
    return(DF)
}

BuildDummy2 <- function(){
    Inter <- 10
    o <- read.csv('DATA/all_profound.csv')
    o <- filter(o, src=='4c',year==1967)
    NN <- floor(sqrt(dim(o)[1]))
    o <- o[1:(NN^2),]
    XY <- expand.grid(seq(-10,10,length.out=NN), seq(-10,10,length.out=NN))
    DF <- mutate(o, X=XY[,1], Y=XY[,2], ClassSize=(1+floor((D_cm-Inter/2)/Inter))*Inter)
    return(DF)
}

#' Create a Plot object
#'
#' This function takes a data.frame and return a Plot object. The data.frame
#' must contain variables: species, D_cm, weight, X, Y 
#' 
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
createPlot <- function(DF, shape="quadrat", coord){
    if (("X" %in% colnames(DF) + "Y" %in% colnames(DF) +
	"species" %in% colnames(DF) + 'D_cm' %in% colnames(DF))<3){
          stop('Need coordinates with names X and Y, and species, and DDBH (D_cm)')
        }
    if (sum(is.na(DF[,c('species','X','Y')])) > 0){
	    stop('Watch for the NA values')
    }
    if (dim(DF)[1]==0){stop('Empty data')}
    if (shape=="quadrat" & length(coord)!=4){stop('Need coord=(xmin, xmax, ymin, ymax) for quadrats')}
    if (shape=="circular" & length(coord)!=4){stop('Need coord=R for circular plots')}
    if (!(shape %in% c('quadrat','circular'))){stop('shape must be quadrat or circular')}
  Plot <- list(DF=DF, shape=shape, coord=coord)
  class(Plot) <- append("Plot", class(DF))
  return(Plot)
}

#' Compute distance between trees in a plot
#'
#' This function takes a Plot object and return a TabDis object 
#' @param Plot object, discribing ONE plot (one year, one site, one source,...)
#' @param Nselec int: number of neighboor retained, 10 by default
#' @return A TabDis object
#' @export
TabDist <- function(Plot, Nselec = 10){
    if (!('Plot' %in% class(Plot))){stop('Need a Plot class arg')}
    Plot <- DisToBorder(Plot)
    DF <- Plot$DF
    N <- dim(DF)[1]
    coo <- cbind(DF$X, DF$Y)
    V <- as.data.table(t(combn(seq(1,dim(coo)[1]),2))) 
    calcd <- function(i1, i2, coo){return(sqrt(apply((coo[i1,]-coo[i2,])^2,1,sum)))}

    V <- cbind(V, Dis=calcd(V$V1,V$V2, coo=coo))
    V <- rbind(V, data.table(V1=V$V2, V2=V$V1, Dis=V$Dis))
    Tdis <- setorder(V, V1, Dis)
    if (N >= Nselec){Tdis <- Tdis[rep(0:(Nselec-1), N)+rep(seq(1,N*(N-1),by=(N-1)),each=Nselec)]}
    Tdis <- mutate(Tdis, sp1=DF$species[Tdis$V1], X1=DF$X[Tdis$V1], Y1=DF$Y[Tdis$V1],
        X2=DF$X[Tdis$V2], Y2=DF$Y[Tdis$V2], sp2=DF$species[Tdis$V2],
        DBH1=DF$D_cm[Tdis$V1], DBH2=DF$D_cm[Tdis$V2], ClassSize1=DF$ClassSize[Tdis$V2],
       	ClassSize2=DF$ClassSize[Tdis$V2], DisToBord=DF$DisToBord[Tdis$V1])
    Tdis <- list(DF=Tdis, shape=Plot$shape, coord=Plot$coord)
    class(Tdis) <- append('DistanceTab', class(Tdis))
    return(Tdis)
}

plot.DistanceTab <- function(Tdis, Nk=4){
  Plot <- group_by(Tdis$DF, V1) %>%
	  summarise(species=sp1[1], DBH1=DBH1[1], ClassSize=ClassSize1[1],
	  X=X1[1], Y=Y1[1], AllIn=(Dis[Nk]<=DisToBord[Nk])) %>% ungroup()
  if (Tdis$shape=='quadrat'){
    pl <- ggplot() + theme(text=element_text(size=20)) + 
          geom_point(data=filter(Plot, AllIn==TRUE), aes(x=X, y=Y, fill=species, size= DBH1),pch=21) + 
	  geom_point(data=filter(Plot, AllIn==FALSE), aes(x=X, y=Y, col=species, size= DBH1), pch=1) +
          geom_rect(aes(xmin=Tdis$coord[1], xmax=Tdis$coord[2], ymin=Tdis$coord[3], ymax=Tdis$coord[4]), fill=NA,col='black')
  }else if (Tdis$shape=='circular'){
    pl <- ggplot() + theme(text=element_text(size=20)) + 
          geom_point(data=filter(Plot, AllIn==TRUE), aes(x=X, y=Y, fill=species, size= DBH1),pch=21) + 
	  geom_point(data=filter(Plot, AllIn==FALSE), aes(x=X, y=Y, col=species, size= DBH1), pch=1)
  }
  print(pl)
}

#' Compute distance between trees and border
#'
#' This function takes a Plot object and return a TabDis object 
#' @param Plot object
#' @return A Plot object with DisToBorder variable in the DF
DisToBorder <- function(Plot){
  if (!('Plot' %in% class(Plot))){stop('Need a Plot class arg')}
  DF <- Plot$DF
  if (Plot$shape=='circular'){
    DisToCenter <- sqrt(DF$X^2+DF$Y^2)
    DF <- mutate(DF, DisToBord=Plot$coord-DisTocenter)
  }else if (Plot$shape=='quadrat'){
    DisToBord <- abs(cbind(DF$X - Plot$coord[1], DF$X - Plot$coord[2],
	      DF$Y - Plot$coord[3], DF$Y - Plot$coord[4]))
    DisT <- apply(DisToBord, 1, min)
    DF <- mutate(DF, DisToBord=DisT)
  }
  return(list(DF=DF, shape=Plot$shape, coord=Plot$coord))
}

#' Compute Species mingling in a plot
#'
#' This function takes a TabDis object and return mingling metrics
#' @param Tabdis object
#' @param Nk int, number of neighboors
#' @param EdgeCorrection str, type of edge correction applied
#' @return data.frame with mingling metrics
#' @export
Compute_mingling <- function(TabDis, Nk=4, EdgeCorrection="NN1"){
	if (Nk<=0|Nk>10){stop('Number of k neighbours Nk must lie between 1 and 10')}
    TabDis <- mutate(TabDis, I1=sp1!=sp2)
    if (TabDis$shape=='circular'){
        Mg <- group_by(TabDis, V1) %>% summarise(k=sum(I1[1:Nk]), Mi=k/Nk, AllIn=(Dis[Nk]<=DisToBord[Nk])>0,
           Fi=(TabDis$coord-Dis[Nk])^2*pi) %>% ungroup()
    }else if (TabDis$shape=='quadrat'){
        Mg <- group_by(TabDis, V1) %>% summarise(k=sum(I1[1:Nk]), Mi=k/Nk, AllIn=(Dis[Nk]<=DisToBord[Nk])>0,
           Fi=(abs(diff(TabDis$coord[c(1,2)]))-Dis[Nk])*(abs(diff(TabDis$coord[c(3,4)]))-Dis[Nk])) %>% ungroup()
    }
    switch(EdgeCorrection,
        NN1={Lhat <- sum(Mg$AllIn/Mg$Fi)
             Mk <- group_by(Mg, k) %>% summarise(mk=(1/Lhat) * sum(AllIn/Fi)) %>%
		ungroup() %>% mutate(M=sum(Mg$Mi*Mg$AllIn/Mg$Fi) / Lhat)
        },
        NN2={Ni <- sum(Mg$AllIn)
             Mk <- group_by(Mg,k) %>% summarise(mk=sum(AllIn)/Ni)  %>%
	        ungroup() %>% mutate(M=sum(Mg$Mi*Mg$AllIn)/sum(Mg$AllIn))
        },
	Exclude={Mg <- filter(Mg, AllIn==TRUE) 
              Mk <- group_by(Mg, k) %>% summarise(mk=n()/dim(Mg)[1]) %>%
	      	ungroup() %>% mutate(M=sum(Mg$Mi) / dim(Mg)[1])
	},
	None={Mk <- group_by(Mg, k) %>% summarise(mk=n()/dim(Mg)[1]) %>%
	      	ungroup() %>% mutate(M=sum(Mg$Mi) / dim(Mg)[1])
	},
        stop("Need a valid edge correction (NN1, NN2, Exclude, None)")
    )
    TT <- group_by(TabDis, V1) %>% summarise(species=sp1[1],ClassSize=ClassSize1[1]) %>% ungroup()
    Ni <- group_by(TT, species) %>% summarise(N=n()) %>% ungroup() %>% mutate(Nt=dim(TT)[1])
    Em <-  sum(Ni$N * (Ni$Nt - Ni$N) / (Ni$Nt * (Ni$Nt-1)))
    return(mutate(Mk,Em=Em))
}

#' Compute Species mingling in a plot
#'
#' This function takes a TabDis object and return mingling metrics
#' @param TabDis object
#' @param Nk int, number of neighboors
#' @param EdgeCorrection str, type of edge correction applied
#' @return data.frame with size differentiation indexes
#' @export
Compute_Size_Diff <- function(TabDis, Nk=4, EdgeCorrection='None'){
	if (Nk<=0|Nk>10){stop('Number of k neighbours Nk must lie between 1 and 10')}
      N <- length(unique(TabDis$V1))
      TabDis <- mutate(TabDis, I=rep(1:10, N)) %>% filter(I<=Nk)
      TabDis <- mutate(TabDis, DBHmax=apply(cbind(TabDis$DBH1, TabDis$DBH2), 1, max),
	DBHmin=apply(cbind(TabDis$DBH1, TabDis$DBH2), 1, min))
    if (TabDis$shape=='circular'){
        SiDi <- group_by(TabDis, V1) %>% summarise(Ti=1-sum(DBHmin/DBHmax) / Nk,
           AllIn=(Dis[Nk]<=DisToBord[Nk])>0,
           Fi=(TabDis$coord-Dis[Nk])^2*pi) %>% ungroup()
    }else if (TabDis$shape=='quadrat'){
        SiDi <- group_by(TabDis, V1) %>% summarise(Ti=1-sum(DBHmin/DBHmax)/Nk,
           AllIn=(Dis[Nk]<=DisToBord[Nk])>0,
           Fi=(abs(diff(TabDis$coord[c(1,2)]))-Dis[Nk])*(abs(diff(TabDis$coord[c(3,4)]))-Dis[Nk])) %>% ungroup()
    }
    switch(EdgeCorrection,
        NN1={Lhat <- sum(SiDi$AllIn/SiDi$Fi)
             T <- sum(SiDi$Ti*SiDi$AllIn/SiDi$Fi) / Lhat
        },
        NN2={T <- sum(SiDi$Ti*SiDi$AllIn) / sum(SiDi$AllIn)
        },
	Exclude={SiDi <- filter(AllIn==TRUE)
	  T <- mean(SiDi$Ti)
	},
	None={T <- mean(SiDi$Ti)
	},
        stop("Need a valid edge correction (NN1, NN2, Exclude, None)")
    )
    TT <- group_by(TabDis, V1) %>% summarise(DBH=DBH1[1]) %>% ungroup() %>% arrange(DBH)
    TT <- cbind(TT,R=c(0, cumsum(TT$DBH)[1:(N-1)]))
    ET <- 1 - 2/(N*(N-1)) * sum(TT$R/TT$DBH)
    return(cbind(T, ET=ET))
} 

#' Compute Uniform angle Index
#'
#' This function takes a TabDis object and return a mean Winkelmass value
#' @param TabDis object
#' @param Nk int, number of neighboors
#' @param EdgeCorrection str, type of edge correction applied
#' @return numeric mean Winkelmass value for the plot
#' @export
Compute_Winkelmass <- function(TabDis, Nk=4){
	if (Nk<=0|Nk>10){stop('Number of k neighbours Nk must lie between 1 and 10')}
   DF <- group_by(TabDis$DF, V1) %>% mutate(N=1:10, AllIn=(Dis[Nk]<=DisToBord)) %>% ungroup() %>% filter(N<=Nk)
   iN <- combn(Nk,2);iN <- cbind(iN,rbind(iN[2,], iN[1,]))

   listAnglesOld <- function(X1, Y1, X2, Y2, iN, Nk, V1, AllIn){
	   X2 <- X2 + runif(length(X2),-1e-4,1e-4) # Add noise to handle perfectly aligned trees
	   Y2 <- Y2 + runif(length(X2),-1e-4,1e-4)
	   if (AllIn[1]==FALSE){return(NA)}
           D1 <- cbind(X2[iN[1,]], Y2[iN[1,]]) - cbind(X1[iN[1,]], Y1[iN[1,]])
           D2 <- cbind(X2[iN[2,]], Y2[iN[2,]]) - cbind(X1[iN[2,]], Y1[iN[2,]])
	   Ags <- (atan2(D2[,2], D2[,1]) - atan2(D1[,2], D1[,1])) * 180 / pi
	   Ags[Ags<0] <- 360 + Ags[Ags<0]
	   T <- as.data.frame(t(iN)) %>% mutate(Ang=Ags)
	   p1 <- filter(T, V1==1) %>% filter(Ang==min(Ang))
	   p2 <- filter(T, V1==p1$V2) %>%  filter(Ang==min(Ang))
	   p3 <- filter(T, V1==p2$V2) %>%  filter(Ang==min(Ang))
	   p4 <- filter(T, V1==p3$V2, V2==1)
	   listAng <- c(p1$Ang, p2$Ang, p3$Ang, p4$Ang)
	   if (round(sum(listAng))!=360){stop(paste0('Trouble with angle
	     calculation with tree Nb ', V1[1]))}
	   listAng <- round(listAng, digits=1)
	   listAng[listAng>180] <- 360 -  listAng[listAng>180]
	   return(sum(listAng<=(360 / Nk)))
   }
	   Tab <- 1 + permutations(3,3)
	   Tab <- cbind(10+Tab[,1],Tab[,1]*10+Tab[,2],Tab[,2]*10+Tab[,3],Tab[,3]*10+1)

  listAngles <- function(X1, Y1, X2, Y2, iN, Nk, V1, AllIn, Tab){
	   if (AllIn[1]==FALSE){return(NA)}
	   X2 <- X2 + runif(length(X2),-1e-4,1e-4) # Add noise to handle perfectly aligned trees
	   Y2 <- Y2 + runif(length(X2),-1e-4,1e-4)
           D1 <- cbind(X2[iN[1,]], Y2[iN[1,]]) - cbind(X1[iN[1,]], Y1[iN[1,]])
           D2 <- cbind(X2[iN[2,]], Y2[iN[2,]]) - cbind(X1[iN[2,]], Y1[iN[2,]])
	   Ags <- (atan2(D2[,2], D2[,1]) - atan2(D1[,2], D1[,1])) * 180 / pi
	   Ags[Ags<0] <- 360 + Ags[Ags<0]
	   T <- as.data.frame(t(iN)) %>% mutate(Ang=Ags, I=10*V1 + V2)
	   Tabnew <- Tab
           Tabnew[] <- T$Ang[match(unlist(Tab), T$I)]
	   In <- which(round(apply(Tabnew,1, sum))==360)
	   if (length(In)!=1){stop(paste0('Trouble with angle
             calculation with tree Nb ', V1[1]))}
	   listAng <- Tabnew[In,]
	   listAng <- round(listAng, digits=1)
	   listAng[listAng>180] <- 360 -  listAng[listAng>180]
	   return(sum(listAng<=(360 / Nk)))
   }
	 return(group_by(DF, V1) %>% summarise(Wink=listAngles(X1, Y1, X2, Y2,
	   iN=iN, Nk=Nk, V1, AllIn, Tab=Tab)/Nk, AllIn=AllIn[1]) %>% ungroup())
}

############ Structural complexity

##' Compute triangle surfaces
#'
#' This function takes deldir object (plot tesselation) and compute sum of triangle surfaces
#' @param deldir object
#' @return data.frame with triangle surfaces
TriangleSurface <- function(TriDF){
  Surface_Proj <- abs(det(t(matrix(c(TriDF$x,TriDF$y,1,1,1),ncol=3)))) / 2

  Surface_3d <- sqrt(det(t(matrix(c(TriDF$x,TriDF$y,1,1,1),ncol=3)))^2 +
		     det(t(matrix(c(TriDF$y,TriDF$z,1,1,1),ncol=3)))^2 +
		     det(t(matrix(c(TriDF$z,TriDF$x,1,1,1),ncol=3)))^2) / 2
  return(data.frame(Surface_Proj, Surface_3d))
}

##' Compute Structural Complexity Index
#'
#' This function takes a Plot object and return structural complexity index
#' @param Plot object
#' @return SCI num
#' @export
StructuralComplexityIndex <- function(Plot){
  if (dim(Plot$DF)[1]<3){stop("Need at least 3 points to build triangles")}
  if (!('H_m' %in% colnames(Plot$DF))){stop('Need tree heights (H_m)')}
  Del <- deldir::deldir(Plot$DF$X, Plot$DF$Y, z=Plot$DF$H_m)
  Triangles <- deldir::triang.list(Del)
  Surf <- do.call(rbind,lapply(Triangles, TriangleSurface))
  SCI <- sum(Surf$Surface_3d)/sum(Surf$Surface_Proj)
  return(SCI)
}

#########################
#Tdis <- TabDist(DF)

#Mg <- Compute_mingling(TabDis, coord=coord)
#SiDi <- Compute_Size_Diff(TabDis, coord=coord)












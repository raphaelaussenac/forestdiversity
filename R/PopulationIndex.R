#' Compute diversity metrics
#'
#' This function compute the Gini index and Hill numbers for a population 
#' @param dataSet, data.frame of population, can gather multiple sites, years and sources
#' @param N_var string, name of the variable used to compute indices
#' @param ClassInter num, Class size (same units as the variable chosen)
#' @param ClassIni num, min dbh of the first class
#' @param type string, define wether to compute proportion in terms of frequencies or relative basal area
#' @param Out string, define wether to compute Hillnb or entropy based indices ('HillNb' or 'Entropy'
#' @return A data.frame containing diversity metrics for each site/year/source
#' @export
CalcDivIndex <- function(dataSet, Nvar = 'D_cm', ClassInter = 10, ClassIni=7.5, type = 'BA', Out='HillNb'){
    if (!(Out %in% c('HillNb', 'Entropy'))){stop("Out should be 'HillNb' or 'Entropy'")}
    dataSet <- dplyr::mutate(dataSet, Var=dataSet[[Nvar]])
    if (Nvar!='species'){
            if(min(dataSet$Var) < ClassIni){
	    warning('Careful some trees have Size below ClassIni, they will be discarded')
	    dataSet <- dplyr::filter(dataSet, Var>=ClassIni)
	}
        ClassInterval <- seq(ClassIni, max(dataSet$Var) + 2* ClassInter, by=ClassInter) 
        ClassName <- ClassInterval + mean(diff(ClassInterval)) / 2
        dataSet <- dplyr::mutate(dataSet, Class=ClassName[findInterval(Var, ClassInterval)])
    }else{
        dataSet <- dplyr::mutate(dataSet, Class=as.character(Var))
    }
    if (is.null(dataSet[["Var"]])){stop('Need to choose a correct variable first')}
    if (!('site' %in% names(dataSet))){dataSet <- dplyr::mutate(dataSet, site='NA')}
    if (!('weight' %in% names(dataSet))){dataSet <- dplyr::mutate(dataSet, weight=1)}
    listNameGrouping <- 'site'
    if (('year' %in% names(dataSet))){listNameGrouping <- c(listNameGrouping, 'year')}
    if (('src' %in% names(dataSet))){listNameGrouping <- c(listNameGrouping, 'src')}
    if (('postThinning' %in% names(dataSet))){listNameGrouping <- c(listNameGrouping, 'postThinning')}
    if (('postDisturbance' %in% names(dataSet))){listNameGrouping <- c(listNameGrouping, 'postDisturbance')}
    if (!('src' %in% names(dataSet))){dataSet <- dplyr::mutate(dataSet, src='NA')}
    dataSet <- dplyr::mutate(dataSet, BA=pi*(D_cm/200)^2*weight)
    dataSet <- data.table::as.data.table(dataSet)
    if (type=='BA'){
        if (Nvar!="species"){
	     DivIndex <- dataSet[, .(Gini=GiniPop(Var, BA, weight),
                 H=DivPop(Class, BA, OutFormat='str', Out=Out)), by=listNameGrouping]
             DivIndex <- dplyr::mutate(DivIndex, 
                        Nclass=as.numeric(do.call(rbind, strsplit(H, '/'))[,1]),
                        Shannon=as.numeric(do.call(rbind, strsplit(H, '/'))[, 2]),
                        Simpson=as.numeric(do.call(rbind, strsplit(H, '/'))[, 3]))
	     DivIndex <- dplyr::select(DivIndex, -H)

	}else{
             DivIndex <- dataSet[, .(H=DivPop(Class, BA, OutFormat='str', Out=Out)), by=listNameGrouping]
             DivIndex <- dplyr::mutate(DivIndex,  
                        Nclass=as.numeric(do.call(rbind, strsplit(H, '/'))[,1]),
                        Shannon=as.numeric(do.call(rbind, strsplit(H, '/'))[, 2]),
                        Simpson=as.numeric(do.call(rbind, strsplit(H, '/'))[, 3]))
	     DivIndex <- dplyr::select(DivIndex, -H)
	}
    }else{
        if (Nvar!='species'){
	     DivIndex <- dataSet[, .(Gini=GiniPop(Var, BA, weight),
                 H=DivPop(Class, weight, OutFormat='str', Out=Out)), by=listNameGrouping]
             DivIndex <- dplyr::mutate(DivIndex,
                        Nclass=as.numeric(do.call(rbind, strsplit(H, '/'))[,1]),
                        Shannon=as.numeric(do.call(rbind, strsplit(H, '/'))[, 2]),
                        Simpson=as.numeric(do.call(rbind, strsplit(H, '/'))[, 3]))
	     DivIndex <- dplyr::select(DivIndex, -H)
	}else{
             DivIndex <- dataSet[, .(H=DivPop(Class, weight, OutFormat='str', Out=Out)), by=listNameGrouping]
             DivIndex <- dplyr::mutate(DivIndex,
                        Nclass=as.numeric(do.call(rbind, strsplit(H, '/'))[,1]),
                        Shannon=as.numeric(do.call(rbind, strsplit(H, '/'))[, 2]),
                        Simpson=as.numeric(do.call(rbind, strsplit(H, '/'))[, 3]))
	     DivIndex <- dplyr::select(DivIndex, -H)
	}
    }
    if (Out=='HillNb'){
        names(DivIndex)[which(names(DivIndex)=="Shannon")] <- "H1"
        names(DivIndex)[which(names(DivIndex)=="Simpson")] <- "H2"
    }
    return(DivIndex)
}

#' Compute Gini index
#'
#' This function compute the Gini index for a population based on the BA proportion 
#' @param Size numeric vector, size of each treee in the population (in cm)
#' @param BA numeric vector, basal area associated with each tree (in m2ha-1)
#' @param Weight, numeric vector, weight associated with each tree
#' @param PLOT, logical, if set to 1 plot the lorenz curve
#' @return The Gini index for the population
#' @export
GiniPop <- function(Size, BA, weight = 1, PLOT=FALSE){
	# We handle cases where several trees have same size
    DF <- dplyr::group_by(data.frame(S=Size, B=BA, W=weight), S)
    DF <- dplyr::summarise(DF, B=sum(B), W=sum(W))
    DF <- dplyr::ungroup(DF)
    Size <- DF$S; BA=DF$B; weight=DF$W
    oSize <- order(Size)
    Size <- Size[oSize]
    BA <- BA[oSize]
    weight <- weight[oSize]
    ## Lorenz Curve
    x <- cumsum(weight)/sum(weight) # CDF share of pop
    y <- cumsum(BA)/sum(BA) # CDF share of  BA
    x <- c(0, x)
    y <- c(0, y)
    if (length(x)==1){
      return(NA)
    }else{
      dx <- diff(x)
    }
    A <- sum(c(x[1],dx)*(y+c(0,y[1:(length(y)-1)]))/2) # Area under the Lorenz Curve
    if (PLOT==TRUE){
        pl <- ggplot2::ggplot(data.frame(x=x,y=y), ggplot2::aes(x=x,y=y)) +
	    ggplot2::geom_area() + ggplot2::geom_abline(col='red') +
	    ggplot2::xlab('Cumulative proportion of size') +
	    ggplot2::ylab('Cumulative proportion of basal area') +
	    ggplot2::theme_bw(base_size=20) +
	    ggplot2::ggtitle(paste0('Gini coefficient : ', round(1-2*A, digits=2)))
        print(pl)
    }
    return(1-2*A)
}

##################################################
###################################################
#' Compute Hill Nb
#'
#' This function compute the Hill numbers index for a population 
#' @param Class numeric vector, class of each treee in the population (in cm)
#' @param Weight, numeric vector, weight associated with each tree, can be the basal area
#' @return The Hill Numbers index for the population
#' @export
HillPop <- function(Class, weight, OutFormat='list'){
    In <- EntropyPop(Class, weight, OutFormat='list')
    Out <- data.frame(Nclass=In$Nclass, H1=exp(In$Shannon), H2=1/(In$Simpson))
    if (OutFormat=='list'){
        return(Out)
    }else if (OutFormat=='str'){
        Out <- paste(Out$Nclass, Out$H1, Out$H2, sep='/') 
        return(Out)
    }else{
	    stop('Outformat must be list or str')
    }
}

##################################################
###################################################
#' Compute Entropy indices
#'
#' This function compute the Hill numbers index for a population 
#' @param Class numeric vector, class of each treee in the population (in cm)
#' @param Weight, numeric vector, weight associated with each tree, can be the basal area
#' @return The Hill Numbers index for the population
#' @export
EntropyPop <- function(Class, weight, OutFormat='list'){
    P <- dplyr::group_by(data.frame(Class, weight), Class)
    P <- dplyr::summarise(P, p = sum(weight))
    P <- dplyr::ungroup(P)
    P <- dplyr::mutate(P, p=p/sum(p))
    Out <- data.frame(Shannon=-sum(P$p*log(P$p)), Simpson=sum(P$p^2), Nclass=dim(P)[1])
    if (OutFormat=='list'){
        return(Out)
    }else if (OutFormat=='str'){
        Out <- paste(Out$Nclass, Out$Shannon, Out$Simpson, sep='/') 
        return(Out)
    }else{
	    stop('Outformat must be list or str')
    }
}

##################################################
###################################################
#' Wrapper function for Div indices
#'

DivPop <- function(Class, weight, OutFormat, Out){
    if (Out=='HillNb'){return(HillPop(Class=Class, weight=weight, OutFormat=OutFormat))}
    if (Out=='Entropy'){return(EntropyPop(Class=Class, weight=weight, OutFormat=OutFormat))}
}

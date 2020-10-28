#' Return Diversity metrics
#'
#' This function loads a according to the evalSite value.
# It return the diversity metrics for each year * site in a data.frame
#' @param dataSet data.frame
#' @param Nvar string, name of the diversity axis (D_cm for size, species for species)
#' @param Inter int, width of each size class category
#' @param path string, path of the data stored
#' @return A data.frame of metrics
#' @export
ReturnDivIndex <- function(evalSite, Inter=10, path='DATA'){
    if (!file.exists(paste0(path,'/','all_',evalSite,'.csv'))){stop('Need to build dataset first')}
    dataTemp <- read.csv(file=paste0(path,'/','all_',evalSite,'.csv'))
    DivIndexSize <- CalcDivIndex(dataTemp, Nvar="D_cm", Inter=Inter) %>%
	    dplyr::mutate(ShSize=Sh, GSSize=GS, SimpSize=Simp, GiSize=GI, NClassSize=NClass) %>%
	    dplyr::select(-Sh, -GS, -Simp, -GI, -NClass)
    DivIndexSp <- CalcDivIndex(dataTemp, Nvar="species", Inter=Inter) %>%
	    dplyr::mutate(ShSp=Sh, GSSp=GS, SimpSp=Simp, NSp=NClass) %>%
	    dplyr::select(-Sh, -GS, -Simp, -GI, -NClass)
    DivIndex <- dplyr::left_join(DivIndexSize, DivIndexSp, by=c('year','site','src'))
    return(as.data.frame(DivIndex))
}


#' Compute diversity metrics
#'
#' This function compute the Gini index and Hill numbers for a population 
#' @param dataSet, data.frame of population, can gather multiple sites, years and sources
#' @param Type string, define wether to compute proportion in terms of frequencies or relative basal area
#' @return A data.frame containing diversity metrics for each site/year/source
#' @export
CalcDivIndex <- function(dataSet, Nvar = 'D_cm', Inter = 10, type = 'BA'){
    dataSet <- dplyr::mutate(dataSet, Var=dataSet[[Nvar]])
    if (Nvar %in% c('D_cm','H_m','V_m3')){
	    dataSet <- dplyr::mutate(dataSet, Class = (1+floor((Var-Inter/2)/Inter))*Inter)
    }else if (Nvar=='species'){
	    dataSet <- dplyr::mutate(dataSet, Class=as.character(Var))
    }else{
	    stop('Nvar specified not in output')
    }

    if (is.null(dataSet[["Var"]])){stop('Need to choose variable first')}
    dataSet <- dplyr::mutate(dataSet, BA=pi*(D_cm/200)^2*weight)
    PClass <- dplyr::group_by(dataSet, year, site, src) %>% dplyr::mutate(N = sum(weight), BAt=sum(BA)) %>%
      dplyr::ungroup() %>% dplyr::group_by(year, Class, site, src) %>% dplyr::summarise(p=(sum(weight)/N[1]),
      pBA=sum(BA)/BAt[1], BAt=BAt[1]) %>% dplyr::ungroup()
    if (type=='BA'){PClass$p <- PClass$pBA}
    HillNB <- dplyr::group_by(PClass, year, site, src) %>% dplyr::summarise(Sh=-sum(p * log(p)),
	    NClass=n(), GS=1-sum(p^2), Simp=sum(p^2)) %>% dplyr::ungroup()
    if (is.numeric(dataSet$Var)){
      GiniIndex <- dplyr::group_by(dataSet, year, site, src) %>% dplyr::summarise(GI=Gini(Var,BA,weight)) %>% dplyr::ungroup()
      DivIndex <- dplyr::left_join(HillNB, GiniIndex, by=c('year','site','src'))
    }else{
      DivIndex <- dplyr::mutate(HillNB, GI=NA) 
    }
    return(DivIndex)
}

#' Compute Gini index
#'
#' This function compute the Gini index for a population 
#' @param Size numeric vector, size of each treee in the population (in cm)
#' @param BA numeric vector, basal area associated with each tree (in m2ha-1)
#' @param Weight, numeric vector, weight associated with each tree
#' @return The Gini index for the population
#' @export
Gini <- function (Size, BA = rep(1, length = length(x)), weight = rep(1, length = length(x))){
    oSize <- order(Size)
    Size <- Size[oSize]
    BA <- BA[oSize]
    weight <- weight[oSize]
    ## Lorenz Curve
    x <- cumsum(weight)/sum(weight) # CDF share of pop
    y <- cumsum(BA)/sum(BA) # CDF share of  BA
    if (length(x)==1){
      dx <- 0
    }else{
      dx <- diff(x)
    }
    A <- sum(c(x[1],dx)*(y+c(0,y[1:(length(y)-1)]))/2) # Area under the Lorenz Curve
    return(1-2*A)
}



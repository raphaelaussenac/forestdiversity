
#' Compute resilience metrics for a dataSet
#'
#' This function compute resilience metrics for a dataset 
#' @param dataSet, data.frame formatted using format_Salem or format_Samsara with multiple simulation conditions
#' @param Nvar string, name of the variable studied
#' @return A data.frame containing resilience metrics for each simulation. A NA value is returned if there is an error
#' -99 means that recovery was not complete and metrics could not be computed
#' @export
EventResilience <- function(dataSet, Nvar='V_m3', RecTime=20){
    if (!(Nvar %in% names(dataSet))){stop('Need a variable from dataSet')}
    if (!('VirtualExperiment' %in% class(dataSet))){stop('Need format first')}
    dataSetini <- dataSet[preDisturbance==TRUE, ]
    dataSetini <- dplyr::select(dataSetini, -postThinning, -postDisturbance, -preDisturbance)
    dataSetini <- dataSetini[, lapply(.SD, mean), by=c('site','src')]
    names(dataSetini)[!(names(dataSetini) %in% c('year','site','src'))] <-
        paste0(names(dataSetini)[!(names(dataSetini) %in% c('year','site','src'))], 'ini')

    dataSet <- dplyr::mutate(dataSet, Var=dataSet[, ..Nvar][[1]])
    TT <- dataSet[,.(RecMet=tryCatch(RecoveryMetrics(Var, year, preDisturbance, RecTime=RecTime, FormatOut='str'),
       error=function(e){return('-1/-1/-1')})), by='site']
    TT <- dplyr::mutate(TT, Theta=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[, 1]),
	TimeRec=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[,2]),
	DegRec=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[,3]))
    TT <- dplyr::select(TT, -RecMet)
    TT[TT==-1] <- NA
    TT <- dplyr::mutate(TT, Nvar=Nvar)
    TT <- merge(TT, dataSetini, by="site")
    return(TT)
}

#' Compute resilience metrics for one simulation
#'
#' This function compute resilience metrics for list of variables
#' @param Var, numeric the variable studied
#' @param Year, numeric, the years corresponding to var
#' @param PreDisturb, logical, whether the corresponding year is before perturbation 
#' @param RecTime, integer, Nb of years to compute recovery degree
#' @param FormatOut, str, form of output
#' @return resilience metrics in a list or a string
#' @export
RecoveryMetrics <- function(Var, Year, PreDisturb, RecTime=20, FormatOut='list'){
    VarBef <- mean(Var[PreDisturb==TRUE])
    YearBef <- max(Year[PreDisturb==TRUE])
    if (length(VarBef)==0){stop('No value before perturbation')}
    if (is.na(VarBef)){stop('Something wrong before perturbation')}
    VarPost <- Var[PreDisturb==FALSE]
    if (length(VarPost)==0){stop('No value after perturbation')}
    Intens <- VarPost[1] / VarBef
    if (Intens >= 1){warning('Perturbation did not reduce the variable')}
    dT <- Year[which(Var >= VarBef & PreDisturb==FALSE)][1] - Year[PreDisturb==TRUE]
    VardT <- Var[which(Year>=(YearBef + RecTime))][1]
    Out <- data.frame(Theta=Intens/dT, TimeRec=dT, DegRec=VardT/VarBef)
    if (is.na(dT)){Out$Theta <- -99;Out$TimeRec <- -99}
    if (is.na(VardT)){Out$DegRec <- -99}
    if (FormatOut=='list'){
	return(Out)
    }else if (FormatOut=='str'){
	return(paste(Out$Theta, Out$TimeRec, Out$DegRec, sep='/'))
    }else{
        stop('FormatOut should be list or str')
    }
}

#' Format Samsara dataset
#'
#' This function format the data for further analysis
#' @param dataRaw, data.frame, dataSet containing the simulation output of Samsara
#' @return dataSet, data.table, formatted outuput
format_samsara_Tree <- function(dataRaw){
    dataSet <- dplyr::select(dataRaw, -Int.Energy..MJ.year.1., -Pot.Energy..MJ.year.1.,
	-quality, -marketValue..euros., -dendroHabitats, -ecologicalScore, -treeId, -speciesCode)
    names(dataSet) <- c('site', 'year', 'eventName', 'species', 'X', 'Y', 'D_cm', 'H_m', 'V_m3')
}

#' Format Samsara dataset2
#'
#' This function format the data for further analysis
#' @param dataRaw, data.frame, dataSet containing the simulation output of Samsara
#' @return dataSet, data.table, formatted outuput
format_samsara_Pop <- function(dataRaw){
    dataSet <- dplyr::select(dataRaw, simulationId, eventName, year, speciesName,
        N_ha, G_ha, V_ha, Dg, Dm, Carbon_AG_ha, Carbon_BG_ha, L4Cover,
	NSaplings_ha, NRecruits_ha_yr, NDead_ha_yr, cut_N_ha, cut_V_ha, deadCut_N_ha, deadCut_V_ha) %>%
        dplyr::filter(speciesName=='AllSpecies')
    names(dataSet) <- c('site', 'Event', 'year', 'species', 
        'N', 'G', 'V_m3', 'Dg', 'Dm', 'Carbon_AG_ha', 'Carbon_BG_ha',
	'L4cover', 'NSaplings_ha', 'NRecruits_ha_yr' , 'NDead_ha_yr',
       	'cut_N_ha', 'cut_V_ha', 'deadCut_N_ha', 'deadCut_V_ha')
    dataSet <- dplyr::filter(dataSet, Event!='SaplingCreation')
    dataSet <- dplyr::mutate(dataSet, EventPerturb=substr(Event, 1, 7)=="Disturb",
	PreDisturb=substr(Event, 1, 8)=="Creation")
    return(dataSet)   
}

#' Format Salem dataset
#'
#' This function format the data for further analysis
#' @param dataRaw, data.frame, dataSet containing the simulation output of Salem
#' @param ClassInter, num, parameter used to compute diversity indices
#' @param ClassIni, num, parameter used to compute diversity indices
#' @param type, str, parameter used to compute diversity indices
#' @return dataSet, data.table, formatted outuput
#' @export
format_salem <- function(dataRaw, ClassInter=10, ClassIni=7.5, Out='HillNb', type='BA'){
    listNameGrouping <- c('year', 'site', 'src', 'postThinning', 'postDisturbance')
    if (!('year' %in% names(dataRaw))){stop('Need at least severak years')}
    if (!('site' %in% names(dataRaw))){dataRaw <- dplyr::mutate(dataRaw, site='NA')}
    if (!('src' %in% names(dataRaw))){dataRaw <- dplyr::mutate(dataRaw, src='salem')}
    if (!('postThinning' %in% names(dataRaw))){dataRaw <- dplyr::mutate(dataRaw, postThinning=FALSE)}
    if (!('postDisturbance' %in% names(dataRaw))){stop('Need logical variable postDisturbance')}
    dataRaw <- data.table::as.data.table(dataRaw)
    dataRaw$postThinning <- as.logical(dataRaw$postThinning)
    dataRaw$postDisturbance <- as.logical(dataRaw$postDisturbance)
    HetIndexSize <- CalcDivIndex(dataRaw, 'D_cm', ClassInter=ClassInter, ClassIni=ClassIni, type=type)
    names(HetIndexSize)[!(names(HetIndexSize) %in% listNameGrouping)] <- 
	   paste0(names(HetIndexSize)[!(names(HetIndexSize) %in% listNameGrouping)], 'Size')
    HetIndexSp <- CalcDivIndex(dataRaw, 'species', ClassInter=10, ClassIni=ClassIni, type=type)
    names(HetIndexSp)[!(names(HetIndexSp) %in% listNameGrouping)] <-
         paste0(names(HetIndexSp)[!(names(HetIndexSp) %in% listNameGrouping)], 'Sp')
    dataSet <- dataRaw[, .(V_m3=sum(weight*V_m3), Dg=sqrt(sum(weight*D_cm^2)/sum(weight)),
	 N=sum(weight), BA=pi*sum(weight*(D_cm/200)^2)), by=listNameGrouping]
    dataSet <- merge(dataSet, HetIndexSize, by=listNameGrouping, all.x=TRUE)
    dataSet <- merge(dataSet, HetIndexSp, by=listNameGrouping, all.x=TRUE)
    dataSet <- dplyr::mutate(dataSet, BAinc=c(diff(BA), NA))
    class(dataSet) <- append('VirtualExperiment', class(dataSet))
    dataSet <- dplyr::mutate(dataSet, preDisturbance=FALSE)
    dataSet <- dataSet[, preDisturbance:=(year<=min(year) & postDisturbance==FALSE), by=list(site, src)]
    return(dataSet)
}

#' Plot virtual experiment
#'
#' This function plot a variable from formatted dataSet
#' @param dataSet, formatted dataset (VirtualExperiment object)
#' @param Nvar, string, dataSet, variable to be plotted
#' @return plot Variable chosen
#' @export
plot.VirtualExperiment <- function(dataSet, Nvar='BA'){
    if (!(Nvar %in% names(dataSet))){stop('Need a variable from dataSet')}
    if (!('VirtualExperiment' %in% class(dataSet))){stop('Need format first')}
    dataSet <- dplyr::mutate(dataSet, Var=dataSet[, ..Nvar][[1]])
    pl <- ggplot2::ggplot(dataSet, ggplot2::aes(x=year, y=Var, col=preDisturbance)) +
        ggplot2::geom_point() + ggplot2::facet_wrap(~site, scales='free') +
       	ggplot2::geom_line(col='black') +
       	ggplot2::ggtitle(paste0('Variable : ', Nvar, ' / Model : ', dataSet$src[1])) +
	ggplot2::theme_bw(base_size=16) + ggplot2::ylab('')
    print(pl)
}

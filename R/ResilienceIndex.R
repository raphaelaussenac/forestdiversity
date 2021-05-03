
#' Compute resilience metrics for a dataSet
#'
#' This function compute resilience metrics for a dataset 
#' @param dataSet, data.frame formatted using format_Salem or format_Samsara with multiple simulation conditions
#' @param Nvar string, name of the variable studied
#' @return A data.frame containing resilience metrics for each simulation. A NA value is returned if there is an error
#' -99 means that recovery was not complete and metrics could not be computed
#' @export
EventResilience <- function(dataSet, Nvar='V_m3', RecTime=20, normalize='baseline'){
    if (!(Nvar %in% names(dataSet))){stop('Need a variable from dataSet')}
    if (!('VirtualExperiment' %in% class(dataSet))){stop('Need format first')}
    dataSetini <- dataSet[preDisturbance==TRUE, ]
    dataSetini <- dplyr::select(dataSetini, -postThinning, -postDisturbance, -preDisturbance, -YearDisturbance, -CONTROL)
    dataSetini <- dataSetini[, lapply(.SD, mean), by=c('site','src')]
    names(dataSetini)[!(names(dataSetini) %in% c('year','site','src'))] <-
        paste0(names(dataSetini)[!(names(dataSetini) %in% c('year','site','src'))], 'ini')
    dataSet <- dplyr::mutate(dataSet, Var=dataSet[, ..Nvar][[1]])
    TT <- dataSet[,.(RecMet=tryCatch(RecoveryMetrics(Var, year, preDisturbance,
        YearDisturbance, RecTime=RecTime, normalize=normalize, FormatOut='str'),
       error=function(e){return('-1/-1/-1/-1')})), by='site']
    TT <- dplyr::mutate(TT, TimeToRecover=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[, 1]),
	DegreeRecovery=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[, 2]),
	ThetaRecovery=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[, 3]),
	Impact=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[, 4]))
#	Phi1=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[,4]),
#	Phi2=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[,5]))
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
RecoveryMetrics <- function(Var,
			    Year,
			    preDisturb,
			    YearDisturbance,
			    RecTime=20,
			    normalize='baseline',
			    FormatOut='list'){
    com <- ''
    YearDisturbance <- YearDisturbance[1]
    Pd <- min(Var[preDisturb==FALSE])
    C0 <- mean(Var[preDisturb==TRUE])
    if (normalize=="baseline"){
        Var <- Var/C0
	Pd <- Pd/C0
        C0 <- 1
    }else if (normalize=="impact"){
        Var <- (Var-Pd)/(C0 - Pd)
        C0 <- 1
	Pd <- 0
    }
    Td <- Year[which(Var==min(Var[preDisturb==FALSE]))]
    Tf <- Year[which(preDisturb==FALSE & Var>=mean(Var[preDisturb==TRUE]))][1]
    T0 <- YearDisturbance[1]
    P0 <- Var[preDisturb==FALSE & Year==min(Year[preDisturb==FALSE])]
    Tx <- YearDisturbance[1] + RecTime
    Px <- Var[which(Year>=YearDisturbance[1]+RecTime)[1]]
    TimeToRecover <- Tf - T0
    DegreeRecovery <- Px / C0
    ThetaRecovery <- (C0 - Pd) / (Tf - T0)
    AdT <- TimeToRecover * C0
    if (max(Year) < (T0 + RecTime)){com <- 'RecTime is after the maximum year'; DegreeRecovery <- -99}
    if (is.na(Tf)){com <- paste0(com, '/ No recovery'); TimeToRecover <- -99; ThetaRecovery <- -99}
    if (sum(preDisturb==FALSE)<1){com <- paste0(com, '/ No data after disturbance'); DegreeRecovery <- -99; ThetaRecovery <- -99; TimeToRecover <- -99}
    if (C0 < Pd){com <- paste0(com, '/ Perturbation did not reduce the variable')}
    if (is.na(C0)){com <- paste0(com, '/ No data before perturbation'); DegreeRecovery <- -99; ThetaRecovery <- -99; TimeToRecover <- -99}
    if (max(Year) < (T0 + RecTime)){com <- 'RecTime is after the maximum year'; DegreeRecovery <- -99}
    Out <- data.frame(TimeToRecover=TimeToRecover, DegreeRecovery=DegreeRecovery,
        ThetaRecovery=ThetaRecovery, RecTime=RecTime, Impact=C0-Pd, com=com)
    if (FormatOut=='list'){
	return(Out)
    }else if (FormatOut=='str'){
	return(paste(Out$TimeToRecover, Out$DegreeRecovery, Out$ThetaRecovery, Out$Impact, sep='/'))
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
#' @param Out, str, 'HillNb' or 'Entropy' 
#' @return dataSet, data.table, formatted outuput
#' @export
format_salem <- function(dataRaw, ClassInter=10, ClassIni=7.5, Out='HillNb', type='BA'){
    listNameGrouping <- c('site', 'src')
    if (('year' %in% names(dataRaw))){listNameGrouping <- c(listNameGrouping, 'year')}
    if (('postThinning' %in% names(dataRaw))){listNameGrouping <- c(listNameGrouping, 'postThinning')}
    if (('postDisturbance' %in% names(dataRaw))){listNameGrouping <- c(listNameGrouping, 'postDisturbance')}
    if (!('src' %in% names(dataRaw))){dataRaw <- dplyr::mutate(dataRaw, src='salem')}
    dataRaw <- data.table::as.data.table(dataRaw)
    dataRaw$postThinning <- as.logical(dataRaw$postThinning)
    dataRaw$postDisturbance <- as.logical(dataRaw$postDisturbance)
    HetIndexSize <- CalcDivIndex(dataRaw, 'D_cm', ClassInter=ClassInter, ClassIni=ClassIni, type=type, Out=Out)
    names(HetIndexSize)[!(names(HetIndexSize) %in% listNameGrouping)] <- 
	   paste0(names(HetIndexSize)[!(names(HetIndexSize) %in% listNameGrouping)], 'Size')
    HetIndexSp <- CalcDivIndex(dataRaw, 'species', ClassInter=10, ClassIni=ClassIni, type=type)
    names(HetIndexSp)[!(names(HetIndexSp) %in% listNameGrouping)] <-
         paste0(names(HetIndexSp)[!(names(HetIndexSp) %in% listNameGrouping)], 'Sp')
    dataSet <- dataRaw[, .(V_m3=sum(weight*V_m3), Dg=sqrt(sum(weight*D_cm^2)/sum(weight)),
         N=sum(weight), BA=pi*sum(weight*(D_cm/200)^2)), by=listNameGrouping]
    dataSet <- merge(dataSet, HetIndexSize, by=listNameGrouping, all.x=TRUE)
    dataSet <- merge(dataSet, HetIndexSp, by=listNameGrouping, all.x=TRUE)
    dataSet <- dataSet[, BAinc:=c(diff(BA), NA), by=listNameGrouping]
    if (sum(dataSet$postDisturbance==TRUE)==0){
        dataSet <- dplyr::mutate(dataSet, YearDisturbance=NA, preDisturbance=FALSE, CONTROL=TRUE)
    }else{
        dataSet <- dataSet[, YearDisturbance:=min(year[postDisturbance==TRUE]), by=list(site, src)]
        dataSet <- dplyr::mutate(dataSet, preDisturbance=(year<=YearDisturbance & postDisturbance==FALSE),
	    CONTROL=FALSE)
    }
    class(dataSet) <- append('VirtualExperiment', class(dataSet))
    return(dataSet)
}

#' Plot virtual experiment
#'
#' This function plot a variable from formatted dataSet
#' @param dataSet, formatted dataset (VirtualExperiment object)
#' @param Nvar, string, dataSet, variable to be plotted
#' @return plot Variable chosen
#' @export
plot.VirtualExperiment <- function(dataSet, Nvar='BA', RecTime=20, normalize='baseline'){
    if (!(Nvar %in% names(dataSet))){stop('Need a variable from dataSet')}
    if (!('VirtualExperiment' %in% class(dataSet))){stop('Need format first')}
    dataSet <- dplyr::mutate(dataSet, Var=dataSet[, ..Nvar][[1]])
    dataSet <- data.table::setDT(dataSet)
    if (dataSet$CONTROL[1]==FALSE){
        dataSet <- dataSet[, C0:=mean(Var[preDisturbance==TRUE]), by=list(site, src)]
        dataSet <- dataSet[, Pd:=min(Var[preDisturbance==FALSE]), by=list(site, src)]
        if (normalize=='baseline'){
	    dataSet <- dplyr::mutate(dataSet, Var=Var/C0)
	}else if (normalize=='impact'){
	    dataSet <- dplyr::mutate(dataSet, Var= (Var-Pd)/(C0 - Pd))
	}
        Points <- dataSet[, .(C0=mean(Var[preDisturbance==TRUE]),  # Mean value of variable  before perturbation
            Pd=min(Var[preDisturbance==FALSE]),
            Td=year[which(Var==min(Var[preDisturbance==FALSE]))],  # Minimum value an time of perturbed 
            Tf=year[which(preDisturbance==FALSE & Var>=mean(Var[preDisturbance==TRUE]))][1], # First year of recovery of initial value
	    T0=YearDisturbance[1],
       	    P0=Var[preDisturbance==FALSE & year==min(year[preDisturbance==FALSE])],
	    Tx=YearDisturbance[1]+RecTime,
            Px=Var[which(year>=YearDisturbance[1]+RecTime)[1]]), by=list(site, src)]
        Area <- merge(dplyr::select(dataSet, -C0), Points, by=c("site", "src"), all.x=TRUE)
        pl <- ggplot2::ggplot(data=Points) + ggplot2::geom_hline(ggplot2::aes(yintercept=C0), linetype=2, col='blue') +
	ggplot2::geom_vline(ggplot2::aes(xintercept=Tf), linetype=2) +
	ggplot2::geom_vline(ggplot2::aes(xintercept=Tx), linetype=2) + 
	ggplot2::geom_vline(ggplot2::aes(xintercept=Td), linetype=2) +
	ggplot2::geom_text(ggplot2::aes(y=0, x=1+Tf, label='Tf'), size=4) +
	ggplot2::geom_text(ggplot2::aes(y=0, x=1+Tx, label='Tx'), size=4) +
	ggplot2::geom_text(ggplot2::aes(y=0, x=Td, label='Td'), size=4) +
	ggplot2::geom_text(ggplot2::aes(y=C0*1.2, x=1+Td, label='C0=Yd'), size=4) +
	ggplot2::geom_text(ggplot2::aes(y=Pd*0.9, x=1+Td, label='Pd'), size=4) +
	ggplot2::geom_text(ggplot2::aes(y=Px*0.9, x=1+Tx, label='Px'), size=4) +
	ggplot2::geom_point(data=dataSet, ggplot2::aes(x=year, y=Var, col=preDisturbance), size=2) +
	ggplot2::geom_line(data=dataSet, ggplot2::aes(x=year, y=Var)) +
	ggplot2::geom_ribbon(data=dplyr::filter(Area, year<=Tf, year>=T0), ggplot2::aes(x=year, ymin=C0, ymax=Var), alpha=0.7, fill='yellow') +
	ggplot2::geom_ribbon(data=dplyr::filter(Area, year>=Tf,year<=Tx), ggplot2::aes(x=year, ymin=C0, ymax=Var), alpha=0.4, fill='yellow') +
	ggplot2::geom_ribbon(data=dplyr::filter(Area, year<=Tf, year>=T0), ggplot2::aes(x=year, ymin=0*C0, ymax=C0), alpha=0.2, fill='red') +
	ggplot2::scale_color_manual(values=c('red', 'blue')) +
	ggplot2::facet_wrap(~site, scales='free', ncol=1) + ggplot2::theme_bw(base_size=16) +
	ggplot2::xlab('Years')
    }else if(dataSet$CONTROL[1]==TRUE){
        pl <- ggplot2::ggplot(data=dataSet, ggplot2::aes(x=year,y=Var, size=2)) +
	    ggplot2::geom_line() +
	    ggplot2::facet_wrap(~site, scales='free', ncol=1) + ggplot2::theme_bw(base_size=16)
    }
    print(pl)
}

buildExpControl <- function(DFexp, DFcontrol){
   DFEXP <- format_salem(DFexp)
   DFCONTROL <- format_salem(DFcontrol)
   DF <- rbind(DFEXP, DFCONTROL)
   class(dataSet) <- append('VirtualExperimentControl', class(dataSet))
}

#' Plot virtual experiment
#'
#' This function plot a variable from formatted dataSet
#' @param dataSet, formatted dataset (VirtualExperiment object)
#' @param Nvar, string, dataSet, variable to be plotted
#' @return plot Variable chosen
#' @export
plot.VirtualExperimentControl <- function(dataSet, Nvar='V_m3', RecTime=20, normalize='baseline'){
    if (!(Nvar %in% names(dataSet))){stop('Need a variable from dataSet')}
    if (!('VirtualExperimentControl' %in% class(dataSet))){stop('Need format first')}
    dataSet <- dplyr::mutate(dataSet, Var=dataSet[, ..Nvar][[1]])
    if (dataSet$CONTROL[1]==FALSE){
        dataSet <- dataSet[, C0:=mean(Var[preDisturbance==TRUE],
            Pd=min(Var[preDisturbance==FALSE])), by=list(site, src)]
        if (normalize=='baseline'){dataSet <- dplyr::mutate(dataSet, Var=Var/C0)}
        Points <- dataSet[, .(C0=mean(Var[preDisturbance==TRUE]),  # Mean value of variable  before perturbation
            Pd=min(Var[preDisturbance==FALSE]),
            Td=year[which(Var==min(Var[preDisturbance==FALSE]))],  # Minimum value an time of perturbed 
            Tf=year[which(preDisturbance==FALSE & Var>=mean(Var[preDisturbance==TRUE]))][1], # First year of recovery of initial value
	    T0=YearDisturbance[1],
       	    P0=Var[preDisturbance==FALSE & year==min(year[preDisturbance==FALSE])],
	    Tx=YearDisturbance[1]+RecTime,
            Px=Var[which(year>=YearDisturbance[1]+RecTime)[1]]), by=list(site, src)]

        Area <- merge(dplyr::select(dataSet, -C0), Points, by=c("site", "src"), all.x=TRUE)
        pl <- ggplot2::ggplot(data=Points) + ggplot2::geom_hline(ggplot2::aes(yintercept=C0), linetype=2, col='blue') +
	ggplot2::geom_vline(ggplot2::aes(xintercept=Tf), linetype=2) +
	ggplot2::geom_vline(ggplot2::aes(xintercept=Tx), linetype=2) + 
	ggplot2::geom_vline(ggplot2::aes(xintercept=Td), linetype=2) +
	ggplot2::geom_text(ggplot2::aes(y=0, x=1+Tf, label='Tf'), size=10) +
	ggplot2::geom_text(ggplot2::aes(y=0, x=1+Tx, label='Tx'), size=10) +
	ggplot2::geom_text(ggplot2::aes(y=0, x=Td, label='Td'), size=10) +
	ggplot2::geom_text(ggplot2::aes(y=C0*1.05, x=1+Td, label='Y0=Yd'), size=10) +
	ggplot2::geom_text(ggplot2::aes(y=Pd*0.95, x=1+Td, label='Pd'), size=10) +
	ggplot2::geom_text(ggplot2::aes(y=Px*0.95, x=1+Tx, label='Px'), size=10) +
	ggplot2::geom_point(data=dataSet, ggplot2::aes(x=year, y=Var, col=preDisturbance), size=2) +
	ggplot2::geom_line(data=dataSet, ggplot2::aes(x=year, y=Var)) +
	ggplot2::geom_ribbon(data=dplyr::filter(Area, year<=Tf, year>=T0), ggplot2::aes(x=year, ymin=C0, ymax=Var), alpha=0.7, fill='yellow') +
	ggplot2::geom_ribbon(data=dplyr::filter(Area, year>=Tf, year<=Tx), ggplot2::aes(x=year, ymin=C0, ymax=Var), alpha=0.4, fill='yellow') +
	ggplot2::geom_ribbon(data=dplyr::filter(Area, year<=Tf, year>=T0), ggplot2::aes(x=year, ymin=0*C0, ymax=C0), alpha=0.2, fill='red') +
	ggplot2::scale_color_manual(values=c('red', 'blue')) +
	ggplot2::facet_wrap(~site, scales='free') + ggplot2::theme_bw(base_size=16)
    }else if(dataSet$CONTROL[1]==TRUE){
        pl <- ggplot2::ggplot(data=dataSet, ggplot2::aes(x=year,y=Var, size=2)) +
	    ggplot2::geom_line() +
	    ggplot2::facet_wrap(~site, scales='free') + ggplot2::theme_bw(base_size=16)
    }
    print(pl)
}













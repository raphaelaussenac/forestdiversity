
#' Compute resilience metrics for a dataSet
#'
#' This function compute resilience metrics for a dataset 
#' @param dataSet, data.frame formatted using format_Salem or format_Samsara with multiple simulation conditions
#' @param Nvar string, name of the variable studied
#' @return A data.frame containing resilience metrics for each simulation
#' @export
EventResilience <- function(dataSet, Nvar='V_m3'){
    dataSet <- data.table(as.data.table(dataSet))
    TT <- dataSet[,.(RecMet=tryCatch(RecoveryMetrics(V_m3, year, PreDisturb, FormatOut='str'), error=function(e){return('-1/-1/-1')})), by='simulation']
    TT <- dplyr::mutate(TT, Theta=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[, 1]),
			TimeRec=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[,2]),
			DegRec=as.numeric(do.call(rbind, strsplit(RecMet, '/'))[,3]))
    TT <- dplyr::select(TT, -RecMet)
    TT[TT==-1] <- NA
    return(TT)
}

#' Compute resilience metrics
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
    VarBef <- Var[PreDisturb==TRUE] 
    if (length(VarBef)==0){stop('No value before perturbation')}
    VarPost <- Var[PreDisturb==FALSE]
    if (length(VarPost)==0){stop('No value after perturbation')}
    Intens <- VarPost[1] / VarBef
    if (Intens >= 1){warning('Perturbation did not reduce the variable')}
    dT <- Year[which(Var >= VarBef & PreDisturb==FALSE)][1] - Year[PreDisturb==TRUE]
    VardT <- Var[which(Year>=Year[PreDisturb==TRUE] + RecTime)[1]]
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
format_Samsara_Tree <- function(dataRaw){
    dataSet <- dplyr::select(dataRaw, -Int.Energy..MJ.year.1., -Pot.Energy..MJ.year.1.,
	-quality, -marketValue..euros., -dendroHabitats, -ecologicalScore, -treeId, -speciesCode)
    names(dataSet) <- c('simulation', 'year', 'eventName', 'species', 'X', 'Y', 'D_cm', 'H_m', 'V_m3')
}

#' Format Samsara dataset2
#'
#' This function format the data for further analysis
#' @param dataRaw, data.frame, dataSet containing the simulation output of Samsara
#' @return dataSet, data.table, formatted outuput
format_Samsara_Pop <- function(dataRaw){
    dataSet <- dplyr::select(dataRaw, simulationId, eventName, year, speciesName,
        N_ha, G_ha, V_ha, Dg, Dm, Carbon_AG_ha, Carbon_BG_ha, L4Cover,
	NSaplings_ha, NRecruits_ha_yr, NDead_ha_yr, cut_N_ha, cut_V_ha, deadCut_N_ha, deadCut_V_ha) %>%
        dplyr::filter(speciesName=='AllSpecies')
    names(dataSet) <- c('simulation', 'Event', 'year', 'species', 
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
#' @return dataSet, data.table, formatted outuput
#' @export
format_Salem <- function(dataRaw){
    dataRaw <- data.table::as.data.table(dataRaw)
    dataSet <- dataRaw[, .(V_m3=sum(weight*V_m3), Dg=sqrt(sum(weight*D_cm^2)/sum(weight))), by=list(site,year)]
    dataSet <- dataSet[, PreDisturb:=(year==min(year)), by='site']
    names(dataSet)[names(dataSet)=="site"] <- 'simulation'
    return(dataSet)
}

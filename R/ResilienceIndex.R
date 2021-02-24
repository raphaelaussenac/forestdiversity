## Format : year, variable name 


EventResilience <- function(dataSet, Nvar='V_m3'){
    dataRaw <- read.csv('DATA/VirtualExperiment/outputDistributions.txt',sep='\t')
    A <- format_SamsaraPop(dataRaw)
    TTR <- dplyr::group_by(A, simulation) %>% dplyr::summarise(RecMet=
        tryCatch(RecoveryMetrics(V_m3, year, PreDisturb), error=function(e){return(NA)})) %>%
        dplyr::ungroup()
    TTR <- cbind(TTR[, 1], TTR$RecMet)
    return(TTR)
}

RecoveryMetrics <- function(Var, Year, PreDisturb, RecTime=20){
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
    return(Out)
}

format_Samsara_Tree <- function(dataRaw){
    dataSet <- dplyr::select(dataRaw, -Int.Energy..MJ.year.1., -Pot.Energy..MJ.year.1.,
	-quality, -marketValue..euros., -dendroHabitats, -ecologicalScore, -treeId, -speciesCode)
    names(dataSet) <- c('simulation', 'year', 'eventName', 'species', 'X', 'Y', 'D_cm', 'H_m', 'V_m3')
}

format_Samsara_Pop <- function(dataRaw){
    dataSet <- dplyr::select(dataRaw, simulationId, eventName, year, speciesName,
        N_ha, G_ha, V_ha, Dg, Dm, Carbon_AG_ha, Carbon_BG_ha, L4Cover,
	NSaplings_ha, NRecruits_ha_yr, NDead_ha_yr, cut_N_ha, cut_V_ha, deadCut_N_ha, deadCut_V_ha) %>%
        dplyr::filter(speciesName=='AllSpecies')
    names(dataSet) <- c('simulation', 'Event', 'year', 'species', 
        'N', 'G', 'V_m3', 'Dg', 'Dm', 'Carbon_AG_ha', 'Carbon_BG_ha',
	'L4cover', 'NSaplings_ha', 'NRecruits_ha_yr' , 'NDead_ha_yr',
       	'cut_N_ha', 'cut_V_ha', 'deadCut_N_ha', 'deadCut_V_ha')
    dataSet <- filter(dataSet, Event!='SaplingCreation')
    dataSet <- dplyr::mutate(dataSet, EventPerturb=substr(Event, 1, 7)=="Disturb",
	PreDisturb=substr(Event, 1, 8)=="Creation")
    return(dataSet)   
}

format_Salem <- function(dataRaw){
    dataSet <- dplyr::group_by(dataRaw, site, year) %>%
       dplyr::summarise(V_m3=sum(weight*V_m3), Dg=sqrt(sum(weight*D_cm^2)/sum(weight))) %>%
       dplyr::group_by(site) %>% dplyr::mutate(PreDisturb=(year==min(year))) %>% dplyr::ungroup()
    names(dataSet)[names(dataSet)=="site"] <- 'simulation'
    return(dataSet)
}

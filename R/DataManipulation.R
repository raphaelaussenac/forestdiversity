GatherDataSim <- function(evalSite, listsrc = c('data', 'landclim', '4c', 'salem')){
# retrieve list of file (simulations and observations)
  fileNames <- Sys.glob(file.path('DATA','Raw',evalSite,'*.csv'))
  if(evalSite == 'profound'){
    listsite <- c('kroof', 'solling-beech', 'solling-spruce')
  }
  alldf <- data.frame()
  for (src in listsrc){
    for (i in grep(pattern = src, x=tolower(basename(fileNames)))){
      temp <- read.csv(fileNames[i])
      if(evalSite == 'profound'){
        site <- names(which(sapply(listsite,grepl, x=basename(fileNames[i]))==TRUE))
        temp <- cbind(site, temp)
        temp$weight <- 1
      }
      temp$src <- src
      alldf <- rbind(alldf, temp)
    }
  }
  out <- write.csv(alldf, paste0('./data/all_', evalSite, '.csv'), row.names = FALSE)
  return(file.exists(paste0('./data/all_', evalSite, '.csv')))
}


#' Calculate stand population metrics'
#' This function takes a data.frame and return a Plot object. The data.frame
#' must contain variables: species, D_cm, weight, X, Y 
#'
#' @param evalSite str, name of the site
#' @return data.frame with # calculate yearly aggregated index (N, Dg, BA, BAI, etc) at sp level and stand level (species=="allsp")
#' @export
standVarCalc <- function(evalSite){
  alldf <- read.csv(file.path('DATA',paste0('all_', evalSite, '.csv')))
  YearlySP <- group_by(alldf, year, src, site, species) %>% 
	dplyr::summarise(N = sum(weight),
		  Dg = sqrt(sum(D_cm^2 * weight)/sum(weight)),
		  H = sum(H_m * weight) /sum(weight),
		  BA = sum((pi * (D_cm/200)^2) * weight),
		  V = sum(V_m3 * weight)) %>% ungroup()
  YearlySt <- group_by(alldf, year, src, site) %>% 
	dplyr::summarise(N = sum(weight),
		  Dg = sqrt(sum(D_cm^2 * weight)/sum(weight)),
		  H = sum(H_m * weight) /sum(weight),
		  BA = sum((pi * (D_cm/200)^2) * weight),
		  V = sum(V_m3 * weight)) %>% ungroup() %>% mutate(species='allsp')
  Yearly <- rbind(YearlySP, YearlySt)
  YearlyObs <- group_by(Yearly, site) %>% dplyr::filter(year %in% unique(year[src=='data'])) %>% ungroup()
  YearlyNonObs <- group_by(Yearly, site) %>% dplyr::filter(!(year %in% unique(year[src=='data']))) %>% ungroup()
  yb <- function(year){
 	if (length(year)>1){
		return(c(NA,year[1:(length(year)-1)]))
	}else{
		return(NA)
	}
  }
  DecYear <- 0
  if (evalSite=='bauges'){DecYear <- 1}
  YearlyObs <- group_by(YearlyObs, site, src, species) %>% 
	dplyr::mutate(yearbefore=yb(year), BAbefore=yb(BA), BAI_yr=(BA-BAbefore)/(year-yearbefore+DecYear)) %>%
       	ungroup() %>% dplyr::select(-yearbefore, -BAbefore)
  return(rbind(YearlyObs, mutate(YearlyNonObs, BAI_yr=NA)))
}




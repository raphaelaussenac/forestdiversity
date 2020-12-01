## Package to compute forest diversity index in IMaestro project

Outputs of model/data must be matrix with a line for each year/species/DBH.

Variables recorded must be at least year, species, DBH and a weight.


### Population (=stand) level metrics :
DBH-Distribution related : 
        Nb of DBH class
        Hill Numbers (Shannon, Simpson) based on DBH Class
        Gini index (species only)
        Skewness of DBH distribution
        Gini of DBH distribution

Species-related :
        Nb of species
        Hill Numbers (Shannon, Simpson) based on species


Metrics are calculated for each year, each site and each data source

### Plot level metrics: spatial metrics can be computed if each tree location are known
DBH-Distribution related : 
        DBH segregation index
        Uniform angle index 
        Structural Complexity index

Species-related :
        Species mingling index

# Functions to compute diversity 
Need first to compute a TabDist object using create_TabDist function using the data (data.frame), the shape of the plot and its coordinates
Metrics can then be computed using the TabDist obejct
Refs for the metrics and correction : 

Pommerening and UriaDiez, 2017 : Do large forest trees tend towards high species mingling ?
Hui and Gadow, 2002 : Characterizing forest spatial structure and diversity
Zenner and Hibbs, 2000 : A new method for modeling the heterogeneity of forest structure
Pommerening and Stoyan, 2006 : Edge-correction needs in estimation indices of spatial forest structure

### Landscape level metrics 

### Resilience metrics 

One event case

	Resistance
	Recovery rate
	Recovery time
	Degree of recovery at time tx=20
	Efficiency at time tx

Regime of disturbances

	Temporal instability
	Temporal autocorrelation
	Non-permanent time
	Non-permanent intensity
	Net change
	


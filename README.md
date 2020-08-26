## Package to compute forest diversity index in IMaestro project

Outputs of model/data must be matrix with a line for each year/species/DBH.


Variables recorded must be at least year, species, DBH and a weight.


### Model ES outputs:

	Timber volume harvested by species and dbh class
	Current annual volume increment per species and dbh class
	Above ground carbon
	Below ground carbon
	Dead wood carbon
	Soil carbon

Biodiversity

	Dead wood volume
	Abundance of large standing dead trees
	Dead wood diversity
	Number of trees with dendromicrohabitats
	Canopy cover
	

Disturbance
	Volume damaged by disturbances per species and dbh class
	



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

### Landscape level metrics *:

### Resilience metrics *:

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
	


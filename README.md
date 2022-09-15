# Bat One Health Functions for the Multi Scale Modeling Effort
A set of functions related to the multi scale modeling work for Bat One Health. These set of functions allow for the attachment of environmental data to coordinates that are used in all of the component models. These environmental data can be currently found [here](https://drive.google.com/drive/folders/1cfwvPG9wID0MgaP332Dt2KXR_zsZXcRR?usp=sharing).

### getAustraliaEnv
A helper function for the attachAustraliaEnv function that is used to access raster data for a particular environmental variable matching characteristics given (lag time, year, month).

### attachAustraliaEnv
A function used to attach environmental data to a given set of points for use in further modeling processes. Requires a two column data.frame of coordinates in the GDA 2020 coordinate reference system (or an sf points object with the same CRS) along with vectors of years and months related to those points.

### estimatePrevalence
This is the main function used to estimate prevalence within our study area. Given a vector of years and months, this function will estimate prevalence within the area based on stress. If wanted, the user can modify the number of repetitions for each year-month calculation of prevalence, change the reference max prevalence (currently based on estimates from field data from [Field et al. 2015](https://doi.org/10.1371/journal.pone.0144055)), and parallelize the calculations. 

### \data

### locations_sites.csv 
csv file of roost location estimates from Hume et al 2015 PLoS One paper

### hendra-virus-transmission-data-east-australia.csv
tested urine pools; downloaded from: https://www.data.qld.gov.au/dataset/hev-infection-flying-foxes-eastern-australia

### To do list
* Fix whatever else breaks

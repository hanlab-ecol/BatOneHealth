# Bat One Health Functions for the Multi Scale Modeling Effort
A set of functions related to the modeling work for Bat One Health.

### getAustraliaEnv
A helper function for the attachAustraliaEnv function that is used to access raster data for a particular environmental variable matching characteristics given (lag time, year, month).

### attachAustraliaEnv
A function used to attach environmental data to a given set of points for use in further modeling processes. Requires a two column data.frame of coordinates in the GDA 2020 coordinate reference system along with vectors of years and months related to those points.

### To do list
* Allow attachAustraliaEnv to work with sf objects and with different CRSs
* Add in other indices to attachAustraliaEnv
* Put error messages regarding month and year arguments, allow for the use of a single year and month value
* Put in ability to choose variables to be log transformed (maybe another helper function)
* Fix whatever else breaks

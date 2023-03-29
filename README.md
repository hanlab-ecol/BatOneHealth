# Bat One Health Functions for the Multi Scale Modeling Effort
Here we present a series of functions and code for facilitating the recreation of our multi-scale modeling work as part of Bat One Health. The functions described below enable the estimation of monthly Hendra virus prevalence within a region of eastern Australia from the years 2008 to 2019 at a spatial resolution of 5 km. Most of the functions are a series of helper functions used to allow for the attachment of environmental data to coordinates that are used in all of the component models. These environmental data can be currently found in a FigShare repository [here](https://figshare.com/s/ddb5a1584609b20f6596).

The main script you will be working with to run this model is `MultiScaleExample.Rmd` which is mainly concerned with small bits of setup for the running of the `estimatePrevalence` function.

The boundaries of our study area are shown in the polygon designated below.

```geojson
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "id": 1,
      "properties": {
        "ID": 0
      },
      "geometry": {
        "type": "Polygon",
        "coordinates": [
          [
              [138.23433,-38.36296],
              [138.23433,-15.27318],
              [153.63063,-15.27318],
              [153.63063,-38.36296],
              [138.23433,-38.36296]
          ]
        ]
      }
    }
  ]
}
```

## Setup
Environmental data should be downloaded from the link above, unzipped, and placed in a folder named AustralianClimateData that is to be placed in the working directory. You should download 15 total items from the FigShare link (bat roost predictions, maximum temperature, minimum temperature, temperature difference, precipitation, precipitation anomaly, NDVI, soil moisture, morning vapor pressure, solar exposure, potential evapotranspiration, land cover, ONI, SOI, and SAM). The folder structure should be as follows:
>Your_Working_Directory<br>
>\NewRoostPredictedSurfacesModel.Rdata<br>
>\RehabModel_01062022.Rdata<br>
>\shortage_reduced_environmental.txt<br>
>\AustralianClimateData<br>
>\\\BatRoostPredictions<br>
>\\\Evapotranspiration<br>
>\\\PrecipitationLag<br>
>\\\\\PrecipitationLag2Months<br>
>\\\\\PrecipitationLag3Months<br>
>\\\\\PrecipitationLag6Months<br>
>...<br>

Make sure that the two model objects and the food shortage model results found in the data folder of this repository are in your working directory as well as the folder containing all of the environmental data. This setup will help to prevent some troubleshooting problems. If you download all of the items in the FigShare and place them in a containing folder (`AustralianClimateData` shown above), the structure described above should be retained.

## Expected Output
The output of the main estimatePrevalence function will be a series of rasters estimating Hendra virus prevalence for all months in the specified years. Depending on the model components chosen (see estimatePrevalence function and MultiScaleExample for more information), you will output rasters for estimated prevalence, rehab stress, new roost stress, and a null expectation of stress derived from roost environmental suitability. The output will always include estimated prevalence and the null expectation, with new roost and rehab stress rasters only created when those model components are included.

Note that function will not produce an object in R. If there are no error messages produced after the function runs, you should see the aforementioned raster files created in your working directory.

## Descriptions
The following are descriptions of each of the items included in this repository and how they work together:

### MultiScaleExample
This is an example script showing all of the packages, how to source the functions and data mentioned further below, and how to implement the estimatePrevalence function to generate Hendra virus prevalence predictions for our study area. This file is commented throughout on what each code chunk is doing and further information on each particular function can be found on their respective pages in this repository.

### getAustraliaEnv
A helper function for the attachAustraliaEnv function that is used to access raster data for a particular environmental variable matching characteristics given (lag time, year, month).

### attachAustraliaEnv
A function used to attach environmental data to a given set of points for use in further modeling processes. Requires a two column data.frame of coordinates in the GDA 2020 coordinate reference system (or an sf points object with the same CRS) along with vectors of years and months related to those points.

### roost.counts
A function that uses a quarterly survey of roost counts with a generalized additive model (GAM) to produce estimates of number of roosts present on a quarterly basis within the time frame of our study. This allows us to account for shifts in the number of roosts based on seasonal changes rather than assuming a static number of roosts at all times. 

### estimatePrevalence
This is the main function used to estimate prevalence within our study area. Given a vector of years and months, this function will estimate prevalence within the area based on stress. If wanted, the user can modify the number of repetitions for each year-month calculation of prevalence, change the reference max prevalence (currently based on estimates from field data from [Field et al. 2015](https://doi.org/10.1371/journal.pone.0144055)), and parallelize the calculations. 

## Data descriptions
These are descriptions of all of the objects placed in \data

### locations_sites.csv 
A CSV file of roost location estimates from [Field et al. 2015 PLoS One paper](https://doi.org/10.1371/journal.pone.0144055).

### hendra-virus-transmission-data-east-australia.csv
Hendra virus prevalence data from tested urine pools by [Field et al. 2015](https://doi.org/10.1371/journal.pone.0144055). These data are accessible from the Queensland government's [Open Data Portal](https://www.data.qld.gov.au/dataset/hev-infection-flying-foxes-eastern-australia).

### bff_one_week_per_quarter_ilya.csv
This is a CSV of black flying fox roost counts summarized from the [Flying Fox Monitoring Program](https://www.data.qld.gov.au/dataset/flying-fox-monitoring-program) by quarter of the year. This dataset is used by the roost.counts function above to estimate the number of roosts for a given month. 

### shortage_reduced_environmental.txt
This is a text file containing the monthly results from a model predicting food shortages (Lagergren and Ruiz-Aravena et al. in prep).

### NewRoostPredictedSurfacesModel
This is a gbm model object representing a model trained on roosts that were determined to be newly formed. This model is used to predict stress deriving from the formation of new roosts when used with data from each month queried by the estimatePrevalence function.

### RehabModel_01062022
This is also a gbm model object representing a model trained on bat rehabilitation data. This model is used to predict stress that would lead to bats ending up at a rehabilitation center (this is not necessarily limited to food stress) when used with data from each month queried by the estimatePrevalence function.

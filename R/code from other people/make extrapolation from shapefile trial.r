library(rgdal)
library(sf)
library(sp)

OWF_area_shapefile <- readOGR(dsn = "data/BOEM_Renewable_Energy_Shapefiles_4_20230830/", layer = "BOEM_Wind_Leases_8_30_2023") # load the shape file
OWF_area_shapefile_sf <- st_as_sf(OWF_area_shapefile) # convert to a sf object for convenience

### Crop the shape file to only keep the mid-Atlantic region
bounding_box <- st_bbox(c(xmin = -76, xmax = -70, ymin = 36, ymax = 41.5), crs = st_crs(OWF_area_shapefile_sf))
OWF_area_shapefile_sf <- st_crop(OWF_area_shapefile_sf, bounding_box)
remove(bounding_box)

### a quick visualization to check
# ggplot() +
#   geom_sf(data = OWF_area_shapefile_sf) +
#   theme_minimal() 

### if all good, convert the sf back to sp file
OWF_area_shapefile <- as(OWF_area_shapefile_sf, "Spatial"); remove(OWF_area_shapefile_sf)
class(OWF_area_shapefile)

### calculate UTM zone based on average longitude 
lon <- sum(bbox(OWF_area_shapefile)[1,])/2 
utmzone <- floor((lon + 180)/6)+1; remove(lon) 

### convert polygon to UTM
# test <- spTransform(OWF_area_shapefile, CRS("+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")) # from Thorson script
OWF_area_shapefile_UTM <- spTransform(OWF_area_shapefile, CRS(paste0("+proj=utm +zone=", utmzone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))); remove(utmzone)



### Construct the extrapolation grid for VAST using sf package
# specify the grid size in meters (since working in UTM) to control the grid resolution. This step is slow at high resolutions 
region_grid <- st_make_grid(OWF_area_shapefile_UTM, cellsize = 2000, what = "centers")  



make_settings( n_x = 100, 
               Region = example$Region, 
               
               purpose = "index2", 
               bias.correct = FALSE )





## Convert region_grid to Spatial Points to SpatialPointsDataFrame
region_grid <- as(region_grid, "Spatial")
region_grid_sp <- as(region_grid, "SpatialPointsDataFrame")
## combine shapefile data (region_polygon) with Spatial Points
## (region_grid_spatial) & place in SpatialPointsDataFrame data
## (this provides you with your strata identifier (here called
## Id) in your data frame))
region_grid_sp@data <- over(region_grid, region_polygon)
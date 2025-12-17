library(sf)
library(rgdal)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(mapdata)

### --------------------------------------------------- ###
# rgdal approach
### --------------------------------------------------- ###

my_spdf <- readOGR( 
  dsn = "data/BOEM_Renewable_Energy_Shapefiles_4_20230830/", 
  layer = "BOEM_Wind_Leases_8_30_2023",
  verbose=FALSE
)

summary(my_spdf)
length(my_spdf)
head(my_spdf)

par(mar=c(0,0,0,0))
plot(my_spdf)

# 'fortify' the data to get a dataframe format required by ggplot2
library(broom)
spdf_fortified <- tidy(my_spdf)

USA <- ne_countries(scale = "medium", returnclass = "sf", country = 'United States of America')

ggplot(data = USA) +
  geom_sf(fill= "cornsilk") +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-65, -80), ylim = c(33, 42.5)) +
  geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), colour = "black", fill = NA) +
  theme(panel.grid.major = element_line(color = gray(.85), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"), legend.position = "none")


USA.state <- map_data('state')
centroid_labels <- usmapdata::centroid_labels("states")

ggplot(data = USA.state) +
  geom_polygon(aes(x = long, y = lat, group = group), fill= "cornsilk", color = "black") +
  # coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-68, -78.2), ylim = c(33, 42.5)) +
  geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), color = "black", fill = "black") +
  theme(panel.grid.major = element_line(color = gray(.85), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"), legend.position = "none")  
  



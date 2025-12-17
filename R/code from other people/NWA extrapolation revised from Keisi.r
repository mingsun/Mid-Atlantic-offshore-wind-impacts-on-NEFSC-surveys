



utils::data(northwest_atlantic_grid, package = "SpatialDeltaGLMM")
Data_Extrap <- northwest_atlantic_grid
Area_km2_x = Data_Extrap[, "Area_in_survey_km2"]


utils::data(northwest_atlantic, package = "VAST")
load_example(northwest_atlantic)

test <- read.csv("data/northwest_atlantic_grid.csv")

ggplot(test) +
  geom_point(aes(x = Lon, y = Lat, color = EPU)) +
  # geom_point(data = full_VAST_df, aes(x = LON, y = LAT))
  geom_point(data = subset(full_df, OWF == "INSIDE"), aes(x = LON, y = LAT))


test <- load_example(data_set = "BC_pacific_cod")
test$strata.limits
remove(test)

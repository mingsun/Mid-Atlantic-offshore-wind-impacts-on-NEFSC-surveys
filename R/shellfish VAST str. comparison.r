library(tidyverse)
library(visreg)


# 0. data set-up ---------------------------------------------------------------------------------------------

catch_by_tow_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") %>%
  select(ID, REGION)

full_df <- read.csv("results/stratified.mean.indices/surfclam/full.info.tow.list.csv")[,-1] %>%
  filter(OWF == "OUTSIDE") %>%
  left_join(catch_by_tow_df)

remove(catch_by_tow_df)


# 5. VAST ---------------------------------------------------------------------------------------------

library(VAST)
library(splines)
library(Metrics)

## prepare data ----
full_VAST_df <- full_df %>%  
  mutate(AreaSwept_km2 = ASWEPT_M2/1000000) %>%
  add_column(Vessel = NA, Pass = 0) %>% # required by the model, don't know why  
  rename(Lat = LAT, Lon = LON) %>% # the exact spelling is required by VAST
  mutate(Year = as.integer(factor(YEAR)))

## prepare settings and covariate & formula ----

### create the setting file
settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[2], # number of knots, the more knots the longer it runs
                           purpose = "index2", # to calculate index
                           Region = "user", 
                           # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
                           # FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0), 
                           FieldConfig = c("Omega1" = 0, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 1),
                           # specify various options for turning on/off spatial and spatial-temporal variation factors in the two linear predictors
                           # omega: spatial variation ; epsilon: spatiotemporal variation
                           # details see here: https://rdrr.io/github/James-Thorson/VAST/man/make_data.html
                           RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0), 
                           # specify whether either intercepts (Beta1 and Beta2) or spatio-temporal variation (Epsilon1 and Epsilon2) is structured among time intervals
                           # untouched default setting: all is zero, so all are fixed effect
                           OverdispersionConfig= c(0,0),
                           # vessel/targeting effects for encounter probability and positive catch rate
                           # untouched default setting (0,0) means no effects
                           ObsModel = c(2,1),
                           # link functions for encounter probability and positive catch rate are available
                           # untouched default setting here
                           bias.correct = TRUE # to save time set as false, epsilon bias-correction
)

full_VAST_df$TEMP_Std <- scale(full_VAST_df$TEMP, center = TRUE, scale = TRUE)[,1]
full_VAST_df$DEPTH_Std <- scale(full_VAST_df$DEPTH, center = TRUE, scale = TRUE)[,1]

### design formulas for the two variables
p1_formula = ~ bs(TEMP_Std, degree = 2) + bs(DEPTH_Std, degree = 2)
p2_formula = ~ bs(TEMP_Std, degree = 2) + bs(DEPTH_Std, degree = 2)



## SVAtoSNE dataset model ----
SVAtoSNE_VAST_df <- subset(full_VAST_df, REGION == "SVAtoSNE") %>%
  mutate(Year = as.numeric(factor(YEAR)))

### standardize the covariate to have mean 0 and standard deviation 1.0 as suggested by VAST
SVAtoSNE_VAST_df$TEMP_Std <- scale(SVAtoSNE_VAST_df$TEMP, center = TRUE, scale = TRUE)[,1]
SVAtoSNE_VAST_df$DEPTH_Std <- scale(SVAtoSNE_VAST_df$DEPTH, center = TRUE, scale = TRUE)[,1]

### loop over multiple error/link options to identify the best model ----

scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(3,4), AIC = NA, RMSE = NA) ##details see "make_data" description
i = 3

settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])

## prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(3)],
                                                DirPath = "results/temporal model/surfclam/models objects/VAST/")

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"



# VAST_temp_SVAtoSNE <- fit_model(settings = settings,
#                                 input_grid = extrapolation_region$Data_Extrap,
#                                 Lat_i = SVAtoSNE_VAST_df$Lat,
#                                 Lon_i = SVAtoSNE_VAST_df$Lon,
#                                 t_i = as.numeric(SVAtoSNE_VAST_df$Year),
#                                 b_i = SVAtoSNE_VAST_df$BIOMASS,
#                                 a_i = as_units(SVAtoSNE_VAST_df$AreaSwept_km2, "km^2"),
#                                 # covariate_data = SVAtoSNE_VAST_df, # use full_VAST_df here to avoid missing year (2020) for GBK
#                                 # X1_formula = p1_formula,
#                                 # X2_formula = p2_formula,
#                                 working_dir = "results/temporal model/surfclam/models objects/VAST/")

if (is.null(VAST_temp_SVAtoSNE)) { next }
if (length(VAST_temp_SVAtoSNE$Report) == 1) {
  if (VAST_temp_SVAtoSNE$Report == "Model is not converged") { next }
}

# save(VAST_temp_SVAtoSNE, file = paste0("results/temporal model/surfclam/models objects/VAST_temp_SVAtoSNE_",
#                                        scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))



load("results/temporal model/surfclam/models objects/VAST_temp_SVAtoSNE.Rdata")


VAST_temp_SVAtoSNE$data_frame$b_i <- strip_units(VAST_temp_SVAtoSNE$data_frame$b_i) 
VAST_temp_SVAtoSNE$data_frame$a_i <- strip_units(VAST_temp_SVAtoSNE$data_frame$a_i)

VAST_comp_full_dataset_df <- full_VAST_df %>%
  select(OWF, Year, YEAR, Lat, Lon, TEMP, DEPTH, AreaSwept_km2, BIOMASS) %>%
  # mutate(Year = as.integer(Year)) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "VAST", DATASET = "FULL") %>%
  mutate(PREDICTED_BIOMASS =  predict(x = VAST_temp_SVAtoSNE,
                                      what = "D_i",
                                      Lat_i = .$Lat,
                                      Lon_i = .$Lon,
                                      t_i = .$Year,
                                      a_i = .$AreaSwept_km2,
                                      do_checks = FALSE))


class(VAST_comp_full_dataset_df$YEAR)
order(unique(VAST_comp_full_dataset_df$Year))

scenario_df$AIC[i] = VAST_temp_SVAtoSNE$parameter_estimates$AIC # -4544.321 vs -4642.019



library(tidyverse)

# 0. data set-up ---------------------------------------------------------------------------------------------

 ## prepare spatial and environmental data, tow level catch data are also here

catch_by_tow_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") %>%
  select(ID, REGION)

full_df <- read.csv("results/stratified.mean.indices/surfclam/tow.list.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>%
  left_join(catch_by_tow_df) 

# full_df <- filter(full_df, REGION == "SVAtoSNE") # only work with the Southern region as the North is not impacted 

## prepare the OWF overlay info

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

full_df <- full_df %>% 
  mutate(OWF = ifelse(ID.temp %in% unique(AS_QQ_overlay_df$ID.temp), "INSIDE", "OUTSIDE")) %>% 
  select(-c(ID.temp)) %>%
  rename(BIOMASS = KGPERTOW) %>%
  filter(!is.na(TEMP), !is.na(LON))

sum(full_df$OWF == "INSIDE") 
sum(full_df$OWF == "OUTSIDE") 
  
remove(AS_QQ_overlay_df)

 ## look at the distribution of NUMBER
hist(full_df$BIOMASS)
sum(full_df$BIOMASS == 0)/nrow(full_df) # 0.2568729

write.csv(full_df, "results/stratified.mean.indices/surfclam/full.info.tow.list.csv")
remove(catch_by_tow_df)

# -------------------------------------------------------------------------------------------------------- #




# 0. identify variable weighting for SDM ---------------------------------------------------------------------------------------------

 ## BRT to calculate variable weighting
# library(dismo)
# library(gbm)

# variable selection following Kleisner et al., 2017 Progress in Oceanography & 2016 PlosOne

 ## determine n.trees with gbm.step
# BRT_0 <-  gbm.step(data = full_df, gbm.x = c(8:12), gbm.y = 16,  
#                    family = "poisson", tree.complexity = 5, learning.rate = 0.03, bag.fraction = 0.75)

BRT_0 <- gbm(BIOMASS ~ factor(YEAR) + factor(SEASON) + LAT + LON + SURFTEMP + TEMP + DEPTH, data = full_df, 
                 distribution = "poisson", # non-negative counts
                 n.trees = 1000, bag.fraction = 0.75) # training data proportion

summary(BRT_0)
gbm.perf(BRT_0) # aiming for lower error  

BRT_final <- gbm(BIOMASS ~ factor(YEAR) + LAT + TEMP + DEPTH, data = full_df, 
             distribution = "poisson", 
             n.trees = 1000, bag.fraction = 0.75)
summary(BRT_final)
gbm.perf(BRT_final) # based on Residual Standard Error (RSE)

 ## Correlation analyses

 ## select correlation between variables and go back to last step
  ### previous SDM works selected YEAR + SEASON + s(SURFTEMP) + s(TEMP) + s(DEPTH) (Kleisner et al., 2017 Progress in Oceanography & 2016 PlosOne)
  ### KK also used delta lognormal distribution, we followed it here
library("PerformanceAnalytics")

cor.p <- chart.Correlation(full_df[,c("YEAR", "LAT", "TEMP", "DEPTH")], histogram = TRUE, pch = 19)

## final predictor list: factor(YEAR) + LAT + TEMP + DEPTH

remove(BRT_0, BRT_final)

# -------------------------------------------------------------------------------------------------------- #




# 1. Tweedie GAM ---------------------------------------------------------------------------------------------

library(mgcv)

full_WG_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))

 ## full dataset model ----

TW.GAM_dist_full <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = tw(), data = full_WG_df)

summary(TW.GAM_dist_full)
AIC(TW.GAM_dist_full)

 ### predict the results and save
TW.GAM_comp_full_dataset_df <- full_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Tweedie-GAM", DATASET = "FULL") %>%
  mutate(PREDICTED_BIOMASS =  predict(TW.GAM_dist_full, newdata = ., type = "response"))
  
## outside dataset model (transfer ability) ----

TW.GAM_dist_outside <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = tw(), data = subset(full_WG_df, OWF == "OUTSIDE"))

summary(TW.GAM_dist_outside)
AIC(TW.GAM_dist_outside)

### predict the results and save
TW.GAM_comp_outside_dataset_df <- full_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Tweedie-GAM", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_BIOMASS =  predict(TW.GAM_dist_outside, newdata = ., type = "response"))

TW.GAM_comp <- rbind(TW.GAM_comp_full_dataset_df, TW.GAM_comp_outside_dataset_df)
write.csv(TW.GAM_comp, "results/spatial model/surfclam/Tweedie.GAM_abundance_comparison.csv", row.names = FALSE)

save(TW.GAM_dist_full, full_WG_df, file = "results/spatial model/surfclam/models objects/TW.GAM_dist_full.Rdata")
save(TW.GAM_dist_outside, full_WG_df, file = "results/spatial model/surfclam/models objects/TW.GAM_dist_outside.Rdata")

rm(list = setdiff(ls(), c("full_df", "data_partition_list")))

# -------------------------------------------------------------------------------------------------------- #



# 2. delta GAM ---------------------------------------------------------------------------------------------

library(mgcv)

 # Grüss et al. 2014 FR, delta  approach is a two stage GAM
 # stage 1: a GAM to model the probability of presence using a binomial GAM with logit link
 # stage 2: a GAM to model abundance if present using a quasi-Poisson GAM with log link
 # and then multiplying the results together to obtain an overall abundance, AIC are added for total AIC

## full dataset model ----
full_DELTA.GAM_df <- full_df %>%
  mutate(PRESENCE = ifelse(BIOMASS > 0, 1, 0)) %>%
  mutate(YEAR = as.factor(YEAR))

DELTA.S1.GAM_dist_full <- gam(PRESENCE ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = binomial(link = "logit"), data = full_DELTA.GAM_df)
DELTA.S2.GAM_dist_full <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = tw(), data = subset(full_DELTA.GAM_df, PRESENCE == 1))

summary(DELTA.S1.GAM_dist_full)
summary(DELTA.S2.GAM_dist_full)
AIC(DELTA.S1.GAM_dist_full)
AIC(DELTA.S2.GAM_dist_full)

### predict the results and save
DELTA.GAM_comp_full_dataset_df <- full_DELTA.GAM_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, BIOMASS, PRESENCE) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Delta-GAM", DATASET = "FULL") %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_dist_full, newdata = ., type = "response"),
         PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_dist_full, newdata = ., type = "response"))

## outside dataset model (transfer ability) ----

DELTA.S1.GAM_dist_outside <- gam(PRESENCE ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, OWF == "OUTSIDE"))
DELTA.S2.GAM_dist_outside <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = tw(), data = subset(full_DELTA.GAM_df, OWF == "OUTSIDE" & PRESENCE == 1))

summary(DELTA.S1.GAM_dist_outside)
summary(DELTA.S2.GAM_dist_outside)
AIC(DELTA.S1.GAM_dist_outside)
AIC(DELTA.S2.GAM_dist_outside)

### predict the results and save
DELTA.GAM_comp_outside_dataset_df <- full_DELTA.GAM_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, BIOMASS, PRESENCE) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Delta-GAM", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_dist_outside, newdata = ., type = "response"),
         PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_dist_outside, newdata = ., type = "response"))

DELTA.GAM_comp <- rbind(DELTA.GAM_comp_full_dataset_df, DELTA.GAM_comp_outside_dataset_df)
write.csv(DELTA.GAM_comp, "results/spatial model/surfclam/Delta.GAM_abundance_comparison.csv", row.names = FALSE)

save(DELTA.S1.GAM_dist_full, full_DELTA.GAM_df, file = "results/spatial model/surfclam/models objects/DELTA.S1.GAM_dist_full.Rdata")
save(DELTA.S2.GAM_dist_full, full_DELTA.GAM_df, file = "results/spatial model/surfclam/models objects/DELTA.S2.GAM_dist_full.Rdata")
save(DELTA.S1.GAM_dist_outside, full_DELTA.GAM_df, file = "results/spatial model/surfclam/models objects/DELTA.S1.GAM_dist_outside.Rdata")
save(DELTA.S2.GAM_dist_outside, full_DELTA.GAM_df, file = "results/spatial model/surfclam/models objects/DELTA.S2.GAM_dist_outside.Rdata")

rm(list = setdiff(ls(), c("full_df", "data_partition_list")))

# -------------------------------------------------------------------------------------------------------- #


# 3. RANDOM FOREST ---------------------------------------------------------------------------------------------

library(randomForest)

full_RF_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))

 ## determine the number of trees to grow (n.tree) and number of variables runadomly samples at each split (mtry)
set.seed(1) # for reproducibility
train_index <- sample(1:nrow(full_df), 0.8 * nrow(full_df)) # 80% for training
train_data <- full_df[train_index, ]
valid_data <- full_df[-train_index, ]

  ### Perform cross-validation to determine optimal mtry
rf_tune <- tuneRF(x = train_data[,c("YEAR", "LAT", "TEMP", "DEPTH")],
                  y = train_data$BIOMASS,
                  ntree = 1000, # change n.tree manually 1000, 1500
                  mtryStart = 2,
                  stepFactor = 1, improve = 0.01, trace = TRUE, plot = TRUE)

remove(train_index, train_data, valid_data, rf_tune)

## full dataset model ----
 
set.seed(2) # ensure reproducibility
RF_dist_full <- randomForest(BIOMASS ~ YEAR + LAT + TEMP + DEPTH, data = full_RF_df, mtry = 2, ntree = 1000)

print(RF_dist_full) # Print the model summary

### predict the results and save
RF_comp_full_dataset_df <- full_RF_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Random Forest",DATASET = "FULL") %>%
  mutate(PREDICTED_BIOMASS =  predict(RF_dist_full, newdata = ., type = "response"))


## outside dataset model (transfer ability) ----

set.seed(2) # ensure reproducibility
RF_dist_outside <- randomForest(BIOMASS ~ YEAR + LAT + TEMP + DEPTH, data = subset(full_RF_df, OWF == "OUTSIDE"), mtry = 2, ntree = 1000)

print(RF_dist_outside) # Print the model summary

### predict the results and save
RF_comp_outside_dataset_df <- full_RF_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Random Forest", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_BIOMASS =  predict(RF_dist_outside, newdata = ., type = "response"))

RF_comp <- rbind(RF_comp_full_dataset_df, RF_comp_outside_dataset_df)
write.csv(RF_comp, "results/spatial model/surfclam/RF_abundance_comparison.csv", row.names = FALSE)

save(RF_dist_full, full_RF_df, file = "results/spatial model/surfclam/models objects/RF_dist_full.Rdata")
save(RF_dist_outside, full_RF_df, file = "results/spatial model/surfclam/models objects/RF_dist_outside.Rdata")


rm(list = setdiff(ls(), c("full_df", "data_partition_list")))

# -------------------------------------------------------------------------------------------------------- #



# 4. Bayesian Additive Regression Trees ---------------------------------------------------------------------------------------------

library(dbarts)

full_BART_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))  %>%
  mutate(PRESENCE = ifelse(BIOMASS > 0, 1, 0))


## 4.1 full dataset model ----

### 4.1.1 fit presence ----

BART_presence_full <- dbarts::bart(x.train = full_BART_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], y.train = full_BART_df[, c("PRESENCE")], keeptrees = TRUE)

invisible(BART_presence_full$fit$state) # very important if you want to use the models later

cutoff <- InformationValue::optimalCutoff(full_BART_df$PRESENCE, fitted(BART_presence_full)) # determine cutoff for presence/absence




### 4.1.2 fit biomass for data that presence = 1 ----

full_BART_PRESENCE_df <- subset(full_BART_df, PRESENCE == 1)

BART_BIOMASS_full <- dbarts::bart(x.train = full_BART_PRESENCE_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], y.train = full_BART_PRESENCE_df[, c("BIOMASS")], keeptrees = TRUE)

invisible(BART_BIOMASS_full$fit$state) # very important if you want to use the models later




### 4.1.3 predict for the entire region dataset ----

#### presence 
BART_presence_full_predict <- dbarts:::predict.bart(BART_presence_full , newdata = full_BART_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_presence_full_predict <- apply(BART_presence_full_predict, 2, median)


#### biomass
temp_presence_df <- full_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_full_predict <- dbarts:::predict.bart(BART_BIOMASS_full , newdata = temp_presence_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_biomass_full_predict <- apply(BART_biomass_full_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_full_predict)



### 4.1.4 combine biomass and presence from prediction ----

BART_comp_full_dataset_df <- full_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Bayesian Additive Regression Trees", DATASET = "FULL") %>%
  mutate(PREDICTED_PRESENCE = BART_presence_full_predict) %>%
  mutate(PREDICTED_PRESENCE = ifelse(PREDICTED_PRESENCE > cutoff, 1, 0)) %>%  # use cutoff to categorize probability into 0 and 1
  left_join(temp_presence_df) %>%
  mutate(PREDICTED_BIOMASS = ifelse(is.na(PREDICTED_BIOMASS), 0, PREDICTED_BIOMASS)) %>%  # fill in the absent data with 0
  mutate(DELTA_PREDICTION = PREDICTED_PRESENCE * PREDICTED_BIOMASS) %>%
  mutate(DELTA_PREDICTION = ifelse(DELTA_PREDICTION < 0, 0, DELTA_PREDICTION)) %>% # post-processing to avoid negative BIOMASS
  select(-c(PREDICTED_BIOMASS)) %>%
  rename(PREDICTED_BIOMASS = DELTA_PREDICTION)




## 4.2 outside dataset model (transfer ability) ----

outside_BART_df <- subset(full_BART_df, OWF == "OUTSIDE")

### 4.2.1 fit presence ----

BART_presence_outside <- dbarts::bart(x.train = outside_BART_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], y.train = outside_BART_df[, c("PRESENCE")], keeptrees = TRUE)

invisible(BART_presence_outside$fit$state) # very important if you want to use the models later

cutoff <- InformationValue::optimalCutoff(outside_BART_df$PRESENCE, fitted(BART_presence_outside)) # determine cutoff for presence/absence




### 4.2.2 fit biomass for data that presence = 1 ----

outside_BART_PRESENCE_df <- subset(outside_BART_df, PRESENCE == 1)

BART_BIOMASS_outside <- dbarts::bart(x.train = outside_BART_PRESENCE_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], y.train = outside_BART_PRESENCE_df[, c("BIOMASS")], keeptrees = TRUE)

invisible(BART_BIOMASS_outside$fit$state) # very important if you want to use the models later




### 4.2.3 predict for the entire region dataset ----

#### presence 
BART_presence_outside_predict <- dbarts:::predict.bart(BART_presence_outside , newdata = full_BART_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_presence_outside_predict <- apply(BART_presence_outside_predict, 2, median)



#### biomass
temp_presence_df <- full_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_outside_predict <- dbarts:::predict.bart(BART_BIOMASS_outside , newdata = temp_presence_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_biomass_outside_predict <- apply(BART_biomass_outside_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_outside_predict)



### 4.2.4 combine biomass and presence from prediction ----

BART_comp_outside_dataset_df <- full_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Bayesian Additive Regression Trees", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_PRESENCE = BART_presence_outside_predict) %>%
  mutate(PREDICTED_PRESENCE = ifelse(PREDICTED_PRESENCE > cutoff, 1, 0)) %>%  # use cutoff to categorize probability into 0 and 1
  left_join(temp_presence_df) %>%
  mutate(PREDICTED_BIOMASS = ifelse(is.na(PREDICTED_BIOMASS), 0, PREDICTED_BIOMASS)) %>%  # fill in the absent data with 0
  mutate(DELTA_PREDICTION = PREDICTED_PRESENCE * PREDICTED_BIOMASS) %>%
  mutate(DELTA_PREDICTION = ifelse(DELTA_PREDICTION < 0, 0, DELTA_PREDICTION)) %>% # post-processing to avoid negative BIOMASS
  select(-c(PREDICTED_BIOMASS)) %>%
  rename(PREDICTED_BIOMASS = DELTA_PREDICTION)


## 4.3 save results ----

BART_comp <- rbind(BART_comp_full_dataset_df, BART_comp_outside_dataset_df)
write.csv(BART_comp, "results/spatial model/surfclam/BART_abundance_comparison.csv", row.names = FALSE)

save(BART_presence_full, full_BART_df, file = "results/spatial model/surfclam/models objects/BART_presence_full.Rdata")
save(BART_BIOMASS_full, full_BART_df, file = "results/spatial model/surfclam/models objects/BART_BIOMASS_full.Rdata")
save(BART_presence_outside, outside_BART_df, file = "results/spatial model/surfclam/models objects/BART_presence_outside.Rdata")
save(BART_BIOMASS_outside, outside_BART_df, file = "results/spatial model/surfclam/models objects/BART_BIOMASS_outside.Rdata")

rm(list = setdiff(ls(), c("full_df", "data_partition_list")))

# -------------------------------------------------------------------------------------------------------- #



# 5. VAST ---------------------------------------------------------------------------------------------

library(VAST)
library(splines)
library(Metrics)

 ## prepare data ----
full_VAST_df <- full_df %>%  
  mutate(AreaSwept_km2 = ASWEPT_M2/1000000) %>%
  add_column(Vessel = NA, Pass = 0) %>% # required by the model, don't know why  
  rename(Lat = LAT, Lon = LON) %>% # the exact spelling is required by VAST
  mutate(Year = as.numeric(factor(YEAR)))

NWA_info <- read.csv("data/northwest_atlantic_grid.csv") # difference in area between the two sources?


 ## prepare extrapolation area ----
  
extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)],
                                                  # if use "Georges_Bank", "Mid_Atlantic_Bight", there will be 436 rows missing for LONGFIN sQUID
                                                  # all outside OWF, not affecting the gap analysis
                                                DirPath = "results/spatial model/longfin.squid/models objects/VAST/")

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

 ### below are some quick check of the generated extrapolation
plot(extrapolation_region)

length(unique(full_df$STRATUM))
length(unique(extrapolation_region$Data_Extrap$stratum_number))
setdiff(full_df$STRATUM, test$stratum_number) # stratum that are excluded in the extrapolation
nrow(subset(full_df, STRATUM %in% setdiff(full_df$STRATUM, test$stratum_number))) # tows that are excluded
unique(subset(full_df, STRATUM %in% setdiff(full_df$STRATUM, test$stratum_number))$OWF) # location of excluded tows


 ## prepare settings and covariate & formula ----

  ### create the setting file
settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[2], # number of knots, the more knots the longer it runs
                           purpose = "index2", # to calculate index
                           Region = "user", 
                            # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
                           # FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0), 
                           FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
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

  ### standardize the covariate to have mean 0 and standard deviation 1.0 as suggested by VAST
full_VAST_df$TEMP_Std <- scale(full_VAST_df$TEMP, center = TRUE, scale = TRUE)[,1]
full_VAST_df$DEPTH_Std <- scale(full_VAST_df$DEPTH, center = TRUE, scale = TRUE)[,1]

mean(full_VAST_df$TEMP_Std); sd(full_VAST_df$TEMP_Std)  
mean(full_VAST_df$DEPTH_Std); sd(full_VAST_df$DEPTH_Std)  

  ### design formulas for the two variables
p1_formula = ~ bs(TEMP_Std, degree = 2) + bs(DEPTH_Std, degree = 2)
p2_formula = ~ bs(TEMP_Std, degree = 2) + bs(DEPTH_Std, degree = 2)

  ## run full dataset model (loop over multiple error/link options to identify the best model) ----

scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description
# scenario_df <- expand.grid(error.dist = c(0, 1, 2, 4, 5, 7, 9, 10, 11, 12, 13, 14), link.func = c(0,1,2,3,4), AIC = NA, RMSE = NA) ##details see "make_data" description


# i = 4 is the optimal model
for (i in 1:nrow(scenario_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_full_model <- tryCatch(
    # attempt to fit the model
    {
      VAST_full_model <- fit_model(settings = settings,
                                   input_grid = extrapolation_region$Data_Extrap, 
                                   Lat_i = full_VAST_df$Lat,
                                   Lon_i = full_VAST_df$Lon,
                                   t_i = as.numeric(full_VAST_df$Year),
                                   b_i = full_VAST_df$BIOMASS,
                                   a_i = as_units(full_VAST_df$AreaSwept_km2, "km^2"),
                                   # covariate_data = full_VAST_df, # data frame of covariate values with columns Lat, Lon, and Year, and other columns matching names in formula
                                   # X1_formula = p1_formula,
                                   # X2_formula = p2_formula,
                                   working_dir = "results/spatial model/surfclam/models objects/VAST/")
      VAST_full_model # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
    
  )

  # skip to the next iteration if NULL is returned
  if (is.null(VAST_full_model)) { next }
  if (length(VAST_full_model$Report) == 1) {
    if (VAST_full_model$Report == "Model is not converged") { next }
  }
      
  # a quick check of the fitting convergence
  # check_fit(VAST_full_model$parameter_estimates) # no message means no issue

  
  ### predict the results and save ----
  
  # Remove units from object
  VAST_full_model$data_frame$b_i <- strip_units(VAST_full_model$data_frame$b_i) 
  VAST_full_model$data_frame$a_i <- strip_units(VAST_full_model$data_frame$a_i)
  
  VAST_comp_full_dataset_df <- full_VAST_df %>%
    select(OWF, Year, Lat, Lon, TEMP, DEPTH, AreaSwept_km2, BIOMASS) %>%
    add_column(SPECIES = "SURFCLAM", MODEL = "VAST", DATASET = "FULL") %>%
    mutate(PREDICTED_BIOMASS =  predict(x = VAST_full_model,
                                       what = "D_i",
                                       Lat_i = .$Lat,
                                       Lon_i = .$Lon,
                                       t_i = .$Year,
                                       a_i = .$AreaSwept_km2,
                                       do_checks = FALSE))
  
  scenario_df$AIC[i] = VAST_full_model$parameter_estimates$AIC
  scenario_df$RMSE[i] = rmse(VAST_comp_full_dataset_df$BIOMASS, VAST_comp_full_dataset_df$PREDICTED_BIOMASS)
  
  save(VAST_full_model, file = paste0("results/spatial model/surfclam/models objects/VAST/VAST_full_",
                                                    scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # remove(VAST_full_model, VAST_comp_full_dataset_df)
}
  
write.csv(scenario_df, "results/spatial model/surfclam/models objects/VAST/VAST_full_model_selection.csv")


  ## run outside dataset model ----
outside_VAST_df <- subset(full_VAST_df, OWF == "OUTSIDE")
scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(3,4), AIC = NA, RMSE = NA) ##details see "make_data" description

# i = 3 is the optimal model
for (i in 1:nrow(scenario_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_outside_model <- tryCatch(
    # attempt to fit the model
    {
      VAST_outside_model <- fit_model(settings = settings,
                                      input_grid = extrapolation_region$Data_Extrap, 
                                      Lat_i = outside_VAST_df$Lat,
                                      Lon_i = outside_VAST_df$Lon,
                                      t_i = as.numeric(outside_VAST_df$Year),
                                      b_i = outside_VAST_df$BIOMASS,
                                      a_i = as_units(outside_VAST_df$AreaSwept_km2, "km^2"),
                                      # covariate_data = full_VAST_df, # data frame of covariate values with columns Lat, Lon, and Year, and other columns matching names in formula
                                      # X1_formula = p1_formula,
                                      # X2_formula = p2_formula,
                                      working_dir = "results/spatial model/surfclam/models objects/VAST/")
      VAST_outside_model # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
    
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_outside_model)) { next }
  if (length(VAST_outside_model$Report) == 1) {
    if (VAST_outside_model$Report == "Model is not converged") { next }
  }
  
  # a quick check of the fitting convergence
  # check_fit(VAST_outside_model$parameter_estimates) # no message means no issue
  
  
  ### predict the results and save ----
  
  # Remove units from object
  VAST_outside_model$data_frame$b_i <- strip_units(VAST_outside_model$data_frame$b_i) 
  VAST_outside_model$data_frame$a_i <- strip_units(VAST_outside_model$data_frame$a_i)
  
  VAST_comp_outside_dataset_df <- full_VAST_df %>%
    select(OWF, Year, Lat, Lon, TEMP, DEPTH, AreaSwept_km2, BIOMASS) %>%
    add_column(SPECIES = "SURFCLAM", MODEL = "VAST", DATASET = "OUTSIDE") %>%
    mutate(PREDICTED_BIOMASS =  predict(x = VAST_outside_model,
                                       what = "D_i",
                                       Lat_i = .$Lat,
                                       Lon_i = .$Lon,
                                       t_i = .$Year,
                                       a_i = .$AreaSwept_km2,
                                       do_checks = FALSE))
  
  scenario_df$AIC[i] = VAST_outside_model$parameter_estimates$AIC
  scenario_df$RMSE[i] = rmse(VAST_comp_outside_dataset_df$BIOMASS, VAST_comp_outside_dataset_df$PREDICTED_BIOMASS)
  
  save(VAST_outside_model, file = paste0("results/spatial model/surfclam/models objects/VAST/VAST_outside_",
                                                          scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # remove(VAST_outside_model, VAST_comp_outside_dataset_df)
}

write.csv(scenario_df, "results/spatial model/surfclam/models objects/VAST/VAST_outside_model_selection.csv", row.names = FALSE)

# combine the best fit results ----

  ## load full model ----
load("results/spatial model/surfclam/models objects/VAST_full_model.Rdata")

VAST_full_model$data_frame$b_i <- strip_units(VAST_full_model$data_frame$b_i) 
VAST_full_model$data_frame$a_i <- strip_units(VAST_full_model$data_frame$a_i)

VAST_comp_full_dataset_df <- full_VAST_df %>%
  select(OWF, Year, Lat, Lon, TEMP, DEPTH, AreaSwept_km2, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "VAST", DATASET = "FULL") %>%
  mutate(PREDICTED_BIOMASS =  predict(x = VAST_full_model,
                                      what = "D_i",
                                      Lat_i = .$Lat,
                                      Lon_i = .$Lon,
                                      t_i = .$Year,
                                      a_i = .$AreaSwept_km2,
                                      do_checks = FALSE))

  ## load outside model ----
load("results/spatial model/surfclam/models objects/VAST_outside_model.Rdata")

VAST_outside_model$data_frame$b_i <- strip_units(VAST_outside_model$data_frame$b_i) 
VAST_outside_model$data_frame$a_i <- strip_units(VAST_outside_model$data_frame$a_i)

VAST_comp_outside_dataset_df <- full_VAST_df %>%
  select(OWF, Year, Lat, Lon, TEMP, DEPTH, AreaSwept_km2, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "VAST", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_BIOMASS =  predict(x = VAST_outside_model,
                                      what = "D_i",
                                      Lat_i = .$Lat,
                                      Lon_i = .$Lon,
                                      t_i = .$Year,
                                      a_i = .$AreaSwept_km2,
                                      do_checks = FALSE))

  ## assemble ----
VAST_comp <- rbind(VAST_comp_full_dataset_df, VAST_comp_outside_dataset_df)
write.csv(VAST_comp, "results/spatial model/surfclam/VAST_abundance_comparison.csv", row.names = FALSE)


rm(list = setdiff(ls(), c("full_df", "data_partition_list")))

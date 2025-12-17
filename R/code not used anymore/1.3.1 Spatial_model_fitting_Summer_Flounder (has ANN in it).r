library(tidyverse)

# 0. data set-up ---------------------------------------------------------------------------------------------

 ## prepare spatial and environmental data

tow_list_df <- read.csv("results/stratified.mean.indices/summer.flounder/tow.list.csv") # the tows used in the summer flounder assessment

station.spring.df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVSTA.csv") 
station.spring.df$SEASON <- "SPRING"

station.fall.df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVSTA.csv") 
station.fall.df$SEASON <- "FALL"

station.df <- rbind(station.spring.df, station.fall.df)

full_df <- station.df %>%
  select(CRUISE6, STRATUM, TOW, STATION, EST_YEAR, SEASON, STATUS_CODE, DECDEG_BEGLAT, DECDEG_BEGLON, SURFTEMP, BOTTEMP, AVGDEPTH) %>%
  mutate(STRATUM = ifelse(nchar(STRATUM) == 5, str_sub(STRATUM, 2, 5), STRATUM),
         ID = paste(CRUISE6, STRATUM, TOW, sep = "."),
         BLOCK = paste(EST_YEAR, SEASON, STRATUM, sep = ".")) %>%
  rename(YEAR = EST_YEAR, LAT = DECDEG_BEGLAT, LON = DECDEG_BEGLON) %>%
  filter(ID %in% unique(tow_list_df$ID))
  
## prepare the OWF overlay info

BTS_overlay_df <- read.csv("results/BTS_tows_overlay_final.csv")

full_df <- full_df %>% mutate(OWF = ifelse(ID %in% unique(BTS_overlay_df$ID), "INSIDE", "OUTSIDE"))

## prepare the tow-level catch data

catch_by_tow_df <- read.csv("results/stratified.mean.indices/summer.flounder/calibrated.total.biomass.by.tow.csv") [,-c(10, 12:15)] # fit biomass, not abundance
# catch_by_tow_df <- read.csv("results/stratified.mean.indices/summer.flounder/uncalibrated.total.catch.by.tow.csv")[,-c(11:15)]
# catch_by_tow_df <- read.csv("results/stratified.mean.indices/summer.flounder/calibrated.total.catch.by.tow.csv")[,-c(11:15)]


full_df <- full_df %>%
  mutate(STRATUM = as.integer(STRATUM)) %>%
  left_join(catch_by_tow_df) %>% # 9049 tows
  # mutate(NUMBER = as.integer(NUMBER)) %>% # convert to integer for possion distribution
  na.omit() # 7549 tows

 ## a quick breakdown of the data amount
sum(full_df$OWF == "INSIDE") 
sum(full_df$OWF == "OUTSIDE") 

remove(tow_list_df, station.spring.df, station.fall.df, BTS_overlay_df, catch_by_tow_df, station.df)

 ## look at the distribution of NUMBER
hist(full_df$BIOMASS)
sum(full_df$BIOMASS == 0)/nrow(full_df)

write.csv(full_df, "results/stratified.mean.indices/summer.flounder/full.info.tow.list.csv")

# -------------------------------------------------------------------------------------------------------- #




# 0. identify variable weighting for SDM ---------------------------------------------------------------------------------------------

 ## BRT to calculate variable weighting
# library(dismo)
# library(gbm)

# variable selection following Kleisner et al., 2017 Progress in Oceanography & 2016 PlosOne

 ## determine n.trees with gbm.step
# BRT_0 <-  gbm.step(data = full_df, gbm.x = c(8:12), gbm.y = 16,  
#                    family = "poisson", tree.complexity = 5, learning.rate = 0.03, bag.fraction = 0.75)

BRT_0 <- gbm(BIOMASS ~ factor(YEAR) + factor(SEASON) + LAT + LON + SURFTEMP + BOTTEMP + AVGDEPTH, data = full_df, 
                 distribution = "poisson", # non-negative counts
                 n.trees = 1000, bag.fraction = 0.75) # training data proportion

summary(BRT_0)
gbm.perf(BRT_0) # aiming for lower error  

BRT_final <- gbm(BIOMASS ~ factor(YEAR) + LAT + BOTTEMP + AVGDEPTH, data = full_df, 
             distribution = "poisson", 
             n.trees = 1000, bag.fraction = 0.75)
summary(BRT_final)
gbm.perf(BRT_final) # based on Residual Standard Error (RSE)

 ## Correlation analyses

 ## select correlation between variables and go back to last step
  ### previous SDM works selected YEAR + SEASON + s(SURFTEMP) + s(BOTTEMP) + s(AVGDEPTH) (Kleisner et al., 2017 Progress in Oceanography & 2016 PlosOne)
  ### KK also used delta lognormal distribution, we followed it here
library("PerformanceAnalytics")

cor.p <- chart.Correlation(full_df[,c("YEAR", "LAT", "LON", "SURFTEMP", "BOTTEMP", "AVGDEPTH")], histogram = TRUE, pch = 19)
cor.p <- chart.Correlation(full_df[,c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], histogram = TRUE, pch = 19)

## final predictor list: factor(YEAR) + LAT + BOTTEMP + AVGDEPTH
## season is excluded due to low weighting
## excluded: SEASON (low weighting); SURFTEMP, LON (high correlation)

remove(BRT_0, BRT_final)

# -------------------------------------------------------------------------------------------------------- #


# 0. reproducible 100 dataset partition for K-fold analyses ---------------------------------------------------------------------------------------------

library(caret)

set.seed(100) # always set at 100 to ensure reproducibility for all models using the same 100 datasets
data_partition_list <- vector("list", 100) # create an empty list object to save the fitted models
sum(full_df$OWF == "INSIDE")/length(full_df$OWF) # 0.09, determine the training set size based on stations in the lab 

for (i in 1:100) {
  split_temp <- createDataPartition(full_df$OWF, p = (1 - 0.09), list = FALSE)
  data_partition_list[[i]] <- split_temp
}; remove(split_temp, i)

# -------------------------------------------------------------------------------------------------------- #


# 1. Tweedie GAM ---------------------------------------------------------------------------------------------

library(mgcv)

full_WG_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))

 ## full dataset model ----

TW.GAM_dist_full <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = tw(), data = full_WG_df)

summary(TW.GAM_dist_full)
AIC(TW.GAM_dist_full)

 ### predict the results and save
TW.GAM_comp_full_dataset_df <- full_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Tweedie-GAM", DATASET = "FULL") %>%
  mutate(PREDICTED_BIOMASS =  predict(TW.GAM_dist_full, newdata = ., type = "response"))
  
## outside dataset model (transfer ability) ----

TW.GAM_dist_outside <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = tw(), data = subset(full_WG_df, OWF == "OUTSIDE"))

summary(TW.GAM_dist_outside)
AIC(TW.GAM_dist_outside)

### predict the results and save
TW.GAM_comp_outside_dataset_df <- full_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Tweedie-GAM", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_BIOMASS =  predict(TW.GAM_dist_outside, newdata = ., type = "response"))

TW.GAM_comp <- rbind(TW.GAM_comp_full_dataset_df, TW.GAM_comp_outside_dataset_df)
write.csv(TW.GAM_comp, "results/spatial model/summer.flounder/Tweedie.GAM_abundance_comparison.csv", row.names = FALSE)

save(TW.GAM_dist_full, full_WG_df, file = "results/spatial model/summer.flounder/models objects/TW.GAM_dist_full.Rdata")
save(TW.GAM_dist_outside, full_WG_df, file = "results/spatial model/summer.flounder/models objects/TW.GAM_dist_outside.Rdata")

# ## random effort reduction model (a K-fold validation approach) ----
# 
# TW.GAM_iter_model_list <- vector("list", 100) # create an empty list object to save the fitted models
# 
# for (i in 1:100) { # 100 iteration will take a few minutes
#   
#   train_temp <- full_WG_df[data_partition_list[[i]],]
#   test_temp <- full_WG_df[-data_partition_list[[i]],]
#   
#   TW.GAM_dist_temp <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = tw(), data = train_temp)
#   TW.GAM_iter_model_list[[i]] <- list(model = TW.GAM_dist_temp, train_data = train_temp, test_data = test_temp)                   
#   
#   print(paste(i, "finished at", Sys.time()))                         
# }; remove(train_temp, test_temp, TW.GAM_dist_temp)
# 
# save(TW.GAM_iter_model_list, file = "results/spatial model/summer.flounder/models objects/Tw.GAM_Kfold_100.Rdata") # note: predictions need to be made in script 1.5.x

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

DELTA.S1.GAM_dist_full <- gam(PRESENCE ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = binomial(link = "logit"), data = full_DELTA.GAM_df)
DELTA.S2.GAM_dist_full <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = gaussian(link = "log"), data = subset(full_DELTA.GAM_df, PRESENCE == 1))

summary(DELTA.S1.GAM_dist_full)
summary(DELTA.S2.GAM_dist_full)
AIC(DELTA.S1.GAM_dist_full)
AIC(DELTA.S2.GAM_dist_full)

### predict the results and save
DELTA.GAM_comp_full_dataset_df <- full_DELTA.GAM_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS, PRESENCE) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Delta-GAM", DATASET = "FULL") %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_dist_full, newdata = ., type = "response"),
         PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_dist_full, newdata = ., type = "response"))

## outside dataset model (transfer ability) ----

DELTA.S1.GAM_dist_outside <- gam(PRESENCE ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, OWF == "OUTSIDE"))
DELTA.S2.GAM_dist_outside <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = gaussian(link = "log"), data = subset(full_DELTA.GAM_df, OWF == "OUTSIDE" & PRESENCE == 1))

summary(DELTA.S1.GAM_dist_outside)
summary(DELTA.S2.GAM_dist_outside)
AIC(DELTA.S1.GAM_dist_outside)
AIC(DELTA.S2.GAM_dist_outside)

### predict the results and save
DELTA.GAM_comp_outside_dataset_df <- full_DELTA.GAM_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS, PRESENCE) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Delta-GAM", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_dist_outside, newdata = ., type = "response"),
         PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_dist_outside, newdata = ., type = "response"))

DELTA.GAM_comp <- rbind(DELTA.GAM_comp_full_dataset_df, DELTA.GAM_comp_outside_dataset_df)
write.csv(DELTA.GAM_comp, "results/spatial model/summer.flounder/Delta.GAM_abundance_comparison.csv", row.names = FALSE)

save(DELTA.S1.GAM_dist_full, full_DELTA.GAM_df, file = "results/spatial model/summer.flounder/models objects/DELTA.S1.GAM_dist_full.Rdata")
save(DELTA.S2.GAM_dist_full, full_DELTA.GAM_df, file = "results/spatial model/summer.flounder/models objects/DELTA.S2.GAM_dist_full.Rdata")
save(DELTA.S1.GAM_dist_outside, full_DELTA.GAM_df, file = "results/spatial model/summer.flounder/models objects/DELTA.S1.GAM_dist_outside.Rdata")
save(DELTA.S2.GAM_dist_outside, full_DELTA.GAM_df, file = "results/spatial model/summer.flounder/models objects/DELTA.S2.GAM_dist_outside.Rdata")

# ## random effort reduction model (a K-fold validation approach) ----
# 
# DELTA.GAM_iter_model_list <- vector("list", 100) # create an empty list object to save the fitted models
# 
# for (i in 1:100) { # 100 iteration will take a few minutes
#   
#   train_temp <- full_DELTA.GAM_df[data_partition_list[[i]],]
#   test_temp <- full_DELTA.GAM_df[-data_partition_list[[i]],]
#   
#   DELTA.S1.GAM_dist_temp <- gam(PRESENCE ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = binomial(link = "logit"), data = subset(train_temp, OWF == "OUTSIDE"))
#   DELTA.S2.GAM_dist_temp <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = gaussian(link = "log"), data = subset(train_temp, OWF == "OUTSIDE" & PRESENCE == 1))
#   
#   DELTA.GAM_iter_model_list[[i]] <- list(model.s1 = DELTA.S1.GAM_dist_temp, model.s2 = DELTA.S2.GAM_dist_temp, train_data = train_temp, test_data = test_temp)                   
#   
#   print(paste(i, "finished at", Sys.time()))                         
# }; remove(train_temp, test_temp, DELTA.S1.GAM_dist_temp, DELTA.S2.GAM_dist_temp)
# 
# save(DELTA.GAM_iter_model_list, file = "results/spatial model/summer.flounder/models objects/DELTA.GAM_Kfold_100.Rdata") # note: predictions need to be made in script 1.5.x

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
rf_tune <- tuneRF(x = train_data[,c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")],
                  y = train_data$BIOMASS,
                  ntree = 1000, # change n.tree manually 1000, 1500
                  mtryStart = 2,
                  stepFactor = 1, improve = 0.01, trace = TRUE, plot = TRUE)

remove(train_index, train_data, valid_data, rf_tune)

## full dataset model ----
 
set.seed(2) # ensure reproducibility
RF_dist_full <- randomForest(BIOMASS ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = full_RF_df, mtry = 2, ntree = 1000)

print(RF_dist_full) # Print the model summary

### predict the results and save
RF_comp_full_dataset_df <- full_RF_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Random Forest",DATASET = "FULL") %>%
  mutate(PREDICTED_BIOMASS =  predict(RF_dist_full, newdata = ., type = "response"))


## outside dataset model (transfer ability) ----

set.seed(2) # ensure reproducibility
RF_dist_outside <- randomForest(BIOMASS ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = subset(full_RF_df, OWF == "OUTSIDE"), mtry = 2, ntree = 1000)

print(RF_dist_outside) # Print the model summary

### predict the results and save
RF_comp_outside_dataset_df <- full_RF_df %>%
  # filter(OWF == "INSIDE") %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Random Forest", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_BIOMASS =  predict(RF_dist_outside, newdata = ., type = "response"))

RF_comp <- rbind(RF_comp_full_dataset_df, RF_comp_outside_dataset_df)
write.csv(RF_comp, "results/spatial model/summer.flounder/RF_abundance_comparison.csv", row.names = FALSE)

save(RF_dist_full, full_RF_df, file = "results/spatial model/summer.flounder/models objects/RF_dist_full.Rdata")
save(RF_dist_outside, full_RF_df, file = "results/spatial model/summer.flounder/models objects/RF_dist_outside.Rdata")


# ## training dataset model (a K-fold validation approach) ----
# 
# for (i in 1:100) { # 100 iteration will take 1-2 hour
#   
#   train_temp <- full_RF_df[data_partition_list[[i]],]
#   test_temp <- full_RF_df[-data_partition_list[[i]],]
# 
#   set.seed(2) # ensure reproducibility
#   RF_dist_temp <- randomForest(BIOMASS ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = train_temp, mtry = 2, ntree = 1000)
#   RF_iter_model_list <- list(model = RF_dist_temp, train_data = train_temp, test_data = test_temp)                   
#   
#   # save each iteration seperately as they are too large
#   save(RF_iter_model_list, file = paste0("results/spatial model/summer.flounder/models objects/RF_Kfold_100/", i, "_RF_model.Rdata")) # note: predictions need to be made in script 1.5.x
#   
#   print(paste(i, "finished at", Sys.time()))   
#   
# }; remove(train_temp, test_temp, RF_iter_model_list)


rm(list = setdiff(ls(), c("full_df", "data_partition_list")))

# -------------------------------------------------------------------------------------------------------- #



# 4. Artificial Neural Network ---------------------------------------------------------------------------------------------

library(nnet)

full_ANN_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))

## full dataset model ----

set.seed(1) # ensure reproducibility
ANN_dist_full <- nnet(BIOMASS ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = full_ANN_df, size = 3, # a heuristic approach to determine the number of hidden nodes = round((4 + 1) / 2)
                      maxit = 1000, linout = TRUE) # iter needs to be set to converge, lineout = TRUE means this is a regression not a classification

print(ANN_dist_full) # Print the model summary

### predict the results and save
ANN_comp_full_dataset_df <- full_ANN_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Artificial Neural Network",DATASET = "FULL") %>%
  mutate(PREDICTED_BIOMASS =  as.numeric(predict(ANN_dist_full, newdata = ., type = "raw"))) # type needs to be raw for continuous output with nnet

## outside dataset model (transfer ability) ----

set.seed(1) # ensure reproducibility
ANN_dist_outside <- nnet(BIOMASS ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = subset(full_ANN_df, OWF == "OUTSIDE"), size = 3, # a heuristic approach to determine the number of hidden nodes = round((4 + 1) / 2)
                         maxit = 1000, linout = TRUE) # iter needs to be set to converge, lineout = TRUE means this is a regression not a classification

print(ANN_dist_outside) # Print the model summary

### predict the results and save
ANN_comp_outside_dataset_df <- full_ANN_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Artificial Neural Network", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_BIOMASS =  as.numeric(predict(ANN_dist_outside, newdata = ., type = "raw")))

ANN_comp <- rbind(ANN_comp_full_dataset_df, ANN_comp_outside_dataset_df)
write.csv(ANN_comp, "results/spatial model/summer.flounder/ANN_abundance_comparison.csv", row.names = FALSE)

save(ANN_dist_full, full_ANN_df, file = "results/spatial model/summer.flounder/models objects/ANN_dist_full.Rdata")
save(ANN_dist_outside, full_ANN_df, file = "results/spatial model/summer.flounder/models objects/ANN_dist_outside.Rdata")


# ## training dataset model (a K-fold validation approach) ----
# 
# ANN_iter_model_list <- vector("list", 100) # create an empty list object to save the fitted models
# 
# for (i in 1:100) { # 100 iteration will take more than 1 hour
#   
#   train_temp <- full_ANN_df[data_partition_list[[i]],]
#   test_temp <- full_ANN_df[-data_partition_list[[i]],]
#   
#   set.seed(1) # ensure reproducibility
#   ANN_dist_temp <- nnet(BIOMASS ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = train_temp, size = 3, # a heuristic approach to determine the BIOMASS of hidden nodes = round((4 + 1) / 2)
#                         maxit = 1000, linout = TRUE)
#   ANN_iter_model_list[[i]] <- list(model = ANN_dist_temp, train_data = train_temp, test_data = test_temp)                   
#   
#   print(paste(i, "finished at", Sys.time()))                         
# }; remove(train_temp, test_temp, ANN_dist_temp)
# 
# save(ANN_iter_model_list, file = "results/spatial model/summer.flounder/models objects/ANN_Kfold_100.Rdata") # note: predictions need to be made in script 1.5.x

rm(list = setdiff(ls(), c("full_df", "data_partition_list")))

# -------------------------------------------------------------------------------------------------------- #



# 5. VAST ---------------------------------------------------------------------------------------------

library(VAST)
library(splines)
library(Metrics)

 ## prepare data ----
full_VAST_df <- full_df %>%   
  add_column(AreaSwept_km2 = 0.024, # 24000 m2 from doi:10.1093/icesjms/fsv166 and BTS protocol
             Vessel = NA, Pass = 0) %>% # required by the model, don't know why  
  # filter(YEAR %in% c(2010:2022)) %>% # just use one year to save time
  rename(Year = YEAR, Lat = LAT, Lon = LON) # the exact spelling is required by VAST

NWA_info <- read.csv("data/northwest_atlantic_grid.csv") # difference in area between the two sources?


 ## prepare extrapolation area ----
  
extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)],
                                                  # if use "Georges_Bank", "Mid_Atlantic_Bight", there will be 436 rows missing for summer flounder
                                                  # all outside OWF, not affecting the gap analysis
                                                DirPath = "results/spatial model/summer.flounder/models objects/VAST/")

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

 ### below are some quick check of the generated extrapolation
plot(extrapolation_region)

length(unique(full_df$STRATUM))
length(unique(extrapolation_region$Data_Extrap$stratum_number))
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
full_VAST_df$BOTTEMP_Std <- scale(full_VAST_df$BOTTEMP, center = TRUE, scale = TRUE)[,1]
full_VAST_df$AVGDEPTH_Std <- scale(full_VAST_df$AVGDEPTH, center = TRUE, scale = TRUE)[,1]

mean(full_VAST_df$BOTTEMP_Std); sd(full_VAST_df$BOTTEMP_Std)  
mean(full_VAST_df$AVGDEPTH_Std); sd(full_VAST_df$AVGDEPTH_Std)  

  ### design formulas for the two variables
p1_formula = ~ bs(BOTTEMP_Std, degree = 2) + bs(AVGDEPTH_Std, degree = 2)
p2_formula = ~ bs(BOTTEMP_Std, degree = 2) + bs(AVGDEPTH_Std, degree = 2)

  ## run full dataset model (loop over multiple error/link options to identify the best model) ----

scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description

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
                                   covariate_data = full_VAST_df, # data frame of covariate values with columns Lat, Lon, and Year, and other columns matching names in formula
                                   X1_formula = p1_formula,
                                   X2_formula = p2_formula,
                                   working_dir = "results/spatial model/summer.flounder/models objects/VAST/")
      VAST_full_model # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
    
  )

  # skip to the next iteration if NULL is returned
  if (is.null(VAST_full_model)) { next }
  if (length(VAST_full_model$Report) == 1) {
    if (VAST_full_model$Report != "Model is not converged") { next }
  }
      
  # a quick check of the fitting convergence
  # check_fit(VAST_full_model$parameter_estimates) # no message means no issue

  
  ### predict the results and save ----
  
  # Remove units from object
  VAST_full_model$data_frame$b_i <- strip_units(VAST_full_model$data_frame$b_i) 
  VAST_full_model$data_frame$a_i <- strip_units(VAST_full_model$data_frame$a_i)
  
  VAST_comp_full_dataset_df <- full_VAST_df %>%
    select(OWF, Year, Lat, Lon, BOTTEMP, AVGDEPTH, AreaSwept_km2, BIOMASS) %>%
    add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "VAST", DATASET = "FULL") %>%
    mutate(PREDICTED_BIOMASS =  predict(x = VAST_full_model,
                                       what = "D_i",
                                       Lat_i = .$Lat,
                                       Lon_i = .$Lon,
                                       t_i = .$Year,
                                       a_i = .$AreaSwept_km2,
                                       do_checks = FALSE))
  
  scenario_df$AIC[i] = VAST_full_model$parameter_estimates$AIC
  scenario_df$RMSE[i] = rmse(VAST_comp_full_dataset_df$BIOMASS, VAST_comp_full_dataset_df$PREDICTED_BIOMASS)
  
  save(VAST_full_model, file = paste0("results/spatial model/summer.flounder/models objects/VAST/VAST_full_",
                                                    scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # remove(VAST_full_model, VAST_comp_full_dataset_df)
}
  
write.csv(scenario_df, "results/spatial model/summer.flounder/models objects/VAST/VAST_full_model_selection.csv")


  ## run outside dataset model ----
outside_VAST_df <- subset(full_VAST_df, OWF == "OUTSIDE")
scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description

# i = 4 is the optimal model
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
                                      covariate_data = full_VAST_df, # data frame of covariate values with columns Lat, Lon, and Year, and other columns matching names in formula
                                      X1_formula = p1_formula,
                                      X2_formula = p2_formula,
                                      working_dir = "results/spatial model/summer.flounder/models objects/VAST/")
      VAST_outside_model # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
    
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_outside_model)) { next }
  if (length(VAST_outside_model$Report) == 1) {
    if (VAST_outside_model$Report != "Model is not converged") { next }
  }
  
  # a quick check of the fitting convergence
  # check_fit(VAST_outside_model$parameter_estimates) # no message means no issue
  
  
  ### predict the results and save ----
  
  # Remove units from object
  VAST_outside_model$data_frame$b_i <- strip_units(VAST_outside_model$data_frame$b_i) 
  VAST_outside_model$data_frame$a_i <- strip_units(VAST_outside_model$data_frame$a_i)
  
  VAST_comp_outside_dataset_df <- full_VAST_df %>%
    select(OWF, Year, Lat, Lon, BOTTEMP, AVGDEPTH, AreaSwept_km2, BIOMASS) %>%
    add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "VAST", DATASET = "OUTSIDE") %>%
    mutate(PREDICTED_BIOMASS =  predict(x = VAST_outside_model,
                                       what = "D_i",
                                       Lat_i = .$Lat,
                                       Lon_i = .$Lon,
                                       t_i = .$Year,
                                       a_i = .$AreaSwept_km2,
                                       do_checks = FALSE))
  
  scenario_df$AIC[i] = VAST_outside_model$parameter_estimates$AIC
  scenario_df$RMSE[i] = rmse(VAST_comp_outside_dataset_df$BIOMASS, VAST_comp_outside_dataset_df$PREDICTED_BIOMASS)
  
  save(VAST_outside_model, file = paste0("results/spatial model/summer.flounder/models objects/VAST/VAST_outside_",
                                                          scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # remove(VAST_outside_model, VAST_comp_outside_dataset_df)
}

write.csv(scenario_df, "results/spatial model/summer.flounder/models objects/VAST/VAST_outside_model_selection.csv", row.names = FALSE)

# combine the best fit results ----

 ## load full model ----
load("results/spatial model/summer.flounder/models objects/VAST_full_model.Rdata")

VAST_full_model$data_frame$b_i <- strip_units(VAST_full_model$data_frame$b_i) 
VAST_full_model$data_frame$a_i <- strip_units(VAST_full_model$data_frame$a_i)

VAST_comp_full_dataset_df <- full_VAST_df %>%
  select(OWF, Year, Lat, Lon, BOTTEMP, AVGDEPTH, AreaSwept_km2, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "VAST", DATASET = "FULL") %>%
  mutate(PREDICTED_BIOMASS =  predict(x = VAST_full_model,
                                      what = "D_i",
                                      Lat_i = .$Lat,
                                      Lon_i = .$Lon,
                                      t_i = .$Year,
                                      a_i = .$AreaSwept_km2,
                                      do_checks = FALSE))

  ## load outside model ----
load("results/spatial model/summer.flounder/models objects/VAST_outside_model.Rdata")

VAST_outside_model$data_frame$b_i <- strip_units(VAST_outside_model$data_frame$b_i) 
VAST_outside_model$data_frame$a_i <- strip_units(VAST_outside_model$data_frame$a_i)

VAST_comp_outside_dataset_df <- full_VAST_df %>%
  select(OWF, Year, Lat, Lon, BOTTEMP, AVGDEPTH, AreaSwept_km2, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "VAST", DATASET = "OUTSIDE") %>%
  mutate(PREDICTED_BIOMASS =  predict(x = VAST_outside_model,
                                      what = "D_i",
                                      Lat_i = .$Lat,
                                      Lon_i = .$Lon,
                                      t_i = .$Year,
                                      a_i = .$AreaSwept_km2,
                                      do_checks = FALSE))

  ## assemble ----

VAST_comp <- rbind(VAST_comp_full_dataset_df, VAST_comp_outside_dataset_df)
write.csv(VAST_comp, "results/spatial model/summer.flounder/VAST_abundance_comparison.csv", row.names = FALSE)


rm(list = setdiff(ls(), c("full_df", "data_partition_list")))

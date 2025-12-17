library(tidyverse)
library(visreg)


# 0. data set-up ---------------------------------------------------------------------------------------------

catch_by_tow_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") %>%
  select(ID, REGION)

full_df <- read.csv("results/stratified.mean.indices/surfclam/full.info.tow.list.csv")[,-1] %>%
  filter(OWF == "OUTSIDE") %>%
  left_join(catch_by_tow_df)

remove(catch_by_tow_df)

# -------------------------------------------------------------------------------------------------------- #



# 1. Tweedie GAM ---------------------------------------------------------------------------------------------

library(mgcv)

full_WG_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))

## GBK dataset model ----
TW.GAM_temp_GBK <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = tw(), data = subset(full_WG_df, REGION == "GBK"))

TW.GAM_GBK_YEAR_df <- visreg::visreg(TW.GAM_temp_GBK, data = subset(full_WG_df, REGION == "GBK"), "YEAR", scale = "response",  plot = FALSE)

TW.GAM_temp_GBK_df <- data.frame(YEAR = as.numeric(as.character(TW.GAM_GBK_YEAR_df$fit$YEAR)),
                                  REGION = "GBK",
                                  MODEL = "Tweedie-GAM",
                                  fit = TW.GAM_GBK_YEAR_df$fit$visregFit,
                                  lwr = TW.GAM_GBK_YEAR_df$fit$visregLwr,
                                  upr = TW.GAM_GBK_YEAR_df$fit$visregUpr) 

## SVAtoSNE dataset model ----
TW.GAM_temp_SVAtoSNE <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = tw(), data = subset(full_WG_df, REGION == "SVAtoSNE"))

TW.GAM_SVAtoSNE_YEAR_df <- visreg::visreg(TW.GAM_temp_SVAtoSNE, data = subset(full_WG_df, REGION == "SVAtoSNE"), "YEAR", scale = "response",  plot = FALSE)

TW.GAM_temp_SVAtoSNE_df <- data.frame(YEAR = as.numeric(as.character(TW.GAM_SVAtoSNE_YEAR_df$fit$YEAR)),
                                 REGION = "SVAtoSNE",
                                 MODEL = "Tweedie-GAM",
                                 fit = TW.GAM_SVAtoSNE_YEAR_df$fit$visregFit,
                                 lwr = TW.GAM_SVAtoSNE_YEAR_df$fit$visregLwr,
                                 upr = TW.GAM_SVAtoSNE_YEAR_df$fit$visregUpr)


  ## combine and save ----
TW.GAM_temp_df <- rbind(TW.GAM_temp_GBK_df, TW.GAM_temp_SVAtoSNE_df) 
write.csv(TW.GAM_temp_df, "results/temporal model/surfclam/Tweedie.GAM_YEAR_Indices.csv", row.names = FALSE)

save(TW.GAM_temp_GBK, full_WG_df, file = "results/temporal model/surfclam/models objects/TW.GAM_temp_GBK.Rdata")
save(TW.GAM_temp_SVAtoSNE, full_WG_df, file = "results/temporal model/surfclam/models objects/TW.GAM_temp_SVAtoSNE.Rdata")

rm(list = setdiff(ls(), "full_df"))

# -------------------------------------------------------------------------------------------------------- #



# 2. delta GAM ---------------------------------------------------------------------------------------------

library(mgcv)

full_DELTA.GAM_df <- full_df %>%
  mutate(PRESENCE = ifelse(BIOMASS > 0, 1, 0)) %>%
  mutate(YEAR = as.factor(YEAR))

## 2.1 GBK dataset model ----
DELTA.S1.GAM_temp_GBK <- gam(PRESENCE ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, REGION == "GBK"))
DELTA.S2.GAM_temp_GBK <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = tw(), data = subset(full_DELTA.GAM_df, REGION == "GBK" & PRESENCE == 1))

AIC(DELTA.S2.GAM_temp_GBK)

# load(file = "results/temporal model/surfclam/models objects/DELTA.s1.GAM_temp_GBK.Rdata")
# load(file = "results/temporal model/surfclam/models objects/DELTA.s2.GAM_temp_GBK.Rdata")

### generate 2-stage predicted values
DELTA.GAM_PREDICTED_GBK <- full_DELTA.GAM_df %>%
  filter(REGION == "GBK") %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, BIOMASS, PRESENCE) %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_temp_GBK, newdata = ., type = "response"),
         PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_temp_GBK, newdata = ., type = "response"))

### generate GAM model based on final prediction for year effect - tw(), Gamma(link = "log"), or gaussian(link = "log")
DELTA.GAM_fit_year_GBK <- gam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), data = DELTA.GAM_PREDICTED_GBK, family = Gamma(link = "log"),
                               nthreads = parallel::detectCores()) 

AIC(DELTA.GAM_fit_year_GBK) # gamma = -7914.249, tw = NA, gaussian = -7138.563

save(DELTA.GAM_fit_year_GBK, DELTA.GAM_PREDICTED_GBK, file = "results/temporal model/surfclam/models objects/DELTA.GAM_year_effect_GAM_GBK.Rdata")
# load("results/temporal model/surfclam/models objects/DELTA.GAM_year_effect_GAM_GBK.Rdata")

DELTA.GAM_year_effect_GBK <- visreg(DELTA.GAM_fit_year_GBK, data = DELTA.GAM_PREDICTED_GBK, "YEAR", scale = "response", plot = FALSE)

### save
DELTA.GAM_temp_GBK_df <- data.frame(YEAR = as.numeric(as.character(DELTA.GAM_year_effect_GBK$fit$YEAR)),
                                     REGION = "GBK",
                                     MODEL = "Delta-GAM",
                                     fit = DELTA.GAM_year_effect_GBK$fit$visregFit,
                                     lwr = DELTA.GAM_year_effect_GBK$fit$visregLwr,
                                     upr = DELTA.GAM_year_effect_GBK$fit$visregUpr)




  ## 2.2 SVAtoSNE dataset model ----

DELTA.S1.GAM_temp_SVAtoSNE <- gam(PRESENCE ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, REGION == "SVAtoSNE"))
DELTA.S2.GAM_temp_SVAtoSNE <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = Gamma(link = "log"), data = subset(full_DELTA.GAM_df, REGION == "SVAtoSNE" & PRESENCE == 1))

AIC(DELTA.S2.GAM_temp_SVAtoSNE)

# load(file = "results/temporal model/surfclam/models objects/DELTA.s1.GAM_temp_SVAtoSNE.Rdata")
# load(file = "results/temporal model/surfclam/models objects/DELTA.s2.GAM_temp_SVAtoSNE.Rdata")


    ### generate 2-stage predicted values
DELTA.GAM_PREDICTED_SVAtoSNE <- full_DELTA.GAM_df %>%
  filter(REGION == "SVAtoSNE") %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, BIOMASS, PRESENCE) %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_temp_SVAtoSNE, newdata = ., type = "response")) %>%
  filter(YEAR != 2014) %>%  #  2014 SVAtoSNE only has one tow of catch, remove it 
  mutate(PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_temp_SVAtoSNE, newdata = ., type = "response"))


### generate GAM model based on final prediction for year effect - tw(), Gamma(link = "log"), or gaussian(link = "log")
DELTA.GAM_fit_year_SVAtoSNE <- gam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), data = DELTA.GAM_PREDICTED_SVAtoSNE, family = Gamma(link = "log"), 
                                 nthreads = parallel::detectCores()) 

AIC(DELTA.GAM_fit_year_SVAtoSNE) # gamma = -12942.08, tw = NA, gaussian = -12433.69

save(DELTA.GAM_fit_year_SVAtoSNE, DELTA.GAM_PREDICTED_SVAtoSNE, file = "results/temporal model/surfclam/models objects/DELTA.GAM_year_effect_GAM_SVAtoSNE.Rdata")
# load("results/temporal model/surfclam/models objects/DELTA.GAM_year_effect_GAM_SVAtoSNE.Rdata")

DELTA.GAM_year_effect_SVAtoSNE <- visreg(DELTA.GAM_fit_year_SVAtoSNE, data = DELTA.GAM_PREDICTED_SVAtoSNE, "YEAR", scale = "response", plot = FALSE)

### save
DELTA.GAM_temp_SVAtoSNE_df <- data.frame(YEAR = as.numeric(as.character(DELTA.GAM_year_effect_SVAtoSNE$fit$YEAR)),
                                       REGION = "SVAtoSNE",
                                       MODEL = "Delta-GAM",
                                       fit = DELTA.GAM_year_effect_SVAtoSNE$fit$visregFit,
                                       lwr = DELTA.GAM_year_effect_SVAtoSNE$fit$visregLwr,
                                       upr = DELTA.GAM_year_effect_SVAtoSNE$fit$visregUpr)

      #### add 2014 as zero
DELTA.GAM_temp_SVAtoSNE_df <- DELTA.GAM_temp_SVAtoSNE_df %>%
  bind_rows(data.frame(YEAR = 2014, REGION = "SVAtoSNE", MODEL = "Delta-GAM", fit = 0, lwr = 0, upr = 0)) %>%
  arrange(YEAR)



## 2.3 annual dataset model ----

load("results/spatial model/surfclam/models objects/DELTA.S1.GAM_dist_outside.Rdata")
load("results/spatial model/surfclam/models objects/DELTA.S2.GAM_dist_outside.Rdata")

### generate 2-stage predicted values
DELTA.GAM_PREDICTED_annual <- full_DELTA.GAM_df %>%
  filter(OWF == "OUTSIDE") %>%
  add_column(REGION = "ANNUAL") %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, BIOMASS, PRESENCE) %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_dist_outside, newdata = ., type = "response"),
         PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_dist_outside, newdata = ., type = "response")) %>%
  mutate(PREDICTED_BIOMASS = as.numeric(PREDICTED_BIOMASS))

### generate GAM model based on final prediction for year effect - tw(), Gamma(link = "log"), or gaussian(link = "log")

DELTA.GAM_fit_year_annual <- bam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), data = DELTA.GAM_PREDICTED_annual, family = Gamma(link = "log"), 
                                 nthreads = parallel::detectCores()) 

AIC(DELTA.GAM_fit_year_annual) # gamma = -20637.47, tw = NA, gaussian = -19848.76

save(DELTA.GAM_fit_year_annual, DELTA.GAM_PREDICTED_annual, file = "results/temporal model/surfclam/models objects/DELTA.GAM_year_effect_GAM_annual.Rdata")
# load("results/temporal model/surfclam/models objects/DELTA.GAM_year_effect_GAM_annual.Rdata")

DELTA.GAM_year_effect_annual <- visreg(DELTA.GAM_fit_year_annual, data = DELTA.GAM_PREDICTED_annual, "YEAR", scale = "response", plot = FALSE)

### save
DELTA.GAM_temp_annual_df <- data.frame(YEAR = as.numeric(as.character(DELTA.GAM_year_effect_annual$fit$YEAR)),
                                       REGION = "ANNUAL",
                                       MODEL = "Delta-GAM",
                                       fit = DELTA.GAM_year_effect_annual$fit$visregFit,
                                       lwr = DELTA.GAM_year_effect_annual$fit$visregLwr,
                                       upr = DELTA.GAM_year_effect_annual$fit$visregUpr)



## 2.4 combine and save ----
DELTA.GAM_temp_df <- rbind(DELTA.GAM_temp_GBK_df, DELTA.GAM_temp_SVAtoSNE_df, DELTA.GAM_temp_annual_df) 
write.csv(DELTA.GAM_temp_df, "results/temporal model/surfclam/Delta.GAM_YEAR_Indices.csv", row.names = FALSE)

save(DELTA.S1.GAM_temp_GBK, full_DELTA.GAM_df, file = "results/temporal model/surfclam/models objects/DELTA.s1.GAM_temp_GBK.Rdata")
save(DELTA.S2.GAM_temp_GBK, full_DELTA.GAM_df, file = "results/temporal model/surfclam/models objects/DELTA.s2.GAM_temp_GBK.Rdata")
save(DELTA.S1.GAM_temp_SVAtoSNE, full_DELTA.GAM_df, file = "results/temporal model/surfclam/models objects/DELTA.s1.GAM_temp_SVAtoSNE.Rdata")
save(DELTA.S2.GAM_temp_SVAtoSNE, full_DELTA.GAM_df, file = "results/temporal model/surfclam/models objects/DELTA.s2.GAM_temp_SVAtoSNE.Rdata")

rm(list = setdiff(ls(), "full_df"))

# -------------------------------------------------------------------------------------------------------- #


# 3. RANDOM FOREST ---------------------------------------------------------------------------------------------

library(randomForest)
library(pdp) # to extract RF model results

full_RF_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))

## GBK dataset model ----

set.seed(2) # ensure reproducibility
RF_temp_GBK <- randomForest(BIOMASS ~ YEAR + LAT + TEMP + DEPTH, data = subset(full_RF_df, REGION == "GBK"), mtry = 2, ntree = 1000)

### Create Partial Dependence Plot Data
RF_GBK_year_df <- partial(RF_temp_GBK, pred.var = "YEAR", grid.resolution = 100, plot = FALSE, train = subset(full_RF_df, REGION == "GBK")) 

RF_temp_GBK_df <- data.frame(YEAR = as.numeric(as.character(RF_GBK_year_df$YEAR)),
                              REGION = "GBK",
                              MODEL = "Random Forest",
                              fit = RF_GBK_year_df$yhat, #colMeans(boot_fits),
                              lwr = NA, #apply(boot_fits, 2, quantile, probs = 0.025),
                              upr = NA) #apply(boot_fits, 2, quantile, probs = 0.975))

## SVAtoSNE dataset model ----

set.seed(2) # ensure reproducibility
RF_temp_SVAtoSNE <- randomForest(BIOMASS ~ YEAR + LAT + TEMP + DEPTH, data = subset(full_RF_df, REGION == "SVAtoSNE"), mtry = 2, ntree = 1000)

### Create Partial Dependence Plot Data
RF_SVAtoSNE_year_df <- partial(RF_temp_SVAtoSNE, pred.var = "YEAR", grid.resolution = 100, plot = FALSE, train = subset(full_RF_df, REGION == "SVAtoSNE")) 

RF_temp_SVAtoSNE_df <- data.frame(YEAR = as.numeric(as.character(RF_SVAtoSNE_year_df$YEAR)),
                                REGION = "SVAtoSNE",
                                MODEL = "Random Forest",
                                fit = RF_SVAtoSNE_year_df$yhat, #colMeans(boot_fits),
                                lwr = NA, #apply(boot_fits, 2, quantile, probs = 0.025),
                                upr = NA) #apply(boot_fits, 2, quantile, probs = 0.975))


## combine and save ----
RF_temp_df <- rbind(RF_temp_GBK_df, RF_temp_SVAtoSNE_df) 
write.csv(RF_temp_df, "results/temporal model/surfclam/RF_YEAR_Indices.csv", row.names = FALSE)

save(RF_temp_GBK, full_RF_df, file = "results/temporal model/surfclam/models objects/RF_temp_GBK.Rdata")
save(RF_temp_SVAtoSNE, full_RF_df, file = "results/temporal model/surfclam/models objects/RF_temp_SVAtoSNE.Rdata")

rm(list = setdiff(ls(), "full_df"))


# -------------------------------------------------------------------------------------------------------- #



# 4. Bayesian Additive Regression Trees ---------------------------------------------------------------------------------------------
library(dbarts)
library(pdp) # to extract BART model results

full_BART_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR)) %>%
  mutate(PRESENCE = ifelse(BIOMASS > 0, 1, 0))

## 4.1 GBK dataset model ----

GBK_BART_df <- subset(full_BART_df, REGION == "GBK")


### 4.1.1 fit presence ----

BART_presence_GBK <- dbarts::bart(x.train = GBK_BART_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], y.train = GBK_BART_df[, c("PRESENCE")], keeptrees = TRUE)

invisible(BART_presence_GBK$fit$state) # very important if you want to use the models later

cutoff <- InformationValue::optimalCutoff(GBK_BART_df$PRESENCE, fitted(BART_presence_GBK)) # determine cutoff for presence/absence


### 4.1.2 fit biomass for data that presence = 1 ----

GBK_BART_PRESENCE_df <- subset(GBK_BART_df, PRESENCE == 1)

BART_BIOMASS_GBK <- dbarts::bart(x.train = GBK_BART_PRESENCE_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], y.train = GBK_BART_PRESENCE_df[, c("BIOMASS")], keeptrees = TRUE)

invisible(BART_BIOMASS_GBK$fit$state) # very important if you want to use the models later


### 4.1.3 predict for the entire region dataset ----

#### presence 
BART_presence_GBK_predict <- dbarts:::predict.bart(BART_presence_GBK , newdata = GBK_BART_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_presence_GBK_predict <- apply(BART_presence_GBK_predict, 2, median)


#### biomass
temp_presence_df <- GBK_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_GBK_predict <- dbarts:::predict.bart(BART_BIOMASS_GBK , newdata = temp_presence_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_biomass_GBK_predict <- apply(BART_biomass_GBK_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_GBK_predict)


### 4.1.4 combine biomass and presence from prediction ----

BART_comp_GBK_dataset_df <- GBK_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Bayesian Additive Regression Trees", DATASET = "GBK") %>%
  mutate(PREDICTED_PRESENCE = BART_presence_GBK_predict) %>%
  mutate(PREDICTED_PRESENCE = ifelse(PREDICTED_PRESENCE > cutoff, 1, 0)) %>%  # use cutoff to categorize probability into 0 and 1
  left_join(temp_presence_df) %>%
  mutate(PREDICTED_BIOMASS = ifelse(is.na(PREDICTED_BIOMASS), 0, PREDICTED_BIOMASS)) %>%  # fill in the absent data with 0
  mutate(DELTA_PREDICTION = PREDICTED_PRESENCE * PREDICTED_BIOMASS) %>%
  mutate(DELTA_PREDICTION = ifelse(DELTA_PREDICTION < 0, 0, DELTA_PREDICTION)) %>% # post-processing to avoid negative BIOMASS
  select(-c(PREDICTED_BIOMASS)) %>%
  rename(PREDICTED_BIOMASS = DELTA_PREDICTION)



### 4.1.5 fit GAM for year effect ---- 

####  - tw(), Gamma(link = "log"), or gaussian(link = "log")
BART.GAM_fit_year_GBK <- bam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), data = BART_comp_GBK_dataset_df, family = tw(), 
                              nthreads = parallel::detectCores()) 

AIC(BART.GAM_fit_year_GBK) # gamma = NA, tw = -787.3621, gaussian = NA

save(BART.GAM_fit_year_GBK, BART_comp_GBK_dataset_df, file = "results/temporal model/surfclam/models objects/BART_year_effect_GAM_GBK.Rdata")
# load("results/temporal model/surfclam/models objects/BART_year_effect_GAM_GBK.Rdata")

BART.GAM_year_effect_GBK <- visreg(BART.GAM_fit_year_GBK, data = BART_comp_GBK_dataset_df, "YEAR", scale = "response", plot = FALSE)

BART_temp_GBK_df <- data.frame(YEAR = as.numeric(as.character(BART.GAM_year_effect_GBK$fit$YEAR)),
                                REGION = "GBK",
                                MODEL = "Bayesian Additive Regression Trees",
                                fit = BART.GAM_year_effect_GBK$fit$visregFit, 
                                lwr = BART.GAM_year_effect_GBK$fit$visregLwr,
                                upr = BART.GAM_year_effect_GBK$fit$visregUpr) 



## 4.2 SVAtoSNE dataset model ----

SVAtoSNE_BART_df <- subset(full_BART_df, REGION == "SVAtoSNE")


### 4.2.1 fit presence ----

BART_presence_SVAtoSNE <- dbarts::bart(x.train = SVAtoSNE_BART_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], y.train = SVAtoSNE_BART_df[, c("PRESENCE")], keeptrees = TRUE)

invisible(BART_presence_SVAtoSNE$fit$state) # very important if you want to use the models later

cutoff <- InformationValue::optimalCutoff(SVAtoSNE_BART_df$PRESENCE, fitted(BART_presence_SVAtoSNE)) # determine cutoff for presence/absence


### 4.2.2 fit biomass for data that presence = 1 ----

SVAtoSNE_BART_PRESENCE_df <- subset(SVAtoSNE_BART_df, PRESENCE == 1)

BART_BIOMASS_SVAtoSNE <- dbarts::bart(x.train = SVAtoSNE_BART_PRESENCE_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], y.train = SVAtoSNE_BART_PRESENCE_df[, c("BIOMASS")], keeptrees = TRUE)

invisible(BART_BIOMASS_SVAtoSNE$fit$state) # very important if you want to use the models later


### 4.2.3 predict for the entire region dataset ----

#### presence 
BART_presence_SVAtoSNE_predict <- dbarts:::predict.bart(BART_presence_SVAtoSNE , newdata = SVAtoSNE_BART_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_presence_SVAtoSNE_predict <- apply(BART_presence_SVAtoSNE_predict, 2, median)


#### biomass
temp_presence_df <- SVAtoSNE_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_SVAtoSNE_predict <- dbarts:::predict.bart(BART_BIOMASS_SVAtoSNE , newdata = temp_presence_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_biomass_SVAtoSNE_predict <- apply(BART_biomass_SVAtoSNE_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_SVAtoSNE_predict)


### 4.2.4 combine biomass and presence from prediction ----

BART_comp_SVAtoSNE_dataset_df <- SVAtoSNE_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Bayesian Additive Regression Trees", DATASET = "SVAtoSNE") %>%
  mutate(PREDICTED_PRESENCE = BART_presence_SVAtoSNE_predict) %>%
  mutate(PREDICTED_PRESENCE = ifelse(PREDICTED_PRESENCE > cutoff, 1, 0)) %>%  # use cutoff to categorize probability into 0 and 1
  left_join(temp_presence_df) %>%
  mutate(PREDICTED_BIOMASS = ifelse(is.na(PREDICTED_BIOMASS), 0, PREDICTED_BIOMASS)) %>%  # fill in the absent data with 0
  mutate(DELTA_PREDICTION = PREDICTED_PRESENCE * PREDICTED_BIOMASS) %>%
  mutate(DELTA_PREDICTION = ifelse(DELTA_PREDICTION < 0, 0, DELTA_PREDICTION)) %>% # post-processing to avoid negative BIOMASS
  select(-c(PREDICTED_BIOMASS)) %>%
  rename(PREDICTED_BIOMASS = DELTA_PREDICTION)



### 4.2.5 fit GAM for year effect ---- 

####  - tw(), Gamma(link = "log"), or gaussian(link = "log")
BART.GAM_fit_year_SVAtoSNE <- bam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), data = BART_comp_SVAtoSNE_dataset_df, family = tw(), 
                                nthreads = parallel::detectCores()) 

AIC(BART.GAM_fit_year_SVAtoSNE) # gamma = NA, tw = -2619.904, gaussian = NA

save(BART.GAM_fit_year_SVAtoSNE, BART_comp_SVAtoSNE_dataset_df, file = "results/temporal model/surfclam/models objects/BART_year_effect_GAM_SVAtoSNE.Rdata")
# load("results/temporal model/surfclam/models objects/BART_year_effect_GAM_SVAtoSNE.Rdata")

BART.GAM_year_effect_SVAtoSNE <- visreg(BART.GAM_fit_year_SVAtoSNE, data = BART_comp_SVAtoSNE_dataset_df, "YEAR", scale = "response", plot = FALSE)

BART_temp_SVAtoSNE_df <- data.frame(YEAR = as.numeric(as.character(BART.GAM_year_effect_SVAtoSNE$fit$YEAR)),
                                  REGION = "SVAtoSNE",
                                  MODEL = "Bayesian Additive Regression Trees",
                                  fit = BART.GAM_year_effect_SVAtoSNE$fit$visregFit, 
                                  lwr = BART.GAM_year_effect_SVAtoSNE$fit$visregLwr,
                                  upr = BART.GAM_year_effect_SVAtoSNE$fit$visregUpr) 



## 4.3 Annual dataset model ----

load("results/spatial model/surfclam/models objects/BART_presence_outside.Rdata")
load("results/spatial model/surfclam/models objects/BART_BIOMASS_outside.Rdata")

invisible(BART_presence_outside$fit$state)
invisible(BART_BIOMASS_outside$fit$state) 

cutoff <- InformationValue::optimalCutoff(outside_BART_df$PRESENCE, fitted(BART_presence_outside))

### 4.3.1 predict for the entire region dataset ----

#### presence 
BART_presence_annual_predict <- dbarts:::predict.bart(BART_presence_outside, newdata = outside_BART_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_presence_annual_predict <- apply(BART_presence_annual_predict, 2, median)

#### biomass
temp_presence_df <- outside_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_annual_predict <- dbarts:::predict.bart(BART_BIOMASS_outside , newdata = temp_presence_df[, c("YEAR", "LAT", "TEMP", "DEPTH")], type = "bart")
BART_biomass_annual_predict <- apply(BART_biomass_annual_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_annual_predict)


### 4.3.2 combine biomass and presence from prediction ----

BART_comp_annual_dataset_df <- outside_BART_df %>%
  select(OWF, YEAR, LAT, TEMP, DEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SURFCLAM", MODEL = "Bayesian Additive Regression Trees", DATASET = "ANNUAL") %>%
  mutate(PREDICTED_PRESENCE = BART_presence_annual_predict) %>%
  mutate(PREDICTED_PRESENCE = ifelse(PREDICTED_PRESENCE > cutoff, 1, 0)) %>%  # use cutoff to categorize probability into 0 and 1
  left_join(temp_presence_df) %>%
  mutate(PREDICTED_BIOMASS = ifelse(is.na(PREDICTED_BIOMASS), 0, PREDICTED_BIOMASS)) %>%  # fill in the absent data with 0
  mutate(DELTA_PREDICTION = PREDICTED_PRESENCE * PREDICTED_BIOMASS) %>%
  mutate(DELTA_PREDICTION = ifelse(DELTA_PREDICTION < 0, 0, DELTA_PREDICTION)) %>% # post-processing to avoid negative BIOMASS
  select(-c(PREDICTED_BIOMASS)) %>%
  rename(PREDICTED_BIOMASS = DELTA_PREDICTION)


### 4.3.3 fit GAM for year effect ---- 

####  - tw(), Gamma(link = "log"), or gaussian(link = "log")
BART.GAM_fit_year_annual <- bam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), data = BART_comp_annual_dataset_df, family = tw(), 
                                nthreads = parallel::detectCores()) 

AIC(BART.GAM_fit_year_annual) # gamma = NA, tw = 14691.9, gaussian = NA

save(BART.GAM_fit_year_annual, BART_comp_annual_dataset_df, file = "results/temporal model/surfclam/models objects/BART_year_effect_GAM_annual.Rdata")
# load("results/temporal model/surfclam/models objects/BART_year_effect_GAM_annual.Rdata")

BART.GAM_year_effect_annual <- visreg(BART.GAM_fit_year_annual, data = BART_comp_annual_dataset_df, "YEAR", scale = "response", plot = FALSE)

BART_temp_annual_df <- data.frame(YEAR = as.numeric(as.character(BART.GAM_year_effect_annual$fit$YEAR)),
                                  REGION = "ANNUAL",
                                  MODEL = "Bayesian Additive Regression Trees",
                                  fit = BART.GAM_year_effect_annual$fit$visregFit, 
                                  lwr = BART.GAM_year_effect_annual$fit$visregLwr,
                                  upr = BART.GAM_year_effect_annual$fit$visregUpr) 


## 4.4 combine and save ----

BART_temp_df <- rbind(BART_temp_GBK_df, BART_temp_SVAtoSNE_df, BART_temp_annual_df) %>%
  mutate(fit = ifelse(fit < 0, 0, fit))  # post-processing to ensure minimal value should be 0

write.csv(BART_temp_df, "results/temporal model/surfclam/BART_YEAR_Indices.csv", row.names = FALSE)

save(BART_presence_GBK, GBK_BART_df, file = "results/temporal model/surfclam/models objects/BART_presence_GBK.Rdata")
save(BART_BIOMASS_GBK, GBK_BART_df, file = "results/temporal model/surfclam/models objects/BART_BIOMASS_GBK.Rdata")
save(BART_presence_SVAtoSNE, SVAtoSNE_BART_df, file = "results/temporal model/surfclam/models objects/BART_presence_SVAtoSNE.Rdata")
save(BART_BIOMASS_SVAtoSNE, SVAtoSNE_BART_df, file = "results/temporal model/surfclam/models objects/BART_BIOMASS_SVAtoSNE.Rdata")

rm(list = setdiff(ls(), "full_df"))



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

full_VAST_df$TEMP_Std <- scale(full_VAST_df$TEMP, center = TRUE, scale = TRUE)[,1]
full_VAST_df$DEPTH_Std <- scale(full_VAST_df$DEPTH, center = TRUE, scale = TRUE)[,1]

    ### design formulas for the two variables
p1_formula = ~ bs(TEMP_Std, degree = 2) + bs(DEPTH_Std, degree = 2)
p2_formula = ~ bs(TEMP_Std, degree = 2) + bs(DEPTH_Std, degree = 2)



    ## GBK dataset model ----
GBK_VAST_df <- subset(full_VAST_df, REGION == "GBK") %>%
  mutate(Year = as.numeric(factor(YEAR)))

      ### standardize the covariate to have mean 0 and standard deviation 1.0 as suggested by VAST
GBK_VAST_df$TEMP_Std <- scale(GBK_VAST_df$TEMP, center = TRUE, scale = TRUE)[,1]
GBK_VAST_df$DEPTH_Std <- scale(GBK_VAST_df$DEPTH, center = TRUE, scale = TRUE)[,1]


      ### loop over multiple error/link options to identify the best model ----

scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(3,4), AIC = NA, RMSE = NA) ##details see "make_data" description

# i = 1 is the optimal model
for (i in 1:nrow(scenario_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  ## prepare extrapolation area ----
  
  extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                  strata.limits = data.frame(STRATA = "EPU"),
                                                  epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2)],
                                                  # if use "Georges_Bank", "Mid_Atlantic_Bight", there will be 436 rows missing for summer flounder
                                                  # all outside OWF, not affecting the gap analysis
                                                  DirPath = "results/temporal model/surfclam/models objects/VAST/")
  
  colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_temp_GBK <- tryCatch(
    # attempt to fit the model
    {
      VAST_temp_GBK <- fit_model(settings = settings,
                                  input_grid = extrapolation_region$Data_Extrap, 
                                  Lat_i = GBK_VAST_df$Lat,
                                  Lon_i = GBK_VAST_df$Lon,
                                  t_i = as.numeric(GBK_VAST_df$Year), 
                                  b_i = GBK_VAST_df$BIOMASS,
                                  a_i = as_units(GBK_VAST_df$AreaSwept_km2, "km^2"),
                                  # covariate_data = GBK_VAST_df, 
                                  # X1_formula = p1_formula,
                                  # X2_formula = p2_formula,
                                  working_dir = "results/temporal model/surfclam/models objects/VAST/")
      VAST_temp_GBK # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_GBK)) { next }
  if (length(VAST_temp_GBK$Report) == 1) {
    if (VAST_temp_GBK$Report == "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_GBK$parameter_estimates$AIC
  
  save(VAST_temp_GBK, file = paste0("results/temporal model/surfclam/models objects/VAST_GBK_",
                                     scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # save(VAST_temp_GBK, full_VAST_df, file = paste0("results/temporal model/surfclam/models objects/VAST_temp_GBK.Rdata"))
  # load(file = paste0("results/temporal model/surfclam/models objects/VAST_temp_GBK.Rdata"))
  
}

write.csv(scenario_df, "results/temporal model/surfclam/models objects/VAST/VAST_GBK_model_selection.csv")

### save the fitted value ----

# load("results/temporal model/surfclam/models objects/VAST_temp_GBK.Rdata")

VAST_GBK_data <- plot(VAST_temp_GBK, working_dir = "results/temporal model/surfclam/models objects/VAST/figure/")

VAST_temp_GBK_df <- data.frame(YEAR = VAST_GBK_data$Index$Table$Time,
                                REGION = "GBK",
                                MODEL = "VAST",
                                fit = VAST_GBK_data$Index$Table$Estimate, 
                                lwr = VAST_GBK_data$Index$Table$Estimate - 1.96 * VAST_GBK_data$Index$Table$`Std. Error for Estimate`, 
                                upr = VAST_GBK_data$Index$Table$Estimate + 1.96 * VAST_GBK_data$Index$Table$`Std. Error for Estimate`) 




## SVAtoSNE dataset model ----
SVAtoSNE_VAST_df <- subset(full_VAST_df, REGION == "SVAtoSNE") %>%
  mutate(Year = as.numeric(factor(YEAR)))

### standardize the covariate to have mean 0 and standard deviation 1.0 as suggested by VAST
SVAtoSNE_VAST_df$TEMP_Std <- scale(SVAtoSNE_VAST_df$TEMP, center = TRUE, scale = TRUE)[,1]
SVAtoSNE_VAST_df$DEPTH_Std <- scale(SVAtoSNE_VAST_df$DEPTH, center = TRUE, scale = TRUE)[,1]

### loop over multiple error/link options to identify the best model ----

scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(3,4), AIC = NA, RMSE = NA) ##details see "make_data" description

# i = 3 is the optimal model
for (i in 1:nrow(scenario_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  ## prepare extrapolation area ----
  
  extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                  strata.limits = data.frame(STRATA = "EPU"),
                                                  epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(3)],
                                                  DirPath = "results/temporal model/surfclam/models objects/VAST/")
  
  colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_temp_SVAtoSNE <- tryCatch(
    # attempt to fit the model
    {
      VAST_temp_SVAtoSNE <- fit_model(settings = settings,
                                    input_grid = extrapolation_region$Data_Extrap, 
                                    Lat_i = SVAtoSNE_VAST_df$Lat,
                                    Lon_i = SVAtoSNE_VAST_df$Lon,
                                    t_i = as.numeric(SVAtoSNE_VAST_df$Year), 
                                    b_i = SVAtoSNE_VAST_df$BIOMASS,
                                    a_i = as_units(SVAtoSNE_VAST_df$AreaSwept_km2, "km^2"),
                                    # covariate_data = SVAtoSNE_VAST_df, # use full_VAST_df here to avoid missing year (2020) for GBK
                                    # X1_formula = p1_formula,
                                    # X2_formula = p2_formula,
                                    working_dir = "results/temporal model/surfclam/models objects/VAST/")
      VAST_temp_SVAtoSNE # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_SVAtoSNE)) { next }
  if (length(VAST_temp_SVAtoSNE$Report) == 1) {
    if (VAST_temp_SVAtoSNE$Report == "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_SVAtoSNE$parameter_estimates$AIC
  
  save(VAST_temp_SVAtoSNE, file = paste0("results/temporal model/surfclam/models objects/VAST_temp_SVAtoSNE_",
                                       scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # save(VAST_temp_SVAtoSNE, file = paste0("results/temporal model/surfclam/models objects/VAST_temp_SVAtoSNE.Rdata"))
  
  # load(file = paste0("results/temporal model/surfclam/models objects/VAST_temp_SVAtoSNE.Rdata"))
}

write.csv(scenario_df, "results/temporal model/surfclam/models objects/VAST/VAST_SVAtoSNE_model_selection.csv")




### save the fitted value ----

  ## GBK
# load("results/temporal model/surfclam/models objects/VAST_temp_GBK.Rdata")

# VAST_GBK_data <- plot(VAST_temp_GBK, working_dir = "results/temporal model/surfclam/models objects/VAST/figure/")

VAST_temp_GBK_df <- data.frame(YEAR = unique(GBK_VAST_df$YEAR),
                                    REGION = "GBK",
                                    MODEL = "VAST",
                                    fit = NA, 
                                    lwr = NA, 
                                    upr = NA) 

  ## SVAtoSNE
load("results/temporal model/surfclam/models objects/VAST_temp_SVAtoSNE.Rdata")

VAST_SVAtoSNE_data <- plot(VAST_temp_SVAtoSNE, working_dir = "results/temporal model/surfclam/models objects/VAST/figure/")

VAST_temp_SVAtoSNE_df <- data.frame(YEAR = unique(SVAtoSNE_VAST_df$YEAR),
                                  REGION = "SVAtoSNE",
                                  MODEL = "VAST",
                                  fit = VAST_SVAtoSNE_data$Index$Table$Estimate, 
                                  lwr = VAST_SVAtoSNE_data$Index$Table$Estimate - 1.96 * VAST_SVAtoSNE_data$Index$Table$`Std. Error for Estimate`, 
                                  upr = VAST_SVAtoSNE_data$Index$Table$Estimate + 1.96 * VAST_SVAtoSNE_data$Index$Table$`Std. Error for Estimate`) 


## load all region outside model ----
load("results/spatial model/surfclam/models objects/VAST_outside_model.Rdata")

VAST_all_data <- plot(VAST_outside_model, working_dir = "results/temporal model/surfclam/models objects/VAST/figure/")

VAST_temp_all_df <- data.frame(YEAR = unique(full_VAST_df$YEAR),
                                  REGION = "ALL",
                                  MODEL = "VAST",
                                  fit = VAST_all_data$Index$Table$Estimate, 
                                  lwr = VAST_all_data$Index$Table$Estimate - 1.96 * VAST_all_data$Index$Table$`Std. Error for Estimate`, 
                                  upr = VAST_all_data$Index$Table$Estimate + 1.96 * VAST_all_data$Index$Table$`Std. Error for Estimate`) 
  ## combine and save ----
VAST_temp_df <- rbind(VAST_temp_GBK_df, VAST_temp_SVAtoSNE_df, VAST_temp_all_df) 
write.csv(VAST_temp_df, "results/temporal model/surfclam/VAST_YEAR_Indices.csv", row.names = FALSE)

rm(list = setdiff(ls(), "full_df"))

# -------------------------------------------------------------------------------------------------------- #




# 6. simple delta GAM ---------------------------------------------------------------------------------------------

library(mgcv)

full_DELTA.GAM_df <- full_df %>%
  mutate(PRESENCE = ifelse(BIOMASS > 0, 1, 0)) %>%
  mutate(YEAR = as.factor(YEAR))

## 2.1 GBK dataset model ----
DELTA.S1.GAM_temp_GBK <- gam(PRESENCE ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, REGION == "GBK"))
DELTA.S2.GAM_temp_GBK <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = tw(), data = subset(full_DELTA.GAM_df, REGION == "GBK" & PRESENCE == 1))

AIC(DELTA.S2.GAM_temp_GBK)

# load(file = "results/temporal model/surfclam/models objects/DELTA.s1.GAM_temp_GBK.Rdata")
# load(file = "results/temporal model/surfclam/models objects/DELTA.s2.GAM_temp_GBK.Rdata")

### generate 2-stage predicted values
DELTA.s1.GAM_GBK_df <- visreg(DELTA.S1.GAM_temp_GBK, data = subset(full_DELTA.GAM_df, REGION == "GBK"), "YEAR", scale = "response", plot = FALSE)
DELTA.s2.GAM_GBK_df <- visreg(DELTA.S2.GAM_temp_GBK, data = subset(full_DELTA.GAM_df, REGION == "GBK" & PRESENCE == 1), "YEAR", scale = "response", plot = FALSE)

DELTA.GAM_temp_GBK_df <- data.frame(YEAR = as.numeric(as.character(DELTA.s1.GAM_GBK_df$fit$YEAR)),
                                    REGION = "GBK",
                                    MODEL = "Delta-GAM",
                                    fit = DELTA.s1.GAM_GBK_df$fit$visregFit * DELTA.s2.GAM_GBK_df$fit$visregFit,
                                    lwr = DELTA.s1.GAM_GBK_df$fit$visregLwr * DELTA.s2.GAM_GBK_df$fit$visregLwr,
                                    upr = DELTA.s1.GAM_GBK_df$fit$visregUpr * DELTA.s2.GAM_GBK_df$fit$visregUpr)




## 2.2 SVAtoSNE dataset model ----

DELTA.S1.GAM_temp_SVAtoSNE <- gam(PRESENCE ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, REGION == "SVAtoSNE"))
DELTA.S2.GAM_temp_SVAtoSNE <- gam(BIOMASS ~ YEAR + s(LAT) + s(TEMP) + s(DEPTH), family = Gamma(link = "log"), data = subset(full_DELTA.GAM_df, REGION == "SVAtoSNE" & PRESENCE == 1))

AIC(DELTA.S2.GAM_temp_SVAtoSNE)

# load(file = "results/temporal model/surfclam/models objects/DELTA.s1.GAM_temp_SVAtoSNE.Rdata")
# load(file = "results/temporal model/surfclam/models objects/DELTA.s2.GAM_temp_SVAtoSNE.Rdata")


### generate 2-stage predicted values
DELTA.s1.GAM_SVAtoSNE_df <- visreg(DELTA.S1.GAM_temp_SVAtoSNE, data = subset(full_DELTA.GAM_df, REGION == "SVAtoSNE"), "YEAR", scale = "response", plot = FALSE)
DELTA.s2.GAM_SVAtoSNE_df <- visreg(DELTA.S2.GAM_temp_SVAtoSNE, data = subset(full_DELTA.GAM_df, REGION == "SVAtoSNE" & PRESENCE == 1), "YEAR", scale = "response", plot = FALSE)

DELTA.GAM_temp_SVAtoSNE_df <- data.frame(YEAR = as.numeric(as.character(DELTA.s2.GAM_SVAtoSNE_df$fit$YEAR)), # 2014 for s1 is 0
                                    REGION = "SVAtoSNE",
                                    MODEL = "Delta-GAM",
                                    fit = DELTA.s1.GAM_SVAtoSNE_df$fit$visregFit[-7] * DELTA.s2.GAM_SVAtoSNE_df$fit$visregFit, # 2014 for s1 is 0
                                    lwr = DELTA.s1.GAM_SVAtoSNE_df$fit$visregLwr[-7] * DELTA.s2.GAM_SVAtoSNE_df$fit$visregLwr,
                                    upr = DELTA.s1.GAM_SVAtoSNE_df$fit$visregUpr[-7] * DELTA.s2.GAM_SVAtoSNE_df$fit$visregUpr)


#### add 2014 as zero
DELTA.GAM_temp_SVAtoSNE_df <- DELTA.GAM_temp_SVAtoSNE_df %>%
  bind_rows(data.frame(YEAR = 2014, REGION = "SVAtoSNE", MODEL = "Delta-GAM", fit = 0, lwr = 0, upr = 0)) %>%
  arrange(YEAR)



## 2.3 annual dataset model ----

load("results/spatial model/surfclam/models objects/DELTA.S1.GAM_dist_outside.Rdata")
load("results/spatial model/surfclam/models objects/DELTA.S2.GAM_dist_outside.Rdata")

### generate 2-stage predicted values
DELTA.s1.GAM_SVAtoSNE_df <- visreg(DELTA.S1.GAM_temp_SVAtoSNE, data = subset(full_DELTA.GAM_df, REGION == "SVAtoSNE"), "YEAR", scale = "response", plot = FALSE)
DELTA.s2.GAM_SVAtoSNE_df <- visreg(DELTA.S2.GAM_temp_SVAtoSNE, data = subset(full_DELTA.GAM_df, REGION == "SVAtoSNE" & PRESENCE == 1), "YEAR", scale = "response", plot = FALSE)


DELTA.GAM_year_effect_annual <- visreg(DELTA.GAM_fit_year_annual, data = DELTA.GAM_PREDICTED_annual, "YEAR", scale = "response", plot = FALSE)

### save
DELTA.GAM_temp_annual_df <- data.frame(YEAR = as.numeric(as.character(DELTA.GAM_year_effect_annual$fit$YEAR)),
                                       REGION = "ANNUAL",
                                       MODEL = "Delta-GAM",
                                       fit = DELTA.GAM_year_effect_annual$fit$visregFit,
                                       lwr = DELTA.GAM_year_effect_annual$fit$visregLwr,
                                       upr = DELTA.GAM_year_effect_annual$fit$visregUpr)



## 2.4 combine and save ----
DELTA.GAM_temp_df <- rbind(DELTA.GAM_temp_GBK_df, DELTA.GAM_temp_SVAtoSNE_df) 
write.csv(DELTA.GAM_temp_df, "results/temporal model/surfclam/simple.Delta.GAM_YEAR_Indices.csv", row.names = FALSE)

save(DELTA.S1.GAM_temp_GBK, full_DELTA.GAM_df, file = "results/temporal model/surfclam/models objects/DELTA.s1.GAM_temp_GBK.Rdata")
save(DELTA.S2.GAM_temp_GBK, full_DELTA.GAM_df, file = "results/temporal model/surfclam/models objects/DELTA.s2.GAM_temp_GBK.Rdata")
save(DELTA.S1.GAM_temp_SVAtoSNE, full_DELTA.GAM_df, file = "results/temporal model/surfclam/models objects/DELTA.s1.GAM_temp_SVAtoSNE.Rdata")
save(DELTA.S2.GAM_temp_SVAtoSNE, full_DELTA.GAM_df, file = "results/temporal model/surfclam/models objects/DELTA.s2.GAM_temp_SVAtoSNE.Rdata")

rm(list = setdiff(ls(), "full_df"))

# -------------------------------------------------------------------------------------------------------- #

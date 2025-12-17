library(tidyverse)
library(visreg)


# 0. data set-up ---------------------------------------------------------------------------------------------

full_df <- read.csv("results/stratified.mean.indices/summer.flounder/full.info.tow.list.csv")[,-1] %>%
  filter(OWF == "OUTSIDE")

# -------------------------------------------------------------------------------------------------------- #



# 1. Tweedie GAM ---------------------------------------------------------------------------------------------

library(mgcv)

full_WG_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))

  ## fall dataset model ----
TW.GAM_temp_fall <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = tw(), data = subset(full_WG_df, SEASON == "FALL"))

TW.GAM_FALL_YEAR_df <- visreg::visreg(TW.GAM_temp_fall, data = subset(full_WG_df, SEASON == "FALL"), "YEAR", scale = "response",  plot = FALSE)

TW.GAM_temp_fall__df <- data.frame(YEAR = as.numeric(as.character(TW.GAM_FALL_YEAR_df$fit$YEAR)),
                                   Season = "FALL",
                                   MODEL = "Tweedie-GAM",
                                   fit = TW.GAM_FALL_YEAR_df$fit$visregFit,
                                   lwr = TW.GAM_FALL_YEAR_df$fit$visregLwr,
                                   upr = TW.GAM_FALL_YEAR_df$fit$visregUpr) 

  ## spring dataset model ----
TW.GAM_temp_spring <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = tw(), data = subset(full_WG_df, SEASON == "SPRING"))

TW.GAM_SPRING_YEAR_df <- visreg::visreg(TW.GAM_temp_spring, data = subset(full_WG_df, SEASON == "SPRING"), "YEAR", scale = "response",  plot = FALSE)

TW.GAM_temp_spring__df <- data.frame(YEAR = as.numeric(as.character(TW.GAM_SPRING_YEAR_df$fit$YEAR)),
                                   Season = "SPRING",
                                   MODEL = "Tweedie-GAM",
                                   fit = TW.GAM_SPRING_YEAR_df$fit$visregFit,
                                   lwr = TW.GAM_SPRING_YEAR_df$fit$visregLwr,
                                   upr = TW.GAM_SPRING_YEAR_df$fit$visregUpr) 

  ## annual dataset model ----
load("results/spatial model/summer.flounder/models objects/TW.GAM_dist_outside.Rdata")

TW.GAM_ANNUAL_YEAR_df <- visreg::visreg(TW.GAM_dist_outside, data = subset(full_WG_df, OWF == "OUTSIDE"), "YEAR", scale = "response",  plot = FALSE)

TW.GAM_temp_annual__df <- data.frame(YEAR = as.numeric(as.character(TW.GAM_ANNUAL_YEAR_df$fit$YEAR)),
                                     Season = "ANNUAL",
                                     MODEL = "Tweedie-GAM",
                                     fit = TW.GAM_ANNUAL_YEAR_df$fit$visregFit,
                                     lwr = TW.GAM_ANNUAL_YEAR_df$fit$visregLwr,
                                     upr = TW.GAM_ANNUAL_YEAR_df$fit$visregUpr) 

  ## combine and save ----
TW.GAM_temp_df <- rbind(TW.GAM_temp_fall__df, TW.GAM_temp_spring__df, TW.GAM_temp_annual__df) 
write.csv(TW.GAM_temp_df, "results/temporal model/summer.flounder/Tweedie.GAM_YEAR_Indices.csv", row.names = FALSE)

save(TW.GAM_temp_fall, full_WG_df, file = "results/temporal model/summer.flounder/models objects/TW.GAM_temp_fall.Rdata")
save(TW.GAM_temp_spring, full_WG_df, file = "results/temporal model/summer.flounder/models objects/TW.GAM_temp_spring.Rdata")

rm(list = setdiff(ls(), "full_df"))

# -------------------------------------------------------------------------------------------------------- #



# 2. delta GAM ---------------------------------------------------------------------------------------------

library(mgcv)

full_DELTA.GAM_df <- full_df %>%
  mutate(PRESENCE = ifelse(BIOMASS > 0, 1, 0)) %>%
  mutate(YEAR = as.factor(YEAR))

  ## 2.1 fall dataset model ----
DELTA.S1.GAM_temp_fall <- gam(PRESENCE ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, SEASON == "FALL"))
DELTA.S2.GAM_temp_fall <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = tw(), data = subset(full_DELTA.GAM_df, SEASON == "FALL" & PRESENCE == 1))

AIC(DELTA.S2.GAM_temp_fall)

# load(file = "results/temporal model/summer.flounder/models objects/DELTA.s1.GAM_temp_fall.Rdata")
# load(file = "results/temporal model/summer.flounder/models objects/DELTA.s2.GAM_temp_fall.Rdata")

    ### generate 2-stage predicted values
DELTA.GAM_PREDICTED_fall <- full_DELTA.GAM_df %>%
  filter(SEASON == "FALL") %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS, PRESENCE) %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_temp_fall, newdata = ., type = "response"),
         PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_temp_fall, newdata = ., type = "response"))

    ### generate GAM model based on final prediction for year effect - tw(), Gamma(link = "log"), or gaussian(link = "log")
DELTA.GAM_fit_year_fall <- gam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), data = DELTA.GAM_PREDICTED_fall, family = gaussian(link = "log"),
                               nthreads = parallel::detectCores()) 

AIC(DELTA.GAM_fit_year_fall) # gamma =  -4.260088, tw = NA, gaussian = -2215.16

save(DELTA.GAM_fit_year_fall, DELTA.GAM_PREDICTED_fall, file = "results/temporal model/summer.flounder/models objects/DELTA.GAM_year_effect_GAM_fall.Rdata")
# load("results/temporal model/summer.flounder/models objects/DELTA.GAM_year_effect_GAM_fall.Rdata")

DELTA.GAM_year_effect_fall <- visreg(DELTA.GAM_fit_year_fall, data = DELTA.GAM_PREDICTED_fall, "YEAR", scale = "response", plot = FALSE)

    ### save
DELTA.GAM_temp_fall_df <- data.frame(YEAR = as.numeric(as.character(DELTA.GAM_year_effect_fall$fit$YEAR)),
                                     Season = "FALL",
                                     MODEL = "Delta-GAM",
                                     fit = DELTA.GAM_year_effect_fall$fit$visregFit,
                                     lwr = DELTA.GAM_year_effect_fall$fit$visregLwr,
                                     upr = DELTA.GAM_year_effect_fall$fit$visregUpr)




  ## 2.2 spring dataset model ----

DELTA.S1.GAM_temp_spring <- gam(PRESENCE ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, SEASON == "SPRING"))
DELTA.S2.GAM_temp_spring <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = tw(), data = subset(full_DELTA.GAM_df, SEASON == "SPRING" & PRESENCE == 1))

AIC(DELTA.S2.GAM_temp_spring)

# load(file = "results/temporal model/summer.flounder/models objects/DELTA.s1.GAM_temp_spring.Rdata")
# load(file = "results/temporal model/summer.flounder/models objects/DELTA.s2.GAM_temp_spring.Rdata")

    ### generate 2-stage predicted values
DELTA.GAM_PREDICTED_spring <- full_DELTA.GAM_df %>%
  filter(SEASON == "SPRING") %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS, PRESENCE) %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_temp_spring, newdata = ., type = "response"),
         PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_temp_spring, newdata = ., type = "response"))

    ### generate GAM model based on final prediction for year effect - tw(), Gamma(link = "log"), or gaussian(link = "log")
DELTA.GAM_fit_year_spring <- gam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), data = DELTA.GAM_PREDICTED_spring, family = Gamma(link = "log"), 
                                 nthreads = parallel::detectCores()) 

AIC(DELTA.GAM_fit_year_spring) # gamma = -3089.262, tw = , gaussian = -2270.313

save(DELTA.GAM_fit_year_spring, DELTA.GAM_PREDICTED_spring, file = "results/temporal model/summer.flounder/models objects/DELTA.GAM_year_effect_GAM_spring.Rdata")
# load("results/temporal model/summer.flounder/models objects/DELTA.GAM_year_effect_GAM_spring.Rdata")

DELTA.GAM_year_effect_spring <- visreg(DELTA.GAM_fit_year_spring, data = DELTA.GAM_PREDICTED_spring, "YEAR", scale = "response", plot = FALSE)

    ### save
DELTA.GAM_temp_spring_df <- data.frame(YEAR = as.numeric(as.character(DELTA.GAM_year_effect_spring$fit$YEAR)),
                                     Season = "SPRING",
                                     MODEL = "Delta-GAM",
                                     fit = DELTA.GAM_year_effect_spring$fit$visregFit,
                                     lwr = DELTA.GAM_year_effect_spring$fit$visregLwr,
                                     upr = DELTA.GAM_year_effect_spring$fit$visregUpr)


  ## 2.3 annual dataset model ----

load("results/spatial model/summer.flounder/models objects/DELTA.S1.GAM_dist_outside.Rdata")
load("results/spatial model/summer.flounder/models objects/DELTA.S2.GAM_dist_outside.Rdata")

    ### generate 2-stage predicted values
DELTA.GAM_PREDICTED_annual <- full_DELTA.GAM_df %>%
  filter(OWF == "OUTSIDE") %>%
  add_column(SEASON = "ANNUAL") %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, BIOMASS, PRESENCE) %>%
  mutate(PREDICTED_PRESENCE = predict(DELTA.S1.GAM_dist_outside, newdata = ., type = "response"),
         PREDICTED_BIOMASS =  PREDICTED_PRESENCE * predict(DELTA.S2.GAM_dist_outside, newdata = ., type = "response")) %>%
  mutate(PREDICTED_BIOMASS = as.numeric(PREDICTED_BIOMASS))

    ### generate GAM model based on final prediction for year effect - tw(), Gamma(link = "log"), or gaussian(link = "log")

DELTA.GAM_fit_year_annual <- bam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), data = DELTA.GAM_PREDICTED_annual, family = gaussian(link = "log"), 
                                 nthreads = parallel::detectCores()) 

AIC(DELTA.GAM_fit_year_annual) # gamma = -3821.365, tw = NA, gaussian = -4546.359

save(DELTA.GAM_fit_year_annual, DELTA.GAM_PREDICTED_annual, file = "results/temporal model/summer.flounder/models objects/DELTA.GAM_year_effect_GAM_annual.Rdata")
# load("results/temporal model/summer.flounder/models objects/DELTA.GAM_year_effect_GAM_annual.Rdata")

DELTA.GAM_year_effect_annual <- visreg(DELTA.GAM_fit_year_annual, data = DELTA.GAM_PREDICTED_annual, "YEAR", scale = "response", plot = FALSE)

    ### save
DELTA.GAM_temp_annual_df <- data.frame(YEAR = as.numeric(as.character(DELTA.GAM_year_effect_annual$fit$YEAR)),
                                       Season = "ANNUAL",
                                       MODEL = "Delta-GAM",
                                       fit = DELTA.GAM_year_effect_annual$fit$visregFit,
                                       lwr = DELTA.GAM_year_effect_annual$fit$visregLwr,
                                       upr = DELTA.GAM_year_effect_annual$fit$visregUpr)



  ## 2.4 combine and save ----
DELTA.GAM_temp_df <- rbind(DELTA.GAM_temp_fall_df, DELTA.GAM_temp_spring_df, DELTA.GAM_temp_annual_df) 
write.csv(DELTA.GAM_temp_df, "results/temporal model/summer.flounder/Delta.GAM_YEAR_Indices.csv", row.names = FALSE)

save(DELTA.S1.GAM_temp_fall, full_DELTA.GAM_df, file = "results/temporal model/summer.flounder/models objects/DELTA.s1.GAM_temp_fall.Rdata")
save(DELTA.S2.GAM_temp_fall, full_DELTA.GAM_df, file = "results/temporal model/summer.flounder/models objects/DELTA.s2.GAM_temp_fall.Rdata")
save(DELTA.S1.GAM_temp_spring, full_DELTA.GAM_df, file = "results/temporal model/summer.flounder/models objects/DELTA.s1.GAM_temp_spring.Rdata")
save(DELTA.S2.GAM_temp_spring, full_DELTA.GAM_df, file = "results/temporal model/summer.flounder/models objects/DELTA.s2.GAM_temp_spring.Rdata")

rm(list = setdiff(ls(), "full_df"))

# -------------------------------------------------------------------------------------------------------- #


# 3. RANDOM FOREST ---------------------------------------------------------------------------------------------

library(randomForest)
library(pdp) # to extract RF model results

full_RF_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))

  ## fall dataset model ----

set.seed(2) # ensure reproducibility
RF_temp_fall <- randomForest(BIOMASS ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = subset(full_RF_df, SEASON == "FALL"), mtry = 2, ntree = 1000)

    ### Create Partial Dependence Plot Data
RF_fall_year_df <- partial(RF_temp_fall, pred.var = "YEAR", grid.resolution = 100, plot = FALSE, train = subset(full_RF_df, SEASON == "FALL")) 

#     ### Bootstrapping for Confidence Intervals
# set.seed(100) # Ensure reproducibility
# boot_fits <- matrix(NA, nrow = 100, ncol = length(unique(full_RF_df$YEAR))) # Initialize matrices to store bootstrap results
# 
# for (i in 1:100) { # Perform bootstrapping
#   boot_indices <- sample(nrow(subset(full_RF_df, SEASON == "FALL")), replace = TRUE)
#   boot_data <- subset(full_RF_df, SEASON == "FALL")[boot_indices, ]
#   
#   boot_model <- randomForest(NUMBER ~ ., data = boot_data, importance = TRUE)
#   boot_pdp <- partial(boot_model, pred.var = "YEAR", grid.resolution = length(year_values), 
#                       plot = FALSE, train = boot_data)
#   boot_fits[i, ] <- boot_pdp$yhat
#   print(paste("iter",i,"finished"))
# }; remove(boot_indices, boot_data, boot_model, boot_pdp,i)

RF_temp_fall_df <- data.frame(YEAR = as.numeric(as.character(RF_fall_year_df$YEAR)),
                                       Season = "FALL",
                                       MODEL = "Random Forest",
                                       fit = RF_fall_year_df$yhat, #colMeans(boot_fits),
                                       lwr = NA, #apply(boot_fits, 2, quantile, probs = 0.025),
                                       upr = NA) #apply(boot_fits, 2, quantile, probs = 0.975))


## spring dataset model ----

set.seed(2) # ensure reproducibility
RF_temp_spring <- randomForest(BIOMASS ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = subset(full_RF_df, SEASON == "SPRING"), mtry = 2, ntree = 1000)

### Create Partial Dependence Plot Data
RF_spring_year_df <- partial(RF_temp_spring, pred.var = "YEAR", grid.resolution = 100, plot = FALSE, train = subset(full_RF_df, SEASON == "SPRING")) 

RF_temp_spring_df <- data.frame(YEAR = as.numeric(as.character(RF_spring_year_df$YEAR)),
                              Season = "SPRING",
                              MODEL = "Random Forest",
                              fit = RF_spring_year_df$yhat, #colMeans(boot_fits),
                              lwr = NA, #apply(boot_fits, 2, quantile, probs = 0.025),
                              upr = NA) #apply(boot_fits, 2, quantile, probs = 0.975))


## annual dataset model ----

load("results/spatial model/summer.flounder/models objects/RF_dist_outside.Rdata")

RF_annual_df <- partial(RF_dist_outside, pred.var = "YEAR", grid.resolution = 100, plot = FALSE, train = subset(full_RF_df, OWF == "OUTSIDE"))

RF_temp_annual_df <- data.frame(YEAR = as.numeric(as.character(RF_annual_df$YEAR)),
                                Season = "ANNUAL",
                                MODEL = "Random Forest",
                                fit = RF_annual_df$yhat, #colMeans(boot_fits),
                                lwr = NA, #apply(boot_fits, 2, quantile, probs = 0.025),
                                upr = NA) #apply(boot_fits, 2, quantile, probs = 0.975))


## combine and save ----
RF_temp_df <- rbind(RF_temp_fall_df, RF_temp_spring_df, RF_temp_annual_df) 
write.csv(RF_temp_df, "results/temporal model/summer.flounder/RF_YEAR_Indices.csv", row.names = FALSE)

save(RF_temp_fall, full_RF_df, file = "results/temporal model/summer.flounder/models objects/RF_temp_fall.Rdata")
save(RF_temp_spring, full_RF_df, file = "results/temporal model/summer.flounder/models objects/RF_temp_spring.Rdata")

rm(list = setdiff(ls(), "full_df"))


# -------------------------------------------------------------------------------------------------------- #



# 4. Bayesian Additive Regression Trees ---------------------------------------------------------------------------------------------
library(dbarts)
library(pdp) # to extract BART model results

full_BART_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR)) %>%
  mutate(PRESENCE = ifelse(BIOMASS > 0, 1, 0))

  ## 4.1 fall dataset model ----

fall_BART_df <- subset(full_BART_df, SEASON == "FALL")


    ### 4.1.1 fit presence ----

BART_presence_fall <- dbarts::bart(x.train = fall_BART_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], y.train = fall_BART_df[, c("PRESENCE")], keeptrees = TRUE)

invisible(BART_presence_fall$fit$state) # very important if you want to use the models later

cutoff <- InformationValue::optimalCutoff(fall_BART_df$PRESENCE, fitted(BART_presence_fall)) # determine cutoff for presence/absence


    ### 4.1.2 fit biomass for data that presence = 1 ----

fall_BART_PRESENCE_df <- subset(fall_BART_df, PRESENCE == 1)

BART_BIOMASS_fall <- dbarts::bart(x.train = fall_BART_PRESENCE_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], y.train = fall_BART_PRESENCE_df[, c("BIOMASS")], keeptrees = TRUE)

invisible(BART_BIOMASS_fall$fit$state) # very important if you want to use the models later


    ### 4.1.3 predict for the entire region dataset ----

      #### presence 
BART_presence_fall_predict <- dbarts:::predict.bart(BART_presence_fall , newdata = fall_BART_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_presence_fall_predict <- apply(BART_presence_fall_predict, 2, median)


      #### biomass
temp_presence_df <- fall_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_fall_predict <- dbarts:::predict.bart(BART_BIOMASS_fall , newdata = temp_presence_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_biomass_fall_predict <- apply(BART_biomass_fall_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_fall_predict)


    ### 4.1.4 combine biomass and presence from prediction ----

BART_comp_fall_dataset_df <- fall_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Bayesian Additive Regression Trees", DATASET = "FALL") %>%
  mutate(PREDICTED_PRESENCE = BART_presence_fall_predict) %>%
  mutate(PREDICTED_PRESENCE = ifelse(PREDICTED_PRESENCE > cutoff, 1, 0)) %>%  # use cutoff to categorize probability into 0 and 1
  left_join(temp_presence_df) %>%
  mutate(PREDICTED_BIOMASS = ifelse(is.na(PREDICTED_BIOMASS), 0, PREDICTED_BIOMASS)) %>%  # fill in the absent data with 0
  mutate(DELTA_PREDICTION = PREDICTED_PRESENCE * PREDICTED_BIOMASS) %>%
  mutate(DELTA_PREDICTION = ifelse(DELTA_PREDICTION < 0, 0, DELTA_PREDICTION)) %>% # post-processing to avoid negative BIOMASS
  select(-c(PREDICTED_BIOMASS)) %>%
  rename(PREDICTED_BIOMASS = DELTA_PREDICTION)



    ### 4.1.5 fit GAM for year effect ---- 

     ####  - tw(), Gamma(link = "log"), or gaussian(link = "log")
BART.GAM_fit_year_fall <- bam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), data = BART_comp_fall_dataset_df, family = tw(), 
                                 nthreads = parallel::detectCores()) 

AIC(BART.GAM_fit_year_fall) # gamma = NA, tw = 8991.962, gaussian = NA

save(BART.GAM_fit_year_fall, BART_comp_fall_dataset_df, file = "results/temporal model/summer.flounder/models objects/BART_year_effect_GAM_fall.Rdata")
# load("results/temporal model/summer.flounder/models objects/BART_year_effect_GAM_fall.Rdata")

BART.GAM_year_effect_fall <- visreg(BART.GAM_fit_year_fall, data = BART_comp_fall_dataset_df, "YEAR", scale = "response", plot = FALSE)

BART_temp_fall_df <- data.frame(YEAR = as.numeric(as.character(BART.GAM_year_effect_fall$fit$YEAR)),
                                Season = "FALL",
                                MODEL = "Bayesian Additive Regression Trees",
                                fit = BART.GAM_year_effect_fall$fit$visregFit, 
                                lwr = BART.GAM_year_effect_fall$fit$visregLwr,
                                upr = BART.GAM_year_effect_fall$fit$visregUpr) 



  ## 4.2 spring dataset model ----

spring_BART_df <- subset(full_BART_df, SEASON == "SPRING")


    ### 4.2.1 fit presence ----

BART_presence_spring <- dbarts::bart(x.train = spring_BART_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], y.train = spring_BART_df[, c("PRESENCE")], keeptrees = TRUE)

invisible(BART_presence_spring$fit$state) # very important if you want to use the models later

cutoff <- InformationValue::optimalCutoff(spring_BART_df$PRESENCE, fitted(BART_presence_spring)) # determine cutoff for presence/absence


    ### 4.2.2 fit biomass for data that presence = 1 ----

spring_BART_PRESENCE_df <- subset(spring_BART_df, PRESENCE == 1)

BART_BIOMASS_spring <- dbarts::bart(x.train = spring_BART_PRESENCE_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], y.train = spring_BART_PRESENCE_df[, c("BIOMASS")], keeptrees = TRUE)

invisible(BART_BIOMASS_spring$fit$state) # very important if you want to use the models later


    ### 4.2.3 predict for the entire region dataset ----

      #### presence 
BART_presence_spring_predict <- dbarts:::predict.bart(BART_presence_spring , newdata = spring_BART_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_presence_spring_predict <- apply(BART_presence_spring_predict, 2, median)


      #### biomass
temp_presence_df <- spring_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_spring_predict <- dbarts:::predict.bart(BART_BIOMASS_spring , newdata = temp_presence_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_biomass_spring_predict <- apply(BART_biomass_spring_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_spring_predict)


    ### 4.2.4 combine biomass and presence from prediction ----

BART_comp_spring_dataset_df <- spring_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Bayesian Additive Regression Trees", DATASET = "SPRING") %>%
  mutate(PREDICTED_PRESENCE = BART_presence_spring_predict) %>%
  mutate(PREDICTED_PRESENCE = ifelse(PREDICTED_PRESENCE > cutoff, 1, 0)) %>%  # use cutoff to categorize probability into 0 and 1
  left_join(temp_presence_df) %>%
  mutate(PREDICTED_BIOMASS = ifelse(is.na(PREDICTED_BIOMASS), 0, PREDICTED_BIOMASS)) %>%  # fill in the absent data with 0
  mutate(DELTA_PREDICTION = PREDICTED_PRESENCE * PREDICTED_BIOMASS) %>%
  mutate(DELTA_PREDICTION = ifelse(DELTA_PREDICTION < 0, 0, DELTA_PREDICTION)) %>% # post-processing to avoid negative BIOMASS
  select(-c(PREDICTED_BIOMASS)) %>%
  rename(PREDICTED_BIOMASS = DELTA_PREDICTION)



    ### 4.2.5 fit GAM for year effect ---- 

        ####  - tw(), Gamma(link = "log"), or gaussian(link = "log")
BART.GAM_fit_year_spring <- bam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), data = BART_comp_spring_dataset_df, family = tw(), 
                                 nthreads = parallel::detectCores()) 

AIC(BART.GAM_fit_year_spring) # gamma = NA, tw = 5967.085, gaussian = NA

save(BART.GAM_fit_year_spring, BART_comp_spring_dataset_df, file = "results/temporal model/summer.flounder/models objects/BART_year_effect_GAM_spring.Rdata")
# load("results/temporal model/summer.flounder/models objects/BART_year_effect_GAM_spring.Rdata")

BART.GAM_year_effect_spring <- visreg(BART.GAM_fit_year_spring, data = BART_comp_spring_dataset_df, "YEAR", scale = "response", plot = FALSE)

BART_temp_spring_df <- data.frame(YEAR = as.numeric(as.character(BART.GAM_year_effect_spring$fit$YEAR)),
                                Season = "SPRING",
                                MODEL = "Bayesian Additive Regression Trees",
                                fit = BART.GAM_year_effect_spring$fit$visregFit, 
                                lwr = BART.GAM_year_effect_spring$fit$visregLwr,
                                upr = BART.GAM_year_effect_spring$fit$visregUpr) 



  ## 4.3 Annual dataset model ----

load("results/spatial model/summer.flounder/models objects/BART_presence_outside.Rdata")
load("results/spatial model/summer.flounder/models objects/BART_BIOMASS_outside.Rdata")

invisible(BART_presence_outside$fit$state)
invisible(BART_BIOMASS_outside$fit$state) 

cutoff <- InformationValue::optimalCutoff(outside_BART_df$PRESENCE, fitted(BART_presence_outside))

    ### 4.3.1 predict for the entire region dataset ----

      #### presence 
BART_presence_annual_predict <- dbarts:::predict.bart(BART_presence_outside, newdata = outside_BART_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_presence_annual_predict <- apply(BART_presence_annual_predict, 2, median)

      #### biomass
temp_presence_df <- outside_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_annual_predict <- dbarts:::predict.bart(BART_BIOMASS_outside , newdata = temp_presence_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_biomass_annual_predict <- apply(BART_biomass_annual_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_annual_predict)


    ### 4.3.2 combine biomass and presence from prediction ----

BART_comp_annual_dataset_df <- outside_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Bayesian Additive Regression Trees", DATASET = "ANNUAL") %>%
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
BART.GAM_fit_year_annual <- bam(PREDICTED_BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), data = BART_comp_annual_dataset_df, family = tw(), 
                                 nthreads = parallel::detectCores()) 

AIC(BART.GAM_fit_year_annual) # gamma = NA, tw = 14691.9, gaussian = NA

save(BART.GAM_fit_year_annual, BART_comp_annual_dataset_df, file = "results/temporal model/summer.flounder/models objects/BART_year_effect_GAM_annual.Rdata")
# load("results/temporal model/summer.flounder/models objects/BART_year_effect_GAM_annual.Rdata")

BART.GAM_year_effect_annual <- visreg(BART.GAM_fit_year_annual, data = BART_comp_annual_dataset_df, "YEAR", scale = "response", plot = FALSE)

BART_temp_annual_df <- data.frame(YEAR = as.numeric(as.character(BART.GAM_year_effect_annual$fit$YEAR)),
                                  Season = "ANNUAL",
                                  MODEL = "Bayesian Additive Regression Trees",
                                  fit = BART.GAM_year_effect_annual$fit$visregFit, 
                                  lwr = BART.GAM_year_effect_annual$fit$visregLwr,
                                  upr = BART.GAM_year_effect_annual$fit$visregUpr) 


  ## 4.4 combine and save ----

BART_temp_df <- rbind(BART_temp_fall_df, BART_temp_spring_df, BART_temp_annual_df) %>%
  mutate(fit = ifelse(fit < 0, 0, fit))  # post-processing to ensure minimal value should be 0

write.csv(BART_temp_df, "results/temporal model/summer.flounder/BART_YEAR_Indices.csv", row.names = FALSE)

save(BART_presence_fall, fall_BART_df, file = "results/temporal model/summer.flounder/models objects/BART_presence_fall.Rdata")
save(BART_BIOMASS_fall, fall_BART_df, file = "results/temporal model/summer.flounder/models objects/BART_BIOMASS_fall.Rdata")
save(BART_presence_spring, spring_BART_df, file = "results/temporal model/summer.flounder/models objects/BART_presence_spring.Rdata")
save(BART_BIOMASS_spring, spring_BART_df, file = "results/temporal model/summer.flounder/models objects/BART_BIOMASS_spring.Rdata")

rm(list = setdiff(ls(), "full_df"))




# -------------------------------------------------------------------------------------------------------- #



# 5. VAST ---------------------------------------------------------------------------------------------

library(VAST)
library(splines)
library(Metrics)

  ## prepare data ----
full_VAST_df <- full_df %>%   # convert YEAR into a factor
  add_column(AreaSwept_km2 = 0.024, # 24000 m2 from doi:10.1093/icesjms/fsv166 and BTS protocol
             Vessel = NA, Pass = 0) %>% # required by the model, don't know why  
  rename(Year = YEAR, Lat = LAT, Lon = LON) # the exact spelling is required by VAST

  ## prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)],
                                                # if use "Georges_Bank", "Mid_Atlantic_Bight", there will be 436 rows missing for summer flounder
                                                # all outside OWF, not affecting the gap analysis
                                                DirPath = "results/spatial model/summer.flounder/models objects/VAST/")

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

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

full_VAST_df$BOTTEMP_Std <- scale(full_VAST_df$BOTTEMP, center = TRUE, scale = TRUE)[,1]
full_VAST_df$AVGDEPTH_Std <- scale(full_VAST_df$AVGDEPTH, center = TRUE, scale = TRUE)[,1]

    ### design formulas for the two variables
p1_formula = ~ bs(BOTTEMP_Std, degree = 2) + bs(AVGDEPTH_Std, degree = 2)
p2_formula = ~ bs(BOTTEMP_Std, degree = 2) + bs(AVGDEPTH_Std, degree = 2)



  ## fall dataset model ----
fall_VAST_df <- subset(full_VAST_df, SEASON == "FALL")

    ### standardize the covariate to have mean 0 and standard deviation 1.0 as suggested by VAST
fall_VAST_df$BOTTEMP_Std <- scale(fall_VAST_df$BOTTEMP, center = TRUE, scale = TRUE)[,1]
fall_VAST_df$AVGDEPTH_Std <- scale(fall_VAST_df$AVGDEPTH, center = TRUE, scale = TRUE)[,1]



    ### loop over multiple error/link options to identify the best model ----

scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description

# i = 6 is the optimal model
for (i in 1:nrow(scenario_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_temp_fall <- tryCatch(
    # attempt to fit the model
    {
      VAST_temp_fall <- fit_model(settings = settings,
                                  input_grid = extrapolation_region$Data_Extrap, 
                                  Lat_i = fall_VAST_df$Lat,
                                  Lon_i = fall_VAST_df$Lon,
                                  t_i = as.numeric(fall_VAST_df$Year), 
                                  b_i = fall_VAST_df$BIOMASS,
                                  a_i = as_units(fall_VAST_df$AreaSwept_km2, "km^2"),
                                  # covariate_data = full_VAST_df, # use full_VAST_df here to avoid missing year (2020) for fall
                                  # X1_formula = p1_formula,
                                  # X2_formula = p2_formula,
                                  working_dir = "results/temporal model/summer.flounder/models objects/VAST/")
      VAST_temp_fall # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_fall)) { next }
  if (length(VAST_temp_fall$Report) == 1) {
    if (VAST_temp_fall$Report == "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_fall$parameter_estimates$AIC
  
  save(VAST_temp_fall, file = paste0("results/temporal model/summer.flounder/models objects/VAST_fall_",
                                                    scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # save(VAST_temp_fall, full_VAST_df, file = paste0("results/temporal model/summer.flounder/models objects/VAST_temp_fall.Rdata"))
  # load(file = paste0("results/temporal model/summer.flounder/models objects/VAST_temp_fall.Rdata"))
  
}

write.csv(scenario_df, "results/temporal model/summer.flounder/models objects/VAST/VAST_fall_model_selection.csv")

 ### save the fitted value ----

load("results/temporal model/summer.flounder/models objects/VAST_fall_model.Rdata")

VAST_fall_data <- plot(VAST_temp_fall, working_dir = "results/temporal model/summer.flounder/models objects/VAST/figure/")

VAST_temp_fall_df <- data.frame(YEAR = VAST_fall_data$Index$Table$Time,
                                Season = "FALL",
                                MODEL = "VAST",
                                fit = VAST_fall_data$Index$Table$Estimate, 
                                lwr = VAST_fall_data$Index$Table$Estimate - 1.96 * VAST_fall_data$Index$Table$`Std. Error for Estimate`, 
                                upr = VAST_fall_data$Index$Table$Estimate + 1.96 * VAST_fall_data$Index$Table$`Std. Error for Estimate`) 


  ## spring dataset model ----
spring_VAST_df <- subset(full_VAST_df, SEASON == "SPRING")

    ### standardize the covariate to have mean 0 and standard deviation 1.0 as suggested by VAST
spring_VAST_df$BOTTEMP_Std <- scale(spring_VAST_df$BOTTEMP, center = TRUE, scale = TRUE)[,1]
spring_VAST_df$AVGDEPTH_Std <- scale(spring_VAST_df$AVGDEPTH, center = TRUE, scale = TRUE)[,1]

    ### loop over multiple error/link options to identify the best model ----

scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description

# i = 6 is the optimal model
for (i in 1:nrow(scenario_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_temp_spring <- tryCatch(
    # attempt to fit the model
    {
      VAST_temp_spring <- fit_model(settings = settings,
                                  input_grid = extrapolation_region$Data_Extrap, 
                                  Lat_i = spring_VAST_df$Lat,
                                  Lon_i = spring_VAST_df$Lon,
                                  t_i = as.numeric(spring_VAST_df$Year), 
                                  b_i = spring_VAST_df$BIOMASS,
                                  a_i = as_units(spring_VAST_df$AreaSwept_km2, "km^2"),
                                  covariate_data = full_VAST_df, # use full_VAST_df here to avoid missing year (2020) for fall
                                  # X1_formula = p1_formula,
                                  # X2_formula = p2_formula,
                                  working_dir = "results/temporal model/summer.flounder/models objects/VAST/")
      VAST_temp_spring # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_spring)) { next }
  if (length(VAST_temp_spring$Report) == 1) {
    if (VAST_temp_spring$Report != "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_spring$parameter_estimates$AIC

  save(VAST_temp_spring, file = paste0("results/temporal model/summer.flounder/models objects/VAST_spring_",
                                                   scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # save(VAST_temp_spring, file = paste0("results/temporal model/summer.flounder/models objects/VAST_temp_spring.Rdata"))
  
  # load(file = paste0("results/temporal model/summer.flounder/models objects/VAST_temp_spring.Rdata"))
}

write.csv(scenario_df, "results/temporal model/summer.flounder/models objects/VAST/VAST_spring_model_selection.csv")


### save the fitted value ----

load("results/temporal model/summer.flounder/models objects/VAST_spring_model.Rdata")

VAST_spring_data <- plot(VAST_temp_spring, working_dir = "results/temporal model/summer.flounder/models objects/VAST/figure/")

VAST_temp_spring_df <- data.frame(YEAR = VAST_spring_data$Index$Table$Time,
                                Season = "SPRING",
                                MODEL = "VAST",
                                fit = VAST_spring_data$Index$Table$Estimate, 
                                lwr = VAST_spring_data$Index$Table$Estimate - 1.96 * VAST_spring_data$Index$Table$`Std. Error for Estimate`, 
                                upr = VAST_spring_data$Index$Table$Estimate + 1.96 * VAST_spring_data$Index$Table$`Std. Error for Estimate`) 

## annual dataset model ----

load("results/spatial model/summer.flounder/models objects/VAST_outside_model.Rdata")

VAST_annual_data <- plot(VAST_outside_model, working_dir = "results/temporal model/summer.flounder/models objects/VAST/figure/")

VAST_temp_annual_df <- data.frame(YEAR = VAST_annual_data$Index$Table$Time,
                                  Season = "ANNUAL",
                                  MODEL = "VAST",
                                  fit = VAST_annual_data$Index$Table$Estimate, 
                                  lwr = VAST_annual_data$Index$Table$Estimate - 1.96 * VAST_annual_data$Index$Table$`Std. Error for Estimate`, 
                                  upr = VAST_annual_data$Index$Table$Estimate + 1.96 * VAST_annual_data$Index$Table$`Std. Error for Estimate`) 

  ## combine and save ----

VAST_temp_df <- rbind(VAST_temp_fall_df, VAST_temp_spring_df, VAST_temp_annual_df) 
write.csv(VAST_temp_df, "results/temporal model/summer.flounder/VAST_YEAR_Indices.csv", row.names = FALSE)

rm(list = setdiff(ls(), "full_df"))



# 7. delta GAM (simple original version) ---------------------------------------------------------------------------------------------

library(mgcv)

full_DELTA.GAM_df <- full_df %>%
  mutate(PRESENCE = ifelse(BIOMASS > 0, 1, 0)) %>%
  mutate(YEAR = as.factor(YEAR))

## 2.1 fall dataset model ----
DELTA.S1.GAM_temp_fall <- gam(PRESENCE ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, SEASON == "FALL"))
DELTA.S2.GAM_temp_fall <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = tw(), data = subset(full_DELTA.GAM_df, SEASON == "FALL" & PRESENCE == 1))

AIC(DELTA.S2.GAM_temp_fall)

# load(file = "results/temporal model/summer.flounder/models objects/DELTA.s1.GAM_temp_fall.Rdata")
# load(file = "results/temporal model/summer.flounder/models objects/DELTA.s2.GAM_temp_fall.Rdata")

### generate 2-stage predicted values
DELTA.s1.GAM_fall_df <- visreg(DELTA.S1.GAM_temp_fall, data = subset(full_DELTA.GAM_df, SEASON == "FALL"), "YEAR", scale = "response", plot = FALSE)
DELTA.s2.GAM_fall_df <- visreg(DELTA.S2.GAM_temp_fall, data = subset(full_DELTA.GAM_df, SEASON == "FALL" & PRESENCE == 1), "YEAR", scale = "response", plot = FALSE)

DELTA.GAM_temp_fall_df <- data.frame(YEAR = as.numeric(as.character(DELTA.s1.GAM_fall_df$fit$YEAR)),
                                     Season = "FALL",
                                     MODEL = "Delta-GAM",
                                     fit = DELTA.s1.GAM_fall_df$fit$visregFit * DELTA.s2.GAM_fall_df$fit$visregFit,
                                     lwr = DELTA.s1.GAM_fall_df$fit$visregLwr * DELTA.s2.GAM_fall_df$fit$visregLwr,
                                     upr = DELTA.s1.GAM_fall_df$fit$visregUpr * DELTA.s2.GAM_fall_df$fit$visregUpr)



## 2.2 spring dataset model ----

DELTA.S1.GAM_temp_spring <- gam(PRESENCE ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = binomial(link = "logit"), data = subset(full_DELTA.GAM_df, SEASON == "SPRING"))
DELTA.S2.GAM_temp_spring <- gam(BIOMASS ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), family = tw(), data = subset(full_DELTA.GAM_df, SEASON == "SPRING" & PRESENCE == 1))

AIC(DELTA.S2.GAM_temp_spring)

# load(file = "results/temporal model/summer.flounder/models objects/DELTA.s1.GAM_temp_spring.Rdata")
# load(file = "results/temporal model/summer.flounder/models objects/DELTA.s2.GAM_temp_spring.Rdata")

### generate 2-stage predicted values
DELTA.s1.GAM_spring_df <- visreg(DELTA.S1.GAM_temp_spring, data = subset(full_DELTA.GAM_df, SEASON == "SPRING"), "YEAR", scale = "response", plot = FALSE)
DELTA.s2.GAM_spring_df <- visreg(DELTA.S2.GAM_temp_spring, data = subset(full_DELTA.GAM_df, SEASON == "SPRING" & PRESENCE == 1), "YEAR", scale = "response", plot = FALSE)

DELTA.GAM_temp_spring_df <- data.frame(YEAR = as.numeric(as.character(DELTA.s1.GAM_spring_df$fit$YEAR)),
                                       Season = "SPRING",
                                       MODEL = "Delta-GAM",
                                       fit = DELTA.s1.GAM_spring_df$fit$visregFit * DELTA.s2.GAM_spring_df$fit$visregFit,
                                       lwr = DELTA.s1.GAM_spring_df$fit$visregLwr * DELTA.s2.GAM_spring_df$fit$visregLwr,
                                       upr = DELTA.s1.GAM_spring_df$fit$visregUpr * DELTA.s2.GAM_spring_df$fit$visregUpr)


## 2.3 annual dataset model ----

load("results/spatial model/summer.flounder/models objects/DELTA.S1.GAM_dist_outside.Rdata")
load("results/spatial model/summer.flounder/models objects/DELTA.S2.GAM_dist_outside.Rdata")

### generate 2-stage predicted values
DELTA.s1.GAM_ANNUAL_df <- visreg::visreg(DELTA.S1.GAM_dist_outside, data = subset(full_DELTA.GAM_df, OWF == "OUTSIDE"), "YEAR", scale = "response",  plot = FALSE)
DELTA.s2.GAM_ANNUAL_df <- visreg::visreg(DELTA.S2.GAM_dist_outside, data = subset(full_DELTA.GAM_df, OWF == "OUTSIDE" & PRESENCE == 1), "YEAR", scale = "response", plot = FALSE)


DELTA.GAM_temp_annual_df <- data.frame(YEAR = as.numeric(as.character(DELTA.s1.GAM_ANNUAL_df$fit$YEAR)),
                                       Season = "ANNUAL",
                                       MODEL = "Delta-GAM",
                                       fit = DELTA.s1.GAM_ANNUAL_df$fit$visregFit * DELTA.s2.GAM_ANNUAL_df$fit$visregFit,
                                       lwr = DELTA.s1.GAM_ANNUAL_df$fit$visregLwr * DELTA.s2.GAM_ANNUAL_df$fit$visregLwr,
                                       upr = DELTA.s1.GAM_ANNUAL_df$fit$visregUpr * DELTA.s2.GAM_ANNUAL_df$fit$visregUpr)


## 2.4 combine and save ----
DELTA.GAM_temp_df <- rbind(DELTA.GAM_temp_fall_df, DELTA.GAM_temp_spring_df, DELTA.GAM_temp_annual_df) 
write.csv(DELTA.GAM_temp_df, "results/temporal model/summer.flounder/simple.Delta.GAM_YEAR_Indices.csv", row.names = FALSE)

save(DELTA.S1.GAM_temp_fall, full_DELTA.GAM_df, file = "results/temporal model/summer.flounder/models objects/DELTA.s1.GAM_temp_fall.Rdata")
save(DELTA.S2.GAM_temp_fall, full_DELTA.GAM_df, file = "results/temporal model/summer.flounder/models objects/DELTA.s2.GAM_temp_fall.Rdata")
save(DELTA.S1.GAM_temp_spring, full_DELTA.GAM_df, file = "results/temporal model/summer.flounder/models objects/DELTA.s1.GAM_temp_spring.Rdata")
save(DELTA.S2.GAM_temp_spring, full_DELTA.GAM_df, file = "results/temporal model/summer.flounder/models objects/DELTA.s2.GAM_temp_spring.Rdata")

rm(list = setdiff(ls(), "full_df"))


library(dbarts)

full_BART_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR))  %>%
  mutate(PRESENCE = ifelse(BIOMASS > 0, 1, 0))


# for spatial modeling ----

## 4.1 full dataset model ----

### 4.1.1 fit presence ----

BART_presence_full <- dbarts::bart(x.train = full_BART_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], y.train = full_BART_df[, c("PRESENCE")], keeptrees = TRUE)

invisible(BART_presence_full$fit$state) # very important if you want to use the models later

cutoff <- InformationValue::optimalCutoff(full_BART_df$PRESENCE, fitted(BART_presence_full)) # determine cutoff for presence/absence




### 4.1.2 fit biomass for data that presence = 1 ----

full_BART_PRESENCE_df <- subset(full_BART_df, PRESENCE == 1)

BART_BIOMASS_full <- dbarts::bart(x.train = full_BART_PRESENCE_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], y.train = full_BART_PRESENCE_df[, c("BIOMASS")], keeptrees = TRUE)

invisible(BART_BIOMASS_full$fit$state) # very important if you want to use the models later




### 4.1.3 predict for the entire region dataset ----

#### presence 
BART_presence_full_predict <- dbarts:::predict.bart(BART_presence_full , newdata = full_BART_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_presence_full_predict <- apply(BART_presence_full_predict, 2, median)


#### biomass
temp_presence_df <- full_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_full_predict <- dbarts:::predict.bart(BART_BIOMASS_full , newdata = temp_presence_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_biomass_full_predict <- apply(BART_biomass_full_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_full_predict)



### 4.1.4 combine biomass and presence from prediction ----

BART_comp_full_dataset_df <- full_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Bayesian Additive Regression Trees", DATASET = "FULL") %>%
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

BART_presence_outside <- dbarts::bart(x.train = outside_BART_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], y.train = outside_BART_df[, c("PRESENCE")], keeptrees = TRUE)

invisible(BART_presence_outside$fit$state) # very important if you want to use the models later

cutoff <- InformationValue::optimalCutoff(outside_BART_df$PRESENCE, fitted(BART_presence_outside)) # determine cutoff for presence/absence




### 4.2.2 fit biomass for data that presence = 1 ----

outside_BART_PRESENCE_df <- subset(outside_BART_df, PRESENCE == 1)

BART_BIOMASS_outside <- dbarts::bart(x.train = outside_BART_PRESENCE_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], y.train = outside_BART_PRESENCE_df[, c("BIOMASS")], keeptrees = TRUE)

invisible(BART_BIOMASS_outside$fit$state) # very important if you want to use the models later




### 4.2.3 predict for the entire region dataset ----

#### presence 
BART_presence_outside_predict <- dbarts:::predict.bart(BART_presence_outside , newdata = full_BART_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_presence_outside_predict <- apply(BART_presence_outside_predict, 2, median)



#### biomass
temp_presence_df <- full_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  filter(PRESENCE == 1)

BART_biomass_outside_predict <- dbarts:::predict.bart(BART_BIOMASS_outside , newdata = temp_presence_df[, c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")], type = "bart")
BART_biomass_outside_predict <- apply(BART_biomass_outside_predict, 2, median)

temp_presence_df <- temp_presence_df %>%
  mutate(PREDICTED_BIOMASS = BART_biomass_outside_predict)



### 4.2.4 combine biomass and presence from prediction ----

BART_comp_outside_dataset_df <- full_BART_df %>%
  select(OWF, YEAR, LAT, BOTTEMP, AVGDEPTH, PRESENCE, BIOMASS) %>%
  add_column(SPECIES = "SUMMER FLOUNDER", MODEL = "Bayesian Additive Regression Trees", DATASET = "OUTSIDE") %>%
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
write.csv(BART_comp, "results/spatial model/summer.flounder/BART_BIOMASS_comparison.csv", row.names = FALSE)

save(BART_presence_full, full_BART_df, file = "results/spatial model/summer.flounder/models objects/BART_presence_full.Rdata")
save(BART_BIOMASS_full, full_BART_df, file = "results/spatial model/summer.flounder/models objects/BART_BIOMASS_full.Rdata")
save(BART_presence_outside, full_BART_df, file = "results/spatial model/summer.flounder/models objects/BART_dist_outside.Rdata")
save(BART_BIOMASS_outside, full_BART_df, file = "results/spatial model/summer.flounder/models objects/BART_BIOMASS_outside.Rdata")

rm(list = setdiff(ls(), c("full_df", "data_partition_list")))



# for temporal modeling (year effect) ----

# Year_presence_df <- partial(BART_presence, pred.var = "YEAR", plot = FALSE, train = full_BART_df) 
# Year_BIOMASS_df <- partial(BART_BIOMASS, pred.var = "YEAR", plot = FALSE, train = full_BART_df) 

# fit extra GAM to get the year effect

combined_gam <- gam(DELTA_PREDICTION ~ YEAR + s(LAT) + s(BOTTEMP) + s(AVGDEPTH), data = full_BART_df, family = tw())

year_effect_df <- visreg(combined_gam, data = full_BART_df, "YEAR", scale = "response", plot = FALSE)


annual_indices_df <- data.frame(YEAR = as.numeric(as.character(year_effect_df$fit$YEAR)),
                                Season = "ANNUAL",
                                MODEL = "Bayesian Additive Regression Trees",
                                fit = year_effect_df$fit$visregFit, 
                                lwr = NA,
                                upr = NA) 





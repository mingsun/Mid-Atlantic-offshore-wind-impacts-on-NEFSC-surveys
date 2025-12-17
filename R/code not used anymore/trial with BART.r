

catch_by_tow_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") %>%
  select(ID, REGION)

full_df <- read.csv("results/stratified.mean.indices/surfclam/full.info.tow.list.csv")[,-1] %>%
  filter(OWF == "OUTSIDE") %>%
  left_join(catch_by_tow_df)

remove(catch_by_tow_df)


library(dbarts)


full_BART_df <- full_df %>%   # convert YEAR into a factor
  mutate(YEAR = as.factor(YEAR)) %>%
  filter(REGION == "SVAtoSNE")


BART_temp <- dbarts::bart(x.train = full_BART_df[, c("YEAR", "LAT", "DEPTH", "TEMP")], y.train = full_BART_df[, c("BIOMASS")], keeptrees = TRUE)

invisible(BART_temp$fit$state) # very important if you want to use the models later

summary(BART_temp)

full_BART_df$fit <- fitted(BART_temp)

test <- dbarts:::predict.bart(BART_temp , newdata = full_BART_df[, c("YEAR", "LAT", "DEPTH", "TEMP")], type = "bart")

full_BART_df$pred <- apply(test, 2, median)
full_BART_dflower_bound <- apply(test, 2, quantile, probs = 0.025)
full_BART_dfupper_bound <- apply(test, 2, quantile, probs = 0.975)

library(pdp)

partial_year <- partial(BART_temp, pred.var = "YEAR", train = full_BART_df)
plotPartial(partial_year)


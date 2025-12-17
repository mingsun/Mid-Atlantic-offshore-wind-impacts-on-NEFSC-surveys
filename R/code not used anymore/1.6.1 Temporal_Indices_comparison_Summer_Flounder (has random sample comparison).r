library(tidyverse)


# load indices data ---------------------------------------------------------------------------------------------

FD_indices_df <- read.csv("results/stratified.mean.indices/summer.flounder/original.indices.csv") %>%
  add_column(MODEL = "Design-based full dataset") %>%
  rename(fit = STRATIFIED_MEAN_N, lwr = lo_CI_95, upr = up_CI_95) %>%
  select(YEAR, SEASON, MODEL, fit, lwr, upr)

OD_indices_df <- read.csv("results/stratified.mean.indices/summer.flounder/impacted.indices.csv") %>%
  add_column(MODEL = "Design-based outside dataset") %>%
  rename(fit = STRATIFIED_MEAN_N, lwr = lo_CI_95, upr = up_CI_95) %>%
  select(YEAR, SEASON, MODEL, fit, lwr, upr)

# random_indices_df <- read.csv("results/stratified.mean.indices/summer.flounder/random.reduction.indices.csv")
# random_indices_df$TS <- "RANDOM"

TW.GAM_df <- read.csv("results/temporal model/summer.flounder/Tweedie.GAM_YEAR_Indices.csv")
DELTA.GAM_df <- read.csv("results/temporal model/summer.flounder/Delta.GAM_YEAR_Indices.csv")
ANN_df <- read.csv("results/temporal model/summer.flounder/ANN_YEAR_Indices.csv")
RF_df <- read.csv("results/temporal model/summer.flounder/RF_YEAR_Indices.csv")
VAST_df <- read.csv("results/temporal model/summer.flounder/VAST_YEAR_Indices.csv")

model_based_df <- rbind(TW.GAM_df, DELTA.GAM_df, ANN_df, RF_df, VAST_df) %>%
  # filter(SEASON == "ANNUAL")  %>%# season specific models are not meaningful
  rename(SEASON = Season) %>%
  mutate(MODEL = factor(MODEL, levels = c("Tweedie-GAM", "Delta-GAM", "Random Forest", "Artificial Neural Network", "VAST"))) 

levels(model_based_df$MODEL)

  ## standardize to the model-based indices value to its first year ----

model_based_df <- model_based_df %>%
  group_by(MODEL, SEASON) %>%
  mutate(fit_1982 = fit[YEAR == 1982],
         fit = fit / fit_1982,
         lwr = lwr / fit_1982,
         upr = upr / fit_1982) %>%
  ungroup() %>%
  select(-fit_1982)

# remove(TW.GAM_df, DELTA.GAM_df, ANN_df, RF_df, VAST_df, FD_indices_wide_df)


# 0. indices trend comparison ---------------------------------------------------------------------------------------------

design_based_indices_df <- rbind(FD_indices_df, OD_indices_df)

# random_abundance_df <- random_indices_df %>%
#   group_by(YEAR, SEASON) %>%
#   summarize(MEDIAN_RANDOM = median(STRATIFIED_MEAN_N), 
#             QUANT_0.025 = quantile(STRATIFIED_MEAN_N, probs = c(0.025)),
#             QUANT_0.975 = quantile(STRATIFIED_MEAN_N, probs = c(0.975)))

  ## design-based indices ----

design_based_p <- ggplot(design_based_indices_df, aes(x = YEAR)) +
  # geom_ribbon(data = random_abundance_df, aes(ymin = QUANT_0.025, ymax = QUANT_0.975), fill = "grey", alpha = 0.3) + # first plot the ribbon for the random estimate
  geom_errorbar(aes(ymin = lwr, ymax = upr, color = MODEL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = fit, color = MODEL), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = fit, color = MODEL)) +
  facet_wrap(.~SEASON, nrow = 2) +
  scale_color_manual(values = c("indianred2", "steelblue3"), labels = c("Full dataset", "Outside dataset"), name = "Time Series") +
  ylab("survey indices [N/tow]") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

png("plots/temporal model/summer flounder/design_based_index.png",  width = 8, height = 6, units = 'in', res = 800)
print(design_based_p)
dev.off()

  ## model-based indices ----

model_based_p <- ggplot(model_based_df, aes(x = YEAR)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr, color = MODEL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = fit, color = MODEL), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = fit, color = MODEL)) +
  facet_wrap(MODEL ~ SEASON, nrow = 5) +
  # scale_color_manual(values = c("indianred2", "steelblue3")) +
  ylab("survey indices [standardized]") +
  theme_minimal() +
  theme(legend.position = "none") +
  # theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
  #       legend.justification = c(0, 1)) +
  coord_cartesian(ylim = c(0,6))

png("plots/temporal model/summer flounder/model_based_index.png",  width = 10, height = 10, units = 'in', res = 800)
print(model_based_p)
dev.off()


# 1. correlation (spearman test as it works better with non-linear relationship) ---------------------------------------------------------------------------------------------

cor_df <- expand.grid(REFERENCE = "ORIGINAL", TS = c("IMPACTED", "RANDOM"), SEASON = c("FALL", "SPRING"), COR_MEDIAN = NA, QUANT_0.025 = NA, QUANT_0.975 = NA)
cor_df <- arrange(cor_df, TS)

## between the original and impacted indices
cor_df$COR_MEDIAN[1] <- cor(subset(original_indices_df, SEASON == "FALL")$STRATIFIED_MEAN_N, subset(impacted_indices_df, SEASON == "FALL")$STRATIFIED_MEAN_N, method = "spearman")
cor_df$COR_MEDIAN[2] <- cor(subset(original_indices_df, SEASON == "SPRING")$STRATIFIED_MEAN_N, subset(impacted_indices_df, SEASON == "SPRING")$STRATIFIED_MEAN_N, method = "spearman")

## between the original and random indices
for (s in c("FALL", "SPRING")) {
  
  cor_by_season <- numeric(0)
  
  for (iter in 1:1000) {
    temp_random_df <- subset(random_indices_df, SEASON == s & ITER == iter)
    cor_temp <- cor(subset(original_indices_df, SEASON == s)$STRATIFIED_MEAN_N, temp_random_df$STRATIFIED_MEAN_N, method = "spearman")
    cor_by_season <- c(cor_by_season, cor_temp)
  }
  
  cor_df[cor_df$TS == "RANDOM" & cor_df$SEASON == s, ]$COR_MEDIAN <- median(cor_by_season)
  cor_df[cor_df$TS == "RANDOM" & cor_df$SEASON == s, ]$QUANT_0.025 <- quantile(cor_by_season, probs = c(0.025))
  cor_df[cor_df$TS == "RANDOM" & cor_df$SEASON == s, ]$QUANT_0.975 <- quantile(cor_by_season, probs = c(0.975))

}; remove(s, iter, temp_random_df, cor_temp, cor_by_season)

ggplot(cor_df) +
  geom_errorbar(aes(x = SEASON, ymin = QUANT_0.025, ymax = QUANT_0.975, color = TS), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(x = SEASON, y = COR_MEDIAN, color = TS), position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("indianred2", "steelblue3")) +
  ylab("Spearman's Correlation Coefficient") +
  theme_minimal() +
  theme(legend.position = "none")

# 2. precision loss (relative change in CV by year) ---------------------------------------------------------------------------------------------

## between the original and impacted indices
combine_cv_df <- original_indices_df[,c(1,2,8)] %>%
  rename(ORIGINAL_CV = CV) %>%
  left_join(impacted_indices_df[,c(1,2,8)]) %>%
  rename(IMPACTED_CV = CV) %>%
  mutate(PROPORTIONAL_CHANGE_CV = (IMPACTED_CV - ORIGINAL_CV)/ORIGINAL_CV) %>%
  group_by(SEASON) %>%
  mutate(MEAN_CHANGE_CV = mean(PROPORTIONAL_CHANGE_CV))

## between the original and random indices
random_cv_df <- original_indices_df[,c(1,2,8)] %>%
  rename(ORIGINAL_CV = CV) %>%
  left_join(random_indices_df[,c(1,2,8,9)]) %>%
  rename(RANDOM_CV = CV) %>%
  mutate(PROPORTIONAL_CHANGE_CV = (RANDOM_CV - ORIGINAL_CV)/ORIGINAL_CV) %>%
  group_by(SEASON) %>%
  mutate(MEAN_CHANGE_CV = mean(PROPORTIONAL_CHANGE_CV)) %>%
  group_by(YEAR, SEASON) %>%
  summarize(QUANT_0.025 = quantile(PROPORTIONAL_CHANGE_CV, probs = c(0.025)),
            QUANT_0.975 = quantile(PROPORTIONAL_CHANGE_CV, probs = c(0.975)))

ggplot(combine_cv_df, aes(x = YEAR)) +
  # first plot the ribbon for the random estimate
  geom_ribbon(data = random_cv_df, aes(ymin = QUANT_0.025, ymax = QUANT_0.975), fill = "grey", alpha = 0.3) +
  # then the original and impacted datag
  geom_line(aes(y = MEAN_CHANGE_CV), color = "red", linetype = 2) +
  geom_point(aes(y = PROPORTIONAL_CHANGE_CV)) + 
  geom_line(aes(y = PROPORTIONAL_CHANGE_CV)) +
  facet_wrap(.~SEASON, nrow = 2, scale = "free_y") +
  ylab("relative change in CV") +
  theme_minimal()

# 3. accuracy loss (annual relative error by year) ---------------------------------------------------------------------------------------------

## between the original and impacted indices
combine_re_df <- indices_df[,c(1,2,3,9)] %>%
  group_by(SEASON, TS) %>%
  mutate(GEOMETRIC_MEAN = exp(mean(log(STRATIFIED_MEAN_N))),
         ADJ_STRATIFIED_MEAN_N = STRATIFIED_MEAN_N/GEOMETRIC_MEAN) %>% # standardize the original and impacted indices by diving them by geometric means across years to avoid bias
  ungroup() %>%
  select(-c(STRATIFIED_MEAN_N, GEOMETRIC_MEAN)) %>%
  pivot_wider(names_from = TS, values_from = ADJ_STRATIFIED_MEAN_N) %>%
  mutate(RE = (IMPACTED - ORIGINAL)/ORIGINAL) %>%
  group_by(SEASON) %>%
  mutate(MEDIAN_RE = median(RE))

## between the original and random indices
original_indices_iter_df <- bind_rows(replicate(1000, original_indices_df, simplify = FALSE), .id = "ITER")
original_indices_iter_df$ITER <- as.integer(original_indices_iter_df$ITER )

random_re_df <- random_indices_df %>%
  bind_rows(original_indices_iter_df) %>%
  select(TS, ITER, YEAR, SEASON, STRATIFIED_MEAN_N) %>%
  group_by(TS, ITER, SEASON) %>%
  mutate(GEOMETRIC_MEAN = exp(mean(log(STRATIFIED_MEAN_N))),
         ADJ_STRATIFIED_MEAN_N = STRATIFIED_MEAN_N/GEOMETRIC_MEAN) %>% # standardize the original and random indices by diving them by geometric means across years to avoid bias
  ungroup() %>%
  select(-c(STRATIFIED_MEAN_N, GEOMETRIC_MEAN)) %>%
  pivot_wider(names_from = TS, values_from = ADJ_STRATIFIED_MEAN_N) %>%
  mutate(RE = (RANDOM - ORIGINAL)/ORIGINAL) %>%
  group_by(YEAR, SEASON) %>%
  summarize(MEDIAN_RE_BY_ITER = median(RE),
            QUANT_0.025 = quantile(RE, probs = c(0.025)),
            QUANT_0.975 = quantile(RE, probs = c(0.975)))

ggplot(combine_re_df, aes(x = YEAR)) +
  # first plot the ribbon for the random estimate
  geom_ribbon(data = random_re_df, aes(ymin = QUANT_0.025, ymax = QUANT_0.975), fill = "grey", alpha = 0.3) +
  # then the original and impacted dataggeom_hline(yintercept = 0, color = "grey", linetype = 2) +
  geom_line(aes(y = MEDIAN_RE), color = "red", linetype = 2) +
  geom_point(aes(y = RE)) + 
  geom_line(aes(y = RE)) +
  facet_wrap(.~SEASON, nrow = 2, scale = "free_y") +
  ylab("RE of standardized survey indices") +
  theme_minimal()

# 4. trend bias (trend in RE from last step) ---------------------------------------------------------------------------------------------

## create a full df with two sets of RE
full_re_df <- combine_re_df %>%
  left_join(random_re_df) %>%
  rename(RE_IMPACTED = RE, RE_RANDOM = MEDIAN_RE_BY_ITER) %>%
  select(-c(MEDIAN_RE, QUANT_0.025, QUANT_0.975))

## a df to store the results
slope_df <- expand.grid(REFERENCE = "ORIGINAL", TS = c("IMPACTED", "RANDOM"), SEASON = c("FALL", "SPRING"), SLOPE_MEDIAN = NA, QUANT_0.025 = NA, QUANT_0.975 = NA)
slope_df <- arrange(slope_df, TS)

## between the original and impacted indices
slope_df$SLOPE_MEDIAN[1] <- coef(lm(RE_IMPACTED ~ YEAR, data = subset(full_re_df, SEASON == "FALL")))["YEAR"] * 10 # slope in units of RE per decade
slope_df$SLOPE_MEDIAN[2] <- coef(lm(RE_IMPACTED ~ YEAR, data = subset(full_re_df, SEASON == "SPRING")))["YEAR"] * 10 # slope in units of RE per decade

## between the original and random indices

random_re_by_iter_df <- random_indices_df %>%
  bind_rows(original_indices_iter_df) %>%
  select(TS, ITER, YEAR, SEASON, STRATIFIED_MEAN_N) %>%
  group_by(TS, ITER, SEASON) %>%
  mutate(GEOMETRIC_MEAN = exp(mean(log(STRATIFIED_MEAN_N))),
         ADJ_STRATIFIED_MEAN_N = STRATIFIED_MEAN_N/GEOMETRIC_MEAN) %>% # standardize the original and random indices by diving them by geometric means across years to avoid bias
  ungroup() %>%
  select(-c(STRATIFIED_MEAN_N, GEOMETRIC_MEAN)) %>%
  pivot_wider(names_from = TS, values_from = ADJ_STRATIFIED_MEAN_N) %>%
  mutate(RE = (RANDOM - ORIGINAL)/ORIGINAL) 

 ## loop to calculate slope for each iteration

for (s in c("FALL", "SPRING")) { 
  
  slope_random_by_season <- numeric()
  
  for (iter in 1:1000) { 
    temp_random_df <- subset(random_re_by_iter_df, SEASON == s & ITER == iter)
    temp_result <- coef(lm(RE ~ YEAR, data = temp_random_df))["YEAR"] * 10
    slope_random_by_season <- c(slope_random_by_season, temp_result)
  }

  slope_df[slope_df$TS == "RANDOM" & slope_df$SEASON == s, ]$SLOPE_MEDIAN <- median(slope_random_by_season)
  slope_df[slope_df$TS == "RANDOM" & slope_df$SEASON == s, ]$QUANT_0.025 <- quantile(slope_random_by_season, probs = c(0.025))
  slope_df[slope_df$TS == "RANDOM" & slope_df$SEASON == s, ]$QUANT_0.975 <- quantile(slope_random_by_season, probs = c(0.975))
  
} ; remove(s, iter, temp_random_df, temp_result)

## between the original and impacted indices
ggplot(full_re_df, aes(x = YEAR)) +
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  geom_point(aes(y = RE_IMPACTED)) + 
  geom_smooth(aes(y = RE_IMPACTED), method = "lm", se = TRUE) +
  facet_wrap(.~SEASON, nrow = 2, scale = "free_y") +
  ylab("Trend in RE") +
  theme_minimal() 

## then compare the impacted slope to the and random slope distribution
ggplot(slope_df) +
  geom_errorbar(aes(x = SEASON, ymin = QUANT_0.025, ymax = QUANT_0.975, color = TS), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(x = SEASON, y = SLOPE_MEDIAN, color = TS), position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("indianred2", "steelblue3")) +
  ylab("Trend in RE") +
  theme_minimal() +
  theme(legend.position = "none")


# 5. correlation  between survey exclusion impacts (trend in RE from last step) ---------------------------------------------------------------------------------------------


# (6). correlation  between survey exclusion impacts (trend in RE from last step) ---------------------------------------------------------------------------------------------

# not viable here as we only have 4 species and it makes sense to include more species

# From Sean Anderson et al. ICES MPA paper
# .precision: annual coefficient of of variation CV
# .accruacy: annual relative error (RE)
# .Trend bias: trend in RE with a linear regression by year and report the slope in units of RE per decade. This is like a diverge from the original dataset.
# .Correlates of survey-restriction impacts: Median.absolute.RE ~ CV.status.quo.survey.index + abundance.prop.in/out.closure.area + trends.in.abundane.prop.in/out.closure.area
# 
# For mitigation strategy:
#   In addition to the ones we have
# .Model-based indices
# .down-sampled: randomly remove corrensponding effort loss from the entire area regardless of OWF
# .up-sampled: redistribute the survey effort outside the ofw area


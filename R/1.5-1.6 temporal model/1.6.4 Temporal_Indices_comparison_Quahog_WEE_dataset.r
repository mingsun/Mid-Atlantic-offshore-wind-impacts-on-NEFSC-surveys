library(tidyverse)


# load indices data ---------------------------------------------------------------------------------------------

FD_indices_df <- read.csv("results/stratified.mean.indices/quahog/original.biomass.indices.csv") %>%
  add_column(MODEL = "Design-based full dataset") %>%
  rename(fit = STRATIFIED_MEAN_N, lwr = lo_CI_95, upr = up_CI_95) %>%
  select(YEAR, REGION, MODEL, fit, lwr, upr) %>%
  filter(REGION == "SVAtoSNE")

OD_indices_df <- read.csv("results/stratified.mean.indices/quahog/impacted.biomass.indices.csv") %>%
  add_column(MODEL = "Design-based outside dataset") %>%
  rename(fit = STRATIFIED_MEAN_N, lwr = lo_CI_95, upr = up_CI_95) %>%
  select(YEAR, REGION, MODEL, fit, lwr, upr) %>%
  filter(REGION == "SVAtoSNE")

TW.GAM_df <- read.csv("results/temporal model/quahog/Tweedie.GAM_YEAR_Indices.csv")
DELTA.GAM_df <- read.csv("results/temporal model/quahog/simple.Delta.GAM_YEAR_Indices.csv")
# DELTA.GAM_df <- read.csv("results/temporal model/quahog/Delta.GAM_YEAR_Indices.csv")
BART_df <- read.csv("results/temporal model/quahog/BART_YEAR_Indices.csv")
RF_df <- read.csv("results/temporal model/quahog/RF_YEAR_Indices.csv")
VAST_df <- read.csv("results/temporal model/quahog/VAST_YEAR_Indices.csv")

model_based_df <- rbind(TW.GAM_df, DELTA.GAM_df, BART_df, RF_df, VAST_df) %>%
  filter(REGION == "SVAtoSNE") %>% 
  # rename(REGION = REGION) %>%
  mutate(MODEL = factor(MODEL, 
                        levels = c("Tweedie-GAM", "Delta-GAM", "Random Forest", "Bayesian Additive Regression Trees", "VAST"),
                        labels = c("Tweedie-GAM", "Delta-GAM", "RF", "BART", "VAST")))  

levels(model_based_df$MODEL)

  ## standardize to the model-based indices value to its first year ----

model_based_df <- model_based_df %>%
  group_by(MODEL, REGION) %>%
  mutate(fit_1982 = fit[YEAR == 2002],
         fit = fit / fit_1982,
         lwr = lwr / fit_1982,
         upr = upr / fit_1982) %>%
  ungroup() %>%
  select(-fit_1982)

write.csv(model_based_df, "results/temporal model/quahog/model_based_indices.csv", row.names = FALSE)

# remove(TW.GAM_df, DELTA.GAM_df, ANN_df, RF_df, VAST_df, FD_indices_wide_df)


# 0. indices trend comparison ---------------------------------------------------------------------------------------------

design_based_indices_df <- rbind(FD_indices_df, OD_indices_df)

write.csv(design_based_indices_df, "results/temporal model/quahog/design_based_indices.csv", row.names = FALSE)


  ## design-based indices ----

design_based_p <- ggplot(design_based_indices_df, aes(x = YEAR)) +
  # geom_ribbon(data = random_abundance_df, aes(ymin = QUANT_0.025, ymax = QUANT_0.975), fill = "grey", alpha = 0.3) + # first plot the ribbon for the random estimate
  geom_errorbar(aes(ymin = lwr, ymax = upr, color = MODEL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = fit, color = MODEL), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = fit, color = MODEL)) +
  facet_wrap(.~REGION, nrow = 2) +
  scale_color_manual(values = c("indianred2", "steelblue3"), labels = c("Full dataset", "Outside dataset"), name = "Time Series") +
  ylab("biomass indices [g/tow]") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

png("plots/temporal model/quahog/design_based_index.png",  width = 8, height = 6, units = 'in', res = 800)
print(design_based_p)
dev.off()

  ## model-based indices ----

model_based_p <- ggplot(model_based_df, aes(x = YEAR)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr, color = MODEL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = fit, color = MODEL), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = fit, color = MODEL)) +
  facet_wrap(MODEL ~ REGION, nrow = 5) +
  # scale_color_manual(values = c("indianred2", "steelblue3")) +
  ylab("biomass indices [standardized]") +
  theme_minimal() +
  theme(legend.position = "none") +
  # theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
  #       legend.justification = c(0, 1)) +
  coord_cartesian(ylim = c(0,6))

png("plots/temporal model/quahog/model_based_index.png",  width = 4, height = 8, units = 'in', res = 800)
print(model_based_p)
dev.off()


# 1. correlation (spearman test as it works better with non-linear relationship) ---------------------------------------------------------------------------------------------

non_FD_df <- rbind(OD_indices_df, model_based_df) 

cor_df <- expand.grid(MODEL = c("Design-based outside dataset", as.character(levels(model_based_df$MODEL))), REGION = c("SVAtoSNE"), Spearman_Cor = NA)

for (i in 1:nrow(cor_df)) {
  
  temp_ind <- subset(non_FD_df, MODEL == cor_df$MODEL[i]) 
  temp_ref <- OD_indices_df
  
  temp_year <- intersect(temp_ind$YEAR, temp_ref$YEAR)
  
  temp_ind <- subset(temp_ind, YEAR %in% temp_year) 
  temp_ref <- subset(temp_ref, YEAR %in% temp_year)
  
  cor_df$Spearman_Cor[i] <- cor(temp_ind$fit, temp_ref$fit, method = "spearman")
  
}; remove(temp_ind, temp_ref, temp_year)


cor_df <- cor_df %>%
  group_by(MODEL) %>%
  mutate(MEAN_COR = mean(Spearman_Cor)) %>%
  ungroup() %>%
  mutate(SCENARIO = factor(REGION, levels = c("SVAtoSNE")))

write.csv(cor_df, "results/temporal model/quahog/correlation.df.csv", row.names = FALSE)
  

library(RColorBrewer)
display.brewer.all()
palette <- brewer.pal(4, "Set1")

cor_p <- ggplot(cor_df) +
  geom_bar(aes(x = SCENARIO, y = Spearman_Cor, fill = SCENARIO), position = "dodge", stat = "identity", linewidth = 1) +
  # geom_hline(aes(yintercept = MEAN_COR), linetype = 2, linewidth = 1) +
  facet_wrap(.~MODEL, nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  ylab("Spearman Correlation Coefficient") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  scale_fill_manual(values = palette, name = "")

png("plots/temporal model/quahog/correlation.png",  width = 12, height = 6, units = 'in', res = 800)
print(cor_p)
dev.off()

rm(list = setdiff(ls(), c("non_FD_df", "FD_indices_df", "OD_indices_df", "MODEL_based_df")))

# -------------------------------------------------------------------------------------------------------- #



# 2. precision loss (relative change in CV by year) ---------------------------------------------------------------------------------------------

 ## prepare reference data
OD_CV_df <- read.csv("results/stratified.mean.indices/quahog/impacted.biomass.indices.csv") %>%
  filter(REGION == "SVAtoSNE") %>%
  select(YEAR, REGION, CV) %>%
  rename(REF_CV = CV, REGION_REF = REGION)

 
 ## assemble the full df
CV_df <- non_FD_df %>%
  filter(upr != 0) %>%
  mutate(CV = (((upr - fit) / qnorm(0.975)) / fit)) %>% # calculate CV based on CI
  left_join(OD_CV_df) %>%
  filter(!is.na(REF_CV), !is.na(CV)) %>%
  mutate(PROPORTIONAL_CHANGE_CV = (CV - REF_CV)/REF_CV) %>%
  group_by(MODEL) %>%
  mutate(MEAN_CHANGE_CV = mean(PROPORTIONAL_CHANGE_CV)) %>%
  ungroup() %>%
  mutate(MODEL = factor(MODEL, levels = c("Design-based outside dataset", "Tweedie-GAM", "Delta-GAM", "RF", "BART", "VAST"))) %>%
  arrange(MODEL)

write.csv(CV_df, "results/temporal model/quahog/CV.df.csv", row.names = FALSE)


pre_loss_p <- ggplot(CV_df, aes(x = YEAR)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(y = MEAN_CHANGE_CV), color = "indianred") +
  geom_point(aes(y = PROPORTIONAL_CHANGE_CV)) + 
  geom_line(aes(y = PROPORTIONAL_CHANGE_CV)) +
  facet_wrap(. ~ MODEL, scale = "free_y", nrow = 5) +
  geom_text(aes(x = min(YEAR), y = Inf, label = paste0("Mean=", round(MEAN_CHANGE_CV,2))), 
            hjust = -0.1, vjust = 1.1, color = "indianred") +
  ylab("relative change in CV") +
  theme_classic()

png("plots/temporal model/quahog/precision_loss_year.png",  width = 8, height = 8, units = 'in', res = 800)
print(pre_loss_p)
dev.off()

rm(list = setdiff(ls(), c("non_FD_df", "FD_indices_df", "OD_indices_df", "MODEL_based_df")))

# -------------------------------------------------------------------------------------------------------- #



# 3. accuracy loss (annual relative error by year) ---------------------------------------------------------------------------------------------

## prepare reference data
OD_RE_df <- OD_indices_df %>%
  mutate(GEOMETRIC_MEAN = exp(mean(log(fit))),
         fit_REF = fit/GEOMETRIC_MEAN) %>% # standardize the values by diving them by geometric means across years to avoid bias
  select(YEAR, REGION, fit_REF) 

# ## prepare annual data
# RE_annual_df <- filter(non_FD_df, REGION == "ANNUAL") %>%
#   bind_rows(filter(non_FD_df, REGION == "ANNUAL")) %>%
#   rename(REGION_IND = REGION) %>%
#   add_column(REGION_REF = c(rep("FALL", nrow(.)/2), c(rep("SPRING", nrow(.)/2)))) %>%
#   select(YEAR, MODEL, REGION_IND, REGION_REF, fit)

## assemble the full df
RE_df <- non_FD_df %>%
  filter(fit != 0) %>%
  group_by(MODEL) %>%
  mutate(GEOMETRIC_MEAN = exp(mean(log(fit))),
       fit = fit/GEOMETRIC_MEAN) %>%  # standardize the values by diving them by geometric means across years to avoid bias
  ungroup() %>%
  select(YEAR, REGION, MODEL, fit) %>%
  left_join(OD_RE_df) %>%
  filter(!is.na(fit) & !is.na(fit_REF)) %>%
  mutate(RE = (fit - fit_REF)/fit_REF) %>%
  group_by(MODEL) %>%
  mutate(MEAN_RE = mean(RE)) %>%
  ungroup()  %>%
  mutate(MODEL = factor(MODEL, levels = c("Design-based outside dataset", "Tweedie-GAM", "Delta-GAM", "RF", "BART", "VAST"))) %>%
  arrange(MODEL) 

write.csv(RE_df, "results/temporal model/quahog/RE.df.csv", row.names = FALSE)


acc_loss_p <- ggplot(RE_df, aes(x = YEAR)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(y = MEAN_RE), color = "indianred") +
  geom_point(aes(y = RE)) + 
  geom_line(aes(y = RE)) +
  facet_wrap(. ~ MODEL, scale = "free_y", nrow = 5) +
  geom_text(aes(x = min(YEAR), y = Inf, label = paste0("Mean=", round(MEAN_RE,2))), 
            hjust = -0.1, vjust = 1.1, color = "indianred") +
  ylab("relative change in RE") +
  theme_classic()


png("plots/temporal model/quahog/accuracy_loss_year.png",  width = 12, height = 8, units = 'in', res = 800)
print(acc_loss_p)
dev.off()

# rm(list = setdiff(ls(), c("non_FD_df", "FD_indices_df", "OD_indices_df", "MODEL_based_df")))

# -------------------------------------------------------------------------------------------------------- #




# 4. trend bias (trend in correlation between RE and YEAR from last step) ---------------------------------------------------------------------------------------------

trend_bias_df <- expand.grid(MODEL = c("Design-based outside dataset", as.character(levels(model_based_df$MODEL))), REGION = c("SVAtoSNE"), DECADAL_SLOPE = NA, R_SQUARE = NA)

# i=1
for (i in 1:nrow(trend_bias_df)) {
  
  temp_df <- subset(RE_df, MODEL == trend_bias_df$MODEL[i])
  temp_df <- subset(temp_df, RE <= 5 & RE >= -0.9)
  temp_lm <- try(lm(RE ~ YEAR, data = temp_df), silent = TRUE)
  
  if(nrow(temp_df) == 1) next

  if (inherits(temp_lm, "try-error")) {
    next  # skip this iteration
  }
  
  
  trend_bias_df$DECADAL_SLOPE[i] <- round(coef(temp_lm)["YEAR"] * 10, 3) # slope in units of RE per decade
  trend_bias_df$R_SQUARE[i] <- round(summary(temp_lm)$r.squared, 2)
  trend_bias_df$P_VALUE[i] <- round(summary(temp_lm)$coefficients[2,4], 3)
  
}; remove(temp_df, temp_lm)


 ## assemble
trend_bias_df <- trend_bias_df %>%
  right_join(RE_df)

write.csv(trend_bias_df, "results/temporal model/quahog//trend_bias_df.csv", row.names = FALSE)


trend_re_p <- ggplot(trend_bias_df, aes(x = YEAR)) +
  geom_point(aes(y = RE)) + 
  geom_smooth(aes(y = RE), method = "lm", se = TRUE) +
  facet_wrap(. ~ MODEL, scale = "free_y", nrow = 4) +
  geom_text(aes(x = min(YEAR), y = Inf, label = paste0("Decadal.Slope=", DECADAL_SLOPE)), 
            hjust = -0.1, vjust = 1.1, color = "indianred") +
  ylab("Trend in RE") +
  theme_classic()

png("plots/temporal model/quahog/trend_in_RE.png",  width = 12, height = 8, units = 'in', res = 800)
print(trend_re_p)
dev.off()





library(tidyverse)


# load full dataset abundance indices ---------------------------------------------------------------------------------------------

TW.GAM_df <- read.csv("results/temporal model/surfclam/full dataset/Tweedie.GAM_YEAR_Indices.csv")
DELTA.GAM_df <- read.csv("results/temporal model/surfclam/full dataset/simple.Delta.GAM_YEAR_Indices.csv")
# DELTA.GAM_df <- read.csv("results/temporal model/surfclam/Delta.GAM_YEAR_Indices.csv")
BART_df <- read.csv("results/temporal model/surfclam/full dataset/BART_YEAR_Indices.csv")
RF_df <- read.csv("results/temporal model/surfclam/full dataset/RF_YEAR_Indices.csv")
VAST_df <- read.csv("results/temporal model/surfclam/full dataset/VAST_YEAR_Indices.csv")

full_ind_df <- rbind(TW.GAM_df, DELTA.GAM_df, BART_df, RF_df, VAST_df) %>%
  filter(REGION == "SVAtoSNE") %>%
  add_column(Dataset = "FULL", .after = 0) %>%
  mutate(MODEL = factor(MODEL, 
                        levels = c("Tweedie-GAM", "Delta-GAM", "Random Forest", "Bayesian Additive Regression Trees", "VAST"),
                        labels = c("Tweedie-GAM", "Delta-GAM", "RF", "BART", "VAST"))) 

remove(TW.GAM_df, DELTA.GAM_df, BART_df, RF_df, VAST_df)




# load WEE dataset abundance indices ---------------------------------------------------------------------------------------------

TW.GAM_df <- read.csv("results/temporal model/surfclam/Tweedie.GAM_YEAR_Indices.csv")
DELTA.GAM_df <- read.csv("results/temporal model/surfclam/simple.Delta.GAM_YEAR_Indices.csv")
# DELTA.GAM_df <- read.csv("results/temporal model/surfclam/Delta.GAM_YEAR_Indices.csv")
BART_df <- read.csv("results/temporal model/surfclam/BART_YEAR_Indices.csv")
RF_df <- read.csv("results/temporal model/surfclam/RF_YEAR_Indices.csv")
VAST_df <- read.csv("results/temporal model/surfclam/VAST_YEAR_Indices.csv")

wee_ind_df <- rbind(TW.GAM_df, DELTA.GAM_df, BART_df, RF_df, VAST_df) %>%
  filter(REGION == "SVAtoSNE") %>%
  add_column(Dataset = "WEE", .after = 0) %>%
  mutate(MODEL = factor(MODEL, 
                        levels = c("Tweedie-GAM", "Delta-GAM", "Random Forest", "Bayesian Additive Regression Trees", "VAST"),
                        labels = c("Tweedie-GAM", "Delta-GAM", "RF", "BART", "VAST"))) 

remove(TW.GAM_df, DELTA.GAM_df, BART_df, RF_df, VAST_df)



# combine and edit ---------------------------------------------------------------------------------------------

ind_df <- rbind(full_ind_df, wee_ind_df); remove(full_ind_df, wee_ind_df)

ind_df <- ind_df %>%       # standardize the model-based indices value to its median value 
  group_by(Dataset, MODEL, REGION) %>%
  mutate(fit_1982 = fit[YEAR == 2002],
         fit = fit / fit_1982,
         lwr = lwr / fit_1982,
         upr = upr / fit_1982) %>%
  ungroup() %>%
  select(-fit_1982) 






# 0. indices trend visual comparison ---------------------------------------------------------------------------------------------

indices_p <- ggplot(ind_df, aes(x = YEAR)) +
  # geom_ribbon(data = random_abundance_df, aes(ymin = QUANT_0.025, ymax = QUANT_0.975), fill = "grey", alpha = 0.3) + # first plot the ribbon for the random estimate
  geom_errorbar(aes(ymin = lwr, ymax = upr, color = Dataset), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = fit, color = Dataset), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = fit, color = Dataset)) +
  facet_grid(MODEL ~ REGION) +
  scale_color_manual(values = c("indianred2", "steelblue3"), labels = c("Full dataset", "WEE dataset"), name = "Time Series") +
  ylab("biomass indices [g/tow]") +
  theme_minimal() +
  theme(legend.position = c(0.02, 0.98), # Place the legend in the top-left corner
        legend.justification = c(0, 1)) +
  coord_cartesian(y = c(0,6))

png("plots/temporal model/surfclam/model_based_indices_comparison.png",  width = 8, height = 6, units = 'in', res = 800)
print(indices_p)
dev.off()

  




# 1. correlation (spearman test as it works better with non-linear relationship) ---------------------------------------------------------------------------------------------

cor_df <- expand.grid(MODEL = as.character(levels(ind_df$MODEL)), REGION = c("SVAtoSNE"), Spearman_Cor = NA)




# i = 1
for (i in 1:nrow(cor_df)) {
  
  temp_ind <- subset(ind_df, Dataset == "WEE"  & MODEL == cor_df$MODEL[i])
  temp_ref <- subset(ind_df, Dataset == "FULL" & MODEL == cor_df$MODEL[i])
  
  temp_year <- intersect(temp_ind$YEAR, temp_ref$YEAR)
  
  temp_ind <- subset(temp_ind, YEAR %in% temp_year) 
  temp_ref <- subset(temp_ref, YEAR %in% temp_year)
  
  cor_df$Spearman_Cor[i] <- cor(temp_ind$fit, temp_ref$fit, method = "spearman")
  
}; remove(temp_ind, temp_ref, temp_year)


write.csv(cor_df, "results/temporal model/surfclam/full dataset/correlation.between.model-based.indices.df.csv", row.names = FALSE)
  

# library(RColorBrewer)
# display.brewer.all()
# palette <- brewer.pal(4, "Set1")
# 
# cor_p <- ggplot(cor_df) +
#   geom_bar(aes(x = SCENARIO, y = Spearman_Cor, fill = SCENARIO), position = "dodge", stat = "identity", linewidth = 1) +
#   # geom_hline(aes(yintercept = MEAN_COR), linetype = 2, linewidth = 1) +
#   facet_wrap(.~MODEL, nrow = 1) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1),
#         legend.position = "none") +
#   ylab("Spearman Correlation Coefficient") +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
#   scale_fill_manual(values = palette, name = "")
# 
# png("plots/temporal model/surfclam/correlation.png",  width = 12, height = 6, units = 'in', res = 800)
# print(cor_p)
# dev.off()
# 
# rm(list = setdiff(ls(), c("non_FD_df", "FD_indices_df", "OD_indices_df", "MODEL_based_df")))

# -------------------------------------------------------------------------------------------------------- #





# 2. precision loss (relative change in CV by year) ---------------------------------------------------------------------------------------------

CV_df <- ind_df %>%
  mutate(CV = (((upr - fit) / qnorm(0.975)) / fit)) %>% # calculate CV based on CI
  select(Dataset, YEAR, MODEL, CV) %>%
  pivot_wider(names_from = Dataset, values_from = CV) %>% # convert from long to wide
  filter(!is.na(FULL) & !is.na(WEE)) %>% # remove years with NA in cv (RF) or with NaN (2020 fall estimate, where mo survey occurred)
  mutate(PROPORTIONAL_CHANGE_CV = (WEE - FULL)/FULL) %>%
  group_by(MODEL) %>%
  mutate(MEAN_CHANGE_CV = mean(PROPORTIONAL_CHANGE_CV))

write.csv(CV_df, "results/temporal model/surfclam/full dataset/CV.between.model-based.indices.df.csv", row.names = FALSE)


# pre_loss_p <- ggplot(CV_df, aes(x = YEAR)) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_line(aes(y = MEAN_CHANGE_CV), color = "indianred") +
#   geom_point(aes(y = PROPORTIONAL_CHANGE_CV)) + 
#   geom_line(aes(y = PROPORTIONAL_CHANGE_CV)) +
#   facet_wrap(. ~ MODEL, scale = "free_y", nrow = 5) +
#   geom_text(aes(x = min(YEAR), y = Inf, label = paste0("Mean=", round(MEAN_CHANGE_CV,2))), 
#             hjust = -0.1, vjust = 1.1, color = "indianred") +
#   ylab("relative change in CV") +
#   theme_classic()
# 
# png("plots/temporal model/surfclam/precision_loss_year.png",  width = 8, height = 8, units = 'in', res = 800)
# print(pre_loss_p)
# dev.off()
 

# -------------------------------------------------------------------------------------------------------- #




# 3. accuracy loss (annual relative error by year) ---------------------------------------------------------------------------------------------

RE_df <- ind_df %>%
  filter(fit != 0) %>% # remove years with 0 (2014 Delta GAM and VAST)
  mutate(GEOMETRIC_MEAN = exp(mean(log(fit))),
         fit_REF = fit/GEOMETRIC_MEAN) %>% # standardize the values by diving them by geometric means across years to avoid bias
  select(Dataset, YEAR, MODEL, fit) %>%
  pivot_wider(names_from = Dataset, values_from = fit) %>% # convert from long to wide
  mutate(RE = (WEE - FULL)/FULL) %>%
  group_by(MODEL) %>%
  mutate(MEAN_RE = mean(RE))

write.csv(RE_df, "results/temporal model/surfclam/full dataset/RE.between.model-based.indices.df.csv", row.names = FALSE)
  
# acc_loss_p <- ggplot(RE_df, aes(x = YEAR)) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_line(aes(y = MEAN_RE), color = "indianred") +
#   geom_point(aes(y = RE)) + 
#   geom_line(aes(y = RE)) +
#   facet_wrap(. ~ MODEL, scale = "free_y", nrow = 5) +
#   geom_text(aes(x = min(YEAR), y = Inf, label = paste0("Mean=", round(MEAN_RE,2))), 
#             hjust = -0.1, vjust = 1.1, color = "indianred") +
#   ylab("relative change in RE") +
#   theme_classic()
# 
# 
# png("plots/temporal model/surfclam/accuracy_loss_year.png",  width = 12, height = 8, units = 'in', res = 800)
# print(acc_loss_p)
# dev.off()

# rm(list = setdiff(ls(), c("non_FD_df", "FD_indices_df", "OD_indices_df", "MODEL_based_df")))

# -------------------------------------------------------------------------------------------------------- #




# 4. trend bias (trend in correlation between RE and YEAR from last step) ---------------------------------------------------------------------------------------------

trend_bias_df <- expand.grid(MODEL = as.character(levels(ind_df$MODEL)), REGION = c("SVAtoSNE"), DECADAL_SLOPE = NA, R_SQUARE = NA)

# i=1
for (i in 1:nrow(trend_bias_df)) {
  
  temp_df <- subset(RE_df, MODEL == trend_bias_df$MODEL[i])
  temp_lm <- lm(RE ~ YEAR, data = temp_df)
  
  trend_bias_df$DECADAL_SLOPE[i] <- round(coef(temp_lm)["YEAR"] * 10, 3) # slope in units of RE per decade
  trend_bias_df$R_SQUARE[i] <- round(summary(temp_lm)$r.squared, 2)
  trend_bias_df$P_VALUE[i] <- round(summary(temp_lm)$coefficients[2,4], 3)
  
}; remove(temp_df, temp_lm)






write.csv(trend_bias_df, "results/temporal model/surfclam/full dataset/trend.bias.between.model-based.indices.df.csv", row.names = FALSE)


# trend_re_p <- ggplot(trend_bias_df, aes(x = YEAR)) +
#   geom_point(aes(y = RE)) + 
#   geom_smooth(aes(y = RE), method = "lm", se = TRUE) +
#   facet_wrap(. ~ MODEL, scale = "free_y", nrow = 4) +
#   geom_text(aes(x = min(YEAR), y = Inf, label = paste0("Decadal.Slope=", DECADAL_SLOPE)), 
#             hjust = -0.1, vjust = 1.1, color = "indianred") +
#   ylab("Trend in RE") +
#   theme_classic()
# 
# png("plots/temporal model/surfclam/trend_in_RE.png",  width = 12, height = 8, units = 'in', res = 800)
# print(trend_re_p)
# dev.off()





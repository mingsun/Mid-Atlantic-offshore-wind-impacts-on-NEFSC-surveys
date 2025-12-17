library(tidyverse)
library(mgcv)
library(VAST)
library(caret)
library(ggplot2)
library(openair)
library(randomForest)
library(dbarts)
detach("package:plyr", unload = TRUE)


# 1. RMSE: accuracy between predictions of “full dataset” and the “outside-OWA dataset” ----

TW.GAM_comp <- read.csv("results/spatial model/surfclam/Tweedie.GAM_abundance_comparison.csv")
DELTA.GAM_comp <- read.csv("results/spatial model/surfclam/Delta.GAM_abundance_comparison.csv")
RF_comp <- read.csv("results/spatial model/surfclam/RF_abundance_comparison.csv")
BART_comp <- read.csv("results/spatial model/surfclam/BART_abundance_comparison.csv")
VAST_comp <- read.csv("results/spatial model/surfclam/VAST_abundance_comparison.csv"); colnames(VAST_comp)[2:3] <- c("YEAR", "LAT")

Pred_df <- bind_rows(TW.GAM_comp, DELTA.GAM_comp[,-c(7,11)], RF_comp, BART_comp[,-c(6,11)], VAST_comp[,-c(4,7)]) %>%
  rename(AREA = OWF) %>%
  mutate(MODEL = factor(MODEL, levels = c("Tweedie-GAM", "Delta-GAM", "Random Forest", "Bayesian Additive Regression Trees", "VAST"),
                        labels = c("Tweedie-GAM", "Delta-GAM", "RF", "BART", "VAST"))) 

write.csv(Pred_df, "results/spatial model/surfclam/Pred_df.csv")

  ## calculate RMSE by OWF area
Pred_RMSE <- Pred_df %>%
  group_by(MODEL, DATASET, AREA) %>%
  summarize(RMSE = caret::RMSE(PREDICTED_BIOMASS, BIOMASS))

  ## calculate RMSE for area-aggregated results
Pred_RMSE_all <- Pred_df %>%
  group_by(MODEL, DATASET) %>%
  summarize(RMSE = caret::RMSE(PREDICTED_BIOMASS, BIOMASS)) %>%
  add_column(AREA = "ALL") %>%
  bind_rows(Pred_RMSE) %>%
  arrange(MODEL, DATASET, AREA) %>%
  select(MODEL, DATASET, AREA, RMSE)

remove(TW.GAM_comp, DELTA.GAM_comp, RF_comp, BART_comp, VAST_comp, Pred_RMSE)

  ## convert the results to show the difference in RMSE between using full and outside dataset for better comparison
RMSE_diff_df <- Pred_RMSE_all %>%
  group_by(MODEL, AREA) %>%
  summarize(Diff_RMSE = RMSE[DATASET == "OUTSIDE"] - RMSE[DATASET == "FULL"])

write.csv(RMSE_diff_df, "results/spatial model/surfclam/RMSE_diff_df.csv")

  ## plot
facet_labels <- c("ALL" = "Entire Region", "INSIDE" = "Inside WEA", "OUTSIDE" = "Outside WEA")

RMSE_diff_p <- ggplot(RMSE_diff_df) +
  geom_bar(aes(x = MODEL, y = Diff_RMSE, fill = factor(sign(Diff_RMSE))), position = "dodge", stat = "identity", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(.~AREA, nrow = 1, labeller = labeller(AREA = facet_labels)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  ylab("RMSE.Outside.Dataset - RMSE.Full.Dataset") +
  # scale_y_continuous(breaks = seq(-2, 2, by = 0.1)) +
  scale_fill_manual(values = c("-1" = "springgreen3", "1" = "indianred2"),
                    labels = c("Decrease", "Increase"),
                    name = "") +
  coord_cartesian(ylim = c(-0.2, 1.2))

png("plots/spatial model/surfclam/RMSE_difference_between_FD_OD.png",  width = 8, height = 8, units = 'in', res = 800)
print(RMSE_diff_p)
dev.off()

# -------------------------------------------------------------------------------------------------------- #




# 2. RMSE vs effort/sample loss by year: build relationship ----

  ## effort and sample loss by year
effort_df <- read.csv("results/AS_OQ_survey_effort_loss.csv")[,2:5] %>%
  rename(YEAR = EST_YEAR, EFFORT_LOSS = ratio) %>%
  filter(YEAR >= 1982) %>% # summer flounder assessment uses survey data since 1982
  group_by(YEAR) %>%
  summarize(EFFORT_LOSS = sum(overlay.n)/sum(total.n))

sample_df <- read.csv("results/AS_OQ_sample_loss_all_species.csv") %>%
  filter(SVSPP == 403,      # (summer flounder: 103, surf clam: 403, quahog 409, longfin squid 503)
         YEAR >= 1982) %>%  # summer flounder assessment uses survey data since 1982
  select(c(3,6)) %>%
  rename(SAMPLE_LOSS = precluded.ratio)  
  
loss_info_df <- merge(effort_df, sample_df); remove(effort_df, sample_df)

  ## calculate RMSE for inside area by year
Pred_RMSE_by_year <- Pred_df %>%
  filter(AREA == "INSIDE") %>%
  group_by(MODEL, DATASET, YEAR) %>%
  summarize(RMSE = caret::RMSE(PREDICTED_BIOMASS, BIOMASS)) %>%
  group_by(MODEL, YEAR) %>%
  summarize(Diff_RMSE = RMSE[DATASET == "OUTSIDE"] - RMSE[DATASET == "FULL"]) %>%
  left_join(loss_info_df) %>%
  mutate(Diff_RMSE_std = scale(Diff_RMSE), # standardize the variables
         EFFORT_LOSS_std = scale(EFFORT_LOSS),
         SAMPLE_LOSS_std = scale(SAMPLE_LOSS)); remove(loss_info_df)

  ## plot the linear relationship
# Pred_RMSE_by_year <- subset(Pred_RMSE_by_year, SAMPLE_LOSS != 0)

library(ggpmisc) # the package to show slope and R-square

RMSE_effort_p <- ggplot(Pred_RMSE_by_year, aes(x = EFFORT_LOSS_std, y = Diff_RMSE_std)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")  +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = "top") +
  facet_wrap(.~MODEL, ncol = 2) +
  labs(title = "Relationship between effort loss and Diff_RMSE in standardized scale") +
  theme_classic()

png("plots/spatial model/summer.flounder/RMSE_effort_relationship.png",  width = 6, height = 9, units = 'in', res = 800)
print(RMSE_effort_p)
dev.off()

RMSE_sample_p <- ggplot(Pred_RMSE_by_year, aes(x = SAMPLE_LOSS_std, y = Diff_RMSE_std)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")  +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = "top") +
  facet_wrap(.~MODEL, ncol = 2) +
  labs(title = "Relationship between sample loss and Diff_RMSE in standardized scale") +
  theme_classic()

png("plots/spatial model/summer.flounder/RMSE_sample_relationship.png",  width = 6, height = 9, units = 'in', res = 800)
print(RMSE_sample_p)
dev.off()

detach("package:ggpmisc", unload = TRUE)

# -------------------------------------------------------------------------------------------------------- #




# 3. Nash-Sutcliffe Efficiency (NSE) of INSIDE data: accuracy of model prediction relative to the average of the observation ----

  ## Function to calculate the reliability index (-Inf to 1, with 1 corresponding to a perfect match)
nse <- function(observation, prediction) {
  if (length(observation) != length(prediction)) {
    stop("The lengths of observation and prediction values must be the same.")
  }
  numerator <- sum((observation - prediction)^2)
  denominator <- sum((observation - mean(observation))^2)
  nse_value <- 1 - (numerator / denominator)
  return(nse_value)
}

  ## Calculate the NSE
NSE_df <- Pred_df %>%
  # filter(AREA == "INSIDE") %>%
  group_by(MODEL, DATASET, AREA) %>%
  summarise(NSE = nse(BIOMASS, PREDICTED_BIOMASS))

## calculate NSE for area-aggregated results
NSE_all_df <- Pred_df %>%
  group_by(MODEL, DATASET) %>%
  summarize(NSE = nse(BIOMASS, PREDICTED_BIOMASS)) %>%
  add_column(AREA = "ALL") %>%
  bind_rows(NSE_df) %>%
  arrange(MODEL, DATASET, AREA) %>%
  select(MODEL, DATASET, AREA, NSE)

## convert the results to show the difference in NSE between using full and outside dataset for better comparison
NSE_diff_df <- NSE_all_df %>%
  group_by(MODEL, AREA) %>%
  summarize(Diff_NSE = NSE[DATASET == "OUTSIDE"] - NSE[DATASET == "FULL"])

remove(NSE_df)

  ## plot
NSE_diff_plot <- ggplot(NSE_diff_df) +
  geom_bar(aes(x = MODEL, y = Diff_NSE, fill = factor(sign(Diff_NSE))), position = "dodge", stat = "identity", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(.~AREA, nrow = 1, labeller = labeller(AREA = facet_labels)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  ylab("NSE.Outside.Dataset - NSE.Full.Dataset") +
  scale_y_continuous(breaks = seq(-2, 2, by = 0.1)) +
  scale_fill_manual(values = c("-1" = "indianred2", "1" = "springgreen3"),
                    labels = c("Decrease", "Increase"),
                    name = "") +
  coord_cartesian(ylim = c(-0.1, 0.02))

png("plots/spatial model/summer.flounder/NSE_difference_between_FD_OD.png",  width = 8, height = 8, units = 'in', res = 800)
print(NSE_diff_plot)
dev.off()



# -------------------------------------------------------------------------------------------------------- #



# 4. Coefficient of Variation (CV) of residual: Prediction precision ----

  ## functions to calculate CV
cv_residuals <- function(observation, prediction) {
  residuals <- observation - prediction
  return(sd(residuals) / mean(abs(residuals)))
}

## Calculate the CV
CV_df <- Pred_df %>%
  # filter(AREA == "INSIDE") %>%
  group_by(MODEL, DATASET, AREA) %>%
  summarise(CV = cv_residuals(BIOMASS, PREDICTED_BIOMASS))

## calculate CV for area-aggregated results
CV_all_df <- Pred_df %>%
  group_by(MODEL, DATASET) %>%
  summarize(CV = cv_residuals(BIOMASS, PREDICTED_BIOMASS)) %>%
  add_column(AREA = "ALL") %>%
  bind_rows(CV_df) %>%
  arrange(MODEL, DATASET, AREA) %>%
  select(MODEL, DATASET, AREA, CV)

## convert the results to show the difference in NSE between using full and outside dataset for better comparison
CV_diff_df <- CV_all_df %>%
  group_by(MODEL, AREA) %>%
  summarize(Diff_CV = CV[DATASET == "OUTSIDE"] - CV[DATASET == "FULL"])

remove(CV_df)

write.csv(CV_diff_df, "results/spatial model/surfclam/CV_diff_df.csv")

## plot
CV_diff_plot <- ggplot(CV_diff_df) +
  geom_bar(aes(x = MODEL, y = Diff_CV, fill = factor(sign(Diff_CV))), position = "dodge", stat = "identity", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(.~AREA, nrow = 1, labeller = labeller(AREA = facet_labels)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  ylab("CV.Outside.Dataset - CV.Full.Dataset") +
  scale_y_continuous(breaks = seq(-2, 2, by = 0.1)) +
  scale_fill_manual(values = c("1" = "indianred2", "-1" = "springgreen3"),
                    # labels = c("Decrease", "Increase"),
                    name = "")

png("plots/spatial model/surfclam/CV_difference_between_FD_OD.png",  width = 8, height = 8, units = 'in', res = 800)
print(CV_diff_plot)
dev.off()

# -------------------------------------------------------------------------------------------------------- #




# 5. Taylor Diagram for model evaluation with conditioning ----

TD_df <- Pred_df %>%
  mutate(MODEL.DATASET = paste(MODEL, DATASET, sep = " * "))

TD_p <- TaylorDiagram(TD_df, obs = "BIOMASS", mod = "PREDICTED_BIOMASS", normalise = TRUE, group = "MODEL.DATASET")

png("plots/spatial model/surfclam/Taylor_Diagram_both.png",  width = 8, height = 8, units = 'in', res = 800)
print(TD_p)
dev.off()

TD_FD_p <- TaylorDiagram(subset(Pred_df, DATASET == "FULL"), obs = "BIOMASS", mod = "PREDICTED_BIOMASS", normalise = TRUE, group = "MODEL")

png("plots/spatial model/surfclam/Taylor_Diagram_FD.png",  width = 5, height = 5, units = 'in', res = 800)
print(TD_FD_p)
dev.off()


TD_OD_p <- TaylorDiagram(subset(Pred_df, DATASET == "OUTSIDE"), obs = "BIOMASS", mod = "PREDICTED_BIOMASS", normalise = TRUE, group = "MODEL")

png("plots/spatial model/surfclam/Taylor_Diagram_OD.png",  width = 5, height = 5, units = 'in', res = 800)
print(TD_OD_p)
dev.off()



# -------------------------------------------------------------------------------------------------------- #


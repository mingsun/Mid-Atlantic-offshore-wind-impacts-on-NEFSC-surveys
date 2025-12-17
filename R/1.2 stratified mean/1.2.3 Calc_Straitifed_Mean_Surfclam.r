library(tidyverse)

# load tow data -----------------------------------------------------------------------

# this is a clean version provided by Dan. Data are already filtered so good to use.
# but separate the data compiling because sensor were introduced in 1997 to measure the accurate tow swept area


  ## 1982 - 1994 data ---- 

tow_1982_df <- read.csv("data/NOAA.stock.data/surfclam/TOW_DATA.CSV") %>% 
  select(CRUISE6, YR, STATION, STRATUM, TOW, LENGRP, NPERTOW, KGPERTOW, ASWEPT_M2, LAT, LON, DEPTH, TEMP) %>%
  rename(YEAR = YR) %>%
  filter(YEAR <= 1994) %>%
  filter(nchar(TOW) != 9) %>% # remove the repeated tows, looks like the data are already filled in to the original record
  mutate(ID = paste(CRUISE6, STRATUM, STATION, sep = ".")) %>%
  group_by(YEAR, STRATUM, STATION) %>%
  mutate(NPERTOW = sum(NPERTOW), KGPERTOW = sum(KGPERTOW)) %>%
  ungroup() %>%
  select(-c(LENGRP)) %>%
  distinct() %>%
  group_by(CRUISE6, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  ungroup() 


  ## 1997 - 2023 data ---- 

tow_1997_df <- read.csv("data/NOAA.stock.data/surfclam/TOW_DATA_length_aggregated.CSV") %>% 
  select(CRUISE6, YR, STATION, STRATUM, TOW, LENGRP, NPERTOW, KGPERTOW, ASWEPT_M2, LAT, LON, DEPTH, TEMP) %>%
  rename(YEAR = YR) %>%
  filter(YEAR >= 1997) %>%
  filter(nchar(TOW) != 9) %>% # remove the repeated tows, looks like the data are already filled in to the original record
  mutate(ID = paste(CRUISE6, STRATUM, STATION, sep = ".")) %>%
  select(-c(LENGRP)) %>%
  distinct() %>%
  group_by(CRUISE6, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  ungroup() 


  ## combine data ----
tow_df <- rbind(tow_1982_df, tow_1997_df) %>%
  filter(!is.na(NPERTOW) & !is.na(KGPERTOW)) %>%
  mutate(LON = - LON)

write.csv(tow_df, "results/stratified.mean.indices/surfclam/tow.list.csv", row.names = FALSE) # this is supposed to be tow header only, but I saved everything here


# -------------------------------------------------------------------------------------------------- #


# survey strata area as weighting factors ---------------------------------------------------------------------------------------------

SOUTH.STRATA <- paste0(1:6, "S")
NORTH.STRATA <- paste0(7:12, "S")

stra_area_df <- read.csv("data/NOAA.stock.data/surfclam/STRATA_MEANS.CSV") %>%
  select(STRATUM, AREASQNM) %>%
  mutate(REGION = ifelse(STRATUM %in% SOUTH.STRATA, "SVAtoSNE", "GBK")) %>%
  distinct()

write.csv(stra_area_df, "results/stratified.mean.indices/surfclam/strata.area.csv", row.names = FALSE) # biomass/abundance both included


station_df <- tow_df %>%
  select(YEAR, STRATUM) %>%
  left_join(stra_area_df) %>%
  distinct() %>%
  group_by(YEAR, REGION) %>%
  mutate(WEIGHT = AREASQNM/sum(AREASQNM))

catch_by_tow_df <- tow_df %>%
  left_join(station_df)

write.csv(catch_by_tow_df, "results/stratified.mean.indices/surfclam/total.catch.by.tow.csv", row.names = FALSE) # biomass/abundance both included


# I. all tows calculation ---------------------------------------------------------------------------------------------

## 1.1 calculate mean numbers by stratum ---------------------------------------------------------------------------------------------

mean_N_stratum_df <- catch_by_tow_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((NPERTOW - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(SE_STRATUM = sqrt(VAR_STRATUM)/sqrt(TOTAL_N_STATION)) %>%
  # mutate(VAR_STRATUM = var(EXPCATCHNUM)) %>% # variance by stratum, same as last line
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)
  

## 1.2 calculate stratified mean numbers ---------------------------------------------------------------------------------------------

stratified_mean_N_df <- mean_N_stratum_df %>%
  # select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()

write.csv(stratified_mean_N_df, "results/stratified.mean.indices/surfclam/original.indices.csv", row.names = FALSE)



## 1.3 calculate mean biomass by stratum ---------------------------------------------------------------------------------------------

mean_B_stratum_df <- catch_by_tow_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_B_STRATUM = mean(KGPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((KGPERTOW - MEAN_B_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(SE_STRATUM = sqrt(VAR_STRATUM)/sqrt(TOTAL_N_STATION)) %>%
  # mutate(VAR_STRATUM = var(EXPCATCHNUM)) %>% # variance by stratum, same as last line
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_B_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


## 1.4 calculate stratified mean biomass ---------------------------------------------------------------------------------------------

stratified_mean_B_df <- mean_B_stratum_df %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_B = weighted.mean(MEAN_B_STRATUM, w = WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_B, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_B, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_B) %>%
  ungroup()

write.csv(stratified_mean_N_df, "results/stratified.mean.indices/surfclam/original.biomass.indices.csv", row.names = FALSE)



# II. now exclude the tows overlapped with OWF ---------------------------------------------------------------------------------------------

## identify overlapped tows ----

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

impacted_tow_df <- catch_by_tow_df %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% # add this because the WEA overlay strata are old BTS strata
  filter(!ID.temp %in% unique(AS_QQ_overlay_df$ID.temp)) %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) # total station No needs to be recalculated for each stratum


## check if any stratum is completely removed ----
original_stratum_list <- catch_by_tow_df %>%
  select(YEAR, REGION, STRATUM) %>%
  distinct() %>%
  add_column(TS = "ORIGINAL")

impacted_stratum_list <- impacted_tow_df %>%
  select(YEAR, REGION, STRATUM) %>%
  distinct() %>%
  add_column(TS = "IMPACTED") 

comp_stratum_No <- original_stratum_list %>%
  bind_rows(impacted_stratum_list) %>%
  group_by(YEAR, REGION, TS) %>%
  summarize(N.STRATUM = length(STRATUM)) %>%
  pivot_wider(names_from = TS, values_from = N.STRATUM) %>%
  mutate(MISSING_STRATUM = ORIGINAL - IMPACTED)

missing_stratum_list <- original_stratum_list[,c(1:3)] %>%
  anti_join(impacted_stratum_list[,c(1:3)]) %>%
  group_by(YEAR, REGION) %>%
  mutate(N_STRATUM_LOST = length(unique(STRATUM)))


## 2.1 mean numbers by stratum ----
impacted_mean_N_stratum_df <- impacted_tow_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((NPERTOW - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(SE_STRATUM = sqrt(VAR_STRATUM)/sqrt(TOTAL_N_STATION)) %>%
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM, AREASQNM, WEIGHT)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


## 2.2 stratified mean numbers ----
impacted_stratified_mean_N_df <- impacted_mean_N_stratum_df %>%
  # select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  group_by(YEAR, REGION) %>%
  mutate(REL_WEIGHT = AREASQNM/sum(AREASQNM, na.rm = TRUE)) %>% # this is necessary in case any stratum fully overlap with WEA
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = REL_WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(REL_WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()

write.csv(impacted_stratified_mean_N_df, "results/stratified.mean.indices/surfclam/impacted.indices.csv", row.names = FALSE)


## 2.3 mean biomass by stratum ----
impacted_mean_B_stratum_df <- impacted_tow_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_B_STRATUM = mean(KGPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((KGPERTOW - MEAN_B_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(SE_STRATUM = sqrt(VAR_STRATUM)/sqrt(TOTAL_N_STATION)) %>%
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_B_STRATUM, VAR_STRATUM, AREASQNM, WEIGHT)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


## 2.4 stratified mean biomass ----
impacted_stratified_mean_B_df <- impacted_mean_B_stratum_df %>%
  group_by(YEAR, REGION) %>%
  mutate(REL_WEIGHT = AREASQNM/sum(AREASQNM, na.rm = TRUE)) %>% # this is necessary in case any stratum fully overlap with WEA
  summarize(STRATIFIED_MEAN_B = weighted.mean(MEAN_B_STRATUM, w = REL_WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(REL_WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_B, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_B, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_B) %>%
  ungroup()

write.csv(impacted_stratified_mean_N_df, "results/stratified.mean.indices/surfclam/impacted.indices.biomass.csv", row.names = FALSE)

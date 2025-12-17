library(tidyverse)

BTS_overlay_df <- read.csv("results/BTS_tows_overlay_final.csv")

# survey strata used in assessment -----------------------------------------------------------------------

fall_BTS_stra_list <- c(01010, 01020, 01030, 01040, 01050, 01060, 01070, 01080, 01090, 01100, 01110, 01120, 01130, 01140, 01150,
                        01160, 01170, 01180, 01190, 01200, 01210, 01220, 01230, 01250, 01260, 01610, 01620, 01630, 01640, 01650,
                        01660, 01670, 01680, 01690, 01700, 01710, 01720, 01730, 01740, 01750, 01760, 03020, 03050, 03080, 03110,
                        03140, 03170, 03200, 03230, 03260, 03290, 03320, 03350, 03380, 03410, 03440, 03450, 03460, 03560, 03590,
                        03600, 03610, 03650, 03660) # 63
fall_BTS_stra_list <- substr(fall_BTS_stra_list, 1,4) # remove the first 0 as strata in other excel files can be broken

spring_BTS_stra_list <- c(01010, 01020, 01030, 01040, 01050, 01060, 01070, 01080, 01090, 01100, 01110, 01120, 01130, 01140, 01150,
                          01160, 01170, 01180, 01190, 01200, 01210, 01220, 01230, 01250, 01260, 01610, 01620, 01630, 01640, 01650,
                          01660, 01670, 01680, 01690, 01700, 01710, 01720, 01730, 01740, 01750, 01760, 03020, 03050, 03080, 03110,
                          03140, 03170, 03200, 03230, 03260, 03290, 03320, 03350, 03380, 03410, 03440, 03450, 03460, 03560, 03590,
                          03600, 03610, 03650, 03660) # 63
spring_BTS_stra_list <- substr(spring_BTS_stra_list, 1,4) # remove the first 0 as strata in other excel files can be broken

# ----------------------------------------------------------------------------------------------- #



# survey strata area as weighting factors ---------------------------------------------------------------------------------------------

fall_stra_area_df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/SVDBS_SupportTables/SVDBS_SVMSTRATA.csv")
fall_stra_area_df <- fall_stra_area_df %>%
  filter(stratum %in% fall_BTS_stra_list) %>%
  select(c(stratum, stratum_area)) %>%
  mutate(stratum = as.character(stratum)) %>%
  # mutate(REL_WEIGHT = stratum_area/sum(stratum_area)) %>%
  rename(STRATUM = stratum, STRATUM_AREA = stratum_area) %>%
  add_column(SEASON = "FALL")

spring_stra_area_df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/SVDBS_SupportTables/SVDBS_SVMSTRATA.csv")
spring_stra_area_df <- spring_stra_area_df %>%
  filter(stratum %in% spring_BTS_stra_list) %>%
  select(c(stratum, stratum_area)) %>%
  mutate(stratum = as.character(stratum)) %>%
  # mutate(REL_WEIGHT = stratum_area/sum(stratum_area)) %>%
  rename(STRATUM = stratum, STRATUM_AREA = stratum_area) %>%
  add_column(SEASON = "SPRING")

stra_area_df <- rbind(fall_stra_area_df, spring_stra_area_df); remove(fall_stra_area_df, spring_stra_area_df)

# ----------------------------------------------------------------------------------------------- #





# full tow information ---------------------------------------------------------------------------------------------

fall_tow_df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVSTA.csv")  
fall_tow_df <- fall_tow_df %>%
  select(CRUISE6, STRATUM, TOW, STATION, ID, SHG, TOGA) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  mutate(STRATUM = substr(STRATUM,2,5)) %>%
  add_column(SEASON = "FALL")  %>%
  filter(STRATUM %in% fall_BTS_stra_list) %>% # filter by the species strata list
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))
  
spring_tow_df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVSTA.csv") 
spring_tow_df <- spring_tow_df %>%
  select(CRUISE6, STRATUM, TOW, STATION, ID, SHG, TOGA) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4)),
         STRATUM = as.character(STRATUM)) %>%
  add_column(SEASON = "SPRING")  %>%
  filter(STRATUM %in% spring_BTS_stra_list) %>% # filter by the species strata list
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))

# only day time tow were used, so add an extra criteria here
day_tow_df <- read.csv("data/NOAA.stock.data/longfin.squid/correct data to use/082372_survey_dist_map_fixed_day_tows_only_20250515_prod.csv") %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) %>%
  select(ID, CATCH_WT_CAL) %>%
  mutate(CATCH_WT_CAL = replace_na(CATCH_WT_CAL, 0))

# combine and final filter
full_tow_df <- rbind(fall_tow_df, spring_tow_df) %>%
  # filter(SHG <= 136) %>% # tow evaluation criteria
  # filter(YEAR >= 1976) %>% # longfin squid assessment uses survey data since 1976
  filter(ID %in% day_tow_df$ID) %>%
  left_join(day_tow_df) %>%
  left_join(stra_area_df) %>%
  mutate(BLOCK = paste(YEAR, SEASON, STRATUM, sep = ".")) %>%
  mutate(STRATUM = as.numeric(STRATUM)) %>%
  select(YEAR, SEASON, STRATUM, BLOCK, everything()) # 21238 rows

write.csv(full_tow_df, "results/stratified.mean.indices/longfin.squid/tow.list.csv", row.names = FALSE)

# ----------------------------------------------------------------------------------------------- #


# catch information ----

  ## calculate the calibrated biomass by tow for modeling ----

load(file = "results/BTS_total_catch_by_tow.Rdata")

cat_biomass_df <- cat.df %>%
  filter(SVSPP == 503) %>% #longfin squid 503
  select(-c(CRUISE, SVSPP, CATCHSEX, CATCH_COMMENT, SCIENTIFIC_NAME, EXPCATCHNUM, overlay)) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  filter(YEAR >= 1976) %>%  # longfin squid assessment uses survey data since 1976
  mutate(STRATUM = ifelse(nchar(STRATUM) == 5, str_sub(STRATUM, 2, 5), STRATUM)) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = "."),
         BLOCK = paste(YEAR, SEASON, STRATUM, sep = ".")) %>%
  mutate(CALIBRATION_FACTOR = case_when(
    SEASON == "SPRING" & YEAR > 2008 ~ 1.5308, # only apply calibration factor to >2008
    SEASON == "FALL" & YEAR > 2008 ~ 1.5099,
    TRUE ~ 1))  %>%
  mutate(BIOMASS = EXPCATCHWT/CALIBRATION_FACTOR) %>%
  # mutate(BIOMASS = EXPCATCHWT) %>%
  select(-c(EXPCATCHWT))%>%
  arrange(BLOCK)

cat_biomass_df <- cat_biomass_df %>%
  right_join(full_tow_df) %>% # right join to add zero catch tows
  mutate(BIOMASS = ifelse(is.na(BIOMASS), 0, BIOMASS)) %>%
  select(YEAR, SEASON, STRATUM, BLOCK, everything()) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>% # this is important because number of tows changed due to tow filtering
  ungroup() %>%
  arrange(BLOCK)

write.csv(cat_biomass_df, "results/stratified.mean.indices/longfin.squid/calibrated.total.biomass.by.tow.csv", row.names = FALSE)

# ----------------------------------------------------------------------------------------------- #


# 1.1 calculate mean numbers by stratum ---------------------------------------------------------------------------------------------
# all methods below follow (https://noaa-edab.github.io/survdat/articles/calc_strat_mean.html), used Chris's data CATCH_WT_CAL

mean_N_stratum_df <- cat_biomass_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  # mutate(MEAN_N_STRATUM = sum(BIOMASS, na.rm = TRUE)/TOTAL_N_STATION) %>% # mean within a strata
  # mutate(VAR_STRATUM = sum((BIOMASS - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(MEAN_N_STRATUM = sum(CATCH_WT_CAL, na.rm = TRUE)/TOTAL_N_STATION) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((CATCH_WT_CAL - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  ungroup()

# ----------------------------------------------------------------------------------------------- #



# 1.2 calculate stratified mean numbers ---------------------------------------------------------------------------------------------

stratified_mean_N_df <- mean_N_stratum_df %>%
  select(c(YEAR, SEASON, STRATUM, TOTAL_N_STATION, STRATUM_AREA, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  # select(c(YEAR, SEASON, STRATUM, N.STATION, REL_WEIGHT, MEAN_N_STRATUM)) %>%
  distinct() %>% # downsize the data frame to a minimal without repetitive rows
  group_by(YEAR, SEASON) %>%
  mutate(REL_WEIGHT = STRATUM_AREA/sum(STRATUM_AREA, na.rm = TRUE)) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = REL_WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(REL_WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()

write.csv(stratified_mean_N_df, "results/stratified.mean.indices/longfin.squid/original.indices.csv", row.names = FALSE)

# ----------------------------------------------------------------------------------------------- #


# quick plot to compare its to NOAA official data ---------------------------------------------------------------------------------------------
ref_df <- read.csv("data/NOAA.stock.data/longfin.squid/correct data to use/STOCKEFF_SV_082372_UNIT_NONE_daynight_stratified_indices_20250515.csv")
ref_df <- subset(ref_df, YEAR >= 1976 & METRIC == "Weight (kg)/tow")

ggplot(stratified_mean_N_df) +
  geom_point(aes(x = YEAR, y = STRATIFIED_MEAN_N)) +
  geom_line(aes(x = YEAR, y = STRATIFIED_MEAN_N)) +
  geom_vline(xintercept = 2008) +
  facet_wrap(.~SEASON, nrow = 2, scale = "free_y") +
  geom_point(data = ref_df, aes(x = YEAR, y = DAY_INDEX), color = "red") +
  geom_line(data = ref_df, aes(x = YEAR, y = DAY_INDEX), color = "red") 

# ----------------------------------------------------------------------------------------------- #





# II. now exclude the tows overlapped with OWF ---------------------------------------------------------------------------------------------

## identify tows
impacted_tow_df <- full_tow_df %>%
  filter(!ID %in% unique(BTS_overlay_df$ID)) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) # total station No needs to be recalculated for each stratum


# ## check if any stratum is completely removed
# original_stratum_list <- full_tow_df %>%
#   select(YEAR, SEASON, STRATUM) %>%
#   distinct() %>%
#   add_column(TS = "ORIGINAL")
# 
# impacted_stratum_list <- impacted_tow_df %>%
#   select(YEAR, SEASON, STRATUM) %>%
#   distinct() %>%
#   add_column(TS = "IMPACTED") 
# 
# comp_stratum_No <- original_stratum_list %>%
#   bind_rows(impacted_stratum_list) %>%
#   group_by(YEAR, SEASON, TS) %>%
#   summarize(N.STRATUM = length(STRATUM)) %>%
#   pivot_wider(names_from = TS, values_from = N.STRATUM) %>%
#   mutate(MISSING_STRATUM = ORIGINAL - IMPACTED)
# 
# missing_stratum_list <- original_stratum_list[,c(1:3)] %>%
#   anti_join(impacted_stratum_list[,c(1:3)]) %>%
#   group_by(YEAR, SEASON) %>%
#   mutate(N_STRATUM_LOST = length(unique(STRATUM)))


## filter catch by excluding the impacted tow
BTS_overlay_df <- read.csv("results/BTS_tows_overlay_final.csv")

impacted_cat_df <- cat_biomass_df %>%
  filter(!ID %in% unique(BTS_overlay_df$ID)) %>%
  group_by(BLOCK) %>%
  mutate(TOTAL_N_STATION = length(STATION))

## 2.1 mean numbers by stratum ----
impacted_mean_N_stratum_df <- impacted_cat_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  # mutate(MEAN_N_STRATUM = sum(BIOMASS, na.rm = TRUE)/TOTAL_N_STATION) %>% # mean within a strata
  # mutate(VAR_STRATUM = sum((BIOMASS - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(MEAN_N_STRATUM = sum(CATCH_WT_CAL, na.rm = TRUE)/TOTAL_N_STATION) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((CATCH_WT_CAL - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  ungroup() 

## 2.2 stratified mean numbers ----
impacted_stratified_mean_N_df <- impacted_mean_N_stratum_df %>%
  select(c(YEAR, SEASON, STRATUM, ID, TOTAL_N_STATION, STRATUM_AREA, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  # select(c(YEAR, SEASON, STRATUM, N.STATION, REL_WEIGHT, MEAN_N_STRATUM)) %>%
  group_by(YEAR, SEASON) %>%
  mutate(REL_WEIGHT = STRATUM_AREA/sum(STRATUM_AREA, na.rm = TRUE)) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = REL_WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(REL_WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()

write.csv(impacted_stratified_mean_N_df, "results/stratified.mean.indices/longfin.squid/impacted.indices.csv", row.names = FALSE)





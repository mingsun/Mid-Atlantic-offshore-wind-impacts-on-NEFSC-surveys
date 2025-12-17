library(tidyverse)

# survey strata used in assessment -----------------------------------------------------------------------

fall_BTS_stra_list <- c(01010, 01050, 01090, 01610, 01650, 01690, 01730, 03010, 03020, 03030, 03040, 03050, 03060, 03070, 
                       03080, 03090, 03100, 03110, 03120, 03130, 03140, 03150, 03160, 03170, 03180, 03190, 03200, 03210, 
                       03220, 03230, 03240, 03250, 03260, 03270, 03280, 03290, 03300, 03310, 03320, 03330, 03340, 03350, 
                       03360, 03370, 03380, 03390, 03400, 03410, 03420, 03430, 03440, 03450, 03460, 03470, 03480, 03490, 
                       03500, 03510, 03520, 03530, 03540, 03550, 03560, 03570, 03580, 03590, 03600, 03610) # 68
fall_BTS_stra_list <- substr(fall_BTS_stra_list, 1,4) # remove the first 0 as strata in other excel files can be broken 

spring_BTS_stra_list <- c(01010, 01020, 01030, 01040, 01050, 01060, 01070, 01080, 01090, 01100, 01110, 01120, 01610, 01620, 
                         01630, 01640, 01650, 01660, 01670, 01680, 01690, 01700, 01710, 01720, 01730, 01740, 01750, 01760) # 28
spring_BTS_stra_list <- substr(spring_BTS_stra_list, 1,4) # remove the first 0 as strata in other excel files can be broken 


# survey strata area as weighting factors ---------------------------------------------------------------------------------------------

fall_stra_area_df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/SVDBS_SupportTables/SVDBS_SVMSTRATA.csv")
fall_stra_area_df <- fall_stra_area_df %>%
  filter(stratum %in% fall_BTS_stra_list) %>%
  add_column(SEASON = "FALL")

spring_stra_area_df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/SVDBS_SupportTables/SVDBS_SVMSTRATA.csv")
spring_stra_area_df <- spring_stra_area_df %>%
  filter(stratum %in% spring_BTS_stra_list) %>%
  add_column(SEASON = "SPRING")

stra_area_df <- fall_stra_area_df %>%
  bind_rows(spring_stra_area_df) %>%
  mutate(STRATUM = as.character(stratum),
         STRATUM_AREA = stratum_area) %>%
  select(c(STRATUM, STRATUM_AREA)) %>%
  # mutate(REL_WEIGHT = stratum_area/sum(stratum_area)) %>%
  distinct()

remove(fall_stra_area_df, spring_stra_area_df)

# full tow information ---------------------------------------------------------------------------------------------

# Fall BTS: SHG <= 136, 1963-2008
# Fall BTS: TOGA <= 132x, 2009-2019
# Spring BTS: SHG <= 136, 1963-2008 
# Spring BTS: TOGA <= 132x, 2009-2019

fall_tow_df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVSTA.csv")  
fall_tow_df <- fall_tow_df %>%
  select(CRUISE6, STRATUM, TOW, STATION, ID, SHG, TOGA, STATUS_CODE) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4)),
         STRATUM = substr(STRATUM,2,5)) %>%
  add_column(SEASON = "FALL")  %>%
  filter(STRATUM %in% fall_BTS_stra_list) %>% # filter by the species strata list
  mutate(ID = paste(CRUISE6, STRATUM, TOW, STATION, sep = ".")) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))

spring_tow_df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVSTA.csv") 
spring_tow_df <- spring_tow_df %>%
  select(CRUISE6, STRATUM, TOW, STATION, ID, SHG, TOGA, STATUS_CODE) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4)),
         STRATUM = as.character(STRATUM)) %>%
  add_column(SEASON = "SPRING")  %>%
  filter(STRATUM %in% spring_BTS_stra_list) %>% # filter by the species strata list
  mutate(ID = paste(CRUISE6, STRATUM, TOW, STATION, sep = ".")) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))

full_tow_df <- bind_rows(fall_tow_df, spring_tow_df)
full_tow_df <- full_tow_df %>% 
  left_join(stra_area_df) %>% 
  filter((SHG <= 136 & YEAR >= 1963 & YEAR <= 2008) |
           (TOGA <= 1330  & YEAR >=2009)) %>% # tow evaluation criteria
  filter(YEAR >= 1982) %>% # summer flounder assessment uses survey data since 1982
  mutate(BLOCK = paste(YEAR, SEASON, STRATUM, sep = ".")) %>%
  select(YEAR, SEASON, STRATUM, BLOCK, everything())

# write.csv(full_tow_df, "results/survey.indices/summer.flounder/tow.list.csv", row.names = FALSE)

# catch information ---------------------------------------------------------------------------------------------

 ## calculate the calibrated biomass by tow for modeling


load(file = "results/BTS_total_catch_by_tow.Rdata")

cat_df <- cat.df %>%
  filter(SVSPP == 103) %>%
  select(-c(CRUISE, SVSPP, CATCHSEX, CATCH_COMMENT, SCIENTIFIC_NAME, EXPCATCHNUM, overlay)) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  filter(YEAR >= 1982) %>% # summer flounder assessment uses survey data since 1982
  mutate(STRATUM = ifelse(nchar(STRATUM) == 5, str_sub(STRATUM, 2, 5), STRATUM)) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, STATION, sep = "."),
         BLOCK = paste(YEAR, SEASON, STRATUM, sep = ".")) %>%
  mutate(CALIBRATION_FACTOR = case_when(
    SEASON == "SPRING" & YEAR > 2008 ~ 3.0657, # only apply calibration factor to >2008
    SEASON == "FALL" & YEAR > 2008 ~ 2.1409,
    TRUE ~ 1))  %>%
  mutate(BIOMASS = EXPCATCHWT/CALIBRATION_FACTOR) %>%
  select(-c(EXPCATCHWT))%>%
  arrange(BLOCK)

cat_df <- cat_df %>%
  right_join(full_tow_df) %>% # right join to add zero catch tows
  mutate(BIOMASS = ifelse(is.na(BIOMASS), 0, BIOMASS)) %>%
  select(YEAR, SEASON, STRATUM, BLOCK, everything()) %>%
  arrange(BLOCK)
  
# write.csv(cat_biomass_df, "results/stratified.mean.indices/summer.flounder/calibrated.total.biomass.by.tow.csv", row.names = FALSE)


## add zero catch tows

cat_df <- cat_df %>%
  right_join(full_tow_df) %>% # right join to add zero catch tows
  mutate(BIOMASS = ifelse(is.na(BIOMASS), 0, BIOMASS)) %>%
  select(YEAR, SEASON, STRATUM, BLOCK, everything()) %>%
  arrange(BLOCK)

write.csv(cat_df, "results/stratified.mean.indices/summer.flounder/calibrated.total.biomass.by.tow.csv", row.names = FALSE)

# remove(cat.df)

# calculate mean numbers by stratum ---------------------------------------------------------------------------------------------
# all methods below follow (https://noaa-edab.github.io/survdat/articles/calc_strat_mean.html)

mean_N_stratum_df <- cat_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(BIOMASS)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((BIOMASS - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  # mutate(VAR_STRATUM = var(EXPCATCHNUM)) %>% # variance by stratum, same as last line
  ungroup() 

# calculate stratified mean numbers ---------------------------------------------------------------------------------------------

stratified_mean_N_df <- mean_N_stratum_df %>%
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

write.csv(stratified_mean_N_df, "results/stratified.mean.indices/summer.flounder/original.biomass.indices.csv", row.names = FALSE)

# quick plot to compare its to NOAA official data ---------------------------------------------------------------------------------------------
ref_df <- read.csv("data/NOAA.stock.data/summer.flounder/STOCKEFF_SV_172735_UNIT_NONE_strat_mean.csv")[1:117, c(8,7,10)]
ref_df <- subset(ref_df, YEAR >= 1982)

ggplot(stratified_mean_N_df) +
  geom_point(aes(x = YEAR, y = STRATIFIED_MEAN_N)) +
  geom_line(aes(x = YEAR, y = STRATIFIED_MEAN_N)) +
  geom_vline(xintercept = 2008) +
  facet_wrap(.~SEASON, nrow = 2) +
  geom_point(data = ref_df, aes(x = YEAR, y = INDEX), color = "red") +
  geom_line(data = ref_df, aes(x = YEAR, y = INDEX), color = "red") 



# II. now exclude the tows overlapped with OWF ---------------------------------------------------------------------------------------------

## identify tows
impacted_tow_df <- full_tow_df %>%
  filter(!ID %in% unique(BTS_overlay_df$ID)) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) # total station No needs to be recalculated for each stratum

## check if any stratum is completely removed
original_stratum_list <- full_tow_df %>%
  select(YEAR, SEASON, STRATUM) %>%
  distinct() %>%
  add_column(TS = "ORIGINAL")

impacted_stratum_list <- impacted_tow_df %>%
  select(YEAR, SEASON, STRATUM) %>%
  distinct() %>%
  add_column(TS = "IMPACTED") 

comp_stratum_No <- original_stratum_list %>%
  bind_rows(impacted_stratum_list) %>%
  group_by(YEAR, SEASON, TS) %>%
  summarize(N.STRATUM = length(STRATUM)) %>%
  pivot_wider(names_from = TS, values_from = N.STRATUM) %>%
  mutate(MISSING_STRATUM = ORIGINAL - IMPACTED)

missing_stratum_list <- original_stratum_list[,c(1:3)] %>%
  anti_join(impacted_stratum_list[,c(1:3)]) %>%
  group_by(YEAR, SEASON) %>%
  mutate(N_STRATUM_LOST = length(unique(STRATUM)))

## filter catch by excluding the impacted tow
BTS_overlay_df <- read.csv("results/BTS_tows_overlay_final.csv")

impacted_cat_df <- cat_df %>%
  filter(!ID %in% unique(BTS_overlay_df$ID)) %>%
  group_by(BLOCK) %>%
  mutate(TOTAL_N_STATION = length(STATION))

## mean numbers by stratum
impacted_mean_N_stratum_df <- impacted_cat_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(BIOMASS)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((BIOMASS - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  # mutate(VAR = var(EXPCATCHNUM)) %>% # variance by stratum
  ungroup() 

## stratified mean numbers
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

write.csv(impacted_stratified_mean_N_df, "results/stratified.mean.indices/summer.flounder/impacted.biomass.indices.csv", row.names = FALSE)

# III. randomly exclude the same amounts of overlapped tows using iterations ---------------------------------------------------------------------------------------------

## first identify the tow sample size for each year and season
BTS_overlay_df <- read.csv("results/BTS_tows_overlay_final.csv")

tow_No_year_df <- full_tow_df %>%
  mutate(OVERLAY = ifelse(ID %in% BTS_overlay_df$ID, "TRUE", "FALSE")) %>%
  group_by(YEAR, SEASON) %>%
  summarize(NO_TOTAL_TOW = length(OVERLAY),
            No_IMPACTED_TOW = sum(OVERLAY == "TRUE"),
            No_REMAINING_TOW = NO_TOTAL_TOW - No_IMPACTED_TOW,
            PROP_REMAINING_TOW = round(No_REMAINING_TOW/NO_TOTAL_TOW, 3))

## 1000 iterations to realize the randomization

set.seed(1)
# iter = 2;y = 1982; s = "FALL"

random_df_list <- list()

for (iter in 1:1000) {
  
  ## generate random tow dataset by year and season
  random_tow_df <- data.frame()
  
  for (y in unique(full_tow_df$YEAR)) {
    for (s in unique(full_tow_df$SEASON)) {
      
      temp_tow_df <- subset(tow_No_year_df, YEAR == y & SEASON == s)
      if (nrow(temp_tow_df) == 0) {next}

      random.tow.row.No <- sort(sample(1:temp_tow_df$NO_TOTAL_TOW, temp_tow_df$No_REMAINING_TOW, replace = FALSE))
      
      temp_dataset_df <- full_tow_df %>%
        filter(YEAR == y, SEASON == s) %>% # the full dataset filtered by y and s
        arrange(BLOCK)
      temp_dataset_df <- temp_dataset_df[random.tow.row.No,]
      
      random_tow_df <- rbind(random_tow_df, temp_dataset_df)
    }
  }
  
  remove(temp_dataset_df, temp_tow_df, random.tow.row.No)
  
  ## identify tows
  random_tow_df <- random_tow_df %>%
    group_by(YEAR, SEASON, STRATUM) %>%
    mutate(TOTAL_N_STATION = length(STATION)) # total station No needs to be recalculated for each stratum

  ## filter catch by keeping the random tow
  random_cat_df <- cat_df %>%
    filter(ID %in% unique(random_tow_df$ID)) %>%
    group_by(BLOCK) %>%
    mutate(TOTAL_N_STATION = length(STATION))
  
  ## mean numbers by stratum
  random_mean_N_stratum_df <- random_cat_df %>%
    group_by(YEAR, SEASON, STRATUM) %>%
    mutate(MEAN_N_STRATUM = mean(BIOMASS)) %>% # mean within a strata
    mutate(VAR_STRATUM = sum((BIOMASS - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
    # mutate(VAR = var(EXPCATCHNUM)) %>% # variance by stratum
    ungroup() 
  
  ## stratified mean numbers
  random_stratified_mean_N_df <- random_mean_N_stratum_df %>%
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
    ungroup() %>%
    add_column(ITER = iter)
  
  random_df_list[[iter]] <- random_stratified_mean_N_df
  
  print(paste("iteration", iter, "finished"))
  remove(random_tow_df, random_cat_df, random_mean_N_stratum_df, random_stratified_mean_N_df)
  
}

random_stratified_mean_N_full_df <- bind_rows(random_df_list)

write.csv(random_stratified_mean_N_full_df, "results/stratified.mean.indices/summer.flounder/random.reduction.indices.csv", row.names = FALSE)

library(tidyverse)
# library(plyr)

# station data load and handling ----
station.spring.df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVSTA.csv") 
station.fall.df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVSTA.csv") 

station.spring.df$SEASON <- "SPRING"
station.fall.df$SEASON <- "FALL"

station.df <- rbind(station.spring.df, station.fall.df); remove(station.spring.df, station.fall.df)

# station distribution 
# the below survey stratum fall in the MAB and SNE regions, GOM and GB are excluded
# offshore: 1-12, 61-76, offshore trawl code is 01, 08 (e.g., 01760)
# inshore:  1-11, 45-56, inshore code is 03, 07
# just use the station code extracted from SVDBS_SVMSTRATA

# !!!!!! note that the station list below is only for spatiotemporal modelign in the mid-atlantic, abundance indices spatial range can be larger than this
station.list <- c(paste0(0, c(1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,1110,1120, 1250)), # SNE offshore
                  paste0(sprintf("01%02d", 61:76), 0), # MAB offshore
                  paste0(sprintf("03%02d", 1:14), 0), # SNE inshore
                  paste0(sprintf("03%02d", 15:44), 0), # MAB inshore
                  paste0(sprintf("03%02d", 45:55), 0), # SNE inshore
                  "03910"  # MAB LONG IS SD
                  )

# station.spring.df %>%
#   filter(STRATUM %in% station.list)
       
station.df <- station.df %>%
  select(CRUISE6, STRATUM, TOW, STATION, AREA, SVVESSEL, SVGEAR, EST_YEAR, DECDEG_BEGLAT, DECDEG_BEGLON, SEASON) %>%
  filter(STRATUM %in% station.list| STRATUM %in% substr(station.list, 2, 5))
  # filter(EST_YEAR == 2017)

write.csv(station.df, "results/BTS_tows_list_spatiotemporal_modeling.csv", row.names = FALSE)

# # plot by stratum
# ggplot() +
#   geom_polygon(data = USA.state, aes(x = long, y = lat, group = group), fill= "cornsilk", color = "black") +
#   labs(x = "Longitude", y = "Latitude") +
#   coord_sf(xlim = c(-68, -78.2), ylim = c(33, 42.5)) +
#   geom_polygon(data = my_spdf, aes(x = long, y = lat, group = group), fill = "grey50", alpha = 0.3) +
#   geom_point(data = station.df, aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = STRATUM), size = 2) +
#   theme(panel.grid.major = element_line(color = gray(.85), linetype = "dashed", size = 0.5), 
#         panel.background = element_rect(fill = "aliceblue"), legend.position = "none")

# -------------------------------------------------------------------------------------------------- #



# identify overlapping stations to OWF area ----

# this step is done in Arcgis
# all survey head 180, 30 mins, 2 knots, resulting in 1 knot range 1.852 kilometers
# so we create a 1.852 km buffer for all station, just to keep a distance from the OWF leased area

overlay.spring.df <- read.csv("results/survey overlay/BTS_spring_station_overlay.csv")
overlay.spring.df <- overlay.spring.df %>%
  select(CRUISE6, STRATUM, TOW, STATION, ID, AREA, SVVESSEL, SVGEAR, EST_YEAR, DECDEG_BEGLAT, DECDEG_BEGLON) %>%
  add_column(SEASON = "SPRING")  %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) %>% # all IDs are bad due to scientific notation, generate a new set of ID
  filter(STRATUM %in% substr(station.list, 2, 5)) # substr of the four last digit in station.list, because the first digit in overlay.df was lost due to data saving

overlay.fall.df <- read.csv("results/survey overlay/BTS_fall_station_overlay.csv")
overlay.fall.df <- overlay.fall.df %>%
  select(CRUISE6, STRATUM, TOW, STATION, ID, AREA, SVVESSEL, SVGEAR, EST_YEAR, DECDEG_BEGLAT, DECDEG_BEGLON) %>%
  add_column(SEASON = "FALL") %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) %>% # all IDs are bad due to scientific notation, generate a new set of ID
  filter(STRATUM %in% substr(station.list, 2, 5)) # substr of the four last digit in station.list, because the first digit in overlay.df was lost due to data saving

overlay.df <- rbind(overlay.spring.df, overlay.fall.df); remove(overlay.spring.df, overlay.fall.df)

write.csv(overlay.df, "results/BTS_tows_overlay_final.csv", row.names = FALSE)

# -------------------------------------------------------------------------------------------------- #




# survey effort loss by year in retrospective ----

overlay.no.df <- overlay.df %>%
  group_by(EST_YEAR, SEASON) %>%
  summarise(overlay.n = length(TOW))
  
station.no.df <- station.df %>%
  group_by(EST_YEAR, SEASON) %>%
  summarise(total.n = length(TOW))

overlay.prop.df <- merge(overlay.no.df, station.no.df); remove(overlay.no.df, station.no.df)
overlay.prop.df$ratio <- overlay.prop.df$overlay.n/overlay.prop.df$total.n

se.plot <- ggplot(overlay.prop.df) +
  geom_line(aes(x = EST_YEAR, y = ratio)) +
  geom_point(aes(x = EST_YEAR, y = ratio)) +
  scale_x_continuous(breaks = c(1963:2022)) +
  # scale_x_continuous(breaks = seq(from = 1963, to = 2022, by = 3)) +
  facet_wrap(.~SEASON, nrow = 2) +
  labs(x = "YEAR", y = "proportion of survey effort loss") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5)) +
  coord_cartesian(xlim = c(1965, 2020.5))

png("plots/BTS_survey_effort_loss.png",  width = 16, height = 10, units = 'in', res = 800)
print(se.plot)
dev.off()

write.csv(overlay.prop.df, "results/BTS_survey_effort_loss.csv")

# -------------------------------------------------------------------------------------------------- #



# sample loss size (abundance) in retrospective ----

spring.cat.df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVCAT.csv")[,1:13]
spring.cat.df <- spring.cat.df %>%
  add_column(SEASON = "SPRING") %>%
  mutate(STRATUM = as.numeric(STRATUM)) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) %>% # all IDs are bad due to scientific notation, generate a new set of ID
  transform(overlay = ifelse(ID %in% as.factor(overlay.df$ID), "precluded", "n"))
  
fall.cat.df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVCAT.csv")[,1:13] 
fall.cat.df <- fall.cat.df %>%
  add_column(SEASON = "FALL") %>%
  mutate(STRATUM = as.numeric(STRATUM)) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) %>% # all IDs are bad due to scientific notation, generate a new set of ID
  transform(overlay = ifelse(ID %in% as.factor(overlay.df$ID), "precluded", "n"))

cat.df <- rbind(spring.cat.df, fall.cat.df); remove(spring.cat.df, fall.cat.df)

save(cat.df, file = "results/BTS_total_catch_by_tow.Rdata")

cat.df <- cat.df %>%
  filter(is.na(ID) == FALSE & is.na(SCIENTIFIC_NAME) == FALSE) %>%
  mutate(YEAR = substr(CRUISE6, 1, 4)) %>%
  group_by(SCIENTIFIC_NAME, SVSPP, YEAR, SEASON, overlay) %>%
  summarise(abun = sum(EXPCATCHNUM, na.rm=TRUE)) %>%
  ungroup() %>%
  # ddply(.(SCIENTIFIC_NAME, SVSPP, YEAR, SEASON, overlay), summarize, abun = sum(EXPCATCHNUM, na.rm=TRUE)) %>%
  dplyr::rename("SPECIES" = "SCIENTIFIC_NAME") %>%
  spread(overlay, abun) %>%
  mutate_at(c('precluded', 'n'), ~replace_na(.,0)) %>%
  mutate(total = precluded + n, ratio = precluded/total)

save(cat.df, file = "results/BTS_sample_loss_all_species_by_season.Rdata")

## for all species (total abundance loss over the year) ----

cat.all.df <- cat.df %>%
  group_by(SPECIES, SVSPP) %>%
  summarise(total.precluded = sum(precluded), total = sum(total), precluded.ratio = round(total.precluded/total, 4)) %>%
  ungroup() %>%
  # ddply(.(SPECIES, SVSPP), summarize, total.precluded = sum(precluded), total = sum(total), precluded.ratio = round(total.precluded/total, 4)) %>%
  filter(total > 100) %>%
  arrange(desc(precluded.ratio)) %>%
  mutate(rank = 1:nrow(.))

write.csv(cat.all.df, "results/BTS_sample_loss_all_species.csv", row.names = FALSE)

## for all species (total abundance loss by year) ----

cat.all.by.year.df <- cat.df %>%
  group_by(SPECIES, SVSPP, YEAR) %>%
  summarise(total.precluded = sum(precluded), total = sum(total), precluded.ratio = round(total.precluded/total, 4)) %>%
  ungroup() %>%
  # ddply(.(SPECIES, SVSPP, YEAR), summarize, total.precluded = sum(precluded), total = sum(total), precluded.ratio = round(total.precluded/total, 4)) %>%
  filter(total > 100) %>%
  arrange(desc(precluded.ratio)) %>%
  mutate(rank = 1:nrow(.))

write.csv(cat.all.by.year.df, "results/BTS_sample_loss_all_species_by_year.csv", row.names = FALSE)

# cat.all.df$rank[cat.all.df$SPECIES == "Arctica islandica (ocean quahog)"]; cat.all.df$precluded.ratio[cat.all.df$SPECIES == "Arctica islandica (ocean quahog)"]
# cat.all.df$rank[cat.all.df$SPECIES == "Loligo pealeii (longfin squid)"]; cat.all.df$precluded.ratio[cat.all.df$SPECIES == "Loligo pealeii (longfin squid)"]
# cat.all.df$rank[cat.all.df$SPECIES == "Paralichthys dentatus (summer flounder)"]; cat.all.df$precluded.ratio[cat.all.df$SPECIES == "Paralichthys dentatus (summer flounder)"]
# cat.all.df$rank[cat.all.df$SPECIES == "Spisula solidissima (Atlantic surfclam)"]; cat.all.df$precluded.ratio[cat.all.df$SPECIES == "Spisula solidissima (Atlantic surfclam)"]

### plot all species ----

cat.all.pl <- ggplot(cat.all.df) +
  geom_bar(aes(x = reorder(SPECIES, precluded.ratio), y = precluded.ratio), stat = "identity") +
  labs(y = "survey abundance loss", x = "species") +
  coord_flip() +
  theme_classic()

png("plots/BTS_sample_loss_all_species.png",  width = 10, height = 10, units = 'in', res = 800)
print(cat.all.pl)
dev.off()

# -------------------------------------------------------------------------------------------------- #


## for the four key species (summer flounder: 103, surf clam: 403, quahog 409, longfin squid 503) ----

# cat.four.df <- cat.all.df %>%
#   filter(SVSPP %in% c(103, 403, 409, 503))

cat.four.df <- cat.df %>%
  filter(SVSPP %in% c(103, 403, 409, 503)) %>%
  filter(SVSPP %in% c(103, 503)) 

cat.four.df[cat.four.df == "Arctica islandica (ocean quahog)"] <- "ocean quahog"
cat.four.df[cat.four.df == "Loligo pealeii (longfin squid)"] <- "longfin squid"
cat.four.df[cat.four.df == "Paralichthys dentatus (summer flounder)"] <- "summer flounder"
cat.four.df[cat.four.df == "Spisula solidissima (Atlantic surfclam)"] <- "Atlantic surfclam"

cf.plot <- ggplot(cat.four.df) +
  geom_line(aes(x = as.numeric(YEAR), y = ratio)) +
  geom_point(aes(x = as.numeric(YEAR), y = ratio)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 2023, color = "grey") +
  scale_x_continuous(breaks = seq(from = 1963, to = 2023, by = 4)) +
  facet_grid(SPECIES~SEASON, scale = "free_y") +
  labs(x = "YEAR", y = "proportion of samples loss") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5)) +
  coord_cartesian(xlim = c(1965, 2020.5))

png("plots/BTS_sample_loss_four_species.png",  width = 12, height = 7, units = 'in', res = 800)
print(cf.plot)
dev.off()
  
# -------------------------------------------------------------------------------------------------- #




# cleaning the catch-at-lentgh by tow ----

spring_catch_at_len_df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVLEN.csv")
spring_catch_at_len_df <- spring_catch_at_len_df %>%
  select(-c(LENGTH_COMMENT, SCIENTIFIC_NAME)) %>%
  filter(SVSPP %in% c(103, 403, 409, 503)) %>%
  add_column(SEASON = "SPRING") %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) # all IDs are bad due to scientific notation, generate a new set of ID

fall_catch_at_len_df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVLEN.csv") 
fall_catch_at_len_df <- fall_catch_at_len_df %>%
  select(-c(LENGTH_COMMENT, SCIENTIFIC_NAME)) %>%
  filter(SVSPP %in% c(103, 403, 409, 503)) %>%
  add_column(SEASON = "FALL") %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) # all IDs are bad due to scientific notation, generate a new set of ID

catch_at_len_df <- rbind(spring_catch_at_len_df, fall_catch_at_len_df); remove(spring_catch_at_len_df, fall_catch_at_len_df)

save(catch_at_len_df, file = "results/BTS_catch_at_len_by_tow.Rdata")

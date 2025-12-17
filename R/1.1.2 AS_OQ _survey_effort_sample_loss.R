library(tidyverse)
library(stringr)
# library(plyr)

# station data load and handling ---- (no spring/fall division)
station.df <- read.csv("data/22565_NEFSCAltanticSurfClamOceanQuahogFisheriesIndependentSurveyData/22565_UNION_FSCS_SVSTA.csv")

# no need to differentiate stratam

# identify overlapping stations to OWF area ----

# this step is done in Arcgis
# dredging tows at 3 knots for 5 minutes, nominal tow distance 154m  
# so we create a 154m buffer for all station, just to keep a distance from the OWF leased area

overlay.df <- read.csv("results/survey overlay/AS_OQ_station_overlay_154m.csv")

# try this later? (NEFSC BTS protocol: all survey head 180, 30 mins, 2 knots, resulting in 1 knot range 1.852 kilometers)
# overlay.df <- read.csv("results/survey overlay/AS_OQ_station_overlay_1nm.csv")
overlay.df <- overlay.df %>%
  select(CRUISE6, STRATUM, TOW, STATION, ID, AREA, SVVESSEL, SVGEAR, EST_YEAR, DECDEG_BEGLAT, DECDEG_BEGLON) %>%
  mutate(ID = paste(CRUISE6, STRATUM, STATION, sep = ".")) # all IDs are bad due to scientific notation, generate a new set of ID

write.csv(overlay.df, "results/AS_OQ_tows_overlay_final.csv", row.names = FALSE)

# survey effort loss by year in retrospective ----
overlay.no.df <- overlay.df %>%
  group_by(EST_YEAR) %>%
  summarise(overlay.n = length(TOW))

station.no.df <- station.df %>%
  group_by(EST_YEAR) %>%
  summarise(total.n = length(TOW))

overlay.prop.df <- merge(overlay.no.df, station.no.df); remove(overlay.no.df, station.no.df)
overlay.prop.df$ratio <- overlay.prop.df$overlay.n/overlay.prop.df$total.n

se.plot <- ggplot(overlay.prop.df) +
  geom_line(aes(x = EST_YEAR, y = ratio)) +
  geom_point(aes(x = EST_YEAR, y = ratio)) +
  scale_x_continuous(breaks = unique(overlay.prop.df$EST_YEAR)) +
  labs(x = "YEAR", y = "proportion of survey effort loss") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5)) +
  coord_cartesian(xlim = c(1983, 2021), ylim = c(0.01,0.2))

png("plots/AS_OQ_survey_effort_loss.png",  width = 16, height = 5, units = 'in', res = 800)
print(se.plot)
dev.off()

write.csv(overlay.prop.df, "results/AS_OQ_survey_effort_loss.csv")

# sample loss size (abundance) in retrospective ----

cat.df <- read.csv("data/22565_NEFSCAltanticSurfClamOceanQuahogFisheriesIndependentSurveyData/22565_UNION_FSCS_SVCAT.csv")

cat.df <- cat.df %>%
  mutate(ID = paste(CRUISE6, STRATUM, STATION, sep = ".")) %>% # all IDs are bad due to scientific notation, generate a new set of ID
  transform(overlay = ifelse(ID %in% as.factor(overlay.df$ID), "precluded", "n")) %>%
  filter(is.na(ID) == FALSE & is.na(SCIENTIFIC_NAME) == FALSE) %>%
  mutate(YEAR = substr(CRUISE6, 1, 4)) %>%
  group_by(SCIENTIFIC_NAME, SVSPP, YEAR, overlay) %>%
  summarize(abun = sum(EXPCATCHNUM, na.rm=TRUE)) %>%
  ungroup() %>%
  # ddply(.(SCIENTIFIC_NAME, SVSPP, YEAR, overlay), summarize, abun = sum(EXPCATCHNUM, na.rm=TRUE)) %>%
  dplyr::rename("SPECIES" = "SCIENTIFIC_NAME") %>%
  spread(overlay, abun) %>%
  mutate_at(c('precluded', 'n'), ~replace_na(.,0)) %>%
  mutate(total = precluded + n, ratio = precluded/total)



## for all species (total abundance loss over the year) ----

cat.all.df <- cat.df %>%
  group_by(SPECIES, SVSPP) %>%
  summarize(total.precluded = sum(precluded), total = sum(total), precluded.ratio = round(total.precluded/total, 4)) %>%
  ungroup() %>%
  # ddply(.(SPECIES, SVSPP), summarize, total.precluded = sum(precluded), total = sum(total), precluded.ratio = round(total.precluded/total, 4)) %>%
  filter(total > 100) %>%
  arrange(desc(precluded.ratio)) %>%
  mutate(rank = 1:nrow(.))

write.csv(cat.all.df, "results/AS_OQ_sample_loss_all_species.csv", row.names = FALSE)

### plot all species ----

cat.all.pl <- ggplot(cat.all.df) +
  geom_bar(aes(x = reorder(SPECIES, precluded.ratio), y = precluded.ratio), stat = "identity") +
  labs(y = "survey abundance loss", x = "species") +
  coord_flip() +
  theme_classic()

png("plots/AS_OQ_sample_loss_all_species.png",  width = 10, height = 10, units = 'in', res = 800)
print(cat.all.pl)
dev.off()

## for the four key species (summer flounder: 103, surf clam: 403, quahog 409, longfin squid 503) ----

cat.four.df <- cat.df %>%
  filter(SVSPP %in% c(403, 409)) 

cat.four.df[cat.four.df == "Arctica islandica (ocean quahog)"] <- "ocean quahog"
# cat.four.df[cat.four.df == "Loligo pealeii (longfin squid)"] <- "longfin squid"
# cat.four.df[cat.four.df == "Paralichthys dentatus (summer flounder)"] <- "summer flounder"
cat.four.df[cat.four.df == "Spisula solidissima (Atlantic surfclam)"] <- "Atlantic surfclam"

cf.plot <- ggplot(cat.four.df) +
  geom_line(aes(x = as.numeric(YEAR), y = ratio)) +
  geom_point(aes(x = as.numeric(YEAR), y = ratio)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 2023, color = "grey") +
  scale_x_continuous(breaks = as.numeric(unique(cat.four.df$YEAR))) +
  facet_wrap(.~SPECIES, ncol = 1) +
  labs(x = "YEAR", y = "proportion of samples loss") +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5)) +
  coord_cartesian(xlim = c(1983, 2022)) 

png("plots/AS_OQ_sample_loss_four_species.png",  width = 10, height = 8, units = 'in', res = 800)
print(cf.plot)
dev.off()

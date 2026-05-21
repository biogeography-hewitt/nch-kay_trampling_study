####packages
library(dplyr)
library(tidyverse)
library(car)
library(broom)


####reading in data
raw_data_2024 <- read.csv(file = "data/nch-kay_trampling_data_2024_raw.csv", skip = 6, nrows = 280)

raw_data_2024[is.na(raw_data_2024)] <- 0


####cleaning data

raw_data_2024_clean <- raw_data_2024 %>% 
  pivot_longer(cols = LUEPEC:MICTOL,
               names_to = "species",
               values_to = "pct_cover") %>% 
  select("site.name", "transect.code", "treat.transplot", "treatcode", "altitude", "slope", "species", "pct_cover")


####pooling data within transect

raw_data_2024_pooled <- raw_data_2024_clean %>% 
  group_by(transect.code, species) %>% 
  summarise(pct_cover_analyzed = mean(pct_cover),
            site.name = unique(site.name),
            treatcode = unique(treatcode)
            
  )  


raw_data_2024_pooled$treatcode <- as.factor(raw_data_2024_pooled$treatcode)       
raw_data_2024_pooled$site.name <- as.factor(raw_data_2024_pooled$site.name)



####log transforming the data
raw_data_2024_pooled_log <- raw_data_2024_pooled %>% 
  mutate(pct_cover_analyzed = log(pct_cover_analyzed + 1)) 


####sqrt transforming the data
raw_data_2024_pooled_sqrt <- raw_data_2024_pooled %>% 
  mutate(pct_cover_analyzed = sqrt(pct_cover_analyzed))

aov(species ~ treatcode, data = raw_data_2024_pooled)


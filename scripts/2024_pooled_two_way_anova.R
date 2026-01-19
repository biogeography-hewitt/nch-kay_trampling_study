####packages
library(dplyr)
library(tidyverse)
library(car)
library(broom)


####reading in data
raw_data_2024 <- read.csv(file = "data/nch-kay_trampling_data_2024.csv", skip = 6, nrows = 280)

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



most_abundant_forbs <- c("LUEPEC", "VACOVA", "CASMER", "PHYEMP", "DIPSIT")
most_abundant_grams <- c("CARSPE", "CARNIG", "VAHATR", "DANINT")



#####FORB ANALYSIS#####



####raw data, interaction model
forb_pooled_int_summary_table <- two_way_anova_int_sp_list(most_abundant_forbs, raw_data_2024_pooled)


####raw data, non-interaction model
forb_pooled_summary_table <- two_way_anova_sp_list(most_abundant_forbs, raw_data_2024_pooled)


####log raw data, interaction model
forb_log_pooled_int_summary_table <- two_way_anova_int_sp_list(most_abundant_forbs, raw_data_2024_pooled_log)


####log raw data, non-interaction model
forb_log_pooled_summary_table <- two_way_anova_sp_list(most_abundant_forbs, raw_data_2024_pooled_log)


####sqrt raw data, interaction model
forb_sqrt_pooled_int_summary_table <- two_way_anova_int_sp_list(most_abundant_forbs, raw_data_2024_pooled_sqrt)


####sqrt raw data, non-interaction model
forb_sqrt_pooled_summary_table <- two_way_anova_sp_list(most_abundant_forbs, raw_data_2024_pooled_sqrt)



#####GRAMINOID ANALYSIS#####



####raw data, interaction model
gram_pooled_int_summary_table <- two_way_anova_int_sp_list(most_abundant_grams, raw_data_2024_pooled)


####raw data, non-interaction model
gram_pooled_summary_table <- two_way_anova_sp_list(most_abundant_grams, raw_data_2024_pooled)


####log raw data, interaction model
gram_log_pooled_int_summary_table <- two_way_anova_int_sp_list(most_abundant_grams, raw_data_2024_pooled_log)


####log raw data, non-interaction model
gram_log_pooled_summary_table <- two_way_anova_sp_list(most_abundant_grams, raw_data_2024_pooled_log)


####sqrt raw data, interaction model
gram_sqrt_pooled_int_summary_table <- two_way_anova_int_sp_list(most_abundant_grams, raw_data_2024_pooled_sqrt)


####sqrt raw data, non-interaction model
gram_sqrt_pooled_summary_table <- two_way_anova_sp_list(most_abundant_grams, raw_data_2024_pooled_sqrt)


####summary tables
forb_pooled_int_summary_table
forb_pooled_summary_table
forb_log_pooled_int_summary_table
forb_log_pooled_summary_table
forb_sqrt_pooled_int_summary_table
forb_sqrt_pooled_summary_table


gram_pooled_int_summary_table
gram_pooled_summary_table
gram_log_pooled_int_summary_table
gram_log_pooled_summary_table
gram_sqrt_pooled_int_summary_table
gram_sqrt_pooled_summary_table


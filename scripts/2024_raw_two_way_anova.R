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


#creating an identical copy to use for raw data analysis
raw_data_2024_unpooled <- raw_data_2024_clean %>% 
  mutate(pct_cover_analyzed = pct_cover)


raw_data_2024_unpooled$treatcode <- as.factor(raw_data_2024_clean$treatcode)
raw_data_2024_unpooled$site.name <- as.factor(raw_data_2024_clean$site.name)


####log transforming the data
raw_data_2024_unpooled_log <- raw_data_2024_unpooled %>% 
  mutate(pct_cover_analyzed = log(pct_cover_analyzed + 1))


####sqrt transforming the data
raw_data_2024_unpooled_sqrt <- raw_data_2024_unpooled %>% 
  mutate(pct_cover_analyzed = sqrt(pct_cover_analyzed)) 



most_abundant_forbs <- c("LUEPEC", "VACOVA", "CASMER", "PHYEMP", "DIPSIT")
most_abundant_grams <- c("CARSPE", "CARNIG", "VAHATR", "DANINT")



#####FORB ANALYSIS#####



####raw data, interaction model
forb_raw_int_summary_table <- two_way_anova_int_sp_list(most_abundant_forbs, raw_data_2024_unpooled)


####raw data, non-interaction model
forb_raw_summary_table <- two_way_anova_sp_list(most_abundant_forbs, raw_data_2024_unpooled)


####log raw data, interaction model
forb_log_raw_int_summary_table <- two_way_anova_int_sp_list(most_abundant_forbs, raw_data_2024_unpooled_log)


####log raw data, non-interaction model
forb_log_raw_summary_table <- two_way_anova_sp_list(most_abundant_forbs, raw_data_2024_unpooled_log)


####sqrt raw data, interaction model
forb_sqrt_raw_int_summary_table <- two_way_anova_int_sp_list(most_abundant_forbs, raw_data_2024_unpooled_sqrt)


####sqrt raw data, non-interaction model
forb_sqrt_raw_summary_table <- two_way_anova_sp_list(most_abundant_forbs, raw_data_2024_unpooled_sqrt)



#####GRAMINOID ANALYSIS#####



####raw data, interaction model
gram_raw_int_summary_table <- two_way_anova_int_sp_list(most_abundant_grams, raw_data_2024_unpooled)


####raw data, non-interaction model
gram_raw_summary_table <- two_way_anova_sp_list(most_abundant_grams, raw_data_2024_unpooled)


####log raw data, interaction model
gram_log_raw_int_summary_table <- two_way_anova_int_sp_list(most_abundant_grams, raw_data_2024_unpooled_log)


####log raw data, non-interaction model
gram_log_raw_summary_table <- two_way_anova_sp_list(most_abundant_grams, raw_data_2024_unpooled_log)


####sqrt raw data, interaction model
gram_sqrt_raw_int_summary_table <- two_way_anova_int_sp_list(most_abundant_grams, raw_data_2024_unpooled_sqrt)


####sqrt raw data, non-interaction model
gram_sqrt_raw_summary_table <- two_way_anova_sp_list(most_abundant_grams, raw_data_2024_unpooled_sqrt)



forb_raw_int_summary_table
forb_raw_summary_table
forb_log_raw_int_summary_table
forb_log_raw_summary_table
forb_sqrt_raw_int_summary_table
forb_sqrt_raw_summary_table


gram_raw_int_summary_table
gram_raw_summary_table
gram_log_raw_int_summary_table
gram_log_raw_summary_table
gram_sqrt_raw_int_summary_table
gram_sqrt_raw_summary_table


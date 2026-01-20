#### packages

library(dplyr)
library(tidyverse)
library(car)
library(broom)



#### function for non-interaction model output

two_way_anova_sp_list <- function(sp_list, dataset){
  
  shapiro_vector = c()
  treatcode_pvalues = c()                                         #create blank vectors to store values
  site.name_pvalues = c()
  
  for (i in sp_list) {
    raw_data <- dataset %>%                                       #filters the dataset for specific species using species code
      filter(species == i)                                                  
                     
    model <- aov(pct_cover_analyzed ~ treatcode + site.name,      #creates ANOVA model for analysis
                 data = raw_data)
    
    tidy_model <- tidy(model)                                     #allows for easy p-value extraction
    
    shapiro <- shapiro.test(model$residuals)                      #testing for normality
    
    
    shapiro_vector[i] <-  shapiro$p.value                          #saving p values to vectors

    treatcode_pvalues[i] <- tidy_model$p.value[1]
    site.name_pvalues[i] <- tidy_model$p.value[2]
    
  }
  
  
  summary_table <- data.frame(shapiro_pvalues = shapiro_vector,          #creates summary table of p values
                              treatcode_pvalues = treatcode_pvalues,
                              site.name_pvalues = site.name_pvalues)
      
  return(summary_table)
}



####function for interaction model output

two_way_anova_int_sp_list <- function(sp_list, dataset){
  
  shapiro_vector = c()
  levene_vector = c()
  treatcode_pvalues = c()                                          #create blank vectors to store values
  site.name_pvalues = c()
  interaction_pvalues = c()
  
  for (i in sp_list) {
    raw_data <- dataset %>%                                        #filters the dataset for specific species
      filter(species == i)                                                  
    
    model <- aov(pct_cover_analyzed ~ treatcode*site.name,         #creates model for analysis
                 data = raw_data)
                                                                  
    tidy_model <- tidy(model)                                      #allows for easy p-value extraction
    
    shapiro <- shapiro.test(model$residuals)                       #testing for normality
    levene <- leveneTest(model)                                    #testing for equality of variance
    
    
    shapiro_vector[i] <-  shapiro$p.value                          #saving p values to vectors
    levene_vector[i] <-  levene$`Pr(>F)`[1]
    
    treatcode_pvalues[i] <- tidy_model$p.value[1]
    site.name_pvalues[i] <- tidy_model$p.value[2]
    interaction_pvalues[i] <- tidy_model$p.value[3]
    
  }
  
  
  summary_table <- data.frame(shapiro_pvalues = shapiro_vector,       #creates summary table of p values
                              levene_pvalues = levene_vector,
                              treatcode_pvalues = treatcode_pvalues,
                              site.name_pvalues = site.name_pvalues,
                              interaction_pvalues = interaction_pvalues)
  
  return(summary_table)
}



########EXAMPLE########


#reading in data

raw_data_2024 <- read.csv(file = "data/nch-kay_trampling_data_2024_raw.csv", skip = 6, nrows = 280)

raw_data_2024[is.na(raw_data_2024)] <- 0       #changes NA values to 0



####pivoting data to long form, selecting columns of interest

raw_data_2024_clean <- raw_data_2024 %>%       
  pivot_longer(cols = LUEPEC:MICTOL,
               names_to = "species",
               values_to = "pct_cover") %>% 
  select("site.name", "transect.code", "treat.transplot", "treatcode", "altitude", "slope", "species", "pct_cover")


####taking mean percent cover for each species within each transect

raw_data_2024_pooled <- raw_data_2024_clean %>% 
  group_by(transect.code, species) %>% 
  summarise(pct_cover_analyzed = mean(pct_cover),
            site.name = unique(site.name),       
            treatcode = unique(treatcode))                      
            
    
raw_data_2024_pooled$treatcode <- as.factor(raw_data_2024_pooled$treatcode)
raw_data_2024_pooled$site.name <- as.factor(raw_data_2024_pooled$site.name)



####Creating function arguments

ex_sp_list <-  c("LUEPEC", "VACOVA", "CASMER", "PHYEMP", "DIPSIT")
ex_dataset <- raw_data_2024_pooled


####testing interaction function

ex_int_output <- two_way_anova_int_sp_list(ex_sp_list, ex_dataset)
ex_int_output


####testing non-interaction function

ex_output <- two_way_anova_sp_list(ex_sp_list, ex_dataset)
ex_output


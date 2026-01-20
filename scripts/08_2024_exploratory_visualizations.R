####packages
library(dplyr)
library(tidyverse)
library(car)
library(broom)
library(gridExtra)

####ALSO REQUIRED - run the two_way_anova_functions script####


####reading in data
raw_data_2024 <- read.csv(file = "data/nch-kay_trampling_data_2024_raw.csv", skip = 6, nrows = 280)


raw_data_2024[is.na(raw_data_2024)] <- 0

####cleaning data

raw_data_2024_clean <- raw_data_2024 %>% 
  pivot_longer(cols = LUEPEC:MICTOL,
               names_to = "species",
               values_to = "pct_cover") %>% 
  select("site.name", "transect.code", "treat.transplot", "treatment", "altitude", "slope", "species", "pct_cover") %>% 
  mutate(treatment = factor(treatment, levels = c("near", "far")))


####visualizations####
sp_list <- c("LUEPEC", "VACOVA", "PHYEMP", "CASMER","HIETRI", "VAHATR", "CARSPE", "CARNIG", "littercover", "soil")   #list of species to show in visualizations

for (i in sp_list){     #loops through list, creates a visualization for each species, assigns it to an object called each species' code
  assign(i,
         ggplot(data = filter(raw_data_2024_clean, species == i), 
                aes(x = site.name, y = pct_cover, fill = treatment)) +
         geom_boxplot() +
         scale_fill_manual(values = c("#b10026", "#0570b0"),
                           breaks = c("near", "far")) +
         labs(x = "Site",
              y = "Percent cover (%)",
              color = "Transect",
              title = i) +
         theme_bw() +
         theme(axis.title.x = element_text(size = 12),
                 axis.title.y = element_text(size = 12),
                 axis.text = element_text(face = "bold", size = 8),
                 legend.title = element_text(face = "bold", size = 10),
                 legend.text = element_text(size = 10))
  )
  
}

LUEPEC

sp_site_treatment_figure <- grid.arrange(LUEPEC, VACOVA, PHYEMP, CASMER, HIETRI, VAHATR, CARSPE, CARNIG, littercover, soil, 
                                          nrow = 5)

ggsave(filename = "sp_site_treatment_figure.png", 
       plot = sp_site_treatment_figure,
       path = "figure",
       width = 10,
       height = 16)



  
  

  
####packages
library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(broom)

####ALSO REQUIRED - run the two_way_anova_functions script####


####reading in data
raw_data_2024 <- read.csv(file = "data/nch-kay_trampling_data_2024_raw.csv", skip = 6, nrows = 280)


raw_data_2024[is.na(raw_data_2024)] <- 0

####cleaning data, logging pct cover

raw_data_2024_clean <- raw_data_2024 %>% 
  pivot_longer(cols = LUEPEC:MICTOL,
               names_to = "species",
               values_to = "pct_cover") %>% 
  select("site.name", "transect.code", "treat.transplot", "treatment","treatcode", "altitude", "slope", "species", "pct_cover") %>% 
  mutate(treatment = factor(treatment, levels = c("near", "far")),
         log_pct_cover = log(pct_cover))

####creating dataframe for combined soil/rock plot



####visualizations####
sp_list <- c("LUEPEC", "VACOVA", "PHYEMP", "CASMER","HIETRI", "VAHATR", "CARSPE", "CARNIG", "littercover", "soil")   #list of species to show in visualizations

for (i in sp_list){     #loops through list, creates a visualization for each species, assigns it to an object called each species' code
  assign(i,
         ggplot(data = filter(raw_data_2024_clean, species == i), 
                aes(x = site.name, y = pct_cover, fill = treatment)) +
         geom_boxplot() +
         scale_fill_manual(values = c("#FF0000", "#ACD6E4"),
                           breaks = c("near", "far")) +
         scale_x_discrete(labels = NULL) +
         labs(x = NULL,
              y = NULL,
              fill = "Treatment") +
         theme_bw() +
         theme(axis.title.x = element_text(size = 14),
               axis.title.y = element_text(size = 14),
               axis.text = element_text(face = "bold", size = 10),
               legend.title = element_text(face = "bold", size = 14),
               legend.text = element_text(size = 14),
               axis.line = element_line(color='black'),
               plot.background = element_blank(),
               plot.title = element_text(face = "italic", size = 16),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
  )
  
}


#adding axis and plot titles
PHYEMP <- PHYEMP +
  ylab("Cover (%)") +
  ggtitle("Phylodoce empetriformis")

CASMER <- CASMER +
  ggtitle("Cassiope mertensiana")

VACOVA <- VACOVA +
  ylab("Cover (%)") +
  ggtitle("Vaccinium ovalifolium")

LUEPEC <- LUEPEC +
  ggtitle("Luetkea pectinata")

HIETRI <- HIETRI +
  ylab("Cover (%)") +
  ggtitle("Hieracium tristis")

VAHATR <- VAHATR +
  ggtitle("Vahlodea atropurpurea")

CARSPE <- CARSPE +
  scale_x_discrete(labels = c("BT", "PR", "TM")) +
  ylab("Cover (%)") +
  ggtitle("Carex spectabilis")

CARNIG <- CARNIG +
  scale_x_discrete(labels = c("BT", "PR", "TM")) +
  ggtitle("Carex nigricans")


#combining into a panel
sp_site_treatment_figure <- ggarrange(PHYEMP, CASMER, VACOVA, LUEPEC, HIETRI, VAHATR, CARSPE, CARNIG, nrow = 4, ncol = 2, common.legend = T, legend = "top")

sp_site_treatment_figure

ggsave(filename = "sp_site_treatment_figure.png", 
       plot = sp_site_treatment_figure,
       path = "figure",
       width = 10,
       height = 16)


#ground cover panel
littercover <- littercover +
  scale_x_discrete(labels = c("BT", "PR", "TM")) +
  ylab("Cover (%)")

soil <- soil +
  scale_x_discrete(labels = c("BT", "PR", "TM")) +
  ylab("Cover (%)")

ground_cover_treatment_figure <- ggarrange(littercover, soil, nrow = 1, ncol = 2, common.legend = T, legend = "bottom")

ground_cover_treatment_figure

ggsave(filename = "ground_cover_treatment_figure.png", 
       plot = ground_cover_treatment_figure,
       path = "figure",
       width = 10,
       height = 4)

#### ANOVAs for NJB paper
littercover_data <- raw_data_2024_clean %>% 
  filter(species == "littercover")

tidy(aov(pct_cover ~ treatcode + site.name, data = littercover_data))
tidy(aov(pct_cover ~ treatcode*site.name, data = littercover_data))

soil_data <- raw_data_2024_clean %>% 
  filter(species == "soil")

tidy(aov(pct_cover ~ treatcode + site.name, data = soil_data))
tidy(aov(pct_cover ~ treatcode*site.name, data = soil_data))
  

####Same as above, but for logged data

for (i in sp_list){     #loops through list, creates a visualization for each species, assigns it to an object called each species' code + 'log'
  assign(paste(i, "log"),
         ggplot(data = filter(raw_data_2024_clean, species == i), 
                aes(x = site.name, y = log_pct_cover, fill = treatment)) +
           geom_boxplot() +
           scale_fill_manual(values = c("#FF0000", "#0570b0"),
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

`LUEPEC log`


sp_site_treatment_figure <- grid.arrange(`LUEPEC log`, `VACOVA log`, `PHYEMP log`, `CASMER log`, `HIETRI log`, `VAHATR log`, `CARSPE log`, `CARNIG log`, `littercover log`, `soil log`, 
                                         nrow = 5)

ggsave(filename = "LOG_sp_site_treatment_figure.png", 
       plot = sp_site_treatment_figure,
       path = "figure",
       width = 10,
       height = 16)
  
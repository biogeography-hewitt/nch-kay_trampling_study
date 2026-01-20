####packages
library(dplyr)
library(tidyverse)
library(gridExtra)


####reading in data
raw_data_2025 <- read.csv(file = "data/nch-kay_trampling_data_2025_raw.csv", skip = 4, nrows = 121)

raw_data_2025[is.na(raw_data_2025)] <- 0 


####cleaning data,

raw_data_2025_clean <- raw_data_2025 %>% 
  pivot_longer(cols = antmed:cypto.crust.w.lichen.present,
               names_to = "species",
               values_to = "pct_cover") %>% 
  select("species", "transect.no",  "site.name", "pct_cover", "transect.code", "quad")


####filtering for only perpendicular transects, calculating distance from trail

raw_data_2025_perp <- raw_data_2025_clean %>% 
  filter(transect.no %in% c(29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5)) %>% 
  mutate(dist_to_trail = case_when(
    transect.no %in% c(29.5, 31.5, 33.5, 35.5) ~ quad - 0.5,
    transect.no %in% c(30.5, 32.5, 34.5, 36.5) ~ 2 * quad + 5
  )
  )


####filtering for 2 plots on near-trail parallel transects closest to perpendicular transects, adding together to get 1m^2 coverage, setting distance to trail equal to 0

raw_data_2025_par <- raw_data_2025_clean %>% 
  filter(transect.no %in% c(29, 31, 33, 35),
         quad %in% c(1, 2)) %>% 
  group_by(species, transect.no) %>% 
  summarise(site.name = unique(site.name),
            pct_cover = sum(pct_cover),
            transect.code = case_when(
              unique(transect.code) == "HC-29" ~ "HC-P-29.5",
              unique(transect.code) == "HC-31" ~ "HC-P-31.5",
              unique(transect.code) == "HC-33" ~ "HC-P-33.5",
              unique(transect.code) == "HC-35" ~ "HC-P-35.5" 
            ),
            quad = 0,
            dist_to_trail = 0)


####calculating mean percent cover

pooled_data_2025_perp <- raw_data_2025_perp %>% 
  group_by(species, dist_to_trail) %>% 
  summarise(mean_pct_cover = mean(pct_cover),
            se_pct_cover = sd(pct_cover)/sqrt(4),
            dist_to_trail = unique(dist_to_trail))



pooled_data_2025_par <- raw_data_2025_par %>% 
  group_by(species) %>% 
  summarise(dist_to_trail = 0,
            mean_pct_cover = mean(pct_cover),
            se_pct_cover = sd(pct_cover)/sqrt(4))


####combining dataframes of perpendicular and parallel to get 0 - 15.5m coverage

# combine by row
raw_data_2025_comb <- rbind(raw_data_2025_par, raw_data_2025_perp)

pooled_data_2025_comb <- rbind(pooled_data_2025_par, pooled_data_2025_perp)

# update ID column
raw_data_2025_comb$ID <- 1:nrow(raw_data_2025_comb) 

pooled_data_2025_comb$ID <- 1:nrow(pooled_data_2025_comb)

#combining transect codes to make lines continuous
raw_data_2025_comb <- raw_data_2025_comb %>% 
  mutate(transect.code = case_when(
    transect.code %in% c("HC-P-29.5", "HC-P-30.5") ~ "HC-P-29.5/30.5",
    transect.code %in% c("HC-P-31.5", "HC-P-32.5") ~ "HC-P-31.5/32.5",
    transect.code %in% c("HC-P-33.5", "HC-P-34.5") ~ "HC-P-33.5/34.5",
    transect.code %in% c("HC-P-35.5", "HC-P-36.5") ~ "HC-P-35.5/36.5"
  ))




####visualizations####

sp_list <- c("luepec", "vacova", "phyemp", "casmer","hietri", "dipsit", "vahatr", "carspe", "carnig", "dicranum.spp.")   #list of species to show in visualizations

for (i in sp_list){     #loops through list, creates a visualization for each species, assigns it to an object called each species' code
  assign(i,
         ggplot(data = filter(raw_data_2025_comb, species == i), 
                aes(x = dist_to_trail, y = pct_cover, colour = transect.code)) +
           geom_line(linewidth = 1, alpha = 0.95) +
           scale_colour_manual(values = c("HC-P-29.5/30.5" = "#0570b0", "HC-P-31.5/32.5" = "#238443", "HC-P-33.5/34.5" = "#6a51a3", "HC-P-35.5/36.5" = "#b10026"),
                               labels = c("29.5/30.5", "31.5/32.5", "33.5/34.5", "35.5/36.5")) +
           scale_x_continuous(breaks = c(0, 0.5, 1.5, 2.5, 3.5, 4.5, 7, 9, 11, 13, 15),
                              labels = c(0, 0.5, 1.5, 2.5, 3.5, 4.5, 7, 9, 11, 13, 15))+
           labs(x = "Distance from trail (m)",
                y = "Percent cover (%)",
                color = "Transect",
                title = i) +
           theme_bw() +
           theme(axis.title.x = element_text(size = 12),
                 axis.title.y = element_text(size = 12),
                 axis.text = element_text(face = "bold", size = 8),
                 legend.title = element_text(face = "bold", size = 10),
                 legend.text = element_text(size = 10),
                 panel.grid = element_blank())
  )
}
luepec

carnig <- carnig +
  scale_y_continuous(breaks = c(0, 0.50, 1.00),
                     labels = c(0, 0.5, 1))
carnig





raw_pct_cover_transects_panel <- grid.arrange(luepec, vacova, phyemp, casmer, hietri, dipsit, vahatr, carspe, carnig, dicranum.spp.,
             nrow = 5)

ggsave(filename = "raw_pct_cover_transects_panel.png", 
       plot = raw_pct_cover_transects_panel,
       path = "figure",
       width = 12,
       height = 10)


for (i in sp_list){   #same as above, for pooled data this time
  assign(paste(i, "1"),
         ggplot(data = filter(pooled_data_2025_comb, species == i), 
                aes(x = dist_to_trail, y = mean_pct_cover)) +
           geom_line(linewidth = 1, alpha = 0.95) +
           geom_errorbar(aes(ymax = mean_pct_cover + se_pct_cover, 
                             ymin = mean_pct_cover - se_pct_cover)) +
           scale_x_continuous(breaks = c(0, 0.5, 1.5, 2.5, 3.5, 4.5, 7, 9, 11, 13, 15),
                              labels = c(0, 0.5, 1.5, 2.5, 3.5, 4.5, 7, 9, 11, 13, 15))+
           labs(x = "Distance from trail (m)",
                y = "Mean percent cover (%)",
                title = i) +
           theme_bw() +
           theme(axis.title.x = element_text(size = 12),
                 axis.title.y = element_text(size = 12),
                 axis.text = element_text(face = "bold", size = 8),
                 legend.title = element_text(face = "bold", size = 10),
                 legend.text = element_text(size = 10),
                 panel.grid = element_blank(),
                 panel.border = element_rect(size = 1.5)))
}


`luepec 1`



pooled_pct_cover_panel <- grid.arrange(`luepec 1`, `vacova 1`, `phyemp 1`, `casmer 1`, `hietri 1`, 
             `dipsit 1`, `vahatr 1`, `carspe 1`, `carnig 1`, `dicranum.spp. 1`,
             nrow = 5)

ggsave(filename = "pooled_pct_cover_panel.png", 
       plot = pooled_pct_cover_panel,
       path = "figure",
       width = 12,
       height = 10)

# Data on trampling effects of near and far from trail
# Purpose to wrangle diversity info
# Sept 2024 NH;
# Updated Oct 1 2024_Go to line 86 end for Oct 17 Script; 
# Updated Mar 13 2025 for ggplot versions of beta diversity NMDS plots, and 
# added plot with points colored by site

# Clears the workspace
rm(list=ls())

# Ask R where it is right now
getwd()

setwd()
## correct format is file
# ground_cover.csv

# install.packages("ggplot2")
# install.packages("tidyverse")
# install.packages("dplyr")

library(ggplot2)
library(dplyr)
library(tidyverse)

## Save pooled data - skip and go to line 94 if you don't want to save pooled data

# BRING IN THE DATA

# load the dataset. You can either use file.choose command
df <- read.csv(file.choose()) # I chose "species_matrix_factorcolumns.csv" from my RScripts working directory
# or use the command read.csv() 
# df = read.csv("___.csv")

# df = read.csv("/Users/ninahewitt/Desktop/Teaching/Geog 374 Statistics in Geography/R Scripts/[filename_here].csv")

## LOOK AT OUR DATA
df          # look at all your data
head(df)    # look at the first few entries
tail(df)    # look at the last few entries
str(df)     # look at the overall structure

# Convert all columns to numeric except 1st 6 - which are as factor 
# Load necessary library
library(dplyr)

df <- df %>%
  mutate(across(1:6, as.factor)) %>%
  mutate(across(7:ncol(df), as.numeric))

# Display the structure of the modified data frame
str(df)

# Write csv
# write.csv(df, file = "name.csv")

## don't need (already acheived) Replace NA with zeros for numeric columns
# df <- mydata %>%
#  mutate(across(where(is.numeric), ~replace_na(., 0)))
## or
#mydata[numeric_vars] <- lapply(mydata[numeric_vars], function(x) ifelse(is.na(x), 0, x))

# Display the modified data frame
head(df)

# Compute mean values for each transect
mean_values <- df %>%
  group_by(transnum) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# works, but next time use transno/transnum as it does not order them properly

# command to AI turor: "Set the first 6 columns to "as factor" and the remaining (from column 7) to numeric. Create boxplots by near vs far for values in columns 7 to x (antmed to vervor

 
write.csv(mean_values, file = "species_means.csv")


# Write csv
# write.csv(df, file = "name.csv")

## don't need (already acheived) Replace NA with zeros for numeric columns
# df <- mydata %>%
#  mutate(across(where(is.numeric), ~replace_na(., 0)))
## or
#mydata[numeric_vars] <- lapply(mydata[numeric_vars], function(x) ifelse(is.na(x), 0, x))

# Display the modified data frame
head(df)




######### Oct 17 2025 workflow
# Steps recorded in Word document:
#"Workflow R code to analyse species diversity composition data.docx"

# Clears the workspace
rm(list=ls())

# Ask R where it is right now
getwd()

setwd()
## correct format is file
# ground_cover.csv

library(ggplot2)
library(dplyr)
# install.packages("tidyverse")
library(tidyverse)

# BRING IN THE DATA

# load the dataset. You can either use file.choose command
df <- read.csv(file = "data/species_matrix_factorcolumns.csv") # I chose "species_matrix_factorcolumns.csv" (from Sept 13, 2025) from my RScripts working directory
# can also choose "species_means.csv" for means only
# for just mosses in Connor's transects: "species_matrix_mosses_CGW_transects_nocrust.csv"
# for just vascular: "species_matrix_factorcolumns_forbs_grasses_shrubs_trees.csv"
# for just vascular, no graminoids: "species_matrix_factorcolumns_forbs_shrubs_trees.csv"
# or use the command read.csv() 
# df = read.csv("___.csv")

# df = read.csv("/Users/ninahewitt/Desktop/Teaching/Geog 374 Statistics in Geography/R Scripts/[filename_here].csv")

## LOOK AT OUR DATA
df          # look at all your data
head(df)    # look at the first few entries
tail(df)    # look at the last few entries
str(df)     # look at the overall structure

# Convert all columns to numeric except 1st 6 - which are as factor 
# Load necessary library
library(dplyr)

df <- df %>%
  mutate(across(1:6, as.factor)) %>%
  mutate(across(7:ncol(df), as.numeric))
##OR if for species means dataset, only first 3 columns contain factor, so write:
# df <- df %>%
#  mutate(across(1:3, as.factor)) %>%
#  mutate(across(4:ncol(df), as.numeric))

# Display the structure of the modified data frame
str(df)

# Load necessary packages
# install.packages("lattice")
# install.packages("permute")
# install.packages("vegan")
# install.packages("dplyr")
# install.packages("ggplot2")

library(lattice)
library(permute)
library(dplyr)
library(ggplot2)
library(vegan)
library(BiodiversityR)

## 2. Alpha Diversity (Shannon, Simpson)
# Alpha Diversity: Measures diversity within each transect, calculated using Shannon and Simpson indices
## calculate common diversity indices for each transect and compare the results.
 
# Calculate Shannon and Simpson diversity indices
df$Shannon <- diversity(df[, 7:ncol(df)], index = "shannon")
df$Simpson <- diversity(df[, 7:ncol(df)], index = "simpson")
# or for species means dataframe:
# df$Shannon <- diversity(df[, 4:ncol(df)], index = "shannon")
# df$Simpson <- diversity(df[, 4:ncol(df)], index = "simpson")

# Summarize by site (Near vs Far)
div_summary <- df %>%
  group_by(treatment) %>%
  summarise(
    Shannon_mean = mean(Shannon),
    Shannon_sd = sd(Shannon),
    Simpson_mean = mean(Simpson),
    Simpson_sd = sd(Simpson)
  )

print(div_summary)

# Results for raw data with 280 rows were:
# A tibble: 2 × 5
#treatment Shannon_mean Shannon_sd Simpson_mean Simpson_sd
# <fct>            <dbl>      <dbl>        <dbl>      <dbl>
# 1 far            1.37      0.467        0.658      0.185
# 2 near           1.10      0.518        0.615      0.228
 
# Results for "species_means.csv" data (means pooled by transect)
# A tibble: 2 × 5
# A tibble: 2 × 5
# treatment Shannon_mean Shannon_sd Simpson_mean Simpson_sd
# <fct>            <dbl>      <dbl>        <dbl>      <dbl>
# 1 far             1.82       0.411        0.761     0.113 
# 2 near            1.77       0.224        0.772     0.0488

# Compare diversity between sites
t_test_shannon <- t.test(Shannon ~ treatment, data = df)
t_test_simpson <- t.test(Simpson ~ treatment, data = df)

print(t_test_shannon) # Results - significantly different t = 4.5484, df = 275.03, p-value = 8.116e-06
print(t_test_simpson) # Results - not significant t = 1.7517, df = 266.45, p-value = 0.08097]

# When species' means per transect are used, results NOT significant 

# 3. Beta Diversity (Bray-Curtis, Jaccard)

# Use dissimilarity measures and ordination (NMDS) to assess and visualize differences in species composition.

# Compute Bray-Curtis dissimilarity matrix
species_dat <- df[, 7:ncol(df)]

# species_dat <- df[, 4:ncol(df)] # if for species means dataframe

bray_curtis <- vegdist(species_dat, method = "bray")
bray_curtis # present results as a matrix of dissimilarity

# Perform NMDS
nmds <- metaMDS(bray_curtis, k = 2) # ran 20 iterations; best solution was not repeated
nmds # see results - 

# Plot NMDS
plot(nmds, type = "t")
ordiplot(nmds, display = "site", type = "n")
points(nmds, display = "site", col = ifelse(df$treatment == "near", "blue", "red"), pch = 19) # doesn't work
legend("bottomleft", legend = c("near", "far"), col = c("blue", "red"), pch = 19)

# Plot NMDS but add envelopes [code added Jan 22, 2026, NH]
#####  TBD later today!
plot <-  plot(nmds, type = "t")
ordiplot(nmds, display = "site", type = "n")
points(nmds, display = "site", col = ifelse(df$treatment == "near", "blue", "red"), pch = 19) # doesn't work
legend("bottomleft", legend = c("near", "far"), col = c("blue", "red"), pch = 19)

#draws polygon around groups
for (i in unique (df$treatment)) ordihull (nmds, groups = df$treatment, show.group = i, col = c("red", "blue"), draw = 'polygon', label = F)


# Plot NMDS but number the points
plot(nmds, type = "t")
ordiplot(nmds, display = "treatment", type = "n") # gives transect number
points(nmds, display = "site", col = ifelse(df$treatment == "near", "blue", "red"), pch = 19) # doesn't work
legend("bottomleft", legend = c("near", "far"), col = c("blue", "red"), pch = 19)


# PERMANOVA (Permutational MANOVA, McArdle and Anderson 2001; Anderson 2001) to test differences in species composition between treatments
# is an Analysis of variance using distance matrices - for partitioning distance matrices among sources of variation and fitting linear models
# (e.g., factors, polynomial regression) to distance matrices; uses a permutation test with pseudo-
# F ratios.
permanova <- adonis2(bray_curtis ~ treatment, data = df) # doesn't make sense
print(permanova)


###### Step: color by site with new color scheme:
## Repeat steps on lines 86 to 184

# Load additional package for better color palettes
library(ggplot2)

# Perform NMDS as before
# Compute Bray-Curtis dissimilarity matrix
species_dat <- df[, 7:ncol(df)]

# species_dat <- df[, 4:ncol(df)] # if for species means dataframe

bray_curtis <- vegdist(species_dat, method = "bray")
bray_curtis # present results as a matrix of dissimilarity
nmds <- metaMDS(bray_curtis, k = 2) # ran 20 iterations; best solution was not repeated
# best was from try 18
nmds # see results - 

# Extract NMDS scores
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$site <- df$site  # Add site info

# Relabel the Site column
nmds_scores$site <- recode(df$site, 
                           `1` = "TM", 
                           `2` = "BT", 
                           `3` = "PR")

# Plot NMDS with ggplot2, colored by site
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = site)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "green", "orange")) +  # Customize colors
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove all gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Add a clean border
    legend.position = "right"
  ) +
  labs(title = "NMDS Plot Colored by Site",
       x = "NMDS1",
       y = "NMDS2",
       color = "site") 

### Repeat Site plot labelled with plot numbers
## Repeat steps in previous lines 223-43 and lines 86 to 184 if haven't

# Relabel the Site column
nmds_scores$site <- recode(df$site, 
                           `1` = "TM", 
                           `2` = "BT", 
                           `3` = "PR")

# Add a column for transect numbers (assuming sequential from 0 to n-1)
nmds_scores$transnum <- 0:(nrow(nmds_scores) - 1)

# Plot NMDS with ggplot2, colored by site and labelle by transect and plot number
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = site)) +
  geom_point(size = 3) +
  geom_text(aes(label = transnum), vjust = -1, size = 3) +  # Add transect numbers
  scale_color_manual(values = c("blue", "green", "orange")) +  # Customize colors
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove all gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Add a clean border
    legend.position = "right"
  ) +
  labs(title = "NMDS Plot Colored by Site",
       x = "NMDS1",
       y = "NMDS2",
       color = "site") 



# PERMANOVA to test differences in species composition between site
permanova2 <- adonis2(bray_curtis ~ site, data = df) 
print(permanova2)

# Try same thing for Treatment - plot is nicer
# Load additional package for better color palettes
library(ggplot2)

# Perform NMDS as before
# Compute Bray-Curtis dissimilarity matrix
species_dat <- df[, 7:ncol(df)]

# species_dat <- df[, 4:ncol(df)] # if for species means dataframe

bray_curtis <- vegdist(species_dat, method = "bray")
bray_curtis # present results as a matrix of dissimilarity
nmds <- metaMDS(bray_curtis, k = 2) # ran 20 iterations; best solution was not repeated
nmds # see results - 

#### Now repeat for treatment with this style of graph:
# Extract NMDS scores and add treatment info
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$treatment <- df$treatment  # Ensure this matches your treatment column name

# Plot NMDS with ggplot2, colored by treatment
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "darkorange")) +  # Customize colors
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove all gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Add a clean border
    legend.position = "right"
  ) +
  labs(title = "NMDS Plot Colored by Treatment",
       x = "NMDS1",
       y = "NMDS2",
       color = "treatment")



# PERMANOVA to test differences in species composition between site
permanova <- adonis2(bray_curtis ~ treatment, data = df) # 
print(permanova)



### Redo with numbers

nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$treatment <- df$treatment  # Ensure this matches your treatment column name

# Add a column for transect numbers (assuming sequential from 0 to n-1)
nmds_scores$transnum <- 0:(nrow(nmds_scores) - 1)

# Plot NMDS with ggplot2, colored by treatment
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = treatment)) +
  geom_point(size = 3) +
  geom_text(aes(label = transnum), vjust = -1, size = 3) +  # Add transect numbers
  scale_color_manual(values = c("blue", "darkorange")) +  # Customize colors
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove all gridlines
    panel.border = element_rect(color = "black", fill = NA),  # Add a clean border
    legend.position = "right"
  ) +
  labs(title = "NMDS Plot Colored by Treatment",
       x = "NMDS1",
       y = "NMDS2",
       color = "treatment")

#  scale_color_manual(values = c("darkorange", "forestgreen")) +  # Customize treatment colors


# 4. Abundance Analysis
# If you want to check if overall species abundance differs between sites:
  
# Sum total abundance for each transect
df$Total_Abundance <- rowSums(df[, 7:ncol(df)])
# df$Total_Abundance <- rowSums(df[, 4:ncol(df)]) # for species means dataframe
# Compare abundance between sites
t_test_abundance <- t.test(Total_Abundance ~ treatment, data = df)
print(t_test_abundance) # Results significant: t = 8.4208, df = 277.95, p-value = 2.018e-15,
# mean in group far ; mean in group near 
# 78.01411             39.29537 

## Result for species means - pooled by transect: t = 3.4997, df = 25.672, p-value = 0.001719
# mean in group far  mean in group near 
# 95 percent confidence interval: 31.76429 to 122.31993
#     157.12813           80.08602 
## Doesn't make sense - this is computing a measure that is not relevant - mean is not 78% so what does that stand for? I think
# it is total abundance, see line 392 - so is the sum of each row - the sume of abundance
# in each plot is about 78% cover in far (undisturbed) and 39% in near (trampled), since
# this includes only live above-ground species cover, not ground cover (rock, soil, litter)

# Explanation of Steps:
# Alpha Diversity: Measures diversity within each transect, calculated using Shannon and Simpson indices.
#	Beta Diversity: Measures how species composition differs between transects/sites, analyzed with Bray-Curtis and visualized via NMDS; each dot is a plot
#	PERMANOVA: Used to test whether there are statistically significant differences in community composition between the sites.
#	Abundance: A simple comparison of total species counts between the sites.
# You can extend this analysis by adding more specific metrics or visualizations as required by your research.



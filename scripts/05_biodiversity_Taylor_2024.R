#########################
# Data on cover of plant species to evaluate trampling effects of near and far from trail
# Goal: to wrangle diversity info; compute species composition NMDS and visualizations; conduct analyeses (PERMANOVA)
# Script: Sept 2024 Nina Hewitt; 
# Sources: BC Parks Github (Cassandra Elphinstone) - for the Shannon/Simpson comps
# Edited: Oct 1 2024_Go to line 86 end for Oct 17 Script;
# Updated Mar 13 2025 for ggplot versions of beta diversity NMDS plots, and 
# added plot with points colored by site
# Updated Feb 22, 23, 2026
# - Added sqrt transformation for percent cov values; collapsed bryophytes (single variable/columnn)
# - Ran with crust (used this, though little change); and crust remove
# - Re-run permanova, nmds and plot nmds in ggplot with hulls 
# - Combined site and treatment into single figure and improved color scheme
# - Finally re-plotted figures with centroids (sites) and ellipses (treatments)
#    representing the 68% standard error region around treatment centroids.
#########################



## A. WRANGLE DATAFRAME TO EXTRACT TRANSECT MEAN VALUES
# To obtain and save pooled data, proceed
# To use plot level data skip and go to line 87 (for Shannon/Simpson's diversity);
# or straight to lines 160-__ for the Bray-Curtisdissimiliarity and NMDS, Permanova]

# clear workspace
rm(list = ls())

## LOAD PACKAGES
library(ggplot2)
library(dplyr)
library(tidyverse)
library(data.table) # faster data.frame

# BRING IN THE DATA

## read the csv data file (ensure you are in the correct working directory)
df <-  fread(file.path("data","species_matrix_factorcolumns.csv")) # this command works if you have loaded the data.table library/package
## otherwise load dataset with file.choose command
# df <- read.csv(file.choose()) # choose "species_matrix_factorcolumns.csv" from my [updated trampling matrixes] or RScripts working [Sept. 2025 version] directory
## otherwise load dataset with file.choose command

# For transect level means only, choose "species_means.csv" (from updated trampling matrixes in rscripts folder)
# for just mosses in Connor's transects: "species_matrix_mosses_CGW_transects_nocrust.csv"
# for just vascular: "species_matrix_factorcolumns_forbs_grasses_shrubs_trees.csv"
# for just vascular, no graminoids: "species_matrix_factorcolumns_forbs_shrubs_trees.csv"

## LOOK AT OUR DATA
df          # look at all your data
head(df)    # look at the first few entries
tail(df)    # look at the last few entries
str(df)     # look at the overall structure

# Convert all columns to numeric except 1st 6 [or 8] - which are as factor 
# Load necessary library
library(dplyr)

df <- df %>%
  mutate(across(1:6, as.factor)) %>%
  mutate(across(7:ncol(df), as.numeric))

# Display the structure of the modified data frame
str(df)

## If dataframe requires it, replace NA with zeros for numeric columns (already acheived)
# df <- mydata %>%
#  mutate(across(where(is.numeric), ~replace_na(., 0)))
## or
#mydata[numeric_vars] <- lapply(mydata[numeric_vars], function(x) ifelse(is.na(x), 0, x))

# Display the modified data frame
head(df)

# skip to Bray Curtis, NMDS (we employed the zero-inflated models in file "01Diversity.R" for diversity analyses)
# Compute mean values for each transect
mean_values <- df %>%
  group_by(transnum) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Write csv
write.csv(mean_values, file = "species_means.csv")
# write.csv(df, file = "name.csv")

# Display the modified data frame
head(df)


######### COMPUTE SHANNON AND SIMPSON'S DIVERSITY INDICES
## We use a zero-inflated Bayes Model in the manuscript

# Clears the workspace
rm(list=ls())

# Ask R where it is right now
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(data.table) # faster data.frame

# BRING IN THE DATA

## read the csv data file (ensure you are in the correct working directory)
df <-  fread(file.path("data","species_matrix_factorcolumns.csv")) # this command works if you have loaded the data.table library/package

## B. Alpha Diversity (Shannon, Simpson)
# Alpha Diversity: Measures diversity within each transect, calculated using Shannon and Simpson indices
## calculate common diversity indices for each transect and compare the results.
 
# Calculate Shannon and Simpson diversity indices
df$Shannon <- diversity(df[, 7:ncol(df)], index = "shannon")
df$Simpson <- diversity(df[, 7:ncol(df)], index = "simpson")
# Or for species means dataframe, adjust accordingly:

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

# When species' transect means are used, results NOT significant 
## * However, a zero-inflated model was more appropriate to this analysis and was used instead (see 01Diversity.R Script)

## If using this dataframe for next steps, REMOVE the last 2 columns from the dataframe (simpsons, shannon) as these 
# # are derived values from previous step - use either:
# species_dat <- species_dat[, !colnames(species_dat) %in% c("Shannon", "Simpson")] # or
# species_dat <- species_dat[, !tolower(colnames(species_dat)) %in% c("shannon", "simps

#######
## C. COMPUTE BRAY-CURTIS DISSIMILARITY AND CONDUCT NMDS ORDINATION, PERMANOVA ANAL
## UPDATED: March 27 2026

## Use dissimilarity measures and ordination (NMDS) to assess and visualize differences 
## in species composition.

# -------------------------------
# 0. Load packages
# -------------------------------
rm(list=ls()) # clears workspace

library(data.table)
library(dplyr)
library(ggplot2)
library(vegan)
library(permute)

# -------------------------------
# 1. Load data
# -------------------------------
df <- fread(file.path("data","species_matrix_factorcolumns.csv")) ## read the csv data file (ensure in the correct wd and data.table library loaded
## otherwise load dataset with file.choose command
# df <- read.csv(file.choose()) # choose "species_matrix_factorcolumns.csv" (from Sept 13, 2025 version from my RScripts working directory [updated trampling matrixes])
# or use the command read.csv() if you are in the correct working directory
# df = read.csv("___.csv")

# Convert columns
df <- df %>%
  mutate(across(c(1:4,6), as.factor)) %>%   # keep elevation numeric (col 5)
  mutate(across(7:ncol(df), as.numeric))

# -------------------------------
# 2. Create plot labels + remove outliers
# -------------------------------
df$plot_label <- paste0(df$transect, "_", df$quad) #  Create plot labels FIRST (needed for filtering)
outliers <- c("PR-17_3", "BT-15_3") # Define outliers to remove
df <- df[!(df$plot_label %in% outliers), ] # Remove outlier plots from dataframe

# -------------------------------
# 3. Split metadata and species
# -------------------------------
meta <- df[, 1:6]   # metadata
species_dat <- df[, 7:ncol(df)]
# Force to a plain data.frame right after creating species_dat:

species_dat <- as.data.frame(df[, 7:ncol(df)])

# REMOVE plot_label if it slipped in
species_dat <- species_dat[, !colnames(species_dat) %in% "plot_label"]

## species_dat should now be 63 variables for 59 native vasc, bryo taxa plus bryosp1, bryosp2; crust; agrcap (invasive grass left in - it indicates compositional diffs so leave in)

# -------------------------------
# 4. Define cryptogam species (just bryophytes)
# -------------------------------
# Place all mosses in a category, since Bray–Curtis and NMDS assume:
# - Each column represents comparable biological units and differences reflect ecology, not data structure
# - Current dataset mixes species-level bryophytes (fine resolution) with “bryosp1" "2” # which artificially inflates dissimilarity.
# - "crust" is another category (coarse resolution) but we will leave as separate

cryptogam_spp <- c(
  "aulpal","braspp","bryosp1","bryosp2","bucsud","cerpur",
  "dicoly","dicpal","ditspp","meilya","nipelo","nipmut",
  "nippyg","phifon","pohnut","polpil","polspp","rhysqu",
  "rhyrob","scioed","strcon","barhat","lopven"
)

# NB: this includes only bryophytes even though is named cryptogam_spp

# -------------------------------
# 5. Build community matrix
# -------------------------------

# Start with species data
comm <- species_dat

# Combine cryptogams (+ crust if desired - I did not combine)
comm$cryptogams <- rowSums(
  comm[, colnames(comm) %in% cryptogam_spp],
  na.rm = TRUE
)

# OPTIONAL: if crust does not overlap with moss species can combine with moss:
# comm$cryptogams <- rowSums(
#   comm[, colnames(comm) %in% cryptogam_spp],
#   na.rm = TRUE
# ) + comm$crust
## I simply left crust in the dataframe as a sep column along with crypto (only bryophytes)

# OPTIONAL: include crust in cryptogams
# comm$cryptogams <- comm$cryptogams + comm$crust

# Remove original bryophyte species columns
comm <- comm[, !colnames(comm) %in% cryptogam_spp] 


# -------------------------------
# 6. Remove empty rows (plots with all zeros)
# -------------------------------
keep_rows <- rowSums(comm, na.rm = TRUE) > 0

comm <- comm[keep_rows, ]
meta <- meta[keep_rows, ] # filter metadata to match

# -------------------------------
# 7. Square-root Transform 
# -------------------------------
comm_sqrt <- sqrt(comm)

# -------------------------------
# 8. Compute Bray-Curtis
# -------------------------------
bc <- vegdist(comm_sqrt, method = "bray")

# -------------------------------
# 9. NMDS
# -------------------------------
nmds <- metaMDS(bc, k = 2, trymax = 200)
nmds$stress

## try a stress plot of the results of the metaMDS function:

stressplot(nmds) # looks good for both plot level data and transect means
## if there is "Large scatter around the line suggests that original 
# dissimilarities are not well preserved in the reduced number of dimensions." 
# (https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/) In this case it looks good 

# Run 20 stress 0.2318025 
# *** Best solution repeated 1 times
# > nmds$stress
# [1] 0.2123354
# > stressplot(nmds)

# -------------------------------
# 10. Extract scores + metadata
# -------------------------------
nmds_scores <- as.data.frame(scores(nmds))

# Use filtered metadata
nmds_scores$treatment <- meta$treatment
nmds_scores$site <- meta$site

## [for ellipses (treatment) with centroids (sites) go to lines 410 below]

# -------------------------------
# 11. Create convex hulls (treatment)
# -------------------------------
get_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

hulls <- nmds_scores %>%
  group_by(treatment) %>%
  group_modify(~ get_hull(.x))

# -------------------------------
# 12. NMDS plot (treatment)
# -------------------------------
ggplot(nmds_scores, aes(NMDS1, NMDS2, color = treatment, fill = treatment)) +
  geom_point(size = 3) +
  geom_polygon(data = hulls, aes(fill = treatment),
               alpha = 0.2, color = "black") +
  scale_color_manual(values = c("blue", "darkorange")) +
  scale_fill_manual(values = c("blue", "darkorange")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(
    title = "NMDS plot of beta diversity by treatment",
    x = "NMDS1",
    y = "NMDS2"
  )

# -------------------------------
# 13. PERMANOVA (CORRECT DATA SOURCE)
# -------------------------------
adonis2(bc ~ treatment * site, data = meta, permutations = 999)

# Results, Mar 27, 2026:
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bc ~ treatment * site, data = meta, permutations = 999)
# Df SumOfSqs     R2      F Pr(>F)    
# Model      5   15.769 0.2332 16.361  0.001 ***
#   Residual 269   51.852 0.7668                  
# Total    274   67.621 1.0000     

# -------------------------------
# 14. Dispersion tests
# -------------------------------
disp_treat <- betadisper(bc, meta$treatment)
anova(disp_treat)
boxplot(disp_treat)

disp_site <- betadisper(bc, meta$site)
anova(disp_site)
boxplot(disp_site)

# -------------------------------
# 15. Site labels for plotting
# -------------------------------
site_labels <- c("1" = "TM", "2" = "BT", "3" = "PR")
nmds_scores$site_name <- site_labels[as.character(nmds_scores$site)]

## Back to my original code chunks:

# -------------------------------
# Create convex hulls for site
# -------------------------------
get_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
hulls_site <- nmds_scores %>%
  group_by(site_name) %>%
  group_modify(~ get_hull(.x))

# -------------------------------
# Plot NMDS by site
# -------------------------------
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = site_name, fill = site_name)) +
  geom_point(size = 3) +
  geom_polygon(data = hulls_site, aes(x = NMDS1, y = NMDS2, fill = site_name),
               alpha = 0.2, color = "black") +
  scale_color_manual(values = c("TM" = "forestgreen", "BT" = "purple", "PR" = "goldenrod")) +
  scale_fill_manual(values = c("TM" = "forestgreen", "BT" = "purple", "PR" = "goldenrod")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  ) +
  labs(
    title = "NMDS Plot of Vascular plant and Bryophyte Beta Diversity by Site",
    x = "NMDS1",
    y = "NMDS2",
    color = "Site",
    fill = "Site"
  )


## combined (treatment and site)

library(ggplot2)
library(dplyr)

# -------------------------------
# Map site labels
# -------------------------------
site_labels <- c("1" = "TM", "2" = "BT", "3" = "PR")
nmds_scores$site_name <- site_labels[as.character(nmds_scores$site)]

# -------------------------------
# Create convex hulls for site
# -------------------------------
get_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
hulls_site <- nmds_scores %>%
  group_by(site_name) %>%
  group_modify(~ get_hull(.x))

# -------------------------------
# Plot NMDS with site hulls and treatment points
# -------------------------------
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  # Hulls for sites
  geom_polygon(data = hulls_site, aes(x = NMDS1, y = NMDS2, fill = site_name),
               alpha = 0.2, color = "black") +
  # Points for treatments
  geom_point(aes(color = treatment, shape = treatment), size = 3) +
  scale_fill_manual(values = c("TM" = "forestgreen", "BT" = "purple", "PR" = "goldenrod")) +
  scale_color_manual(values = c("trampled" = "blue", "reference" = "darkorange")) +
  scale_shape_manual(values = c("trampled" = 19, "reference" = 17)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  ) +
  labs(
    title = "NMDS of Beta Diversity by Site and Trample Condition",
    x = "NMDS1",
    y = "NMDS2",
    fill = "Site",
    color = "Condition",
    shape = "Condition"
  )

## better for color blind: Okabe & Ito palette (or modified),
## Note that PR-17, plot 3 and BT-15, plot 3 are both outliers

library(ggplot2)
library(dplyr)

# Map site labels
site_labels <- c("1" = "TM", "2" = "BT", "3" = "PR")
nmds_scores$site_name <- site_labels[as.character(nmds_scores$site)]

# Convex hulls by site
get_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
hulls_site <- nmds_scores %>%
  group_by(site_name) %>%
  group_modify(~ get_hull(.x))

# Colorblind-friendly palette for 3 sites
site_palette <- c(
  "TM" = "#E69F00",  # orange
  "BT" = "#009E73",  # green
  "PR" = "#0072B2"   # blue
)

# Shape palette for treatment
treat_shapes <- c("trampled" = 19, "reference" = 17)
treat_colors <- c("trampled" = "#56B4E9", "reference" = "#D55E00") # optional, different from site

#### nice fig, but use below with centroids

# Plot NMDS
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  # Hulls for sites
  geom_polygon(data = hulls_site,
               aes(x = NMDS1, y = NMDS2, fill = site_name),
               alpha = 0.2, color = "black") +
  # Points colored/shaped by treatment
  geom_point(aes(color = treatment, shape = treatment), size = 3) +
  # Fill for hulls
  scale_fill_manual(values = site_palette) +
  # Points
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = treat_shapes) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  ) +
  labs(
    title = "NMDS of vascular plant and bryophyte Beta Diversity",
    x = "NMDS1",
    y = "NMDS2",
    fill = "Site",
    color = "Condition",
    shape = "Condition"
  )


## add plot numbers - create a new df called nmds_scores with new column 

library(ggplot2)
library(dplyr)
# install.packages("ggrepel")
library(ggrepel)  # for non-overlapping labels

# -------------------------------
# Create plot labels
# -------------------------------
nmds_scores$plot_label <- paste0(meta$transect, "_", meta$quad)


# -------------------------------
# NMDS plot: hulls by site, points by treatment, labels by plot
# -------------------------------

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  # Hulls for sites
  geom_polygon(data = hulls_site,
               aes(x = NMDS1, y = NMDS2, fill = site_name),
               alpha = 0.2, color = "black") +
  # Points colored/shaped by treatment
  geom_point(aes(color = treatment, shape = treatment), size = 3) +
  # Labels for each point (transect + quad)
  geom_text_repel(aes(label = plot_label),
                  size = 3,
                  max.overlaps = 20) +  # adjust max-overlaps as needed
  # Color palettes
  scale_fill_manual(values = site_palette) +
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = treat_shapes) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  ) +
  labs(
    title = "NMDS of Beta Diversity with Plot Labels",
    x = "NMDS1",
    y = "NMDS2",
    fill = "Site",
    color = "Condition",
    shape = "Condition"
  )


## Centroids instead - USE THIS

# -------------------------------
# Centroids by Site
# -------------------------------
centroids_site <- nmds_scores %>%
  group_by(site_name) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2)
  )


# -------------------------------
# Centroids by Treatment
# -------------------------------
centroids_treat <- nmds_scores %>%
  group_by(treatment) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2)
  )


# NMDS plot with centroids (no hulls)
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  
  # Points
  geom_point(aes(color = treatment, shape = treatment), size = 3) +
  
  # Site centroids
  geom_point(data = centroids_site,
             aes(x = NMDS1, y = NMDS2, fill = site_name),
             shape = 21,          # filled circle with border
             color = "black",
             size = 5,
             stroke = 1.2) +
  
  # Labels for centroids
  geom_text(data = centroids_site,
            aes(label = site_name),
            vjust = -1,
            fontface = "bold") +
  
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = treat_shapes) +
  scale_fill_manual(values = site_palette) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  ) +
  labs(
    title = "NMDS of Vascular Plant and Bryophyte Beta Diversity",
    x = "NMDS1",
    y = "NMDS2",
    color = "Condition",
    shape = "Condition",
    fill = "Site"
  )

## Add dispersion structure using SE ellipses
# stat_ellipse(aes(color = treatment), level = 0.68)

(plot_nmds_export <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  
  # SE ellipses (draw FIRST so points sit on top)
  stat_ellipse(aes(color = treatment),
               level = 0.68,      # ~ multivariate SE analogue
               linewidth = 1,
               linetype = 1) +
  
  # Points
  geom_point(aes(color = treatment, shape = treatment), size = 3) +
  
  # Site centroids
  geom_point(data = centroids_site,
             aes(x = NMDS1, y = NMDS2, fill = site_name),
             shape = 21,
             color = "black",
             size = 5,
             stroke = 1.2) +
  
  # Centroid labels
  geom_text(data = centroids_site,
            aes(label = site_name),
            vjust = -1,
            fontface = "bold") +
  
  scale_color_manual(values =  c("#27A81E", "#DEBF50")) +
  scale_shape_manual(values = treat_shapes) +
  scale_fill_manual(values = site_palette) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  ) +
  labs(
    #title = "NMDS of vascular plant, bryophyte and crust beta diversity",
    x = "NMDS1",
    y = "NMDS2",
    color = "Condition",
    shape = "Condition",
    fill = "Site"
  ))



ggsave(file.path("figure","Fig_NMDS_colors.jpg"),plot_nmds_export,
       unit = "mm", width = 180,height = 150,dpi = 250, scale =0.9)  

#### END

### Old steps, 
# 4. Abundance Analysis
# to check if overall species abundance differs between sites:

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


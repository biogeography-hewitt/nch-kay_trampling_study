## Description: this script creates stacked bar charts of cover of 
## lifeform and ground cover classes by treatment and site, and analyse differences
## Created: Mar 7, 2025, Nina Hewitt
## Most recent update: March 24, Nina Hewitt to recreate as a project
## in RProjects_2026 folder (in Documents)
## Later will upload to Github Repo for this project and publication

##
# Note - the raw data in google has species that were not encountered in plots
# incl. AntRos, AntUmb, EryGra, PenPro, RanEsc, FesBra, VacMem, SalCom, SorSit

# Clears the workspace
rm(list=ls())

# Ask R where it is right now
getwd()

# BRING IN THE DATA

##### packages ####
library(data.table) # faster data.frame
library(ggplot2) # plots

## read the csv data file (ensure you are in the correct working directory
df <-  fread(file.path("data","species_means_bylifeform.csv")) # this command works if you have loaded the data.table library/package

## otherwise load dataset with file.choose command
# df <- read.csv(file.choose()) # I chose "species_means_bylifeform_summed.csv" Feb 5, 2025, from my RScripts working directory
# then use "species_means_bylifeform_summed_b.csv" (slightly diff transect labels in "_c..")

## LOOK AT OUR DATA
df          # look at all your data
head(df)    # look at the first few entries
tail(df)    # look at the last few entries
str(df)     # look at the overall structure

# Convert all columns to numeric except the first 5 which are as factor 
# Load necessary library [check dataframe and adjust accordingly]
library(dplyr) 

df <- df %>%
  mutate(across(1:5, as.factor)) %>%
  mutate(across(6:ncol(df), as.numeric))

## Display the structure of the modified data frame
str(df)

## Visualizations
# load libraries
## to run analyses, proceed to lines 165 onward

library(ggplot2)
library(dplyr)
library(tidyr)

# Reshape data from wide to long format
df_long <- df %>%
  pivot_longer(cols = -c(transnum, treatment, site, transectid, siteno),  # Keep factors, convert lifeform columns
               names_to = "Lifeform",
               values_to = "Cover") %>%
  filter(!is.na(Cover) & Cover > 0)  # Remove NAs and zeros if needed

# Summarize total cover by Lifeform, Treatment, and Site
df_summary <- df_long %>%
  group_by(treatment, site, Lifeform) %>%
  summarise(Total_Cover = sum(Cover, na.rm = TRUE), .groups = "drop")

# Create Stacked Bar Chart with white background
######  Custom colors scheme for color impaired 
## - USED THIS VERSION, but with values adjusted to 100%, as in lines 212-261

library(ggplot2)
library(dplyr)
library(tidyr)

# Reorder Site factor for logical flow
df_summary$site <- factor(df_summary$Site, levels = c("Taylor meadows", "Black tusk", "Panorama ridge"))

lifeform_colors <- c(
  "bryophyte" = "#009E73",   # Medium green (colorblind-friendly)
  "conifer" = "#E69F00",     # Deep orange
  "dshrub" = "#56B4E9",      # Sky blue (light contrast)
  "eshrub" = "#CC79A7",      # Soft magenta-pink (differentiates from red)
  "forbs" = "#F0E442",       # Bright yellow (high visibility)
  "graminoid" = "#0072B2",   # Darker blue (contrasts with forbs)
  "litter" = "#999999",      # Medium gray (neutral)
  "rockcrust" = "#D55E00",   # Rust orange (good contrast with blue)
  "soil" = "#660066"         # Deep purple (contrasts with all other colors)
)

ggplot(df_summary, aes(x = site, y = Total_Cover, fill = Lifeform)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black") +  # Add black outlines
  facet_wrap(~treatment) +
  scale_fill_manual(values = lifeform_colors) +  # Apply custom colors
  labs(title = "Lifeform Composition by Site and Treatment",
       x = "Site",
       y = "Total cover",
       fill = "Functional group") +
  theme_minimal(base_size = 14) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))

## Redo using proportion cover (out of 100 by site) and
## Custom colors for visually color impaired

library(tidyr)
library(dplyr)
library(data.table) # faster data.frame


## reading the csv data file if you are in the correct working directory
## (RProjects_2026)
df_summ_2 <-  fread(file.path("data","fun_group_summed_prop_cov.csv")) # this command works if you have loaded the data.table library/package

## otherwise load dataset with file.choose command:
# df_summ_2 <- read.csv(file.choose()) # I chose "fun_group_summed_prop_cov.csv" 
## working excel file is "fun_group_adjust_tot_cov.xlsx" in my RScripts working directory

## make sure vars 1-3 are as factor 4- are as numeric
head(df_summ_2)
tail(df_summ_2)

df_summ_2 <- df_summ_2 %>%
  mutate(across(1:3, as.factor)) %>%
  mutate(across(7:ncol(df), as.numeric))

# Display the structure of the modified data frame
str(df_summ_2)


library(ggplot2)
library(dplyr)
library(tidyr)

# Reorder Site factor for logical flow
df_summ_2$site <- factor(df_summary_2$site, levels = c("Taylor meadows", "Black tusk", "Panorama ridge"))

lifeform_colors <- c(
  "bryophyte" = "#009E73",   # Medium green (colorblind-friendly)
  "conifer" = "#E69F00",     # Deep orange
  "dshrub" = "#56B4E9",      # Sky blue (light contrast)
  "eshrub" = "#CC79A7",      # Soft magenta-pink (differentiates from red)
  "forbs" = "#F0E442",       # Bright yellow (high visibility)
  "graminoid" = "#0072B2",   # Darker blue (contrasts with forbs)
  "litter" = "#999999",      # Medium gray (neutral)
  "rockcrust" = "#D55E00",   # Rust orange (good contrast with blue)
  "soil" = "#660066"         # Deep purple (contrasts with all other colors)
)

ggplot(df_summ_2, aes(x = site, y = Prop_Cover, fill = lifeform)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black") +  # Add black outlines
  facet_wrap(~treatment) +
  scale_fill_manual(values = lifeform_colors) +  # Apply custom colors
  labs(title = "Lifeform Composition by Site and Treatment",
       x = "Site",
       y = "Cover (percent of site total)",
       fill = "Functional group") +
  theme_minimal(base_size = 14) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))


## Analyses of cover by lifeform and treatment

# Use a linear mixed-effects model
# First take dataframe (same as uploaded in row 21) called df
# and make it long so that each row is one transect and one lifeform
# call is df_long2, since df_long, above, has 0s removed (is only 180 obs)
# while this includes 0s that will go into analysis (is 252 obs).

library(tidyr)
library(dplyr)

df_long2 <- df %>%
  pivot_longer(
    cols = forb:rockcrust,
    names_to = "lifeform",
    values_to = "cover"
  )

head(df_long2)
tail(df_long2)

## Fit mixed model.
## Model tests: "Is cover significantly different by lifeform and treatment?"

# Core model:

library(lme4)
# install.packages("lmerTest")
help(lmerTest)
library(lmerTest)

m1 <- lmer(
  cover ~ lifeform * treatment + site +
    (1 | transnum),
  data = df_long2
)
anova(m1)

# RESULTS
# > anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                     Sum Sq Mean Sq NumDF DenDF F value Pr(>F)    
# lifeform           28851.4  3606.4     8   232 27.0563 <2e-16 ***
# treatment              8.4     8.4     1   232  0.0632 0.8018    
# site                  22.4    11.2     2   232  0.0840 0.9195    
# lifeform:treatment 23990.9  2998.9     8   232 22.4982 <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# So lifeforms are sig different and the interaction between lifeform and treatment are sig

m1 # Linear mixed model fit by REML ['lmerModLmerTest']
# Formula: cover ~ lifeform * treatment + site + (1 | transnum)


## Next: Follow-up (if interaction is significant)

#install.packages("emmeans")
library(emmeans)
emmeans(m1, pairwise ~ treatment | lifeform)

# These results indicate that “Trampling significantly reduced eshrub cover;
# although forb cover was greater in undisturbed, differences were not 
# significant, likely due to Black Tusk; Graminoid cover was not 
# significantly different.
# Bare soil was greater and litter significantly lower in trampled.

# or try [doesn't make sense]:
(lsm <- ls_means(m1))
ls_means(m1, which = "treatment", pairwise = TRUE)

## ls_means also have plot and as.data.frame methods:
## Not run: 
plot(lsm, which=c("treatment", "site"))
as.data.frame(lsm)
## Inspect the LS-means contrasts:
show_tests(lsm, fractions=TRUE)$treatment


## Combine with Permanova

library(vegan)

cover_mat <- df %>%
  select(forb:rockcrust)

adonis2(
  cover_mat ~ site * treatment,
  data = df,
  method = "bray"
)

# Permutation test for adonis under reduced model
## ***USE THIS PERMANOVA *** - it deals with zeros by using transect averages for each functional group

### NB: "This is a permutation method used in multivariate analysis 
### to test the significance of ecological community structure
## evaluates factors by shuffling data while respecting constraints 
## (strata), testing for significance against the reduced, null model 
## rather than the full model....does not rely on 
## traditional parametric assumptions like multivariate normality,
## making it suitable for complex community datasets
## ...  It is widely used to analyze the effect of treatments,
## environmental variables, or temporal factors on microbial or 
## ecological community structures." Github ; https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/

## ** Results ** use this

# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = cover_mat ~ site * treatment, data = df, method = "bray")
# Df SumOfSqs     R2           F     Pr(>F)    
# Model     5   2.6135 0.5344 5.0503  0.001 ***
# Residual 22   2.2770 0.4656                  
# Total    27   4.8905 1.0000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## REDO with log values
# library(lme4)
# install.packages("lmerTest")
help(lmerTest)
library(lmerTest)

# create a logged variable, add 1 due to zeros

df_long2$log_cover <- log(df_long2$cover + 1)

m2 <- lmer(
  log_cover ~ lifeform * treatment + site +
    (1 | transnum),
  data = df_long2
)
anova(m2)

# > anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                     Sum Sq Mean Sq NumDF DenDF F value Pr(>F)    
# lifeform           180.660 22.5824     8   232 25.2600 <2e-16 ***
# treatment            1.863  1.8629     1   232  2.0837 0.1502    
# site                 0.026  0.0130     2   232  0.0146 0.9855    
# lifeform:treatment 107.816 13.4770     8   232 15.0749 <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Shows that lifeform and lifeform*treatment are sig

m2 # Linear mixed model fit by REML ['lmerModLmerTest']
# Formula: logcover ~ lifeform * treatment + site + (1 | transnum)


## Next: Follow-up (if interaction is significant)

# install.packages("emmeans")
library(emmeans)
emmeans(m2, pairwise ~ treatment | lifeform)

# These results indicate that “Trampling significantly reduced bryophyte and eshrub cover;
# although cover for all lifeforms except graminoids was greater in undisturbed, transect mean cover differences were 
# significant only for bryophytes and forbs, likely due to Black Tusk; Graminoid cover was not 
# significantly different.
# Bare soil was greater and litter significantly lower in trampled.


## Combine with Permanova

library(vegan)

cover_mat <- df %>%
  select(forb:rockcrust)

adonis2(
  cover_mat ~ site * treatment,
  data = df,
  method = "bray"
)

# Permutation test for adonis under reduced model

# 
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = cover_mat ~ site * treatment, data = df, method = "bray")
# Df SumOfSqs     R2      F Pr(>F)    
# Model     5   2.6135 0.5344 5.0503  0.001 ***
#   Residual 22   2.2770 0.4656                  
# Total    27   4.8905 1.0000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### *** Use this analysis - see explanation in Word doc (Nina's)

## REDO with plot level (rather than summed transect level means) for each group
## upload new file "lifeform_matrix_plotlevel.csv" in Rscripts, in folder "Updated Trampling Matrixes_Nov2024_2026"

## START OVER with plot level dataframe of lifeform/functional types:
# Clears the workspace
rm(list=ls())

## read the csv data file (ensure you are in the correct working directory
df_plot <-  fread(file.path("data","lifeform_matrix_plotlevel.csv")) # this command works if you have loaded the data.table library/package
## will need to fix fact that each lifeform column is numbered consecutively
# see lines 405 onward.


## otherwise load dataset with file.choose command
# df_plot <- read.csv(file.choose()) # new file "lifeform_matrix_plotlevel.csv" (local file in folder "Updated Trampling Matrixes_Nov2024_2026")


## LOOK AT OUR DATA
df_plot          # look at all your data
head(df_plot)    # look at the first few entries
tail(df_plot)    # look at the last few entries
str(df_plot)     # look at the overall structure

# Convert all columns to numeric except 1st 6 - which are as factor 
# Load necessary library
library(dplyr)

# df_plot <- df_plot %>% e.g., [check how many columns are factor and adjust accordingly]

df_plot <- df_plot %>%
  mutate(across(1:8, as.factor)) %>%
  mutate(across(9:ncol(df_plot), as.numeric))

# Display the structure of the modified data frame
str(df_plot)
## need to fix fact that each lifeform column is numbered consecutively per l. 405 onward.


## Analyses of cover by lifeform and treatment

# Use a linear mixed-effects model

## Model 1: vegetation lifeforms only

# First take dataframe (same as uploaded in row 21) called df
# and make it long so that each row is one transect and one lifeform
# call is df_long2, since df_long, above, has 0s removed (is only 180 obs)
# while this includes 0s that will go into analysis (is 252 obs).
library(tidyverse)

df_long_plot <- df_plot %>%
  pivot_longer(
    cols = starts_with(c("forb", "graminoid", "bryophyte", "eshrub", "dshrub", "conifer", "litter", "soil", "rockcrust" )),
    names_to = "lifeform",
    values_to = "cover"
  ) %>%
  mutate(
    lifeform = case_when(
      str_detect(lifeform, "forb") ~ "forb",
      str_detect(lifeform, "graminoid") ~ "graminoid",
      str_detect(lifeform, "bryophyte") ~ "bryophyte",
      str_detect(lifeform, "eshrub") ~ "eshrub",
      str_detect(lifeform, "dshrub") ~ "dshrub",
      str_detect(lifeform, "conifer") ~ "conifer",
      str_detect(lifeform, "litter") ~ "litter",
      str_detect(lifeform, "soil") ~ "soil",
      str_detect(lifeform, "rockcrust") ~ "rockcrust",
    )
  )

## filter just lifeforms - vegetation [for substrate jump to line 558]
df_veg <- df_long_plot %>%
  filter(lifeform %in% c("forb", "graminoid", "bryophyte", "eshrub", "dshrub", "conifer"))

## Fit mixed model.
# A: MINIMUM MODEL
## Model tests: "Is cover significantly different by lifeform and treatment?"

## load library
library(lme4)
# install.packages("lmerTest")
# help(lmerTest)
library(lmerTest)

##log values
m_veg <- lmer(
  log(cover + 1) ~ lifeform * treatment + site +
    (1 | transnum),
  data = df_veg
)
anova(m_veg)

## Results - *** use this ***

# > anova(m_veg)
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# lifeform           650.48 130.097     5 17042.0 407.536 < 2.2e-16 ***
# treatment            7.56   7.559     1    36.4  23.680 2.206e-05 ***
# site                 1.23   0.616     2    24.0   1.929    0.1672    
# lifeform:treatment  74.00  14.800     5 17042.0  46.363 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > 

# So lifeform and treatment are sig different; and the interaction between lifeform and treatment is sig

## run pairwise comparisons on substrate
# install.packages("emmeans")
# help(emmeans)
library(emmeans)

emmeans(m_veg, pairwise ~ treatment | lifeform)

# see results in ms. 
## Ignore message the D.f. calcs disabled due to > 3000 observations - the model
## still calculates test stats based on asymptotic (large-sample) theory; 
## estimates and contrasts are unchanged.


### SUBSTRATE USE THIS - ground values linear mixed effect model 

library(lme4)
# install.packages("lmerTest")
# help(lmerTest)
library(lmerTest)

df_sub <- df_long_plot %>%
  filter(lifeform %in% c("litter", "soil", "rockcrust"))

m_sub <- lmer(
  log(cover + 1) ~ lifeform * treatment + site +
    (1 | transnum),
  data = df_sub
)

anova(m_sub)

# Substrate ***Use this*** - note that lifeform is just the label for the column
# in reporting it will be labelled substrate class

## anova(m_sub)
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
# lifeform[substrate] 182.64   91.32     2   808  73.7536 < 2.2e-16 ***
# treatment           11.81   11.81     1    24   9.5365  0.005029 ** 
# site                 3.67    1.84     2    24   1.4823  0.247172    
# lifeform:treatment 928.38  464.19     2   808 374.8938 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > 

## run pairwise comparisons on substrate
# install.packages("emmeans")
# help(emmeans)
library(emmeans)

emmeans(m_sub, pairwise ~ treatment | lifeform)


# $contrasts
# lifeform = litter:
#   contrast   estimate    SE   df t.ratio p.value
# far - near   1.7148 0.214 43.4   8.010 <0.0001
# 
# lifeform = rockcrust:
#   contrast   estimate    SE   df t.ratio p.value
# far - near  -0.0638 0.214 43.4  -0.298  0.7671
# 
# lifeform = soil:
#   contrast   estimate    SE   df t.ratio p.value
# far - near  -3.3603 0.214 43.4 -15.696 <0.0001
# 
# Results are averaged over the levels of: site 
# Note: contrasts are still on the log(mu + 1) scale. Consider using
# regrid() if you want contrasts of back-transformed estimates. 
# Degrees-of-freedom method: kenward-roger 

# Results: Cover values for litter were sig greater in undisturbed, and bare soil in trampled; rock cover did not differ significantly


## B: 2nd MODEL option (if plots are uniquely identified) ## Try with the advanced model

## log values:

## load library
library(lme4)
# install.packages("lmerTest")
# help(lmerTest)
library(lmerTest)

##log values
m_veg2 <- lmer(
  log(cover + 1) ~ lifeform * treatment + site +
    (1 | transnum/quad),
  data = df_veg
)
anova(m_veg2)

# > anova(m_veg2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# lifeform           650.48 130.097     5 17042.0 407.536 < 2.2e-16 ***
# treatment            7.56   7.559     1    36.4  23.680 2.206e-05 ***
# site                 1.23   0.616     2    24.0   1.929    0.1672    
# lifeform:treatment  74.00  14.800     5 17042.0  46.363 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > 


## Next: Follow-up (if interaction is significant)
## plot level data paired comparisons:

# install.packages("emmeans")
library(emmeans)
emmeans(m_veg2, pairwise ~ treatment | lifeform)

## Results: Although all taxa except graminoids and dshrubs had 
## greater cover values in undisturbed, only eshrub, conifer and bryophyte were significantly different.



## Combine with Permanova

## Start from plot-level, vegetation-only data. already (correctly) to exclude substrate.

# Clears the workspace
rm(list=ls())

library(tidyverse)
# install.packages("vegan")
library(vegan)
library(data.table)


## read the csv data file (ensure you are in the correct working directory
df <-  fread(file.path("data","lifeform_matrix_plotlevel.csv")) # command works if you have loaded the data.table library/package

# alternatively, use file.choose command
# df <- read.csv(file.choose()) # choose "lifeform_matrix_plotlevel" from my RScripts working directory

df <- df %>%
  mutate(across(1:8, as.factor)) %>%
  mutate(across(9:ncol(df), as.numeric))

# Display the structure of the modified data frame
str(df)


# 2. Keep only vegetation lifeforms + metadata
lifeforms <- c("forb", "graminoid", "bryophyte",
               "eshrub", "dshrub", "conifer")

df_veg <- df %>%
  select(site, treatment, transnum, quad, all_of(lifeforms))

# 3. Create a unique plot ID (recommended)
# This avoids accidental duplication later.
# plot_id is a new variable

df_veg <- df_veg %>%
  unite(plot_id, site, transnum, quad, remove = FALSE)

# 4. Split community matrix and metadata
# Community matrix (one row per plot, one column per lifeform)

comm_mat <- df_veg %>%
  select(all_of(lifeforms))

# Metadata
# this matrix just contains the info about plot id
meta <- df_veg %>%
  select(plot_id, site, treatment, transnum)

# Step 1: Diagnose (confirm the issue) 
# Run this once to see it explicitly:

rowSums(comm_mat)

# or:

which(rowSums(comm_mat) == 0)

# # Filter out empty plots (the fix) (add this before Step 5)
nonzero <- rowSums(comm_mat) > 0

comm_mat2 <- comm_mat[nonzero, ]
meta2 <- meta[nonzero, ]

any(rowSums(comm_mat2) == 0)
# should return FALSE

# Permutation test for adonis under reduced model

# 5. Run PERMANOVA (with restricted permutations)

# Because plots are nested within transects, permutations must be constrained.
library(vegan)

set.seed(123)

perm <- adonis2(
  comm_mat2 ~ treatment * site,
  data = meta2,
  permutations = 999,
  method = "bray",
  strata = meta2$transnum
)

perm


# ## Not significant
# adonis2(formula = comm_mat2 ~ treatment * site, data = meta2, permutations = 999, method = "bray", strata = meta2$transnum)
# Df SumOfSqs      R2      F Pr(>F)
# Model      5     8.19 0.12146 6.4148      1
# Residual 232    59.24 0.87854              
# Total    237    67.43 1.00000    
# Which is odd because the transect level one was; and llm for plot level is
# However, the LMM answers:
## LMM / GLMM answers “Do individual lifeforms differ in cover by treatment/site?”
# whereas PERMANOVA answers “Does overall lifeform composition differ by treatment/site?”
## conclude that, while lifeforms do differ in cover (abundance) by treatment,
## overall lifeform composition does not differ by treatment/site.

## Did not include the following in manuscript:

# Option to visualize: 
nmds <- metaMDS(comm_mat2, distance = "bray", k = 2, trymax = 100)

plot(nmds, type = "n")
points(nmds, display = "sites",
       col = as.factor(meta2$treatment), pch = 19)


### B: 2nd MODEL option (if plots are uniquely identified)
## load library

m_veg2 <- lmer(
  log(cover + 1) ~ lifeform * treatment + site +
    (1 | transnum/quad),
  data = df_veg
)

anova(m_veg2)


#### END OF SCRIPT -- Extra material below
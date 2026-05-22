######################
# Compaction data
# Goal: to analyse (with 2 Factor ANOVA) and visualize (boxplots, violinplots) the compaction values 
# by treatment and site obtained in 2023 with soil penetrometer. Methods - see manuscript.
# Script: Nina Hewitt, Nov 30 2025
# Sources: 
# 1. Hewitt, N. Geog 374 materials 2017-2025; 
# 2. Soetewey, Antoined (2026) Stats and R, https://statsandr.com/blog/two-way-anova-in-r/ 
# Edited: March 2026 to update treatment/condition labels (NH)

#########################
# clear workspace
rm(list = ls())

##### packages ####

## load libraries
library(data.table) # faster data.frame
library(ggplot2) # plots
library(dplyr)
library(tidyverse)

# BRING IN THE DATA

## read the csv data file (ensure you are in the correct working directory
mydata <-  fread(file.path("data","compaction_2023.csv")) # this command works if you have loaded the data.table library/package

## otherwise load dataset with file.choose command
# df <- read.csv(file.choose()) # compaction_2023.csv" ["Compaction_2023_v2.csv" in local drive]
# Measurements are in columns; sites are rows, UB = Unbrushed, 
# B = Brushed (using a paintbrush to remove top layer of any litter; or disturbed, loose soil if present)
# compact_B_else_UB = contains brushed values except where only UB values were taken/required, since we 
# took both where brushed (loose/unconsolidated debris) was present.
## Note that this datafram has remove the 3 'NA' values (BT-11, plot 6 (big rock); and BT 15,16 (rocky))
# where no measurements were taken.

## LOOK AT OUR DATA
mydata          # look at all your data
head(mydata)    # look at the first few entries
tail(mydata)    # look at the last few entries
str(mydata)     # look at the overall structure

# change columns to AsFactor or numeric accordingly:
mydata <- mydata %>%
  mutate(across(2:5, as.factor)) %>%
  mutate(across(6:8, as.numeric))

# Display the structure of the modified data frame
str(mydata)

# Display the modified data frame
head(mydata)
tail(mydata)

## 2 Factor ANOVA
# Added Jan 2026

## Although subgroups have large sample size (>29), so don't really need to test for normality
table(mydata$treatment, mydata$site)
## We will check nonetheless 

# First, create and save model

# Two-way ANOVA with interaction
mod <- aov(compact_B_else_UB ~ treatment * site,
           data = mydata
)

mod # display results

## Test for Normality with QQ-plot of residuals.
plot(mod, which = 2) # looks fairly normal points follow straight line - diagonal. Some deviation at ends, but is expected

# or QQ-plots the confidence interval around the reference line (Henry’s line).
library(car)
qqPlot(mod$residuals,
       id = FALSE # remove point identification
)
# points follow line and fall within the confidence band, so assume normality

## Histogram of residuals to check normality:
hist(mod$residuals) # relatively normal, some left skew, but minimal

# Shapiro-Wilk normality test:
shapiro.test(mod$residuals) 
# do not reject the null H that the residuals follow a normal distribution p=0.2894

# verify homogeneity of variances or homoscedasticity visually with plot() function:
plot(mod, which = 3) # Spread of residuals is relatively constant
# but red line is somewhat sloping (rather than completely horizontal and flat)
# so while it looks like constant variance assumption reasonably satisfied can test 
# more formally with the Levene’s test (also in the {car} package)

leveneTest(mod)
# Result: must reject the null hypothesis that the variances are equal (p-value = 3.488e-05 ***).

## Now run the model on log values. 

## First, create and save model

# Two-way ANOVA with interaction
mod2 <- aov(log(compact_B_else_UB) ~ treatment * site,
           data = mydata
)

mod2 

# run Levene's test again:
leveneTest(mod2)
# after logging, is normal (p-value = 0.1385)

## Therefore, use log values 

## you have already run the 2 way ANOVA procedure above; run again with log values.
# print results
summary(mod2)
# Both treatment and site are significantly different; and there is significant interaction
# > summary(mod) [log values of compaction]
#                  Df Sum Sq  Mean Sq F value   Pr(>F)    
# treatment         1 14.789  14.789  151.483  < 2e-16 ***
#  site             2  4.572   2.286   23.415 4.65e-10 ***
#  treatment:site   2  1.821   0.911    9.327 0.000123 ***
#  Residuals      254 24.797   0.098                     

# But, since this may be an unbalanced design (there are unequal numbers of subjects in each subgroup -- 
# different sample sizes for site), and in which interaction is significant, run Type 3 Anova:

### Type 3 anova on existing model, where "mod2" is the model
Anova(mod2, type = "III") 
# or same result for:
Anova(mod2, type = 3)

## USED THIS**
# Anova Table (Type III tests)
# Response: log(compact_B_else_UB)
#                 Sum Sq  Df  F value    Pr(>F)    
# (Intercept)     8.2722  1 84.7344  < 2.2e-16 ***
#  treatment      0.7294  1  7.4718  0.0067080 ** 
#  site           0.8108  2  4.1527  0.0168026 *  
#  treatment:site 1.8211  2  9.3271  0.0001234 ***
#  Residuals      24.7969 254 

Anova(mod, type = "III") # unlogged values - since qq plot looks good and large sample (if we ignore Levene's test)

# > Anova(mod, type = "III") # unlogged values
# Anova Table (Type III tests)
# Response: compact_B_else_UB
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)     99.919   1 193.1150 < 2.2e-16 ***
#  treatment        3.113   1   6.0169   0.01484 *  
#  site             2.052   2   1.9827   0.13983    
#   treatment:site  16.463   2  15.9093 3.093e-07 ***
#   Residuals      131.421 254                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### Pairwise comparisons [log values]

#unlogged values:
TukeyHSD(mod,
         which = "site"
)

##Results: Fit: aov(formula = compact_B_else_UB ~ treatment * site, data = mydata)
# $site
# diff     lwr        upr     p adj
# PR-BT 0.67277950  0.3817679 0.9637911 0.0000004
# TM-BT 0.76890752  0.4996354 1.0381797 0.0000000
# TM-PR 0.09612802 -0.1482343 0.3404903 0.6234658

#logged values:

TukeyHSD(mod2,
         which = "site"
)

#Results: Fit: aov(formula = log(compact_B_else_UB) ~ treatment * site, data = mydata)
#$site
# diff         lwr       upr     p adj
# PR-BT 0.31026050  0.18385196 0.4366690 0.0000001
# TM-BT 0.32047326  0.20350781 0.4374387 0.0000000
# TM-PR 0.01021276 -0.09593244 0.1163580 0.9720377


TukeyHSD(mod2,
         which = "treatment:site"
)

#Results:
# > TukeyHSD(mod2,
#            +          which = "treatment:site"
#            + )
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# Fit: aov(formula = log(compact_B_else_UB) ~ treatment * site, data = mydata)
# 
# $`treatment:site`

# diff                                          lwr          upr     p adj
# trampled:BT-reference:BT   0.222412694 -0.0112411489  0.45606654 0.0722703
# reference:PR-reference:BT  0.217057375  0.0003543636  0.43376039 0.0493588
# trampled:PR-reference:BT   0.630191135  0.4134881237  0.84689415 0.0000000
# reference:TM-reference:BT  0.114164474 -0.0864635295  0.31479248 0.5767752
# trampled:TM-reference:BT   0.752178756  0.5520996671  0.95225784 0.0000000
# reference:PR-trampled:BT  -0.005355319 -0.2241829206  0.21347228 0.9999998
# trampled:PR-trampled:BT    0.407778441  0.1889508395  0.62660604 0.0000029
# reference:TM-trampled:BT  -0.108248220 -0.3111691895  0.09467275 0.6440216
# trampled:TM-trampled:BT    0.529766062  0.3273877877  0.73214434 0.0000000
# trampled:PR-reference:PR   0.413133760  0.2125057568  0.61376176 0.0000002
# reference:TM-reference:PR -0.102892901 -0.2860403727  0.08025457 0.5905357
# trampled:TM-reference:PR   0.535121381  0.3525753802  0.71766738 0.0000000
# reference:TM-trampled:PR  -0.516026661 -0.6991741329 -0.33287919 0.0000000
# trampled:TM-trampled:PR    0.121987621 -0.0605583800  0.30453362 0.3929451
# trampled:TM-reference:TM   0.638014282  0.4748749453  0.80115362 0.0000000

# ANOVA different method 
mod3 <- lm(compact_B_else_UB ~ treatment * site,
           data = mydata
)

# print results
summary(mod3) # this may be more appropriate to your data due to unbalanced design


# Visualizations

## Updated, Nov, 2025

## create boxplots of plot values by transect and treatment
## LOOK AT OUR DATA with BOXPLOTS
# with base R:
boxplot(compact_B_else_UB ~ treatment, 
        data = mydata, col = "lightgray",
        varwidth = TRUE, notch = FALSE, 
        main = "Soil compaction levels by trampling condition",
        ylab = "Penetration resistance (kg/cm^2)")

# with dotplot
stripchart(data = mydata, compact_B_else_UB ~ treatment, vertical = TRUE, method = "jitter", pch = 16, col = 'purple', add = TRUE)

ggplot(data = mydata, aes(x = treatment, y = compact_B_else_UB)) +
  geom_boxplot() +
  xlab("Condition") +
  ylab(expression(Penetration~resistance~(kg/cm^{2}))) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "steelblue", size = 1.5) +
  theme(
    panel.background = element_blank(), # Optional: removes the background color
    axis.line = element_line(colour = "black")) # Optional: adds a black axis line)

site_palette <- c(
  "TM" = "#E69F00",  # orange
  "BT" = "#009E73",  # green
  "PR" = "#0072B2"   # blue
)

## Color code points by site and using viridis palette (for color-blindness): 
(plot_compact_export <- ggplot(data = mydata, aes(x = treatment, y = compact_B_else_UB)) +
  geom_boxplot(outliers = F) +
  xlab("Condition") +
  ylab(expression(Penetration~resistance~(kg/cm^{2}))) +
  geom_jitter(
    aes(color = site),
    width = 0.2,
    alpha = 0.7,
    size = 1.75
  ) +
  scale_color_manual(breaks = c("TM","BT","PR"),values =site_palette)+   # viridis palette + legend title
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ))
ggsave(file.path("figure","Fig_Compaction_boxplot.jpg"),plot_compact_export,
       unit = "mm",width = 180,height = 120,dpi = 250,scale = 0.8)
###END####
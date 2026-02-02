######################
# Compaction data
# Goal: to analyse (with 2 Factor ANOVA) and visualize (boxplots, violinplots) the compaction values 
# by treatment and site obtained in 2023 with soil penetrometer. Methods - see manuscript.
# Script by Nina Hewitt, Nov 30 2025
# Sources: 
# 1. Hewitt, N. Geog 374 materials 2017-2025; 
# 2. Soetewey, Antoined (2026) Stats and R, https://statsandr.com/blog/two-way-anova-in-r/ 
# Edited: Jan 2026 (NH)

# note to self - to comment or uncomment: highlight chunk and press Cmd + Shift + C 

#########################
# clear workspace
rm(list = ls())

## load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)

# load data
## Either use command (and adjust file location and name)
# compaction.data <- read.csv(file = "./data/raw_data/data_2023 - Compaction475.csv") ## old file - adjust
# site.data <- read.csv(file = "./data/raw_data/data_2023 - Site_data.csv")

# or use command:
mydata <- read.csv(file.choose())
# read csv "Compaction_2023_v2.csv" into R (has new single column for "compact_B_else_UB"); 
# Measurements are in columns; sites are rows, UB = Unbrushed,
# B = Brushed (using a paintbrush to remove top layer of any litter; or disturbed, loose soil if present)

## LOOK AT OUR DATA
mydata          # look at all your data
head(mydata)    # look at the first few entries
tail(mydata)    # look at the last few entries
str(mydata)     # look at the overall structure

# change columns to AsFactor or numeric accordingly:
mydata <- mydata %>%
  mutate(across(2:5, as.factor)) %>%
  mutate(across(6:10, as.numeric))

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

mod

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

# verifyomogeneity of variances or homoscedasticity visually with plot() function:
plot(mod, which = 3) # Spread of the residuals is relatively constant
# but the red line is somewhat sloping (rather than perfectly horizontal and flat)
# so while looks like constant variance assumption satisfied can test more formally with the Levene’s test
# (also from the {car} package)

leveneTest(mod)
# We must reject the null hypothesis that the variances are equal (p-value = 3.488e-05 ***).

# Now run the model on log values. 
# First, create and save model

# Two-way ANOVA with interaction
mod2 <- aov(log(compact_B_else_UB) ~ treatment * site,
           data = mydata
)

mod2 # should overwrite original

# run Levene's test again:
leveneTest(mod2)
# after logging, is normal (p-value = 0.1385)

## Use log values 

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
# different sample sizes for site), and in which interaction is significant, should use Type 3 Anova:

### Run a type 3 anova on the existing model, where "mod2" is the model
Anova(mod2, type = "III") 
# or same result for:
Anova(mod2, type = 3)

#Anova Table (Type III tests)
#Response: log(compact_B_else_UB)
#Sum Sq  Df  F value    Pr(>F)    
# (Intercept)    16.2049   1 165.9909 < 2.2e-16 ***
#  treatment       0.7294   1   7.4718 0.0067080 ** 
#  site            5.5821   2  28.5893 6.334e-12 ***
#  treatment:site  1.8211   2   9.3271 0.0001234 ***
#  Residuals      24.7969 254 

Anova(mod, type = "III") # unlogged values - since qq plot looks good and large sample (ignore Levene's test)

# > Anova(mod, type = "III") # unlogged values
# Anova Table (Type III tests)
# Response: compact_B_else_UB
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)    151.347   1 292.5114 < 2.2e-16 ***
#   treatment        3.113   1   6.0169   0.01484 *  
#   site            39.208   2  37.8893 3.978e-15 ***
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
# diff        lwr       upr     p adj
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

## Nicer, but need to remove the 3 'NA' values from dataframe (BT-11, plot 6 (big rock); and BT 15,16 (rocky))
# upload and start again with file saved as "Compaction_2023_v2.csv"
# with ggplot
ggplot(data = mydata, aes(x = treatment, y = compact_B_else_UB)) +
  geom_boxplot() +
  xlab("Condition") +
  ylab(expression(Penetration~resistance~(kg/cm^{2}))) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "steelblue", size = 1.5) +
  theme(
    panel.background = element_blank(), # Optional: removes the background color
    axis.line = element_line(colour = "black")) # Optional: adds a black axis line)

## Need to use dataframe without 'NA' values, _vw

## Update, Dec 4, 2025 color code points by site and using viridis palette (for color-blindness): 
ggplot(data = mydata, aes(x = treatment, y = compact_B_else_UB)) +
  geom_boxplot() +
  xlab("Condition") +
  ylab(expression(Penetration~resistance~(kg/cm^{2}))) +
  geom_jitter(
    aes(color = site),
    width = 0.2,
    alpha = 0.7,
    size = 1.5
  ) +
  scale_color_viridis_d(name = "Site") +   # viridis palette + legend title
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

# Stripcharts only

# Create plot that contains one strip chart per variable:
# Look at, e.g., compaction by treatment

stripchart(compact_B_else_UB ~ treatment,
           data = mydata,
           main = 'Soil compaction levels by treatment',
           xlab = 'Treatment', 
           ylab = 'Penetration resistance (kg/cm^2)',
           col = c('steelblue', 'purple'),
           pch = 16,
           method = 'jitter',
           vertical = TRUE)


## Stopped here.

###END####

## additional notes for creating violinplots and stripcharts

# View documentation on stripcharts

?stripchart

# a strip chart for each variable in a single plot:

#create list of variables
x <- list('Treatment' = mydata$treatment, 'Compaction' = mydata$compact_B_else_UB)

# show x
x

# I got the full list. What could I have typed instead? Hint: see above

# Create plot that contains one strip chart per variable

# Look at, e.g., compaction by treatment, plot vertically

stripchart(compact_B_else_UB ~ treatment,
           data = mydata,
           main = 'Compaction by condition',
           xlab = 'Treatment', 
           ylab = 'Compaction (kg/cm^sq)',
           col = c('steelblue', 'purple'),
           pch = 16,
           method = 'jitter',
           vertical = TRUE)

# View documentation on violin plots

?geom_violin

## Create Basic violin plots

library(ggplot2)
# Basic violin plot
q <- ggplot(mydata, aes(x=treatment, y=compact_B_else_UB)) + 
  geom_violin()
q

# Rotate the violin plot
q + coord_flip()

# Set trim argument to FALSE
ggplot(dat, aes(x=treatment, y=compact_B_else_UB)) + 
  geom_violin(trim=FALSE)


# violin plot with mean points
q + stat_summary(fun=mean, geom="point", shape=23, size=2)

# violin plot with median points
q + stat_summary(fun=median, geom="point", size=2, color="red")

# Is the data left or right skewed? Hint: if the median is lower
# on the scale than the mean, they are right skewed, and visa versa

## Customize legend, etc

# Basic violin plot
ggplot(mydata, aes(x=treatment, y=compact_B_else_UB)) + 
  geom_violin(trim=FALSE, fill="gray") +
  labs(title="Compaction by condition",x="Treatment", y = "Compaction (kg/cm^sq)") +
  geom_boxplot(width=0.1) +
  theme_classic()

# Change color by groups
dp <- ggplot(mydata, aes(x=treatment, y=compact_B_else_UB)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  labs(title="Compaction by condition",x="Treatment", y = "Compaction (kg/cm^sq)")
dp + theme_classic()

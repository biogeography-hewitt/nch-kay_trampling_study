##### packages ####
library(data.table) # faster data.frame
library(brms) # bayesian model fitting using stan
library(ggplot2) # plots
library(stringr) # string manipulation
library(vegan) # diversity indices computation
#### how does trampling affect species richness ? ####

## reading the csv, cleaning some values
data_t <- fread(file.path("nch-kay_trampling_study", "data","Trampling_Data_2024_Clean_nozeros.csv"))
data_long <- melt(data_t,id.vars = colnames(data_t)[1:9]) 
data_long[,value:=ifelse(value == "*","",value)]
data_long[,value:=ifelse(value == "","0",value)]
data_long[,value:=as.numeric(value)]
data_long[,value:=ifelse(is.na(value),0 ,value)]
data_long[,.N,by = variable]

## removing "species" not of interest
data_long <- data_long[!variable%in% c("soil","litter","rock","crust","bryosp2","bryosp1"),]

table(data_long$value)
table(data_long$variable)

#### Richness ####
## compute richness as species present in a subplot, and shannon diversity
data_richness <- data_long[,.(richness = sum(value!=0), 
                              shannon = diversity(value,"shannon")), by = eval(colnames(data_long)[1:9])]
data_richness[,altitude_scaled := scale(altitude)] # center and scale the predictor incase we want to use it
data_richness[,site := substr(transect,1,2)]## getting the site variable back with string manipulation
data_richness[,transect_pair := paste0(site,"_",str_pad(trans.pair, width=2, side="left", pad="0"))]# naming the transect for plots

table(data_richness$transect)
table(data_richness$treatment)
table(data_richness$commun)
table(data_richness$trans.pair)

## the distribution look like poisson or negbinomial
ggplot(data_richness,aes(x = richness,fill = transect))+
  facet_wrap(~treatment)+
  geom_histogram(binwidth = 1)+
  theme_classic()

## Bayesian model 

# setting priors, knowledge we know already to aid the fit
prior_rich <- set_prior("normal(2,0.25)",class = "Intercept") #I already know the mean species richness 
# in control plot is around exp(2)
prior_rich <- c(prior_rich,set_prior("normal(500,20)", class = "shape" ))#I fix the shape paramter very hight to help the 
# covnergence, as this parameter can go to infinity without any improvment of the fit

# 1 means varying intercept(control richness) across site / pair / transect.
# 1 + treatment means varying effect of treatment across pairs
model_rich <-  brm(richness ~ treatment  +(1|site) +  (1 + treatment|trans.pair) +(1|transect) ,
                                  data = data_richness,
                                  family = negbinomial(), # distribution family
                                  backend = "cmdstanr", # Need to be installed, use backend = "rstan" works aswell
                                  iter = 2000, # number of iteration of the algo per chains (models)
                                  prior = prior_rich, # the priors
                                  warmup = 500, # number of iterations discarded 
                                  chains = 3, # number of model to fit, they need to converge to a unique solution
                                  cores = 3, # using a CPU core per model to fasten the computation
                                  control = list(adapt_delta = 0.99), # increase the computation time to better estimate parameters, need if divergent transition
                                  threads = threading(3), # this number can go up or down to fasten computation,
                                  init = 0, # lower the risk of crash
                                  file = file.path("model","richness_model_4")) # where the model is stored
# will load the model instead of fitting one if the file already exists

pp_check(model_rich) # posterio predictive check to make sur we fit the real distribution
summary(model_rich) # summary of fixed and random effects
plot(model_rich) # checking the convergence: fuzzy caterpillar, the chains have converged

sjPlot::plot_model(model_rich,type = "pred",terms = "treatment")

## getting the model prediction back for every transect
preds_rich <- fitted(model_rich,
                     newdata = data_richness[,.N, by = .(treatment,site,trans.pair,transect,altitude_scaled)],
                     re_formula = ~  (1 |site) + (1+ treatment|trans.pair)  +(1|transect)) # I specifiy here the heirarchy I want in the preds


preds_rich <- cbind(preds_rich,data_richness[,.N, by = .(treatment,site,trans.pair,transect_pair,transect,altitude_scaled)])

# plot of pred vs real div
(rich_pred_plot <- ggplot(data_richness,aes (x = transect_pair, y = richness , fill = treatment,group =treatment))+
  theme_classic()+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  geom_jitter(pch = 21, size = 2,alpha = 0.5,
              position = position_jitterdodge(0.2,0.1,0.6))+
  geom_segment(data = preds_rich, aes(y = Q2.5, yend = Q97.5),position = position_dodge(0.6))+
  geom_point(position = position_dodge(0.6),data = preds_rich, aes(y = Estimate),pch = 21, size = 4)+
  scale_color_manual(values = trail_color <- c("#27A81E", "#DEBF50"), labels = trail_labs <- c("Far (<5m)","Trampled"))+
  scale_fill_manual(values = trail_color)+
  labs(x = "Transect pair", y = "Species richness", color = lab_trt <- "Position along\nthe trail",fill = lab_trt)
)

# export
ggsave(file.path("figure","Fig_richness_transect.jpg"),rich_pred_plot,unit = "cm",width = 18,height = 11,dpi = 300)

ggplot(data_richness,aes( x = treatment, y = richness,fill = site))+
  geom_violin()+
  theme_classic()



#### how does trampling affect Shannon diversity? ####

ggplot(data_richness,aes(x = shannon,fill = transect))+
  facet_wrap(~treatment)+
  geom_histogram(binwidth = 0.1)+
  theme_classic()

prior_shannon <- set_prior("normal(1,0.25)",class = "Intercept")
prior_shannon <- c(prior_shannon,set_prior("student_t(3, 0.5, 0.5)", class = "sd", group = "trans.pair"))
prior_shannon <- c(prior_shannon,set_prior("student_t(3, 0.5, 0.5)", class = "sd", group = "site"))
prior_shannon <- c(prior_shannon,set_prior("student_t(3, 0.5, 0.5)", class = "sd", group = "transect"))
prior_shannon <- c(prior_shannon,set_prior("normal(500,300)", class = "shape" ))

model_shannon <-  brm(shannon ~ treatment  +(1 |site) +  (1+ treatment|trans.pair) +(1|transect) ,
                   data = data_richness,
                   family = student(),
                   backend = "cmdstanr", 
                   iter = 2000,
                   prior = prior_shannon,
                   warmup = 500,
                   chains = 3,
                   cores = 3,
                   control = list(adapt_delta = 0.99),
                   threads = threading(3),
                   init = 0,
                   file = file.path("model","shannon_model"))

conditional_effects(model_shannon,"treatment")
pp_check(model_shannon)+coord_cartesian(xlim = c(0,3))
summary(model_shannon)
plot(model_shannon)

preds_shannon <- fitted(model_shannon,
                     newdata = data_richness[,.N, by = .(treatment,site,trans.pair,transect,altitude_scaled)],
                     re_formula = ~  (1 |site) + (1+ treatment|trans.pair)  +(1|transect))


preds_shannon <- cbind(preds_shannon,data_richness[,.N, by = .(treatment,site,trans.pair,transect_pair,transect,altitude_scaled)])
preds_shannon[,shannon := Estimate]

(shannon_pred_plot <- ggplot(data_richness,aes (x = transect_pair, y = shannon , fill = treatment,group =treatment))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
    geom_jitter(pch = 21, size = 2,alpha = 0.5,
                position = position_jitterdodge(0.2,0.1,0.6))+
    geom_segment(data = preds_shannon, aes(y = Q2.5, yend = Q97.5),position = position_dodge(0.6))+
    geom_point(position = position_dodge(0.6),data = preds_shannon, aes(y = Estimate),pch = 21, size = 4)+
    scale_color_manual(values = trail_color <- c("#27A81E", "#DEBF50"), labels = trail_labs <- c("Far (<5m)","Trampled"))+
    scale_fill_manual(values = trail_color)+
    labs(x = "Transect pair", y = "Shannon diversity index", color = lab_trt <- "Position along\nthe trail",fill = lab_trt)
)

ggsave(file.path("figure","Fig_Shannon_transect.jpg"),shannon_pred_plot,unit = "cm",width = 18,height = 11,dpi = 300)


library(data.table)
library(brms)
library(ggplot2)
library(stringr)
library(vegan)
## how does trampling affect species richness ?

data_t <- fread(file.path("data","Trampling_Data_2024_Clean_nozeros.csv"))
data_long <- melt(data_t,id.vars = colnames(data_t)[1:9]) 
data_long[,value:=ifelse(value == "*","",value)]
data_long[,value:=ifelse(value == "","0",value)]
data_long[,value:=as.numeric(value)]
data_long[,value:=ifelse(is.na(value),0 ,value)]
data_long[,.N,by = variable]

data_long <- data_long[!variable%in% c("soil","litter","rock","crust","bryosp2","bryosp1"),]


ggplot(data_long,aes(x = value,fill = treatment))+
  geom_histogram(position = "dodge")+
  facet_wrap(~ variable,scales = "free_y")+
  theme_classic()


table(data_long$value)

#### Richness ####
data_richness <- data_long[,.(richness = sum(value!=0), shannon = diversity(value,"shannon")), by = eval(colnames(data_long)[1:9])]
data_richness[,altitude_scaled := scale (altitude)]
data_richness[,site := substr(transect,1,2)]

table(data_richness$transect)
table(data_richness$treatment)
table(data_richness$commun)
table(data_richness$trans.pair)

ggplot(data_richness,aes(x = richness,fill = transect))+
  facet_wrap(~treatment)+
  geom_histogram(binwidth = 1)+
  theme_classic()

library(brms)

prior_rich <- set_prior("normal(2,0.25)",class = "Intercept")
prior_rich <- c(prior_rich,set_prior("student_t(3, 0.5, 0.5)", class = "sd", group = "trans.pair"))
prior_rich <- c(prior_rich,set_prior("student_t(3, 0.5, 0.5)", class = "sd", group = "site"))
prior_rich <- c(prior_rich,set_prior("student_t(3, 0.5, 0.5)", class = "sd", group = "transect"))
prior_rich <- c(prior_rich,set_prior("normal(500,300)", class = "shape" ))


data_richness[,.N, by = .(treatment,site,trans.pair,transect,altitude_scaled)]
data_richness[,transect_pair := paste0(site,"_",str_pad(trans.pair, width=2, side="left", pad="0"))]

default_prior(model_rich)

model_rich <-  brm(richness ~ treatment  +(1+ |site) +  (1+ treatment|trans.pair) +(1|transect) ,
                                  data = data_richness,
                                  family = negbinomial(),
                                  backend = "cmdstanr", 
                                  iter = 2000,
                                  prior = prior_rich,
                                  warmup = 500,
                                  chains = 3,
                                  cores = 3,
                                  control = list(adapt_delta = 0.98),
                                  threads = threading(3),
                                  init = 0,
                                  file = file.path("model","richness_model")  )

pp_check(model_rich)
summary(model_rich)
plot(model_rich)

sjPlot::plot_model(model_rich,type = "pred",terms = "treatment")

preds_rich <- fitted(model_rich,
                     newdata = data_richness[,.N, by = .(treatment,site,trans.pair,transect,altitude_scaled)],
                     re_formula = ~  (1 |site) + (1+ treatment|trans.pair)  +(1|transect))


preds_rich <- cbind(preds_rich,data_richness[,.N, by = .(treatment,site,trans.pair,transect_pair,transect,altitude_scaled)])

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

ggsave(file.path("figure","Fig_richness_transect.jpg"),rich_pred_plot,unit = "cm",width = 18,height = 11,dpi = 300)
ggplot(data_richness,aes( x = treatment, y = richness,fill = site))+
  geom_violin()+
  theme_classic()



#### Shanon div ####
data_richness[,transect_pair := paste0(site,"_",str_pad(trans.pair, width=2, side="left", pad="0"))]

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
                   file = file.path("model","shannon_model")  )
data_richness
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
    labs(x = "Transect pair", y = "Hannon diversity index", color = lab_trt <- "Position along\nthe trail",fill = lab_trt)
)

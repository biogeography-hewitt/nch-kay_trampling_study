library(data.table)
library(brms)
library(ggplot2)
library(stringr)

## how does trampling affect species cover ?
# some species has never been recorded, only 0, removing them
data_cover <- data_long
data_cover[,species := variable];data_cover[,variable := NULL]
data_cover[,cover := value];data_cover[,value := NULL]

data_cover[,never_found := all(cover == 0),by = species]
data_cover <- data_cover[never_found == F,]; data_cover[,never_found := NULL]
data_cover[,altitude_scaled := scale (altitude)]
data_cover[,site := substr(transect,1,2)]
data_cover[,transect_pair := paste0(site,"_",str_pad(trans.pair, width=2, side="left", pad="0"))]


data_cover[cover != 0 ,length(unique(species)),by = .(treatment,transect_pair)]
data_cover[cover != 0 ,length(unique(species)),by = .(transect_pair)]

sp_pool_of_site <- data_cover[cover != 0 ,.(species = unique(species)),by = .(transect_pair)]

## remove from the sites species that has never been found in the site / transect
data_cover_trim <- merge(data_cover,sp_pool_of_site)
data_cover_trim[,species := droplevels(species)]
data_cover_trim[,cover_1 := ceiling(cover)]
data_cover_trim[,cover_1 := cover / 100]


ggplot(data_cover_trim,aes( x = cover))+geom_histogram(binwidth = 1)+theme_classic()
ggplot(data_cover_trim,aes( x = cover, fill = species))+geom_histogram(binwidth = 1)+theme_classic()+facet_wrap(~ treatment)

prior_cover <- set_prior("normal(-3,2)",class = "Intercept",dpar = "mu")
prior_cover <- c(prior_cover,set_prior("normal(1,1)",class = "Intercept",dpar = "zi"))
prior_cover <- c(prior_cover,set_prior("student_t(3, 0.5, 0.5)", class = "sd", group = "trans.pair"))
#prior_cover <- c(prior_cover,set_prior("student_t(3, 0.5, 0.5)", class = "sd", group = "site"))
prior_cover <- c(prior_cover,set_prior("student_t(3, 0.5, 0.5)", class = "sd", group = "transect"))
prior_cover <- c(prior_cover,set_prior("normal(10,20)", class = "phi" ))

formula_model <- bf(cover_1 ~ treatment + (1|trans.pair) + (1|transect) + (treatment|species) + (1|transect:species)+ (1|trans.pair:species),
   zi ~   treatment + (1|trans.pair) + (1|transect) + (treatment|species) + (1|transect:species)+ (1|trans.pair:species))

model_cover <- brm(formula_model,
                   data = data_cover_trim,
                   family = zero_inflated_beta(),
                   backend = "cmdstanr", 
                   iter = 1000,
                   prior = prior_cover ,
                   warmup = 200,
                   chains = 3,
                   cores = 3,
                   threads = threading(4),
                   file = file.path("model","ZIB_species_cover_model"),
                   init = 0)



pp_check(model_cover)+coord_cartesian(xlim = c(0,0.5))
summary(model_cover)
plot(model_cover)
data_cover[,sum(cover != 0),by = .(treatment,species)]

sjPlot::plot_model(model_cover,type = "pred",terms = c("species","treatment"))
sjPlot::plot_model(model_cover,type = "pred",terms = c("treatment"), dpar = 'mu')
sjPlot::plot_model(model_cover,type = "pred",terms = c("treatment"), dpar = 'hu')

conditional_effects(model_cover,effects = c("species:treatment"),re_formula = ~ (1|species))
conditional_effects(model_cover,effects = c("species:treatment"),re_formula = ~ (1|species),dpar = "mu")
conditional_effects(model_cover,effects = c("species:treatment"),re_formula = ~ (1|species),dpar = "hu")

conditional_effects(model_cover,effects = c("treatment"),re_formula = NA)
conditional_effects(model_cover,effects = c("treatment"),re_formula = NA,dpar = "mu")
conditional_effects(model_cover,effects = c("treatment"),re_formula = NA,dpar = "zi")

preds_cover <- fitted(model_cover,
                      newdata =data_cover_trim[,.N, by = .(treatment,site,trans.pair,transect,transect_pair,species)],
                      re_formula = ~ (1|trans.pair) + (1|transect) + (treatment|species) + (1|transect:species)+ (1|trans.pair:species))

preds_cover <- cbind(data_cover_trim[,.N, by = .(treatment,site,trans.pair,transect,transect_pair,species)],preds_cover)

ggplot(data_cover_trim[cover!=0,],aes( x = species, y = cover,group = treatment,fill = treatment))+
  geom_point(pch = 21, size = 2,alpha = 0.5,
              position = position_jitterdodge(0.2,0.1,0.6))+
  #facet_wrap(~site, scale = "free_x",nrow = 3)+
  theme_classic()

data_cover_trim[,species_order := as.character(species)]
new_frame <- data_cover_trim[,.N ,by = .(treatment,species,species_order),]
pred_sp <- cbind(new_frame,fitted(model_cover,
                  newdata =new_frame,
                  re_formula = ~ (1|species) ) * 100)


(figure_pred_sp <- ggplot(data_cover_trim[,],aes( x = species_order, y = cover,group = treatment,fill = treatment, color = treatment))+
  geom_point( size = 1,alpha = 0.5, position = position_jitterdodge(0.2,0.1,0.6),fill = NA)+
  geom_pointrange(aes(y = Estimate,ymin = Q2.5,ymax = Q97.5 ),pred_sp,
                  color = "grey5",pch = 21,size = 0.5,lwd = 0.5,lineend = "round",position = position_dodge(0.6))+
 theme_classic()+
   theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),legend.position = c(0.2,0.8))+
  
  scale_color_manual(values = trail_color <- c("#27A81E", "#DEBF50"), labels = trail_labs <- c("Far (<5m)","Trampled"))+
  scale_fill_manual(values = trail_color, labels = trail_labs )+
    labs(x = "", y = "Cover (%)",color = lab_trt,fill = lab_trt ))

ggsave(file.path("figure","Fig_cover_sp.jpg"),figure_pred_sp,unit = "cm",width = 18,height = 14,dpi = 300,scale = 1.25)

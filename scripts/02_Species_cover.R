##### packages ####
library(data.table) # faster data.frame
library(brms) # bayesian model fitting using stan
library(ggplot2) # plots
library(stringr) # string manipulation
library(vegan) # diversity indices computation
#### how does trampling affect species cover ? ####
data_cover <- data_long
data_cover[,species := variable];data_cover[,variable := NULL] # renaming
data_cover[,cover := value];data_cover[,value := NULL]


# some species has never been recorded, only 0, removing them
data_cover[,never_found := all(cover == 0),by = species]
data_cover <- data_cover[never_found == F,]; data_cover[,never_found := NULL]

## site naming and scaling altitude 
data_cover[,altitude_scaled := scale (altitude)]
data_cover[,site := substr(transect,1,2)]
data_cover[,transect_pair := paste0(site,"_",str_pad(trans.pair, width=2, side="left", pad="0"))]


data_cover[cover != 0 ,length(unique(species)),by = .(treatment,transect_pair)]
data_cover[cover != 0 ,length(unique(species)),by = .(transect_pair)]

## ASSUMPTION: A species is in the local pool of a pair of transect if it has been recorded once
## CONSEQUENCE: species that has not been recorded once at a pair are REMOVED from the data
## CONSEQUENCE: Species recorded once in a pair will be 0 (true absence) when not recorded in other plots of the pair
sp_pool_of_site <- data_cover[cover != 0 ,.(species = unique(species)),by = .(transect_pair)]

## remove from the sites species that has never been found in the site / transect
data_cover_trim <- merge(data_cover,sp_pool_of_site)
data_cover_trim[,species := droplevels(species)]
data_cover_trim[,cover_1 := cover / 100] # beta distrib needs cover between 0 and 1 

##looks like zero inflated beta
ggplot(data_cover_trim,aes( x = cover_1, fill = species))+geom_histogram(binwidth = 1)+theme_classic()+facet_wrap(~ treatment)

#### Modelling 

## setting the prior knowledge: 

prior_cover <- set_prior("normal(-3,1)",class = "Intercept",dpar = "mu") # I expect the average plant to be tiny : cover of inv_logit(-3): 0.05 , 5 % cover
prior_cover <- c(prior_cover,set_prior("normal(10,10)", class = "phi" )) # Setting the shape parameter to be wide to accommodate for the data dispertion
# those are wide prior, so nuch much informations given here

formula_model <- bf(cover_1 ~ treatment + (1|trans.pair) + (1|transect) + (1 + treatment|species) + (1|transect:species)+ (1|trans.pair:species),
   zi ~   treatment + (1|trans.pair) + (1|transect) + (1 + treatment|species) + (1|transect:species)+ (1|trans.pair:species))

model_cover <- brm(formula_model, # takes a while to compute 
                   data = data_cover_trim,
                   family = zero_inflated_beta(),
                   backend = "cmdstanr", # "rstan" if not working
                   iter = 1000,
                   prior = prior_cover ,
                   warmup = 200,
                   chains = 3,
                   cores = 3,
                   threads = threading(4),
                   file = file.path("model","ZIB_species_cover_model"),
                   init = 0)


## distribution is fitting as well as possible, but the not-continuous cover is making it harder to see
pp_check(model_cover)+coord_cartesian(xlim = c(0,0.5))
summary(model_cover)
plot(model_cover)## parameters converged because the caterpillar are fuzzy

## rough plots to look at the effect of trampling per species
conditional_effects(model_cover,effects = c("species:treatment"),re_formula = ~ (1 + treatment|species))
conditional_effects(model_cover,effects = c("species:treatment"),re_formula = ~ (1 + treatment|species),dpar = "mu")
conditional_effects(model_cover,effects = c("species:treatment"),re_formula = ~ (1 + treatment|species),dpar = "hu")

## rough plots to look at the general effect of trampling
conditional_effects(model_cover,effects = c("treatment"),re_formula = NA) # average of prob * cover
conditional_effects(model_cover,effects = c("treatment"),re_formula = NA,dpar = "mu") # cover
conditional_effects(model_cover,effects = c("treatment"),re_formula = NA,dpar = "zi")  # prob of absence

## trick for alphabetical order
data_cover_trim[,species_order := as.character(species)]

## getting prediction from the model, only asking for species specific random effects, everything else assumed average
new_data_for_pred <- data_cover_trim[,.N ,by = .(treatment,species,species_order),]

pred_sp <- cbind(new_frame,fitted(model_cover,
                  newdata =new_data_for_pred, # no dpar argument : the whole model prob * cover, dpar = "mu" cover, dpar= "zi" probability of absence 
                  re_formula = ~ (1+treatment|species) ) * 100)

## plots! can be adapted to show full model, just zi, just mu
(figure_pred_sp <- ggplot(data_cover_trim[,],aes( x = species_order, y = cover,group = treatment,fill = treatment, color = treatment))+
  geom_point( size = 1,alpha = 0.5, position = position_jitterdodge(0.2,0.1,0.6),fill = NA)+
  geom_pointrange(aes(y = Estimate,ymin = Q2.5,ymax = Q97.5 ),pred_sp,
                  color = "grey5",pch = 21,size = 0.5,lwd = 0.5,lineend = "round",position = position_dodge(0.6))+
  theme_classic()+
    scale_y_continuous(limits =  c(0,50),oob = scales::squish)+ ##squishing the very high cover, can be removed
   theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
         legend.position = c(0.1,0.8),
         axis.title.x = element_blank())+
  
  scale_color_manual(values = trail_color <- c("#27A81E", "#DEBF50"), labels = trail_labs <- c("Far (<5m)","Trampled"))+
  scale_fill_manual(values = trail_color, labels = trail_labs )+
    labs(x = "", y = "Cover (%)",color = lab_trt,fill = lab_trt ))

# export
ggsave(file.path("figure","Fig_cover_sp.jpg"),figure_pred_sp,unit = "cm",width = 18,height = 14,dpi = 300,scale = 1.25)

## getting the raw draws from the model to compute delta predictions
new_data_for_pred <- data_cover_trim[,.N ,by = .(treatment,species,species_order),]

pred_sp <- data.table(fitted(model_cover,
       newdata =new_data_for_pred, # no dpar argument : the whole model prob * cover, dpar = "mu" cover, dpar= "zi" probability of absence 
       re_formula = ~ (1+treatment|species), summary = F) * 100)
pred_sp <- t(pred_sp)
pred_sp <- cbind(new_frame[,],pred_sp)

pred_sp <- melt(pred_sp,id.vars = c("treatment","species","species_order","N"))
pred_sp[,draws_id := variable ] ; pred_sp[,variable := NULL ]  
pred_sp[,pred_cover := value ] ; pred_sp[,value := NULL ]  

pred_sp <- pred_sp[,.(delta_cover = diff(-pred_cover)),by = .(species,species_order,draws_id,N)]
pred_sp_summarized <- pred_sp[,.(Estimate = mean(delta_cover),
                                 Q2.5 = quantile(delta_cover,probs = 0.025),
                                 Q97.5  = quantile(delta_cover,probs = 0.975),
                                 prob_of_being_signif = mean(delta_cover<0)),by = .(species,species_order,N)]


pred_sp_summarized[,signif := ifelse(sign(Q2.5) == sign(Q97.5),"Signif","Non")]
pred_sp_summarized[,signif_2way := ifelse(sign(Q2.5) == sign(Q97.5) & sign(Estimate) == -1,"Decline",
                                             ifelse(sign(Q2.5) == sign(Q97.5) & sign(Estimate) == 1,"Increase","Non signif"))]

(ggplot(pred_sp_summarized,aes( x= reorder(species_order,Estimate), y= Estimate, ymin = Q2.5, ymax  = Q97.5,color = signif))+
    geom_hline(yintercept = 0, lty = 2, color = "grey20")+
    geom_pointrange()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          legend.position = c(0.7,0.3),
          axis.title.x = element_blank())+
    scale_color_manual(values = c("grey60","grey10"))+
    labs( y= "Distribution of Δcover")
  )



s
####alternative probabilistic representation
(ggplot(pred_sp_summarized,aes( x= reorder(species_order,prob_of_being_signif), y= prob_of_being_signif,color = signif))+
    geom_hline(yintercept = 0, lty = 2, color = "grey20")+
    geom_point()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          legend.position = c(0.7,0.3),
          axis.title.x = element_blank())+
    scale_color_manual(values = c("grey60","grey10"))+
    labs( y= "Distribution of Δcover")
)


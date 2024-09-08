#######
# Field Population
# IT 
#  Diet work
######

rm(list=ls())

library(tidyverse) 
library(mgcv)
library(MuMIn)
library(multcomp)
library(gratia)
library(ggpubr)
library(patchwork)
library(broom)
library(knitr)

#Function
std <- function(x) sd(x)/sqrt(length(x))

dat <- read.csv("2015_Australia/data/raw/IT_data/IT_dat.csv") #Replace with wherever you put your data
str(dat)
dat$Population <- as.factor(dat$Population)
dat$Sex <- as.factor(dat$Sex)
dat$Diet.Pair <- as.factor(dat$Diet.Pair)

dat <- dat %>% drop_na("Protein.consumed..g.")


names(dat)
#raw data viz

ggplot(dat,aes(x=(Protein.consumed..g.),y=(Carb.consumed..g.),color=Population)) + geom_point() +
  coord_equal(ratio=1) + xlim(0,.05) + ylim(0,.05)



# non-linear relationship with initial GH mass
mod1 <- gam(list( 
  Carb.consumed..g. ~ Population + Diet.Pair + Sex + s(Initial.mass..g.,k=30),
  Protein.consumed..g. ~ Population + Diet.Pair + Sex + s(Initial.mass..g.,k=30)),
  family=mvn(d=2),select=TRUE, data=dat
)

summary(mod1) #checking results
plot(mod1, all.terms=TRUE) #plotting ALL graphs graphs with a '.1' signify protein consumption
k.check(mod1) #These p-values should NOT be significant
gam.check(mod1)

#linear relationship with initial GH mass
mod2 <- gam(list(
  Carb.consumed..g. ~ Population + Diet.Pair + Sex + (Initial.mass..g.),
  Protein.consumed..g. ~ Population + Diet.Pair + Sex + (Initial.mass..g.)),
  family=mvn(d=2),select=TRUE, data=dat
)

summary(mod1)
plot(mod2, all.terms=TRUE)

#model without initial mass
mod3 <- gam(list(
  Carb.consumed..g. ~ Population + Sex + Diet.Pair ,
  Protein.consumed..g. ~ Population + Sex + Diet.Pair),
  family=mvn(d=2),select=TRUE, data=dat
)

mod3$edf1
lot(mod3, all.terms=TRUE)

#NULL model (this should be the worst model)
mod4 <- gam(list(
  Carb.consumed..g. ~ 1,
  Protein.consumed..g. ~ 1),
  family=mvn(d=2),select=TRUE, data=dat
)

summary(mod1)
plot(mod4, all.terms=TRUE)

#Below are three ways to select the models (see here to understand the difference: https://stats.stackexchange.com/questions/577/is-there-any-reason-to-prefer-the-aic-or-bic-over-the-other)
BIC(mod1,mod2,mod3,mod4) 
AIC(mod1,mod2,mod3,mod4)
AICc(mod1,mod2,mod3,mod4) #this is AIC adjusted for small sample size.

#Model 3 seems to be the best fit and most parsimonous
out <- broom::tidy(mod3)
kable(out)
write.csv(out,file="~/australian_plague_locust_model/Spatial_Heterogeneity/Model_fit/model_summaries/APL_modGI_summary.csv")
?tidy


dat$Pred_carb <- predict(mod3,type="response")[,1]
dat$Pred_prot <- predict(mod3,type="response")[,2]


population_mean <- dat %>%
  group_by(Population) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb),
            std.p = std(Pred_prot),std.c = std(Pred_carb))

diet_mean <- dat %>%
  group_by(Diet.Pair) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb),
            std.p = std(Pred_prot),std.c = std(Pred_carb))

sex_mean <- dat %>%
  group_by(Sex) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb),
            std.p = std(Pred_prot),std.c = std(Pred_carb))

library(MetBrewer)
MetBrewer::colorblind_palettes

pop_graph <- ggplot(dat,aes(x=Protein.consumed..g.,y=Carb.consumed..g.,fill=Population)) + 
  geom_point(shape=21) +  geom_point(data=population_mean, aes(x = mean.p, y = mean.c,fill=Population), size=5,shape=21) +
  coord_equal(ratio=1) + xlim(0,.055) + ylim(0,.055) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") + scale_fill_met_d(name="Demuth") +
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2)

diet_graph <- ggplot(dat,aes(x=Protein.consumed..g.,y=Carb.consumed..g.,fill=Diet.Pair)) + 
  geom_point(shape=21) +  geom_point(data=diet_mean, aes(x = mean.p, y = mean.c,fill=Diet.Pair), size=5,shape=21) +
  coord_equal(ratio=1) + xlim(0,.055) + ylim(0,.055) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") + scale_fill_manual(values=c("#4daf4a","#984ea3")) +
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0)

sex_graph <- ggplot(dat,aes(x=Protein.consumed..g.,y=Carb.consumed..g.,fill=Sex)) + 
  geom_point(shape=21) +  geom_point(data=sex_mean, aes(x = mean.p, y = mean.c,fill=Sex), size=5,shape=21) +
  coord_equal(ratio=1) + xlim(0,.055) + ylim(0,.055) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") + scale_fill_manual(values=c("#ff7f00","#ffff33")) +
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2)

layout <- "
AA
BC
"

fig1 <- pop_graph + SGR_plot + DevTime_plot + 
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')


Sup_fig1 <- diet_graph + sex_graph + plot_annotation(tag_levels = 'A')

ggsave(pop_graph,file="presentation/field_pop_it.png",dpi=600,height=3.95,width=3.95,units="in")


gam.check(mod3)

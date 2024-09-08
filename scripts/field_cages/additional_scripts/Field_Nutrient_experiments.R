#######
# Field Cage
#  Work
#   2015 Australia
#######

rm(list=ls())
library(tidyverse)
library(magrittr)
library(mgcv)
library(gratia)
library(ggpubr)
library(patchwork)
library(viridis)
library(DHARMa)
library(scales)

cage_dat <- read.csv("2015_Australia/data/raw/Field Cage/LocustCageData.csv") %>%
  mutate(Treatment = fct_recode(Treatment, "None" = "0N", "Medium" = "MedH", "High" = "HighN"))

cage_IT <- read.csv("2015_Australia/data/raw/Field Cage/Field cage IT.csv")

cage_over_time <- read.csv("2015_Australia/data/raw/Field Cage/Field_Cage_Nutrients_Through_Time.csv") %>%
  mutate(Treatment = fct_recode(Treatment, "None" = "0", "Medium" = "M", "High" = "H",
                                "None" = '1', "Medium" = "2", "High" = "3")) %>% as_tibble()


survival <- read.csv("2015_Australia/data/raw/Field Cage/LocustCageData_v02.xlsx - by cage.csv")
#formatting data
str(cage_dat)
cols <- c("CageNum", "Color", "Instar", "Sex","Treatment","FieldPlot")
cage_dat[cols] <- lapply(cage_dat[cols], factor)

names(cage_dat)

cage_dat <- cage_dat %>% drop_na("enterpogam....of.veg.") %>% dplyr::select(!c(24:28))

tab<-cage_dat %>% group_by(Treatment,CageNum) %>% tally()
view(tab)

flevels <-c("None","Medium","High")
cage_dat$Treatment <- factor(cage_dat$Treatment, levels=flevels)

cage_dat <- cage_dat %>% drop_na(Treatment)


names(cage_dat)

ggplot(cage_dat,aes(x=C.mg.mg,y=Protein.mg.mg)) + geom_point() + geom_smooth(method="lm") 

cor(cage_dat$C.mg.mg,cage_dat$Protein.mg.mg)

cage_dat$CageNum <- as.factor(cage_dat$CageNum )

names(cage_dat)
ggplot(cage_dat,aes(x=Treatment,y=Mass..g.,fill=Sex))+geom_boxplot()

mod <- gam(
  Mass..g. ~ Treatment + Sex + s(Protein.mg.mg,Carb.mg.mg) + s((CageNum),bs="re"),data=cage_dat,select=TRUE,family=scat()
)

plant_dat <- cage_dat %>%
  group_by(CageNum) %>%
  summarize(
    Carb.mg.mg_mean = mean(Carb.mg.mg),
    Protein.mg.mg_mean = mean(Protein.mg.mg),
    Treatment = Treatment
  ) %>% as.matrix() %>% as.data.frame() %>% drop_na(Treatment,Carb.mg.mg_mean)

plant_dat$Carb.mg.mg_mean <- as.numeric(plant_dat$Carb.mg.mg_mean)
plant_dat$Protein.mg.mg_mean <- as.numeric(plant_dat$Protein.mg.mg_mean)
plant_dat$CageNum <- as.factor(plant_dat$CageNum)
plant_dat$Treatment <- factor(plant_dat$Treatment, levels=flevels)
str(plant_dat)

mod_plant <- gam(list(
  Carb.mg.mg_mean ~ Treatment + s(CageNum,bs="re") ,
  Protein.mg.mg_mean ~ Treatment+ s(CageNum,bs="re")),
  family=mvn(d=2),select=TRUE, data=plant_dat
)

draw(mod_plant)
plot(mod_plant)
summary(mod_plant)

plant_dat$Pred_carb <- predict(mod_plant,type="response")[,1]
plant_dat$Pred_prot <- predict(mod_plant,type="response")[,2]

draw(mod)
appraise(mod)
summary(mod)
concurvity(mod,full=FALSE)


cage_dat$Pred <- predict(mod,type="response",newdata=cage_dat)


ggplot(cage_dat,aes(x=Protein.mg.mg,y=Carb.mg.mg,z=Pred)) + stat_summary_hex(bins=50) + scale_fill_viridis()

aovmod<-aov(Pred ~ Treatment,data=cage_dat)
TukeyHSD(aovmod)

treatment_mean_plant <- plant_dat %>%
  group_by(Treatment) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb))

treatment_mean_GHweight <- cage_dat %>%
  group_by(Treatment) %>%
  summarise(mean.p = mean(Pred,na.rm=TRUE))

plant_graph <- ggplot(plant_dat,aes(x=Protein.mg.mg_mean,y=Carb.mg.mg_mean,fill=Treatment)) + geom_point(shape=21,color="black")  + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(data=treatment_mean_plant, aes(x = mean.p, y = mean.c,fill=Treatment), size=5,shape=21,color="black") +
  coord_equal(ratio=1) + xlim(0,.25) + ylim(0,.25) + theme_pubr() + ylab("Carbohydrate content (mg/mg)") +
  xlab("Protein content (mg/mg)") + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) + ggtitle("Plant Nutrients by Treatment")


cage_dat <- cage_dat %>% filter(Sex %in% c("F","M"))

 GH_weight_graph <- 
  ggplot(cage_dat,aes(x=Treatment,y=Mass..g.,color=Sex)) + geom_jitter(width = .2) +
  scale_color_manual(values=c("#1b9e77","#d95f02")) +
  geom_point(data=treatment_mean_GHweight, aes(x = Treatment, y = mean.p), size=5,shape=21,color="black",fill="black") +
  theme_pubr() + ylab("Mass (g)") +
  xlab("Treatment") + ggtitle("Final grasshopper weight") +
   theme(text = element_text(size=15))


#working on surivial/adult proportions dat

str(survival)

cols <- c("cageNum", "N.level")
survival[cols] <- lapply(survival[cols], factor)

survival <- survival %>% mutate(
  Treatment = case_when(
    N.level == "175" ~ "High",
    N.level == "87.5" ~ "Medium",
    N.level == "0" ~ "None"
  )
)

flevels <-c("None","Medium","High")
survival$Treatment <- factor(survival$Treatment, levels=flevels)

Survival_graph <- ggplot(survival,aes(x=Treatment,y=SurviveProportion)) +
  geom_boxplot() + geom_jitter(width=.2) + theme_pubr() + ggtitle("Survival proportion") +
  ylab("Proportion") + scale_y_continuous(limits=c(0,1)) +
  theme(text = element_text(size=15))


Adultprop_graph <- ggplot(survival,aes(x=Treatment,y=AdultProportion)) +
  geom_boxplot() + geom_jitter(width=.2) + theme_pubr() + ggtitle("Adult proportion") +
  ylab("Proportion") +
  theme(text = element_text(size=15))


Survival_mod <- gam(SurviveProportion ~ Treatment + s(Protein.mg.mg,Carb.mg.mg), data=survival,select=TRUE,family=betar())
Adult_mod <- gam(AdultProportion ~ Treatment + s(Protein.mg.mg,Carb.mg.mg), data=survival,select=TRUE,family=betar())

summary(Survival_mod)
summary(Adult_mod)
appraise(Survival_mod)
appraise(Adult_mod)



flevels 

plant_graph + GH_weight_graph + Survival_graph + Adultprop_graph + plot_annotation(tag_levels = 'A')


plant_dat$Treatment

plant_graph_none <- ggplot((plant_dat),aes(x=Protein.mg.mg_mean,y=Carb.mg.mg_mean,fill=Treatment))  + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  coord_equal(ratio=1) + xlim(0,.25) + ylim(0,.25) + theme_pubr() + ylab("Carbohydrate content (mg/mg)") +
  xlab("Protein content (mg/mg)") + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) + ggtitle("Plant Nutrients by Treatment")



plant_graph_0N <- ggplot((plant_dat %>% filter(Treatment == "0N")),aes(x=Protein.mg.mg_mean,y=Carb.mg.mg_mean,fill=Treatment)) + geom_point(shape=21,color="black")  + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(data=(treatment_mean_plant%>% filter(Treatment == "0N")), aes(x = mean.p, y = mean.c,fill=Treatment), size=5,shape=21,color="black") +
  coord_equal(ratio=1) + xlim(0,.25) + ylim(0,.25) + theme_pubr() + ylab("Carbohydrate content (mg/mg)") +
  xlab("Protein content (mg/mg)") + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) + ggtitle("Plant Nutrients by Treatment")

plant_graph_0N_MedH <- ggplot((plant_dat %>% filter(Treatment != "HighN")),aes(x=Protein.mg.mg_mean,y=Carb.mg.mg_mean,fill=Treatment)) + geom_point(shape=21,color="black")  + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(data=(treatment_mean_plant  %>% filter(Treatment != "HighN")), aes(x = mean.p, y = mean.c,fill=Treatment), size=5,shape=21,color="black") +
  coord_equal(ratio=1) + xlim(0,.25) + ylim(0,.25) + theme_pubr() + ylab("Carbohydrate content (mg/mg)") +
  xlab("Protein content (mg/mg)") + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) + ggtitle("Plant Nutrients by Treatment")

plant_graph <- ggplot(plant_dat,aes(x=Protein.mg.mg_mean,y=Carb.mg.mg_mean,fill=Treatment)) + geom_point(shape=21,color="black")  + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(data=treatment_mean_plant, aes(x = mean.p, y = mean.c,fill=Treatment), size=5,shape=21,color="black") +
  coord_equal(ratio=1) + xlim(0,.25) + ylim(0,.25) + theme_pubr() + ylab("Carbohydrate content (mg/mg)") +
  xlab("Protein content (mg/mg)") + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) + ggtitle("Plant Nutrients by Treatment")



## Graphing plant nutrients over time

head(cage_over_time)
str(cage_over_time)
names(cage_over_time)

cols <- c("X", "Plot", "Cage","Treatment","Species","Position")
cage_over_time[cols] <- lapply(cage_over_time[cols], factor)

flevels <- c("None","Medium","High")

cage_over_time <- cage_over_time %>% select(1:6,8:14)

#cage_over_time$Treatment <- substring(cage_over_time$Plot, 1,1)
cage_over_time$Block <- substring(cage_over_time$Position, 1,1)
cage_over_time$Treatment<- factor(cage_over_time$Treatment ,levels=flevels)

cage_over_time$Block <- factor(cage_over_time$Block)

# By treatment

carbon <- ggplot(cage_over_time,aes(x=X,y=C.mg.mg, fill=Treatment)) + geom_boxplot(outlier.size = -1) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(position=position_jitterdodge(jitter.width=0.1),aes(fill=Treatment), color="black",pch=21) +
  ylab("Carbon (mg/mg)") + xlab("Date") + theme_pubr() + ggtitle("Carbon content by treatment and time")

Nitrogen <- ggplot(cage_over_time,aes(x=X,y=N.mg.mg, fill=Treatment)) + geom_boxplot(outlier.size = -1) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(position=position_jitterdodge(jitter.width=0.1),aes(fill=Treatment), color="black",pch=21) +
  ylab("Nitrogen (mg/mg)") + xlab("Date") + theme_pubr() + ggtitle("Nitrogen content by treatment and time")

Protein <- ggplot(cage_over_time,aes(x=X,y=Protein.mg.mg, fill=Treatment)) + geom_boxplot(outlier.size = -1) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(position=position_jitterdodge(jitter.width=0.1),aes(fill=Treatment), color="black",pch=21) +
  ylab("Protein (mg/mg)") + xlab("Date") + theme_pubr() + ggtitle("Protein content by treatment and time")

Carbohydrates <- ggplot(cage_over_time,aes(x=X,y=Carb.mg.mg, fill=Treatment)) + geom_boxplot(outlier.size = -1) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(position=position_jitterdodge(jitter.width=0.1),aes(fill=Treatment), color="black",pch=21) +
  ylab("Carbohydrate (mg/mg)") + xlab("Date") + theme_pubr() + ggtitle("Carbohydrate content by treatment and time")

Prot_carb_over_time <- ggplot(cage_over_time,aes(x=Protein.mg.mg,y=Carb.mg.mg, shape=X, color=Treatment)) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point() +
  ylab("Carbohydrate (mg/mg)") + xlab("Date") + theme_pubr()

carbon + Nitrogen + Protein + Carbohydrates

str(cage_over_time)


mod_plant_overtime_CarbProt <- gam(list(
  Carb.mg.mg ~ Treatment * X + s(Species,bs="re") + s(Block,bs="re") + s(Cage,bs="re") ,
  Protein.mg.mg ~ Treatment * X + s(Species,bs="re") + s(Block,bs="re") + s(Cage,bs="re")),
  family=mvn(d=2),select=TRUE, data=cage_over_time
)

str(cage_over_time)
summary(mod_plant_overtime_CarbProt)
gam.check(mod_plant_overtime_CarbProt)
plot(mod_plant_overtime_CarbProt, all.terms = T)
?broom::tidy

paramteric <- broom::tidy(mod_plant_overtime_CarbProt,parametric = TRUE)
nonpara <- broom::tidy(mod_plant_overtime_CarbProt,parametric = FALSE)

#write.csv(paramteric,"2015_Australia/output/field_parametric_results.csv")
#write.csv(nonpara,"2015_Australia/output/field_nonparametric_results.csv")

knitr::kable(out)



mod_plant_overtime_CarbonNitro <- gam(list(
  C.mg.mg ~ Treatment * X + s(Species,bs="re") + s(Block,bs="re") + s(Cage,bs="re") ,
  N.mg.mg ~ Treatment * X +  s(Species,bs="re") + s(Block,bs="re") + s(Cage,bs="re")),
  family=mvn(d=2),select=TRUE, data=cage_over_time)

str(cage_over_time)
summary(mod_plant_overtime_CarbonNitro)
gam.check(mod_plant_overtime_CarbonNitro)
plot(mod_plant_overtime_CarbProt, all.terms = T)


cage_over_time$Pred_carb <- predict(mod_plant_overtime_CarbProt,newdata=cage_over_time,type="response")[,1]
cage_over_time$Pred_prot <- predict(mod_plant_overtime_CarbProt,newdata=cage_over_time,type="response")[,2]

cage_over_time$Pred_carbon <- predict(mod_plant_overtime_CarbonNitro,newdata=cage_over_time,type="response")[,1]
cage_over_time$Pred_nitro <- predict(mod_plant_overtime_CarbonNitro,newdata=cage_over_time,type="response")[,2]


treatment_mean_Treatment_time <- cage_over_time %>%
  group_by(Treatment,X) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb)) %>% drop_na()

treatment_mean_treatment<- cage_over_time %>%
  group_by(Treatment) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb)) %>% drop_na()

treatment_mean_Treatment_time_C_N <- cage_over_time %>%
  group_by(Treatment,X) %>%
  summarise(mean.carbon = mean(Pred_carbon),mean.nitro = mean(Pred_nitro)) %>% drop_na()

treatment_mean_treatment <- cage_over_time %>%
  group_by(Treatment) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb)) %>% drop_na()


Prot_carb_over_timeA <- ggplot(cage_over_time,aes(x=Protein.mg.mg,y=Carb.mg.mg,  color=Treatment)) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point() +
  geom_point(data=(treatment_mean_treatment), aes(x = mean.p, y = mean.c,fill=Treatment,), size=5) +
  coord_equal(ratio=1) + xlim(0,.3) + ylim(0,.3) + 
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) +
  ylab("Carbohydrate (mg/mg)") + xlab("Protein (mg/mg)") + theme_pubr() + ggtitle("Nutrients by fertilization")

Prot_carb_over_time1 <- ggplot((cage_over_time %>% filter(X == "11/11/2015")),aes(x=Protein.mg.mg,y=Carb.mg.mg, fill=Treatment)) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(color="black",pch=21) +
  geom_point(data=(treatment_mean_Treatment_time %>% filter(X == "11/11/2015")), aes(x = mean.p, y = mean.c,fill=Treatment), size=5, color="black",pch=21) +
  #facet_wrap(~X) +
  coord_equal(ratio=1) + xlim(0,.3) + ylim(0,.3) + 
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) +
  ylab("Carbohydrate (mg/mg)") + xlab("Protein (mg/mg)") + theme_pubr() + ggtitle("11/11/2015")+
  theme(legend.title = element_blank())+
  theme(text = element_text(size=15))


Prot_carb_over_time2 <- ggplot((cage_over_time %>% filter(X == "11/25/2015")),aes(x=Protein.mg.mg,y=Carb.mg.mg, fill=Treatment)) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(color="black",pch=21) + ylab("") +
  geom_point(data=(treatment_mean_Treatment_time %>% filter(X == "11/25/2015")), aes(x = mean.p, y = mean.c,fill=Treatment), size=5, color="black",pch=21) +
  #facet_wrap(~X) +
  coord_equal(ratio=1) + xlim(0,.3) + ylim(0,.3) + 
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) +
  ylab("Carbohydrate (mg/mg)") + xlab("Protein (mg/mg)") + theme_pubr() + ggtitle("11/25/2015")+
  theme(legend.title = element_blank())+
  theme(text = element_text(size=15))


Prot_carb_over_time3 <- ggplot((cage_over_time %>% filter(X == "12/1/2015")),aes(x=Protein.mg.mg,y=Carb.mg.mg, fill=Treatment)) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(color="black",pch=21) +  ylab("") +
  geom_point(data=(treatment_mean_Treatment_time %>% filter(X == "12/1/2015")), aes(x = mean.p, y = mean.c,fill=Treatment), size=5, color="black",pch=21) +
  #facet_wrap(~X) +
  coord_equal(ratio=1) + xlim(0,.3) + ylim(0,.3) + 
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) +
  ylab("Carbohydrate (mg/mg)") + xlab("Protein (mg/mg)") + theme_pubr() + ggtitle("12/1/2015") +
  theme(legend.title = element_blank()) +
  theme(text = element_text(size=15))

?theme
Nitro_carb_over_time1 <- ggplot((cage_over_time %>% filter(X == "11/11/2015")),aes(x=N.mg.mg,y=C.mg.mg, fill=Treatment)) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(color="black",pch=21) +
  geom_point(data=(treatment_mean_Treatment_time_C_N %>% filter(X == "11/11/2015")), aes(x = mean.nitro, y = mean.carbon,fill=Treatment), size=5, color="black",pch=21) +
  #facet_wrap(~X) +
  coord_equal(ratio=1) + xlim(0,.5) + ylim(0,.5) + 
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) +
  ylab("Carbohydrate (mg/mg)") + xlab("Protein (mg/mg)") + theme_pubr() + ggtitle("11/11/2015")

Nitro_carb_over_time2 <- ggplot((cage_over_time %>% filter(X == "11/25/2015")),aes(x=N.mg.mg,y=C.mg.mg, fill=Treatment)) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(color="black",pch=21) +
  geom_point(data=(treatment_mean_Treatment_time_C_N %>% filter(X == "11/25/2015")), aes(x = mean.nitro, y = mean.carbon,fill=Treatment), size=5, color="black",pch=21) +
  #facet_wrap(~X) +
  coord_equal(ratio=1) + xlim(0,.5) + ylim(0,.5) + 
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) +
  ylab("Carbohydrate (mg/mg)") + xlab("Protein (mg/mg)") + theme_pubr() + ggtitle("11/25/2015")

Nitro_carb_over_time3 <- ggplot((cage_over_time %>% filter(X == "12/1/2015")),aes(x=Protein.mg.mg,y=Carb.mg.mg, fill=Treatment)) + 
  scale_color_manual(values=c("#addd8e","#31a354","#005a32")) +
  scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(color="black",pch=21) +
  geom_point(data=(treatment_mean_Treatment_time_C_N %>% filter(X == "12/1/2015")), aes(x = mean.nitro, y = mean.carbon,fill=Treatment), size=5, color="black",pch=21) +
  #facet_wrap(~X) +
  coord_equal(ratio=1) + xlim(0,.5) + ylim(0,.5) + 
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) +
  ylab("Carbohydrate (mg/mg)") + xlab("Protein (mg/mg)") + theme_pubr() + ggtitle("Plant nutrients - 12/1/2015")


Field_Trials_result <- ( Prot_carb_over_time1 | Prot_carb_over_time2 | Prot_carb_over_time3) / 
( GH_weight_graph | Survival_graph | Adultprop_graph) + plot_annotation(tag_levels = 'A')

ggsave(Field_Trials_result,file="Figures/Figure3.png",height=10,width=15,units="in", dpi=600)


# looking at Nitrogen early to see pattern

Nitrogen_M <- ggplot((cage_over_time),aes(x=X,y=N.mg.mg, fill=Species)) + geom_boxplot(outlier.size = -1) + 
  #scale_color_manual(values=c("#addd8e","#31a354","#005a32")) + scale_fill_manual(values=c("#addd8e","#31a354","#005a32")) +
  geom_point(position=position_jitterdodge(jitter.width=.5),aes(fill=Species), color="black",pch=21) +
  facet_wrap(~Treatment) +
  #scale_fill_viridis(discrete = TRUE) +
  ylab("Nitrogen (mg/mg)") + xlab("Date") + theme_pubr() + ggtitle("Nitrogen content by time and species split by treatment")



## protein by Nitrgon

PvN <- ggplot(cage_over_time,aes(x=N.mg.mg,y=Protein.mg.mgs))+geom_point() + geom_smooth(method="gam") + theme_pubr() + 
  annotate("text", x = 0.01, y = .25, label = "R^2 = 0.364")


cor(cage_over_time$N.mg.mg,cage_over_time$Protein.mg.mg)


Nitrogen_M / PvN

#######
# Field Cage
#  Work Intake Target
#   2015 Australia
#######

rm(list=ls())
library(tidyverse)
library(magrittr)
library(mgcv)
library(gratia)
library(ggpubr)
library(patchwork)
library(gganimate)
library(gifski)
library(MANOVA.RM)

std <- function(x) sd(x)/sqrt(length(x))

dat <- read.csv("2015_Australia/data/raw/Field Cage/Field_Cage_IT.csv", skip = 1)
str(dat)

cols <- c("LocustIDNum", "FieldPlot", "FieldCageNum", "Color","Sex")
dat[cols] <- lapply(dat[cols], factor)
names(dat)

dat$Sex[dat$Sex == "B"] <- "M"



names(dat)
dat <- dat %>%
  drop_na(Locusts.final.mass..g.) %>% mutate(
    TotalC_step1 = (HighCarb.initial.g. - HighCarb.final.g.) * 0.42,
    TotalP_step1 = (HighProt.initial.g. - HighProt.final.g.) * 0.42,
    TotalC_step2 = (HighCarb.initial.g..1 - HighCarb.final.g..1) * 0.42,
    TotalP_step2 = (HighProt.initial.g..1 - HighProt.final.g..1) * 0.42,
    TotalC_step3 = (HighCarb.initial.g..2 - HighCarb.final.g..2) * 0.42,
    TotalP_step3 = (HighProt.initial.g..2 - HighProt.final.g..2) * 0.42,
    TotalC_step4 = (HighCarb.initial.g..3 - HighCarb.final.g..3) * 0.42,
    TotalP_step4 = (HighProt.initial.g..3 - HighProt.final.g..3) * 0.42,
    TotalC_step5 = (HighCarb.initial.g..4 - HighCarb.final.g..4) * 0.42,
    TotalP_step5 = (HighProt.initial.g..4 - HighProt.final.g..4) * 0.42
  ) %>% mutate(Treatment = case_when(
    FieldCageNum == 1 ~ "High",
    FieldCageNum == 2 ~ "High",
    FieldCageNum == 3 ~ "High",
    FieldCageNum == 4 ~ "High",
    FieldCageNum == 21 ~ "High",
    FieldCageNum == 22 ~ "High",
    FieldCageNum == 23 ~ "High",
    FieldCageNum == 24 ~ "High",
    FieldCageNum == 29 ~ "High",
    FieldCageNum == 30 ~ "High",
    FieldCageNum == 31 ~ "High",
    FieldCageNum == 32 ~ "High",
    FieldCageNum == 9 ~ "Med",
    FieldCageNum == 10 ~ "Med",
    FieldCageNum == 11 ~ "Med",
    FieldCageNum == 12 ~ "Med",
    FieldCageNum == 17 ~ "Med",
    FieldCageNum == 18 ~ "Med",
    FieldCageNum == 19 ~ "Med",
    FieldCageNum == 20 ~ "Med",
    FieldCageNum == 25 ~ "Med",
    FieldCageNum == 26 ~ "Med",
    FieldCageNum == 27 ~ "Med",
    FieldCageNum == 28 ~ "Med",
    FieldCageNum == 5 ~ "none",
    FieldCageNum == 6 ~ "none",
    FieldCageNum == 7 ~ "none",
    FieldCageNum == 8 ~ "none",
    FieldCageNum == 13 ~ "none",
    FieldCageNum == 14 ~ "none",
    FieldCageNum == 15 ~ "none",
    FieldCageNum == 16 ~ "none",
    FieldCageNum == 33 ~ "none",
    FieldCageNum == 34 ~ "none",
    FieldCageNum == 35 ~ "none",
    FieldCageNum == 36 ~ "none")) %>% filter(
      TotalC_step1 >= 0,TotalP_step1 >= 0,
      TotalC_step2 >= 0,TotalP_step2 >= 0,
      TotalC_step3 >= 0,TotalP_step3 >= 0,
      TotalC_step4 >= 0,TotalP_step4 >= 0,
      TotalC_step5 >= 0,TotalP_step5 >= 0
      ) %>% mutate(ID = row_number())

dat2 <- dat %>% 
  tidyr::pivot_longer(cols = c(starts_with("TotalC_"), starts_with("TotalP_")),
                      names_to = c(".value", "interval"), 
                      names_pattern = "(.*)_(.*)") %>% dplyr::select(2:3,6:8,54:58) %>%
  mutate(interval = recode(interval, step1 = "1")) %>% mutate(interval = recode(interval, step2 = "2")) %>%
  mutate(interval = recode(interval, step3 = "3")) %>% mutate(interval = recode(interval, step4 = "4")) %>%
  mutate(interval = recode(interval, step5 = "5")) 
str(dat2)
dat2$ID <- as.factor(dat2$ID)
dat2$Treatment <- as.factor(dat2$Treatment)
dat2$interval <- as.factor(dat2$interval)


plot <- ggplot(dat2,aes(x=(TotalP),y=(TotalC),color=Treatment)) + geom_point() +
  coord_equal(ratio=1) + xlim(0,.03) + ylim(0,.03)  +
  xlab("Protein consumed (g)") + ylab("Carbohydrates consumed (g)") + theme_pubr()+
  geom_abline(slope = 1,linetype=1) + geom_abline(slope = 2,linetype=2) +
  facet_wrap(~interval) + scale_color_manual(values=c("#4daf4a","#984ea3"))

str(dat2)

mod1 <- gam(list( 
  TotalC ~   Sex + (interval) * Treatment + s(Locust.initialmass.g.),
  TotalP ~    Sex + (interval) * Treatment +  s(Locust.initialmass.g.)),
  family=mvn(d=2),select=TRUE, data=dat2
)

mod2 <- gam(list( 
  TotalC ~  Sex + (interval) * Treatment,   
  TotalP ~  Sex + (interval) * Treatment),
  family=mvn(d=2),select=TRUE, data=dat2
)

dat2$interval <- as.integer(dat2$interval)

mod2.5 <- gam(list( 
  TotalC ~  Sex + s(ID, bs="re") +  s((interval)) + Treatment,  #this is the correct way to specify repeated measures....but I dont have enough sample size....
  TotalP ~  Sex + s(ID, bs="re") +  s((interval)) + Treatment),
  family=mvn(d=2),select=TRUE, data=dat2
)

mod3 <- gam(list( 
  TotalC ~  Sex ,
  TotalP ~  Sex),
  family=mvn(d=2),select=TRUE, data=dat2
)

mod4 <- gam(list( 
  TotalC ~   1,
  TotalP ~   1),
  family=mvn(d=2),select=TRUE, data=dat2
)

mod5 <- gam(list( 
  TotalC ~   Treatment,
  TotalP ~   Treatment),
  family=mvn(d=2),select=TRUE, data=dat2
)
summary(mod5)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)

plot(mod2,all.terms = TRUE)



#Trying to use repeated measures MANOVA


summary_dat <- dat2 %>% group_by(Treatment,interval,Sex) %>%
  tally() %>% pivot_wider(names_from = interval,values_from = n)

view(summary_dat)
data(EEG)
eeg <- spread(EEG, feature, resp)
fit <- multRM(cbind(brainrate, complexity) ~ sex * region, data = eeg,
              subject = "id", within = "region", iter = 200)
str(dat2)


dat2$interval <- as.factor(dat2$interval)
dat2$ID <- as.integer(dat2$ID)
dat2 <- as.data.frame(dat2)


fit <- multRM(cbind(TotalC,TotalP) ~ Sex * interval, data=dat2,
              subject="ID",within="interval",iter=200, resampling = "paramBS"
)

summary(fit)
simCI(fit, contrast = "pairwise", type = "Tukey")


AICc(mod1,mod2,mod3,mod4,mod5) #mod2, mod1, mod3, mod4
AIC(mod1,mod2,mod3,mod4,mod5) #mod1, mod2, mod3, mod4
BIC(mod1,mod2,mod3,mod4,mod5) #mod4, mod3, mod2, mod1

dat2$Pred_carb <- predict(mod2,type="response")[,1]
dat2$Pred_prot <- predict(mod2,type="response")[,2]

GHSPP_mean <- dat2 %>%
  group_by(interval) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb))

sex_mean <- dat2 %>%
  group_by(Sex,interval) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb))

treatment_mean <- dat2 %>%
  group_by(Treatment,interval) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb), ratio = mean.p / mean.c)


dat2 %>%
  filter(interval == "1") %>%
  group_by(Treatment) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb), ratio = mean.p / mean.c)



str(dat2)
consumption_1 <- ggplot( (dat2 %>% filter(interval =="1")),aes(x=TotalP,y=TotalC)) + 
  geom_point() +  geom_point(data=(GHSPP_mean%>% filter(interval =="1")), aes(x = mean.p, y = mean.c), size=5,shape=21,fill="blue") +
  coord_equal(ratio=1) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0)+ xlim(0,.04) + ylim(0,.04) +
  ggtitle("Day 1")

consumption_2 <- ggplot( (dat2 %>% filter(interval =="2")),aes(x=TotalP,y=TotalC)) + 
  geom_point() +  geom_point(data=(GHSPP_mean%>% filter(interval =="2")), aes(x = mean.p, y = mean.c), size=5,shape=21,fill="blue") +
  coord_equal(ratio=1) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0)+ xlim(0,.04) + ylim(0,.04) +
  ggtitle("Day 2")

consumption_4 <-ggplot( (dat2 %>% filter(interval =="3")),aes(x=TotalP,y=TotalC)) + 
  geom_point() +  geom_point(data=(GHSPP_mean %>% filter(interval =="3")), aes(x = mean.p, y = mean.c), size=5,shape=21,fill="blue") +
  coord_equal(ratio=1) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0) + xlim(0,.04) + ylim(0,.04) +
  ggtitle("Day 3-4")


consumption_6 <- ggplot( (dat2 %>% filter(interval =="4")),aes(x=TotalP,y=TotalC)) + 
  geom_point() +  geom_point(data=(GHSPP_mean%>% filter(interval =="4")), aes(x = mean.p, y = mean.c), size=5,shape=21,fill="blue") +
  coord_equal(ratio=1) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0) + xlim(0,.04) + ylim(0,.04) +
  ggtitle("Day 5-6")

consumption_9 <-ggplot( (dat2 %>% filter(interval =="5")),aes(x=TotalP,y=TotalC)) + 
  geom_point() +  geom_point(data=(GHSPP_mean%>% filter(interval =="5")), aes(x = mean.p, y = mean.c), size=5,shape=21,fill="blue") +
  coord_equal(ratio=1) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0) + xlim(0,.04) + ylim(0,.04) +
  ggtitle("Day 7-9")


### Ploting individual lines
names(dat3)
dat3 <- dat %>% mutate(
  day1_C = TotalC_step1,
  day1_P = TotalP_step1,
  day1.2_C = TotalC_step1 + TotalC_step2,
  day1.2_P = TotalP_step1 + TotalP_step2,
  day1.3_C = day1.2_C + TotalC_step3,
  day1.3_P = day1.2_P + TotalP_step3,
  day1.4_C = day1.3_C + TotalC_step4,
  day1.4_P = day1.3_P + TotalP_step4,
  day1.5_C = day1.4_C + TotalC_step5,
  day1.5_P = day1.4_P + TotalP_step5) %>%
  tidyr::pivot_longer(cols = c(starts_with("day1")), names_to = c("interval",".value"),names_sep="_")

dat3_summary <- dat3 %>% group_by(interval) %>% summarize(C = mean(C), P=mean(P), ID=ID)

dat3_summary_treatment <- dat3 %>% group_by(interval,Treatment) %>% summarize(C = mean(C), P=mean(P))
dat3_summary_treatment_id <- dat3 %>% group_by(interval,Treatment,ID) %>% summarize(C = mean(C), P=mean(P))

dat3$Treatment
names(dat3)
Individual_lines <- ggplot(dat3,aes(y=C,x=P,group=as.factor(ID),color=as.factor(ID))) + geom_point() + geom_line(alpha=.5) +
  viridis::scale_color_viridis(discrete=TRUE) +
  geom_line(dat3_summary,mapping=aes(y=C,x=P),color="black",linetype=3) + geom_point(dat3_summary,mapping=aes(y=C,x=P),color="black",size=5) +
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0) + xlim(0,.1) + ylim(0,.1) + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  theme_pubr() + theme(legend.position = "none")


blank <- ggplot(dat3_summary_treatment_id,aes(y=C,x=P,fill=Treatment,shape=Treatment,group=ID, color=as.factor(Treatment)))  +
  geom_abline(slope = 1,linetype=1,intercept=0) + 
  geom_abline(slope = 2,linetype=6,intercept=0) + 
  geom_abline(slope = .5,intercept=0, color="dark green")  +
  ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  theme_pubr() + 
  scale_x_continuous(limits = c(0, .1)) + 
  scale_y_continuous(limits = c(0, 0.1))+ guides(color = FALSE) +
  theme(text = element_text(size=20))  + theme(legend.position="none")



plot <- ggplot(dat3_summary_treatment_id,aes(y=C,x=P,fill=Treatment,shape=Treatment,group=ID, color=as.factor(Treatment))) + geom_point(alpha=.5) + geom_line(alpha=.5) + 
  geom_line(data=dat3_summary_treatment,mapping=aes(y=C,x=P, group=Treatment,color=as.factor(Treatment)),size=1.5,inherit.aes = FALSE,linetype=1) +
  geom_point(data=dat3_summary_treatment,mapping=aes(y=C,x=P,fill=as.factor(Treatment),shape=Treatment),inherit.aes = FALSE,size=5) +
  scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values=c("#31a354","#a1d99b")) +
  scale_color_manual(values=c("#31a354","#a1d99b")) +
  ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  geom_abline(slope = 1,linetype=1,intercept=0) + 
  geom_abline(slope = 2,linetype=6,intercept=0) + 
  geom_abline(slope = .5,intercept=0, color="dark green") + 
  scale_x_continuous(limits = c(0, .1)) + 
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_pubr() + guides(color = FALSE) +
  theme(text = element_text(size=20))  + theme(legend.position="none")

plot

ggsave(blank,file="presentation/consumption_time_blank.png",dpi=700,width=10,height=10,units="in")
ggsave(plot,file="presentation/consumption_time.png",dpi=700,width=10,height=10,units="in")

ggsave(plot,file="figures/consumption_time.png",dpi=700,width=5,height=5,units="in")


Treatment_lines <- ggplot(dat3_summary_treatment_id,aes(y=C,x=P,group=as.factor(ID))) + geom_point() + geom_line(alpha=.5) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(dat3_summary_treatment,mapping=aes(y=C,x=P,color=Treament),linetype=3) + 
  geom_point(dat3_summary_treatment,mapping=aes(y=C,x=P,shape=as.factor(Treatment)),color="black",size=5) +
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0) + xlim(0,.1) + ylim(0,.1) +
  ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  theme_pubr() + guides(color = FALSE) + guides(shape=guide_legend(title="Treatment"))




Treatment_lines


layout <- "
AB#
CDE
"

supfig <- consumption_1 + consumption_2 + consumption_4 + consumption_6 + consumption_9 + 
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')



consumption <- ggplot(dat2,aes(x=TotalP,y=TotalC)) + 
  coord_equal(ratio=1) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0) #+ xlim(0,.05) + ylim(0,.05) 

sex_consumption <- ggplot(dat2,aes(x=TotalP,y=TotalC,fill=Sex)) + 
  geom_point(shape=21) +  geom_point(data=sex_mean, aes(x = mean.p, y = mean.c,fill=Sex), size=5,shape=21) +
  coord_equal(ratio=1)  + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") + scale_fill_manual(values=c("#ff7f00","#ffff33"))  +
  facet_wrap(~interval)+
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0)  +
  theme(panel.spacing.x = unit(6, "mm")) #+ xlim(0,.05) + ylim(0,.05)

Treatment_consumption <- ggplot(dat2,aes(x=TotalP,y=TotalC,fill=Treatment)) + 
  geom_point(shape=21) +  geom_point(data=treatment_mean, aes(x = mean.p, y = mean.c,fill=Treatment), size=5,shape=21) +
  coord_equal(ratio=1)  + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") + scale_fill_manual(values=c("#4daf4a","#984ea3"))  +
  facet_wrap(~interval)+
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0)  +
  theme(panel.spacing.x = unit(6, "mm")) #+ xlim(0,.05) + ylim(0,.05)

consumption + sex_consumption + Treatment_consumption









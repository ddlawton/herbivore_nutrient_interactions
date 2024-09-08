#######
# OA and APL
#  2013
#   Diet work
######
rm(list=ls())
library(tidyverse) 
library(mgcv)
library(MuMIn)
library(multcomp)
library(gratia)
library(ggpubr)
library(patchwork)
std <- function(x) sd(x)/sqrt(length(x))

dat <- read.csv("data/raw/IT_data/FieldCaughtOaCt_TargetExp.xls - rawdata.csv")
names(dat)

dat2 <- dat %>%
  drop_na(GhMb_final.g.)

dat2 <- dat2 %>%
  mutate(
    High_P_day1_eaten =  HproteinMass.g._day1initial - HproteinMass.g._day1final,
    Low_P_day1_eaten =  LproteinMass.g._day1initial - LproteinMass.g._day1final,
    High_P_day2_eaten =  HproteinMass.g._day2initial - HproteinMass.g._day2final,
    Low_P_day2_eaten =  LproteinMass.g._day2initial - LproteinMass.g._day2final,
    High_P_day4_eaten =  HproteinMass.g._day4initial - HproteinMass.g._day4final,
    Low_P_day4_eaten =  LproteinMass.g._day4initial - LproteinMass.g._day4final
  ) %>% mutate(ID = row_number())
names(dat2)





dietA <- dat2 %>% filter(DietPair == 'A')
dietB <- dat2 %>% filter(DietPair == 'B')

dietA <- dietA %>% mutate(
  HighP_Carb_day1 = High_P_day1_eaten*0.16666666667,
  HighP_Prot_day1 = High_P_day1_eaten*0.83333333333,
  LowP_Carb_day1 = Low_P_day1_eaten*0.83333333333,
  LowP_Prot_day1 = Low_P_day1_eaten*0.16666666667,
  HighP_Carb_day2 = High_P_day2_eaten*0.16666666667,
  HighP_Prot_day2 = High_P_day2_eaten*0.83333333333,
  LowP_Carb_day2 = Low_P_day2_eaten*0.83333333333,
  LowP_Prot_day2 = Low_P_day2_eaten*0.16666666667,
  HighP_Carb_day4 = High_P_day4_eaten*0.16666666667,
  HighP_Prot_day4 = High_P_day4_eaten*0.83333333333,
  LowP_Carb_day4 = Low_P_day4_eaten*0.83333333333,
  LowP_Prot_day4 = Low_P_day4_eaten*0.16666666667,
  TotalC_day1 = HighP_Carb_day1 + LowP_Carb_day1,
  TotalP_day1 = HighP_Prot_day1 + LowP_Prot_day1,
  TotalC_day2 = HighP_Carb_day2 + LowP_Carb_day2,
  TotalP_day2 = HighP_Prot_day2 + LowP_Prot_day2,
  TotalC_day4 = HighP_Carb_day4 + LowP_Carb_day4,
  TotalP_day4 = HighP_Prot_day4 + LowP_Prot_day4
)



dietB <- dietB %>% mutate(
  HighP_Carb_day1 = High_P_day1_eaten*0.33333333333,
  HighP_Prot_day1 = High_P_day1_eaten*0.66666666666,
  LowP_Carb_day1 = Low_P_day1_eaten*0.83333333333,
  LowP_Prot_day1 = Low_P_day1_eaten*0.16666666667,
  HighP_Carb_day2 = High_P_day2_eaten*0.33333333333,
  HighP_Prot_day2 = High_P_day2_eaten*0.66666666666,
  LowP_Carb_day2 = Low_P_day2_eaten*0.83333333333,
  LowP_Prot_day2 = Low_P_day2_eaten*0.16666666667,
  HighP_Carb_day4 = High_P_day4_eaten*0.33333333333,
  HighP_Prot_day4 = High_P_day4_eaten*0.66666666666,
  LowP_Carb_day4 = Low_P_day4_eaten*0.83333333333,
  LowP_Prot_day4 = Low_P_day4_eaten*0.16666666667,
  TotalC_day1 = HighP_Carb_day1 + LowP_Carb_day1,
  TotalP_day1 = HighP_Prot_day1 + LowP_Prot_day1,
  TotalC_day2 = HighP_Carb_day2 + LowP_Carb_day2,
  TotalP_day2 = HighP_Prot_day2 + LowP_Prot_day2,
  TotalC_day4 = HighP_Carb_day4 + LowP_Carb_day4,
  TotalP_day4 = HighP_Prot_day4 + LowP_Prot_day4)


dat <- rbind(dietA,dietB)  %>% dplyr::select(1:5,12,28,36,49:54)
names(dat)

dat2 <- dat %>% 
  tidyr::pivot_longer(cols = c(starts_with("TotalC_"), starts_with("TotalP_")),
                                      names_to = c(".value", "interval"), 
                                      names_pattern = "(.*)_(.*)") %>%
  mutate(interval = recode(interval, day1 = "1")) %>% mutate(interval = recode(interval, day2 = "2")) %>%
  mutate(interval = recode(interval, day4 = "4")) %>%
  filter(TotalC > 0) %>% filter(TotalP > 0)


str(dat2)

m <- gam(y ~ s(V1) + V2 + s(time) + s(time, id, bs = 'fs'),
         family=gaussian, data=dat, method = "REML")



cols <- c("CageNum", "GhNum", "GhSpp", "Sex","DietPair","ID","interval")
dat2[cols] <- lapply(dat2[cols], factor)
dat2$interval <- as.numeric(dat2$interval)
dat2$GhSpp <- as.factor(dat2$GhSpp)
dat2$Sex <- as.factor(dat2$Sex)
dat2$DietPair <- as.factor(dat2$DietPair)
ggplot(dat2,aes(x=TotalC,y=TotalP,color=interval,shape=GhSpp)) + geom_point()+
  coord_equal(ratio=1) + xlim(0,.22) + ylim(0,.22) + theme_pubr()
str(dat2)




mod2 <- gam(list( 
  TotalC ~  DietPair + Sex + s(GhMb.g.,k=15) + s(interval, ID, bs = 'fs') + s(GhSpp,bs="re"),
  TotalP ~  DietPair + Sex + s(GhMb.g.,k=15) + s(interval, ID, bs = 'fs') + s(GhSpp,bs="re")),
  family=mvn(d=2),select=TRUE, data=dat2
)

mod3 <- gam(list( 
  TotalC ~  DietPair + Sex + s(interval, ID, bs = 'fs') + s(GhSpp,bs="re"),
  TotalP ~  DietPair + Sex + s(interval, ID, bs = 'fs') + s(GhSpp,bs="re")),
  family=mvn(d=2),select=TRUE, data=dat2
)

mod4 <- gam(list( 
  TotalC ~  Sex + s(interval, ID, bs = 'fs') + s(GhSpp,bs="re"),
  TotalP ~  Sex + s(interval, ID, bs = 'fs') + s(GhSpp,bs="re")),
  family=mvn(d=2),select=TRUE, data=dat2
)

mod5 <- gam(list( 
  TotalC ~   s(interval, ID, bs = 'fs') + s(GhSpp,bs="re") + s(GhSpp,bs="re"),
  TotalP ~   s(interval, ID, bs = 'fs') + s(GhSpp,bs="re") + s(GhSpp,bs="re")),
  family=mvn(d=2),select=TRUE, data=dat2
)

AICc(mod2,mod3,mod4,mod5)
AIC(mod2,mod3,mod4,mod5)
BIC(mod2,mod3,mod4,mod5) #model 4 it is


gam.check(mod1)

summary(mod1)
summary(mod2)
summary(mod4) 



str(CT)
str(dat2)
plot(mod4,all.terms = TRUE)



dat2$Pred_prot <- predict(mod4,type="response")[,2]
dat2$Pred_carb <- predict(mod4,type="response")[,1]

GHSPP_mean <- dat2 %>%
  group_by(GhSpp,interval) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb),
            std.p = std(Pred_prot),std.c = std(Pred_carb))

diet_mean <- dat2 %>%
  group_by(GhSpp,DietPair,interval) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb),
            std.p = std(Pred_prot),std.c = std(Pred_carb))

sex_mean <- dat2 %>%
  group_by(GhSpp,Sex,interval) %>%
  summarise(mean.p = mean(Pred_prot),mean.c = mean(Pred_carb),
            std.p = std(Pred_prot),std.c = std(Pred_carb))

Species_consumption <- ggplot(dat2,aes(x=TotalP,y=TotalC)) + 
  geom_point() +  geom_point(data=GHSPP_mean, aes(x = mean.p, y = mean.c), size=5,shape=21,fill="blue") +
  coord_equal(ratio=1) + xlim(0,.2) + ylim(0,.2) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") +
  facet_wrap(~GhSpp+interval)+
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0)

Species_sex_consumption <- ggplot(dat2,aes(x=TotalP,y=TotalC,fill=Sex)) + 
  geom_point(shape=21) +  geom_point(data=sex_mean, aes(x = mean.p, y = mean.c,fill=Sex), size=5,shape=21) +
  coord_equal(ratio=1) + xlim(0,.2) + ylim(0,.2) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") + scale_fill_manual(values=c("#ff7f00","#ffff33"))  +
  facet_wrap(~GhSpp+interval)+
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0)  +
  theme(panel.spacing.x = unit(6, "mm"))

Species_diet_pair_consumption <- ggplot(dat2,aes(x=TotalP,y=TotalC,fill=DietPair)) + 
  geom_point(shape=21) +  geom_point(data=diet_mean, aes(x = mean.p, y = mean.c,fill=DietPair), size=5,shape=21) +
  coord_equal(ratio=1) + xlim(0,.2) + ylim(0,.2) + theme_pubr() + ylab("Carbohydrate consumed (g)") +
  xlab("Protein consumed (g)") + scale_fill_manual(values=c("#4daf4a","#984ea3")) +
  facet_wrap(~GhSpp+interval)+
  geom_abline(slope = 1,linetype=1,intercept=0) + geom_abline(slope = 2,linetype=2,intercept=0)

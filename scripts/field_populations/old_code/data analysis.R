##############################
#    Australia
#          2k15
##############################
rm(list=ls())

library(tidyverse)
library(mgcv)
library(piecewiseSEM)
library(nlme)
library(lme4)
library(lavaan)
library(piecewiseSEM)
library(gratia)

dat <- read.csv('data/raw/Survey Data/complete_survey_2015.csv')
correct_cover <- read.csv('data/raw/Survey Data/correct_veg_ground_cover.csv')
str(correct_cover)


summarized_cover <-
  correct_cover %>%
  group_by(Site,Transect) %>%
  summarize(dung.median = median(dung),rock.median = median(rock),vegetation.median = median(vegetation),
            grass.avg = (mean(grass)/100), forb.avg = (mean(forb)/100),shrub.avg = (mean(shrubs)/100))%>%
  unite(Site_ID,"Site","Transect",sep="")


final_dat <- merge(dat,summarized_cover,by="Site_ID")
str(final_dat)

cols <- c("dung.median", "rock.median", "vegetation.median")
final_dat[cols] <- lapply(final_dat[cols], factor)  

str(as.integer(final_dat$dung.median))


final_dat2 <-
  final_dat %>%
  mutate(
    dung.precentage = case_when(dung.median == "0" ~ 0,
                                dung.median == "1" ~0.03,
                                dung.median == "2" ~0.08,
                                dung.median == "3" ~0.18)
  )

final_dat %>%
  group_by(dung.median)
  summarize(length(.))



ggplot(final_dat,aes(x=TotalN,y=Grass_N)) + geom_smooth(method="gam")




mod <- gam((
  Grass_Protein ~ te(Grass_N,Grass_C)
), data=final_dat, select=TRUE)


k.check(mod)
draw(mod)
appraise(mod)
summary(mod)

summary(final_dat$dung.median)
final_dat$grass_protein_proportion <- 
  final_dat$Grass_Protein / (final_dat$Grass_Carb+final_dat$Grass_Protein)
final_dat$forb_protein_proportion <- 
  final_dat$Forb_Protein / (final_dat$Forb_Carb+final_dat$Forb_Protein)
final_dat$shrub_protein_proportion <- 
  final_dat$Shrub_Protein / (final_dat$Shrub_Carb+final_dat$Shrub_Protein)


binary_dat <- final_dat
final_dat$X2015_CT
binary_dat$X2015_CT[binary_dat$X2015_CT == 1] <- 0
binary_dat$X2015_CT[binary_dat$X2015_CT == 2] <- 0
binary_dat$X2015_CT[binary_dat$X2015_CT == 3] <- 0
binary_dat$X2015_CT[binary_dat$X2015_CT == 4] <- 0
binary_dat$X2015_CT[binary_dat$X2015_CT == 5] <- 1
binary_dat$X2015_CT[binary_dat$X2015_CT == 6] <- 1
binary_dat$X2015_CT[binary_dat$X2015_CT == 7] <- 1

str(as.factor(binary_dat$X2015_CT))

names(binary_dat)
binary_mod <- glm(X2015_CT~(grass_protein_proportion)+
                    Av_Grass+Av_Forb, family=binomial(),data=binary_dat)

summary(binary_mod)

plot(binary_mod)

ggplot(final_dat,aes(x=grass_protein_proportion,y=X2015_CT))+geom_point()+geom_smooth(method="gam",formula = y~s(x))
str(as.factor(final_dat$X2015_CT))
names(final_dat)
?gam
lm_mod <- lm(X2015_CT~(grass_protein_proportion)+(forb_protein_proportion)+
               Av_Grass+Av_Forb+dung.median+,dat=final_dat)

names(final_dat)
summary(lm_mod)

lmer_mod <- lmer(X2015_CT~(grass_protein_proportion)+(forb_protein_proportion)+
                   Av_Grass+Av_Forb+(1|Site_ID),dat=final_dat)


mod <- gam(X2015_CT~(grass_protein_proportion)+(forb_protein_proportion)+
              Av_Grass+Av_Forb+s(Lat,Long)+s(Site_ID,bs="re"),dat=final_dat)



summary(mod)
drop(mod)
plot(mod)
anova(mod)

manova_dat <-
  final_dat %>%
  drop_na(Grass_Protein,Grass_Carb,Grass_C,Grass_N,TotalN)

x <- cbind(final_dat$Grass_Protein,final_dat$Grass_Carb)
mod <- manova(x~Grass_C+Grass_N+TotalN,dat=final_dat)
summary(mod, test="Pillai")
summary.aov(mod)

ggplot(final_dat,aes(x=Grass_Protein,y=Grass_Carb,color=as.factor(Site_NO)))+geom_point()



length(manova_dat$Grass_C)


depedents <- c(grass$Carb..,grass$Protein..)
mod<-manova(depedents~Field+Location*Species*Year,data=grass)
mod <- manova(depedents~Field,data=grass)



?lme4


?glmer

dat_NA_removed <-
  final_dat %>%
  drop_na(X2015_CT,TotalN,Net_Nit,Net_Min,Av_Grass,Av_Forb,Av_Shrub,grass_protein_proportion,
         forb_protein_proportion,Grass_N,Forb_N,grass_protein_proportion,forb_protein_proportion, 
       dung.median,rock.median,grass.avg,forb.avg,shrub.avg) 
names(dat)

?drop_na 
  
str(dat)


cor(dat$bareground,dat$Av_Shrub)

names <- c('Bioreg' ,'Site_NO','X2015_CT','Norm_Dens','CT_count')
dat[names] <- lapply(dat[names] , factor)

ggplot(dat_NA_removed,aes(y=grass_protein_proportion,x=forb_protein_proportion)) + geom_point() +geom_smooth(method="lm")
str(dat_NA_removed)
qqnorm(log(dat$TotalN))
dat$X2015_CT

mod <- gam(as.numeric(X2015_CT) ~ (Grass_Carb),data=dat,family=ocat(R=7))
summary(mod)
?binomial

?gam
dat$X2015_CT

X2015_CT ~ Nutritiona_landscape + dung.median + rock.median + grass.avg + forb.avg + shrub.avg  
str(dat_NA_removed)

dat_NA_removed$Norm_Dens



?psem
library(psem)

names(dat_NA_removed)

fitted <- psem(
  lme(grass_protein_proportion~Grass_N+Grass_C,random=~1|Site_ID,data=dat_NA_removed),
  lme(Grass_N ~ TotalN,random=~1|Site_ID,data=dat_NA_removed)
)
coefs(fitted, standardize = "scale")

summary(fitted)
explore <-
  dat_NA_removed %>%
  select(X2015_CT,grass_protein_proportion,forb_protein_proportion,dung.median,rock.median,grass.avg,
         shrub.avg,Grass_N,Grass_C,TotalN,Forb_N,Forb_C)

library(CauseAndCorrelation)
Exploratory.path.analysis(explore)

fitted <- update(fitted, grass_protein_proportion %~~% forb_protein_proportion)
fitted <- update(fitted, Forb_N %~~% Grass_N)

summary(fitted)
?psem

str(dat)
summary(lme(forb_protein_proportion~Forb_N+Forb_C,random=~1|Site_ID,data=dat_NA_removed))



coefs(fitted, standardize = "scale", standardize.type = "latent.linear",
intercepts = FALSE)




mod <- glm((X2015_CT) ~ grass_protein_proportion + forb_protein_proportion +
       dung.median+rock.median+grass.avg+shrub.avg,family=gaussian(),data=dat_NA_removed)
summary(mod)
summary(dat_NA_removed$shrub.avg)


dat_NA_removed$sum_cov <- dat_NA_removed$grass.avg + dat_NA_removed$forb.avg + dat_NA_removed$shrub.avg

names


Exploratory.path.analysis(SEM_explore)
?DAG
vanishing.tetrads(SEM_explore)
summary(fitted)

library(lavaan)
?factor
Norm_Dens
dat_NA_removed$X2015_CT <- as.integer(dat_NA_removed$X2015_CT)
dat_NA_removed$Norm_Dens <- as.integer(dat_NA_removed$Norm_Dens)
dat_NA_removed$dung.median <- as.integer(dat_NA_removed$dung.median)
dat_NA_removed$rock.median <- as.integer(dat_NA_removed$rock.median)

fit.mod1<-"
#latent
nutritional_landscape =~ grass_protein_proportion + forb_protein_proportion

#Measured
X2015_CT ~ nutritional_landscape + dung.median + rock.median + grass.avg + shrub.avg

grass_protein_proportion ~ Grass_N + Grass_C

Grass_N ~ TotalN

forb_protein_proportion ~ Forb_N + Forb_C

Forb_N ~ TotalN
#covariances
TotalN ~~0*Forb_C
TotalN ~~0*Grass_C
Grass_C ~~0*Forb_C

"


mod1 <- sem(fit.mod1,data=dat_NA_removed,estimator="MLM",fixed.x=TRUE)
parTable(mod1)
summary(mod1)
mod3@SampleStats@cov
varTable(mod1)
var(dat_NA_removed$Norm_Dens)

coef(mod1)




fit.mod2<-"
#latent
nutritional_landscape =~ start(0.6)*grass_protein_proportion + start(0.3)*forb_protein_proportion

#Measured
X2015_CT ~ nutritional_landscape + dung.median + rock.median + start(.8)*grass.avg + start(.4)*shrub.avg

grass_protein_proportion ~ start(0.04)*Grass_N + start(0.8)*Grass_C

Grass_N ~ start(0.05)*TotalN

forb_protein_proportion ~ start(0.03)*Forb_N + start(0.7)*Forb_C

Forb_N ~ start(0.05)*TotalN
"


mod2 <- sem(fit.mod2,data=dat_NA_removed,estimator="WLSMV",test = "bootstrap")

summary(mod2)
summary(mod1, standardized=TRUE)
mod2@SampleStats@cov
varTable(mod2)



fit.mod3<-"
#latent
nutritional_landscape =~ start(0.6)*grass_protein_proportion + start(0.3)*forb_protein_proportion

#Measured
Norm_Dens ~ nutritional_landscape + dung.median + rock.median + start(.8)*grass.avg + start(.4)*shrub.avg

grass_protein_proportion ~ start(0.04)*Grass_N + start(0.8)*Grass_C

Grass_N ~ start(0.05)*TotalN

forb_protein_proportion ~ start(0.03)*Forb_N + start(0.7)*Forb_C

Forb_N ~ start(0.05)*TotalN
"

library(lavaan)
mod3 <- sem(fit.mod3,data=dat_NA_removed,estimator="WLSMV",test = "bootstrap")

summary(mod3)
summary(mod3, standardized=TRUE)
varTable(mod3)

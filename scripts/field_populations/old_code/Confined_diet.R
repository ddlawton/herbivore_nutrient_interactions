#######
# Field Population
# Confined
#  IT 
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
library(multcomp)
#Function
std <- function(x) sd(x)/sqrt(length(x))

dat <- read.csv("2015_Australia/data/raw/IT_data/confined_diet.csv") 

dat$Sex <- as.factor(dat$Sex)

flevels <- c("7p:35c","14p:28c","21p:21c","35p:7c")
dat$Diet <- factor(dat$Diet,levels = flevels)

dat <- dat %>%
  drop_na("Specific.growth.rate...ln.Mb.final.Mb.initial...time.")


names(dat)

ggplot(dat,aes(x=Diet,y=development.time.with.days.in.5th.instar.for.those.not.yet.molted..days. )) + geom_boxplot() + theme_pubr() +
  ylab("Specific growth rate")

names(dat)

SGR_mod <- gam(
  Specific.growth.rate...ln.Mb.final.Mb.initial...time. ~
    Diet + Sex + s(Locust.initial.weight..g.) , data=dat,select=TRUE, family=gaussian()
)

DevTime_mod <- gam(
  development.time.with.days.in.5th.instar.for.those.not.yet.molted..days. ~
    Diet + Sex + s(Locust.initial.weight..g.) , data=dat,select=TRUE, family=gaussian()
)



#Multiple comparisons
## this gets a little hairy as we have to code a work around to get a GAM to work with this function

#this basically creates a matrix to show what comparisons to test
contr <- matrix(0, nrow = 6, ncol = length(coef(SGR_mod))) #change 'mod2' to what ever model you want to look at
colnames(contr) <- names(coef(SGR_mod)) #change 'mod2' to what ever model you want to look at
rownames(contr) <- c("7p:35c - 14p:28c", "7p:35c - 21p:21c", "7p:35c - 35p:7c","14p:28c - 21p:21c","14p:28c - 35p:7c","21p:21c - 35p:7c") #these are the comparisons we are trying to test
contr[, 2:4] <- rbind(c(1, 0), c(0, 1), c(-1, 1))

multicomp <- glht(SGR_mod, linfct = contr) #This is the actual command that does the test. Again, change 'mod2' to what model you want
summary(multicomp)



summary(SGR_mod)
summary(DevTime_mod)
plot(SGR_mod,all.terms=T)
appraise(Time_mod)


summary(Time_mod)
plot(SGR_mod, all.terms = TRUE)

aov <- aov(Pred_SGR~Diet,data = dat)
sgr_mod_tidy <- tidy(TukeyHSD(aov))

aov <- aov(Pred_DevTime~Diet,data = dat)
DevTime_mod_tidy <- tidy(TukeyHSD(aov))


write.csv(sgr_mod_tidy,file="output/model_result_tables.csv")
write.csv(DevTime_mod_tidy,file="output/DevTime_model_result_tables.csv")

dat$Pred_SGR <- predict(SGR_mod,type="response")
dat$Pred_DevTime <- predict(DevTime_mod,type="response")

ggplot(dat,aes(x=Diet,y=Pred_SGR)) + geom_point()

flevel <- c("1p:5c","1p:2c","1p:1c","5p:1c")

dat_average_SGR <- dat %>%
  group_by(Diet) %>%
  summarise(mean = mean(Pred_SGR),std = std(Pred_SGR)) %>%
  mutate(Diet = case_when(
    Diet == "7p:35c" ~ "1p:5c",
    Diet == "14p:28c" ~ "1p:2c",
    Diet == "21p:21c" ~ "1p:1c",
    Diet == "35p:7c" ~ "5p:1c"), Diet = factor(Diet,levels=flevel)) 

dat_average_DevTime <- dat %>%
  group_by(Diet) %>%
  summarise(mean = mean(Pred_DevTime),std = std(Pred_DevTime)) %>%
  mutate(Diet = case_when(
    Diet == "7p:35c" ~ "1p:5c",
    Diet == "14p:28c" ~ "1p:2c",
    Diet == "21p:21c" ~ "1p:1c",
    Diet == "35p:7c" ~ "5p:1c"), Diet = factor(Diet,levels=flevel)) 

unique(levels(dat$Diet))

dat2 <- dat %>%
  mutate(Diet = case_when(
    Diet == "7p:35c" ~ "1p:5c",
    Diet == "14p:28c" ~ "1p:2c",
    Diet == "21p:21c" ~ "1p:1c",
    Diet == "35p:7c" ~ "5p:1c"), Diet = factor(Diet,levels=flevel)) 



SGR_plot <- 
  ggplot(dat2,aes(x=Diet,y=Pred_SGR)) +
  geom_point(aes(y=Specific.growth.rate...ln.Mb.final.Mb.initial...time.)) +
  geom_point(data=dat_average_SGR,aes(y=mean),shape=23,size=5,fill="#208cc0") +
  ylab("Specific Growth Rate") + xlab("Diet") + theme_pubr() +
  scale_y_continuous(limits=c(0,0.1))
  
SGR_plot
DevTime_plot <- ggplot(dat2,aes(x=Diet,y=Pred_DevTime)) +
  geom_point(aes(y=development.time.with.days.in.5th.instar.for.those.not.yet.molted..days.)) +
  geom_point(data=dat_average_DevTime,aes(y=mean),shape=22,size=5,fill="#f1af3a") +
  ylab("Development Time (days)") + xlab("Diet") + theme_pubr() 
DevTime_plot

ggsave(SGR_plot,file="presentation/SGR_plot.png",dpi=600,height=3.95,width=3.95,units="in")
ggsave(DevTime_plot,file="presentation/DevTime_plot.png",dpi=600,height=3.95,width=3.95,units="in")



ggsave()


SGR_plot +DevTime_plot
Lines_graph <- ggplot(data=dat) + 
  geom_point(aes(x=))


rm(list=ls())
library(tidyverse)
library(magrittr)
library(mgcv)
library(gratia)
library(ggpubr)
library(patchwork)
library(viridis)

cage_dat <- read.csv("data/raw/Field Cage/Field_Cage_Nutrients_Through_Time.csv")

str(cage_dat)
cage_dat$Species <- as.factor(cage_dat$Species)
ggplot(cage_dat,aes(x=N.mg.mg,y=Protein.mg.mg, color=Species)) + geom_point() + theme_pubr() + geom_smooth(method="gam")


cage_dat %>% group_by(Species) %>%
  summarize(cor = cor(N.mg.mg,Protein.mg.mg), method = c("spearman"))
?cor

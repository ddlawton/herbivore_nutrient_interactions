####
# Summarizing ground cover stuff
#
####
rm(list=ls())

library(tidyverse)
library(janitor)

dat <- read_csv("2015_Australia/data/raw/Field Cage/FieldCagePlantSurveys.xlsx - data.csv") %>%
  clean_names() %>%
  mutate(across(everything(), na_if, "<5"),
         across(6:28,as.numeric),
         across(1:5, as.factor)),
         cage_id = as.character(cage_id),
         cage_id = factor(case_when(
           cage_id != "outside" ~ "inside",
           TRUE ~ cage_id
         )))

summary(dat)

names(dat)



outside <- dat %>% filter(cage_id == "outside")

inside <- dat %>% filter(cage_id != "outside") %>%
  remove_empty(which = "cols") %>%
  select(2,4:11) %>%
  mutate_if(is.numeric, ~replace_na(., 0))

    


summary(inside)


dat2 <- dat %>%
  #pivot_longer(cols=c(6:28)) %>%
  mutate(value = replace(value,is.na(value), 0)) %>%
  group_by(n_treatment_kg_n_per_ha,plot_replicate,cage_id,name) %>%
  summarize(mean = round(mean(value,na.rm=TRUE),2)) %>%
  pivot_wider(names_from = name,values_from = mean) %>%
  select(where(~ any(. != 0)))






dat2 <- dat %>%
  pivot_longer(cols=c(6:28)) %>%
  mutate(value = replace(value,is.na(value), 0)) %>%
  group_by(n_treatment_kg_n_per_ha,plot_replicate,cage_id,name) %>%
  summarize(mean = round(mean(value,na.rm=TRUE),2)) %>%
  pivot_wider(names_from = name,values_from = mean) %>%
  select(where(~ any(. != 0)))

write.csv(inside,file="2015_Australia/output/ground_cover_field_trials.csv")

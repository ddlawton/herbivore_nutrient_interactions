#######
# iButton analysis
#  
#  
#######
rm(list=ls())
library(tidyverse)
library(lubridate)

std <- function(x) sd(x)/sqrt(length(x))

lab1 <- as_tibble(read.csv("2015_Australia/data/raw/iButton_data/T_Lab1.csv", skip = 19)) %>%
  mutate(rep = "Lab1")

lab2 <- as_tibble(read.csv("2015_Australia/data/raw/iButton_data/T_Lab2.csv", skip = 19)) %>%
  mutate(rep = "Lab2")

lab3 <- as_tibble(read.csv("2015_Australia/data/raw/iButton_data/T_Lab3.csv", skip = 19)) %>%
  mutate(rep = "Lab3")

lab4 <- as_tibble(read.csv("2015_Australia/data/raw/iButton_data/T_Lab4.csv", skip = 19)) %>%
  mutate(rep = "Lab4")

lab5 <- as_tibble(read.csv("2015_Australia/data/raw/iButton_data/T_Lab5.csv", skip = 19)) %>%
  mutate(rep = "Lab5")

lab_temps <- lab1 %>% rbind(lab2,lab3,lab4,lab5) %>%
  select(Date.Time,Value,rep) %>%
  mutate(Date.Time = parse_date_time(Date.Time, "mdy hms p")) %>%
  filter(Date.Time > "2015-11-11" & Date.Time <"2015-11-29") %>%
  mutate(hour = hour(Date.Time),
         day.night = case_when(
               between(hour,6,19) ~ "day",
               TRUE ~ "night"
  ))

lab_temps %>% group_by(day.night) %>%
  summarize(temperature = mean(Value), temperature_std = std(Value))



lab_temps %>% ggplot(aes(x=day.night,y=Value)) + geom_boxplot()

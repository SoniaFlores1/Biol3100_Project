#Here is an example for how to clean raw data in a project
#This script cleans the penguins data 
library(tidyverse)
library(palmerpenguins)
library(janitor)


df<- penguins_raw %>% clean_names()


df<- df %>% 
  select(species, island, date_egg, culmen_length_mm, culmen_depth_mm,
         flipper_length_mm, body_mass_g, sex)


df$species %>% unique() %>% 
  str_split(" ") %>% 
  map_chr(1)

clean<- df %>% 
  mutate(species = species %>% 
           str_split(" ") %>% 
           map_chr(1),
         year = date_egg %>% 
           str_split("-") %>% 
           map_chr(1) %>% 
           as.numeric()) %>% 
  select(-date_egg) %>% 
  rename(bill_length_mm = culmen_length_mm,
         bill_depth_mm = culmen_depth_mm)
#Let's save the cleaned data in a separate file

write_csv(clean, "./Data/Cleaned/clean_penguin_data_example.csv")

# Step 2:
# Format raw data into functional groups
# Identify differences among clickers (Sam, Corinne and Jonny)
# est. March 2024


library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)


# Import raw compiled data exported from script 1
raw <- read.csv("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/Final Products/RAW_03-06-2024.csv")



# Exploratory visuals:





# Compare differences between Corinne, Sam and Jonny
d <- merge %>% group_by(class, clickr) %>% 
  summarise(totcount = count(count)) 


summarise(totsite = count(unique(site_year)))


summarise(meancount = totcount/count(site_year))





# Export raw data
write.csv(merge,"C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/Final Products/RAW_03-06-2024.csv")

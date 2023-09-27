# Load in all VPI data sets, QA/AC and merge
# Corinne Amir 
# Sept 2023

library(tidyverse)


# Import data
setwd("E:/csv exports/") # External hard drive with 2020-2022 raw data

filenames = list.files(pattern="\\.csv$")
raw <- do.call("rbind", sapply(filenames, read.csv, simplify = FALSE))


raw<- tibble::rownames_to_column(raw, "row_names")

colnames(raw) <- c("site", "id","class","count","percent")

raw$year <- str_sub(raw$site,5,8)
raw$site <- str_sub(raw$site,13,15)

raw[is.na(raw)] <- 0


# Remove classes that don't show up throughout entire dataset

a <- raw %>% group_by(class) %>% 
             summarise(totcount = sum(count)) %>%
             filter(totcount > 0)

raw <- raw %>% filter(class != "Unk_soft_coral" & class %in% a$class) %>%
               filter(class %in% a$class)




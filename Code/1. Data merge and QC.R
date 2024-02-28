# Load in all VPI data sets, QA/AC and merge
# Corinne Amir 
# Sept 2023

library(tidyverse)
library(dplyr)
library(reshape2)


# Import 2020-2022 data
# setwd("E:/csv exports/") # External hard drive with 2020-2022 raw data
setwd("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/2020-2022 VPI exports/")

filenames = list.files(pattern="\\.csv$")
raw <- do.call("rbind", sapply(filenames, read.csv, simplify = FALSE))


raw<- tibble::rownames_to_column(raw, "row_names")

colnames(raw) <- c("site", "id","class","count","percent")

raw$year <- str_sub(raw$site,5,8)
raw$site <- str_sub(raw$site,13,15)

raw[is.na(raw)] <- 0




# Import json files (2) sent in from Sam and merge with raw
setwd("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/Sam exports/")

filenames = list.files(pattern="\\.csv$")
sam <- do.call("rbind", sapply(filenames, read.csv, simplify = FALSE))

sam<- tibble::rownames_to_column(sam, "row_names")

sam$year <- str_sub(sam$row_names,5,8)
sam$site <- str_sub(sam$row_names,13,15)

sam <- sam %>% group_by(site) %>% 
  mutate(percent = 100*(count/sum(count))) %>%
  select(-row_names)

sam <- left_join(sam, distinct(raw[,c("id","class")]))

sam[is.na(sam)] <- 0

col_order <- colnames(raw)
sam <- sam[,col_order]

a <- rbind(raw, sam)





# Import 2019 data
setwd("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/")
raw19 <- read.csv("Jonny_2019.csv")

# raw19 <- read.csv("2019_raw_wide.csv") #import data
# key <- read.csv("VPI_key.csv")
# 
# raw19 <- melt(raw19, id = "Site") # Switch from wide to long format

# raw19 <- select(raw19, -variable) # Merge key and raw data
# key <- rename(key, "value"="Number")
# raw19 <- left_join(raw19,key)
# raw19 <- raw19 %>% group_by(Site,Name) %>% summarise(count = n())

# sitepts <- raw19 %>% group_by(Site) %>% summarise(tot = sum(count)) # Add rows to match other raw data 
# raw19 <- left_join(raw19,sitepts)
# raw19$percent <- raw19$count/raw19$tot
raw19$year <- "2019"
raw19 <- select(raw19,-X)
raw19$site <- as.character(raw19$site)
raw19$site <- str_sub(raw19$site,1,3)
colnames(raw19) <- c("site","class","count","percent","year")


# Import 2014-2018 data and format to merge with 2019-2022 dataset
setwd("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/")
raw14 <- read.csv("1408-1809_Point Count_RAW.csv")

raw14 <- select(raw14,-Grand.Total)
raw14 <- melt(raw14, id = c("Site","Year","TRT")) # Switch from wide to long format

colnames(raw14) <- c("site","year","TRT","class","count")

# sitepts <- raw14 %>% group_by(site,year) %>% summarise(tot = sum(count)) # Add rows to match other raw data 
# raw14 <- left_join(raw14,sitepts)
# raw14$percent <- raw14$count/raw14$tot
# raw14 <- select(raw14,-tot)
raw14$site <- as.character(raw14$site)
raw14$year <- as.character(raw14$year)
raw14$class <- as.character(raw14$class)

raw14 <- raw14 %>% mutate(class = replace(class, class == "A..acuminata", "A. acuminata"))
raw14 <- raw14 %>% mutate(class = replace(class, class == "M..capitata", "M. capitata"))
raw14 <- raw14 %>% mutate(class = replace(class, class == "P..damicornis", "P. damicornis"))
raw14 <- raw14 %>% mutate(class = replace(class, class == "P..lobata", "Porites massive"))
raw14 <- raw14 %>% mutate(class = replace(class, class == "Pocillopora.sp.", "Pocillopora spp."))
raw14 <- raw14 %>% mutate(class = replace(class, class == "Psammocora","Psammocora spp."))
raw14 <- raw14 %>% mutate(class = replace(class, class == "Monitpora.sp.","Montipora spp."))



# Merge and QC 2019 and 2020-2022 datasets

merge <- full_join(raw, raw19)


  # Take a look at the data
a <- merge %>% group_by(site,year,class) %>% 
             summarise(totcount = sum(count)) %>%
             filter(totcount > 0)

b <- merge %>% group_by(year) %>% 
  summarise(totcount = n_distinct(site))

  # Remove classes that don't show up throughout entire dataset
merge <- merge %>% filter(class != "Unk_soft_coral" & class != "Review" & class != "Brillo spike pad"
                          & class != "None" & class != "Sponge" & count != 0)

  # QC names to match with 2014-2018 data
merge <- merge %>% mutate(class = replace(class, class == "CCA_sick" | class == "CCA_rubble", "CCA")) 
merge <- merge %>% mutate(class = replace(class, class == "Turf_rubble","Turf"))
merge <- merge %>% mutate(class = replace(class, class == "Acro_branch","A. acuminata"))
merge <- merge %>% mutate(class = replace(class, class == "Montipora capitata branching" |
                                                 class == "Montipora capitata encrusting" |
                                                 class == "Montipora capitata foliose" ,"M. capitata"))
merge <- merge %>% mutate(class = replace(class, class == "Monti_crust" | class == "Montipora spp" | class == "Monti_foliose" |
                                                 class == "Montipora patula" ,"Montipora spp."))
merge <- merge %>% mutate(class = replace(class, class == "Pocillopora damicornis","P. damicornis"))   
merge <- merge %>% mutate(class = replace(class, class == "Pocillopora meandrina/verrucosa","Pocillopora spp.")) 
merge <- merge %>% mutate(class = replace(class, class == "Fine sediment","Sediment")) 
merge <- merge %>% mutate(class = replace(class, class == "Porites spp_massive","Porites massive"))
merge <- merge %>% mutate(class = replace(class, class == "Psammocora","Psammocora sp."))
merge <- merge %>% mutate(class = replace(class, class == "Lobophora variegata","Lobophora"))
merge <- merge %>% mutate(class = replace(class, class == "Pavona spp","Pavona spp."))
merge <- merge %>% mutate(class = replace(class, class == "Porites superfusa","P. superfusa"))



# Merge 2014-2018 with 2019-2022 dataset

merge1 <- full_join(raw14,merge)
merge1 <- select(merge1,-id)

# Take a look at the data
b <- merge1 %>% group_by(year) %>% 
  summarise(totcount = n_distinct(site))
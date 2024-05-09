# Step 1:
# Load in all VPI data sets, QA/AC and merge
# est. Sept 2023

library(tidyverse)
library(dplyr)
library(reshape2)

#### Import datasets ####
# Import 2020-2022 data
# setwd("E:/csv exports/") # External hard drive with 2020-2022 raw data
setwd("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/2020-2022 VPI exports/")

filenames = list.files(pattern="\\.csv$")
raw <- do.call("rbind", sapply(filenames, read.csv, simplify = FALSE))


raw<- tibble::rownames_to_column(raw, "row_names")

colnames(raw) <- c("site", "id","class","count","percent")

raw$year <- str_sub(raw$site,5,8)
raw$site <- str_sub(raw$site,13,15)

raw <- select(raw, -percent)

raw[is.na(raw)] <- 0

raw$clickr <- rep("CA",times = nrow(raw))



# Import files (8) sent in from Sam and merge with raw
setwd("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/Sam exports/")

filenames = list.files(pattern="\\.csv$")
sam <- do.call("rbind", sapply(filenames, read.csv, simplify = FALSE))

sam<- tibble::rownames_to_column(sam, "row_names")

sam$year <- str_sub(sam$row_names,5,8)
sam$site <- str_sub(sam$row_names,13,15)

sam <- sam %>% group_by(site) %>%
  # mutate(percent = 100*(count/sum(count))) %>%
  select(-row_names)

sam$clickr <- rep("SC",times = nrow(sam))

sam <- left_join(sam, distinct(raw[,c("id","class")]))

sam[is.na(sam)] <- 0

col_order <- colnames(raw)
sam <- sam[,col_order]

sam <- as.data.frame(sam)

sam <- select(sam, -id)
raw <- select(raw, -id)



# Import 2019 data
setwd("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/")
raw19 <- read.csv("Jonny_2019.csv")

raw19$year <- "2019"
raw19 <- select(raw19,-X)
raw19$site <- as.character(raw19$site)
raw19$site <- str_sub(raw19$site,1,3)
colnames(raw19) <- c("site","class","count","percent","year")
raw19 <- select(raw19, -percent)

raw19$clickr <- rep("JC",times = nrow(raw19))



# Import 2014-2018 data and format to merge with 2019-2022 dataset
setwd("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/")
raw14 <- read.csv("1408-1809_Point Count_RAW.csv")

raw14 <- select(raw14,-Grand.Total)
raw14 <- melt(raw14, id = c("Site","Year","TRT")) # Switch from wide to long format

colnames(raw14) <- c("site","year","TRT","class","count")

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

raw14$clickr <- rep("CA",times = nrow(raw14))


#### Merge datasets ####
# Merge 2019 and 2020-2022 datasets
merge <- full_join(raw, raw19)


# Merge Sam's data with 2014-2022 dataset
merge <- full_join(sam, merge)


# Merge 2014-2018 with 2019-2022 dataset
merge <- full_join(raw14,merge)




# Take a look at the data
a <- merge %>% group_by(site,year,class) %>% 
             summarise(totcount = sum(count)) %>%
             filter(totcount > 0)

b <- merge %>% group_by(year) %>% 
  summarise(totcount = n_distinct(site))

c <- merge %>% group_by(class)%>% 
  summarise(totcount = sum(count))



# Clean up dataset
  # Remove classes that don't show up throughout entire dataset
merge <- merge %>% filter(class != "Unk_soft_coral" & class != "Review" & class != "Brillo spike pad"
                          & class != "None" & class != "Sponge" & count != 0 & class != "Instrument" & count != 0)

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

# raw14 %>% group_by(TRT, site) %>% summarise(ncount=unique(site))
merge <- merge %>% mutate(TRT = case_when(site == 101 ~ "REM", site == 103 ~ "RAN", site == 105 ~ "BYSP",
                                          site == 107 ~ "XSP", site == 109 ~ "RAN", site == 111 ~ "REM",
                                          site == 113 ~ "XSP", site == 115 ~ "BYSP", site == 117 ~ "RAN",
                                          site == 119 ~ "CON", site == 121 ~ "CON", site == 123 ~ "XSP",
                                          site == 125 ~ "CON", site == 127 ~ "BYSP", site == 129 ~ "REM"))

merge$site_year <- paste0(merge$site,"_",merge$year)

sapply(merge,unique)



# Add percent cover
merge <- merge %>% group_by(site_year) %>%
            mutate(percent = 100*(count/sum(count))) 


#### Group classes into functional groups ####
sapply(merge, unique)


# Tier 1 funcitonal groups:
calcifying <- c("A. acuminata","Hydnophora","Leptastrea","Montipora.sp.","M. capitata","P. damicornis","Monti_branch","P. superfusa",
                "Porites spp_foliose","Montipora spp.","Psammocora sp.","Unk_Pavona_foliose","Pocillopora spp.","Porites massive",
                "Pavona spp.","Fungia","Psammocora spp.", "CCA","Peyssonnelia","Red_calc_mac","Dictyota","Galaxaura","Halimeda")

noncalcifying <- c("Caulerpa","Cyano","Turf","Lobophora","Dictyosphaeria","Green_macro","Red fleshy macro", "Corallimorph")

nonbiological<- c("Sediment","Limestone","Sand")


merge$tier1 <- ifelse(merge$class %in% calcifying, "calcifying","noncalcifying")
merge$tier1 <- ifelse(merge$class %in% nonbiological, "nonbiological", as.character(merge$tier1))


# Tier 2 functional groups:
coral <- c("A. acuminata","Hydnophora","Leptastrea","Montipora.sp.","M. capitata","P. damicornis","Monti_branch","P. superfusa","Fungia","Pavona spp.",
           "Porites spp_foliose","Montipora spp.","Psammocora sp.","Unk_Pavona_foliose","Pocillopora spp.","Porites massive","Psammocora spp.")

calcifying.algae <- c("CCA","Peyssonnelia","Red_calc_mac","Dictyota","Galaxaura","Halimeda")

other.noncalcifying <- c("Caulerpa","Cyano","Turf","Lobophora","Dictyosphaeria","Green_macro","Red fleshy macro") 

merge$tier2 <- ifelse(merge$class %in% coral, "hard.coral","other.noncalcifying")
merge$tier2 <- ifelse(merge$class %in% calcifying.algae, "calcifying.algae",as.character(merge$tier2))
merge$tier2 <- ifelse(merge$class %in% nonbiological, "nonbiological",as.character(merge$tier2))
merge$tier2 <- ifelse(merge$class == "Corallimorph", "Corallimorph",as.character(merge$tier2))


# Tier 3 functional groups:
other.coral <- c("Hydnophora","Leptastrea","Montipora.sp.","Monti_branch","P. superfusa","Fungia","Pavona spp.","Porites spp_foliose",
                 "Montipora spp.","Psammocora sp.","Unk_Pavona_foliose","Pocillopora spp.","Porites massive","Psammocora spp.")

other.calcifying.algae <- c("Peyssonnelia","Red_calc_mac","Dictyota","Galaxaura","Halimeda")

merge$tier3 <- ifelse(merge$class == "A. acuminata", "A.acuminata", as.character(merge$tier3))
merge$tier3 <- ifelse(merge$class == "M. capitata", "M.capitata",as.character(merge$tier3))
merge$tier3 <- ifelse(merge$class == "P. damicornis", "P.damicornis",as.character(merge$tier3))
merge$tier3 <- ifelse(merge$class == "Corallimorph", "Corallimorph",as.character(merge$tier3))
merge$tier3 <- ifelse(merge$class %in% other.coral, "other.coral",as.character(merge$tier3))
merge$tier3 <- ifelse(merge$class == "CCA", "CCA",as.character(merge$tier3))
merge$tier3 <- ifelse(merge$class %in% other.calcifying.algae, "other.calcifying.algae",as.character(merge$tier3))
merge$tier3 <- ifelse(merge$class %in% other.noncalcifying, "other.noncalcifying",as.character(merge$tier3))
merge$tier3 <- ifelse(merge$class %in% nonbiological, "nonbiological",as.character(merge$tier3))



##### Export merge data ####
write.csv(merge,"C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/Final Products/RAW_03-11-2024.csv",row.names = F)
  
  
  
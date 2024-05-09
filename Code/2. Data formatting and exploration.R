# Step 2:
# Format raw data into functional groups
# Identify differences among clickers (Sam, Corinne and Jonny)
# est. March 2024


library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)


# Import raw compiled data exported from script 1
raw <- read.csv("C:/Users/Corinne.Amir/Documents/GitHub/Cmorph_Removal/CSV files/Final Products/RAW_03-11-2024.csv")
raw[is.na(raw)] <- 0


# Calculate average %cover and SE
avg <- raw %>% group_by(TRT, year, tier3) %>% summarise(avgpercent = mean(percent))
se <- raw %>% group_by(TRT, year, tier3) %>% summarise(se = sd(percent)/sqrt(length(percent)))
avg <- left_join(avg,se)


# Calculate final-initial % cover and SE
initial <- avg %>% filter(year == 1409) %>% rename(i_percent=avgpercent) %>% ungroup() %>% select(c(-year, -se))
final <- avg %>% filter(year == 2022) %>% rename(f_percent=avgpercent)  %>% ungroup() %>% select(c(-year, -se))
change <- full_join(initial,final)
change$i_percent[is.na(change$i_percent)] <- 0
change$f_percent[is.na(change$f_percent)] <- 0
change$change <- change$f_percent-change$i_percent



# Compare differences between Corinne, Sam and Jonny
d <- merge %>% group_by(class, clickr) %>% 
  summarise(totcount = count(count)) 


summarise(totsite = count(unique(site_year)))


summarise(meancount = totcount/count(site_year))




# Final percent covers:

ggplot(avg %>% filter(year == "2022"), aes(TRT, avgpercent), color = tier3, fill = tier3) +
  geom_bar(stat = "identity") + 
  facet_wrap(~tier3) +
  geom_errorbar(aes(ymin=avgpercent-se, ymax=avgpercent+se), width=.5,
                position=position_dodge(0.05), size = .5) 


# Final - Initial percent covers
  
ggplot(change, aes(TRT, change), color = tier3, fill = tier3) +
  geom_bar(stat = "identity") + 
  facet_wrap(~tier3)




#### ANOVAs ####
# Question: Does benthic substrate differ between treatments? Does it differ at different points in the study?
# Answer: One way anovas separated by time point for taxa of interest, permanova for overall benthic

# CCA
cca <- raw %>% filter(tier3 == "CCA" & year != 1408) 

#Assumptions
bartlett.test(percent ~ TRT, data=cca %>% filter(year == 2022)) #1409 only year that doesn't pass (only three plots with CCA)
tapply(sqrt(cca0$percent), list(cca0$TRT), shapiro.test) # 1409: error, 1509: RAN = 0.04, 1609: CON = 0.05, 1809: XSP = 0.04, 2019: XSP = 0.01, REM = 0.02, CON = 0.04 (only XSP sig with sqrt transform)



# Turf

# 


#Test interactions using a 2-way anova?

#cca 1409-2022
cca <- raw %>% filter(tier3 == "CCA" & year != 1408 & year != 1409 & year != 1509) 

#Assumptions
bartlett.test(percent ~ interaction(TRT, as.factor(year)), data=cca) #Homogenous without transform
tapply(cca$percent, list(cca$TRT), shapiro.test) # REM = 0.04, CON = 0.005, BYSP= 0.02 (worse off with sqrt transofrm)
tapply(sqrt(cca$percent), list(cca$year), shapiro.test) # 2021 = 0.05, 1910 = 0.0001 (2021 ok with sqrt transform)

#two-way anova (not really what we're interested in)
summary(aov(CCA~as.factor(Year)*TRT, data = cca)) #trt and year are significant but interaction is not
summary(aov(CCA~as.factor(Year), data = cca)) # f = 7.848, df = 5, p < 0.001
summary(aov(CCA~TRT, data = cca)) # f = 5.649, df = 4, p < 0.001
TukeyHSD(aov(CCA~as.factor(Year), data = cca)) 
# 1709 > 1510, 1606, 1609. p = 0.003, 0.0009, 0.023
# 1910 > 1510, 1606, 1609, 1809. p = 0.0007, 0.0002, 0.007, 0.041
TukeyHSD(aov(CCA~TRT, data = cca)) 
# REM > BYSP, RAN, XSP. p = 0.018, 0.005, 0.0003

#Holmes p adjustment
#Year
pval <- c(0.9992248,0.9857520,0.0029375,0.7963619,0.0007060,0.9114170,0.0008566,0.5791235,0.0001895,0.0229511,0.9894632,0.0066518,0.1126887,0.9984678,0.0407925)
p.adjust(pval, method="holm", n=15)
#1709-1510 = 0.0353, 1910-1510 = 0.0099, 1709-1606 = 0.0111, 1910-1606 = 0.0028

#Treatment
pval<- c(0.8711207,0.9931063,0.0177827,0.7548821, 0.6377295,0.1881712,0.2065407,0.0048761,0.9386729,0.0003428)
p.adjust(pval, method="holm", n=10)
#XSP-REM = 0.0034, REM-RAN = 0.0439


#corallimorph 1510-1809
cmorph <- select(perm.prop.abr, c(Site, TRT, Year, Corallimorph))

#Not parametric so use kruskal-wallis
kruskal.test(cmorph$Corallimorph, cmorph$Year) # NS
kruskal.test(cmorph$Corallimorph, cmorph$TRT) # p = 0.011
dunn.test(cmorph$Corallimorph, cmorph$TRT, method = "holm")
#CON > RAN, REM, XSP, BYSP. RAN different than BYSP
#CON > XSP (p = 0.0463), CON > RAN (p = 0.0025) with holm correction 


#split dataframe by trt
rem <- filter(cmorph, TRT == "REM") # NOT parametric
ran <- filter(cmorph, TRT == "RAN") # NOT parametric
xsp <- filter(cmorph, TRT == "XSP") # NOT parametric
bysp <- filter(cmorph, TRT == "BYSP") # NOT parametric
con <- filter(cmorph, TRT == "CON") # NOT parametric

bartlett.test(Corallimorph ~ as.factor(Year), data=con) #Homogenous without transform
tapply(sqrt(con$Corallimorph), list(con$Year), shapiro.test) # normal

kruskal.test(rem$Corallimorph, rem$Year) # NS
kruskal.test(ran$Corallimorph, ran$Year) # NS
kruskal.test(xsp$Corallimorph, xsp$Year) # NS
kruskal.test(bysp$Corallimorph, bysp$Year) # NS
kruskal.test(con$Corallimorph, con$Year) # p = 0.014
dunn.test(con$Corallimorph, con$Year, method = "holm")
#1510 different than 1910
dunn.test(con$Corallimorph, con$Year, method = "none")
#1510 > 1709, 1809, 1910
#1606 > 1709, 1809, 1910


#split dataframe by year
one <- filter(cmorph, Year == "1510")
two <- filter(cmorph, Year == "1606")
three <-filter(cmorph, Year == "1609")
four <- filter(cmorph, Year == "1709")
five <- filter(cmorph, Year == "1809")


kruskal.test(four$Corallimorph, four$TRT) # NS
kruskal.test(five$Corallimorph, five$TRT) # NS
kruskal.test(one$Corallimorph, one$TRT) # NS

kruskal.test(two$Corallimorph, two$TRT) # p = 0.018
dunn.test(two$Corallimorph, two$TRT, method = "bonferroni") 
#CON has more than RAN
dunn.test(two$Corallimorph, two$TRT, method = "none") 
#CON has more than RAN, REM, XSP. BYSP and RAN are different

kruskal.test(three$Corallimorph, three$TRT) # p = 0.037
dunn.test(three$Corallimorph, three$TRT, method = "bonferroni") 
#CON has more than REM and XSP
dunn.test(three$Corallimorph, three$TRT, method = "none") 
#CON has more than RAN, REM, and XSP



#Turf 1510-1910
turf <- droplevels(select(perm.prop.abr, c(Site, TRT, Year,Turf)))
turf$sqrt <- sqrt(turf$Turf)

#Assumptions
bartlett.test(sqrt ~ interaction(TRT, as.factor(Year)), data=turf) #Homogenous without transform
tapply(turf$sqrt, list(turf$TRT), shapiro.test) # XSP = 0.03 w/o sqrt transform (sqrt transform makes it worse)
tapply(turf$Turf, list(turf$Year), shapiro.test) # 1809 = 0.03 w/o sqrt transform

#two-way anova
summary(aov(Turf~as.factor(Year)*TRT, data = turf)) #trt and year are significant but interaction is not
summary(aov(Turf~as.factor(Year), data = turf)) # f = 6.588, p = 0.00015
summary(aov(Turf~TRT, data = turf)) # F = 2.83 p = 0.0295 
TukeyHSD(aov(Turf~as.factor(Year), data = turf)) 
# 1709 < 1510, 1606, 1609, 1809. p = 0.002, 0.0002, 0.0003, 0.019
#	1910 < 1510, 1606, 1609 1809. p > 0.0001
TukeyHSD(aov(Turf~TRT, data = turf)) # NS

#Add holmes adjustment to tukeyHSD
#Year
pval <- c(0.9867392,0.9932033,0.0022311,0.9837151,0.0000002,0.9999992,0.0002160,0.7678792,0.0000000,0.0003005,0.8152028,0.0000000,0.0192208,0.1722568,0.0000031)
p.adjust(pval,"holm",n=15)
#1709 < 1510, 1606, 1609(p = 0.0201, 0.0024, 0.0030)
#1910 < 1510, 1606, 1609, 1809 (p < 0.001 for all)


# Total Coral
coral <- droplevels(select(perm.prop.abr, c(Site, TRT, Year,Coral)))
coral$sqrt <- sqrt(coral$Coral)
coral$cube <- (coral$Coral)^0.25

#Assumptions
bartlett.test(cube ~ interaction(TRT, as.factor(Year)), data=coral) #NOT Homogenous with any transform
tapply(coral$sqrt, list(coral$TRT), shapiro.test) # multiple NOT normal with sqrt transform
tapply(coral$sqrt, list(coral$Year), shapiro.test) # multiple NOT normal with no transform, sqrt transform. 1510 also not normal with sqrt transform

#Not parametric so use kruskal-wallis
kruskal.test(coral$Coral, coral$Year) # p = 0.0002
kruskal.test(coral$Coral, coral$TRT) # p = 0.0002
dunn.test(coral$Coral, coral$Year, method = "none") 
#1809 > 1510, 1606, 1609. 1709 > 1510
dunn.test(coral$Coral, coral$TRT, method = "none") 
# CON < RAN, XSP, BYSP. REM < RAN, BYSP, XSP


# #split dataframe by year
# one <- filter(coral, Year == "1510") # NOT parametric
# two <- filter(coral, Year == "1606") # NOT parametric
# three <-filter(coral, Year == "1609") # NOT parametric
# four <- filter(coral, Year == "1709") # parametric
# five <- filter(coral, Year == "1809") # parametric
# 
# kruskal.test(one$Coral, one$TRT) # p = 0.027
# dunn.test(one$Coral, one$TRT, method = "bonferroni") 
# #REM different than XSP p = 0.026
# dunn.test(one$Coral, one$TRT, method = "none") 
# #REM and is less than BYSP, XSP, RAN. CON less than BYSP, XSP
# 
# kruskal.test(two$Coral, two$TRT) # p = 0.024
# dunn.test(two$Coral, two$TRT, method = "bonferroni")
# #REM different than XSP p = 0.035
# dunn.test(two$Coral, two$TRT, method = "none")
# #REM different than XSP, RAN, BYSP. CON less than XSP, BYSP, RAN
# 
# kruskal.test(three$Coral, three$TRT) # p = 0.026
# dunn.test(three$Coral, three$TRT, method = "bonferroni")
# #None show sig difference
# dunn.test(three$Coral, three$TRT, method = "none")
# #REM less than BYSP, RAN, XSP. CON different than RAN and XSP
# 
# summary(aov(Coral~TRT, data = four)) #TRT is significant
# TukeyHSD(aov(Coral~TRT, data = four)) # CON and REM < BYSP, RAN, XSP. XSP > RAN
# 
# summary(aov(Coral~TRT, data = five)) #TRT is significant
# TukeyHSD(aov(Coral~TRT, data = five)) # CON and REM < BYSP, RAN, XSP. 
# 
# 
# #Test just transplant treatments
# five.transp <- droplevels(five[-c(4:6,10:12),])
# kruskal.test(five.transp$Coral, five.transp$TRT) #ns


# Which time points are sig different from each other (transplant treatments only)?
coralsp.trts <- coralsp.comb %>% filter(TRT != "CON") %>%  filter(TRT != "REM") %>% filter(Year != "1408") %>% filter(Year != "1409") %>% filter(Year != "1509") 
coralsp.trts <- droplevels(coralsp.trts)


#Total Coral
coralsp.trts$sqrt.total <- sqrt(coralsp.trts$total)
#Assumptions
bartlett.test(sqrt.total~interaction(Year,TRT), data=coralsp.trts) # HOMOGENOUS with sqrt
tapply(coralsp.trts$sqrt.total, list(coralsp.trts$TRT), shapiro.test) # Mostly normal with sqrt
tapply(coralsp.trts$total, list(coralsp.trts$Year), shapiro.test) # NORMAL
#ANOVA
summary(aov(sqrt.total ~ TRT*Year, data = coralsp.trts)) # TRT and year are significant = separate
summary(aov(sqrt.total ~ Year, data = coralsp.trts)) # SIGNIFICANT p < 0.001
summary(aov(sqrt.total ~ TRT, data = coralsp.trts)) # NOT SIGNIFICANT
#Post-hoc
TukeyHSD(aov(sqrt.total ~ Year, data = coralsp.trts)) #
pval <- c(0.4283253,0.3261938,0.0000000,0.0000000,0.0000000, 0.9999712, 0.0000136, 0.0000000,0.0000000, 0.0000255,0.0000000,0.0000000,0.0003934,0.0000000,0.0006821)
p.adjust(pval,"holm",n=15) #1510 < 1910, 1809, 1709 | 1606 < 1910, 1809, 1709 | 1609 < 1910, 1809, 1709 | 1709 < 1809, 1910 | 1809 < 1910


#Acropora
coralsp.trts$sqrt.acro <- sqrt(coralsp.trts$acro)
#Assumptions
bartlett.test(sqrt.acro ~ interaction(TRT, as.factor(Year)), data=coralsp.trts) #Homogenous with sqrt
tapply(coralsp.trts$sqrt.acro, list(coralsp.trts$TRT), shapiro.test) # NORMAL
tapply(coralsp.trts$sqrt.acro, list(coralsp.trts$Year), shapiro.test) # NORMAL (except 16-6 = 0.03 with sqrt)
#ANOVA
summary(aov(sqrt.acro~TRT*Year, data = coralsp.trts)) # treatment and year are significant
summary(aov(sqrt.acro~TRT, data = coralsp.trts)) #NS
summary(aov(sqrt.acro~Year, data = coralsp.trts)) # SIGNIFICANT
#Post hoc
TukeyHSD(aov(sqrt.acro~TRT, data = coralsp.trts)) # none significant
TukeyHSD(aov(sqrt.acro~Year, data = coralsp.trts)) # 
pval <- c(0.8584891,0.6908274,0.0000207,0.0000000,0.0000000, 0.9995828,0.0009506, 0.0000000,0.0000000,0.0025164, 0.0000001,0.0000000,0.0366962,0.0000005,0.0126845)
p.adjust(pval,"holm",n=15) # 1510 < 1910, 1809, 1709 | 1606 < 1910, 1809, 1709 | 1609 < 1910, 1809, 1709 | 1709 < 1910



#Pocillopora
coralsp.trts$sqrt.poc <- sqrt(coralsp.trts$poc)
#Assumptions
bartlett.test(sqrt.poc ~ interaction(TRT, as.factor(Year)), data=coralsp.trts) #Homogenous with sqrt
tapply(coralsp.trts$sqrt.poc, list(coralsp.trts$TRT), shapiro.test) # NORMAL with sqrt
tapply(coralsp.trts$sqrt.poc, list(coralsp.trts$Year), shapiro.test) # NORMAL with sqrt

#ANOVA
summary(aov(sqrt.poc~TRT*Year, data = coralsp.trts)) # year is significant
summary(aov(sqrt.poc~Year, data = coralsp.trts)) 
#Post hoc
TukeyHSD(aov(sqrt.poc~Year, data = coralsp.trts)) # 1809 and 1910 > 1510, 1606, 1609, 1709. 1910 > 1809
pval <- c(0.9449378,0.9399936, 0.4203703,0.0000105, 0.0000009,1.0000000,0.0766322,0.0000004,0.0000000,0.0732876,0.0000004, 0.0000000,0.0050768,0.0005995,0.9824572)
p.adjust(pval,"holm",n=15) # 1510 < 1910, 1809 | 1606 < 1910, 1809 | 1609 < 1910, 1809 | 1709 < 1910, 1809


#Montipora
coralsp.trts$mont.sqrt <- sqrt(coralsp.trts$mont)
#Assumptions
bartlett.test(mont.sqrt ~ interaction(TRT, as.factor(Year)), data=coralsp.trts) #Homogenous w/o sqrt
tapply(coralsp.trts$mont.sqrt, list(coralsp.trts$TRT), shapiro.test) # xsp = 0.01 with sqrt transform
tapply(coralsp.trts$mont.sqrt, list(coralsp.trts$Year), shapiro.test) # NORMAL with sqrt transform
#ANOVA
summary(aov(mont.sqrt~TRT*Year, data = coralsp.trts)) # only year is significant
summary(aov(mont.sqrt~Year, data = coralsp.trts)) # F = 6.41, p = 0.000123 DF = 4,48
#Post hoc
TukeyHSD(aov(mont.sqrt~Year, data = coralsp.trts)) # 1709, 1809, 1910 > 1510
pval <- c(0.0968979,0.2024925,0.0040160,0.0012431,0.0001039,0.9991973,0.8353553,0.6157011,0.1989730, 0.6313474,0.3937780,0.0949226,0.9988916,0.8619797,0.9729259)
p.adjust(pval,"holm",n=15) # 1510 < 1910, 1809  


# When did each treatment experience increases in coral cover?
#split dataframe by trt
rem <- filter(coral, TRT == "REM") # NOT parametric
ran <- filter(coral, TRT == "RAN") #  parametric
xsp <- filter(coral, TRT == "XSP") #  parametric
bysp <- filter(coral, TRT == "BYSP") # parametric
con <- filter(coral, TRT == "CON") # NOT parametric

summary(aov(Coral~as.factor(Year), data = ran)) #year is significant
TukeyHSD(aov(Coral~as.factor(Year), data = ran)) # 1910 > 1510, 1606, 1609, 1709, 1809. 1809 > 1709, 1609, 1606, 1510.

summary(aov(Coral~as.factor(Year), data = xsp)) #year is significant
TukeyHSD(aov(Coral~as.factor(Year), data = xsp)) # 1910 > 1510, 1606, 1609, 1709. 1809 > 1609, 1606, 1510.

summary(aov(Coral~as.factor(Year), data = bysp)) #year is significant
TukeyHSD(aov(Coral~as.factor(Year), data = bysp)) # 1910 > 1510, 1606, 1609, 1709, 1809. 1809 > 1709, 1609, 1606, 1510. 1709 > 1609, 1510

kruskal.test(con$Coral, con$Year) # ns

kruskal.test(rem$Coral, rem$Year) # p = 0.017
dunn.test(rem$Coral, rem$Year, method = "bonferroni")
#1809 > 1510
dunn.test(rem$Coral, rem$Year, method = "none")
#1910 > 1609, 1606, 1510 | 1809 > 1510, 1606, 1609 


#CONCLUSIONS:
#Coral increased in all treatments except for CON. Largest increases in coral cover were 2017 and 2018
#In 2017 turf was replaced by CCA, particularly in plots without coral transplants
#REM had sig greater CCA than transplant treatments and the most CCA occurred in 1709
#Cmorph decreased significantly in CON and did not have sig increase in other treatments. 
#By 2017 cmorph decreased levels similar to other treatments

#non-parametric: coral, cmorph
#parametric: turf, CCA


#### end ANOVAs ####


#### Percent Cover Line Graphs ####

ggplot(avg %>% filter(year != 1408 & year != 1606), aes(x = as.factor(year), y = avgpercent, group = TRT, color = TRT)) +
  geom_line(stat = "identity") +
  geom_point(size = 3) +
  facet_wrap(vars(tier3), scales = "free") +
  geom_errorbar(aes(ymin=avgpercent-se, ymax=avgpercent+se), width=.5,
                position=position_dodge(0.05), size = .5) 



# Line graph by coral species
coralsp <- avg %>% filter(tier3 == "A.acuminata"|tier3 == "P.damicornis" | tier3 == "M.capitata"| tier3 == "other.coral")

# coralsp <- coralsp.comb %>% gather(Group, Percent, acro:total) # switch from wide to long format
# 
# coralsp  <- coralsp  %>% # ungroup(Year) %>% # this step is sometimes necessary to change the year names
#   mutate(Year = factor(Year, levels = c(1408,1409,1509,1510,1606,1609,1709,1809,2019),labels = c("Aug. 2014","Sept. 2014",
#                                                                                                  "Sept. 2015","Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017", "Sept. 2018","Oct. 2019"))) 
# 
# coralsp$TRT <- as.factor(coralsp$TRT)
# levels(coralsp$TRT)[levels(coralsp$TRT)=="BYSP"] <- "SSP" # Fix TRT names
# levels(coralsp$TRT)[levels(coralsp$TRT)=="RAN"] <- "NAG"

# Order levels for legend
coralsp$TRT <- factor(coralsp$TRT, levels = c("XSP", "SSP", "NAG","REM","CON")) 

coralsp$Month <- ifelse(coralsp$Year == "Aug. 2014", -1.4, as.numeric(coralsp$Month)) # need numeric values for proper spacing
coralsp$Month <- ifelse(coralsp$Year == "Sept. 2014", 0.3, as.numeric(coralsp$Month))
coralsp$Month <- ifelse(coralsp$Year == "Sept. 2015", 11.5, as.numeric(coralsp$Month))
coralsp$Month <- ifelse(coralsp$Year == "Oct. 2015", 13.5, as.numeric(coralsp$Month))
coralsp$Month <- ifelse(coralsp$Year == "Jun. 2016", 21, as.numeric(coralsp$Month))
coralsp$Month <- ifelse(coralsp$Year == "Sept. 2016", 24, as.numeric(coralsp$Month))
coralsp$Month <- ifelse(coralsp$Year == "Sept. 2017", 36, as.numeric(coralsp$Month))
coralsp$Month <- ifelse(coralsp$Year == "Sept. 2018", 48, as.numeric(coralsp$Month))
coralsp$Month <- ifelse(coralsp$Year == "Oct. 2019", 61, as.numeric(coralsp$Month))

coralsp <- coralsp[!(coralsp$TRT=="CON" & coralsp$Year == "Sept. 2015"),] # Remove CON treatment from 1510

coral_timeseries_mean <- aggregate(Percent~TRT+Month+Group, coralsp, mean) # summarize to mean values
#coral_timeseries_mean <- aggregate(Percent~TRT+Year+Group, coralsp, mean) # month wasn't working


# Separate dataframe by species
acro <- coral_timeseries_mean %>% filter(Group == "acro")
mont <- coral_timeseries_mean %>% filter(Group == "mont")
poc <- coral_timeseries_mean %>% filter(Group == "poc")
total <- coral_timeseries_mean %>% filter(Group == "total")


#Create Line Graphs of each species 
acroplot <- ggplot(acro, aes(x=Month, y = Percent, group = TRT, color = TRT)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  ylim(0,27) +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=7,
                position=position_dodge(0.05), size = 1) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1,color="black", size = 10),
        axis.text.y = element_text(color = "black",size=10),
        text = element_text(family = "Times New Roman",color="black", size = 10),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "none") +
  scale_x_continuous(breaks=c(-1.43,0.3,11.5,13.49,21,24,36,48,61), labels = scales::number_format(accuracy = 1)) +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4')) +
  #geom_vline(xintercept = 12.5, linetype = "dashed", color = "black", size = .8) +
  geom_segment(aes(x = 12.5, y = -Inf, xend = 12.5, yend = 26.3), linetype = "dashed", col = "black", size = 1) +
  annotate(geom = "text", x = -2.7, y = 27, label = "b)", color = "black", fontface = "bold",size=4,family="Times New Roman") +
  annotate(geom = "text", x = 6.7, y = 27, label = "A. acuminata", color = "black", fontface = "bold.italic", size=4,family="Times New Roman")

montplot <- ggplot(mont, aes(x=Month, y = Percent, group = TRT, color = TRT)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  ylim(-.1,5.3) +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=7,
                position=position_dodge(0.0), size = 1) +
  #labs(y = "Percent", x = "Months Post-Transplantation") +
  theme(axis.text.x = element_text(color="black",angle = 25, hjust = 1,size = 10),
        axis.text.y = element_text(color = "black",size=10),
        text = element_text(family = "Times New Roman",color="black", size= 10),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "none") +
  scale_x_continuous(breaks=c(-1.43,0.3,11.5,13.49,21,24,36,48,61), labels = scales::number_format(accuracy = 1)) +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4')) +
  geom_segment(aes(x = 12.5, y = -Inf, xend = 12.5, yend = 5), linetype = "dashed", col = "black", size = 1) +
  annotate(geom = "text", x = -2.7, y = 5.3, label = "c)", color = "black", fontface = "bold",size=4,family="Times New Roman") +
  annotate(geom = "text", x = 5.3, y = 5.3, label = "M. capitata", color = "black", fontface = "bold.italic", size=4,family="Times New Roman")

pocplot <- ggplot(poc, aes(x=Month, y = Percent, group = TRT, color = TRT)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  ylim(-.1,5.3) +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=7,
                position=position_dodge(0.0), size = 1) +
  #labs(y = "Percent", x = "Months Post-Transplantation") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1,color="black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        text = element_text(family = "Times New Roman",color="black", size = 10),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "none") +
  scale_x_continuous(breaks=c(-1.43,0.3,11.5,13.49,21,24,36,48,61), labels = scales::number_format(accuracy = 1)) +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4')) +
  geom_segment(aes(x = 12.5, y = -Inf, xend = 12.5, yend = 5), linetype = "dashed", col = "black", size = 1) +
  annotate(geom = "text", x = -2.7, y = 5.3, label = "d)", color = "black", fontface = "bold",size=4,family="Times New Roman") +
  annotate(geom = "text", x = 6.8, y = 5.3, label = "P. damicornis", color = "black", fontface = "bold.italic", size=4,family="Times New Roman")

totalplot <- ggplot(total, aes(x=Month, y = Percent, group = TRT, color = TRT)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  ylim(-.1,27) +
  guides(fill=guide_legend(title="Key")) +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=7,
                position=position_dodge(0.0), size = 1) +
  labs(y = "Percent", x = "Months Post-Transplantation") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1,color="black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        text = element_text(family = "Times New Roman",color="black",size=10),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.title = element_text(face="bold",family="Times New Roman",size=10),
        legend.position = c(.068, 0.69),
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=11),
        legend.key.size = unit(.5, "cm"),
        legend.key.width = unit(0.25,"cm")) +
  scale_x_continuous(breaks=c(-1.43,0.3,11.5,13.49,21,24,36,48,61), labels = scales::number_format(accuracy = 1)) +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4')) +
  geom_segment(aes(x = 12.5, y = -Inf, xend = 12.5, yend = 26.3), linetype = "dashed", col = "black", size = 1) +
  # annotate(geom = "line", xmin=-Inf, xmax=12.5, ymin=-Inf, ymax=Inf, color = "grey",alpha = 0.7) +
  annotate(geom = "text",  x = -2.7, y = 27, label = "a)", color = "black", fontface = "bold",size=4, family="Times New Roman") +
  annotate(geom = "text", x = 6, y = 27, label = "Total Coral", color = "black", fontface = "bold", size=4,family="Times New Roman")


# Compile bar plots into one figure
compile <- plot_grid(totalplot, acroplot,montplot, pocplot, ncol = 2, scale = 1, align="v");compile


#Add Axis titles....still working on it 
axislabels <- ggplot(total, aes(x=Month, y = Percent, group = TRT, color = TRT)) + #current workaround 
  geom_line(size = 1) +
  geom_point()+
  ylim(-.1,18) +
  labs(y = "Percent Cover", x = "Months Post-Transplantation") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Times New Roman",color="black", face = "bold", size = 16))

ggdraw(add_sub(compile, "Months Post-Transplantation", vpadding=grid::unit(0,"lines"),y=2, x=1, vjust=4.5, angle=0))


#Export figures
ggsave(plot=compile,path="C:/Users/Corinne.Amir/Documents/Miscellaneous/Cmorph Removal Manuscript/Figures",
       width=9.1,height=8,filename=paste0("Coralsp Linegraph_11-18-2020.jpg"))

ggsave(plot=axislabels,path="C:/Users/Corinne.Amir/Documents/Miscellaneous/Cmorph Removal Manuscript/Figures",
       width=8.5,height=8,filename=paste0("Coralsp Linegraph_axis labels_11-14-2020.jpg"))


#### end ####


#### Stacked Community Bar Plots ####

#Calculate mean value for each functional group in each TRT x Year pairing
mn.perm.prop <- perm.prop.comb %>% group_by(Year,TRT) %>% mutate(mean(CCA)) %>% mutate(mean(Turf)) %>% 
  mutate(mean(Corallimorph)) %>% mutate(mean(Coral)) %>%
  mutate(mean(Fleshy.Algae)) %>%
  mutate(mean(Calc.Algae)) %>% mutate(mean(Nonbiological)) %>%
  slice(which(row_number() %% 3 == 1)) %>%
  select(TRT,Year, 'mean(CCA)','mean(Turf)','mean(Corallimorph)','mean(Coral)','mean(Fleshy.Algae)','mean(Calc.Algae)','mean(Nonbiological)')      


colnames(mn.perm.prop) <- c("TRT","Year" ,"CCA" , "Turf Algae" , "Corallimorph",'Nonbiological', "Coral", 'Non-Calcified Algae' ,'Other Calcified Algae')

t <- mn.perm.prop[,3:9]
rowSums(t) # OK

mn.perm.prop <- mn.perm.prop %>% gather(Group, Percent, CCA:`Other Calcified Algae`)


# Make bar graph of each treatment separately
#Prep data
dummyrows <- data.frame(TRT = "CON", Year = c("1408", "1509")) #, Group = "Corallimorph", Percent = 0) 
mn.perm.prop <- full_join(mn.perm.prop, dummyrows) # Make spacing correct

mn.perm.prop$Year <- as.vector(mn.perm.prop$Year)
mn.perm.prop <- mn.perm.prop %>% ungroup(Year) %>% # Fix x axis names
  mutate(Year = factor(Year, levels = c(1408,1409,1509,1510,1606,1609,1709,1809,1910),
                       labels = c("Aug. 2014", "Sept. 2014", "Sept. 2015", "Oct. 2015", "Jun. 2016", "Sept. 2016", "Sept. 2017", "Sept. 2018","Oct. 2019")))

bysp <- mn.perm.prop %>% filter(TRT == "BYSP")
con <- mn.perm.prop %>% filter(TRT == "CON")
ran <- mn.perm.prop %>% filter(TRT == "RAN")
rem <- mn.perm.prop %>% filter(TRT == "REM")
xsp <- mn.perm.prop %>% filter(TRT == "XSP")



#Plots
conplot <- ggplot(con, aes(x=Year, y=Percent, fill=factor(Group, levels = c("Nonbiological","Non-Calcified Algae","Other Calcified Algae","Coral","Corallimorph","CCA","Turf Algae")))) +
  geom_bar(stat='identity', position='fill',width = .88) + 
  theme_bw() +
  theme(legend.position = "none",aspect.ratio = .95, 
        axis.title.y = element_blank(), axis.title.x = element_blank(),axis.text.y = element_blank(),
        plot.margin = unit(c(1,0,0,0),"cm"),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(labels = c("0", "25", "50","75","100")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black",size = 13)) +
  scale_fill_manual(values=c('grey',"orange3",'yellow','dodgerblue1','coral4','hotpink1', 'chartreuse3')) +
  theme(axis.title = element_blank())

remplot <- ggplot(rem, aes(x=Year, y=Percent, fill=factor(Group, levels = c("Nonbiological","Non-Calcified Algae","Other Calcified Algae","Coral","Corallimorph","CCA","Turf Algae")))) +
  geom_bar(stat='identity', position='fill' ,width = .88) + 
  theme_bw()+
  theme(legend.position = "none", aspect.ratio = .95) +
  scale_y_continuous(labels = c("0", "25", "50","75","100")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 ,color="black",size = 13), 
        axis.title.y = element_blank(),#axis.title.y = element_text(color = "black", "face" = "bold"), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,0,0,0),"cm"),
        text = element_text(family = "Times New Roman")) +
  scale_fill_manual(values=c('grey',"orange3",'yellow','dodgerblue1','coral4','hotpink1', 'chartreuse3')) +
  theme(text=element_text(family="Times New Roman",color = "black",size = 13)) 

ranplot <- ggplot(ran, aes(x=Year, y=Percent, fill=factor(Group, levels = c("Nonbiological","Non-Calcified Algae","Other Calcified Algae","Coral","Corallimorph","CCA","Turf Algae")))) +
  geom_bar(stat='identity', position='fill', width = .88) + 
  theme_bw()+
  theme(legend.position = "none", aspect.ratio = .95,
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),
        plot.margin = unit(c(1,0,0,0),"cm")) + # top, right, bottom, left
  scale_y_continuous(labels = c("0", "25", "50","75","100")) +
  scale_fill_manual(values=c('grey',"orange3",'yellow','dodgerblue1','coral4','hotpink1', 'chartreuse3')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black",size = 13)) 

xspplot <- ggplot(xsp, aes(x=Year, y=Percent, fill=factor(Group, levels = c("Nonbiological","Non-Calcified Algae","Other Calcified Algae","Coral","Corallimorph","CCA","Turf Algae")))) +
  geom_bar(stat='identity', position='fill' ,width = .88) + 
  theme_bw()+
  theme(legend.position = "none", aspect.ratio = .95) +
  scale_y_continuous(labels = c("0", "25", "50","75","100")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 ,color="black", size = 13), 
        axis.title.y = element_blank(),# axis.title.y = element_text(color = "black", face = "bold"), 
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman"),
        axis.text.y = element_text(color = "black", size = 13),
        plot.margin = unit(c(1,0,0,0),"cm")) +
  scale_fill_manual(values=c('grey',"orange3",'yellow','dodgerblue1','coral4','hotpink1', 'chartreuse3')) +
  theme(text=element_text(family="Times New Roman",color = "black")) 


byspplot <- ggplot(bysp, aes(x=Year, y=Percent, fill=factor(Group, levels = c("Nonbiological","Non-Calcified Algae","Other Calcified Algae","Coral","Corallimorph","CCA","Turf Algae")))) +
  geom_bar(stat='identity', position='fill', width = .88) + 
  theme_bw()+
  theme(legend.position = "none", aspect.ratio = .95, 
        axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,0,0,0),"cm")) + # top, right, bottom, left
  scale_y_continuous(labels = c("0", "25", "50","75","100")) +
  scale_fill_manual(values=c('grey',"orange3",'yellow','dodgerblue1','coral4','hotpink1', 'chartreuse3')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black", size = 13)) 

#Make figure for legend (amd y-axos title)
legend<-ggplot(mn.perm.prop, aes(x=Year, y=Percent, fill=factor(Group, levels = c("Nonbiological","Non-Calcified Algae","Other Calcified Algae","Coral","Corallimorph","CCA","Turf Algae")))) +
  geom_bar(stat='identity', position='fill') + 
  scale_fill_manual(values=c('grey',"orange3",'yellow','dodgerblue1','coral4','hotpink1', 'chartreuse3')) +
  guides(fill=guide_legend(title="Key")) +
  theme(legend.title = element_text(face="bold",family="Times New Roman",size = 16),
        legend.text = element_text(size = 16),
        legend.position = c(.2, .5),
        text = element_text(family = "Times New Roman"),
        axis.title.y = element_text(size=12, vjust=1),
        plot.margin = unit(c(0,0,0,1),"cm"))

legend.add <- get_legend(legend) 

# Compile bar plots into one figure
compile <- plot_grid(xspplot,byspplot,ranplot,remplot,conplot, legend.add, ncol = 3, scale = 1
                     ,rel_widths = c(1.,1,1)
                     ,labels = c("a) XSP", "b) SSP", "c) NAG"," d) REM", "e) CON"), align = "v"
                     ,label_fontfamily = "Times New Roman", label_size = 14
                     ,label_x = c(.07,.02,-.00,.03,-.0)  
                     ,label_y = c(0.97,0.97,0.97,0.97,0.97,0.97)); compile


#ggdraw(add_sub(compile, "Label", vpadding=grid::unit(0,"lines"),y=2, x=-0.1, vjust=4.5, angle=90))
#   draw_plot(compile) +
#   draw_label("Percent", x=0.02, y=.55, size = 12, angle=90, fontface="bold", colour="black") 
# draw_label("Year", x=0.02, y=.55, size = 12, angle=90, fontface="bold", colour="black") 

ggsave(plot=compile,path="C:/Users/Corinne.Amir/Documents/Miscellaneous/Cmorph Removal Manuscript/Figures",
       width=11,height=8,filename=paste0("Community Barplot_11-03-2020.jpg"))

#### end ####

#### Grouped Community Bar Plots ####

#Calculate mean value for each functional group in each TRT x Year pairing
perm.prop.long <- perm.prop.comb %>% gather(Group, Percent, CCA:Calc.Alg)
mn.perm.prop <- perm.prop.long %>% group_by(Year,TRT,Group) %>% mutate(mean(Percent)) %>% select(-c(Site, Percent)) %>% distinct()


#Calculate SE
coral.st.err <- function(x) {
  sd(x)/sqrt(length(x))
}
SE_benthic <- aggregate(Percent~Year+TRT+Group,perm.prop.long,coral.st.err)


#Add SE values into data frame
colnames(SE_benthic) <- c("Year"   , "TRT", "Group"  , "SE")
colnames(mn.perm.prop) <- c("TRT","Year" ,"Group" , "Percent")
benthic_timeseries_mean <- left_join(mn.perm.prop,SE_benthic)
benthic_timeseries_mean <- benthic_timeseries_mean %>% mutate(Group = recode(Group, Turf = "Turf Algae", Calc.Algae = "Other Calcified Algae",Fleshy.Algae = "Non-Calcified Algae",))


#Prep data
dummyrows <- data.frame(TRT = "CON", Year = c("1408", "1509"), Percent = 0, SE = 0)
dummygroups <- data.frame(TRT = "CON",Group = c("Turf Algae", "Other Calcified Algae", "CCA", "Corallimorph", "Coral", "Nonbiological", "Non-Calcified Algae"))
dummyadd <- left_join(dummyrows,dummygroups)
benthic_timeseries_mean <- full_join(benthic_timeseries_mean, dummyadd) # Make spacing correct

benthic_timeseries_mean$Year <- as.vector(benthic_timeseries_mean$Year)
benthic_timeseries_mean <- benthic_timeseries_mean %>% ungroup(Year) %>% # Fix x axis names
  mutate(Year = factor(Year, levels = c(1408,1409,1509,1510,1606,1609,1709,1809,1910),
                       labels = c("Aug. 2014", "Sept. 2014", "Sept. 2015", "Oct. 2015", "Jun. 2016", "Sept. 2016", "Sept. 2017", "Sept. 2018","Oct. 2019")))

#Subset dataset
simper.groups <- c("Corallimorph", "Coral","CCA","Turf Algae")
bysp <- benthic_timeseries_mean %>% filter(TRT == "BYSP") %>% filter(Group %in% c("Corallimorph", "Coral","CCA","Turf Algae"))
con <- benthic_timeseries_mean %>% filter(TRT == "CON") %>% filter(Group %in% c("Corallimorph", "Coral","CCA","Turf Algae"))
ran <- benthic_timeseries_mean %>% filter(TRT == "RAN") %>% filter(Group %in% c("Corallimorph", "Coral","CCA","Turf Algae"))
rem <- benthic_timeseries_mean %>% filter(TRT == "REM") %>% filter(Group %in% c("Corallimorph", "Coral","CCA","Turf Algae"))
xsp <- benthic_timeseries_mean %>% filter(TRT == "XSP") %>% filter(Group %in% c("Corallimorph", "Coral","CCA","Turf Algae"))
transp1 <- benthic_timeseries_mean %>% filter(TRT != "CON") %>% filter(TRT != "REM") %>% group_by(Year,Group) %>% select(-TRT) %>%
  summarize(mn.percent = mean(Percent)) %>% filter(Group %in% c("Corallimorph", "Coral","CCA","Turf Algae"))
transp2 <- benthic_timeseries_mean %>% filter(TRT != "CON") %>% filter(TRT != "REM") %>% group_by(Year,Group) %>% select(-TRT) %>%
  summarize(mn.se = mean(SE)) %>% filter(Group %in% c("Corallimorph", "Coral","CCA","Turf Algae"))
transp <- left_join(transp1,transp2)

#Plot
by_group <- ggplot(benthic_timeseries_mean, aes(x = factor(Group, levels = c("Nonbiological","Non-Calcified Algae","Other Calcified Algae","Coral","Corallimorph","CCA","Turf Algae")), y = Percent, fill = Year)) + 
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black",size = 13)) +
  xlab("Benthic Functional Group") +
  theme(text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

con_grouped <-  ggplot(con, aes(x = reorder(Group, -Percent), y = Percent, fill = Year)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.85, color= "black") +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=0.05, position=position_dodge(0.85), size = 1) +
  ylim(0,100) +
  theme_bw() +
  theme(legend.position = "none", aspect.ratio = .95, 
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust = 1,color="black",size = 13),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

rem_grouped <-  ggplot(rem, aes(x = reorder(Group, -Percent), y = Percent, fill = Year)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.85, color= "black") +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=0.05, position=position_dodge(0.85), size = 1) +
  ylim(0,100) +
  theme_bw() +
  theme(legend.position = "none", aspect.ratio = .95, 
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1,color="black",size = 13),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ran_grouped <-  ggplot(ran, aes(x = reorder(Group, -Percent), y = Percent, fill = Year)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.85, color= "black") +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=0.05, position=position_dodge(0.85), size = 1) +
  theme_bw() +
  ylim(0,100) +
  theme(legend.position = "none", aspect.ratio = .95, 
        axis.text.x = element_text(angle = 45, hjust = 1,color="black",size = 13),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),) 

xsp_grouped <-  ggplot(xsp, aes(x = reorder(Group, -Percent), y = Percent, fill = Year)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.85, color= "black") +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=0.05, position=position_dodge(0.85), size = 1) +
  theme_bw() +
  ylim(0,100) +
  theme(legend.position = "none", aspect.ratio = .95, 
        axis.text.x = element_text(angle = 45, hjust = 1,color="black",size = 13),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

bysp_grouped <-  ggplot(bysp, aes(x = reorder(Group, -Percent), y = Percent, fill = Year)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.85, color= "black") +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=0.05, position=position_dodge(.85), size = 1) +
  theme_bw() +
  ylim(0,100) +
  theme(legend.position = "none", aspect.ratio = .95, 
        axis.text.x = element_text(angle = 45, hjust = 1,color="black",size = 13),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

transp_grouped <- ggplot(transp, aes(x = reorder(Group, -mn.percent), y = mn.percent, fill = Year)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.85, color= "black") +
  geom_errorbar(aes(ymin=mn.percent-mn.se, ymax=mn.percent+mn.se), width=0.05, position=position_dodge(.85), size = 1) +
  theme_bw() +
  ylim(0,100) +
  theme(legend.position = "none", aspect.ratio = .95, 
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1,color="black",size = 13),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

legend<-ggplot(bysp, aes(x = reorder(Group, -Percent), y = Percent, fill = Year)) +
  geom_bar(position = "dodge", stat = "identity", width = 1, color= "black") +
  # scale_fill_manual(values=c('grey',"orange3",'yellow','dodgerblue1','coral4','hotpink1', 'chartreuse3')) +
  guides(fill=guide_legend(title="Key")) +
  xlab("Benthic Functional Group") +
  ylab("Percent Cover")+
  theme(legend.title = element_text(face="bold",family="Times New Roman",size = 14),
        legend.text = element_text(size = 14),
        legend.position = c(.2, .52),
        axis.title.x=element_text(size=16,face = "bold"),
        axis.title.y=element_text(size=16,face = "bold"),
        text = element_text(family = "Times New Roman"), 
        plot.margin = unit(c(0,0,0,1),"cm"))

legend.add <- get_legend(legend) 

#Compile graphs
compile1 <- cowplot::plot_grid(con_grouped, rem_grouped, legend.add,
                               ncol = 3, scale = 1,rel_widths = c(1,1,.5),
                               labels = c("a) CON", "b) REM"), align = "v",
                               label_fontfamily = "Times New Roman", label_size = 14.5,
                               label_y = c(0.72,0.72),label_x = c(0.03,0.03));compile1
compile2 <- cowplot::plot_grid(ran_grouped, bysp_grouped, xsp_grouped,
                               ncol = 3, scale = 1, rel_widths = c(1,1,1),
                               labels = c("c) NAG", "d) SSP", "e) XSP"), align = "v",
                               label_fontfamily = "Times New Roman", label_size = 14,
                               label_y = c(0.735,0.735,0.73),label_x = c(0.04,0.05,0.05));compile2     
compile3 <- cowplot::plot_grid(con_grouped, rem_grouped, transp_grouped,
                               ncol = 1, scale = 1, 
                               labels = c("a) CON", "b) REM", "c) Transplant Treatments"), align = "v",
                               label_fontfamily = "Times New Roman", label_size = 12,
                               label_y = c(0.95,0.95,0.95),label_x = c(0.38,0.38,0.32));compile3




ggsave(plot=compile1,path="C:/Users/Corinne.Amir/Documents/Miscellaneous/Cmorph Removal Manuscript/Figures",
       width=11,height=8,filename=paste0("Grouped Community Barplot_toprow_11-17-2020.jpg"))
ggsave(plot=compile2,path="C:/Users/Corinne.Amir/Documents/Miscellaneous/Cmorph Removal Manuscript/Figures",
       width=11,height=8,filename=paste0("Grouped Community Barplot_bottomrow_11-17-2020.jpg")) 
ggsave(plot=compile3,path="C:/Users/Corinne.Amir/Documents/Miscellaneous/Cmorph Removal Manuscript/Figures",
       width=11,height=8,filename=paste0("Grouped Community Barplot_CON_REM_Transp_11-17-2020.jpg")) 
ggsave(plot=legend,path="C:/Users/Corinne.Amir/Documents/Miscellaneous/Cmorph Removal Manuscript/Figures",
       width=11,height=8,filename=paste0("Grouped Community Barplot_axis title_11-17-2020.jpg")) 

#### end ####
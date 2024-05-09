#This r script is for all multivariate analyses of community composition [of the cmorph removal plots across treatment and time]
#Created on 1/12/2020

#Steps performed for multivariate analysis:
#Step 1: Make b-c dissimilarity matrix 
  #make it with untransformed data and ensure that all values are between 1-10
  #if values are not 1-10 sqrt transform

#Step 2: Perform a permanova (9999 permutations) to determine how community comp varied across treatment and year

#Step 2: Perform pairwise comparisons for all variables that are significant in the permanova

#Step 3: Perform simper on factor levels that were significantly different via pairwise comparisons to determine which functional groups are responsible
#        for variation among factor levels

#Step 4: Run an anova (or kruskal-wallis) on the functional groups contributing to the variation in factor levels as determied via simper

#Step 5: Creat an nmds

#### Install packages and upload dataframe ####
library("extrafont") # fonts have to be loaded before ggplot: https://stackoverflow.com/questions/14733732/cant-change-fonts-in-ggplot-geom-text
library(extrafont)
library(gplots)
font_import(paths = "C:/Windows/Fonts", prompt = F)
#extrafont::fonttable(device = "win")
install.packages(c("vegan","tidyr","dunn.test","devtools","cowplot"))
library(vegan)
library(tidyr) 
library(dplyr) 
library(ggplot2)
library(devtools)
library(dunn.test)
library(cowplot)
library(reshape2)

          
raw <- read.csv("C:/Users/corinne.amir/Documents/GitHub/Cmorph_Removal/CSV files/Final Products/RAW_03-11-2024.csv")


#### end ####


#### Create Matrix ####

#Read data hotdog style
a <- dcast(raw, site*year ~ class, sum, value.var = "percent")

total <- rowSums(a[,3:36]) # Calculate total points per site

a <- filter(a, year != 1408) 

a <- a[order(a$year),] # maintain list order 

# b$site <- a$site #Add site row back in
# b$year <- a$year
trt <- c("REM", "RAN", "BYSP", "XSP", "RAN", "REM", "XSP", "BYSP", "RAN", "CON", "CON", "XSP", "CON", "BYSP", "REM")
# b$trt <- rep(trt, times = 11)


# Create same matrix but with tier 3 data
a3 <- dcast(raw, site*year ~ tier3, sum, value.var = "count")
a3$total <- rowSums(a3[,3:11]) # Calculate total points per site
a3 <- filter(a3, year != 1408) 
a3 <- a3[order(a3$year),] # maintain list order 

b3 <- 100*(a3[,3:11]/a3$total)#calculate percentage

a3 <- a3[order(a3$year),] # add metadata back in 
b3$site <- a3$site 
b3$year <- as.factor(a3$year)
b3$trt <- rep(trt, times = 11)

#### end ####


#### Subset and QC dataframe into proper permanova layout ####

#create two dataframes: 1) a species matrix of benthic cover and 2) our grouping variables (treatment and year)
spp.matrix <- as.matrix(b3[,1:9], row.names =NULL)
taxa.matrix <- b3 %>% select(A.acuminata:P.damicornis, row.names=NULL)
taxa.matrix.a <- lapply(taxa.matrix,as.numeric)
View(spp.matrix)

factors.matrix <- b3[,10:12]
factors.matrix$year <- as.factor(factors.matrix$year)


# Create matrix not including 1409
spp.matrix.abr <- b3 %>% filter(year != 1409) %>% select(A.acuminata:P.damicornis, row.names=NULL)
factors.matrix.abr <- filter(factors.matrix, year != 1409)
factors.matrix$year <- as.factor(factors.matrix$year)


#want 'relative abundance' and a table of numerical values
str(spp.matrix) #check! must be numeric = F
str(taxa.matrix) # numeric
str(factors.matrix) # factor


# Want to compare different treatments within the same time point = subset dataset by years
t0 <- filter(b3, year == 1409)
spp.matrix.t0 <- as.matrix(t0[,1:9], row.names =NULL)

t1 <- filter(b3, year == 1509)
spp.matrix.t1 <- as.matrix(t1[,1:9], row.names =NULL)

t2 <- filter(b3, year == 1510)
spp.matrix.t2 <- as.matrix(t2[,1:9], row.names =NULL)

t3 <- filter(b3, year == 1606)
spp.matrix.t3 <- as.matrix(t3[,1:9], row.names =NULL)

t4 <- filter(b3, year == 1609)
spp.matrix.t4 <- as.matrix(t4[,1:9], row.names =NULL)

t5 <- filter(b3, year == 1709)
spp.matrix.t5 <- as.matrix(t5[,1:9], row.names =NULL)

t6 <- filter(b3, year == 1809)
spp.matrix.t6 <- as.matrix(t6[,1:9], row.names =NULL)

t7 <- filter(b3, year == 2019)
spp.matrix.t7 <- as.matrix(t7[,1:9], row.names =NULL)

t8 <- filter(b3, year == 2020)
spp.matrix.t8 <- as.matrix(t8[,1:9], row.names =NULL)

t9 <- filter(b3, year == 2021)
spp.matrix.t9 <- as.matrix(t9[,1:9], row.names =NULL)

t10 <- filter(b3, year == 2022)
spp.matrix.t10 <- as.matrix(t10[,1:9], row.names =NULL)

factors.matrix.t <- b3 %>% filter(year == 1409) %>% select(c(trt,year,site))
factors.matrix.t$year <- as.factor(factors.matrix.t$year)


# Want to compare different individual treatments across time
con <- filter(b3, trt == "CON")
spp.matrix.con <- as.matrix(con[,1:9], row.names = NULL)

rem <- filter(b3, trt == "REM")
spp.matrix.rem <- as.matrix(rem[,1:9], row.names = NULL)

bysp <- filter(b3, trt == "BYSP")
spp.matrix.bysp <- as.matrix(bysp[,1:9], row.names = NULL)

xsp <- filter(b3, trt == "XSP")
spp.matrix.xsp <- as.matrix(xsp[,1:9], row.names = NULL)

ran <- filter(b3, trt == "RAN")
spp.matrix.ran <- as.matrix(ran[,1:9], row.names = NULL)


factors.matrix.trt <- b3 %>% filter(trt == "CON") %>% select(c(year,site)) 


#### end ####


#### Perform Bray-Curtis and run PERMANOVA ####
#Check range of values. If values range is > 0-10 transform the data...correct?
range(spp.matrix) # 0-99 = also ok?
range(sqrt(spp.matrix)) # 0-9.95 = OK


# create b-c matrix
perm.bcmatrix.notrans <- vegdist(spp.matrix, method="bray")
perm.bcmatrix.sqrt <- vegdist(sqrt(spp.matrix), method="bray")
perm.bcmatrix.notrans <- vegdist(taxa.matrix, method="bray")
perm.bcmatrix.sqrt <- vegdist(sqrt(taxa.matrix), method="bray")
perm.bcmatrix.abr <- vegdist(sqrt(spp.matrix.abr), method="bray")

perm.bcmatrix.sqrt.t0 <- vegdist(sqrt(spp.matrix.t0), method="bray")
perm.bcmatrix.sqrt.t1 <- vegdist(sqrt(spp.matrix.t1), method="bray")
perm.bcmatrix.sqrt.t2 <- vegdist(sqrt(spp.matrix.t2), method="bray")
perm.bcmatrix.sqrt.t3 <- vegdist(sqrt(spp.matrix.t3), method="bray")
perm.bcmatrix.sqrt.t4 <- vegdist(sqrt(spp.matrix.t4), method="bray")
perm.bcmatrix.sqrt.t5 <- vegdist(sqrt(spp.matrix.t5), method="bray")
perm.bcmatrix.sqrt.t6 <- vegdist(sqrt(spp.matrix.t6), method="bray")
perm.bcmatrix.sqrt.t7 <- vegdist(sqrt(spp.matrix.t7), method="bray")
perm.bcmatrix.sqrt.t8 <- vegdist(sqrt(spp.matrix.t8), method="bray")
perm.bcmatrix.sqrt.t9 <- vegdist(sqrt(spp.matrix.t9), method="bray")
perm.bcmatrix.sqrt.t10 <- vegdist(sqrt(spp.matrix.t10), method="bray")


perm.bcmatrix.sqrt.con <- vegdist(sqrt(spp.matrix.con), method="bray")
perm.bcmatrix.sqrt.ran <- vegdist(sqrt(spp.matrix.ran), method="bray")
perm.bcmatrix.sqrt.xsp <- vegdist(sqrt(spp.matrix.xsp), method="bray")
perm.bcmatrix.sqrt.rem <- vegdist(sqrt(spp.matrix.rem), method="bray")
perm.bcmatrix.sqrt.bysp <- vegdist(sqrt(spp.matrix.bysp), method="bray")


#Make a heat map
#red cells = low/no change/difference?
#red cells/closer to 0 = less correlation?
heatmap.2(spp.matrix, key = T) 


# check PERMANOVA assumption of equal dispersion among groups 
# if the dispersion is not significantly different among treatments go straight ahead with the permanova
dispersion <- betadisper(perm.bcmatrix.notrans, factors.matrix$trt)
dispersion
boxplot(dispersion)
anova(dispersion) #No sig differences = continue onto permanova!

dispersion <- betadisper(perm.bcmatrix.sqrt, factors.matrix$trt)
dispersion
boxplot(dispersion)
anova(dispersion) #No sig differences = continue onto permanova!
TukeyHSD(dispersion, ordered = FALSE, conf.level = 0.95) # none significant 

dispersion <- betadisper(perm.bcmatrix.abr, factors.matrix.abr$trt)
TukeyHSD(dispersion, ordered = FALSE, conf.level = 0.95) # REM-CON p = 0.02
anova(dispersion) # p = 0.046
boxplot(dispersion)
dispersion


# By time point
dispersion <- betadisper(perm.bcmatrix.sqrt.t0, factors.matrix.t$trt); anova(dispersion) # all insignificant
dispersion <- betadisper(perm.bcmatrix.sqrt.t1, factors.matrix.t$trt); anova(dispersion) 
dispersion <- betadisper(perm.bcmatrix.sqrt.t2, factors.matrix.t$trt); anova(dispersion)
dispersion <- betadisper(perm.bcmatrix.sqrt.t3, factors.matrix.t$trt); anova(dispersion)
dispersion <- betadisper(perm.bcmatrix.sqrt.t4, factors.matrix.t$trt); anova(dispersion)
dispersion <- betadisper(perm.bcmatrix.sqrt.t5, factors.matrix.t$trt); anova(dispersion) 
dispersion <- betadisper(perm.bcmatrix.sqrt.t6, factors.matrix.t$trt); anova(dispersion)
dispersion <- betadisper(perm.bcmatrix.sqrt.t7, factors.matrix.t$trt); anova(dispersion) 
dispersion <- betadisper(perm.bcmatrix.sqrt.t8, factors.matrix.t$trt); anova(dispersion)
dispersion <- betadisper(perm.bcmatrix.sqrt.t9, factors.matrix.t$trt); anova(dispersion)
dispersion <- betadisper(perm.bcmatrix.sqrt.t10, factors.matrix.t$trt)
boxplot(dispersion)
anova(dispersion) #No sig differences = continue onto permanova!


# By treatment
dispersion <- betadisper(perm.bcmatrix.sqrt.bysp, factors.matrix.trt$year); anova(dispersion) # all insignificant
dispersion <- betadisper(perm.bcmatrix.sqrt.xsp, factors.matrix.trt$year); anova(dispersion) 
dispersion <- betadisper(perm.bcmatrix.sqrt.ran, factors.matrix.trt$year); anova(dispersion)
dispersion <- betadisper(perm.bcmatrix.sqrt.rem, factors.matrix.trt$year); anova(dispersion) 
dispersion <- betadisper(perm.bcmatrix.sqrt.con, factors.matrix.trt$year); anova(dispersion) 
TukeyHSD(dispersion, ordered = FALSE, conf.level = 0.95)


#Do the same for year with two factor permanova
dispersion.yr <- betadisper(perm.bcmatrix.notrans, factors.matrix$year) 
dispersion.yr
boxplot(dispersion.yr)
anova(dispersion.yr) #No sig differences = continue onto permanova

dispersion.yr <- betadisper(perm.bcmatrix.sqrt, factors.matrix$year) # warning: some sq distance negative but changed to zero
dispersion.yr
boxplot(dispersion.yr)
anova(dispersion.yr) # p = 0.03

dispersion.yr <- betadisper(perm.bcmatrix.abr, factors.matrix.abr$year)
TukeyHSD(dispersion.yr, ordered = FALSE, conf.level = 0.95) # none significant




#Run PERMANOVAs
  #Using taxa or spp.matrix yields the same results
  #Variables and interaction are super significant regardless of factor order
  #Add a holms correction to the following p values: 0.05/number of tests - rank number of pair + 1
    # 0.05/1-20+1 = 

set.seed(30) #reproducible results

#With sqrt transform and all years
perm.output.trans <-  adonis2(sqrt(spp.matrix)~trt*year, data=factors.matrix, permutations = 9999, method='bray') 
perm.output.trans # all significant 
                  # F-statistic: year = 77, trt = 49, interaction = 4 (bigger = more variation)
                  # R2: year = 63%, trt = 16%, interaction = 12%, total = 91% (percentage of variance explained)

#Without sqrt transform and with all years
perm.output.notrans <- adonis2(taxa.matrix~trt*year, data=factors.matrix, permutations = 9999, method='bray') 
View(perm.output.notrans) # all significant  
                          # F-statistic: year = 70, trt = 27, interaction = 4
                          # R2: year = 65%, trt = 10%, interaction = 15%, total = 90% 

#With sqrt transform and starting at 1509
perm.output.abr <- adonis2(spp.matrix.abr~trt*year, data=factors.matrix.abr, permutations = 9999, method='bray') 
perm.output.abr # all significant  
# F-statistic: year = 19.1, trt = 20.6, interaction = 1.8
# R2: year = 41%, trt = 20%, interaction = 15%, total = 90% 

# Separated by time point 
perm.output.t0 <- adonis2(spp.matrix.t0~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t0 # p = 0.009, F = 155, R2 = 98%

perm.output.t1 <- adonis2(spp.matrix.t1~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t1 # p = 0.005, F = 3.7, R2 = 59%

perm.output.t2 <- adonis2(spp.matrix.t2~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t2 # p = 0.019, F = 5.1, R2 = 67%

perm.output.t3 <- adonis2(spp.matrix.t3~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t3 # p = 0.051, F = 2.7, R2 = 52%

perm.output.t4 <- adonis2(spp.matrix.t4~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t4 # p = NS, F = 1.1, R2 = 31%

perm.output.t5 <- adonis2(spp.matrix.t5~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t5 # p = 0.003, F = 5.2, R2 = 67%

perm.output.t6 <- adonis2(spp.matrix.t6~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t6 # p = 0.024, F = 2.8, R2 = 53%

perm.output.t7 <- adonis2(spp.matrix.t7~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t7 # p = NS, F = 2.4, R2 = 49%

perm.output.t8 <- adonis2(spp.matrix.t8~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t8 # p = 0.010, F = 3.1, R2 = 56%

perm.output.t9 <- adonis2(spp.matrix.t9~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t9 # p = 0.003, F = 4.7, R2 = 65%

perm.output.t10 <- adonis2(spp.matrix.t10~trt, data=factors.matrix.t, permutations = 9999, method='bray') 
perm.output.t10 # p = 0.001, F = 6.5, R2 = 72%


perm.output.bysp <- adonis2(spp.matrix.bysp~year, data=factors.matrix.trt, permutations = 9999, method='bray') 
perm.output.bysp # p < 0.001, F = 19.0, R2 = 90% 
perm.output.xsp <- adonis2(spp.matrix.xsp~year, data=factors.matrix.trt, permutations = 9999, method='bray') 
perm.output.xsp # p < 0.001, F = 16.6, R2 = 88% 
perm.output.ran <- adonis2(spp.matrix.ran~year, data=factors.matrix.trt, permutations = 9999, method='bray') 
perm.output.ran # p < 0.001, F = 22.7, R2 = 91% 
perm.output.rem <- adonis2(spp.matrix.rem~year, data=factors.matrix.trt, permutations = 9999, method='bray') 
perm.output.rem # p < 0.001, F = 20.6, R2 = 90% 
perm.output.con <- adonis2(spp.matrix.con~year, data=factors.matrix.trt, permutations = 9999, method='bray') 
perm.output.con # p < 0.001, F = 11.8, R2 = 84%


# look at f-values from the 999 permuations
densityplot(permustats(perm.output.t10))


#### end ####


#### Format data for NMDS plots ####

# Calculate mean value for each trt x year
cols <- c("A.acuminata" , "CCA" , "Corallimorph" ,"M.capitata","nonbiological", "other.calcifying.algae" 
          ,"other.coral" ,"other.noncalcifying", "P.damicornis")       
spp.matrix.mean<- b3 %>% group_by(trt, year) %>% summarise(across(all_of(cols),mean))
spp.matrix.mean.abr<- b3 %>% filter(year != 1409) %>% group_by(trt, year) %>% summarise(across(all_of(cols),mean))


#ordination by NMDS, try 2 dimensions (k=2)
nmds <- metaMDS(spp.matrix.mean[,3:11], distance = "bray", k = 2, maxit = 999, trymax = 250, wascores = TRUE)
nmds$stress # 0.13 (0.09 with k = 3)

nmds.sqrt <- metaMDS(sqrt(spp.matrix), distance = "bray", k = 2, maxit = 999, trymax = 250, wascores = TRUE)
nmds$stress # 0.13 

nmds.abr <- metaMDS(spp.matrix.mean.abr[,3:11], distance = "bray", k = 2, maxit = 999, trymax = 250, wascores = TRUE)
nmds.abr$stress # 0.17 (0.10 with k = 3, 0.17 without sqrt)


# Base r plot
plot(nmds.abr$points, col = factors.matrix.abr$year, pch = factors.matrix.abr$trt, key =T)
ordihull(nmds.abr, groups = factors.matrix.abr$trt)



#### Plot NMDS successional trajectory ####
# Plot with ggplot
#add points (note that the factor you want to group by must be a 'factor' in the dataframe)
# Year=bytrtyear$Year
# TRT=bytrtyear$TRT
# Year.abr=bytrtyear.abr$Year
# TRT.abr=bytrtyear.abr$TRT

#Extract xy coordinates for NMDS
MDS1 = nmds$points[,1]  # all data
MDS2 = nmds$points[,1] 
MDS3 = nmds$points[,1] # sqrt transfrom
MDS4 = nmds$points[,2]
MDS1.abr = nmds.abr$points[,1]  # all data
MDS2.abr = nmds.abr$points[,1]
MDS3.abr = nmds.abr$points[,1] # sqrt transfrom
MDS4.abr = nmds.abr$points[,2]


#Create dataframe for ggplot
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3, MDS4 = MDS4, 
                  trt = factors.matrix$trt, year = factors.matrix$year)

NMDS.abr = data.frame(MDS1.abr = MDS1.abr, MDS2.abr = MDS2.abr,# MDS3.abr = MDS3.abr, MDS4.abr = MDS4.abr,
                  trt.abr = spp.matrix.mean.abr$trt, year.abr = spp.matrix.mean.abr$year)


# Order levels for legend......HAVE TO DO THIS WITH EVERY OTHER FIGURE
NMDS$year <- as.factor(NMDS$year)
NMDS$trt <- factor(NMDS$trt, levels = c("XSP", "BYSP", "RAN","REM","CON")) 

NMDS.abr$year.abr <- as.factor(NMDS.abr$year.abr)
NMDS.abr$trt.abr <- factor(NMDS.abr$trt.abr, levels = c("XSP", "BYSP", "RAN","REM","CON")) 


 
#Plot nmds with sqrt transform 1408-1809 (successional trajectory)
nmds.sqrt <- ggplot(NMDS, aes(x=MDS3, y=MDS4, color=trt, shape=year)) +
            geom_point(data = NMDS[!NMDS$year==1910,], size=4) +
            theme_bw() +
            # scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4')) +
            scale_shape_manual(values = c("1408"=9, '1409'=8,"1509"=17,"1510" = 18,'1606'=19,"1609"=17,'1709'=15,'1809'=3,'2019'=1,'2020'=2,'2021'=4,'2022'=5),
                               limits=c("1408", '1409',"1509","1510", '1606',"1609",'1709','1809','2019','2020',"2021",'2022'),
                               labels = c("Aug. 2014", "Sept. 2014","Sept. 2015","Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017","Sept. 2018","Oct. 2019","2020","2021","2022")) +
            # scale_shape_manual(values = c("1408" = 19, '1409'=19,"1509"=19,"1510" = 19, '1606'=19,"1609"=19,'1709'=19,'1809'=19,'1910'=3),
            #           limits=c("1408", '1409',"1509","1510", '1606',"1609",'1709','1809'),#guide=F,
            #           labels = c("Aug. 2014", "Sept. 2014","Sept. 2015","Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017","Sept. 2018","Oct. 2019")) +
            # geom_path(aes(x=MDS1, y=MDS2, group=trt), size=1.2, show.legend = F,
            #           arrow = arrow(type = "closed", angle = 30, length = unit(0.1, "inches"))) +
            theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  text = element_blank(),
                  axis.ticks = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=2)) +
                  geom_text(x=-1.015, y=.78, label="a)", color="black", size=11, family="Times New Roman", fontface="bold")
            #annotate(geom = "rect", xmin=-.1, xmax=.4, ymin=-.31, ymax=.17, color = "red", fill = NA, size = 2) 


#Plot years 1509-1809 with sqrt transformation (successional trajectory)
nmds.sqrt.abr <- ggplot(NMDS.abr, aes(x=MDS3.abr, y=MDS4.abr, colour=trt.abr, shape=year.abr)) +
  geom_point( size=4) + 
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'))    +
  scale_shape_manual(values = c("1509"=17,"1510" = 18, '1606'=19,"1609"=17,'1709'=15,'1809'=3,'1910'=1,'2020'=2,'2021'=4,'2022'=5), 
                     limits=c("1510", '1606',"1609",'1709','1809',"1910",'2020',"2021",'2022'),
                      labels = c("Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017","Sept. 2018","Oct. 2019","2020","2021","2022")) +
  geom_path(aes(x=MDS3.abr, y=MDS4.abr, group=trt.abr), size=1.3, show.legend = F,
            arrow = arrow(type = "closed", angle = 30, length = unit(0.15, "inches"))) +
  ylim(-0.2,0.33) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
        geom_text(x=-0.42, y=0.32, label="b)", color="black", size=11, family="Times New Roman", fontface="bold")
        

# Make the legend      
legend <- ggplot(NMDS, aes(x=MDS1, y=MDS2, colour=TRT, shape=Year)) +
  geom_point(stroke = 1.5,size = 7) + 
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'))    +
  scale_shape_manual(values = c("1510" = 19, '1606'=19,"1609"=19,'1709'=19,'1809'=19,'1910'=3), limits=c("1510", '1606',"1609",'1709','1809','1910'),
                     labels = c("Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017","Sept. 2018", "Oct. 2019"), guide = F) +
  labs(color = "Treatment") +
  theme(legend.title = element_text(face="bold",family="Times New Roman",size=21),
        legend.position = c(0.13, 0.5),
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=21,family="Times New Roman"),
        legend.key.size = unit(.5, "cm"),
        legend.key.width = unit(0.2,"cm"))     
legend.add <- get_legend(legend)         



#Plot figures together
compile <- plot_grid(nmds.sqrt, nmds.sqrt.abr, legend.add, ncol = 3, scale = 1,align = "v" 
                     ,rel_widths = c(1,1,.45)
                     ,labels = c("stress = 0.09","stress = 0.11")
                     ,label_fontfamily = "Times New Roman", label_size = 17
                     ,label_x = c(0.58,0.58),label_y = c(0.98,0.98)); compile
                     

#Export figures
ggsave(plot=compile,path="C:/Users/Corinne.Amir/Documents/Miscellaneous/Cmorph Removal Manuscript/Figures",
       width=12,height=6,filename=paste0("NMDS_compile_11-17-2020.jpg"))




#### Plot NMDS with correlation vectors ####
# Calculate correlation vectors
vec.groups <- envfit(nmds$points, spp.matrix.mean[,3:11], perm = 1000)
vec.groups.abr <- envfit(nmds.abr$points, spp.matrix.mean.abr[,3:11], perm = 1000)

# Make dataframe containing distance values and factors
scrs <- as.data.frame(scores(nmds,display="sites")) 
scrs <- cbind(scrs, Treatment = spp.matrix.mean$trt, Year = spp.matrix.mean$year)

scrs.abr <- as.data.frame(scores(nmds.abr,display="sites")) 
scrs.abr <- cbind(scrs.abr, Treatment = spp.matrix.mean.abr$trt, Year = spp.matrix.mean.abr$year)

#This SHOULD scale the envfit values with the R2 to create the correlation vectors
spp.scrs <- as.data.frame(scores(vec.groups, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))

spp.scrs.abr <- as.data.frame(scores(vec.groups.abr, display = "vectors"))
spp.scrs.abr <- cbind(spp.scrs.abr, Species = rownames(spp.scrs.abr))

#Change name of year and groups
scrs <- scrs %>% #ungroup(Year) %>% # Fix x axis names
  mutate(Year = factor(Year, levels = c(1408,1409,1509,1510,1606,1609,1709,1809,1910,2020,2021,2022),
                       labels = c("Aug. 2014", "Sept. 2014", "Sept. 2015", "Oct. 2015", "Jun. 2016", "Sept. 2016", "Sept. 2017", "Sept. 2018","Oct. 2019","2020","2021","2022")))
  
spp.scrs <- spp.scrs %>%
  mutate(Species = factor(Species, levels = cols,
                       labels = cols))

scrs.abr <- scrs.abr %>% #ungroup(Year) %>% # Fix x axis names
  mutate(Year = factor(Year, levels = c(1510,1606,1609,1709,1809,1910,2020,2021,2022),
                       labels = c("Oct. 2015", "Jun. 2016", "Sept. 2016", "Sept. 2017", "Sept. 2018","Oct. 2019","2020","2021","2022")))

spp.scrs.abr <- spp.scrs
 


# Order levels for legend......HAVE TO DO THIS WITH EVERY OTHER FIGURE
scrs$Treatment <- factor(scrs$Treatment, levels = c("XSP", "BYSP", "RAN","REM","CON")) 
scrs.abr$Treatment <- factor(scrs.abr$Treatment, levels = c("XSP", "BYSP", "RAN","REM","CON")) 

NMDS.abr$Year.abr <- as.factor(NMDS.abr$Year.abr)
NMDS.abr$TRT.abr <- factor(NMDS.abr$TRT.abr, levels = c("XSP", "SSP", "NAG","REM","CON")) 
scale_shape_manual(values = c("1408"=9, '1409'=8,"1509"=17,"1510" = 18,'1606'=19,"1609"=17,'1709'=15,'1809'=3,'2019'=1,'2020'=2,'2021'=4,'2022'=5),
                   limits=c("1408", '1409',"1509","1510", '1606',"1609",'1709','1809','2019','2020',"2021",'2022'),
                   labels = c("Aug. 2014", "Sept. 2014","Sept. 2015","Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017","Sept. 2018","Oct. 2019","2020","2021","2022")) +

#Plot 2014-2018 (sqrt transform)
nmds.vec.sqrt <- ggplot(scrs) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Treatment, shape = Year),stroke =1.5,size = 4.5) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), size = 1,
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  # xlim(-1,1.05) + ylim(-0.5,1) +
  geom_text(data = spp.scrs, aes(x = MDS1, y = MDS2, label = Species),
            size = 4, family="Times New Roman") +
  #labs(shape = "Sampling Period") +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'))    +
  scale_shape_manual(values=c(0,21,22,15,19,17,18,13,3,2,4,5)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(family = "Times New Roman",color="black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) 


#Plot 2015-2018 (sqrt transform)
nmds.vec.sqrt.abr <- ggplot(scrs.abr) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Treatment, shape = Year),stroke = 1.5,size = 4.5) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs.abr,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), size = 1,
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  xlim(-0.9,1) + ylim(-0.5,1) +
  geom_text(data = spp.scrs.abr, aes(x = MDS1, y = MDS2, label = Species),
            size = 4.5, family="Times New Roman") +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'),guide=F)    +
  scale_shape_manual(values=c(15,19,17,18,13,3),guide=F) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman",color="black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2))


#Add legend for plots without successional trajectory
legend.add <- ggplot(scrs) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Treatment, shape = Year),stroke = 1.5,size = 4.5) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), size = 1,
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs, aes(x = MDS1, y = MDS2, label = Species),
            size = 4.5, family="Times New Roman") +
  labs(shape = "Sampling Period") +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'))    +
  scale_shape_manual(values=c(0,21,22,15,19,17,18,13)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(face="bold",family="Times New Roman",size=14),
        legend.position = c(0.2, 0.5),
        legend.text = element_text(size=14,family="Times New Roman"),
        legend.key.size = unit(.5, "cm"),
        legend.key.width = unit(0.5,"cm"))  
legend.add <- get_legend(legend.add) 



#Plot 2014-2018 with successional trajectory and benthic group vectors(sqrt transform)
nmds.combine.sqrt.abr <-  ggplot(scrs.abr) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Treatment, size = Year)) +
  geom_point(data = scrs.abr[!scrs.abr$Year.abr==1910,]) + 
  scale_size_manual(values=c(3,3,3,3,3,0.001), guide = F) +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'),guide=F)  +
  coord_fixed() + ## need aspect ratio of 1!
  xlim(-.95,.8) + ylim(-.74,.733) +
  geom_segment(data = spp.scrs.abr,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), size = 1,
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs.abr, aes(x = MDS1,y = MDS2,label = Species),
            size = 4.5, family="Times New Roman") +
  geom_path(aes(x = NMDS1, y = NMDS2, group=Treatment, color = Treatment), size=1., show.legend = F,
            arrow = arrow(type = "closed", angle = 30, length = unit(0.11, "inches"))) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman",color="black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2))


# Plot 2014-2019 using successional trajectory and vectors
nmds.combine.sqrt <- ggplot(scrs) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Treatment, size = Year)) +
  scale_size_manual(values=c(3,3,3,3,3,3,3,3,0.001), guide = F) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_path(aes(x = NMDS1, y = NMDS2, group=Treatment, colour = Treatment), size=1., show.legend = F,
            arrow = arrow(type = "closed", angle = 30, length = unit(0.11, "inches"))) +
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), size = 1,
               arrow = arrow(length = unit(0.22, "cm")), colour = "grey") +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'),guide=F) +
  xlim(-1.15,1.1) + ylim(-.71,1.06) +
  geom_text(data = spp.scrs, aes(x = MDS1, y = MDS2, label = Species),
            size = 4.5, family="Times New Roman") +
  # scale_shape_manual(values=c(15,19,17,18,13,3),guide=F) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman",color="black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2))



#Add legend for plots with successional trajectory
legend.add <- ggplot(scrs) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Treatment, size = Year)) +
  scale_size_manual(values=c(4,4,4,4,4,4,4,4,0.001), guide = F) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_path(aes(x = NMDS1, y = NMDS2, group=Treatment, colour = Treatment), size=1.3, show.legend = T,
            arrow = arrow(type = "closed", angle = 30, length = unit(0.1, "inches"))) +
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), size = 1,
               arrow = arrow(length = unit(0.22, "cm")), colour = "grey") +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4')) +
  geom_text(data = spp.scrs, aes(x = MDS1, y = MDS2, label = Species),
            size = 4.5, family="Times New Roman") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(face="bold",family="Times New Roman",size=18),
        legend.position = c(0.1, 0.5),
        legend.background = element_blank(),
        legend.text = element_text(size=18,family="Times New Roman"),
        legend.key.size = unit(.6, "cm"),
        legend.key.width = unit(0.5,"cm"))  
legend.add <- get_legend(legend.add) 


#Plot figures together
compile <- plot_grid(nmds.combine.sqrt, nmds.combine.sqrt.abr, legend.add, ncol = 3, scale = 1
                     ,rel_widths = c(1.065,1,.45)
                     ,labels = c("a)", "b)"),label_x = c(0.03,0.03),label_y = c(0.85,0.85)
                     ,label_fontfamily = "Times New Roman", label_size = 18);compile
                     #,labels = c("stress = 0.08","stress = 0.13")
                     #,label_x = c(0.61,0.59),label_y = c(0.91,0.91)); compile


#Export figures
ggsave(plot=compile,path="C:/Users/Corinne.Amir/Documents/Miscellaneous/Cmorph Removal Manuscript/Figures",
       width=9,height=4,filename=paste0("NMDS_wVectors_compile_11-18-2020.jpg"))



#### end ####


#### Perform Pairwise Comparisons ####

# Pairwise function  
pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


# pairwise adonis utilizing non-subsetted matrix
bencov<-pairwise.adonis(sqrt(spp.matrix), factors.matrix$year,'bray','bonferroni')
bencov # none statistically different

bencov<-pairwise.adonis(sqrt(spp.matrix), factors.matrix$trt,'bray','bonferroni')
bencov # CON and REM statistically different from transplant treatments (F = 8.5-16.2, R2 = 11.7-20.2, p = 0.01)
       # CON and REM statistically different from each other (F = 5.0, R2 = 7.2, p = 0.03)


# Subset by year
bencov<-pairwise.adonis(spp.matrix.t0, factors.matrix.t$trt,'bray','bonferroni')
bencov # none significant 

bencov<-pairwise.adonis(spp.matrix.t10, factors.matrix.t$trt,'bray','bonferroni')
bencov # none significant


# Subset by treatment
bencov<-pairwise.adonis(spp.matrix.bysp, factors.matrix.trt$year,'bray','bonferroni')
bencov # none significant


#### end pairwise tests ####


#### SIMPER ####
#Can run simper for my personal benefit and then say in the apper that I ran two-way anovas on [all] benthic functional groups and 
#then just report the ones that were significant/showed up in the simper

#Using sqrt abridged data
sim.trt.trans <- simper(spp.matrix.abr, group = factors.matrix.abr$TRT)
sim.trt.trans 
#BYSP-CON: turf (30), CCA (26), Corallimorph (22)
#RAN-CON: turf (32), cca, corallimorph
#XSP-CON: turf (32), CCA, corallimorph

#CON-REM: cca (32), turf, corallimorph

#BYSP-REM: turf (36), cca (34)
#RAN-REM: turf (38), CCA
#XSP-REM: cca (36), turf
summary(sim.trt.trans)
#BYSP-CON: BYSP has more turf, less cca, and less cmorph than CON
#RAN-CON: RAN " " than CON
#XSP-CON: CON " " than CON

#CON-REM: REM has more CCA, same turf, and less cmorph than CON

#BYSP-REM: BYSP has more turf, less cca and less calc.algae than REM
#RAN-REM: RAN " " 
#XSP-REM: XSP has less CCA and calc.algae, more turf and coral

sim.yr.trans <- simper(spp.matrix.abr, group = factors.matrix.abr$Year)
sim.yr.trans 
#1709-1510: turf (34), cca, cmorph
#1709-1606: turf (39), cca
#1709-1609: turf (41), cca

#1709:1809: turf (34), cca, coral

#1809-1510: turf (29), cca, cmorph
#1809-1606: turf (32), cca, coral
#1809-1609: turf (34), cca, coral

summary(sim.yr.trans)
#1709-1510: 1510 had more turf and cmorph, less CCA
#1709-1606: 1606 had more turf, less cca and calc.algae
#1709-1609: 1609 had more turf, less CCA and calc.algae

#1709:1809: 1709 had less turf, more cca, less coral 

#1809-1510: 1510 had more turf and cmorph, less CCA and coral
#1809-1606: 1606 had more turf, less cca and coral
#1809-1609: 1609 had more turf, less cca and coral

## Conclusion:
# 1709 was different from other years due to greater amounts of cca (due to decreases in turf caused by decreases in cmorph)
# In 1809, cca decreased from 1709 (replaced by turf), but was still greater than time points prior to 1709. 
# 1809 has more coral than all years too
# CON and REM had more CCA than all other treatments (turf greater in transplant plots instead)
# RUN ANOVAS ON: cca, turf, coral, cmorph


#unabridged and untransfromed
sim.yr <- simper(spp.matrix, group = factors.matrix$Year, ordered=T)
sim.yr # these are the functional groups to use for testing for univariate differences
#isted in order by most influential to the differnces between year.  
#Note that "coral" does not contribute to differences between 2016-2017. 
summary(sim) #ava, avb are the average abundances per group.
#

sim.trt <- simper(taxa.matrix, group = factors.matrix$TRT)
sim.trt
summary(sim.trt) #ava, avb are the average abundances per group. average = average contribution to overall dissimilarity. 
#CON different than trts because of Turf (30) and cmorph (29) and cca (16)
#REM different than CON because of corallimorph (30), turf (18), and CCA (21)
#REM different than treatments becuase of turf (31), nonbio (21), and cca (20)



#SIMPER using sqrt transformed and unabridged data
sim.trt.trans <- simper(comm = sqrt(taxa.matrix), group = factors.matrix$TRT)
sim.trt.trans 
summary(sim.trt.trans)
#answers are cumulative!!! adds functional groups until you reach about 70% contribution
#CON sig different than other plots because of cmorph
#RAN/XSP/BYSP different than CON because of corallimorph, turf, then cca, then nonbio
#REM different than CON because of corallimorph (27%), cca (20%), turf (19%), nonbio (15%)

sim.yr.trans <- simper(comm = sqrt(taxa.matrix), group = factors.matrix$Year)
sim.yr.trans 
summary(sim.yr.trans)
#overall contributions are very low b/c of transform
#seems to give same conclusions and non-transformed SIMPER
#1408 different than all years because of cmorph (~37)
#1408 different than all years because of nonbio (~37)
#1509 dfferent than all years because of CCA (~25). can see corallimorph still high as well. Except for 1809 comparison fleshy algae (24%)
#1709 different than all years (except 1809) because of CCA. Different than 1809 because of coral
#1809 different than all years because of coral (~23)

#### end SIMPER ####


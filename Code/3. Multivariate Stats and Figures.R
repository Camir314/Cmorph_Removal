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
install.packages("extrafont") # fonts have to be loaded before ggplot: https://stackoverflow.com/questions/14733732/cant-change-fonts-in-ggplot-geom-text
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

          
raw <- read.csv("C:/Users/corinne.amir/Documents/GitHub/Cmorph_Removal/CSV files/Final Products/RAW_03-06-2024.csv")


#### end ####


#### Create Matrix ####

#Read data hotdog style
a <- dcast(raw, site*year ~ class, sum, value.var = "count")

a$total <- rowSums(a[,3:36]) # Calculate total points per site

a <- filter(a, year != 1408) 

a <- a[order(a$year),] # maintain list order 

#calculate percentage
b <- a[,3:36] / a[,37]
b <- b[1:34]*100

total <- rowSums(b[,1:34]) # QC check


a <- a[order(a$year),]
b$site <- a$site #Add site row back in
b$year <- a$year
trt <- c("REM", "RAN", "BYSP", "XSP", "RAN", "REM", "XSP", "BYSP", "RAN", "CON", "CON", "XSP", "CON", "BYSP", "REM")
b$trt <- rep(trt, times = 11)


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


# By time point
dispersion <- betadisper(perm.bcmatrix.sqrt.t0, factors.matrix.t$trt)
dispersion <- betadisper(perm.bcmatrix.sqrt.t1, factors.matrix.t$trt); anova(dispersion) # all insignificant
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


#https://www.researchgate.net/publication/324004081_Drought_and_increased_CO_2_alter_floral_visual_and_olfactory_traits_with_context-dependent_effects_on_pollinator_visitation
# 1. run dispersion test and anovas on each variable tested (choose between distance-to-centroid or other option)
# 2. for variables with significant dispersion, run Tukey to parse out where the problem is
# 3. can talk about how two [treatments] had different magnitudes of variability and why that is ecologically interesting
  #low dispersion could indicate lower turnover and/or more robustness of the community

# file:///C:/Users/corin/Downloads/peerj-2430.pdf
# However, the PERMANOVA results could be affected by non-homogeneous dispersion of the data... Nevertheless, the ordination
# clearly supports a clustering pattern that confirms the importance of habitat and host species factors combined.

#Do the same for year with two factor permanova
dispersion.yr <- betadisper(perm.bcmatrix.notrans, factors.matrix$year) 
dispersion.yr
boxplot(dispersion.yr)
anova(dispersion.yr) #No sig differences = continue onto permanova

dispersion.yr <- betadisper(perm.bcmatrix.sqrt, factors.matrix$year) # warning: some sq distance negative but changed to zero
dispersion.yr
boxplot(dispersion.yr)
anova(dispersion.yr) # p = 0.03
TukeyHSD(dispersion, ordered = FALSE, conf.level = 0.95) # none are significant




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



# look at f-values from the 999 permuations
densityplot(permustats(perm.output.trans))


#### end ####


#### Make nMDS Plots ####
#Make dataframe with mean values for each treatment within each year
CCA_mean <- perm.prop.comb %>% group_by(TRT,Year) %>% dplyr::summarize(CCA = mean(CCA)) 
Turf_mean <- perm.prop.comb %>% group_by(TRT,Year) %>% dplyr::summarize(Turf = mean(Turf))
Corallimorph_mean <- perm.prop.comb %>% group_by(TRT,Year) %>% dplyr::summarize(Corallimorph = mean(Corallimorph)) 
Nonbiological_mean <- perm.prop.comb %>% group_by(TRT,Year) %>% dplyr::summarize(Nonbiological = mean(Nonbiological))
Coral_mean <- perm.prop.comb %>% group_by(TRT,Year) %>% dplyr::summarize(Coral = mean(Coral)) 
Fleshy.Algae_mean <- perm.prop.comb %>% group_by(TRT,Year) %>% dplyr::summarize(Fleshy.Algae = mean(Fleshy.Algae)) 
Calc.Algae_mean <- perm.prop.comb %>% group_by(TRT,Year) %>% dplyr::summarize(Calc.Algae = mean(Calc.Algae))

bytrtyear <- left_join(CCA_mean, Turf_mean) %>%
             left_join(., Corallimorph_mean) %>%
             left_join(., Nonbiological_mean) %>%
             left_join(., Coral_mean) %>%
             left_join(.,Fleshy.Algae_mean) %>%
             left_join(., Calc.Algae_mean)

levels(bytrtyear$TRT)[levels(bytrtyear$TRT)=="BYSP"] <- "SSP" # Fix TRT names
levels(bytrtyear$TRT)[levels(bytrtyear$TRT)=="RAN"] <- "NAG"

#Turn into a matrix
spp.matrix.mean <- bytrtyear[,3:9]
spp.matrix.mean.sqrt <- sqrt(bytrtyear[,3:9])

#Creat abridged version of mean value matrix
bytrtyear.abr <- subset(bytrtyear,Year !="1408" & Year != "1409" & Year != "1509") 
spp.matrix.mean.abr <- bytrtyear.abr[,3:9]
spp.matrix.mean.sqrt.abr <- sqrt(bytrtyear.abr[,3:9])

set.seed(70)
#medaMDS automatically sqrt transforms data and applies bray dissilimarity matrix, which means that we can use the original matrix
NMDS.matrix.mean=metaMDS(spp.matrix.mean, distance = "bray", k=2,trymax = 1000) 
NMDS.matrix.mean$stress #0.09 = good
NMDS.matrix.mean.sqrt=metaMDS(spp.matrix.mean.sqrt, distance = "bray", k=2,trymax = 1000) 
NMDS.matrix.mean.sqrt$stress # 0.09 = good
# NMDS.matrix.mean.asin=metaMDS(asin(0.01*spp.matrix.mean.sqrt), distance = "bray", k=2,trymax = 100) 
# NMDS.matrix.mean.asin$stress # 0.05 = great
# NMDS.matrix.mean.fourth=metaMDS(spp.matrix.mean.sqrt^0.25, distance = "bray", k=2,trymax = 100) 
# NMDS.matrix.mean.fourth$stress # 0.07 = good

#Calculate metaMDS for abridged data:
NMDS.matrix.mean.abr=metaMDS(spp.matrix.mean.abr, distance = "bray", k=2,trymax = 1000) 
NMDS.matrix.mean.abr$stress #0.16 = OK
NMDS.matrix.mean.sqrt.abr=metaMDS(spp.matrix.mean.sqrt.abr, distance = "bray", k=2,trymax = 1000) 
NMDS.matrix.mean.sqrt.abr$stress #0.11 = OK


# Creat NMDS dataframe for plotting
#add points (note that the factor you want to group by must be a 'factor' in the dataframe)
Year=bytrtyear$Year
TRT=bytrtyear$TRT
Year.abr=bytrtyear.abr$Year
TRT.abr=bytrtyear.abr$TRT

#Extract xy coordinates for NMDS
MDS1 = NMDS.matrix.mean$points[,1] # no transform
MDS2 = NMDS.matrix.mean$points[,2]
MDS3 = NMDS.matrix.mean.sqrt$points[,1] # sqrt transfrom
MDS4 = NMDS.matrix.mean.sqrt$points[,2]
# MDS5 = NMDS.matrix.mean.asin$points[,1] # asin sqrt transform 
# MDS6 = NMDS.matrix.mean.asin$points[,2]
# MDS7 = NMDS.matrix.mean.fourth$points[,1] # fourth root transform 
# MDS8 = NMDS.matrix.mean.fourth$points[,2]

MDS1.abr = NMDS.matrix.mean.abr$points[,1] # no transform
MDS2.abr = NMDS.matrix.mean.abr$points[,2]
MDS3.abr = NMDS.matrix.mean.sqrt.abr$points[,1] # sqrt transfrom
MDS4.abr = NMDS.matrix.mean.sqrt.abr$points[,2]

#Create dataframe for ggplot
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3, MDS4 = MDS4, 
                  TRT = bytrtyear$TRT, Year = bytrtyear$Year)
                #  MDS5 = MDS5, MDS6 = MDS6, MDS7 = MDS7, MDS8 = MDS8,) 

NMDS.abr = data.frame(MDS1.abr = MDS1.abr, MDS2.abr = MDS2.abr, MDS3.abr = MDS3.abr, MDS4.abr = MDS4.abr, 
                  TRT.abr = bytrtyear.abr$TRT, Year.abr = bytrtyear.abr$Year)

# Order levels for legend......HAVE TO DO THIS WITH EVERY OTHER FIGURE
NMDS$Year <- as.factor(NMDS$Year)
NMDS$TRT <- factor(NMDS$TRT, levels = c("XSP", "SSP", "NAG","REM","CON")) 

NMDS.abr$Year.abr <- as.factor(NMDS.abr$Year.abr)
NMDS.abr$TRT.abr <- factor(NMDS.abr$TRT.abr, levels = c("XSP", "SSP", "NAG","REM","CON")) 


 
#Plot nmds with sqrt transform 1408-1809 (successional trajectory)
nmds.sqrt <- ggplot(NMDS, aes(x=MDS3, y=MDS4, color=TRT, shape=Year)) +
  geom_point(data = NMDS[!NMDS$Year==1910,], size=4) +
  theme_bw() +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'),guide=F)+
  # scale_shape_manual(values = c("1408" = 9, '1409'=8,"1509"=17,"1510" = 18, '1606'=19,"1609"=17,'1709'=15,'1809'=3),
  #                    limits=c("1408", '1409',"1509","1510", '1606',"1609",'1709','1809'),
  #                    labels = c("Aug. 2014", "Sept. 2014","Sept. 2015","Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017","Sept. 2018")) +
  scale_shape_manual(values = c("1408" = 19, '1409'=19,"1509"=19,"1510" = 19, '1606'=19,"1609"=19,'1709'=19,'1809'=19,'1910'=3),
            limits=c("1408", '1409',"1509","1510", '1606',"1609",'1709','1809'),guide=F,
            labels = c("Aug. 2014", "Sept. 2014","Sept. 2015","Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017","Sept. 2018","Oct. 2019")) +
  geom_path(aes(x=MDS1, y=MDS2, group=TRT), size=1.2, show.legend = F,
            arrow = arrow(type = "closed", angle = 30, length = unit(0.1, "inches"))) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
        geom_text(x=-1.015, y=.78, label="a)", color="black", size=11, family="Times New Roman", fontface="bold")
  #annotate(geom = "rect", xmin=-.1, xmax=.4, ymin=-.31, ymax=.17, color = "red", fill = NA, size = 2) 


#Plot years 1510-1809 with sqrt transformation (successional trajectory)
nmds.sqrt.abr <- ggplot(NMDS.abr, aes(x=MDS3.abr, y=MDS4.abr, colour=TRT.abr, shape=Year.abr)) +
  geom_point(data = NMDS.abr[!NMDS.abr$Year.abr==1910,], size=4) + 
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'),guide=F)    +
  scale_shape_manual(values = c("1510" = 19, '1606'=19,"1609"=19,'1709'=19,'1809'=19,'1910'=3), limits=c("1510", '1606',"1609",'1709','1809',"1910"),
                      labels = c("Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017","Sept. 2018","Oct. 2019"), guide = F) +
  geom_path(aes(x=MDS3.abr, y=MDS4.abr, group=TRT.abr), size=1.3, show.legend = F,
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




# Plot NMDS with correlation vectors
# Calculate correlation vectors
vec.groups <- envfit(NMDS.matrix.mean.sqrt$points, spp.matrix.mean.sqrt, perm = 1000)
vec.groups.abr <- envfit(NMDS.matrix.mean.sqrt.abr$points, spp.matrix.mean.sqrt.abr, perm = 1000)

# Make dataframe containing distance values and factors
scrs <- as.data.frame(scores(NMDS.matrix.mean.sqrt,display="sites")) 
scrs <- cbind(scrs, Treatment = bytrtyear$TRT, Year = bytrtyear$Year)

scrs.abr <- as.data.frame(scores(NMDS.matrix.mean.sqrt.abr,display="sites")) 
scrs.abr <- cbind(scrs.abr, Treatment = bytrtyear.abr$TRT, Year = bytrtyear.abr$Year)

#This SHOULD scale the envfit values with the R2 to create the correlation vectors
spp.scrs <- as.data.frame(scores(vec.groups, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))

spp.scrs.abr <- as.data.frame(scores(vec.groups.abr, display = "vectors"))
spp.scrs.abr <- cbind(spp.scrs.abr, Species = rownames(spp.scrs.abr))

#Change name of year and groups
scrs <- scrs %>% #ungroup(Year) %>% # Fix x axis names
  mutate(Year = factor(Year, levels = c(1408,1409,1509,1510,1606,1609,1709,1809,1910),
                       labels = c("Aug. 2014", "Sept. 2014", "Sept. 2015", "Oct. 2015", "Jun. 2016", "Sept. 2016", "Sept. 2017", "Sept. 2018","Oct. 2019")))
  
spp.scrs <- spp.scrs %>%
  mutate(Species = factor(Species, levels = c("CCA" ,"Turf","Corallimorph","Nonbiological", "Coral", "Fleshy.Algae" ,"Calc.Algae"),
                       labels = c("CCA" ,"Turf Algae","Corallimorph","Nonbiological", "Coral", "Non-Calcified Algae" ,"Calcified Algae")))

scrs.abr <- scrs.abr %>% #ungroup(Year) %>% # Fix x axis names
  mutate(Year = factor(Year, levels = c(1510,1606,1609,1709,1809,1910),
                       labels = c("Oct. 2015", "Jun. 2016", "Sept. 2016", "Sept. 2017", "Sept. 2018","Oct. 2019"))) 

spp.scrs.abr <- spp.scrs.abr %>%
  mutate(Species = factor(Species, levels = c("CCA" ,"Turf","Corallimorph","Nonbiological", "Coral", "Fleshy.Algae" ,"Calc.Algae"),
                          labels = c("CCA" ,"Turf Algae","Corallimorph","Nonbiological", "Coral", "Non-Calcified Algae" ,"Calcified Algae")))


# Order levels for legend......HAVE TO DO THIS WITH EVERY OTHER FIGURE
scrs$Treatment <- factor(scrs$Treatment, levels = c("XSP", "SSP", "NAG","REM","CON")) 
scrs.abr$Treatment <- factor(scrs.abr$Treatment, levels = c("XSP", "SSP", "NAG","REM","CON")) 

NMDS.abr$Year.abr <- as.factor(NMDS.abr$Year.abr)
NMDS.abr$TRT.abr <- factor(NMDS.abr$TRT.abr, levels = c("XSP", "SSP", "NAG","REM","CON")) 


#Plot 2014-2018 (sqrt transform)
nmds.vec.sqrt <- ggplot(scrs) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Treatment, shape = Year),stroke =1.5,size = 4.5) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), size = 1,
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  xlim(-1,1.05) + ylim(-0.5,1) +
  geom_text(data = spp.scrs, aes(x = MDS1, y = MDS2, label = Species),
            size = 4, family="Times New Roman") +
  #labs(shape = "Sampling Period") +
  scale_color_manual(values=c("gold2", "#56B4E9",'purple','hotpink2', 'green4'),guide=F)    +
  scale_shape_manual(values=c(0,21,22,15,19,17,18,13,3),guide=F) +
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

# Using the 1510-1809 data:
# The next step after obtaining a significant interaction is to do pair-wise comparisons for the factor of interest 
# separately within each level of the other factor (i.e., to do pair-wise comparisons among levels of factor A within 
# each level of factor B) and vice versa (if necessary).

# pairwise comparisons: con x 1510,1606,1609,1709,1809  rem x 1510,1606,1609,1709,1809
#                       bysp x 1510,1606,1609,1709,1809  ran x 1510,1606,1609,1709,1809
#                       xsp x 1510,1606,1609,1709,1809  
# and vice versa?
# separate the spp matrix by treatment and then run pairwise.adonis on each dataframe to examine the interactions


#separate dataframes and turn into [sqrt] matrices:
perm.prop.abr <- arrange(perm.prop.abr,perm.prop.abr$TRT)
#still insignificant if you use sqrt transform()
bysp <- droplevels(perm.prop.abr[1:18,])
  bysp.spp <- as.matrix(sqrt(bysp[,4:10]))
  factors.yr <- bysp[,2:3]
con <- perm.prop.abr[19:36,]
  con.spp <- as.matrix(sqrt(con[,4:10]))
ran <- perm.prop.abr[37:54,]
  ran.spp <- as.matrix(sqrt(ran[,4:10]))
rem <- perm.prop.abr[55:72,]
  rem.spp <- as.matrix(sqrt(rem[,4:10]))
xsp <- perm.prop.abr[73:90,]
  xsp.spp <- as.matrix(sqrt(xsp[,4:10]))
  

perm.prop.abr <- arrange(perm.prop.abr,perm.prop.abr$Year)
one <- perm.prop.abr[1:15,]
factors.trt <- one[,2:3]
  one.spp <- as.matrix(sqrt(one[,4:10]))
two <- perm.prop.abr[16:30,]
  two.spp <- as.matrix(sqrt(two[,4:10]))
three <- perm.prop.abr[31:45,]
  three.spp <- as.matrix(sqrt(three[,4:10]))
four <- perm.prop.abr[46:60,]
  four.spp <- as.matrix(sqrt(four[,4:10]))
five <- perm.prop.abr[61:75,]
  five.spp <- as.matrix(sqrt(five[,4:10]))
six <- perm.prop.abr[76:90,]
  sox.spp <- as.matrix(sqrt(five[,4:10]))
  

  #Attempt 1
  #pairwise.adonis() function from 
  pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray',p.adjust.m='holm') 
  {
    library(vegan)
    
    co = combn(unique(as.character(factors)),2)
    pairs = c()
    F.Model =c()
    R2 = c()
    p.value = c()
    
    
    for(elem in 1:ncol(co)){
      if(sim.function == 'daisy'){
        library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
      
      ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
      pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
      F.Model =c(F.Model,ad$aov.tab[1,4]);
      R2 = c(R2,ad$aov.tab[1,5]);
      p.value = c(p.value,ad$aov.tab[1,6])
    }
    p.adjusted = p.adjust(p.value,method=p.adjust.m)
    sig = c(rep('',length(p.adjusted)))
    sig[p.adjusted <= 0.05] <-'.'
    sig[p.adjusted <= 0.01] <-'*'
    sig[p.adjusted <= 0.001] <-'**'
    sig[p.adjusted <= 0.0001] <-'***'
    
    pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
    print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
    return(pairw.res)
    
  } 
  

#run pairwise.adonis...all p = 1
pairwise.adonis(x=bysp.spp,factors=factors.yr$Year)
pairwise.adonis(x=con.spp,factors=factors.yr$Year)
pairwise.adonis(x=ran.spp,factors=factors.yr$Year,sim.function='vegdist',sim.method='bray')
pairwise.adonis(x=rem.spp,factors=factors.yr$Year,sim.function='vegdist',sim.method='bray')
pairwise.adonis(x=xsp.spp,factors=factors.yr$Year,sim.function='vegdist',sim.method='bray')

pairwise.adonis(x=one.spp,factors=factors.trt$TRT,sim.function='vegdist',sim.method='bray') 
pairwise.adonis(x=two.spp,factors=factors.trt$TRT,sim.function='vegdist',sim.method='bray')
pairwise.adonis(x=three.spp,factors=factors.trt$TRT,sim.function='vegdist',sim.method='bray')
pairwise.adonis(x=four.spp,factors=factors.trt$TRT,sim.function='vegdist',sim.method='bray')
pairwise.adonis(x=five.spp,factors=factors.trt$TRT,sim.function='vegdist',sim.method='bray')


#Bonferroni correction would require 0.05/5 = 0.01 significance for trt variable

# #Pair-wise comparisons on year and TRT (although I should be looking at interaction, but this function can't do that)
# pairwise.adonis(x=spp.matrix,factors=factors.matrix$Year,sim.function='vegdist',sim.method='bray',p.adjust.m='bonferroni')
# #1408 and 1409 sig different than all other years, 1709 sig different than 1510 and 1809
# pairwise.adonis(x=spp.matrix),factors=factors.matrix$TRT,sim.function='vegdist',sim.method='bray',p.adjust.m='bonferroni')
# #no TRTs sig different :(
# pairwise.adonis(x=spp.matrix,factors=factors.matrix$TRT,sim.function='vegdist',sim.method='bray',p.adjust.m='fdr')
# #still no TRTs sig different :(

#Running pairwise.adonis, but with sqrt transformed data and abridged time series:
pairwise.adonis(x=sqrt(spp.matrix.abr),factors=factors.matrix.abr$TRT,sim.function='vegdist',sim.method='bray',p.adjust.m='bonferroni') 
#CON different than BYSP, RAN, XSP, REM
#REM different than BYSP, RAN, XSP, CON
pairwise.adonis(x=sqrt(spp.matrix.abr),factors=factors.matrix.abr$Year,sim.function='vegdist',sim.method='bray',p.adjust.m='bonferroni')
#1408 and 1409 sig different than all years (from non abrdged run)
#1709 sig different than 1510, 1606, 1609,and 1809
#1809 also sig different than 1510, 1606, 1609, and 1709

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


#### ANOVAs ####
#Q: Since permanova and pairwise were done with sqrt transform, do all anova have to be done with sqrt?
#A: No....unless coral only shows sig difference if sqrt transformed

#Test interactions 
#cca 1510-1809
cca <- droplevels(select(perm.prop.abr, c(Site, TRT, Year, CCA)))

#Assumptions
bartlett.test(CCA ~ interaction(TRT, as.factor(Year)), data=cca) #Homogenous without transform
tapply(cca$CCA, list(cca$TRT), shapiro.test) # REM = 0.04
tapply(cca$CCA, list(cca$Year), shapiro.test) # 1910 = 0.04

#two-way anova
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


#### Coral Percent Cover Line Graphs ####

ggplot(raw, aes(x = year, y = count, group = TRT, color = TRT)) +
  geom_line(stat = "identity") +
  geom_point(size = 3) +
  facet_wrap(vars(tier3))


# turn data into correct format
coralsp <- coralsp.comb %>% gather(Group, Percent, acro:total) # switch from wide to long format

coralsp  <- coralsp  %>% # ungroup(Year) %>% # this step is sometimes necessary to change the year names
                        mutate(Year = factor(Year, levels = c(1408,1409,1509,1510,1606,1609,1709,1809,2019),labels = c("Aug. 2014","Sept. 2014",
                       "Sept. 2015","Oct. 2015","Jun. 2016","Sept. 2016","Sept. 2017", "Sept. 2018","Oct. 2019"))) 

coralsp$TRT <- as.factor(coralsp$TRT)
levels(coralsp$TRT)[levels(coralsp$TRT)=="BYSP"] <- "SSP" # Fix TRT names
levels(coralsp$TRT)[levels(coralsp$TRT)=="RAN"] <- "NAG"

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


#Calculate SE
coral.st.err <- function(x) {
  sd(x)/sqrt(length(x))
}
SE_coral <- aggregate(Percent~TRT+Month+Group,coralsp,coral.st.err)
#SE_coral <- aggregate(Percent~TRT+Year+Group,coralsp,coral.st.err) # month wasn't working

#Add SE values into data frame
colnames(SE_coral) <- c("TRT"   ,  "Month"  ,  "Group" , "SE")
coral_timeseries_mean <- left_join(coral_timeseries_mean,SE_coral)


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
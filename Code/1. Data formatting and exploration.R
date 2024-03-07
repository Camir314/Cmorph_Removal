# Step 2:
# Format raw data into functional groups
# Identify differences among clickers (Sam, Corinne and Jonny)
# est. March 2024

# Group classes into functional groups

# Compare differences between Corinne, Sam and Jonny
d <- merge %>% group_by(class, clickr) %>% 
  summarise(totcount = count(count)) 


summarise(totsite = count(unique(site_year)))


summarise(meancount = totcount/count(site_year))

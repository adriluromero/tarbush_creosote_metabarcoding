#Jornada_Creosote_Tarbush_Fungal_Metabarcoding

#Protocol for processing Metabarcoding of soil
#Section 7
#File name: 7_Format_table_for_stat_analyses
#Step: Use RStudio to format the output table from FUNGuild for statistical analyses.

#In this script, you will: 
# 1. Split the 'taxonomy' column into 7 columns (Kingdom, Phylum, Class, Order, Family, Genus & Species).
# 2. Move the 7 columns before the taxonomy column.
# 3. Remove the taxonomy column.
# 4. Write the file as csv. 


#Clear environment
rm(list = ls(all.names = TRUE))

#Load the packages
library(stringr)
library(dplyr)

#Import the FUNGuild output file 'FUNGuildCreosoteTarbush.2023.guilds.csv'.
Funguild <- read.csv("/project/egcc/metabarcoding_tarbush_creosote/OTU-tables/2023-OTUs/funguild/FUNGuildCreosoteTarbush.2023.guilds.csv")

#1. Split the 'taxonomy' column into 7 columns (Kingdom, Phylum, Class, Order, Family, Genus & Species).
Funguild[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- str_split_fixed(Funguild$taxonomy, ";", 7)

#2. Move the 7 columns before the taxonomy column.
Funguild <- Funguild %>% relocate(Kingdom:Species, .before = taxonomy)

#3. Remove the taxonomy column. 
Funguild <- Funguild %>% select(-c(taxonomy))

#4. Write the file as csv.
write.csv(Funguild, "/project/egcc/metabarcoding_tarbush_creosote/OTU-tables/2023-OTUs/ForStats/creosote.tarbush.otu.table.csv", row.names = FALSE)

#The next step is to manually curate this dataset. 
#You do so by: 
# A) deleting the samples that have no sequences. For this dataset, refer to the 'Excluded_Samples.txt' file.
# B) removing repeated samples from same plot to avoid pseudoreplication issues. 
# C) removing ASVs with zero counts.  

#Use this for statistical analyses and data visualization. 
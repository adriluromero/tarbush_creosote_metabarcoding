#Jornada_Creosote_Tarbush_Fungal_Metabarcoding

#Protocol for processing Metabarcoding of soil
#Section 5
#File name: 5_Format_ASV_table_for_FUNGuild
#Step: Use RStudio to format the 'otu.tax.table.csv' for input in FUNGuild.

#In this script, you will: 
# 1. Combine the taxonomic information into one column called 'Taxonomy'.
# 2. Finish formatting the table for input in FUNGuild.
# 3. Write a csv 


#Clear environment
rm(list = ls(all.names = TRUE))

#1. Combine the taxonomic information into one column called 'Taxonomy'.

#Import the file. 
ASV <- read.csv("/project/egcc/metabarcoding_tarbush_creosote/OTU-tables/2023-OTUs/otu.tax.table.csv")

#FUNGuild requires one column named "Taxonomy" that follows this format. 
#k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Corticiales;f__Corticiaceae;g__Laetisaria;s__Laetisaria_fuciformis

#Do this by pasting the information from the taxonomic ranks into a column separated by ';'
ASV$taxonomy <- paste(ASV$Kingdom,ASV$Phylum, ASV$Class, ASV$Order, ASV$Family, ASV$Genus, ASV$Species, sep=";")

#2. Finish formatting the table for input in FUNGuild. 

#Rename the first column to "OTU_ID"
colnames(ASV)[1]<-"OTU_ID"

#Since we don't need this for FUNGuild, remove the Kingdom, Phylum, Class, Order, Family, Genus & Species columns
ASV2 <- ASV[, -c(62:68)] 

#3. Write a csv 

write.csv(ASV2, "/project/egcc/metabarcoding_tarbush_creosote/OTU-tables/2023-OTUs/funguild/FUNGUildCreosoteTarbush.csv")

#Proceed to running FUNGuild on command line using the 'FUNGUildCreosoteTarbush.csv' file as .txt
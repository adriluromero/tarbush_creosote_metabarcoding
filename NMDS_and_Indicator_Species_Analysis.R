#Jornada_Creosote_Tarbush_Fungal_Metabarcoding

#Protocol for processing Metabarcoding of soil
#File name: NMDS_and_Indicator_Species_Analysis

#In this script, you will do a:
#1. Non-metric Multi-dimensional Scaling (NMDS) Analysis
#2. Indicator Species Analysis by harvest and shrub type.  


#Clear environment
rm(list = ls(all.names = TRUE))

#Load packages
library(vegan) #For the NMDS
library(funrar) #For calculating relative abundance
require(indicspecies) #For indicator species analysis
library(ggplot2) #For plotting
require(tidyverse)
library(dplyr)

# 1. Non-metric Multi-dimensional Scaling (NMDS) Analysis

#Import the 'creosote.tarbush.otu.table.csv' file. You can find it in the repository.
otu.table <- read.csv("/project/egcc/metabarcoding_tarbush_creosote/OTU-tables/2023-OTUs/ForStats/creosote.tarbush.otu.table.csv", stringsAsFactors = FALSE)

#Select the columns with ASV data needed for this step.
otu.table2 <- dplyr::select(otu.table, c(1:61))

#Transpose the data to have sample names on rows.
abund_table<-as.data.frame(t(otu.table2))

#Assign names to columns.
names(abund_table) <- abund_table[1,]

#Remove first row. 
abund_table <- abund_table[-1, ] # Delete row 1.
str(abund_table) #Values are character but need numeric.

#Make the values into numeric.
abund_table2 <- mutate_all(abund_table, function(x) as.numeric(as.character(x)))
str(abund_table2) #Now, values are numeric.

#Import the samples data. You can find it in the repository as 'samples_file_CT_LJ.csv'. 
samples <- read.csv("/fs1/project/egcc/metabarcoding_tarbush_creosote/OTU-tables/2023-OTUs/ForStats/samples_file_CT_LJ.csv")

#Format to have sample as rows.
sdata2 <- samples %>%
  column_to_rownames(var = "Sample")

#Merge by row names. 
abund_table_3 <- merge(sdata2, abund_table2, by =0, all=TRUE)

#Reformat so samples are rownames.
abund_table_4 <- abund_table_3 %>%
  column_to_rownames(var = "Row.names")

#Run the NMDS

#Select the columns with the ASV data and make it a matrix.
com2 <- abund_table_4[,5:8568]
str(com2)

m_com2 <- as.matrix(com2)

str(m_com2)

#Calculate relative abundances
rel_com2 <- make_relative(m_com2)

#do the NMDS
set.seed(100)
nmds <- metaMDS(rel_com2, distance = "bray", trymax = 2000, k=3, autotransform = F, noshare = TRUE)

nmds

#Call:
#metaMDS(comm = rel_com2, distance = "bray", k = 3, trymax = 2000, autotransform = F, noshare = TRUE) 

#global Multidimensional Scaling using monoMDS

#Data:     rel_com2 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1380031 
#Stress type 1, weak ties
#Two convergent solutions found after 386 tries
#Scaling: centring, PC rotation 
#Species: expanded scores based on 'rel_com2'  

plot(nmds)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))

#add columns to data frame 
data.scores$Shrub = abund_table_4$Shrub
data.scores$Harvest = abund_table_4$Harvest

head(data.scores)

#Re-level the Harvest levels into June, September, February
data.scores$Harvest <- as.factor(data.scores$Harvest)
str(data.scores)
levels(data.scores$Harvest)
#[1] "Feb-22" "Jun-21" "Sep-21"

#Re-level the harvest levels chronologically
data.scores$Harvest <- factor(data.scores$Harvest, levels = c("Jun-21", "Sep-21", "Feb-22"))
levels(data.scores$Harvest)
#[1] "Jun-21" "Sep-21" "Feb-22"

#Plot the NMDS
NMDS.plot <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Harvest, colour = Shrub))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Shrub Type", y = "NMDS2", shape = "Harvest")  + 
  scale_colour_manual(values=c("#386641", "#a7c957")) 

NMDS.plot

jpeg("/fs1/project/egcc/metabarcoding_tarbush_creosote/Results-Figures/NMDS/NMDS_CT.jpeg",
     height = 4, width = 6, units = 'in', res = 600)
plot(NMDS.plot2)
dev.off()

# 2. Indicator Species Analysis by Harvest and Shrub type.

#For this step, use the 'rel_com2' matrix created above since it has relative abundances.

#Store the values for the factors harvest and shrub type.
harvest <- abund_table_4$Harvest
shrub <- abund_table_4$Shrub

#Run analysis by harvest.
indi_sp_harvest <- multipatt(rel_com2, harvest, func = "r.g", control = how(nperm = 9999))
summary(indi_sp_harvest)

#Multilevel pattern analysis
#---------------------------

#Association function: r.g
#Significance level (alpha): 0.05

#Total number of species: 8564
#Selected number of species: 31 
#Number of species associated to 1 group: 31 
#Number of species associated to 2 groups: 0 

#List of species associated to each combination: 

#Group Feb-22  #sps.  5 
#         stat p.value   
#ASV26   0.352  0.0135 * 
#ASV25   0.337  0.0139 * 
#ASV65   0.317  0.0271 * 
#ASV4258 0.261  0.0244 * 
#ASV851  0.219  0.0061 **

#Group  Jun-21  #sps.  11 
#         stat p.value    
#ASV335  0.453  0.0001 ***
#ASV246  0.390  0.0058 ** 
#ASV857  0.338  0.0140 *  
#ASV571  0.330  0.0068 ** 
#ASV113  0.324  0.0206 *  
#ASV6722 0.307  0.0267 *  
#ASV309  0.295  0.0440 *  
#ASV1139 0.293  0.0334 *  
#ASV2466 0.262  0.0151 *  
#ASV1908 0.256  0.0489 *  
#ASV1450 0.229  0.0316 *  

#Group  Sep-21  #sps.  15 
#          stat p.value    
#ASV2    0.422  0.0011 ** 
#ASV13   0.418  0.0027 ** 
#ASV2096 0.399  0.0003 ***
#ASV122  0.399  0.0019 ** 
#ASV1788 0.395  0.0027 ** 
#ASV1352 0.394  0.0010 ***
#ASV5951 0.390  0.0113 *  
#ASV1067 0.375  0.0108 *  
#ASV1309 0.359  0.0136 *  
#ASV1967 0.343  0.0167 *  
#ASV872  0.327  0.0270 *  
#ASV731  0.319  0.0314 *  
#ASV4431 0.307  0.0366 *  
#ASV2802 0.295  0.0476 *  
#ASV494  0.291  0.0486 *  
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Run analysis by shrub type.
indi_sp_shrub <- multipatt(rel_com2, shrub, func = "r.g", control = how(nperm = 9999))
summary(indi_sp_shrub)

#Multilevel pattern analysis
#---------------------------

#Association function: r.g
#Significance level (alpha): 0.05

#Total number of species: 8564
#Selected number of species: 34 
#Number of species associated to 1 group: 34 

#List of species associated to each combination: 

#Group Creosote  #sps.  15 
#          stat p.value    
#ASV250  0.520  0.0001 ***
#ASV321  0.457  0.0001 ***
#ASV731  0.391  0.0005 ***
#ASV104  0.343  0.0043 ** 
#ASV309  0.340  0.0001 ***
#ASV153  0.329  0.0027 ** 
#ASV598  0.315  0.0074 ** 
#ASV473  0.309  0.0117 *  
#ASV872  0.297  0.0023 ** 
#ASV74   0.289  0.0210 *  
#ASV836  0.279  0.0074 ** 
#ASV987  0.271  0.0144 *  
#ASV274  0.254  0.0350 *  
#ASV1319 0.240  0.0389 *  
#ASV1436 0.188  0.0044 ** 

#Group Tarbush  #sps.  19 
#          stat p.value    
#ASV438  0.410  0.0002 ***
#ASV296  0.361  0.0007 ***
#ASV246  0.359  0.0010 ***
#ASV14   0.352  0.0019 ** 
#ASV705  0.352  0.0017 ** 
#ASV111  0.316  0.0112 *  
#ASV791  0.311  0.0032 ** 
#ASV2527 0.307  0.0023 ** 
#ASV631  0.296  0.0133 *  
#ASV318  0.282  0.0147 *  
#ASV494  0.278  0.0231 *  
#ASV718  0.269  0.0046 ** 
#ASV757  0.254  0.0372 *  
#ASV416  0.253  0.0414 *  
#ASV1004 0.253  0.0423 *  
#ASV244  0.250  0.0322 *  
#ASV1085 0.249  0.0191 *  
#ASV94   0.235  0.0281 *  
#ASV822  0.228  0.0458 * 
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
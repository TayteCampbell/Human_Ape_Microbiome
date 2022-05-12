##########################################################################
### R script used to compare 16S sequencing data of gut microbiomes of ###
### humans chimpanzees and gorillas.                                   ###
### See https://www.nature.com/articles/s41396-020-0634-2              ###
##########################################################################

rm(list=ls())


library("ggplot2")
library("vegan")
library("ape")
library("dplyr")
library("reshape2")

theme_tc <- function() {  # this for all the elements common across plots
  theme_bw() %+replace%
    theme(
      legend.key=element_blank(),
      legend.background = element_rect(color = 'black'),
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      legend.key.size = unit(1.5, 'lines'),
      panel.border = element_rect(color="black",size=1.5, fill = NA),
      
      plot.title = element_text(hjust = 0.05, size = 14),
      axis.text = element_text(size = 20,  color = "black"),
      axis.title = element_text(size = 20, face = "bold", color = "black"),
      
      # formatting for facets
      panel.background = element_blank(),
      strip.background = element_rect(colour="black",size =1.5, fill="grey"), #facet formatting
      panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
      panel.spacing.y = unit(1.5, "lines"), #facet spacing for x axis
      strip.text.x = element_text(size=12, face="bold"), #facet labels
      strip.text.y = element_text(size=12, face="bold", angle = 270) #facet labels
    )
}



species = read.table("relabund_16S_rfile.txt", sep="\t", header=TRUE, row.names=1)

#### Separate sample metadata from sample counts

rownames(species) = species$Captivity
species_data_matrix = species[,4:ncol(species)]
rownames(species_data_matrix) = species$Captivity

NAMES = species$Captivity
species_sample = species[,1:3]
rownames(species_sample)=NAMES


#### PCoA analysis 
bray_distance = vegdist(species_data_matrix, method="bray")
principal_coordinates = pcoa(bray_distance)

pcoa_plot = data.frame(principal_coordinates$vectors[,])
pcoa_plot_merged = merge(pcoa_plot, species_sample, by="row.names")

PCOA = principal_coordinates$values$Eigenvalues


PC1 <- 100*(PCOA[1]/sum(PCOA))
PC2 <- 100*(PCOA[2]/sum(PCOA))

##### Plot PCoA


combined_labels = c("dCongo_Human"="Human - Congo","Waptive_Chimpanzee"="Captive Chimpanzee","Waptive_Gorilla"="Captive Gorilla","Wild_Chimpanzee"="Wild Chimpanzee",
           "Wild_Gorilla"="Wild Gorilla","cElSalvador_Human"="Human - El Salvador","gHadza_Human"="Human - Hadza","bLima_Human"="Human - Lima",
           "aUSA_Human"="Human - USA","Matses_Human"="Human - Matses","fTunapuco_Human"="Human - Tunapuco","eMalawi_Human"="Human - Malawi","fAmazon_Human"="Human - Amerindian",
           "Waptive_Zexternalgorilla"="Captive Gorilla - International")
colors = c("dCongo_Human"="#1f78b4","Waptive_Chimpanzee"="darkorchid1","Waptive_Gorilla"="#33a02c","Wild_Chimpanzee"="#ff7f00","Wild_Gorilla"="#e31a1c",   
           "cElSalvador_Human"="skyblue2","gHadza_Human"="saddlebrown","bLima_Human"="gray28","aUSA_Human"="salmon3","Matses_Human"="aquamarine3","fTunapuco_Human"="darkgoldenrod1",
           "eMalawi_Human"="#c7c7c7","fAmazon_Human"="#ffff3f","Waptive_Zexternalgorilla"="#00b2aa")

ggplot(data=pcoa_plot_merged,aes(x=Axis.1,y=Axis.2)) + 
  geom_point(aes(fill=factor(Combined)),shape=21, color="black", size=6) + 
  theme_tc()+
  scale_fill_manual(values = colors, 
                    labels = combined_labels)+
  labs(x = paste("PC1 - Variation Explained", round(PC1,2),"%"), y = paste("PC2 - Variation Explained", round(PC2,2),"%"))

adonis(species_data_matrix~species_sample$Combined, method="bray",permuations = 999)



#### Distance boxplots - calculate distances
bray_summary = melt(as.matrix(bray_distance))


bray_summary = bray_summary %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)



sd = species %>%
  select(Combined,Captivity) %>%
  mutate_if(is.factor,as.character)

colnames(sd) = c("Combined","Var1")
bray_summarySD = left_join(bray_summary, sd, by = "Var1")

colnames(sd) = c("Combined","Var2")
bray_summarySD = left_join(bray_summarySD, sd, by = "Var2")


#write.csv(df,"all_bray_distances.csv",row.names=FALSE)

#### Plot boxplots

ggplot(subset(bray_summarySD,Combined.x %in% c("aUSA_Human","dCongo_Human","Waptive_Chimpanzee","Waptive_Gorilla","Wild_Chimpanzee","Wild_Gorilla")), aes(Combined.y,value,fill=Combined.y)) +
  geom_boxplot(outlier.color = "black",alpha = 0.9) +
  theme_tc() +
  labs(y="Bray-Curtis Distance",x="") + ylim(0,1)+
  facet_wrap(~Combined.x, labeller = labeller(Combined.x=combined_labels))+
  theme(strip.text.x = element_text(size=9.25))+
  theme(axis.text.x = element_blank())+theme(axis.text=element_text(size=17,color="black"),axis.title=element_text(size=23),
                                             legend.background=element_rect(colour="black"),legend.text = element_text(size=10),legend.title=element_text(size=10))+labs(fill=combined_labels)+theme(legend.title=element_blank())+
  scale_fill_manual(values = colors,
                    labels = combined_labels)


ggplot(bray_summarySD, aes(Combined.y,value,fill=Combined.y)) +
  geom_boxplot(outlier.color = "black",alpha = 0.9) +
  theme_tc() +
  labs(y="Bray-Curtis Distance",x="") + ylim(0,1)+
  facet_wrap(~Combined.x, labeller = labeller(Combined.x=combined_labels))+
  theme(strip.text.x = element_text(size=9.25))+
  theme(axis.text.x = element_blank())+theme(axis.text=element_text(size=17,color="black"),axis.title=element_text(size=23),
                                             legend.background=element_rect(colour="black"),legend.text = element_text(size=10),legend.title=element_text(size=10))+labs(fill=combined_labels)+theme(legend.title=element_blank())+
  scale_fill_manual(values = colors,
                    labels = combined_labels)
 
#### Perform and export statistics

bray_matrix = as.matrix(bray_distance)

bray_matrix[lower.tri(bray_matrix)] = NA


bray_summary = melt(bray_matrix)

bray_summary = na.omit(bray_summary)

bray_summary = bray_summary %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)



sd = species %>%
  select(Combined,Captivity) %>%
  mutate_if(is.factor,as.character)

colnames(sd) = c("Combined","Var1")
bray_summarySD = left_join(bray_summary, sd, by = "Var1")

colnames(sd) = c("Combined","Var2")
bray_summarySD = left_join(bray_summarySD, sd, by = "Var2")


bray_summarySD$comparison = paste(bray_summarySD$Combined.x,bray_summarySD$Combined.y,sep="_")
bray_summarySD

kruskal.test(value~comparison, data=bray_summarySD)

output = pairwise.wilcox.test(bray_summarySD$value, bray_summarySD$comparison, p.adjust.method="BH")


melted = melt(output[[3]])


#write.table(bray_summarySD,"bray_summarySD.txt",sep="\t")



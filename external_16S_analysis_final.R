###### NEW external 16S analysis


#Create PCoA plot of data
rm(list=ls())


library("ggplot2")
library("vegan")
library("ape")
library("dplyr")
library("reshape2")


species = read.table("relabund_16S_rfile.txt", sep="\t", header=TRUE, row.names=1)


#dataPA <- (species[,5:1409] > 0)*1 

#write.table(dataPA, "presence_absence_final.txt", sep = "\t")

rownames(species) = species$Captivity
species_data_matrix = species[,4:592]
rownames(species_data_matrix) = species$Captivity

# 2. PCoA analysis 

adonis(species_data_matrix ~ species$Species*species$Captivity, method="bray", permutations=999)



bray_distance = vegdist(species_data_matrix, method="bray")



#write.table(bray_summary, "melted_bray_distance_matrix.txt", sep = "\t")
hist(bray_distance)
averages = meandist(bray_distance, species$Combined)
heatmap(averages)
#write.table(averages,"bray_distance_averages.txt", sep ="\t")


principal_coordinates = pcoa(bray_distance)

pcoa_plot = data.frame(principal_coordinates$vectors[,])
pcoa_plot_merged = merge(pcoa_plot, species, by="row.names")

# 3. Calculate percent variation explained by PC1, PC2

PC1 <- 100*(principal_coordinates$values$Eigenvalues[1]/sum(principal_coordinates$values$Eigenvalues))
PC2 <- 100*(principal_coordinates$values$Eigenvalues[2]/sum(principal_coordinates$values$Eigenvalues))
PC3 <- 100*(principal_coordinates$values$Eigenvalues[3]/sum(principal_coordinates$values$Eigenvalues))

# 4. Plot PCoA


#legend.position = c(0.857,0.89 ),
#0.14,0.105
ggplot(data=pcoa_plot_merged,aes(x=Axis.1,y=Axis.2)) + geom_point(aes(fill=factor(Combined)),shape=21, colour="black", size=6) + theme_bw()  +
  theme_bw(base_size=20) + 
  theme(axis.text=element_text(size=20,color="black"),axis.title=element_text(size=20),legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=18), legend.title=element_text(size=20)) + labs(fill = "Group")+theme(legend.title=element_blank())+
  labs(x = paste("PC1 - Variation Explained", round(PC1,2),"%"), y = paste("PC2 - Variation Explained", round(PC2,2),"%")) +
scale_fill_manual(values = c("dCongo_Human"="#1f78b4","Waptive_Chimpanzee"="darkorchid1","Waptive_Gorilla"="#33a02c","Wild_Chimpanzee"="#ff7f00","Wild_Gorilla"="#e31a1c",   
                               "cElSalvador_Human"="skyblue2","gHadza_Human"="saddlebrown","bLima_Human"="gray28","aUSA_Human"="salmon3","Matses_Human"="aquamarine3","fTunapuco_Human"="darkgoldenrod1",
                               "other"="black","eMalawi_Human"="#c7c7c7","fAmazon_Human"="#ffff3f","Waptive_Zexternalgorilla"="#00b2aa"), 
                    labels = c("dCongo_Human"="Human - Congo","Waptive_Chimpanzee"="Captive Chimpanzee","Waptive_Gorilla"="Captive Gorilla","Wild_Chimpanzee"="Wild Chimpanzee",
                               "Wild_Gorilla"="Wild Gorilla","cElSalvador_Human"="Human - El Salvador","gHadza_Human"="Human - Hadza","bLima_Human"="Human - Lima",
                               "aUSA_Human"="Human - USA","Matses_Human"="Human - Matses","fTunapuco_Human"="Human - Tunapuco","eMalawi_Human"="Human - Malawi","fAmazon_Human"="Human - Amerindian",
                               "Waptive_Zexternalgorilla"="Captive Gorilla - International"))

  

########Distance boxplots

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


y.labels=c("dCongo_Human"="Human - Congo","Waptive_Chimpanzee"="Captive Chimpanzee","Waptive_Gorilla"="Captive Gorilla","Wild_Chimpanzee"="Wild Chimpanzee",
           "Wild_Gorilla"="Wild Gorilla","cElSalvador_Human"="Human - El Salvador","gHadza_Human"="Human - Hadza","bLima_Human"="Human - Lima",
           "aUSA_Human"="Human - USA","Matses_Human"="Human - Matses","fTunapuco_Human"="Human - Tunapuco","eMalawi_Human"="Human - Malawi","fAmazon_Human"="Human - Amerindian",
           "Waptive_Zexternalgorilla"="Captive Gorilla - International")

ggplot(subset(bray_summarySD,Combined.x %in% c("aUSA_Human","dCongo_Human","Waptive_Chimpanzee","Waptive_Gorilla","Wild_Chimpanzee","Wild_Gorilla")), aes(Combined.y,value,fill=Combined.y)) +
  geom_boxplot(outlier.color = "black",alpha = 0.9) +
  theme_bw() +
  labs(y="Bray-Curtis Distance",x="") + ylim(0,1)+
  facet_wrap(~Combined.x, labeller = labeller(Combined.x=y.labels))+
  theme(strip.text.x = element_text(size=9.25))+
  theme(axis.text.x = element_blank())+theme(axis.text=element_text(size=17,color="black"),axis.title=element_text(size=23),
                                             legend.background=element_rect(colour="black"),legend.text = element_text(size=10),legend.title=element_text(size=10))+labs(fill=y.labels)+theme(legend.title=element_blank())+
  scale_fill_manual(values = c("dCongo_Human"="#1f78b4","Waptive_Chimpanzee"="darkorchid1","Waptive_Gorilla"="#33a02c","Wild_Chimpanzee"="#ff7f00","Wild_Gorilla"="#e31a1c",   
                               "cElSalvador_Human"="skyblue2","gHadza_Human"="saddlebrown","bLima_Human"="gray28","aUSA_Human"="salmon3","Matses_Human"="aquamarine3","fTunapuco_Human"="darkgoldenrod1",
                               "other"="black","eMalawi_Human"="#c7c7c7","fAmazon_Human"="#ffff3f","Waptive_Zexternalgorilla"="#00b2aa"),
                    labels = c("dCongo_Human"="Human - Congo","Waptive_Chimpanzee"="Captive Chimpanzee","Waptive_Gorilla"="Captive Gorilla","Wild_Chimpanzee"="Wild Chimpanzee",
                               "Wild_Gorilla"="Wild Gorilla","cElSalvador_Human"="Human - El Salvador","gHadza_Human"="Human - Hadza","bLima_Human"="Human - Lima",
                               "aUSA_Human"="Human - USA","Matses_Human"="Human - Matses","fTunapuco_Human"="Human - Tunapuco","eMalawi_Human"="Human - Malawi","fAmazon_Human"="Human - Amerindian",
                               "Waptive_Zexternalgorilla"="Captive Gorilla - International"))+
  scale_color_manual(values = c("dCongo_Human"="#1f78b4","Waptive_Chimpanzee"="darkorchid1","Waptive_Gorilla"="#33a02c","Wild_Chimpanzee"="#ff7f00","Wild_Gorilla"="#e31a1c",   
                                "cElSalvador_Human"="skyblue2","gHadza_Human"="saddlebrown","bLima_Human"="gray28","aUSA_Human"="salmon3","Matses_Human"="aquamarine3","fTunapuco_Human"="darkgoldenrod1",
                                "other"="black","eMalawi_Human"="#c7c7c7","fAmazon_Human"="#ffff3f","Waptive_Zexternalgorilla"="#00b2aa"))



ggplot(bray_summarySD, aes(Combined.y,value,fill=Combined.y)) +
  geom_boxplot(outlier.color = "black",alpha = 0.9) +
  theme_bw() +
  labs(y="Bray-Curtis Distance",x="") + ylim(0,1)+
  facet_wrap(~Combined.x, labeller = labeller(Combined.x=y.labels))+
  theme(strip.text.x = element_text(size=9.25))+
  theme(axis.text.x = element_blank())+theme(axis.text=element_text(size=17,color="black"),axis.title=element_text(size=23),
                                             legend.background=element_rect(colour="black"),legend.text = element_text(size=10),legend.title=element_text(size=10))+labs(fill=y.labels)+theme(legend.title=element_blank())+
  scale_fill_manual(values = c("dCongo_Human"="#1f78b4","Waptive_Chimpanzee"="darkorchid1","Waptive_Gorilla"="#33a02c","Wild_Chimpanzee"="#ff7f00","Wild_Gorilla"="#e31a1c",   
                               "cElSalvador_Human"="skyblue2","gHadza_Human"="saddlebrown","bLima_Human"="gray28","aUSA_Human"="salmon3","Matses_Human"="aquamarine3","fTunapuco_Human"="darkgoldenrod1",
                               "other"="black","eMalawi_Human"="#c7c7c7","fAmazon_Human"="#ffff3f","Waptive_Zexternalgorilla"="#00b2aa"),
                    labels = c("dCongo_Human"="Human - Congo","Waptive_Chimpanzee"="Captive Chimpanzee","Waptive_Gorilla"="Captive Gorilla","Wild_Chimpanzee"="Wild Chimpanzee",
                               "Wild_Gorilla"="Wild Gorilla","cElSalvador_Human"="Human - El Salvador","gHadza_Human"="Human - Hadza","bLima_Human"="Human - Lima",
                               "aUSA_Human"="Human - USA","Matses_Human"="Human - Matses","fTunapuco_Human"="Human - Tunapuco","eMalawi_Human"="Human - Malawi","fAmazon_Human"="Human - Amerindian",
                               "Waptive_Zexternalgorilla"="Captive Gorilla - International"))+
  scale_color_manual(values = c("dCongo_Human"="#1f78b4","Waptive_Chimpanzee"="darkorchid1","Waptive_Gorilla"="#33a02c","Wild_Chimpanzee"="#ff7f00","Wild_Gorilla"="#e31a1c",   
                                "cElSalvador_Human"="skyblue2","gHadza_Human"="saddlebrown","bLima_Human"="gray28","aUSA_Human"="salmon3","Matses_Human"="aquamarine3","fTunapuco_Human"="darkgoldenrod1",
                                "other"="black","eMalawi_Human"="#c7c7c7","fAmazon_Human"="#ffff3f","Waptive_Zexternalgorilla"="#00b2aa"))






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



#write.table(bray_summarySD,"bray_summarySD.txt",sep="\t")
#write.table(melted,"statistics_bray_curtis.txt", sep ="\t")



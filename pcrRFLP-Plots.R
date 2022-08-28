#---
#title: "PCR-RFLP"
#author: "Ricardo Ramiro and Augusto Franco"
#date: "16/3/2021"
#---

#Libraries
library(tidyverse)
library(cowplot)
library(ungeviz)
library(vegan)


dat<-read_csv("PCR-RFLP/Final_results-PCR-RFLP.csv")
#this is to remove an unneeded column
dat<-dat %>% select(-X1)
#this is just to make the names smaller
dat$Organism<-str_remove(dat$Organism,"BACTERIASP |Methylorubrum ")
dat$Organism<-str_remove(dat$Organism,"Alphaproteobacteria ")
dat$Organism<-str_replace(dat$Organism,"uncultured","unc.")


dat$Organism_Accession<-paste0(dat$Organism," ",dat$Accession)
# once you have columns for both Accession and organismName, you can do x=paste0(organismName," ",Accession) inside aes()

# probably you need to replace BACTERIASP and Met

ggplot(dat %>% subset(Enzyme=="AluI"),aes(x=Organism_Accession,y=fragments_Size))+
  geom_hpline(width = 0.8, size = 1)+facet_wrap(~Enzyme)+
  theme_cowplot()+
  background_grid(major="xy",color.major = "#f0f0f0")+
  scale_y_continuous(limits=c(0,650),breaks=c(0,seq(100,600,100)))+
  theme(axis.text.x = element_text(angle=90,size=10))+
  xlab("Species Accession")+
  ylab("Fragment size")



ggplot(dat %>% subset(Enzyme=="BstUI"),aes(x=Organism_Accession,y=fragments_Size))+
  geom_hpline(width = 0.8, size = 1)+facet_wrap(~Enzyme)+
  theme_cowplot()+
  background_grid(major="xy",color.major = "#f0f0f0")+
  scale_y_continuous(limits=c(0,400),breaks=c(0,seq(100,400,100)))+
  theme(axis.text.x = element_text(angle=90,size=10))+
  xlab("Species Accession")+
  ylab("Fragment size")

ggplot(dat %>% subset(Enzyme=="CfoI"),aes(x=Organism_Accession,y=fragments_Size))+
  geom_hpline(width = 0.8, size = 1)+facet_wrap(~Enzyme)+
  theme_cowplot()+
  background_grid(major="xy",color.major = "#f0f0f0")+
  scale_y_continuous(limits=c(0,400),breaks=c(0,seq(100,400,100)))+
  theme(axis.text.x = element_text(angle=90,size=10))+
  xlab("Species Accession")+
  ylab("Fragment size")

ggplot(dat %>% subset(Enzyme=="RsaI"),aes(x=Organism_Accession,y=fragments_Size))+
  geom_hpline(width = 0.8, size = 1)+facet_wrap(~Enzyme)+
  theme_cowplot()+
  background_grid(major="xy",color.major = "#f0f0f0")+
  scale_y_continuous(limits=c(0,650),breaks=c(0,seq(100,600,100)))+
  theme(axis.text.x = element_text(angle=90,size=10))+
  xlab("Species Accession")+
  ylab("Fragment size")

ggplot(dat %>% select(Organism_Accession,fragments_Size) %>% distinct(),aes(x=Organism_Accession,y=fragments_Size))+
  geom_hpline(width = 0.8, size = 1)+facet_wrap(~"Pooled enzymes")+
  theme_cowplot()+
  background_grid(major="xy",color.major = "#f0f0f0")+
  scale_y_continuous(limits=c(0,650),breaks=c(0,seq(100,600,100)))+
  theme(axis.text.x = element_text(angle=90,size=10))+
  xlab("Species Accession")+
  ylab("Fragment size")

## Clustering  
#I will be using the jaccard distance as this is simply a presence-absence distance and we don't really care if there are 2 or three bands. I will first do this per enzyme and then for all enzymes at once. When doing this for all enzymes, we are assuming a separate restriction reaction would be run per enzyme.


# AluI
mat_AluI<-dat %>% subset(Enzyme=="AluI") %>% select(Organism_Accession, fragments_Size) %>%
  group_by(Organism_Accession,fragments_Size) %>% summarize(counts=n()) %>%
  pivot_wider(names_from = fragments_Size, values_from = counts,values_fill = 0)

mat_AluI<-as.data.frame(mat_AluI)

rownames(mat_AluI)<-mat_AluI$Organism_Accession
mat_AluI$Organism_Accession<-NULL

clust.res<-hclust(vegdist(mat_AluI,method="jaccard"),method="average")
plot(clust.res,main = "AluI")


# BstUI
mat_BstUI<-dat %>% subset(Enzyme=="BstUI") %>% select(Organism_Accession, fragments_Size) %>%
  group_by(Organism_Accession,fragments_Size) %>% summarize(counts=n()) %>%
  pivot_wider(names_from = fragments_Size, values_from = counts,values_fill = 0)

mat_BstUI<-as.data.frame(mat_BstUI)

rownames(mat_BstUI)<-mat_BstUI$Organism_Accession
mat_BstUI$Organism_Accession<-NULL

clust.res<-hclust(vegdist(mat_BstUI,method="jaccard"),method="average")
plot(clust.res,main = "BstUI")

# CfoI
mat_CfoI<-dat %>% subset(Enzyme=="CfoI") %>% select(Organism_Accession, fragments_Size) %>%
  group_by(Organism_Accession,fragments_Size) %>% summarize(counts=n()) %>%
  pivot_wider(names_from = fragments_Size, values_from = counts,values_fill = 0)

mat_CfoI<-as.data.frame(mat_CfoI)

rownames(mat_CfoI)<-mat_CfoI$Organism_Accession
mat_CfoI$Organism_Accession<-NULL

clust.res<-hclust(vegdist(mat_CfoI,method="jaccard"),method="average")
plot(clust.res,main = "CfoI")

# RsaI
mat_RsaI<-dat %>% subset(Enzyme=="RsaI") %>% select(Organism_Accession, fragments_Size) %>%
  group_by(Organism_Accession,fragments_Size) %>% summarize(counts=n()) %>%
  pivot_wider(names_from = fragments_Size, values_from = counts,values_fill = 0)

mat_RsaI<-as.data.frame(mat_RsaI)

rownames(mat_RsaI)<-mat_RsaI$Organism_Accession
mat_RsaI$Organism_Accession<-NULL

clust.res<-hclust(vegdist(mat_RsaI,method="jaccard"),method="average")
plot(clust.res,main = "RsaI")


# clustering with the 4 enzymes

dat$enzyme_fragSize<-paste0(dat$Enzyme,"_",dat$fragments_Size)


mat_all<-dat %>% select(Organism_Accession, enzyme_fragSize) %>%
  group_by(Organism_Accession,enzyme_fragSize) %>% summarize(counts=n()) %>%
  pivot_wider(names_from = enzyme_fragSize, values_from = counts,values_fill = 0)

mat_all<-as.data.frame(mat_all)

rownames(mat_all)<-mat_all$Organism_Accession
mat_all$Organism_Accession<-NULL

clust.res<-hclust(vegdist(mat_all,method="jaccard"),method="average")
plot(clust.res,main = "all enzymes - separate")

# clustering with the 4 enzymes as if they were all pooled in the same gel
mat_all_pooled<-dat %>% select(Organism_Accession,fragments_Size) %>% distinct() %>%
  group_by(Organism_Accession,fragments_Size) %>% summarize(counts=n()) %>%
  pivot_wider(names_from = fragments_Size, values_from = counts,values_fill = 0)


mat_all_pooled<-as.data.frame(mat_all_pooled)

rownames(mat_all_pooled)<-mat_all_pooled$Organism_Accession
mat_all_pooled$Organism_Accession<-NULL

clust.res<-hclust(vegdist(mat_all_pooled,method="jaccard"),method="average")
plot(clust.res,main = "all enzymes - pooled")




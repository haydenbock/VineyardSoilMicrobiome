#This code computes shannon diversity and tests it using a mixed model
#project repository: https://github.com/SuzanneFleishman/rhizoscale
#publication code created for: (under review, will be updated)

#Code compiled by Suzanne Fleishman; this version saved May 11, 2021
#Based on Borcard, D., Gillet, F., and Legendre, P. (2018). Numerical Ecology with R (Springer). (https://doi.org/10.1007/978-1-4419-7976-6)


#Inputs are Phyloseq objects created in "1PhyloseqObjectCreation"

#To complete this code, run it separately for all phyloseq objectss; select at L41
#to compare bacterial and fungal diversity, see L137
#code needs attention at each step



####STATISTICAL RESULTS SAVED IN THIS FILE####

#### Setup and Data Import ####

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs")

### Clear workspace ###
rm(list=ls())

## load packages we may need
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("gridExtra","nlme","lme4", "effects","corncob","magrittr","car", "interactions","RVAideMemoire","tibble","doParallel","phyloseq", "ggplot2","vegan","microbiome","microbiomeSeq",
              "dplyr","Biostrings","ShortRead","ape","ade4","adegraphics","adespatial","tmap","dada2","jtools","rcompanion","ggnewscale","RColorBrewer","emmeans")
ipak(packages)



##load augmented PS and spatial data ##
load("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/Data&Inputs/Inputs/PS_ccxrs.Rdata")



####Select PS objects for this iteration ####
ps.b<-ps.16s.2019.phy #bacteria
ps.f<-ps.its.2019.phy #fungi

#Subset soil only compartment
ps.soil.b <-subset_samples(ps.b, type=="soil")
ps.soil.f <-subset_samples(ps.f, type=="soil")

    #RAREFY
    sample_sums(ps.soil.b)
    
    ps.soil.b<-prune_samples(sample_sums(ps.soil.b)>0, ps.b)
    
    
    set.seed(500)
    ps.soil.rare.b<-rarefy_even_depth(ps.soil.b)
    sample_sums(ps.soil.rare.b)

#reassign PS object, if it was not rarefied
#ps.rare<-ps

####breakdown ps rare into components####

counts.b<-as.data.frame(ps.soil.rare.b@otu_table)

meta.b<-as.data.frame(ps.soil.rare.b@sam_data)

meta.b$depth<-as.factor(meta.b$depth)

#### Calculate alphadiversity ####
rich.b <- estimate_richness(ps.soil.rare.b,measures="Shannon")
richtable.b<-cbind(meta.b,rich.b)
richtable.b$interact<-paste(richtable.b$gc,richtable.b$rs)

#### Examine for normality ####

quartz()
qqnorm(richtable.b$Shannon, pch = 1, frame = FALSE)
qqline(richtable.b$Shannon, col = "steelblue", lwd = 2)


##############################################################################################################################

#### SECOND APPROACH COMPARTMENT TEST ####

###singularity issues must come from random effect that includes type
###singularity issues arise from interaction of depth*type - assuming bc fescue doesn't have deepest depth

ag_type.b<-aggregate(richtable.b, list(by = meta.b$type,meta.b$depth,meta.b$blk), FUN=mean)

ag_type.b<-cbind(ag_type.b[,1:3],ag_type.b$Shannon)

colnames(ag_type.b)<-c("type","depth","blk","Shannon")

Shan_comp_ag_inter.b=lme(Shannon~ depth, random=~1|blk, data=ag_type.b)
Anova(Shan_comp_ag_inter.b)
###GET ERROR: SINGULARITY IN BACKSOLVE

Shan_comp_ag.b=lme(Shannon~ depth+type, random=~1|blk, data=ag_type.b)
Anova(Shan_comp_ag.b)

emmeans(Shan_comp_ag.b, pairwise~type, adjust="tukey")
emmeans(Shan_comp_ag.b, pairwise~depth, adjust="tukey")


#### grapevine soil and rhizosphere only test ####


richtable_rhizo.b<-richtable.b[richtable.b$type=="grape",]

richtable_soil.b<-richtable.b[richtable.b$type=="soil",]




####Iterate and test grape rhizo

emmeans(Shan_rhizo_3rd.b, pairwise~gc*depth, adjust="tukey")


Shan_rhizo_3rd.b=lme(Shannon~ depth*gc*rs, 
                  random=~1|blk/rs, data=richtable_rhizo.b)
Anova(Shan_rhizo_3rd.b)


Shan_rhizo_2nd.b=lme(Shannon~ depth*gc+depth*rs+gc*rs, 
                  random=~1|blk/rs, data=richtable_rhizo.b)
Anova(Shan_rhizo_2nd.b)

#2019 ASV
#Shan_rhizo_2nd.b=lme(Shannon~ depth*rs+gc*rs, 
                     random=~1|blk/rs, data=richtable_rhizo.b)
#Anova(Shan_rhizo_2nd.b)
#emmeans(Shan_rhizo_2nd.b, pairwise~rs*depth, adjust="tukey")



Shan_rhizo_1st.b=lme(Shannon~ depth+gc+rs, 
                  random=~1|blk/rs, data=richtable_rhizo.b)
Anova(Shan_rhizo_1st.b)

Shan_rhizo_gc.b=lme(Shannon~ gc, 
                 random=~1|blk/rs, data=richtable_rhizo.b)
Anova(Shan_rhizo_gc.b)


####Iterate and test gc soil

Shan_soil_3rd.b=lme(Shannon~ depth*gc*rs, 
               random=~1|blk/rs, data=richtable_soil.b)
Anova(Shan_soil_3rd.b)
emmeans(Shan_soil_3rd.b, pairwise~depth, adjust="tukey")


Shan_soil_2nd.b=lme(Shannon~ depth*gc+depth*rs+gc*rs, 
              random=~1|blk/rs, data=richtable_soil.b)
Anova(Shan_soil_2nd.b)
#emmeans(Shan_soil_2nd.b, pairwise~depth*gc, adjust="tukey")
#emmeans(Shan_soil_2nd.b, pairwise~depth*rs, adjust="tukey")



Shan_soil_1st.b=lme(Shannon~ depth+gc+rs, 
                  random=~1|blk/rs, data=richtable_soil.b)
Anova(Shan_soil_1st.b)
emmeans(Shan_soil_1st.b, pairwise~depth, adjust="tukey")

Shan_soil_gc.b=lme(Shannon~ gc, 
                  random=~1|blk/rs, data=richtable_soil.b)
Anova(Shan_soil_gc.b)




################################################################################
################################################################################
#####REPEAT FOR FUNGI####

#RAREFY
sample_sums(ps.f)

ps.f<-prune_samples(sample_sums(ps.f)>0, ps.f)


set.seed(500)
ps.rare.f<-rarefy_even_depth(ps.f)
sample_sums(ps.rare.f)

#reassign PS object, if it was not rarefied
#ps.rare<-ps

####breakdown ps rare into components####

counts.f<-as.data.frame(ps.rare.f@otu_table)

meta.f<-as.data.frame(ps.rare.f@sam_data)

meta.f$depth<-as.factor(meta.f$depth)


#### Calculate alphadiversity ####
rich.f <- estimate_richness(ps.rare.f,measures="Shannon")
richtable.f<-cbind(meta.f,rich.f)

richtable.f$interact<-paste(richtable.f$gc,richtable.f$rs)


#### Examine for normality ####

quartz()
qqnorm(richtable.f$Shannon, pch = 1, frame = FALSE)
qqline(richtable.f$Shannon, col = "steelblue", lwd = 2)


##############################################################################################################################


#### SECOND APPROACH COMPARTMENT TEST ####

###singularity issues must come from random effect that includes type
###singularity issues arise from interaction of depth*type - assuming bc fescue doesn't have deepest depth

ag_type.f<-aggregate(richtable.f, list(by = meta.f$type,meta.f$depth,meta.f$blk), FUN=mean)

ag_type.f<-cbind(ag_type.f[,1:3],ag_type.f$Shannon)

colnames(ag_type.f)<-c("type","depth","blk","Shannon")

Shan_comp_ag_inter.f=lme(Shannon~ depth*type, random=~1|blk, data=ag_type.f)
Anova(Shan_comp_ag_inter.f)
emmeans(Shan_comp_ag_inter.f, pairwise~depth*type, adjust="tukey")


Shan_comp_ag.f=lme(Shannon~ depth+type, random=~1|blk, data=ag_type.f)
Anova(Shan_comp_ag.f)

emmeans(Shan_comp_ag.f, pairwise~type, adjust="tukey")
emmeans(Shan_comp_ag.f, pairwise~depth, adjust="tukey")


#### grapevine soil and rhizosphere only test ####


richtable_rhizo.f<-richtable.f[richtable.f$type=="grape",]

richtable_soil.f<-richtable.f[richtable.f$type=="soil",]

####Iterate and test grape rhizo


Shan_rhizo_3rd.f=lme(Shannon~ depth*gc*rs, 
                   random=~1|blk/rs, data=richtable_rhizo.f)
Anova(Shan_rhizo_3rd.f)


Shan_rhizo_2nd.f=lme(Shannon~ depth*gc+gc*rs, 
                   random=~1|blk/rs, data=richtable_rhizo.f)
Anova(Shan_rhizo_2nd.f)
emmeans(Shan_rhizo_2nd.f, pairwise~depth, adjust="tukey")



Shan_rhizo_1st.f=lme(Shannon~ depth+gc+rs, 
                   random=~1|blk/rs, data=richtable_rhizo.f)
Anova(Shan_rhizo_1st.f)

Shan_rhizo_gc.f=lme(Shannon~ gc, 
                  random=~1|blk/rs, data=richtable_rhizo.f)
Anova(Shan_rhizo_gc.f)


####Iterate and test gc soil

Shan_soil_3rd.f=lme(Shannon~ depth*gc*rs, 
                  random=~1|blk/rs, data=richtable_soil.f)
Anova(Shan_soil_3rd.f)

#emmeans(Shan_soil_3rd.f, pairwise~depth*gc*rs, adjust="tukey")


emmeans(Shan_soil_3rd.f, pairwise~depth*gc, adjust="tukey")




Shan_soil_2nd.f=lme(Shannon~ depth*gc+gc*rs+rs*depth, 
                  random=~1|blk/rs, data=richtable_soil.f)
Anova(Shan_soil_2nd.f)

emmeans(Shan_soil_2nd.f, pairwise~depth*gc, adjust="tukey")


Shan_soil_1st.f=lme(Shannon~ depth+gc+rs, 
                  random=~1|blk/rs, data=richtable_soil.f)
Anova(Shan_soil_1st.f)

Shan_soil_gc.f=lme(Shannon~ gc, 
                 random=~1|blk/rs, data=richtable_soil.f)
Anova(Shan_soil_gc.f)








################################################################################
################################################################################

####plotsignificant factors####


comp.b<-ggplot(richtable.b,aes(type, Shannon, fill=depth, color=type))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03,
               position=position_dodge(.75),dotsize=.3)+
  scale_fill_manual(values=c("white","snow2","darkgrey"))+
  scale_color_manual(values=c("black", "black","black"))+
  ggtitle("Compartment (16s)")


rhizo.b<-ggplot(richtable_rhizo.b,aes(depth, Shannon, fill=interact,color=interact))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03,
               position=position_dodge(.75),dotsize=.3)+
  scale_fill_manual(values=c("palegreen4","palegreen3","snow4","snow3"))+
  scale_color_manual(values=c("black", "black","black","black"))+
  ggtitle("Rhizosphere (16s)")

soil.b<-ggplot(richtable_soil.b,aes(depth, Shannon, fill=interact,color=interact))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03,
               position=position_dodge(.75),dotsize=.3)+
  scale_fill_manual(values=c("palegreen4","palegreen3","snow4","snow3"))+
  scale_color_manual(values=c("black", "black","black","black"))+
  ggtitle("Soil (16s)") 


comp.f<-ggplot(richtable.f,aes(type, Shannon, fill=depth, color=type))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03,
               position=position_dodge(.75),dotsize=.3)+
  scale_fill_manual(values=c("white","snow2","darkgrey"))+
  scale_color_manual(values=c("black", "black","black"))+
  ggtitle("Compartment (ITS)")


rhizo.f<-ggplot(richtable_rhizo.f,aes(depth, Shannon, fill=interact,color=interact))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03,
               position=position_dodge(.75),dotsize=.3)+
  scale_fill_manual(values=c("palegreen4","palegreen3","snow4","snow3"))+
  scale_color_manual(values=c("black", "black","black","black"))+
  ggtitle("Rhizosphere (ITS)")

soil.f<-ggplot(richtable_soil.f,aes(depth, Shannon, fill=interact,color=interact))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03,
               position=position_dodge(.75),dotsize=.3)+
  scale_fill_manual(values=c("palegreen4","palegreen3","snow4","snow3"))+
  scale_color_manual(values=c("black", "black","black","black"))+
  ggtitle("Soil (ITS)") 




quartz(width = 15, height = 10)
grid.arrange(comp.b, soil.b, rhizo.b, comp.f,soil.f,rhizo.f, ncol=3,nrow=2)





############################
############################
##########################


#boxplot for rhizo and soil
quartz()
ggplot(richtable,aes(type, Shannon, fill=type, color=depth))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03,
               position=position_dodge(.75),dotsize=.3)+
  scale_fill_manual(values=c("white","snow2","darkgrey"))+
  scale_color_manual(values=c("black", "black","black"))




#depthxtype

richtable$typegc<-paste(meta$type,meta$gc)
richtable$typegc[1:35]<-"soil"
quartz()
ggplot(richtable,aes(depth, Shannon, fill=typegc))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03,
               position=position_dodge(.75))+
  scale_fill_manual(values=c("palegreen4", "darkmagenta","snow4","white"))































####fungi vs bacteria####
##run code 2x separately to get a table for its and 16s and assign "richtable" to new dataframe name

richtable.its.gen<-richtable
richtable.16s.gen<-richtable

##merge tables
x<-as.data.frame(richtable.its.gen)#fungi
y<-as.data.frame(richtable.16s.gen) #bacteria

richtable.combined<-merge(x, y,
                          by.x = 'traceid', by.y = 'traceid', all = FALSE)

ggplot(richtable.combined,aes(Shannon.x, Shannon.y, col=cluster.x))+
  geom_point()

ggplot(richtable.combined,aes(agecat.x, Shannon.x, col=colorcat.x))+
  geom_boxplot()

richtable.combined2<-richtable.combined[richtable.combined$Shannon.y!=0,]

reg.sh=lme(Shannon.x~ Shannon.y, random=~1|cluster.x, data = richtable.combined2)
sh.pval<-as.data.frame(anova.lme(reg.sh))#type 3 sums of squares

e<-coef(summary(reg.sh))
e

pi.e<-allEffects(reg.sh,xlevels=20)
pi.e.df<-as.data.frame(pi.e[[1]]) #0.95 default confidence

quartz()
ggplot(pi.e.df,aes(Shannon.y, fit)) + 
  geom_line()+
  geom_point(data=richtable.combined2, aes(Shannon.y, Shannon.x),alpha = 0.6,size=2) +
  geom_ribbon(data=pi.e.df,aes(Shannon.y, ymin = lower, ymax = upper), alpha = .2,show.legend=FALSE,colour = NA)+
  xlab("Bacterial Diversity")+
  ylab("Fungal Diversity")

####test response to ratio ####
richtable.combined$ratio<-richtable.combined$Shannon.y/richtable.combined$Shannon.x
  #higher value = higher bacteria and lower fungi


# load packages 

library(phyloseq)
library(ggplot2)
library(ranacapa)
library(gridExtra)
library(gamm4)


library(plyr)
library(dplyr)

library(ggsci)
library(tidyverse)

library(expss)
library(metagMisc)
library(forcats)
library(ggConvexHull)

library(effects)

library(viridis)
library(ggsci)
library(vegan)

library(ggrepel)

library(glmmTMB)
library("MuMIn")

packageVersion('phyloseq')
packageVersion('ggplot2')

##import data

meerkat_final<-readRDS("geigerdata1.rds")

############# fix some errors

str(sample_data(meerkat_final))

sample_data(meerkat_final)$Symtomatic<-as.character(sample_data(meerkat_final)$Symtomatic)
sample_data(meerkat_final)$Symtomatic<- ifelse(sample_data(meerkat_final)$Symtomatic == "unknown", "Non-symptomatic", sample_data(meerkat_final)$Symtomatic)

table(sample_data(meerkat_final)$Symtomatic)

names(sample_data(meerkat_final))[18]<-"Symptomatic"

names(sample_data(meerkat_final))

###### remove samples with low sequencing depth

sample_data(meerkat_final)$Seq_depth<-phyloseq::sample_sums(meerkat_final)

#### add storage info

available_samples <- read.csv("C:/Users/DFMartinP/Desktop/Bachelorarbeit/ANALYSIS/available_samples.csv", sep=";")

head(available_samples)

sample_data(meerkat_final)$Storage<- vlookup(sample_data(meerkat_final)$SampleNo, available_samples, lookup_column = "Sample", result_column = "STORAGE")

table(sample_data(meerkat_final)$Storage)

sample_data(meerkat_final)$Storage<- ifelse(sample_data(meerkat_final)$Storage == "FREEZEDRIED_BUFFER", "FREEZEDRIED", sample_data(meerkat_final)$Storage)

###### InfectionStatus_Condition_Boxplot ############
###### InfectionStatus_Condition_Boxplot ############
###### InfectionStatus_Condition_Boxplot ############
###### InfectionStatus_Condition_Boxplot ############
###### InfectionStatus_Condition_Boxplot ############

ggplot(sample_data(meerkat_final), aes(x = PCR_pos, y = Condition))+
  geom_boxplot(aes(fill = PCR_pos))+
  geom_jitter(width = 0.2, size = 3, alpha = 0.6)+
  scale_fill_manual(values = c("#00FA9A", "red"))

##################################################### rarefaction curves ########################
##################################################### rarefaction curves ########################
##################################################### rarefaction curves ########################
##################################################### rarefaction curves ########################
##################################################### rarefaction curves ########################

ranacapa::ggrare(meerkat_final, step = 500, se = FALSE)+ 
  # xlim(c(0,50000))+ 
  # facet_wrap(~Species, scales="free_y", ncol=4)+
  geom_vline(xintercept=8000)+
  theme_bw()+
  theme(legend.position = "none")+
  geom_line()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(text=element_text(size=14))+
  ylab("Number of ASVs")+
  xlab("Number of reads")+
  #scale_color_brewer(palette = "Dark2")+
  ggtitle("Rarefaction curves with sequencing depth")

##################################rarefy
##################################rarefy
##################################rarefy
##################################rarefy
##################################rarefy

meerkat_rare<- phyloseq::rarefy_even_depth(meerkat_final, sample.size = 8000, rngseed= 100, replace = TRUE, trimOTUs=TRUE,verbose=TRUE)
meerkat_rare1<- subset_samples(meerkat_rare, IndividID != "3116")
meerkat_rare<-meerkat_rare1
meerkat_rare

sample_sums(meerkat_rare)
sum(sample_sums(meerkat_rare))

############## add points for birth and death dates for each individual
############## add points for birth and death dates for each individual

metadata<-data.frame(sample_data(meerkat_rare)) #extract metadata from phyloseq object as a dataframe
str(metadata)

##for plotting purposes

death_date<-metadata[, c("IndividID","DeathDate")]
names(death_date)[2]<-"SampleDate"

birth_date<-metadata[, c("IndividID","BirthDate")]
names(birth_date)[2]<-"SampleDate"

#SampleTimeline
ggplot(sample_data(meerkat_rare), aes(x = SampleDate, y = IndividID))+
  geom_point(aes(fill=PCR_pos), shape = 21, size = 3)+
  scale_fill_manual(values = c( "#00FA9A", "red"))+
  geom_line(aes(group = IndividID))+
  #theme(text=element_text(size=14))+
  geom_point(data = death_date, col = 'black', shape = 4)+
  geom_point(data = birth_date, col = 'black', shape = 8)+
  ggtitle("Sample timeline, coloured by TB infection state")


###################################### stacked barplots ###########################################################
###################################### stacked barplots ###########################################################
###################################### stacked barplots ###########################################################
###################################### stacked barplots ###########################################################
###################################### stacked barplots ###########################################################
###################################### stacked barplots ###########################################################

#phyloseq::plot_bar(meerkat_final, fill="Phylum")+facet_grid(~PCR_pos, scales = "free")

##make stacked barplot at family level


### fix families with weird double slash

tax_table<- data.frame(meerkat_rare@tax_table@.Data)

tax_table$Family<-gsub("\\\\", "", tax_table$Family)

TAX<- tax_table(tax_table)

row.names(TAX)<- row.names(tax_table)

tax_table(meerkat_rare)<-TAX

head(tax_table(meerkat_rare))

colnames(tax_table(meerkat_rare))<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


#make dataframe that groups ASvs by family

barplot_family <- meerkat_rare %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Family)                                      # Sort data frame alphabetically by family

head(barplot_family)

###########################  edit aesthetic parameters for the object 'barplot_family' ###
###########################  edit aesthetic parameters for the object 'barplot_family' ###
###########################  edit aesthetic parameters for the object 'barplot_family' ###


####

top_families<-ddply(barplot_family, .(Family), summarize, total_abundance=mean(Abundance))
top_families<-top_families[order(-top_families$total_abundance),]
head(top_families, 20)
top_20<-as.character(top_families$Family[1:20])
barplot_family$Family_barplot<-NA


barplot_family$Family<-as.character(barplot_family$Family)


#loop to categorise any rare taxa as 'other' (otherwise way too many to plot)
# this is quite advanced, don't worry about it

for (i in 1:length(top_20)){
  
  barplot_family$Family_barplot<-ifelse(barplot_family$Family == top_20[i], top_20[i], barplot_family$Family_barplot)
}

barplot_family$Family_barplot<-ifelse(is.na(barplot_family$Family_barplot), "Other", barplot_family$Family_barplot)

unique(barplot_family$Family_barplot)

str(barplot_family)

################### colours

#change level order from most abundant to least abundant CHECK

barplot_family$Family_barplot<-factor(barplot_family$Family_barplot, level = c( "Lachnospiraceae", "Clostridiaceae 1",  
                                                                                "Peptostreptococcaceae" , "Enterococcaceae", "Bacteroidaceae", 
                                                                                "Geodermatophilaceae", "Enterobacteriaceae", "Coriobacteriaceae", "Ruminococcaceae", "Erysipelotrichaceae",   
                                                                                "Fusobacteriaceae", "Acidaminococcaceae", "Streptococcaceae", "Bacillaceae", "Beijerinckiaceae", 
                                                                                "Peptococcaceae", "Family XIII", "Planococcaceae",
                                                                                "Eggerthellaceae", "Christensenellaceae", "Other"))
unique(barplot_family$Family_barplot)
## make colour set 

colors <- c("dodgerblue3", "gold", "dimgray", "firebrick", "lightskyblue", "midnightblue",
            "khaki4", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise",
            "mediumvioletred", "mediumorchid4", "turquoise", "turquoise4","violetred1", "violetred4", 
            "wheat", "wheat4", "yellow", "white")

#assign each family a colour

names(colors) <- levels(barplot_family$Family_barplot)
colScale <- scale_fill_manual(name = "Family_barplot",values = colors)

################################# plot
################################# plot
################################# plot
################################# plot

# Reorder following the value of another column:

barplot1<-barplot_family %>%
  mutate(Sample = fct_reorder(Sample, SampleDate)) %>% #reorder samples by date
  ggplot(aes(x = Sample, y = Abundance, fill = Family_barplot)) + 
  geom_bar(stat = "identity", position="fill") +
  theme(axis.text.x = element_blank())+
  ggtitle("Family abundance per sample split by PCR result")+
  #facet_wrap(~IndividID, scales="free")+
  colScale+
  facet_wrap(~PCR_pos, scales ="free")+
  theme(legend.position = "none")


barplot2<-barplot_family %>%
  mutate(Sample = fct_reorder(Sample, SampleNumber)) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Family_barplot)) +
  facet_wrap(~IndividID, scales="free")+
  geom_bar(stat = "identity", position="fill") +
  theme(strip.background = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank()) +
  theme(strip.text = element_blank())+
 ggtitle("Family abundance per sample split by individual")+
 colScale+
  theme(legend.position = "none")

grid.arrange(barplot1, barplot2, ncol =1, heights=1:2)

######################### alpha diversity  ###################################
######################### alpha diversity  ###################################
######################### alpha diversity  ###################################
######################### alpha diversity  ###################################
######################### alpha diversity  ###################################

sample_data(meerkat_rare)$Observed<-phyloseq::estimate_richness(meerkat_rare, measures="Observed")
sample_data(meerkat_rare)$Shannon<-phyloseq::estimate_richness(meerkat_rare, measures=c("Shannon"))
sample_data(meerkat_rare)$Faiths<-metagMisc::phyloseq_phylo_div(meerkat_rare, measures = c("PD")) #another method using metagMisc

str(sample_data(meerkat_rare))

sample_data(meerkat_rare)$Observed<-sample_data(meerkat_rare)$Observed$Observed
sample_data(meerkat_rare)$Shannon<-sample_data(meerkat_rare)$Shannon$Shannon
sample_data(meerkat_rare)$Faiths<-sample_data(meerkat_rare)$Faiths$PD

mean(sample_data(meerkat_rare)$Observed)
max(sample_data(meerkat_rare)$Observed)
min(sample_data(meerkat_rare)$Observed)

####################### visualizing #############
####################### visualizing #############
####################### visualizing #############
####################### visualizing #############

metadata_rare<-data.frame(sample_data(meerkat_rare)) #extract metadata into dataframe

####################### CONDITION 
####################### CONDITION 
####################### CONDITION 
####################### CONDITION 

condition_stats<- glmmTMB(Condition ~
                            
                            PCR_pos+
                            Symptomatic+
                            
                            (1|IndividID), 
                          data =metadata_rare, 
                          family = gaussian)

summary(condition_stats)

plot(allEffects(condition_stats))

########################### boxplots alpha diversity ####################################
########################### boxplots alpha diversity ####################################
########################### boxplots alpha diversity ####################################
########################### boxplots alpha diversity ####################################

str(metadata_rare)

P1<-ggplot(metadata_rare, aes(x = PCR_pos, y = Observed))+
  geom_boxplot(aes(fill = PCR_pos), show.legend = FALSE)+
  geom_jitter(width = 0.2, alpha =0.5)+
 # ggtitle("Observed richness")+
  theme_bw(base_size = 16)+
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("#00FA9A", "red"))

P2<-ggplot(metadata_rare, aes(x = PCR_pos, y = Shannon))+
  geom_boxplot(aes(fill = PCR_pos), show.legend = FALSE)+
  geom_jitter(width = 0.2, alpha =0.5)+
  #ggtitle("Shannon diversity")+
  theme_bw(base_size = 16)+
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("#00FA9A", "red"))

P3<-ggplot(metadata_rare, aes(x = PCR_pos, y = Faiths))+
  geom_boxplot(aes(fill = PCR_pos), show.legend = FALSE)+
  geom_jitter(width = 0.2, alpha =0.5)+
 # ggtitle("Phylogenetic diversity")+
  theme_bw(base_size = 16)+
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("#00FA9A", "red"))


P4<-ggplot(metadata_rare, aes(x = Symptomatic, y = Observed))+
  geom_boxplot(aes(fill = Symptomatic), show.legend = FALSE)+
  geom_jitter(width = 0.2, alpha =0.5)+
 # ggtitle("Observed richness")+
  theme_bw(base_size = 16)+
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("skyblue", "red"))

P5<-ggplot(metadata_rare, aes(x = Symptomatic, y = Shannon))+
  geom_boxplot(aes(fill = Symptomatic), show.legend = FALSE)+
  geom_jitter(width = 0.2, alpha =0.5)+
 # ggtitle("Shannon diversity")+
  theme_bw(base_size = 16)+
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("skyblue", "red"))

P6<-ggplot(metadata_rare, aes(x = Symptomatic, y = Faiths))+
  geom_boxplot(aes(fill = Symptomatic), show.legend = FALSE)+
  geom_jitter(width = 0.2, alpha =0.5)+
  #ggtitle("Phylogenetic diversity")+
  theme_bw(base_size = 16)+
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("skyblue", "red"))

grid.arrange(P1, P2, P3, P4, P5, P6, ncol = 3)

##################### statistical model to predict alpha diversity ##########
##################### statistical model to predict alpha diversity ##########
##################### statistical model to predict alpha diversity ##########
##################### statistical model to predict alpha diversity ##########

####### scale continuous variables around zero to help model converge

metadata_rare<- data.frame(sample_data(meerkat_rare))
str(metadata_rare)

metadata_rare$Seq_depth<-as.numeric(scale(metadata_rare$Seq_depth))
metadata_rare$Faiths<-as.numeric(scale(metadata_rare$Faiths))
metadata_rare$Observed<-as.numeric(scale(metadata_rare$Observed))
metadata_rare$Shannon<-as.numeric(scale(metadata_rare$Shannon))

metadata_rare$DaysTilDeath<-as.numeric(scale(metadata_rare$DaysTilDeath))
metadata_rare$AgeAtSampling<-as.numeric(scale(metadata_rare$AgeAtSampling))
metadata_rare$Weight<-as.numeric(scale(metadata_rare$Weight))
metadata_rare$Condition<-as.numeric(scale(metadata_rare$Condition))

# fixed effect variables = models estimate exact estimates per group
# random effect variables = models estimate an overall effect of group but NOT per group

# generalized linear mixed model 

############# observed richness
############# observed richness
############# observed richness

observed_global<- glmmTMB(Observed ~
                           Seq_depth+
                           Storage+
                           Condition+
                           Sex+
                           AgeAtSampling+
                           DaysTilDeath+
                           PCR_pos+
                          Symptomatic+
                          (1|IndividID), 
                         data =metadata_rare, 
                         family = gaussian)

summary(observed_global)

cor(metadata_rare$Observed, predict(observed_global, type="response"))

AIC(observed_global)

AIC_table_observed<-MuMIn::dredge(observed_global)

head(AIC_table_observed, 10)

#### model averaging because no best model

summary(model.avg(AIC_table_observed), subset = delta < 4)

observed_optimum <- glmmTMB(Observed ~
                              Seq_depth+
                              DaysTilDeath+
                              PCR_pos+
                              Symptomatic+
                              Storage+
                              (1|IndividID), 
                            data =metadata_rare, 
                            family = gaussian)

summary(observed_optimum)
summary(observed_global)

plot(allEffects(observed_optimum))
plot(allEffects(observed_global))

#######Shannon
#######Shannon
#######Shannon

shannon_global<- glmmTMB(Shannon ~
                           Seq_depth+
                           Storage+
                           Condition+
                           Sex+
                           AgeAtSampling+
                           DaysTilDeath+
                           PCR_pos+
                           Symptomatic+
                           (1|IndividID), 
                         data =metadata_rare, 
                         family = gaussian)

summary(shannon_global)

cor(metadata_rare$Shannon, predict(shannon_global, type="response"))

AIC(shannon_global)

AIC_table<-MuMIn::dredge(shannon_global)

head(AIC_table, 10)

shannon_optimum <- glmmTMB(Shannon ~
                             Seq_depth+
                             DaysTilDeath+
                             PCR_pos+
                             Symptomatic+
                             Storage+
                             (1|IndividID), 
                           data =metadata_rare, 
                           family = gaussian)

summary(shannon_optimum)
summary(shannon_global)

plot(allEffects(shannon_optimum))

##########################faiths 
##########################faiths 
##########################faiths 

faiths_global<- glmmTMB(Faiths ~
                            Seq_depth+
                            Storage+
                            Condition+
                            Sex+
                            AgeAtSampling+
                            DaysTilDeath+
                            PCR_pos+
                            Symptomatic+
                            (1|IndividID), 
                          data =metadata_rare, 
                          family = gaussian)

summary(faiths_global)

cor(metadata_rare$Faiths, predict(faiths_global, type="response"))

AIC(faiths_global)

AIC_table<-MuMIn::dredge(faiths_global)

head(AIC_table, 10)

faiths_optimum <- glmmTMB(Faiths ~
                            Seq_depth+
                            DaysTilDeath+
                            PCR_pos+
                            Symptomatic+
                            Storage+
                              (1|IndividID), 
                            data =metadata_rare, 
                            family = gaussian)

summary(faiths_optimum)
summary(faiths_global)

plot(allEffects(faiths_optimum))

######## looking at data ########
######## looking at data ########
######## looking at data ########

sum(sample_sums(meerkat_rare))

sample_variables(meerkat_rare)

table(sample_data(meerkat_rare)$PCR_pos)
table(sample_data(meerkat_rare)$Sex)
min(sample_data(meerkat_rare)$SampleDate)
max(sample_data(meerkat_rare)$SampleDate)
table(sample_data(meerkat_rare)$GroupAtSampling)
table(sample_data(meerkat_rare)$Symptomatic)

ntaxa(meerkat_rare)

# total and mean number of reads in final, rarefied dataset

sum(sample_data(meerkat_rare)$Seq_depth)
mean(sample_data(meerkat_rare)$Seq_depth)
min(sample_data(meerkat_rare)$Seq_depth)
max(sample_data(meerkat_rare)$Seq_depth)

######### table on prevalence and abundance using rarefied dataset
######### table on prevalence and abundance using rarefied dataset
######### table on prevalence and abundance using rarefied dataset
######### table on prevalence and abundance using rarefied dataset
######### table on prevalence and abundance using rarefied dataset
######### table on prevalence and abundance using rarefied dataset

prev0 = apply(X = otu_table(meerkat_rare),
              MARGIN = ifelse(taxa_are_rows(meerkat_rare), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(meerkat_rare),
                    tax_table(meerkat_rare))

prevdf$Prevalence<-(prevdf$Prevalence/362)*100
prevdf$rel_abund<-(prevdf$TotalAbundance/(sum(prevdf$TotalAbundance)))*100

head(prevdf)

#order by prevalence and then abundance

prevdf<-prevdf[order(-prevdf$Prevalence,-prevdf$rel_abund),] #sort by prevalence

prevdf<-prevdf[,c(1,2,10,4,5,6,7,8)] #reorder columns

prevdf<-prevdf[order(-prevdf$rel_abund),] #sort by abundance

head(prevdf, 20)

ggplot(prevdf, aes(x = Prevalence, y = rel_abund))+geom_point()

prevdf25<-subset(prevdf, Prevalence > 25) #only keep taxa over 25% prevalence

## FIGURE
install.packages("ggrepel")

ggplot(prevdf25, aes(x = Prevalence, y = rel_abund, label = Family))+
  geom_point(aes(size = rel_abund, fill = Phylum), shape = 21)+
  geom_text_repel( data  = subset(prevdf25, rel_abund > 1.5 | Prevalence > 85))+ #only label points with high abundance or prevalence
  theme_bw(base_size = 16)+
  scale_size(range = c(3,11), guide = F)+
  theme(legend.position = "right")+
  ggtitle("ASV abundance and prevalence")+
  ylab("Relative abundance (%)")+
  xlab("Prevalence (%)")+
  guides(fill = guide_legend(override.aes = list(size = 7)))+
  scale_fill_locuszoom() #color palette from ggsci

#prevalence
ddply(prevdf, .(Phylum), summarize, MeanPrevalence=mean(Prevalence))

#relative abundance
#phylum
sumabund<-ddply(prevdf, .(Phylum), summarize,SumAbundance=sum(rel_abund))
sumabund<-sumabund[order(-sumabund$SumAbundance),]
head(sumabund, 20)
#family
sumabund<-ddply(prevdf, .(Family), summarize,SumAbundance=sum(rel_abund))
sumabund<-sumabund[order(-sumabund$SumAbundance),]
head(sumabund, 20)

################################################################# beta diversity ################
################################################################# beta diversity ################
################################################################# beta diversity ################
################################################################# beta diversity ################
################################################################# beta diversity ################
################################################################# beta diversity ################
################################################################# beta diversity ################
################################################################# beta diversity ################
################################################################# beta diversity ################

## bray curtis
## unweighted unifrac
## weighted unifrac

#################### plots ######################
#################### plots ######################
#################### plots ######################

############## PCR ###################
############## PCR ###################
############## PCR ###################

###################### Bray-Curtis
###################### Bray-Curtis
###################### Bray-Curtis

meerkat_bray <- phyloseq::distance(meerkat_rare, method = "bray")
betadisperser_PCR <- betadisper(meerkat_bray, sample_data(meerkat_rare)$PCR_pos,
                                bias.adjust = F)
permutest(betadisperser_PCR, pairwise = TRUE, permutations = 10000)

plot(betadisperser_PCR, main = "TB PCR effects Bray-Curtis", cex.lab = 1, cex.axis = 1.3, label = T,
     lwd = 4, ellipse = T, conf = 0.95, hull = F, col = c("#00FA9A", "red"), pch =
       c(19,19), cex=1, seg.lwd = 1)

###################### UniFrac
###################### UniFrac
###################### UniFrac

meerkat_unifrac <- phyloseq::distance(meerkat_rare, method = "unifrac")
betadisperser_PCR <- betadisper(meerkat_unifrac, sample_data(meerkat_rare)$PCR_pos,
                                bias.adjust = F)
permutest(betadisperser_PCR, pairwise = TRUE, permutations = 10000)

plot(betadisperser_PCR, main = "TB PCR effects unweighted UniFrac", cex.lab = 1, cex.axis = 1.3, label = T,
     lwd = 4, ellipse = T, conf = 0.95, hull = F, col = c("#00FA9A", "red"), pch =
       c(19,19), cex=1, seg.lwd = 1)

###################### weighted UniFrac
###################### weighted UniFrac
###################### weighted UniFrac

meerkat_wunifrac <- phyloseq::distance(meerkat_rare, method = "wunifrac")
betadisperser_PCR <- betadisper(meerkat_wunifrac, sample_data(meerkat_rare)$PCR_pos,
                                bias.adjust = F)
permutest(betadisperser_PCR, pairwise = TRUE, permutations = 10000)

plot(betadisperser_PCR, main = "TB PCR effects weighted UniFrac", cex.lab = 1, cex.axis = 1.3, label = T,
     lwd = 4, ellipse = T, conf = 0.95, hull = F, col = c("#00FA9A", "red"), pch =
       c(19,19), cex=1, seg.lwd = 1)

############## Symptoms ###################
############## Symptoms ###################
############## Symptoms ###################

###################### Bray-Curtis
###################### Bray-Curtis
###################### Bray-Curtis

meerkat_bray <- phyloseq::distance(meerkat_rare, method = "bray")
betadisperser_Symptoms <- betadisper(meerkat_bray, sample_data(meerkat_rare)$Symptomatic,
                                bias.adjust = F)
permutest(betadisperser_Symptoms, pairwise = TRUE, permutations = 10000)

plot(betadisperser_Symptoms, main = "TB Symptom effects Bray-Curtis", cex.lab = 1, cex.axis = 1.3, label = T,
     lwd = 4, ellipse = T, conf = 0.95, hull = F, col = c("skyblue", "red"), pch =
       c(19,19), cex=1, seg.lwd = 1)

###################### UniFrac
###################### UniFrac
###################### UniFrac

meerkat_unifrac <- phyloseq::distance(meerkat_rare, method = "unifrac")
betadisperser_Symptoms <- betadisper(meerkat_unifrac, sample_data(meerkat_rare)$Symptomatic,
                                bias.adjust = F)
permutest(betadisperser_Symptoms, pairwise = TRUE, permutations = 10000)

plot(betadisperser_Symptoms, main = "TB Symptom effects unweighted UniFrac", cex.lab = 1, cex.axis = 1.3, label = T,
     lwd = 4, ellipse = T, conf = 0.95, hull = F, col = c("skyblue", "red"), pch =
       c(19,19), cex=1, seg.lwd = 1)

###################### weighted UniFrac
###################### weighted UniFrac
###################### weighted UniFrac

meerkat_wunifrac <- phyloseq::distance(meerkat_rare, method = "wunifrac")
betadisperser_Symptoms <- betadisper(meerkat_wunifrac, sample_data(meerkat_rare)$Symptomatic,
                                bias.adjust = F)
permutest(betadisperser_Symptoms, pairwise = TRUE, permutations = 10000)

plot(betadisperser_Symptoms, main = "TB Symptom effects weighted UniFrac", cex.lab = 1, cex.axis = 1.3, label = T,
     lwd = 4, ellipse = T, conf = 0.95, hull = F, col = c("skyblue", "red"), pch =
       c(19,19), cex=1, seg.lwd = 1)

###################### statistics for beta diversity ############
###################### statistics for beta diversity ############
###################### statistics for beta diversity ############
###################### statistics for beta diversity ############
###################### statistics for beta diversity ############
###################### statistics for beta diversity ############

### stats for beta diversity == PERMANOVA
# TYPE OF NON PARAMETRIC TEST

###################### Bray-Curtis
###################### Bray-Curtis
###################### Bray-Curtis

######### need a distance object, NOT an ordinated object

meerkat_bray <- phyloseq::distance(meerkat_rare, method = "bray")

# make a data frame from the sample_data
metadata_rare_df <- data.frame(sample_data(meerkat_rare))

bray_global<-adonis(meerkat_bray ~ 
                     Seq_depth+
                     Storage+
                     Condition+
                     Sex+
                     AgeAtSampling+
                     DaysTilDeath+
                     PCR_pos+
                     Symptomatic, 
                   data = metadata_rare_df)

bray_global

str(metadata_rare_df)
str(meerkat_unifrac)

###################### UniFrac
###################### UniFrac
###################### UniFrac

meerkat_unifrac <- phyloseq::distance(meerkat_rare, method = "unifrac")

unifrac_global<-adonis(meerkat_unifrac ~
                         Seq_depth+
                         Storage+
                         Condition+
                         Sex+
                         AgeAtSampling+
                         DaysTilDeath+
                         PCR_pos+
                         Symptomatic, 
                    data = metadata_rare_df)
unifrac_global

###################### weighted UniFrac
###################### weighted UniFrac
###################### weighted UniFrac

meerkat_wunifrac <- phyloseq::distance(meerkat_rare, method = "wunifrac")

wunifrac_global<-adonis(meerkat_wunifrac ~
                          Seq_depth+
                          Storage+
                          Condition+
                          Sex+
                          AgeAtSampling+
                          DaysTilDeath+
                          PCR_pos+
                          Symptomatic, 
                       data = metadata_rare_df)
wunifrac_global


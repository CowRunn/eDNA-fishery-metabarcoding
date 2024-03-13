
library(tidyverse)
library(ggpubr)
## uploading the final_data dataset used for the analyses (after threshold filtering)

data_species<-read.csv('final_species.csv')%>% select(!1)%>%
  mutate(primer= ifelse(primer== 'mifish', 'MiFish',
      ifelse(primer== 'coi_NS', 'COI', 
              ifelse(primer== 'fish_NS', 'fish-2kb', 
                ifelse(primer== 'metazoan_NS', 'Metazoan' , NA)))))

## filtering for species detected by at least 3 reads (reads>2)
data_species_4<-data_species%>% filter(reads> 2)%>% select(target, sample, species_name, primer, reads, level_taxon, prop_reads_ass)

## uploading the CSV with species landed/observed on-board
observed<-read.csv('../../../samples_diversity/samples_diversity_logbook_replicates.csv')
obs_catch<-observed%>% select(5, 6, 4, 7)%>% mutate (primer= 'observed_landed') %>% rename(obs = logbook_kg)%>% rename(species_name = sscinames) %>% rename(target = Target_species)
obs_catch<-obs_catch%>% mutate (reads= ifelse(obs== 'observed', 1000, paste(obs)))%>% select(1:3, 5, 6)%>% mutate (level_taxon= 'species')
obs_catch$reads<- as.integer(obs_catch$reads)
str(obs_catch)
obs_catch<-obs_catch%>% mutate(prop_reads_ass= reads)
## merge with observed
data_original_catch<-rbind(data_species_4, obs_catch)

## visualising the species detected by eDNA and landed/observed onboard
ggplot((data_original_catch%>% filter(target== 'Nephrops')), aes(x = primer, y = species_name, fill =log(reads, 10)))+
  geom_tile(color = "black", show.legend = FALSE) +
  facet_grid(~sample)+
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlBu"),na.value = 'light grey')+
  geom_text(aes(label = reads), color = "black", size = 3)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90))+ggtitle('Nephrops samples, reads per species and observed or landed')

ggplot((data_original_catch%>% filter(target!= 'Nephrops')), aes(x = primer, y = species_name, fill =log(reads, 10)))+
  geom_tile(color = "black", show.legend = FALSE) +
  facet_grid(~sample)+
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlBu"),na.value = 'light grey') +
  geom_text(aes(label = reads), color = "black", size = 3)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90))+ggtitle('Roundfish samples, reads per species and observed or landed')

### tidying the data to build barplots and boxplot for data visualisation
final_data_species<- data_original_catch%>%filter(level_taxon== 'species')%>% select(!prop_reads_ass) ## already species

#final_data_species<-final_data_species%>% filter(reads >= 5) ## already filtered reads>3
final_data_pa<-final_data_species%>% pivot_wider(values_from = reads, names_from = primer)%>% mutate(pa_fish=`fish-2kb`, pa_mifish= MiFish,pa_coi=COI, pa_metazoan= Metazoan )
final_data_pa$pa_fish[is.na(final_data_pa$pa_fish)]<-0
final_data_pa$pa_fish[final_data_pa$pa_fish>0]<-1

final_data_pa$pa_metazoan[is.na(final_data_pa$pa_metazoan)]<-0
final_data_pa$pa_metazoan[final_data_pa$pa_metazoan>0]<-1

final_data_pa$pa_mifish[is.na(final_data_pa$pa_mifish)]<-0
final_data_pa$pa_mifish[final_data_pa$pa_mifish>0]<-1

final_data_pa$pa_coi[is.na(final_data_pa$pa_coi)]<-0
final_data_pa$pa_coi[final_data_pa$pa_coi>0]<-1

final_data_pa<-final_data_pa%>% mutate(detected_n= pa_fish+ pa_mifish+ pa_metazoan+ pa_coi)%>% select(1:9, 14)

final_data_pa<-final_data_pa%>%
  mutate(detection_method= ifelse(is.na(observed_landed)& detected_n >0, 'eDNA only',
                                  ifelse(observed_landed>0 & detected_n>0, 'obs onboard + eDNA', 
                                         ifelse(observed_landed>0 & detected_n==0, 'obs onboard only','other'))))

final_data_pa_onboard<-final_data_pa%>% filter(detection_method!= 'eDNA only')%>%mutate(detection= ifelse(detection_method== 'obs onboard + eDNA', 'detected', 'not detected'))

## tidying the data
final_data_pa_prop<-final_data_pa%>%group_by(target, sample, detection_method)%>% count()%>% group_by(sample)%>%mutate(species= sum(n))%>%mutate(prop.species= round(n/species, digits = 4))%>% mutate(perc.species= prop.species*100)
final_data_onboard_prop<-final_data_pa_onboard%>%group_by(target, sample, detection)%>% count()%>% group_by(sample)%>%mutate(species= sum(n))%>%mutate(prop.species= round(n/species, digits = 4))%>% mutate(perc.species= prop.species*100)

## uploading the haul names (see Table S1, for hauls and sample name)
haul_name<-read.csv('../analysis_mock/plots_final/haul_name.csv')%>%select(2:8)
haul_name$Haul_n<-factor(haul_name$Haul_n, levels=haul_name$Haul_n)
haul_name$Haul_l<-factor(haul_name$Haul_l, levels=haul_name$Haul_l)## keep the order
final_data_onboard_prop_haul<-left_join(final_data_onboard_prop, haul_name)
Figure2a<-ggplot(final_data_onboard_prop_haul, aes(x=Haul_n,y= perc.species, fill= detection)) +
  geom_bar (stat = 'identity', width = 0.9, colour='white')+
  scale_fill_manual(values= c( '#BBE189', "#aec3f3"))+#'#2eb8a7'
  facet_grid(( ~ target), scales = 'free_x')+ theme_minimal()+ 
  theme(strip.background = element_rect(colour = 'black', fill = 'white'),
        strip.text.x = element_text(size=10, color = 'black', face = 'italic'),
        strip.text.y = element_text(size=10, color = 'black'),
        axis.text.x = element_text(angle = 90, size=10),
        legend.position="top",
        legend.text = element_text(size=10.2),
        legend.background = element_rect(color= 'lightblue', linewidth =0.5, linetype="solid"))+
  ylab('Species observed onboard (%)')+ xlab(NULL)+ labs(fill= 'eDNA:')

final_data_pa_prop_2<-final_data_pa_prop%>%mutate(detection_method_2= ifelse(detection_method== 'obs onboard + eDNA', 'observed onboard + eDNA',
                                                                             ifelse(detection_method== 'obs onboard only', 'observed onboard only', detection_method)))
final_data_pa_prop_2<-left_join(final_data_pa_prop_2, haul_name)
final_data_pa_prop_2<-final_data_pa_prop_2%>%mutate(y_sum= 100, species_label= paste('(', species, ')'))
Figure2b_label<-final_data_pa_prop_2%>%select(sample, species)%>%unique()%>%mutate(tot= 100)

Figure2b<-ggplot(final_data_pa_prop_2, aes(x=Haul_n,y= perc.species, fill= detection_method_2)) +
  geom_bar (stat = 'identity', width = 0.9, colour='white')+
  scale_fill_manual(values= c( '#df397f','#BBE189', "#aec3f3"))+ #'#2eb8a7'
  facet_grid(( ~ target), scales = 'free_x')+ theme_minimal()+ 
  theme(strip.background = element_rect(colour = 'black', fill = 'white'),
        strip.text.x = element_text(size=10, color = 'black', face = 'italic'),
        strip.text.y = element_text(size=10, color = 'black'),
        axis.text.x = element_text(angle = 90, size=10),
        legend.text = element_text(size=10.2),
        legend.position="top",
        legend.background = element_rect(color= 'lightblue', linewidth =0.5, linetype="solid"))+
  labs(y='Percentage of species', x=NULL, fill= 'Detection:')


figure_2<-ggarrange(Figure2a, Figure2b, ncol = 1, nrow = 2, labels = c('(A)', '(B)'))

ggsave(filename = "Figure2.png", 
       plot = figure_2,
       width = 16, 
       height = 13.2,
       units = "cm",
       dpi = 400)

## Figure S5
FigureS5_suppl<-ggplot(final_data_onboard_prop_haul, aes(x=Haul_n,y= n, fill= detection)) +
  geom_bar (stat = 'identity', width = 0.9, colour='white')+
  scale_fill_manual(values= c( '#BBE189', "#aec3f3"))+#'#2eb8a7'
  facet_grid(( ~ target), scales = 'free_x')+ theme_minimal()+ 
  theme(strip.background = element_rect(colour = 'black', fill = 'white'),
        strip.text.x = element_text(size=10, color = 'black', face = 'italic'),
        strip.text.y = element_text(size=10, color = 'black'),
        axis.text.x = element_text(angle = 90, size=10),
        legend.position="top",
        legend.text = element_text(size=10.2),
        legend.background = element_rect(color= 'lightblue', linewidth =0.5, linetype="solid"))+
  ylab('Number of species observed onboard')+ xlab(NULL)+ labs(fill= 'eDNA:')

FigureS5b_suppl<-ggplot(final_data_pa_prop_2, aes(x=Haul_n,y= n, fill= detection_method_2)) +
  geom_bar (stat = 'identity', width = 0.9, colour='white')+
  scale_fill_manual(values= c( '#df397f','#BBE189', "#aec3f3"))+ #'#2eb8a7'
  facet_grid(( ~ target), scales = 'free_x')+ theme_minimal()+ 
  theme(strip.background = element_rect(colour = 'black', fill = 'white'),
        strip.text.x = element_text(size=10, color = 'black', face = 'italic'),
        strip.text.y = element_text(size=10, color = 'black'),
        axis.text.x = element_text(angle = 90, size=10),
        legend.text = element_text(size=10.2),
        legend.position="top",
        legend.background = element_rect(color= 'lightblue', linewidth =0.5, linetype="solid"))+
  labs(y='Number of species', x=NULL, fill= 'Detection:')

figure_S5<-ggarrange(FigureS5_suppl, FigureS5b_suppl, ncol = 1, nrow = 2, labels = c('(A)', '(B)'))

ggsave(filename = "FigureS5.png", 
       plot = figure_S5,
       width = 16, 
       height = 16,
       units = "cm",
       dpi = 400)

######
# figure 3 PCoA (in beta_diversity_onboard)
## figure 4 (heat map) in figures (4-5)
## figure 5 correlations in figures (4-5)


## Wilcoxon tests
### comparisons between number of species (that were observed on-board) detected in Nephrops and Roundfish samples 
final_data_onboard_prop_detected<-final_data_onboard_prop%>% filter(detection== 'detected')
final_data_onboard_prop_detected$target<- as.factor(final_data_onboard_prop_detected$target)
Wilcox_test_N_R<-wilcox.test(final_data_onboard_prop_detected$n~final_data_onboard_prop_detected$target)
Wilcox_test_N_R

## Species detected by eDNA
## comparison of total number of species detected in the eDNA (not only restricted to those species observed on-board and also detected by eDNA)
n_species_across_samples<-taxa_final_data_species%>% filter(primer!='observed_landed')%>% group_by(sample)%>%distinct(species_name)%>%count %>% rename(species= n)
min(n_species_across_samples$species)
max(n_species_across_samples$species)
mean(n_species_across_samples$species)
sd(n_species_across_samples$species)
n_species_across_samples_target<-left_join(n_species_across_samples, target_species)%>% group_by(target)%>% 
  mutate(mean_sp= mean(species), sd_species= sd(species))

test_sp_target<- wilcox.test(n_species_across_samples_target$species ~ n_species_across_samples_target$target)
test_sp_target

## species observed on-board comparison between target fisheries
obs_n_species_across_samples<-taxa_final_data_species%>% filter(primer=='observed_landed')%>% group_by(sample)%>%distinct(species_name)%>%count %>% rename(species= n)
obs_n_species_across_samples<-left_join(obs_n_species_across_samples, target_species)

obs_n_species_across_samples$sample<-gsub("_.*", "", obs_n_species_across_samples$sample)
obs_n_species_across_samples<-obs_n_species_across_samples%>% distinct()
obs_n_species_across_samples_target<-obs_n_species_across_samples%>% group_by(target)%>% mutate(mean_sp= mean(species), sd_species= sd(species))

test_sp_target_obs<- wilcox.test(obs_n_species_across_samples_target$species ~ obs_n_species_across_samples_target$target)
test_sp_target_obs


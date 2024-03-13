## supplementary figures
library(tidyverse)
library(FSA)
## they are based on datasets created in the script figures_manuscript

## Figure S2
taxonomy<-read.csv('../analysis_mock/midori_NS_database/taxonomy.csv')%>% select(!1)
mullet<- c('Chordata', 'Actinopteri', 'Mulliformes', 'Mullidae', 'Mullus', 'Mullus barbatus')
taxonomy_2<- rbind(taxonomy, mullet)
taxa_final_data_species<-left_join(final_data_species, taxonomy_2)%>% filter(species_name!= 'had-whit',  species_name!='dab-Ampla')%>%droplevels()
taxa_final_data_species_primers<- taxa_final_data_species%>% filter(primer!= 'observed_landed')

unique(taxa_final_data_species_primers$Phylum)
unique(taxa_final_data_species_primers$Class)
unique(taxa_final_data_species_primers$species_name)
unique(taxa_final_data_species_primers$Genus)
unique(taxa_final_data_species_primers$Family)
unique(taxa_final_data_species_primers$Order)

n_species<-taxa_final_data_species%>% group_by(sample, primer)%>% count %>% rename(species= n)
n_genus<-taxa_final_data_species%>% group_by(sample, primer)%>% count(Genus)%>% count()%>% rename(genus= n)
n_family<-taxa_final_data_species%>% group_by(sample, primer)%>% count(Family)%>% count()%>% rename(family= n)
n_Order<-taxa_final_data_species%>% group_by(sample, primer)%>% count(Order)%>% count()%>% rename(order= n)
n_Class<-taxa_final_data_species%>% group_by(sample, primer)%>% count(Class)%>% count()%>% rename(class= n)

n_taxonomy<-left_join(n_species, n_genus)
n_taxonomy<-left_join(n_taxonomy, n_family)
n_taxonomy<-left_join(n_taxonomy, n_Order)
n_taxonomy<-left_join(n_taxonomy, n_Class)

n_taxonomy_long<-pivot_longer(n_taxonomy, 3:7, names_to = 'level_taxonomy', values_to = 'numbers') 
n_taxonomy_long$level_taxonomy<- factor(n_taxonomy_long$level_taxonomy, levels = c("class", "order","family","genus", "species"))
FigureS2<-ggplot(n_taxonomy_long, aes(x= primer, y= numbers, fill=level_taxonomy))+
  geom_boxplot()+
  scale_fill_brewer(palette = 'Set2')+
  facet_grid(~primer, scales = "free_x", space = 'free')+
  theme_bw()+
  theme(strip.background = element_rect(colour = 'black', fill = 'white'),
        strip.text.x = element_text(size=10, color = 'black'),
        strip.text.y = element_text(size=10, face="bold", color = 'black'),
        axis.text.x = element_blank(),
        legend.position="top",
        legend.background = element_rect(color= 'lightblue', linewidth =0.5, linetype="solid"))+
  labs (y='Taxa detected across samples', fill= 'Taxonomyc rank', x=NULL)

ggsave(filename = "FigureS2.png", 
       plot = FigureS2,
       width = 22, 
       height = 8,
       units = "cm",
       dpi = 900)

## Figure S3
final_data<-data_original_catch
samples_fish<-c("sample_538_2", "sample_538_3", "sample_544_1", "sample_544_2", "sample_544_3", "sample_549_1", "sample_549_2", "sample_549_3", "sample_569_1", "sample_569_2", "sample_569_3", "sample_575_1",
                "sample_580_1", "sample_585_1")

fish_no_final_data<- final_data%>% filter(primer!= 'observed_landed', primer!= 'fish-2kb')%>% droplevels()%>% group_by(sample)%>% distinct(species_name)
no_sp_no_fish<-as.data.frame(table(fish_no_final_data$sample))%>% rename(sample= Var1, Species= Freq)%>% mutate(sample= paste('sample', sample, sep = '_'))%>% mutate(dataset= 'fish')

fish_final_data<- final_data%>% filter(primer!= 'observed_landed')%>% droplevels()%>% group_by(sample)%>% distinct(species_name)
no_sp_fish<-as.data.frame(table(fish_final_data$sample))%>% rename(sample= Var1, Species= Freq)%>% mutate(sample= paste('sample', sample, sep = '_'))%>% mutate(dataset= 'fish_included')


######################## new analyses based on comments (AP)
no_sp_fish_new<-no_sp_fish%>% rename(species_all= Species)%>% select(1:2)
no_sp_no_fish_new<-no_sp_no_fish%>% rename(Species_detected= Species)%>% select(1:2)

fish_dataset_new<-left_join(no_sp_fish_new, no_sp_no_fish_new)
fish_dataset_new<-fish_dataset_new%>% mutate(difference= species_all - Species_detected)%>% mutate(dataset= 'fish-2kb')
fish_dataset_new<-fish_dataset_new%>% filter(sample%in%samples_fish)

## mean of unique species detected across samples by Fish-2kb marker
mean(fish_dataset_new$difference)
sd(fish_dataset_new$difference)

mean(fish_dataset_new$species_all)
sd(fish_dataset_new$species_all)

mean(fish_dataset_new$Species_detected)
sd(fish_dataset_new$Species_detected)

##unique species detected by metazoan
metazoan_no_final_data<-data_original_catch%>% filter(primer!= 'observed_landed', primer!= 'Metazoan')%>% droplevels()%>% group_by(sample)%>% distinct(species_name)
metazoan_no_final_data_2<-as.data.frame(table(metazoan_no_final_data$sample))%>% rename(sample= Var1, Species_detected= Freq)%>% mutate(sample= paste('sample', sample, sep = '_'))%>% mutate(dataset= 'Metazoan')
unique_metazoan<- left_join(metazoan_no_final_data_2, no_sp_fish_new)%>% mutate(difference= species_all - Species_detected)
mean(unique_metazoan$difference)
sd(unique_metazoan$difference)

## unique species detected by COI
coi_no_final_data<-final_data%>% filter(primer!= 'observed_landed', primer!= 'COI')%>% droplevels()%>% group_by(sample)%>% distinct(species_name)
coi_no_final_data_2<-as.data.frame(table(coi_no_final_data$sample))%>% rename(sample= Var1, Species_detected= Freq)%>% mutate(sample= paste('sample', sample, sep = '_'))%>% mutate(dataset= 'COI')
unique_coi<- left_join(coi_no_final_data_2, no_sp_fish_new)%>% mutate(difference= species_all - Species_detected)
mean(unique_coi$difference)
sd(unique_coi$difference)

##unique species detected by mifish
mifish_no_final_data<-final_data%>% filter(primer!= 'observed_landed', primer!= 'MiFish')%>% droplevels()%>% group_by(sample)%>% distinct(species_name)
mifish_no_final_data_2<-as.data.frame(table(mifish_no_final_data$sample))%>% rename(sample= Var1, Species_detected= Freq)%>% mutate(sample= paste('sample', sample, sep = '_'))%>% mutate(dataset= 'MiFish')
unique_mifish<- left_join(mifish_no_final_data_2, no_sp_fish_new)%>% mutate(difference= species_all - Species_detected)
mean(unique_mifish$difference)
sd(unique_mifish$difference)

## merging the datasets with unique species detected by each marker
unique_data<-rbind( fish_dataset_new, unique_coi, unique_metazoan, unique_mifish)

figureS3<-ggplot(unique_data, aes(x=dataset, y=difference)) +
  geom_boxplot(fill="#aec3f3") +
  #facet_grid(~ dataset)+
  labs(title= '', x= 'Primer set ', y= 'Number of unique species detected') +
  geom_point()+ theme_light()+theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size = 12))

ggsave(filename = "figureS3.png", 
       plot = figureS3,
       width = 16.6, 
       height = 14.4,
       units = "cm",
       dpi = 400)

### Figure S4
final_data_pa_longer<- pivot_longer(final_data_pa, 5:9, names_to = 'primer', values_to = 'reads')%>% drop_na(reads)%>% filter(primer!= 'observed_landed')

final_data_pa_longer$detection_method<- factor(final_data_pa_longer$detection_method, levels=c('eDNA only', 'obs onboard + eDNA'))
figureS4<-ggplot(final_data_pa_longer, aes(x= sample, fill= detection_method))+
  geom_bar(stat = 'count', width = 1, colour='white', position = position_stack())+
  #scale_fill_brewer(palette = 'Set2')+
  #scale_fill_manual(values = c("#E7B800","")colors are here #a6c6d7,#6793e6)+ 
  scale_fill_manual(values = c("#aec3f3","#0146ab"))+ 
  facet_grid(primer~sample, scales = "free_x", space = 'free')+
  theme_bw()+
  theme(strip.background = element_rect(colour = 'black', fill = 'white'),
        strip.text.x = element_text(size=10, color = 'black'),
        strip.text.y = element_text(size=10, face="bold", color = 'black'),
        axis.text.x = element_blank(),
        legend.position="top",
        legend.background = element_rect(color= 'lightblue', linewidth =0.5, linetype="solid"))+
  labs (y='Number of Species', fill= 'Detection', x='Sample')

ggsave(filename = "FigureS4.png", 
       plot = figureS4,
       width = 24.8, 
       height = 14,
       units = "cm",
       dpi = 400)

species_detected<-table(final_data_pa_longer$sample, final_data_pa_longer$primer, final_data_pa_longer$detection_method)

## mean and SD of number of species detected in the eDNA only (not observed on-board) for the 18 samples
final_data_pa_longer_eDNAonly<-final_data_pa_longer%>% filter(detection_method== 'eDNA only')%>% group_by(sample, primer)%>% count()#%>% pivot_wider(values_from = n, names_from = primer)
Summarize(n~ primer, data=final_data_pa_longer_eDNAonly)

## mean and SD of number of species detected in eDNA and observed on board (detected and observed)
final_data_pa_longer_eDNAonboard<-final_data_pa_longer%>% filter(detection_method!= 'eDNA only')%>% group_by(sample, primer)%>% count()#%>% #pivot_wider(values_from = n, names_from = primer)
Summarize(n~ primer, data=final_data_pa_longer_eDNAonboard)

# figure S5 in the script " figures_manuscript" at line 119

## figure S6 (script= replicates_analysis.R)
## figure S7 (script= replicates_analysis.R)

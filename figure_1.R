## figure 1
library(tidyverse)

observed_2<-obs_catch%>% mutate(obs_landed= ifelse(reads== '1000', '*', reads))%>%
  mutate(sample= paste('sample', sample, sep = '_'))%>% select(2,3,8)
ph_data_species_w<-ph_data_otu %>% select(!OTU)%>% pivot_wider(names_from = species_name, values_from = reads_sum)
ph_data_species_l<-pivot_longer(ph_data_species_w, 6:57,  names_to = 'species_name', values_to = 'reads_sum')
obs_ph_data_otu<- left_join(ph_data_species_l, observed_2)
obs_ph_data_otu$species_name<-as.character(obs_ph_data_otu$species_name) 
obs_ph_data_otu$species_name<-factor(obs_ph_data_otu$species_name, levels = (unique(rev(sort(obs_ph_data_otu$species_name)))))
obs_ph_data_otu<-obs_ph_data_otu%>% mutate (sample_id= sample)
obs_ph_data_otu$sample_id<- str_remove(obs_ph_data_otu$sample_id, 'sample_')
hauls<- read.csv('../analysis_mock/plots_final/haul_modified_sample_logbook_depths2.csv')%>% select(2,18)
obs_ph_data_otu<-left_join(obs_ph_data_otu, hauls)
obs_ph_data_otu$Haul_n<-factor(obs_ph_data_otu$Haul_n, levels = c("haul_2","haul_3-1","haul_3-2", "haul_3-3", "haul_4-1", "haul_4-2", "haul_4-3", "haul_5-1",
                                                                  "haul_5-2" ,"haul_5-3" ,"haul_6",   "haul_7" ,  "haul_9-1", "haul_9-2", "haul_9-3", "haul_10",  "haul_11" , "haul_12"))

## plot 1
obs_ph_data_otu_taxonomy<-left_join(obs_ph_data_otu, taxonomy_2)
species_lifestyle<-obs_ph_data_otu_taxonomy%>%drop_na(reads_sum)%>%mutate(detection= ifelse(is.na(obs_landed), 'eDNA-only' , 'eDNA + observed'))
life<-read.csv('../analysis_mock/network_analyses/fishbase_lifestyle2.csv')%>%select(species_name, life)
species_lifestyle<-left_join(species_lifestyle, life)

species_lifestyle%>%filter(detection== 'eDNA-only')%>%group_by(life)%>%count()
## 226 total species eDNA only, 12 pelagic (bathypelagic, pelagic nerictic and pelagic oceanic)
## 454-226= 228 eDNA and observed

12/226 ## 5% of species detected by eDNA only are pelagic

figure1_lifestyle_plot<-ggplot(species_lifestyle, aes(x=Haul_n,  fill= life)) +
  geom_bar (stat = 'count', width = 0.9, colour='black')+
  scale_fill_brewer(palette = 'Set2')+
  facet_grid((detection~target), scales = 'free_x')+ theme_minimal()+ 
  theme(legend.position="bottom",
        legend.background = element_rect(color= 'lightblue', linewidth =0.5, linetype="solid"),
        strip.background = element_rect(colour = 'lightblue',linewidth =0.2, fill = 'white'),
        strip.text.x = element_text(size=10, color = 'black', face = 'italic'),
        strip.text.y = element_text(size=10, color = 'black'),
        axis.text.x = element_text(angle = 90))+
  ggtitle('')+labs(y='Number of species', x=NULL, fill= 'Category')


### Map and final figure in the Script map.R

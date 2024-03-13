## libraries

library(tidyverse)

## Figure 4 with reads sum (across primer sets)
Figure4<-ggplot(obs_ph_data_otu_taxonomy, aes(x = Haul_n, y = species_name, fill = log(reads_sum, 10)))+
  geom_tile(color = "black") + #, show.legend = FALSE
  facet_grid(Order~target, scales='free', space = 'free')+
  #scale_fill_gradientn(colors = hcl.colors(20, "RdYlBu"),na.value = 'white') +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, na.value = 'white' )+
  geom_text(aes(label = obs_landed), color = "black", size = 2.4)+xlab(NULL)+ylab(NULL)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(strip.background = element_rect(colour = 'black', fill = 'white', linewidth = 0.2),
        strip.text.x = element_text(size=10,  color = 'black', face= 'italic'),
        strip.text.y = element_text(size=8, angle = 0,color = 'black'),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size= 10, face = 'italic'),
        plot.tag.position = 'bottom',
        #legend.position="bottom",
        plot.caption = element_text(hjust = 0.25))+ labs(fill= 'log(reads)', caption  = '*=observed onboard; n= kg landed')

# bottom= width 18, new= width 20
ggsave(filename = "figure4.png", 
       plot = Figure4,
       width = 18, 
       height = 30,
       units = "cm",
       dpi = 500)

## adding proportion of reads for each species detected in the eDNA (sum of abundances across primer sets)
data_original_w<-data_original_catch %>%select(!reads)%>%filter(primer!= 'observed_landed')%>%pivot_wider(names_from = primer, values_from = prop_reads_ass)%>%
  mutate_at(c('MiFish', 'Metazoan', 'fish-2kb', 'COI'), ~replace_na(.,0))
data_original_l<-data_original_w%>%group_by(target, sample, species_name)%>% mutate(sum_abundance= (MiFish+ Metazoan+ `fish-2kb`+ COI))
data_original_l<-data_original_l%>%select(sample, target, species_name, sum_abundance)%>% pivot_wider(names_from = species_name, values_from = sum_abundance)
data_original_l<-pivot_longer(data_original_l, 3:54,  names_to = 'species_name', values_to = 'sum_abundance')


observed_abundance<- observed_2
observed_abundance$sample<-str_remove(observed_2$sample, 'sample_') # remove string from a character

data_original_abundance<-left_join(data_original_l, observed_abundance)
data_original_abundance<-left_join(data_original_abundance, haul_name)
data_original_abundance<-left_join(data_original_abundance, taxonomy_2)

figure6_abundance<-ggplot(data_original_abundance, aes(x = Haul_n, y = species_name, fill = log(sum_abundance, 10)))+#
  geom_tile(color = "black") + #, show.legend = FALSE
  facet_grid(Order~target, scales='free', space = 'free')+
  #scale_fill_gradientn(colors = hcl.colors(20, "RdYlBu"),na.value = 'white') +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, na.value = 'white' )+
  geom_text(aes(label = obs_landed), color = "black", size = 2.4)+xlab(NULL)+ylab(NULL)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(strip.background = element_rect(colour = 'black', fill = 'white'),
        strip.text.x = element_text(size=10,  color = 'black', face= 'italic'),
        strip.text.y = element_text(size=8, angle = 0,color = 'black'),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(face = 'italic'),
        plot.tag.position = 'bottom',
        #legend.position="bottom",
        plot.caption = element_text(hjust = 0.25))+ labs(fill= 'log(abundance)', caption  = '*=observed onboard; n= kg landed')

ggsave(filename = "figure6_abundance_rigth_legend.png", 
       plot = figure6_abundance,
       width = 20, 
       height = 30,
       units = "cm",
       dpi = 900)


### Figure 5
summary_landed_species<-read.csv('../analysis_mock/plots_final/subsampling_summary.csv')%>%filter(origin== 'landed')
summary_landed_species<-left_join(summary_landed_species, taxonomy)
#summary_landed_species<-summary_landed_species%>%select(!sample)%>%rename(sample= sample_id)

quantitative_analyses<-data_original_catch%>% select (2,3, 4, 7)%>% pivot_wider(names_from = primer, values_from = prop_reads_ass)
quantitative_analyses<-quantitative_analyses%>%pivot_longer(3:6, names_to = 'primer', values_to = 'prop_reads')
quantitative_analyses<-quantitative_analyses%>%mutate(obs=ifelse(observed_landed== '1000', 'obs', 'landed'))
quantitative_analyses$obs[is.na(quantitative_analyses$obs)]<- 'not obs'

quantitative_landed<-left_join(summary_landed_species, quantitative_analyses)
quantitative_landed<-left_join(quantitative_landed, taxonomy)
label_quantitative<-read.csv('../analysis_mock/plots_final/quantitative_R.csv')
quantitative_landed_h<-quantitative_landed%>%filter(species_name %in% c('Gadus morhua', 'Lophius piscatorius', 'Melanogrammus aeglefinus', 'Glyptocephalus cynoglossus', 'Nephrops norvegicus'))
figure5<-ggplot((quantitative_landed_h), aes(x = n, y = prop_reads))+ #, color=species_name
  geom_point(size=2) + 
  geom_smooth(method = 'lm')+#, show.legend = FALSE
  facet_grid(species_name~primer, scales='free')+
  geom_text(data = label_quantitative,
            mapping= aes(x = Inf, y = Inf, label = paste('p=', p.value)), size=3,
            hjust= 1.16,
            vjust= 2.8)+
  geom_text(data = label_quantitative,
            mapping= aes(x = Inf, y = Inf, label = paste('R^2=', R2)), size=3,
            hjust= 1.16,
            vjust= 1.5)+
  #scale_color_brewer(palette = "Set3")+
  theme_light()+
  theme(strip.background = element_rect(colour = 'black', fill = 'white'),
        strip.text.x = element_text(size=12, color = 'black'),
        strip.text.y = element_text(size=10,angle = 0, color = 'black', face= 'italic'),
        axis.text.y = element_text(size= 12),
  )+labs(x= 'kg landed', y='Relative Read Abundance')

figure5
ggsave(filename = "figure_5.png", 
       plot = figure5,
       width = 26, 
       height = 15,
       units = "cm",
       dpi = 400)

# correlation used to make the label_quantitative dataset
### correlation (COI, not significant)
cod_coi<-quantitative_landed%>%filter(species_name== 'Gadus morhua', primer== 'COI')
coi_coi_lm<-lm(cod_coi$prop_reads ~ cod_coi$n)
summary(coi_coi_lm)

cod_coi<-quantitative_landed%>%filter(species_name== 'Gadus morhua', primer== 'Metazoan')
coi_coi_lm<-lm(cod_coi$prop_reads ~ cod_coi$n)
summary(coi_coi_lm)

cod_coi<-quantitative_landed%>%filter(species_name== 'Gadus morhua', primer== 'MiFish')
coi_coi_lm<-lm(cod_coi$prop_reads ~ cod_coi$n)
summary(coi_coi_lm)


# Lophius all significant
Lophius_coi<-quantitative_landed%>%filter(species_name== 'Lophius piscatorius', primer== 'COI')
Lophius_coi_lm<-lm(Lophius_coi$prop_reads ~ Lophius_coi$n)
summary(Lophius_coi_lm)
par(mfrow = c(2, 2))
plot(Lophius_coi_lm) ## no 

Lophius_fish<-quantitative_landed%>%filter(species_name== 'Lophius piscatorius', primer== 'fish-2kb')%>% drop_na(prop_reads)
Lophius_fish_lm<-lm((Lophius_fish$prop_reads) ~ Lophius_fish$n)
summary(Lophius_fish_lm)
par(mfrow = c(2, 2))
plot(Lophius_fish_lm)


Lophius_mifish<-quantitative_landed%>%filter(species_name== 'Lophius piscatorius', primer== 'MiFish')
Lophius_mifish_lm<-lm(Lophius_mifish$prop_reads ~ Lophius_mifish$n)
summary(Lophius_mifish_lm)
par(mfrow = c(2, 2))
plot(Lophius_mifish_lm) ##no

Lophius_metazoan<-quantitative_landed%>%filter(species_name== 'Lophius piscatorius', primer== 'Metazoan')
Lophius_metazoan_lm<-lm(log(Lophius_metazoan$prop_reads) ~ Lophius_metazoan$n)
summary(Lophius_metazoan_lm)
par(mfrow = c(2, 2))
plot(Lophius_metazoan_lm) #no

## witch (mifish significant)
witch_coi<-quantitative_landed%>%filter(species_name== 'Glyptocephalus cynoglossus', primer== 'COI')
witch_coi_lm<-lm(witch_coi$prop_reads ~ witch_coi$n)
summary(witch_coi_lm)


witch_mifish<-quantitative_landed%>%filter(species_name== 'Glyptocephalus cynoglossus', primer== 'MiFish')
witch_mifish_lm<-lm(witch_mifish$prop_reads ~ witch_mifish$n)
summary(witch_mifish_lm)
par(mfrow = c(2, 2))
plot(witch_mifish_lm)

witch_metazoan<-quantitative_landed%>%filter(species_name== 'Glyptocephalus cynoglossus', primer== 'Metazoan')
witch_metazoan_lm<-lm(witch_metazoan$prop_reads ~ witch_metazoan$n)
summary(witch_metazoan_lm)

# haddock COI, metazoan (significant)
haddock_coi<-quantitative_landed%>%filter(species_name== 'Melanogrammus aeglefinus', primer== 'COI')
haddock_coi_lm<-lm(haddock_coi$prop_reads ~ haddock_coi$n)
summary(haddock_coi_lm)

haddock_metazoan<-quantitative_landed%>%filter(species_name== 'Melanogrammus aeglefinus', primer== 'Metazoan')%>% drop_na(prop_reads)
haddock_metazoan_lm<-lm(haddock_metazoan$prop_reads ~ haddock_metazoan$n)
summary(haddock_metazoan_lm)


lobster_coi<-quantitative_landed%>%filter(species_name== 'Nephrops norvegicus', primer== 'COI')%>% drop_na(prop_reads)
lobster_coi_lm<-lm(log(lobster_coi$prop_reads) ~ lobster_coi$n)
summary(lobster_coi_lm)
hist(lobster_coi$prop_reads)
par(mfrow = c(2, 2))
plot(lobster_coi_lm)

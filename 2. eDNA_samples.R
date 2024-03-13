### eDNA samples analyses 

## libraries
library(tidyverse)
### here we (1) upload the output of the taxonomic assignment; (2) filter the species detected to avoid false positives
NS_coi_samples<-read.csv('../analysis_mock/data_final/level-6 (coi).csv')%>% pivot_longer(2:49, names_to = 'species_original', values_to = 'reads')%>% mutate(species= species_original)%>% mutate(amplicon= 'short', id= '94', primer='coi_NS')
NS_coi_samples<-NS_coi_samples%>% separate(species, c('Phylum',  'Class', 'Order', 'Family', 'Genus', 'genus', 'Species'))
NS_coi_samples[NS_coi_samples == '']<- NA 

NS_fish_samples<-read.csv('../analysis_mock/data_final/level-6 (fish).csv')%>% pivot_longer(2:59, names_to = 'species_original', values_to = 'reads')%>% mutate(species= species_original)%>% mutate(amplicon= 'short', id= '94', primer='fish_NS')
NS_fish_samples<-NS_fish_samples%>% separate(species, c('Phylum',  'Class', 'Order', 'Family', 'Genus', 'genus', 'Species'))
NS_fish_samples[NS_fish_samples == '']<- NA 

NS_metazoan_samples<-read.csv('../analysis_mock/data_final/level-6 (metazoan).csv')%>% pivot_longer(2:78, names_to = 'species_original', values_to = 'reads')%>% mutate(species= species_original)%>% mutate(amplicon= 'short', id= '94', primer='metazoan_NS')
NS_metazoan_samples<-NS_metazoan_samples%>% separate(species, c('Phylum',  'Class', 'Order', 'Family', 'Genus', 'genus', 'Species'))
NS_metazoan_samples[NS_metazoan_samples == '']<- NA 

NS_mifish_samples<-read.csv('../analysis_mock/data_final/level-6 (mifish).csv')%>% pivot_longer(2:60, names_to = 'species_original', values_to = 'reads')%>% mutate(species= species_original)%>% mutate(amplicon= 'short', id= '94')
NS_mifish_samples<-NS_mifish_samples%>% separate(species, c('Phylum',  'Class', 'Order', 'Family', 'Genus', 'genus', 'Species'))
NS_mifish_samples[NS_mifish_samples == '']<- NA 

## manually removing species not distinguishable using MiFish
mifish_unique_species<-NS_mifish_samples%>% mutate(species_res= ifelse(Species== 'gurnardus' | Species== 'cuculus', 'Triglidae',
                                                                       ifelse(Species== 'merlangus' | Species== 'aeglefinus', 'had-whiting',
                                                                              ifelse(Species== 'lucerna', 'Triglidae',
                                                                                     ifelse(Species==  'lanceolatus', 'Ammodytidae',
                                                                                            ifelse(Species== 'limanda' | Species== 'platessoides','dab-Aplaice',
                                                                                                   ifelse(Species== 'laterna', 'A.lant-megrim', 'no')))))))%>% filter(species_res!= 'no')%>% mutate(Genus= species_res, genus= species_res, Species= species_res, species_original= species_res)%>% group_by(sample, species_original)%>% mutate(reads_sum= sum(reads))%>% mutate(reads= reads_sum)%>% distinct()%>% select(1:14)

## add C. lucerna, L. wiff, H. lanceolatus
NS_mifish_samples_2<-NS_mifish_samples%>% mutate(species_res= ifelse(is.na(Species), 'no', 
                                                                     ifelse(Species== 'gurnardus' | Species== 'cuculus', 'Triglidae',
                                                                            ifelse(Species== 'lucerna', 'Triglidae',
                                                                                   ifelse(Species== 'merlangus' | Species== 'aeglefinus', 'had-whiting',
                                                                                          ifelse(Species==  'lanceolatus', 'Ammodytidae',
                                                                                                 ifelse(Species== 'limanda' | Species== 'platessoides','dab-Aplaice',
                                                                                                        ifelse(Species== 'laterna', 'A.lant-megrim', 'no'))))))))%>%filter(species_res== 'no')%>% select(!15)
NS_mifish_samples_3<-rbind(NS_mifish_samples_2, mifish_unique_species)

## uploading information on the fishing haul (target species)
target_species<-observed<-read.csv('samples_diversity_logbook_replicates.csv')%>% select(sample, Target_species)%>%distinct()%>% rename(target=Target_species)

#### calculating the relative abundance
reads_number<- rbind(NS_metazoan_samples, NS_mifish_samples_3, NS_coi_samples, NS_fish_samples)
reads_number_tot<- reads_number%>% group_by(sample, primer)%>% summarise(all_reads=sum(reads))
unassigned<- reads_number%>% filter(is.na(Species))%>% group_by(sample, primer)%>% summarise(reads_unass=sum(reads)) ## unassigned reads includes all reads unassigned and assigned at higher taxonomic level than species (e.g., genus, etc..)

samples_reads_un<-left_join(reads_number_tot, unassigned)%>%mutate(assigned_reads= all_reads - reads_unass)


all_samples<- reads_number%>% mutate (species_name = ifelse(species_original == 'Unassigned.__.__.__.__.__', 'Unassigned',
                                                            ifelse(is.na(Family), paste(Order),
                                                                   ifelse(is.na(Genus), paste(Family),
                                                                          ifelse(is.na(Species), paste(Genus), 
                                                                                 ifelse(Phylum== 'Unassigned', 'Unassigned', paste(Genus, Species)))))))%>% filter(species_name!= 'NA')

all_samples_short<-all_samples %>% mutate (level_taxon= ifelse(is.na(Species), 'higher', 'species'))%>%select(!1)
all_samples_short<- left_join(all_samples_short, target_species)%>%select(!3)

samples_abundance<- left_join(all_samples_short, samples_reads_un) %>%
  mutate(prop_reads_ass= reads/assigned_reads)%>% 
  mutate(prop_tot_reads= reads/all_reads)

## treshold for COI_NS (based on the thresholds calculated for the mock community)
coi_abundance<-samples_abundance%>% filter(primer== 'coi_NS')%>% filter(reads> 0)%>% filter(prop_reads_ass > 8.727374e-05)%>% filter(!is.na(Species))%>% filter(level_taxon== 'species')


### Metazoan remove species below abundance threshold, and for closely related species use the 1.5% reads abundance for filter

metazoan_abundance<- samples_abundance%>% filter(primer== 'metazoan_NS')%>% filter(reads> 0)%>% filter(prop_reads_ass > 1.681300e-03 )%>% filter(!is.na(Species))  ##previously 3.183902e-04
## 268 >0, 207 > prop_assigned

metazoan_abundance_fam<-metazoan_abundance%>% group_by(sample, Family)%>% mutate(family_reads=sum(reads))%>% mutate(family_reads_02= 0.015*family_reads)

pleuronecitdae<-metazoan_abundance_fam%>% filter(Family== 'Pleuronectidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)
pleuronecitdae_sp<-unique(pleuronecitdae$species_name)

lophius<-metazoan_abundance_fam%>% filter(Family== 'Lophiidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

triglidae<-metazoan_abundance_fam%>% filter(Family== 'Triglidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

pollachius<-metazoan_abundance_fam%>% filter(Genus== 'Pollachius')%>% group_by(sample)%>% mutate(family_reads=sum(reads))%>% mutate(family_reads_02= 0.015*family_reads)%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

genus_had_cod_wh<- c('Gadus', 'Merlangius', 'Melanogrammus')
had_cod_whiting<-metazoan_abundance_fam%>% filter (Genus%in% genus_had_cod_wh)%>% group_by(sample)%>% mutate(family_reads=sum(reads))%>% mutate(family_reads_02= 0.015*family_reads)%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)


metazoan_abundance_fam_no_flat<-metazoan_abundance_fam%>% filter(Family!= 'Pleuronectidae')%>% filter(Family!= 'Lophiidae')%>%filter(Family!= 'Triglidae')%>%filter(Genus!= 'Pollachius')%>%filter(Genus!= 'Gadus')%>%filter(Genus!= 'Merlangius')%>%filter(Genus!= 'Melanogrammus')
metazoan_abundance_final<- rbind(metazoan_abundance_fam_no_flat, pleuronecitdae, lophius, triglidae, pollachius, had_cod_whiting)%>% select(1:20)

## Mifish abundance and closely relates species filtering
mifish_abundance<- samples_abundance%>% filter(primer== 'mifish')%>% filter(reads> 0)%>% filter(prop_reads_ass > 1.710430e-04 )%>% filter(!is.na(Species)) ##1.712111e-04

mifish_abundance_fam<-mifish_abundance%>% group_by(sample, Family)%>% mutate(family_reads=sum(reads))%>% mutate(family_reads_02= 0.015*family_reads)

pleuronecitdae<-mifish_abundance_fam%>% filter(Family== 'Pleuronectidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)
pleuronecitdae_sp<-unique(pleuronecitdae$species_name)

lophius<-mifish_abundance_fam%>% filter(Family== 'Lophiidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

triglidae<-mifish_abundance_fam%>% filter(Family== 'Triglidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

clupeidae<-mifish_abundance_fam%>% filter(Family== 'Clupeidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

pollachius<-mifish_abundance_fam%>% filter(Genus== 'Pollachius')%>% group_by(sample)%>% mutate(family_reads=sum(reads))%>% mutate(family_reads_02= 0.015*family_reads)%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

genus_had_cod_wh<- c('Gadus', 'had-whiting')
had_cod_whiting<-mifish_abundance_fam%>% filter (Genus%in% genus_had_cod_wh)%>% group_by(sample)%>% mutate(family_reads=sum(reads))%>% mutate(family_reads_02= 0.015*family_reads)%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

mifish_abundance_fam_no_flat<-mifish_abundance_fam%>% filter(Family!= 'Pleuronectidae')%>% filter(Family!= 'Lophiidae')%>%filter(Family!= 'Clupeidae', Family!= 'Triglidae') %>%filter(Genus!= 'Pollachius')%>% filter (Genus!= 'Gadus', Genus!= 'had-whiting')
mifish_abundance_final<- rbind(mifish_abundance_fam_no_flat, pleuronecitdae, lophius, pollachius, had_cod_whiting, clupeidae, triglidae)%>% select(1:20)


## fish-2kb filtering

fish_abundance<- samples_abundance%>% filter(primer== 'fish_NS')%>% filter(reads> 0)%>% filter(prop_reads_ass > 1.334315e-03 
)%>% filter(!is.na(Species))## previous 4.509176e-04
## 268 >0, 207 > prop_assigned

fish_abundance_fam<-fish_abundance%>% group_by(sample, Family)%>% mutate(family_reads=sum(reads))%>% mutate(family_reads_02= 0.015*family_reads)

pleuronecitdae<-fish_abundance_fam%>% filter(Family== 'Pleuronectidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)
pleuronecitdae_sp<-unique(pleuronecitdae$species_name)

lophius<-fish_abundance_fam%>% filter(Family== 'Lophiidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

triglidae<-fish_abundance_fam%>% filter(Family== 'Triglidae')%>% mutate(above0.02 = ifelse(reads> family_reads_02, 'yes', 'no'))%>% filter(above0.02== 'yes')%>% select(1:22)

## no Sprat reads for fish marker (no need to filter for Clupeidae)
## Pollachius no reads for fish
## H. lanceolatus 1 read, that will be filtered out when inclduing spp identified with at least 3 reads

fish_abundance_fam_no_flat<-fish_abundance_fam%>% filter(Family!= 'Pleuronectidae')%>% filter(Family!= 'Lophiidae')%>%filter(Family!= 'Triglidae')
fish_abundance_final<- rbind(fish_abundance_fam_no_flat, pleuronecitdae, lophius, triglidae)%>% select(1:20)

## final datasets
final_species<-rbind(mifish_abundance_final, metazoan_abundance_final, fish_abundance_final, coi_abundance)

final_species<-final_species%>%filter(Species!= 'Triglidae', Species!='dab-Aplaice', Species!='had-whiting')
write.csv (final_species, 'final_species.csv')

## Taxonomic assignment of the mock reads analyses
### Uploading the libraries
library(tidyverse)
### importing the output of the taxonomic assignment into R (please note the output are csv files downloaded visualizing the table in QIIME)
## Each marker has a CSV file with the results of the taxonomic assignment of the reads
## For the four markers, we used the custom-made MIDORI databases, which include COI, 12S and 16S genes, whole mitogenome sequences filtered for North Sea fish species(including elasmobranchs and myxine)
## COI reads were also assigned using the MIDORI COI database without the ulterior filtering for NS species

# We import the results of the taxonomic assignment and tidy the datasets

NS_coi_mock<-read.csv('../analysis_mock/data_final/level-6 (mock_coi).csv')%>% pivot_longer(2:15, names_to = 'species_original', values_to = 'reads')%>% mutate(species= species_original)%>% mutate(amplicon= 'short', id= '94', primer='coi_NS')
NS_coi_mock<-NS_coi_mock%>% separate(species, c('Phylum',  'Class', 'Order', 'Family', 'Genus', 'genus', 'Species'))
NS_coi_mock[NS_coi_mock == '']<- NA 

midori_coi_mock<-read.csv('../analysis_mock/level-6 (mock_coi_midori).csv')%>% pivot_longer(2:28, names_to = 'species_original', values_to = 'reads')%>% mutate(species= species_original)%>% mutate(amplicon= 'short', id= '94', primer= 'midori_coi')
midori_coi_mock<-midori_coi_mock%>% separate(species, c('Phylum',  'Class', 'Order', 'Family', 'Genus', 'genus', 'Species'))
midori_coi_mock[midori_coi_mock == '']<- NA 

NS_mifish_mock<-read.csv('../analysis_mock/data_final/level-6 (mock_mifish).csv')%>% pivot_longer(2:25, names_to = 'species_original', values_to = 'reads')%>% mutate(species= species_original)%>% mutate(amplicon= 'short',id= '94', primer='mifish_NS')
NS_mifish_mock<-NS_mifish_mock%>% separate(species, c('Phylum',  'Class', 'Order', 'Family', 'Genus', 'genus', 'Species'))
NS_mifish_mock[NS_mifish_mock == '']<- NA 
#NS_mifish_mock<-NS_mifish_mock%>% mutate(mifish_species= paste(Genus, Species))%>% mutate(mifish_res= ifelse(mifish_species%in% species_no_res$species_cannot_distinguish, 'no species',
#ifelse(Genus %in% genus_no_res$genus_cannot_distinguish, 'no genus', 'yes'))) %>% filter(mifish_res== 'yes')%>% select(1:14)

NS_fish_mock<-read.csv('../analysis_mock/data_final/level-6 (mock_fish).csv')%>% pivot_longer(2:21, names_to = 'species_original', values_to = 'reads')%>% mutate(species= species_original)%>% mutate(amplicon= 'long', id= '94', primer='fish_NS')
NS_fish_mock<-NS_fish_mock%>% separate(species, c('Phylum',  'Class', 'Order', 'Family', 'Genus', 'genus', 'Species'))
NS_fish_mock[NS_fish_mock == '']<- NA 

NS_metazoan_mock<-read.csv('../analysis_mock/data_final/level-6 (mock_metazoan).csv')%>% pivot_longer(2:24, names_to = 'species_original', values_to = 'reads')%>% mutate(species= species_original)%>% mutate(amplicon= 'long', id= '94', primer='metazoan_NS')
NS_metazoan_mock<-NS_metazoan_mock%>% separate(species, c('Phylum',  'Class', 'Order', 'Family', 'Genus', 'genus', 'Species'))
NS_metazoan_mock[NS_metazoan_mock == '']<- NA 

## the datasets are merged
NS_mock_samples<-rbind(NS_coi_mock, NS_metazoan_mock,  NS_mifish_mock, NS_fish_mock, midori_coi_mock) 
NS_mock_samples<-NS_mock_samples%>% mutate (species_name = ifelse(species_original == 'Unassigned.__.__.__.__.__', 'Unassigned',
                                                                  ifelse(is.na(Family), paste(Order),
                                                                         ifelse(is.na(Genus), paste(Family),
                                                                                ifelse(is.na(Species), paste(Genus), 
                                                                                       ifelse(Phylum== 'Unassigned', 'Unassigned', paste(Genus, Species)))))))

NS_mock_samples<- NS_mock_samples%>% filter(species_name!= 'NA')%>% mutate (level_taxon= ifelse(is.na(Species), 'higher', 'species'))%>% select(2, 4, 14, 15)

## Filter for species detected with at least 3 reads
NS_mock_samples<-NS_mock_samples%>%filter(reads>2)
## Uploading a data set with the list of species added in the mock community
mock_species<- read.csv('../analysis_mock/mock_species.csv')%>% select(2:4)%>% mutate(level_taxon= 'species')
species_mock<- mock_species$species_name
## Merging the data sets with the species found in the mock by sequencing and the species added in the mock sample
NS_mock_final<-rbind(mock_species, NS_mock_samples)
NS_mock_final<-pivot_wider(NS_mock_final, names_from = primer, values_from = reads)

## Tidy the data set to make the plot (Figure S1)
NS_mock_final_lon<-pivot_longer(NS_mock_final, 4:8, names_to = 'primer', values_to = 'reads')%>% mutate(species_name_2=ifelse(species_name %in% species_mock, paste(species_name,'**'), species_name))
NS_mock_final_lon$reads<-as.numeric(NS_mock_final_lon$reads)
NS_mock_final_lon<-NS_mock_final_lon%>%mutate(database= ifelse(primer== 'midori_coi', 'Midori' , 'North Sea species'))

## Figure S1
data_plot1<-NS_mock_final_lon%>% filter(level_taxon== 'species')
data_plot1$database[data_plot1$database=='Midori']<- 'Midori'
data_plot1$database[data_plot1$database=='North Sea species']<- 'North Sea restricted database'
data_plot1$primer[data_plot1$primer=='midori_coi']<- 'COI'
data_plot1$primer[data_plot1$primer=='coi_NS']<- 'COI'
data_plot1$primer[data_plot1$primer=='fish_NS']<- 'Fish-2kb'
data_plot1$primer[data_plot1$primer=='metazoan_NS']<- 'Metazoan-2kb'
data_plot1$primer[data_plot1$primer=='mifish_NS']<- 'MiFish'
figureS1<-ggplot(data_plot1, aes(x = primer, y = species_name_2, fill=log(reads, 10)))+
  geom_tile(color = "black") +
  facet_grid(~database, space = 'free', scales = 'free')+
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlBu"),na.value = 'white') +
  geom_text(aes(label = reads), color = "black", size = 3)+
  theme_minimal()+
  labs(x=NULL, y= NULL)+
  theme(axis.text.y = element_text(face = 'italic'), legend.position = 'bottom')

ggsave(filename = "FigureS1.png", 
       plot = figureS1,
       width = 16, 
       height = 16,
       units = "cm",
       dpi = 400)

## Calculating the primer-specific thresholds for false positives based on the mock community data

#### sum of the reads assigned at the species level for each primer
number_reads<- NS_mock_final_lon%>% group_by(primer)%>% summarise(sum_reads= sum(reads, na.rm = TRUE))
number_unassigned_reads<- NS_mock_final_lon %>% filter(species_name== 'Unassigned')%>% group_by(primer) %>% mutate(unass_reads= sum(reads))%>% select(primer, unass_reads)%>% distinct()

## Merging the data sets with numbers of assigned reads and unassigned reads
number_reads_3<- left_join(number_reads, number_unassigned_reads)%>% mutate(assigned_reads= sum_reads-unass_reads)

## Looking at Figure S1 I can identify the maximum reads assigned as a false positive in the mock samples 
# NOTE, we exclude closely related species e.g. Clupea harengus and Sprattus sprattus; Lophius spp. 

reads_limit<- c(5,  29, 27, 23, 31)
## COI_NS (S. stellaris 5 reads used); Fish-2kb_NS (29 reads D. oxyrinchus); Metazoan-2kb_NS (27 reads D. oxyrinchus); Mifish_NS (31 reads, T.esmarkii) 

number_reads_4<- cbind(number_reads_3, reads_limit)
number_reads_4<- number_reads_4%>% mutate(abundance_sum= reads_limit/sum_reads) %>% mutate (abundance_assigned= reads_limit/assigned_reads)
write.csv(number_reads_4, 'primer-specificTresholds.csv')


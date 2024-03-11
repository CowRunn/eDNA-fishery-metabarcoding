## phyloseq analyses on consensus replicates
## here, we (1) convert our data into a phyloseq object, (2) calculated alpha diversity, (3) beta diversity using Jaccard distances, (4) PERMANOVA to see if target fishery and depth of the fishing haul are affecting the diversity observed in the eDNA
# uploading the libraries
library(tidyverse)
library(phyloseq)

## preparing the dataset
ph_final_data<- data_original_catch%>% filter(primer!= 'observed_landed')%>% droplevels()%>% select(!reads)%>%pivot_wider(values_from = prop_reads_ass, names_from = primer)%>% mutate(sample= paste('sample', sample, sep = '_'))
# deleted reads (column) and substituted total reads with prop_reads_ass
ph_final_data<-ph_final_data  %>%rowwise() %>% mutate(reads_sum= sum(COI, `fish-2kb`, MiFish, Metazoan, na.rm= TRUE))
unique(ph_final_data$species_name) ##unique species detected 52(excluding dab-amplaice had-wit)

replicates<- read.csv('../analysis_mock/shared_analyses/target.csv')%>% mutate(sample= paste('sample', sample, sep = '_'))
ph_final_data<-left_join(ph_final_data, replicates) 

ph_data<- ph_final_data %>% select(1:4, 9:11)

## preparing the data to convert it to a phyloseq object
ph_final_data_con<- data_original_catch%>% filter(primer!= 'observed_landed')%>% droplevels()%>% select(!prop_reads_ass)%>%pivot_wider(values_from = reads, names_from = primer)%>% mutate(sample= paste('sample', sample, sep = '_'))
ph_final_data_con<-ph_final_data_con  %>%rowwise() %>% mutate(reads_sum= sum(COI, `fish-2kb`, MiFish, Metazoan, na.rm= TRUE))

replicates<- read.csv('../analysis_mock/shared_analyses/target.csv')%>% mutate(sample= paste('sample', sample, sep = '_'))
ph_final_data_con<-left_join(ph_final_data_con, replicates)
ph_data_con<- ph_final_data_con %>% select(1:4, 9:11)

####### building the replicate consensus with species found across sample_replicate
unique_species_538<- ph_data %>% filter(sample_replicates == '538')%>% droplevels()%>% distinct(species_name)%>% mutate(sample= 'sample_538_c', level_taxon= 'species', reads_sum=1, sample_replicates= '538*', target= 'Roundfish', replicates= 'consensus*')
unique_species_544<- ph_data %>% filter(sample_replicates == '544')%>% droplevels()%>% distinct(species_name)%>% mutate(sample= 'sample_544_c', level_taxon= 'species', reads_sum=1, sample_replicates= '544*', target= 'Roundfish', replicates= 'consensus*')
unique_species_549<- ph_data %>% filter(sample_replicates == '549')%>% droplevels()%>% distinct(species_name)%>% mutate(sample= 'sample_549_c', level_taxon= 'species', reads_sum=1, sample_replicates= '549*', target= 'Nephrops', replicates= 'consensus*')
unique_species_569<- ph_data %>% filter(sample_replicates == '569')%>% droplevels()%>% distinct(species_name)%>% mutate(sample= 'sample_569_c', level_taxon= 'species', reads_sum=1, sample_replicates= '569*', target= 'Nephrops', replicates= 'consensus*')

ph_data_consensus_replicate<- rbind ((ph_data_con %>% filter(replicates== 'no')), unique_species_538, unique_species_544, unique_species_549, unique_species_569)

## creating the matrix for thhe phyloseq object
## creating matrix otu
#### matrix otu reads
final_data_primer<-data_original_catch%>%filter(primer!= 'observed_landed')%>% droplevels()
ph_otu_species<- final_data_primer%>% select(species_name)%>%droplevels()%>% distinct()%>% mutate (OTU= paste0("OTU", 1:52)) ##number of unique species detected across samples, primers

## use ph_data_consensus_replicate to replace replicates with a consensus replicate sample
ph_data_otu_con<- left_join(ph_data_consensus_replicate, ph_otu_species) 
ph_data_otu_con_w<- ph_data_otu_con %>% select(OTU, reads_sum, sample)%>% pivot_wider(names_from = sample, values_from = reads_sum)%>% column_to_rownames(var= 'OTU')

ph_data_otu_con_w[(ph_data_otu_con_w>0)]<-1
ph_data_otu_con_w[(is.na(ph_data_otu_con_w))]<-0

ph_data_otu_con_w_mat<-as.matrix(ph_data_otu_con_w) ##matrix OTU reads
head(ph_data_otu_con_w_mat)

## creating otu-taxa matrix
### taxmat
taxonomy<- read.csv('../analysis_mock/midori_NS_database/taxonomy.csv')%>% select(2:7)%>% mutate(Species= species_name)
ph_otu_species_taxonomy<- left_join(ph_otu_species, taxonomy)%>% column_to_rownames(var= 'OTU')

ph_otu_species_taxonomy_final<-as.matrix(ph_otu_species_taxonomy)
head(ph_otu_species_taxonomy_final)

## sample metadata
### sample
sample_data_df_con<- read.csv('../analysis_mock/shared_analyses/sample_logbook_depths2_consensus_repl.csv')%>%column_to_rownames(var= 'sample')%>% mutate(mean_depth= (depth_start+ depth_stop)/2)

sample_data_df_con$sample_replicates<- as.factor(sample_data_df_con$sample_replicates)

str(sample_data_df_con)

## phyloseq object
library(phyloseq)
OTU_con= otu_table(ph_data_otu_con_w_mat, taxa_are_rows = TRUE)
TAX_con= tax_table(ph_otu_species_taxonomy_final)
samples_con= sample_data(sample_data_df_con)

taxa_names(OTU_con)
taxa_names(TAX_con)
sample_names(OTU_con)
sample_names(samples_con)

phy_data_consensus<-phyloseq(OTU_con, TAX_con, samples_con)
plot_bar(phy_data_consensus, fill= 'Species')

### alpha diversity (plot)
plot_richness(phy_data_consensus , measures="Observed", color= ('sample_replicates'))+ scale_color_brewer(palette = 'Set3')+ geom_point(aes( shape= target), size= 4)
estimate_richness(phy_data_consensus, measures="Observed")

pcoa_data_J_con = ordinate(phy_data_consensus, "PCoA", "jaccard")

plot_ordination(phy_data_consensus, pcoa_data_J_con, color="mean_depth", shape="target") + 
  geom_point(size = 4)+
  geom_text(aes(label = sample_replicates), size = 4, vjust = 1.5)+
  #scale_color_brewer(palette = 'Set3')+
  theme_bw()+ ggtitle("PCoA using distance method Jaccard")

## PERMANOVA using the vegan package
library(vegan)
ph_data_df<-data.frame(otu_table(phy_data_consensus, taxa_are_rows = TRUE))
ph_data_df_t<-t(ph_data_df)
ph_data_df_t_matrix<- as.matrix(ph_data_df_t)
# calculating the Jaccard distances (Presence/absence)
distance_jac <- vegdist(ph_data_df_t_matrix, method = "jaccard")
ph_data_df_2<- as.data.frame(ph_data_df_t)%>% rownames_to_column('sample')
sample_data_df_rowcol<- sample_data_df_con%>% rownames_to_column('sample')
ph_data_df_2<- left_join(ph_data_df_2, sample_data_df_rowcol) ##sample_data_df with replicates
ph_data_df_2<- ph_data_df_2%>% mutate(depth_factor= ifelse(mean_depth< -100, 'deeper', 'shallower'))

## permanova using the target fisheries (Nephrops or Roundfish)
sample_div<- adonis2(distance_jac~target  , data=ph_data_df_2, permutations = 9999 )
sample_div

# using depth (coded as shallow and deeper based on 100 m threshold)
sample_div_depth<- adonis2(distance_jac~depth_factor  , data=ph_data_df_2, permutations = 9999 )
sample_div_depth

## PCoA plot of the species detected in the eDNA
cons_pcoa_haul<-plot_ordination(phy_data_consensus, pcoa_data_J_con, color="mean_depth", shape="target") + 
  geom_point(size = 4.2)+
  geom_text_repel(aes(label = Haul_l), size = 4, vjust = 1.5)+
  #scale_color_brewer(palette = 'Set3')+
  theme_bw()+ ggtitle("(B) Diversity detected in eDNA")+ theme(legend.position = "none")

## the PCoA plot of the species observed on-board is in beta_diversity_onboard.R
## figure 3 in beta_diversity_onboard.R

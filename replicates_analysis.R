## analyses of the sample replicates
### four four fishing hauls we analysed three sampling replicates of catch water (i.e., the three catch water samples collected were analysed)
## samples

## libraries
library(phyloseq)
library(vegan)

## building a phyloseq object from the original datasets with all the 18 sequenced samples
final_data_primer<-data_original_catch%>%filter(primer!= 'observed_landed')%>% droplevels()
ph_otu_species<- final_data_primer%>% select(species_name)%>%droplevels()%>% distinct()%>% mutate (OTU= paste0("OTU", 1:52)) ##number of unique species detected across samples, primers

ph_data_otu<- left_join(ph_data, ph_otu_species) 

ph_data_otu_w<- ph_data_otu %>% select(OTU, reads_sum, sample)%>% pivot_wider(names_from = sample, values_from = reads_sum)%>% column_to_rownames(var= 'OTU')

ph_data_otu_w[(ph_data_otu_w>0)]<-1
ph_data_otu_w[(is.na(ph_data_otu_w))]<-0

ph_data_otu_w_mat<-as.matrix(ph_data_otu_w) ##matrix OTU reads
head(ph_data_otu_w_mat)

## creating otu-taxa matrix
#### taxmat
taxonomy<- read.csv('../analysis_mock/midori_NS_database/taxonomy.csv')%>% select(2:7)%>% mutate(Species= species_name)
ph_otu_species_taxonomy<- left_join(ph_otu_species, taxonomy)%>% column_to_rownames(var= 'OTU')

ph_otu_species_taxonomy_final<-as.matrix(ph_otu_species_taxonomy)
head(ph_otu_species_taxonomy_final)

## sample metadata
#### sample
sample_data_df<- read.csv('../analysis_mock/shared_analyses/sample_logbook_depths2.csv')
hauls<- read.csv('../analysis_mock/plots_final/haul_modified_sample_logbook_depths2.csv')%>% select(2,18)
sample_data_df<-left_join(sample_data_df, hauls)

sample_data_df<-sample_data_df%>% column_to_rownames(var= 'sample')%>% mutate(mean_depth= (depth_start+ depth_stop)/2)
sample_data_df$sample_replicates<- as.factor(sample_data_df$sample_replicates)

str(sample_data_df)

## phyloseq object
library(phyloseq)
OTU= otu_table(ph_data_otu_w_mat, taxa_are_rows = TRUE)
TAX= tax_table(ph_otu_species_taxonomy_final)
samples= sample_data(sample_data_df)

taxa_names(OTU)
taxa_names(TAX)
sample_names(OTU)
sample_names(samples)

phy_data<-phyloseq(OTU, TAX, samples)

## creating a new dataset with only the samples for which replicates were analysed from the phyloseq object
phy_repl<-subset_samples(phy_data, replicates== 'yes')
plot_bar(phy_repl, fill= 'Species')

## species richness
plot_sp_repl<-plot_richness(phy_repl , x= 'Haul_n', measures="Observed", color=("sample_replicates"))+ scale_color_brewer(palette = 'Set2')+ geom_point(aes( shape= target), size= 3)+ theme_light()+
  theme(axis.text.x = element_text(angle = 90))+ ylab('Number of species detected')+ theme(legend.position = "none")
#plot_heatmap(phy_repl, na.value="white")

estimate_richness(phy_repl, measures="Observed")

## beta-diversity and PCoA
pcoa_replicates_J = ordinate(phy_repl, "PCoA", "jaccard") 

plot_PCoA_repl<-plot_ordination(phy_repl, pcoa_replicates_J, color="sample_replicates", shape="target") + #sample_replicates used before Haul
  geom_point(size = 3.2)+
  #scale_color_discrete(name = "Fishing haul", labels = c("haul 3", "haul 4", "haul 5", "haul 9"))+
  scale_color_brewer(palette = 'Set2', name = "Fishing haul", labels = c("haul 3", "haul 4", "haul 5", "haul 9"))+
  theme_bw()+ labs(color= 'Sample replicate')#+ ggtitle("PCoA using distance method Jaccard")

## Figure S6
library(cowplot)

FigureS6<-plot_grid(plot_sp_repl, plot_PCoA_repl, rel_widths = c(4, 6))
## added fishing hauls
FigureS6_haul<-plot_grid(plot_sp_repl, plot_PCoA_repl, rel_widths = c(4, 6))

ggsave(filename = "FigureS6_hauls.png", 
       plot = FigureS6_haul,
       width = 25, 
       height = 10,
       units = "cm",
       dpi = 900)

## PERMANOVA
mat_phy_repl<-as.matrix(t(data.frame(otu_table(phy_repl, taxa_are_rows = TRUE))))

distance_jac_repl <- vegdist(mat_phy_repl, method = "jaccard")
repl_phy_data<- as.data.frame(mat_phy_repl)%>% rownames_to_column('sample')
sample_data_df_rowcol_2<-sample_data_df%>% rownames_to_column('sample')
repl_phy_data<- left_join(repl_phy_data, sample_data_df_rowcol_2)##sample_data_df with replicates

# using sample_replicates as a factor
sample_div_repl<- adonis2(distance_jac_repl~sample_replicates  , data= repl_phy_data, permutations = 9999) ## lat_start significant ph_data_df_2
sample_div_repl
## how many more species detected?
species_replicates<-n_species_across_samples%>% mutate(sample= paste('sample', sample, sep = '_'))
species_replicates<-left_join(species_replicates, replicates)%>% filter(replicates== 'yes')

## figure S7 Venn's diagram
#libraries 
library(RColorBrewer)
library(eulerr)
myCol <- brewer.pal(3, "Pastel2")
###sample 538 Haul 3
sample_538_1<-ph_data %>% filter(sample == 'sample_538_1')
sample_538_1<-sample_538_1$species_name
sample_538_2<-ph_data %>% filter(sample == 'sample_538_2')
sample_538_2<-sample_538_2$species_name
sample_538_3<-ph_data %>% filter(sample == 'sample_538_3')
sample_538_3<-sample_538_3$species_name
observed_538<-obs_catch%>% filter(sample== '538_1')
observed_538<-observed_538$species_name
s538<-list(`538_1`= sample_538_1, `538_2`= sample_538_2, `538_3`= sample_538_3, 'observed onboard'= observed_538)
plot(euler(s538, shape = "ellipse"), quantities = TRUE, fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Sample 538')
plot_538<-plot(venn(s538), fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Sample 538')
## adding fishing haul
s538_haul<-list(`haul 3.1`= sample_538_1, `haul 3.2`= sample_538_2, `haul 3.3`= sample_538_3, 'observed onboard'= observed_538)
plot_538_haul<-plot(euler(s538_haul, shape = "ellipse"), quantities = TRUE, fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Haul 3')
plot_538_haul<-plot(venn(s538_haul), fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Haul 3')

## sample 544 Haul 4
sample_544_1<-ph_data %>% filter(sample == 'sample_544_1')
sample_544_1<-sample_544_1$species_name
sample_544_2<-ph_data %>% filter(sample == 'sample_544_2')
sample_544_2<-sample_544_2$species_name
sample_544_3<-ph_data %>% filter(sample == 'sample_544_3')
sample_544_3<-sample_544_3$species_name
observed_544<-obs_catch%>% filter(sample== '544_1')
observed_544<-observed_544$species_name
s544<-list(`544_1`= sample_544_1, `544_2`= sample_544_2, `544_3`= sample_544_3, 'observed onboard'= observed_544)
plot(euler(s544, shape = "ellipse"), quantities = TRUE, fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Sample 544')
plot_544<-plot(venn(s544), fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Sample 544')
## adding fishing haul
s544_haul<-list(`haul 4.1`= sample_544_1, `haul 4.2`= sample_544_2, `haul 4.3`= sample_544_3, 'observed onboard'= observed_544)
plot_544_haul<-plot(euler(s544_haul, shape = "ellipse"), quantities = TRUE, fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Haul 4')
plot_544_haul<-plot(venn(s544_haul), fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Haul 4')

## sample 549 Haul 5
sample_549_1<-ph_data %>% filter(sample == 'sample_549_1')
sample_549_1<-sample_549_1$species_name
sample_549_2<-ph_data %>% filter(sample == 'sample_549_2')
sample_549_2<-sample_549_2$species_name
sample_549_3<-ph_data %>% filter(sample == 'sample_549_3')
sample_549_3<-sample_549_3$species_name

observed_549<-obs_catch%>% filter(sample== '549_1')
observed_549<-observed_549$species_name

s4<-list(`549_1`= sample_549_1, `549_2`= sample_549_2, `549_3`= sample_549_3, 'observed onboard'= observed_549)
plot(euler(s4, shape = "ellipse"), quantities = TRUE, fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Sample 549')
plot_549<-plot(venn(s4), fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Sample 549')

## addign fihsing haul
s549_haul<-list(`haul 5.1`= sample_549_1, `haul 5.2`= sample_549_2, `haul 5.3`= sample_549_3, 'observed onboard'= observed_549)
plot_s549_haul<-plot(euler(s549_haul, shape = "ellipse"), quantities = TRUE, fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Haul 5')
plot_549_haul<-plot(venn(s549_haul), fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Haul 5')

## sample 569 Haul 9
sample_569_1<-ph_data %>% filter(sample == 'sample_569_1')
sample_569_1<-sample_569_1$species_name
sample_569_2<-ph_data %>% filter(sample == 'sample_569_2')
sample_569_2<-sample_569_2$species_name
sample_569_3<-ph_data %>% filter(sample == 'sample_569_3')
sample_569_3<-sample_569_3$species_name

observed_569<-obs_catch%>% filter(sample== '569_1')
observed_569<-observed_569$species_name

s569<-list(`569_1`= sample_569_1, `569_2`= sample_569_2, `569_3`= sample_569_3, 'observed onboard'= observed_569)
plot(euler(s569, shape = "ellipse"), quantities = TRUE, fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Sample 569')
plot_569<-plot(venn(s569), fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Sample 569')

## fishing haul
s569_haul<-list(`haul 9.1`= sample_569_1, `haul 9.2`= sample_569_2, `haul 9.3`= sample_569_3, 'observed onboard'= observed_569)
plot_s569_haul<-plot(euler(s569_haul, shape = "ellipse"), quantities = TRUE, fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Haul 9')
plot_569_haul<-plot(venn(s569_haul), fills = c("#B3E2CD", "#FDCDAC", "#CBD5E8"), main = 'Haul 9')


FigureS5_haul<-plot_grid(plot_538_haul,plot_544_haul, plot_549_haul, plot_569_haul)
## saved Venn diagram
ggsave(filename = "FigureS7.png", 
       plot = FigureS7,
       width = 24, 
       height = 24,
       units = "cm",
       dpi = 400)

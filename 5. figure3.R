library(phyloseq)
library(vegan)
library(cowplot)
## beta diversity of the species observed on-board 
observed_b<-observed%>% rename(species_name= sscinames)%>%
  mutate(reads= 1)%>% mutate(sample= paste('sample', sample, sep = '_'))

sample_unique<-c("sample_534_1", "sample_538_1",  "sample_544_1", "sample_549_1", "sample_559_1", "sample_564_1", "sample_569_1", "sample_575_1","sample_580_1", "sample_585_1")

## matrix of observed species

otu_observed<-observed_b%>% select(species_name)%>% distinct()%>%mutate (OTU= paste0("OTU", 1:35))

observed_b<- left_join(observed_b, otu_observed)
observed_w<-observed_b%>%filter(sample %in% sample_unique)%>% select(OTU, sample, reads)%>% pivot_wider(names_from = sample, values_from = reads)%>% column_to_rownames(var= 'OTU')
observed_w[is.na(observed_w)]<-0
observed_w<- as.matrix(observed_w)

#### otu of observed species
otu_observed_m<-otu_observed%>% column_to_rownames(var = 'OTU')
otu_observed_m<-as.matrix(otu_observed_m)

#### sample
sample_data_df_repl<- read.csv('../analysis_mock/shared_analyses/sample_logbook_depths2.csv')%>% filter(sample %in% sample_unique)%>% droplevels()%>% column_to_rownames(var= 'sample')%>% mutate(mean_depth= (depth_start+ depth_stop)/2)
sample_data_df_repl$sample_replicates<- as.factor(sample_data_df_repl$sample_replicates)

sample_data_df_repl<-sample_data_df_repl%>%mutate(haul_n= paste('haul.',Haul, sep='' ))

## phyloseq object
library(phyloseq)
OTU= otu_table(observed_w, taxa_are_rows = TRUE)
TAX= tax_table(otu_observed_m)
samples= sample_data(sample_data_df_repl)

## object
phy_data_obs<-phyloseq(OTU, TAX, samples)
plot_bar(phy_data_obs, fill= 'species_name')

pcoa_data_J = ordinate(phy_data_obs, "PCoA", "jaccard") ## PCoA script

obs_pcoa<-plot_ordination(phy_data_obs, pcoa_data_J, color="mean_depth", shape="target") + 
  geom_point(size = 4.2)+
  geom_text_repel(aes(label = haul_n), size = 4, vjust = 1.5)+
  #scale_colour_brewer(palette = 'Set3')+
  theme_bw()+ ggtitle("(A) Diversity observed on-board")+ theme(legend.position = "none")

## PERMANOVA
library(vegan)

ph_data_df<-data.frame(otu_table(phy_data_obs, taxa_are_rows = TRUE))
ph_data_df_t<-t(ph_data_df)
ph_data_df_t_matrix<- as.matrix(ph_data_df_t)
distance_jac <- vegdist(ph_data_df_t_matrix, method = "jaccard")
ph_data_df_2<- as.data.frame(ph_data_df_t)%>% rownames_to_column('sample')
sample_data_df_rowcol<- sample_data_df_repl%>% rownames_to_column('sample')
ph_data_df_2<- left_join(ph_data_df_2, sample_data_df_rowcol) ##sample_data_df with replicates

sample_div_target<- adonis2(distance_jac~target, data=ph_data_df_2) 
sample_div_target

sample_div_depth<- adonis2(distance_jac~mean_depth, data=ph_data_df_2) 
sample_div_depth

### plot the PCoA of species detected by eDNA and species observed on-board
obs_pcoa<-plot_ordination(phy_data_obs, pcoa_data_J, color="mean_depth", shape="target") + 
  geom_point(size = 4.2)+
  geom_text_repel(aes(label = haul_n), size = 4, vjust = 1.5)+
  #scale_colour_brewer(palette = 'Set3')+
  theme_bw()+ ggtitle("(A) Diversity observed on-board")+ theme(legend.position = "none")


## build a plot to get the legend
obs_pcoa_legend<-plot_ordination(phy_data_obs, pcoa_data_J, color="mean_depth", shape="target") + 
  geom_point(size = 4.2)+
  geom_text(aes(label = haul_n), size = 4, vjust = 1.5)+
  #scale_colour_brewer(palette = 'Set3')+
  theme_bw()+ ggtitle("Species observed onboard")+ labs(col= 'Mean depth', shape= 'Target fishery')

library(cowplot)
legend <- get_legend(obs_pcoa_legend)
## cons_pcoa_haul is the PCoA plot of the species detected by the eDNA, it is in phy_data_consensus
pcoa_haul<-plot_grid(obs_pcoa, cons_pcoa_haul)
Figure3<-plot_grid(pcoa_haul, legend, rel_widths = c(3, .4))
ggsave(filename = "Figure3.png", 
       plot = Figure3,
       width = 26, 
       height = 12,
       units = "cm",
       dpi = 400)

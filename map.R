## libraries
library(tidyverse)
library(marmap)
library(mapdata)
library(ggrepel)
library("gridExtra")

## downloading the bathymetry of the Skagerrak and North Sea
# defining the latitude and longitude limits
ylim = c(30., 80.)
xlim = c(-30, 45)

xlim <- c(-28, -10)
ylim <- c(62.5, 67.5)
## downloading from NOAA
depth <- 
  getNOAA.bathy(lon1 = xlim[1], lon2 = xlim[2],
                lat1 = ylim[1], lat2 = ylim[2],
                resolution = 1) %>% 
  fortify() %>%  # turn the object into a data.frame
  filter(z <= 0)
glimpse(depth)

depth_Skagerrak<-depth%>% filter(x>4 & x<15)%>% filter(y>55 & y< 60) 
summary(depth_Skagerrak)

## building the map
w <- map_data("worldHires", ylim = c(30, 80), xlim = c(-30, 45))
depth_Skagerrak<- read.csv('../diversity_samples/depth_Skagerrak.csv')
## uploading hauls coordinates
samples_logb_edited<- read.csv('../analysis_mock/shared_analyses/sample_logbook_depths2_consensus_repl.csv')%>% rename(Target=target)
## defining the colors
library(RColorBrewer)
blues_colors_4<-c("#08306B", "#083D7F", "#084B94", "#0E59A2", "#1966AD", "#2474B6",  "#3282BD", "#4090C5", "#519CCB", "#62A8D2", "#75B3D8" ,"#8BBFDC", "#C1D9ED", "#DEEBF7", "#F7FBFF")

map_depth_3<-ggplot() +
  geom_raster(data = depth_Skagerrak,
              aes(x = x, y = y, fill = z), show.legend = FALSE)+# run this as TRUE when you want to get the legend
  geom_polygon(data = w, aes(long, lat, group = group), fill = "darkgrey", color = NA)+ 
  coord_quickmap(xlim = c(8, 12), ylim = c(57, 58.5))+
  geom_contour(data = depth_Skagerrak, 
               aes(x=x, y=y, z=z),
               breaks=c(-100),
               size=c(0.3),
               linetype='dashed',
               colour="white")+
  geom_contour(data = depth_Skagerrak, 
               aes(x=x, y=y, z=z),
               breaks=c(-250),
               size=c(0.6),
               colour="white")+
  geom_contour(data = depth_Skagerrak, 
               aes(x=x, y=y, z=z),
               breaks=c(-50),
               size=c(0.6),
               linetype='dotted',
               colour="blue")+
  scale_fill_gradientn(name= 'Depth', colors= blues_colors_4)+ ##labels=c(0,-100, -200, -300, -400, -500, -600
  geom_point(data= samples_logb_edited,
             aes(x= Longitude_start, y= Latitude_start, group=NULL, fill=NULL, shape= Target, color= Target), size=4, show.legend = FALSE)+# set to TRUE when you have to save the legend
  scale_color_manual(values = c('darkred','orange'))+
  geom_text_repel(data= samples_logb_edited, aes(x= Longitude_start, y= Latitude_start, label = Haul_label), size=4.8)+ labs(x= 'Longitude', y= 'Latitude')+ #sample_label
  theme(panel.background = element_rect(fill = "#F7FBFF",#"#C1D9ED",
                                        colour = "#F7FBFF"), #"#C1D9ED"
        legend.text = element_text(size = 12),
        legend.title = element_text(size=12),
        axis.title = element_text(size=8),
        axis.text = element_text(size= 8),
        plot.margin = unit(c(0.4, 0, 0.2, 0), 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),#)+
        legend.position="bottom", legend.box = 'vertical')+ #'legend.justification = "left"'
  annotate(geom= 'text', y=57.4, x=8.2, label="50m", size= 4.3)+
  annotate(geom= 'text', y=57.5, x=8.1, label="100m", size= 4.3)+
  annotate(geom= 'text', y=57.73, x=8.1, label="250m", size= 4.3)

## map of Europe
Europe<-ggplot() +
  geom_polygon(data = w, aes(long, lat, group = group), fill='darkgrey')+ 
  coord_quickmap(xlim = c(-10, 26), ylim = c(39, 66))+
  geom_rect(aes(xmin = 8, xmax = 12, ymin = 57, ymax = 59), 
            fill = "blue", alpha = 0.2, color = "black")+ 
  theme_bw()+theme(panel.background = element_rect(fill = 'white', #"#8BBFDC",
                                                   colour = 'white'), #"#8BBFDC"), ##C1D9ED
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.title=element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank(), 
                   plot.margin=grid::unit(c(0,0,0,0), "mm"))+labs(x=NULL, y=NULL)

## lifestyle species plot
lifestyle_plot<-ggplot(species_lifestyle, aes(x=Haul_n,  fill= life)) +
  geom_bar (stat = 'count', width = 0.9, colour='white', show.legend = FALSE)+ ## run this as TRUE when you want to get the legend
  scale_fill_brewer(palette = 'Set2')+
  facet_grid((detection~target), scales = 'free_x')+ theme_minimal()+ 
  theme(legend.position="bottom",
        legend.background = element_rect(color= 'lightblue', linewidth =0.5, linetype="solid"),
        legend.text = element_text(size = 8),
        strip.background = element_rect(colour = 'lightblue',linewidth =0.2, fill = 'white'),
        strip.text.x = element_text(size=10, color = 'black', face = 'italic'),
        strip.text.y = element_text(size=10, color = 'black'),
        axis.text.x = element_text(angle = 90))+
  ggtitle('')+labs(y='Number of species', x=NULL, fill= 'Category')

## function to get the legend out of the plots
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
## to get the legend, run the ggplot code for each plot with the option in the second line show.legend=TRUE, then run these to lines to get the legends of the plots
## after getting the legend, set the show.legend= FALSE and run again the ggplot codes to make the plots without legend this time.
lifestyle_legend<- get_legend(lifestyle_plot)
legend_map_center <- get_legend(map_depth_3)

## 

figure_1<-grid.arrange(map_depth_3, lifestyle_plot, legend_map_center, Europe, lifestyle_legend,
                       widths=c(2.6, 1.4 ,1.6, 1.6), heights=c(4, 1.4, 1), layout_matrix= rbind(c(1,1, 2,2), c(1, 1, 5, 5), c(3, 4, 5,5)))

ggsave(filename = "Figure1_grid_new.png", 
       plot = figure_1,
       width = 26, 
       height = 16,
       units = "cm",
       dpi = 900)

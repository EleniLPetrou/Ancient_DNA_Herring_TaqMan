# The purpose of this script is to create a map using complex Lat/Long data

# Load libraries
library(tidyverse)
library(rgdal)
library(ggmap)
library(viridis)
library(ggrepel)
library(cowplot)

# Specify the directories containing data files
DATADIR <- ("./Ancient DNA and archaeology/mapping_data")

# Specify the names of the data files needed for the analysis:
# Text file containing information about spawning beaches
FILE1 <- "WDFW_herring_spawning_location_shapefile_df.txt"

# Text file containing information about archaeological sites
FILE2 <- "mapping_arch_sites.txt"

# text file containing information about modern sampling sites
FILE3 <-  "mapping_modern_pops.txt"

# output file
OUTFILE <- "Figure1_map_manuscript.pdf"

#############################################################
# Read the files needed for the analysis
setwd(DATADIR)

plotting_df <- read.delim(FILE1)
arch_sites <- read.delim(FILE2)
modern_pops <- read.delim(FILE3)

# Check that the spawning areas look right by plotting

myLocation <- c(left = -124, bottom = 47, right = -121.5, top = 49.1)

myMap <- get_stamenmap(bbox = myLocation,
                       source = "stamen", 
                       maptype = "terrain-background")

# Plot the map

plot1 <- ggmap(myMap) +
  
  # plot spawning beaches
  geom_path(data = plotting_df, 
            aes(x = long, y = lat, group = group, color = avg_DOY), size = 0.75) +
  geom_polygon(data = plotting_df, aes(x = long, y = lat, group = group, fill = avg_DOY)) +
 
  # plot archaeological sites
  geom_point(data = arch_sites, 
             aes(x = longitude, y = latitude), 
             size = 2, alpha = 1, shape = 16, color = "red") +
  
  # Specify some other details
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude")

plot1

## Annotate the plot
mybreaks = c(0, 30, 60, 90, 120)
mylabels = c("January", "February", "March", "April", "May")

plot2 <- plot1 +
  
  #Axis labels
  labs(x = "Longitude", y = "Latitude") +
  
  # Color scheme for points
  scale_fill_viridis(option = "plasma",
                     name = "Spawning date",
                     breaks = mybreaks, 
                     labels = mylabels,
                     begin = 0, end = 0.9) +
  
  # Color scheme for paths
  scale_color_viridis(option = "plasma",
                      name = "Sampling date",
                      breaks = mybreaks, 
                      labels = mylabels,
                      begin = 0, end = 0.9,
                      guide = FALSE) +
  
  theme(legend.position = "right") +
  
  # annotate archaeological sites
  geom_label_repel(data = arch_sites, aes(x = longitude, y = latitude, label = site), 
                  size = 3,
                  min.segment.length = 0,
                  segment.size = 0.6,
                  point.padding = 0.2,
                  nudge_x = 0.80, 
                  segment.color = 'black') 

plot2

# save the legend of plot 2 as an object

mylegend <- get_legend(plot2)

## Insert map of North America

myLocation2 <- c(left = -135, bottom = 22, right = -65, top = 60)

myMap2 <- get_stamenmap(bbox = myLocation2,
                       source = "stamen", maptype = "terrain-background", zoom = 5)

# Make a small map that annotates the study area with rectangle 

insert_map <- ggmap(myMap2) +
    annotate("rect", xmin = -121, xmax = -125, 
           ymin = 47, ymax = 50, color = "black", alpha = 0, size = 0.7) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

insert_map

# combine the plots

final_plot <- ggdraw() +
  draw_plot(plot2 + theme(legend.position = 'none'),  x = -0.1, y = 0, width = 0.9, height = 0.9) +
  draw_plot(insert_map, x = 0.70, y = 0.61, width = 0.3, height = 0.3) +
  draw_plot(mylegend, x = 0.32, y = -0.12) +
draw_plot_label(label = c("A", "B"), size = 12, x = c(0.1, 0.70), y = c(0.96, 0.96))

final_plot


# Save plot to file
ggsave(OUTFILE, final_plot)

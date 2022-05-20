# The purpose of this script is to format shapefile data from WDFW 
# into a format that is appropriate for plotting maps in R

# Load libraries
library(tidyverse)
library(rgdal)
library(ggmap)


# Specify the directories containing data files
DATADIR <- ("./Ancient DNA and archaeology/mapping_data")

# Specify the names of the data files needed for the analysis:

# Tab-delimited text file containing info about spawning activity at each site. Made by WDFW.
SPAWN_FILE <- "WDFW_data_herring_WOY.txt" 

# Folder containing shapefile made by WDFW showing all of the herring spawning beaches in WA
WDFW_FOLDER <- "WDFW_Herring_Spawn_Areas"
WDFW_LAYER <- "Herring" #name of the shapefile (without the .shp extension)

# Specify name of output file
OUTFILE <- "WDFW_herring_spawning_location_shapefile_df.txt"

#############################################################
# Read the files needed for the analysis
setwd(DATADIR)

# .shp file containing outlines of the herring spawning beaches. 
shp <- readOGR(dsn = WDFW_FOLDER, layer = WDFW_LAYER, stringsAsFactors = F)

# text file containing info about the Julian date of spawning
df <- read.delim(SPAWN_FILE)

############################################################
# Process the data in preparation of plotting

# Make some dataframes from the shapefile components
shapefile.dt <- data.table(shp@data)

# Use the broom package to turn the shapefile into a dataframe
shapefile_df <- broom::tidy(shp, region = "Name")
lapply(shapefile_df, class)
head(shapefile_df)

# Calculate the mean spawning date (DOY) for each spawning location in the df

avg_spwn_time_df <-  df %>% 
  filter(PeakWeek == 1) %>% 
  group_by(Area) %>% 
  summarise(avg_DOY = mean(WOY)*7)

head(avg_spwn_time_df)

plotting_df <- left_join(shapefile_df, 
                         avg_spwn_time_df, 
                         by = c("id" = "Area")) %>%
  filter(avg_DOY != "NA") #filter missing data, if needed

head(plotting_df)


# Check that the spawning areas look right by plotting

myLocation <- c(left = -125, bottom = 47, right = -121.5, top = 49.1)

myMap <- get_stamenmap(bbox = myLocation,
          source = "stamen", 
          maptype = "terrain-background")


ggmap(myMap) +
  geom_path(data = plotting_df, 
            aes(x = long, y = lat, group = group, color = avg_DOY), 
            size = 1, alpha = 0.8) +
  theme_classic()

# Save the spawn data as a tab-delimited text file for more detailed plotting

write.table(plotting_df, OUTFILE, quote = FALSE, sep = "\t", row.names = FALSE)

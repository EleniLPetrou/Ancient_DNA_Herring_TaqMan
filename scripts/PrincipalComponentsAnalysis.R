# This script uses genotype data (in genepop format) to conduct a Principal Components Analysis with 
# adegenet. The PCA is color-coded by reading in sample metadata from another dataframe. 

# Load the necessary libraries
library(adegenet)
library(hierfstat)
library(tidyverse)
library(cowplot)
library(viridis)

# Specify directories and file names
DATADIR <- "./analyses_PCA"
DATAFILE <- "genepop_ancient_modern_7loci_20190515.gen"
METADATAFILE <- "sample_metadata_20211206.txt"
OUTFILE <- "plot_PCA_manuscript_20211206.pdf"

# setwd
setwd(DATADIR)

####################################################################################
# Read in the files needed for the analysis. First, read in your data as a genepop file.
# The file can be delimited by tabs or spaces but there must abe a comma after each individual. 
# Specify how many characters code each allele with ncode argument. 

Adata <- read.genepop(DATAFILE, ncode = 3)
(summary(Adata))

# read in the sample metadata for plotting
metadata <- read.delim(METADATAFILE, sep = "\t", header = TRUE)

# To replace missing data information with the mean
Adata_scaled <- scaleGen(Adata, NA.method = "mean")

# To conduct the PCA
pca_A <- dudi.pca(Adata_scaled, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
summary(pca_A) 
# The Projected inertia columns store information on the percent variance explained by each axis

A_df <- pca_A$li

# Combine the PCA results df with the sample metadata to make some nice plots
PCA_df <- dplyr::bind_cols(A_df, metadata)

####################################################################################
# Plot the data
# First, partition the PCA results, so you can plot them as individual layers
head(PCA_df)

modern_samples <- dplyr::filter(PCA_df, era == "modern")
Burton_samples <- dplyr::filter(PCA_df, site == "Burton Acres")
Bay_samples <- dplyr::filter(PCA_df, site == "Bay Street")
ancient_samples <- dplyr::filter(PCA_df, era == "ancient")

# Set the breaks for your color scale (continuous)
mybreaks <-  c(0, 30, 60, 90, 120, 150)

# Specify labels for color scale
mylabels = c("January", "February", "March", "April", "May", "June")

# Make plots
plot_burton <- ggplot() +
  geom_point(data = modern_samples, aes(x = Axis1, y = Axis2, color = julian_date), 
             size = 1.5, alpha = 0.6, shape = 16) +
  geom_point(data = Burton_samples, aes(x = Axis1, y = Axis2, shape = layer), 
             size = 3, alpha = 0.8, color = "grey27") +
  scale_color_viridis(option = "plasma", name = "Sampling date",
                      breaks = mybreaks, labels = mylabels, begin = 0, end = 0.9) +
  scale_shape_manual(name = "Age", values = c(2, 0)) + 
  theme_classic() + 
  theme(axis.ticks = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.title = element_text(size = 11)) +
  ylim(-6,3) +
  xlab("35% of variation") +
  ylab("15% of variation")

plot_burton


plot_bay <- ggplot() +
  geom_point(data = modern_samples, aes(x = Axis1, y = Axis2, color = julian_date), 
             size = 1.5, alpha = 0.6, shape = 16) +
  geom_point(data = Bay_samples, aes(x = Axis1, y = Axis2, shape = layer), 
             size = 3, alpha = 0.8, color = "grey27") +
  scale_color_viridis(option = "plasma", name = "sampling date (doy)",
                      breaks = mybreaks, begin = 0, end = 0.9, guide = 'none') +
  scale_shape_manual(name = "Age", values = c(2, 0)) + 
  theme_classic() + 
  theme(axis.ticks = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.title = element_text(size = 11)) +
  ylim(-6,3) +
  xlab("35% of variation") +
  ylab("15% of variation")

plot_bay

# Save legend from first plot
burton_legend <- get_legend(plot_burton)

# Combine plots into one pdf

final_plot <- ggdraw() +
  draw_plot(plot_burton + theme(legend.position = "none"), x = 0, y = 0.5, width = 0.55, height = 0.45) +
  draw_plot(plot_bay, x = 0.0, y = 0.00, width = .8, height = .45) +
  draw_plot(burton_legend, x = 0.30, y = 0.40, width = .8, height = .45) +
  draw_plot_label(label = c("A: Burton Acres", "B: Bay Street"), size = 12,
                  x = c(0, 0), y = c(1, 0.5))
final_plot

# Save final_plot as pdf
ggsave(OUTFILE, final_plot)

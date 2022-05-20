# The purpose of this script is to conduct simulations of individual assignment using the rubias package

library(rubias)
library(tidyverse)
library(ggdist)

# Turn off scientific notation
options(scipen=999)

# Specify the input file name(containing reference samples in rubias format)
ref_file <- "input_ref_pops_3reportinggroups.txt"
################################################################################

# Read in the data (genotypes in rubias format)
ref_df <- read.delim(ref_file, header = TRUE, sep = "\t",  
                     colClasses = c(rep("character", 4), rep("integer", 14)), na.strings = "NA")

# Do leave-one-out self-assignment of individuals
individual_simulations <- self_assign(
  reference = ref_df, 
  gen_start_col = 5)

# Notes from: https://github.com/eriqande/rubias#self-assigning-fish-from-the-reference
# The output of this simulation is a data frame containing likelihood scores. 
# The log_likelihood is the log probability of the fish's genotype given it is 
# from the inferred_collection computed using leave-one-out. 
# The scaled_likelihood is the posterior prob of assigning the fish to the 
# inferred_collection given an equal prior on every collection in the reference. 
# Other columns are as in the output for infer_mixture(). 
# Note that the z_score computed here can be used to assess the distribution of 
# the z_score statistic for fish from known, reference populations. 
# This can be used to compare to values obtained in mixed fisheries

# The output can be summarized by repunit as was done above:
  
sa_to_repu <- individual_simulations %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))


# Take a peek at the highest likelihoods

highest_likelihoods <- individual_simulations %>%
  group_by(indiv) %>%
  filter(scaled_likelihood == max(scaled_likelihood))

# Now let's summarize the data in a different way; let us sum the likelihoods within specific reporting groups.

group_likelihoods <- individual_simulations %>%
  group_by(indiv, inferred_repunit) %>%
  summarise(summed_likelihood = round(sum(scaled_likelihood), 3)) %>%
  filter(summed_likelihood > 0.80)

# Let's add some metadata to these results

meta_df <- ref_df %>%
  select(indiv, repunit)

group_likelihoods <- left_join(meta_df, group_likelihoods, by = "indiv")

# Plot the results
my_cols <- c("#2b83ba","#abdda4",  "#fdae61")

ggplot(group_likelihoods, aes(x = repunit, y = summed_likelihood))+
  stat_dots(aes(color = inferred_repunit), size = 3, alpha = 0.8, side = "both")+
  ylim(0.7,1) +
  ylab("Scaled likelihood") +
  xlab("True reporting group") +
  scale_color_manual(name = "Inferred reporting group",labels = c("Jan-Feb", "Mar-Apr", "May-Jun"), values = my_cols)+
  scale_x_discrete(labels = c("Jan-Feb", "Mar-Apr","May" ))+
  theme_bw()


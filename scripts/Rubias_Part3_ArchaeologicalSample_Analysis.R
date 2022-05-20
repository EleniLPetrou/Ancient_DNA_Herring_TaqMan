# The purpose of this script is to conduct a mixed stock analysis using rubias.
# Modern herring samples are used as the genetic baseline while ancient 
# samples are analyzed as the mixtures.

# Load libraries
library(tidyverse)
library(rubias)

# Specify data directories and file names
#DATADIR <- "./analyses_rubias"
ref_file <- "input_ref_pops_3reportinggroups.txt" # modern herring genotypes
mix_file <- "input_ancient_mix_pops.txt" # ancient herring genotypes

# Specify output file names
output_file1 <- "results_rubias_analysis_ancient_3repgroups.txt" #results
output_file2 <- "results_rubias_posterior_probabilties_ancient_3repgroups.txt"
output_file3 <- "plot_results_rubias_ancient_3repgroups.pdf"

# Read in data
# The data frames have 4 columns of metadata (character format) 
# followed by 14 columns of genotypes (integer format)

ref_df <- read.delim(ref_file, header = TRUE, sep = "\t",  colClasses = c(rep("character", 4), rep("integer", 14)), na.strings = "NA")
mix_df <- read.delim(mix_file, header = TRUE, sep = "\t", colClasses = c(rep("character", 4), rep("integer", 14)),  na.strings = "NA")

# metadata columns:
#1. sample_type: a column telling whether the sample is a reference sample or a mixture sample
#2. repunit: the reporting unit that an individual/collection belongs to. This is required if sample_type is reference. And if sample_type is mixture then repunit must be NA. This must be a character vector.
#3. collection: for reference samples, the name of the population that the individual is from. For mixture samples, this is the name of the particular sample (i.e. stratum or port that is to be treated together in space and time). This must be a character, not a factor.
#4. indiv: a character vector with the ID of the fish. These must be unique.

##############################

#Performing a Genetic Mixture Analysis. This is done with the infer_mixture function. 

mix_est <- infer_mixture(reference = ref_df, 
                         mixture = mix_df, 
                         gen_start_col = 5,
                         method = "MCMC",
                         reps = 10000, burn_in = 1000)

#The result comes back as a list of four tidy data frames.
lapply(mix_est, head)

# Explanation of terms:
# 1. mixing_proportions: the mixing proportions. The column pi holds the estimated mixing proportion for each collection.
# 2. indiv_posteriors: this holds, for each individual, the posterior means of group membership in each collection. Column PofZ holds those values. Column log_likelihood holds the log of the probability of the individuals genotype given it is from the collection. Also included are n_non_miss_loci and n_miss_loci which are the number of observed loci and the number of missing loci at the individual. A list column missing_loci contains vectors with the indices (and the names) of the loci that are missing in that individual. It also includes a column z_score which can be used to diagnose fish that don't belong to any samples in the reference data base (see below).

# Create a dataframe summarizing the most likely reporting unit of origin for each individual by summing across the posterior probabilities across collections in a reporting group. 

options(scipen = 999) # turn off scientific notation

post_probs_df <- mix_est$indiv_posteriors %>%
  select(-missing_loci)  %>%
  group_by(indiv, repunit) %>%
  mutate(sum_post = sum(PofZ)) %>%
  ungroup() %>%
  group_by(indiv) %>%
  top_n(1, PofZ) # this is just to get the top reporting group for each sample

######################################################################
# Step 2: Assess whether individuals are not from any of the reference populations

#In this case, it is useful to look at the raw log-likelihood values computed for the individual, rather than the scaled posterior probabilities. Because aberrantly low values of the genotype log-likelihood can indicate that there is something wrong. However, the raw likelihood that you get will depend on the number of missing loci, etc. rubias deals with this by computing a z-score for each fish. The Z-score is the Z statistic obtained from the fish's log-likelihood (by subtracting from it the expected log-likelihood and dividing by the expected standard deviation). rubias's implementation of the z-score accounts for the pattern of missing data, but it does this without all the simulation that gsi_sim does.This makes it much, much, faster-fast enough that we can compute it by default for every fish and every population.
#Here, we will look at the z-score computed for each fish to the population with the highest posterior. (It is worth noting that you would never want to use the z-score to assign fish to different populations-it is only there to decide whether it looks like it might not have actually come from the population that it was assigned to, or any other population in the reference data set.)

# Get the maximum a-posteriori population for each individual

map_rows <- mix_est$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()

head(map_rows)

# If everything is good, then we expect that the z-scores we see will be roughly normally distributed. 
# We can compare the distribution of z-scores we see with a bunch of simulated normal random variables.

normo <- tibble(z_score = rnorm(1e06))

ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")

#######################################################
# Compute Credible Intervals for the mixture proportion estimates using the $mix_pop_traces df

# First, check how many MCMC sweeps were done:
nsweeps <- max(mix_est$mix_prop_traces$sweep)
nsweeps

# Discard the first 200 sweeps as burn-in
test <- mix_est$mix_prop_traces
head(test)

trace_subset <- mix_est$mix_prop_traces %>%
  filter(sweep > 200) %>%
  group_by(mixture_collection, repunit, sweep) %>%
  summarise(repprop = sum(pi))

head(trace_subset)

# Now we can plot those traces:
ggplot(trace_subset, aes(x = repprop, colour = repunit)) +
  geom_density() +
  facet_wrap(~mixture_collection) +
  theme_bw()

# Next, let's take that information and calculate the 95% credible intervals from it
credible_interval_df <- trace_subset %>%
  group_by(mixture_collection, repunit) %>%
  summarise(loCI = round(quantile(repprop, probs = 0.05), 2),
            hiCI = round(quantile(repprop, probs = 0.95),2))

head(credible_interval_df)

# Next, add the credible interval df to the mean_pi df, to facilitate plotting

plotting_df <- left_join(rep_mix_est, credible_interval_df, by = c("mixture_collection", "repunit"))
head(plotting_df)

# Round the results 2 figures 
plotting_df$repprop_round <- round(plotting_df$repprop, 2)

plotting_df %>%
  group_by(mixture_collection) %>%
  summarise(tot_mix = sum(repprop_round))

##########################################################################
# Plot the results
 
# Pick colors from Rcolorbewer "spectral"
my_cols <- c("#2b83ba","#abdda4",  "#fdae61")

# specify some custom facet labels
my_facet_labels <- c(bay_new = "Bay Street: 400-100 cal BP", 
                     bay_old = "Bay Street: 800-550 cal BP", 
                     burton_new = "Burton Acres: post-contact", 
                     burton_old = "Burton Acres: 910-685 cal BP")

my_rg_labels <- c("Jan-Feb", "Mar-Apr", "May")

# Simple bar plot
bar_plot <- ggplot() +
  geom_bar(data = plotting_df, 
           aes(x = "", y = repprop_round, fill = repunit),
           width = 0.4, 
           stat = "identity", 
           alpha = 0.8) +
  ylab("Estimated mixing proportion") +
  xlab("") +
  facet_wrap(~mixture_collection, 
             labeller = labeller(mixture_collection = my_facet_labels)) +
  scale_fill_manual(name = "Reporting group",labels = my_rg_labels, values = my_cols) +
  theme_classic()

bar_plot

# Bar plot with credible intervals
my_size = 0.5
my_width = 0.2


CI_plot <- bar_plot + 
  geom_errorbar(data = filter(plotting_df, repunit == "May_Jun"),  
                aes(x = "", ymin = loCI, ymax = hiCI), 
                colour = "black", size = my_size, width = my_width) +
  geom_errorbar(data = filter(plotting_df, repunit == "Mar_Apr", mixture_collection != "burton_old"),  aes(x = "", ymin = loCI, ymax = hiCI), 
                colour = "black",  size = my_size , width = my_width) +
  geom_errorbar(data = filter(plotting_df, repunit == "Jan_Feb", mixture_collection == "burton_old"),  aes(x = "", ymin = 1 - hiCI , ymax = 1 - loCI ), 
                colour = "black", size = my_size, width = my_width)

CI_plot

################################################################################
# save outputs of analyses to files

write.table(plotting_df, output_file1, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(post_probs_df, output_file1, quote = FALSE, sep = "\t", row.names = FALSE)
ggsave(output_file3, CI_plot)

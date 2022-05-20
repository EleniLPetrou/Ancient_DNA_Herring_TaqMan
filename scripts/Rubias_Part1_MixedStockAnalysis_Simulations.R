# The purpose of this script is to conduct simulations of mixed stock analysis 
# using a set of reference (aka "baseline" samples) and the R package rubias. 
# Data analyzed using "R version 3.6.1"

################################################################################
# Load libraries
library(rubias)
library(tidyverse)

################################################################################
# Specify the input file name (containing reference samples in rubias format)
ref_file <- "input_ref_pops_3reportinggroups.txt"

# Specify the output file names
output_file1 <- "results_rubias_analysis_self-ass_mixsim_3repgroups.txt"
output_file2 <- "results_rubias_analysis_self-ass_100sim_3repgroups.txt"

################################################################################
# Read in the data
ref_df <- read.delim(ref_file, header = TRUE, sep = "\t",  colClasses = c(rep("character", 4), rep("integer", 14)), na.strings = "NA")

################################################################################
# Simulated mixtures using a leave-one-out type of approach

# If you want to know how much accuracy you can expect given a set of genetic markers and a 
# grouping of populations (collections) into reporting units (repunits), use function
# assess_reference_loo(): This function carries out simulation of mixtures using the 
# leave-one-out approach of Anderson, Waples, and Kalinowski (2008).

herr_sims <- assess_reference_loo(reference = ref_df, 
                                  gen_start_col = 5, #the number of the column in which the genetic data start
                                  reps = 50, # simulate 50 mixture samples
                                  mixsize = 48) # each mixture sample contains 48 fish 


# Sum up the simulation results over reporting units
tmp <- herr_sims %>%
  group_by(iter, repunit) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n)) %>%
  mutate(repu_n_prop = repu_n / sum(repu_n))

################################################################################
# Plot the data

my_cols <- c("#2b83ba","#abdda4",  "#fdae61")

ggplot(tmp, aes(x = true_repprop, y = reprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Simulated mixing proportion") +
  ylab("Estimated mixing proportion")+
  facet_wrap(~ repunit) +
  scale_color_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr", "May-Jun"), values = my_cols)+
  theme_bw()


################################################################################
# Linear regression of the simulated data
head(tmp)

tmp_Jan_Feb <- filter(tmp, repunit == "Jan_Feb" )
linearmod_JanFeb <- lm(reprop_posterior_mean ~ true_repprop, data = tmp_Jan_Feb)
print(linearmod_JanFeb)
summary(linearmod_JanFeb)

tmp_Mar_Apr <- filter(tmp, repunit == "Mar_Apr" )
linearmod_MarApr <- lm(reprop_posterior_mean ~ true_repprop, data = tmp_Mar_Apr)
print(linearmod_MarApr)
summary(linearmod_MarApr)

tmp_May_Jun <- filter(tmp, repunit == "May_Jun" )
linearmod_MayJun <- lm(reprop_posterior_mean ~ true_repprop, data = tmp_May_Jun)
print(linearmod_MayJun)
summary(linearmod_MayJun)


############################################################################
# 100% SIMULATIONS
# Mixtures are simulated such that 100% of the individuals are from one reporting unit.

# Specify reporting group categories
my_reporting_groups <- unique(ref_df$repunit) 

# Specify mixture proportion - 100%
my_scenarios <- lapply(my_reporting_groups, function(x) tibble(repunit = x, ppn = 1.0))
my_scenarios #check

# Specify some scenario names
names(my_scenarios) <- paste("All", my_reporting_groups, sep = "-")
my_scenarios #check

# Do the 100% Simulations
scenario_results <- assess_reference_loo(reference = ref_df, 
                                         gen_start_col = 5, 
                                         reps = 100, 
                                         mixsize = 48,
                                         alpha_repunit = my_scenarios,
                                         alpha_collection = 10)

# Look at the results
head(scenario_results)

# Explanation of results:
  #repunit_scenario: and integer that gives that repunit simulation parameters 
  #collections_scenario: an integer that gives that collection simulation paramters 
  #iter: the simulation number (1 up to reps)
  #repunit: the reporting unit
  #collection: the collection
  #true_pi: the true simulated mixing proportion
  #n: the actual number of fish from the collection in the simulated mixture.
  #post_mean_pi: the posterior mean of mixing proportion.
  #mle_pi: the maximum likelihood of pi obtained using an EM-algorithm.


# Summarize these results for each reporting unit and scenario by summing 
# over all iterations and taking the mean reproductive posterior
tmp2 <- scenario_results %>%
  group_by(repunit, repunit_scenario, iter) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n))

head(tmp2)

# Make a boxplot of 100% simulation assignments
tmp2_100 <- filter(tmp2,true_repprop == 1) #subset the data

ggplot(tmp2_100, aes(x = repunit_scenario, y = reprop_posterior_mean))+
  geom_boxplot(aes(color = repunit))+
  ylim(0,1) +
  theme_bw()+
  xlab("Simulated mixing proportion") +
  ylab("Estimated mixing proportion")+
  scale_x_discrete(labels = c("100% Jan-Feb spawners", "100% Mar-Apr spawners","100% May-Jun spawners" ))+
  scale_color_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr", "May-Jun"), values = my_cols)
  
# Summarize the percent assignments per reporting group
result_df <- tmp2_100 %>%
  group_by(repunit) %>%
  summarise(mean_pi = mean(reprop_posterior_mean), sd_pi = sd(reprop_posterior_mean))

result_df

# Summarize the the percent assignments per reporting group,  include mis-assignment information
hundysim_results_by_group <- scenario_results %>%
  group_by(repunit, repunit_scenario) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n))


# Write the results to a file
write.table(tmp, output_file1, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(hundysim_results_by_group, output_file2, quote = FALSE, sep = "\t", row.names = FALSE)


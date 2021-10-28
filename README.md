# This R script imports a dataset from Exposito-Alonso et al. 2020 measuring survival and seed production of 517 genotypes of Arabidopsis thaliana grown in 8 different environmental condition.

# We separate the dataset by environmental condition and then analyze the dataset to count how many deaths of seeds and seedlings are selective, as well as how many 'deaths' from reduced fecundity as adults.
# We also calculate reproductive excess in each environmental condition.
# Selective deaths are determined by comparing the average deaths (or fecundity) of each genotype to the best genotype at that life history stage in that environmental condition.

# A substantial portion of the code is devoted to correcting for the extreme-value bias introduced by comparing to the best observed genotype in each life history stage and environmental condition.
# This correction is accomplished by fitting an appropriate distribution to average genotype survival or fecundity and then simulating redoing the experiment 10,00 times.
# In each simulation, we resample 517 genotypes from this distribution, then resample from a binomial (for survival) or Poisson (for fecundity) for every genotype, and compare the 'observed' best genotype to the 'true' best genotype.
# The average of this difference across 10,000 simulations is our estimate of the extreme-value bias introduced.

# We then adjust our selective death calculations to incorporate this bias, although it has only a small effect.

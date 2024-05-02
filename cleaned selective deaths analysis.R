library('ggplot2')
library('ggpubr')
library('egg')
library('MASS')
library('scales')
install.packages("scales")
install.packages("ggpubr")
getwd()
setwd('..')

options(scipen = 10000)

#FUNCTIONS

#This function is taken from the internet, it may not be modular.
#Source: https://groups.google.com/g/ggplot2/c/a_xhMoQyxZ4
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2)
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

#This function requires a dataset that is NOT aggregated
#and has replicates with no surviving seeds OMITTED
CalculateVarianceStDevandMeaninSeedProduction <- function(nozeroesdata) {
  nozeroesdata$Variance <- ave(nozeroesdata$Seeds_by_ind, nozeroesdata$Genotype_id, FUN = var)
  nozeroesdata$StDev <- sqrt(nozeroesdata$Variance)
  nozeroesdata$Mean <- ave(nozeroesdata$Seeds_by_ind, nozeroesdata$Genotype_id, FUN = mean)
  nozeroesdata[is.na(nozeroesdata)] = 0.0
  return(nozeroesdata)
}

#I need to add some new columns to the dataset for the main calculation loop.
#I'm also going to set the Replicates column in this function to INCLUDE pots with no surviving seeds.
AddColumns <- function(environmentalconditiondataframe, nreplicatestable) {
  environmentalconditiondataframe$Replicates <- c(1:nrow(environmentalconditiondataframe))
  environmentalconditiondataframe$Selective_deaths_survival <- c(1:nrow(environmentalconditiondataframe))
  environmentalconditiondataframe$Selective_deaths_births <- c(1:nrow(environmentalconditiondataframe))
  environmentalconditiondataframe$Genotype_death_rate <- c(1:nrow(environmentalconditiondataframe))
  environmentalconditiondataframe$Seeds_for_next_generation <- c(1:nrow(environmentalconditiondataframe))
  for (i in 1:nrow(environmentalconditiondataframe)) {
    environmentalconditiondataframe$Replicates[i] <- nreplicatestable[i]
  }
  return(environmentalconditiondataframe)
}

Boxcoxtransform <- function(environmentalconditiondataframe, lambda) {
  environmentalconditiondataframe$boxcox_transform_fecundity <- c(1:nrow(environmentalconditiondataframe))
  for (i in 1:nrow(environmentalconditiondataframe)) {
    environmentalconditiondataframe$boxcox_transform_fecundity[i] <- (environmentalconditiondataframe$Seeds_by_ind[i]^lambda - 1) / lambda
  }
  return(environmentalconditiondataframe)
}

#This function should take in an aggregated dataset where replicates with
#no surviving seeds are INCLUDED
CalculateSelectiveDeathsSurvival <- function(environmentalconditiondataframe, nreplicatestable, isPopulation) {
  for (i in 1:nrow(environmentalconditiondataframe)) {
    environmentalconditiondataframe$Surviving_individuals_pot[i] <- (environmentalconditiondataframe$Surviving_individuals_pot[i] * environmentalconditiondataframe$Replicates[i])
    if (isPopulation == 1) {
      environmentalconditiondataframe$Genotype_death_rate[i] <- (30*environmentalconditiondataframe$Replicates[i] - environmentalconditiondataframe$Surviving_individuals_pot[i]) / (30*environmentalconditiondataframe$Replicates[i])
    } else {
      environmentalconditiondataframe$Genotype_death_rate[i] <- (environmentalconditiondataframe$Replicates[i] - environmentalconditiondataframe$Surviving_individuals_pot[i]) / (environmentalconditiondataframe$Replicates[i])
    }
  }
  min_death_rate <- min(environmentalconditiondataframe$Genotype_death_rate)
  for (i in 1:nrow(environmentalconditiondataframe)) {
    if (isPopulation == 1) {
      environmentalconditiondataframe$Selective_deaths_survival[i] <- (30*environmentalconditiondataframe$Replicates[i]) * (environmentalconditiondataframe$Genotype_death_rate[i] - min_death_rate)
    } else {
      environmentalconditiondataframe$Selective_deaths_survival[i] <- (environmentalconditiondataframe$Replicates[i]) * (environmentalconditiondataframe$Genotype_death_rate[i] - min_death_rate)
    }
  }
  colnames(environmentalconditiondataframe)[colnames(environmentalconditiondataframe) == "Surviving_individuals_pot"] <- "Total_surviving_seeds"
  return(environmentalconditiondataframe)
}

CalculateSelectiveDeathsSurvivalWithBias <- function(environmentalconditiondataframe, nreplicatestable, meanbias, isPopulation) {
  for (i in 1:nrow(environmentalconditiondataframe)) {
    environmentalconditiondataframe$Surviving_individuals_pot[i] <- (environmentalconditiondataframe$Surviving_individuals_pot[i] * environmentalconditiondataframe$Replicates[i])
    if (isPopulation == 1) {
      environmentalconditiondataframe$Genotype_death_rate[i] <- (30*environmentalconditiondataframe$Replicates[i] - environmentalconditiondataframe$Surviving_individuals_pot[i]) / (30*environmentalconditiondataframe$Replicates[i])
    } else {
      environmentalconditiondataframe$Genotype_death_rate[i] <- (environmentalconditiondataframe$Replicates[i] - environmentalconditiondataframe$Surviving_individuals_pot[i]) / (environmentalconditiondataframe$Replicates[i])
    }
  }
  min_death_rate <- min(environmentalconditiondataframe$Genotype_death_rate)
  for (i in 1:nrow(environmentalconditiondataframe)) {
    if (isPopulation == 1) {
      environmentalconditiondataframe$Selective_deaths_survival[i] <- (30*environmentalconditiondataframe$Replicates[i]) * (environmentalconditiondataframe$Genotype_death_rate[i] - (min_death_rate + meanbias))
    } else {
      environmentalconditiondataframe$Selective_deaths_survival[i] <- (environmentalconditiondataframe$Replicates[i]) * (environmentalconditiondataframe$Genotype_death_rate[i] - (min_death_rate + meanbias))
    }
    if (environmentalconditiondataframe$Selective_deaths_survival[i] < 0.0) {
      environmentalconditiondataframe$Selective_deaths_survival[i] <- 0.0
    }
  }
  colnames(environmentalconditiondataframe)[colnames(environmentalconditiondataframe) == "Surviving_individuals_pot"] <- "Total_surviving_seeds"
  return(environmentalconditiondataframe)
}

#This function should take in an aggregated dataset where replicates with
#no surviving seeds are OMITTED (MUST HAVE BEEN OMITTED FIRST THEN AGGREGATED)
#It also needs a replicates table that has OMITTED replicates with no surviving seeds.
#Note that the $Replicates column INCLUDES these replicates.
CalculateSelectiveDeathsBirths <- function(environmentalconditiondataframe, nreplicatestablewithoutpotswithnosurvivingseeds, isPopulation) {
  for (i in 1:nrow(environmentalconditiondataframe)) {
    environmentalconditiondataframe$Surviving_individuals_pot[i] <- (environmentalconditiondataframe$Surviving_individuals_pot[i] * nreplicatestablewithoutpotswithnosurvivingseeds[i])
    environmentalconditiondataframe$Seeds_total_pot[i] <- (environmentalconditiondataframe$Seeds_total_pot[i] * nreplicatestablewithoutpotswithnosurvivingseeds[i])
    if (environmentalconditiondataframe$Surviving_individuals_pot[i] > 0.0) {
      environmentalconditiondataframe$Seeds_by_ind[i] <- environmentalconditiondataframe$Seeds_total_pot[i] / environmentalconditiondataframe$Surviving_individuals_pot[i]
    } else {
      environmentalconditiondataframe$Seeds_by_ind[i] <- 0.0
    }
    if (isPopulation == 1) {
      if (environmentalconditiondataframe$Seeds_total_pot[i] > 30.0*environmentalconditiondataframe$Replicates[i]) {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- 30.0*environmentalconditiondataframe$Replicates[i]
      } else {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- environmentalconditiondataframe$Seeds_total_pot[i]
      }
    } else {
      if (environmentalconditiondataframe$Seeds_total_pot[i] > 10.0*environmentalconditiondataframe$Replicates[i]) {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- 10.0*environmentalconditiondataframe$Replicates[i]
      } else {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- environmentalconditiondataframe$Seeds_total_pot[i]
      }
    }
  }
  max_birth_rate <- max(environmentalconditiondataframe$Seeds_by_ind)
  for (i in 1:nrow(environmentalconditiondataframe)) {
    environmentalconditiondataframe$Selective_deaths_births[i] <- (max_birth_rate - environmentalconditiondataframe$Seeds_by_ind[i]) * environmentalconditiondataframe$Surviving_individuals_pot[i]
  }
  return(environmentalconditiondataframe)
}


CalculateSelectiveDeathsBirthsWithBias <- function(environmentalconditiondataframe, nreplicatestablewithoutpotswithnosurvivingseeds, meanbias, isPopulation) {
  for (i in 1:nrow(environmentalconditiondataframe)) {
    environmentalconditiondataframe$Surviving_individuals_pot[i] <- (environmentalconditiondataframe$Surviving_individuals_pot[i] * nreplicatestablewithoutpotswithnosurvivingseeds[i])
    environmentalconditiondataframe$Seeds_total_pot[i] <- (environmentalconditiondataframe$Seeds_total_pot[i] * nreplicatestablewithoutpotswithnosurvivingseeds[i])
    if (environmentalconditiondataframe$Surviving_individuals_pot[i] > 0.0) {
      environmentalconditiondataframe$Seeds_by_ind[i] <- environmentalconditiondataframe$Seeds_total_pot[i] / environmentalconditiondataframe$Surviving_individuals_pot[i]
    } else {
      environmentalconditiondataframe$Seeds_by_ind[i] <- 0.0
    }
    if (isPopulation == 1) {
      if (environmentalconditiondataframe$Seeds_total_pot[i] > 30.0*environmentalconditiondataframe$Replicates[i]) {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- 30.0*environmentalconditiondataframe$Replicates[i]
      } else {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- environmentalconditiondataframe$Seeds_total_pot[i]
      }
    } else {
      if (environmentalconditiondataframe$Seeds_total_pot[i] > 10.0*environmentalconditiondataframe$Replicates[i]) {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- 10.0*environmentalconditiondataframe$Replicates[i]
      } else {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- environmentalconditiondataframe$Seeds_total_pot[i]
      }
    }
  }
  max_birth_rate <- max(environmentalconditiondataframe$Seeds_by_ind)
  for (i in 1:nrow(environmentalconditiondataframe)) {
    environmentalconditiondataframe$Selective_deaths_births[i] <- ((max_birth_rate - meanbias) - environmentalconditiondataframe$Seeds_by_ind[i]) * environmentalconditiondataframe$Surviving_individuals_pot[i]
    if (environmentalconditiondataframe$Selective_deaths_births[i] < 0.0) {
      environmentalconditiondataframe$Selective_deaths_births[i] <- 0.0
    }
  }
  return(environmentalconditiondataframe)
}


ExtremeValueBiasSurvival <- function(environmentalconditiondataframe, isPopulation) {
  simulationdataframe <- data.frame(
    genotypesurvivalrates = environmentalconditiondataframe$Genotype_death_rate,
    genotypeID = environmentalconditiondataframe$Genotype_id,
    resampledgenotypesurvivalrates = environmentalconditiondataframe$Genotype_death_rate
  )
  resultsdataframe <- data.frame(
    vectorofmaximumsurvivalrates = c(1.0:10000.0),
    vectoroftruesurvivalrates = c(1.0:10000.0),
    vectorofwithingenotypevariances = c(1.0:10000.0),
    vectorofbestgenotypeID = c(1.0:10000.0)
  )
  for (i in 1:10000) {
    for(j in 1:nrow(environmentalconditiondataframe)) {
      simulationdataframe$genotypesurvivalrates[j] <- (1.0 - environmentalconditiondataframe$Genotype_death_rate[j])
    }
    if(isPopulation == 0) {
      for(j in 1:nrow(environmentalconditiondataframe)) {
        simulationdataframe$resampledgenotypesurvivalrates[j] <- (rbinom(1, (environmentalconditiondataframe$Replicates[j]), simulationdataframe$genotypesurvivalrates[j]) / (environmentalconditiondataframe$Replicates[j]))
      }
    } else {
      for(j in 1:nrow(environmentalconditiondataframe)) {
        simulationdataframe$resampledgenotypesurvivalrates[j] <- (rbinom(1, (environmentalconditiondataframe$Replicates[j])*30, simulationdataframe$genotypesurvivalrates[j]) / (environmentalconditiondataframe$Replicates[j] * 30))
      }
    }
    resultsdataframe$vectorofmaximumsurvivalrates[i] <- max(simulationdataframe$resampledgenotypesurvivalrates)
    resultsdataframe$vectoroftruesurvivalrates[i] <- max(simulationdataframe$genotypesurvivalrates)
    resultsdataframe$vectorofbestgenotypeID[i] <- simulationdataframe$genotypeID[which.max(simulationdataframe$genotypesurvivalrates)]
    if (isPopulation == 0) {
      resultsdataframe$vectorofwithingenotypevariances[i] <- resultsdataframe$vectoroftruesurvivalrates[i]*(1 - resultsdataframe$vectoroftruesurvivalrates[i]) / (environmentalconditiondataframe$Replicates[which(environmentalconditiondataframe$Genotype_id == resultsdataframe$vectorofbestgenotypeID[i], arr.ind = FALSE)])
    } else {
      resultsdataframe$vectorofwithingenotypevariances[i] <- resultsdataframe$vectoroftruesurvivalrates[i]*(1 - resultsdataframe$vectoroftruesurvivalrates[i]) / (environmentalconditiondataframe$Replicates[which(environmentalconditiondataframe$Genotype_id == resultsdataframe$vectorofbestgenotypeID[i], arr.ind = FALSE)]*30)
    }
  }
  biasvector <- c(1:10000)
  for (i in 1:10000) {
    biasvector[i] = (resultsdataframe$vectorofmaximumsurvivalrates[i] - resultsdataframe$vectoroftruesurvivalrates[i])
  }
  meanbias <- mean(biasvector)
  return(meanbias)
}


#DATASET IMPORT AND ORGANIZING

#imports full dataset
fullselectivedeathsdata = read.table("Selective deaths/Dataset_Arabidopsis_Exposito-Alonso.txt", header = TRUE)

#sets NA values to zero
#since all NA values are in seed number data entries, they always mean zero seeds found
fullselectivedeathsdata[is.na(fullselectivedeathsdata)] = 0

#there might be a more elegant way to subsection the eight treatment types,
#but since it's a 2x2x2 design with eight combinations of three different variables,
#subsectioning for each variable is the first way I came up with

#subsectioning by seed density treatment (1 seed per pot vs. 20 seeds per pot)
oneseedselectivedeathsdata <- subset(fullselectivedeathsdata, Density_treatment == 'i', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
twentyseedsselectivedeathsdata <- subset(fullselectivedeathsdata, Density_treatment == 'p', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))

#subsectioning seed density datasets by watering treatment (low vs. high)
lowwateroneseedselectivedeathsdata <- subset(oneseedselectivedeathsdata, Watering_treatment == 'l', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
highwateroneseedselectivedeathsdata <- subset(oneseedselectivedeathsdata, Watering_treatment == 'h', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))

lowwatertwentyseedsselectivedeathsdata <- subset(twentyseedsselectivedeathsdata, Watering_treatment == 'l', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
highwatertwentyseedsselectivedeathsdata <- subset(twentyseedsselectivedeathsdata, Watering_treatment == 'h', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))

#subsectioning waterxdensity datasets by climate (madrid vs. tuebingen)
madridlowwateroneseedselectivedeathsdata <- subset(lowwateroneseedselectivedeathsdata, Field_site == 'madrid', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
madridhighwateroneseedselectivedeathsdata <- subset(highwateroneseedselectivedeathsdata, Field_site == 'madrid', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))

madridlowwatertwentyseedsselectivedeathsdata <- subset(lowwatertwentyseedsselectivedeathsdata, Field_site == 'madrid', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
madridhighwatertwentyseedsselectivedeathsdata <- subset(highwatertwentyseedsselectivedeathsdata, Field_site == 'madrid', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))

tuebingenlowwateroneseedselectivedeathsdata <- subset(lowwateroneseedselectivedeathsdata, Field_site == 'tuebingen', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
tuebingenhighwateroneseedselectivedeathsdata <- subset(highwateroneseedselectivedeathsdata, Field_site == 'tuebingen', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))

tuebingenlowwatertwentyseedsselectivedeathsdata <- subset(lowwatertwentyseedsselectivedeathsdata, Field_site == 'tuebingen', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
tuebingenhighwatertwentyseedsselectivedeathsdata <- subset(highwatertwentyseedsselectivedeathsdata, Field_site == 'tuebingen', select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))

#I need to calculate variance in seed production per individual for each genotype
#for the purpose of correcting for extreme value bias later.
#This means I need to drop all replicates with no surviving individuals.
onlypotswithsurvivingseedsMLO <- subset(madridlowwateroneseedselectivedeathsdata, Surviving_individuals_pot > 0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
onlypotswithsurvivingseedsMLP <- subset(madridlowwatertwentyseedsselectivedeathsdata, Surviving_individuals_pot > 0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
onlypotswithsurvivingseedsMHO <- subset(madridhighwateroneseedselectivedeathsdata, Surviving_individuals_pot > 0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
onlypotswithsurvivingseedsMHP <- subset(madridhighwatertwentyseedsselectivedeathsdata, Surviving_individuals_pot > 0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
onlypotswithsurvivingseedsTLO <- subset(tuebingenlowwateroneseedselectivedeathsdata, Surviving_individuals_pot > 0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
onlypotswithsurvivingseedsTLP <- subset(tuebingenlowwatertwentyseedsselectivedeathsdata, Surviving_individuals_pot > 0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
onlypotswithsurvivingseedsTHO <- subset(tuebingenhighwateroneseedselectivedeathsdata, Surviving_individuals_pot > 0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))
onlypotswithsurvivingseedsTHP <- subset(tuebingenhighwatertwentyseedsselectivedeathsdata, Surviving_individuals_pot > 0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind"))


#Note that the datasets above have replicates for each genotype.
#In order to calculate selective deaths, need to know the optimal genotype at each LH stage
#This requires aggregating data from replicates of each genotype
#INCLUDES POTS WITH NO SURVIVING SEEDS IN SURVIVAL AND BIRTHS
aggregatedpotswithnosurvivingseedsMLO <- aggregate(madridlowwateroneseedselectivedeathsdata, by = list(madridlowwateroneseedselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedpotswithnosurvivingseedsMLP <- aggregate(madridlowwatertwentyseedsselectivedeathsdata, by = list(madridlowwatertwentyseedsselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedpotswithnosurvivingseedsMHO <- aggregate(madridhighwateroneseedselectivedeathsdata, by = list(madridhighwateroneseedselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedpotswithnosurvivingseedsMHP <- aggregate(madridhighwatertwentyseedsselectivedeathsdata, by = list(madridhighwatertwentyseedsselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedpotswithnosurvivingseedsTLO <- aggregate(tuebingenlowwateroneseedselectivedeathsdata, by = list(tuebingenlowwateroneseedselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedpotswithnosurvivingseedsTLP <- aggregate(tuebingenlowwatertwentyseedsselectivedeathsdata, by = list(tuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedpotswithnosurvivingseedsTHO <- aggregate(tuebingenhighwateroneseedselectivedeathsdata, by = list(tuebingenhighwateroneseedselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedpotswithnosurvivingseedsTHP <- aggregate(tuebingenhighwatertwentyseedsselectivedeathsdata, by = list(tuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_id), FUN = mean)

#Same aggregation procedure for data where replicates with no surviving seeds
#are OMITTED.
aggregatedonlypotswithsurvivingseedsMLO <- aggregate(onlypotswithsurvivingseedsMLO, by = list(onlypotswithsurvivingseedsMLO$Genotype_id), FUN = mean)
aggregatedonlypotswithsurvivingseedsMLP <- aggregate(onlypotswithsurvivingseedsMLP, by = list(onlypotswithsurvivingseedsMLP$Genotype_id), FUN = mean)
aggregatedonlypotswithsurvivingseedsMHO <- aggregate(onlypotswithsurvivingseedsMHO, by = list(onlypotswithsurvivingseedsMHO$Genotype_id), FUN = mean)
aggregatedonlypotswithsurvivingseedsMHP <- aggregate(onlypotswithsurvivingseedsMHP, by = list(onlypotswithsurvivingseedsMHP$Genotype_id), FUN = mean)
aggregatedonlypotswithsurvivingseedsTLO <- aggregate(onlypotswithsurvivingseedsTLO, by = list(onlypotswithsurvivingseedsTLO$Genotype_id), FUN = mean)
aggregatedonlypotswithsurvivingseedsTLP <- aggregate(onlypotswithsurvivingseedsTLP, by = list(onlypotswithsurvivingseedsTLP$Genotype_id), FUN = mean)
aggregatedonlypotswithsurvivingseedsTHO <- aggregate(onlypotswithsurvivingseedsTHO, by = list(onlypotswithsurvivingseedsTHO$Genotype_id), FUN = mean)
aggregatedonlypotswithsurvivingseedsTHP <- aggregate(onlypotswithsurvivingseedsTHP, by = list(onlypotswithsurvivingseedsTHP$Genotype_id), FUN = mean)


#After the data has been aggregated into averages by replicate,
#I still need to know the number of replicates to use in the equation for selective deaths.
#This is easily accomplished in a frequency table called on the non-aggregated data
nreplicatestableMLO <- table(madridlowwateroneseedselectivedeathsdata$Genotype_id)
nreplicatestableMHO <- table(madridhighwateroneseedselectivedeathsdata$Genotype_id)
nreplicatestableMLP <- table(madridlowwatertwentyseedsselectivedeathsdata$Genotype_id)
nreplicatestableMHP <- table(madridhighwatertwentyseedsselectivedeathsdata$Genotype_id)
nreplicatestableTLO <- table(tuebingenlowwateroneseedselectivedeathsdata$Genotype_id)
nreplicatestableTHO <- table(tuebingenhighwateroneseedselectivedeathsdata$Genotype_id)
nreplicatestableTLP <- table(tuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_id)
nreplicatestableTHP <- table(tuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_id)

#I also need a table with the number of replicates AFTER OMITTING replicates with no surviving seeds
nreplicatestableomittingreplicateswithnosurvivingseedsMLO <- table(onlypotswithsurvivingseedsMLO$Genotype_id)
nreplicatestableomittingreplicateswithnosurvivingseedsMLP <- table(onlypotswithsurvivingseedsMLP$Genotype_id)
nreplicatestableomittingreplicateswithnosurvivingseedsMHO <- table(onlypotswithsurvivingseedsMHO$Genotype_id)
nreplicatestableomittingreplicateswithnosurvivingseedsMHP <- table(onlypotswithsurvivingseedsMHP$Genotype_id)
nreplicatestableomittingreplicateswithnosurvivingseedsTLO <- table(onlypotswithsurvivingseedsTLO$Genotype_id)
nreplicatestableomittingreplicateswithnosurvivingseedsTLP <- table(onlypotswithsurvivingseedsTLP$Genotype_id)
nreplicatestableomittingreplicateswithnosurvivingseedsTHO <- table(onlypotswithsurvivingseedsTHO$Genotype_id)
nreplicatestableomittingreplicateswithnosurvivingseedsTHP <- table(onlypotswithsurvivingseedsTHP$Genotype_id)

#I also need to pull out the total number of replicates (including those with no surviving seeds)
#but ONLY for genotypes with at least one replicate with surviving seeds.
onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsMLO <- subset(aggregatedwithcolumnsnosurvivingseedsMLO, Surviving_individuals_pot > 0.0, select = c("Replicates"))
onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsMLP <- subset(aggregatedwithcolumnsnosurvivingseedsMLP, Surviving_individuals_pot > 0.0, select = c("Replicates"))
onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsMHO <- subset(aggregatedwithcolumnsnosurvivingseedsMHO, Surviving_individuals_pot > 0.0, select = c("Replicates"))
onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsMHP <- subset(aggregatedwithcolumnsnosurvivingseedsMHP, Surviving_individuals_pot > 0.0, select = c("Replicates"))
onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsTLO <- subset(aggregatedwithcolumnsnosurvivingseedsTLO, Surviving_individuals_pot > 0.0, select = c("Replicates"))
onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsTLP <- subset(aggregatedwithcolumnsnosurvivingseedsTLP, Surviving_individuals_pot > 0.0, select = c("Replicates"))
onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsTHO <- subset(aggregatedwithcolumnsnosurvivingseedsTHO, Surviving_individuals_pot > 0.0, select = c("Replicates"))
onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsTHP <- subset(aggregatedwithcolumnsnosurvivingseedsTHP, Surviving_individuals_pot > 0.0, select = c("Replicates"))


#Now add columns to the aggregated data.
#This function only takes the aggregated data including replicates with no surviving seeds
aggregatedwithcolumnsnosurvivingseedsMLO <- AddColumns(aggregatedpotswithnosurvivingseedsMLO, nreplicatestableMLO)
aggregatedwithcolumnsnosurvivingseedsMLP <- AddColumns(aggregatedpotswithnosurvivingseedsMLP, nreplicatestableMLP)
aggregatedwithcolumnsnosurvivingseedsMHO <- AddColumns(aggregatedpotswithnosurvivingseedsMHO, nreplicatestableMHO)
aggregatedwithcolumnsnosurvivingseedsMHP <- AddColumns(aggregatedpotswithnosurvivingseedsMHP, nreplicatestableMHP)
aggregatedwithcolumnsnosurvivingseedsTLO <- AddColumns(aggregatedpotswithnosurvivingseedsTLO, nreplicatestableTLO)
aggregatedwithcolumnsnosurvivingseedsTLP <- AddColumns(aggregatedpotswithnosurvivingseedsTLP, nreplicatestableTLP)
aggregatedwithcolumnsnosurvivingseedsTHO <- AddColumns(aggregatedpotswithnosurvivingseedsTHO, nreplicatestableTHO)
aggregatedwithcolumnsnosurvivingseedsTHP <- AddColumns(aggregatedpotswithnosurvivingseedsTHP, nreplicatestableTHP)

#Now add columns to aggregated data which OMITS replicates with no surviving seeds.
#The number in the $Replicates column will be the total number of replicates for that genotype,
#INCLUDING replicates with no surviving seeds.
aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMLO <- AddColumns(aggregatedonlypotswithsurvivingseedsMLO, onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsMLO$Replicates)
aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMLP <- AddColumns(aggregatedonlypotswithsurvivingseedsMLP, onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsMLP$Replicates)
aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMHO <- AddColumns(aggregatedonlypotswithsurvivingseedsMHO, onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsMHO$Replicates)
aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMHP <- AddColumns(aggregatedonlypotswithsurvivingseedsMHP, onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsMHP$Replicates)
aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTLO <- AddColumns(aggregatedonlypotswithsurvivingseedsTLO, onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsTLO$Replicates)
aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTLP <- AddColumns(aggregatedonlypotswithsurvivingseedsTLP, onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsTLP$Replicates)
aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTHO <- AddColumns(aggregatedonlypotswithsurvivingseedsTHO, onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsTHO$Replicates)
aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTHP <- AddColumns(aggregatedonlypotswithsurvivingseedsTHP, onlygenotypeswithsomesurvivingseedsbutincludingreplicateswithnosurvivingseedsTHP$Replicates)

#Following lines create datasets to compare the standard deviation within genotypes.
onlypotswithsurvivingseedsstdevforbirthsMLO <- CalculateVarianceStDevandMeaninSeedProduction(onlypotswithsurvivingseedsMLO)
onlypotswithsurvivingseedsstdevforbirthsMLP <- CalculateVarianceStDevandMeaninSeedProduction(onlypotswithsurvivingseedsMLP)
onlypotswithsurvivingseedsstdevforbirthsMHO <- CalculateVarianceStDevandMeaninSeedProduction(onlypotswithsurvivingseedsMHO)
onlypotswithsurvivingseedsstdevforbirthsMHP <- CalculateVarianceStDevandMeaninSeedProduction(onlypotswithsurvivingseedsMHP)
onlypotswithsurvivingseedsstdevforbirthsTLO <- CalculateVarianceStDevandMeaninSeedProduction(onlypotswithsurvivingseedsTLO)
onlypotswithsurvivingseedsstdevforbirthsTLP <- CalculateVarianceStDevandMeaninSeedProduction(onlypotswithsurvivingseedsTLP)
onlypotswithsurvivingseedsstdevforbirthsTHO <- CalculateVarianceStDevandMeaninSeedProduction(onlypotswithsurvivingseedsTHO)
onlypotswithsurvivingseedsstdevforbirthsTHP <- CalculateVarianceStDevandMeaninSeedProduction(onlypotswithsurvivingseedsTHP)

#The following lines remove genotypes with a variance of zero.
onlysurvivingseedsnozeroesinvarianceMLO <- subset(onlypotswithsurvivingseedsstdevforbirthsMLO, Variance > 0.0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Variance", "Mean", "StDev"))
onlysurvivingseedsnozeroesinvarianceMLP <- subset(onlypotswithsurvivingseedsstdevforbirthsMLP, Variance > 0.0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Variance", "Mean", "StDev"))
onlysurvivingseedsnozeroesinvarianceMHO <- subset(onlypotswithsurvivingseedsstdevforbirthsMHO, Variance > 0.0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Variance", "Mean", "StDev"))
onlysurvivingseedsnozeroesinvarianceMHP <- subset(onlypotswithsurvivingseedsstdevforbirthsMHP, Variance > 0.0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Variance", "Mean", "StDev"))
onlysurvivingseedsnozeroesinvarianceTLO <- subset(onlypotswithsurvivingseedsstdevforbirthsTLO, Variance > 0.0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Variance", "Mean", "StDev"))
onlysurvivingseedsnozeroesinvarianceTLP <- subset(onlypotswithsurvivingseedsstdevforbirthsTLP, Variance > 0.0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Variance", "Mean", "StDev"))
onlysurvivingseedsnozeroesinvarianceTHO <- subset(onlypotswithsurvivingseedsstdevforbirthsTHO, Variance > 0.0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Variance", "Mean", "StDev"))
onlysurvivingseedsnozeroesinvarianceTHP <- subset(onlypotswithsurvivingseedsstdevforbirthsTHP, Variance > 0.0, select = c("Genotype_id", "Field_site", "Watering_treatment", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Variance", "Mean", "StDev"))

#Now aggregate the resulting dataset to ensure that genotypes with more replicates don't get more represented in the regression.
aggregatedonlysurvivingseedsnozeroesinvarianceMLO <- aggregate(onlysurvivingseedsnozeroesinvarianceMLO, by = list(onlysurvivingseedsnozeroesinvarianceMLO$Genotype_id), FUN = mean)
aggregatedonlysurvivingseedsnozeroesinvarianceMLP <- aggregate(onlysurvivingseedsnozeroesinvarianceMLP, by = list(onlysurvivingseedsnozeroesinvarianceMLP$Genotype_id), FUN = mean)
aggregatedonlysurvivingseedsnozeroesinvarianceMHO <- aggregate(onlysurvivingseedsnozeroesinvarianceMHO, by = list(onlysurvivingseedsnozeroesinvarianceMHO$Genotype_id), FUN = mean)
aggregatedonlysurvivingseedsnozeroesinvarianceMHP <- aggregate(onlysurvivingseedsnozeroesinvarianceMHP, by = list(onlysurvivingseedsnozeroesinvarianceMHP$Genotype_id), FUN = mean)
aggregatedonlysurvivingseedsnozeroesinvarianceTLO <- aggregate(onlysurvivingseedsnozeroesinvarianceTLO, by = list(onlysurvivingseedsnozeroesinvarianceTLO$Genotype_id), FUN = mean)
aggregatedonlysurvivingseedsnozeroesinvarianceTLP <- aggregate(onlysurvivingseedsnozeroesinvarianceTLP, by = list(onlysurvivingseedsnozeroesinvarianceTLP$Genotype_id), FUN = mean)
aggregatedonlysurvivingseedsnozeroesinvarianceTHO <- aggregate(onlysurvivingseedsnozeroesinvarianceTHO, by = list(onlysurvivingseedsnozeroesinvarianceTHO$Genotype_id), FUN = mean)
aggregatedonlysurvivingseedsnozeroesinvarianceTHP <- aggregate(onlysurvivingseedsnozeroesinvarianceTHP, by = list(onlysurvivingseedsnozeroesinvarianceTHP$Genotype_id), FUN = mean)

max(aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Seeds_by_ind)

#For the extreme-value bias analysis on the life history stage of seed production,
#I need to set the standard deviation for within-genotype variation as a function
#of the mean. I calculate that slope by doing a simple regression.
#I require that the line passes through 0,0 -- this simplifies the later math
#and also seems like a more realistic assumption for the regression.
#OMITTED POTS WITH NO SURVIVING SEEDS AND GENOTYPES WITH NO VARIANCE
#MAKE SURE THE DATASET IS THE AGGREGATED VERSION
lm(StDev ~ Mean + 0, data = aggregatedonlysurvivingseedsnozeroesinvarianceMLO)
lm(StDev ~ Mean + 0, data = aggregatedonlysurvivingseedsnozeroesinvarianceMLP)
lm(StDev ~ Mean + 0, data = aggregatedonlysurvivingseedsnozeroesinvarianceMHO)
lm(StDev ~ Mean + 0, data = aggregatedonlysurvivingseedsnozeroesinvarianceMHP)
lm(StDev ~ Mean + 0, data = aggregatedonlysurvivingseedsnozeroesinvarianceTLO)
lm(StDev ~ Mean + 0, data = aggregatedonlysurvivingseedsnozeroesinvarianceTLP)
lm(StDev ~ Mean + 0, data = aggregatedonlysurvivingseedsnozeroesinvarianceTHO)
lm(StDev ~ Mean + 0, data = aggregatedonlysurvivingseedsnozeroesinvarianceTHP)


#ANALYSES

#Following are the basic selective death calculations,
#without accounting for extreme-value bias.
survivalselectivedeathsMLO <- CalculateSelectiveDeathsSurvival(aggregatedwithcolumnsnosurvivingseedsMLO, nreplicatestableMLO, 0)
survivalselectivedeathsMLP <- CalculateSelectiveDeathsSurvival(aggregatedwithcolumnsnosurvivingseedsMLP, nreplicatestableMLP, 1)
survivalselectivedeathsMHO <- CalculateSelectiveDeathsSurvival(aggregatedwithcolumnsnosurvivingseedsMHO, nreplicatestableMHO, 0)
survivalselectivedeathsMHP <- CalculateSelectiveDeathsSurvival(aggregatedwithcolumnsnosurvivingseedsMHP, nreplicatestableMHP, 1)
survivalselectivedeathsTLO <- CalculateSelectiveDeathsSurvival(aggregatedwithcolumnsnosurvivingseedsTLO, nreplicatestableTLO, 0)
survivalselectivedeathsTLP <- CalculateSelectiveDeathsSurvival(aggregatedwithcolumnsnosurvivingseedsTLP, nreplicatestableTLP, 1)
survivalselectivedeathsTHO <- CalculateSelectiveDeathsSurvival(aggregatedwithcolumnsnosurvivingseedsTHO, nreplicatestableTHO, 0)
survivalselectivedeathsTHP <- CalculateSelectiveDeathsSurvival(aggregatedwithcolumnsnosurvivingseedsTHP, nreplicatestableTHP, 1)

birthsselectivedeathsMLO <- CalculateSelectiveDeathsBirths(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMLO, nreplicatestableomittingreplicateswithnosurvivingseedsMLO, 0)
birthsselectivedeathsMLP <- CalculateSelectiveDeathsBirths(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMLP, nreplicatestableomittingreplicateswithnosurvivingseedsMLP, 1)
birthsselectivedeathsMHO <- CalculateSelectiveDeathsBirths(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMHO, nreplicatestableomittingreplicateswithnosurvivingseedsMHO, 0)
birthsselectivedeathsMHP <- CalculateSelectiveDeathsBirths(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMHP, nreplicatestableomittingreplicateswithnosurvivingseedsMHP, 1)
birthsselectivedeathsTLO <- CalculateSelectiveDeathsBirths(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTLO, nreplicatestableomittingreplicateswithnosurvivingseedsTLO, 0)
birthsselectivedeathsTLP <- CalculateSelectiveDeathsBirths(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTLP, nreplicatestableomittingreplicateswithnosurvivingseedsTLP, 1)
birthsselectivedeathsTHO <- CalculateSelectiveDeathsBirths(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTHO, nreplicatestableomittingreplicateswithnosurvivingseedsTHO, 0)
birthsselectivedeathsTHP <- CalculateSelectiveDeathsBirths(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTHP, nreplicatestableomittingreplicateswithnosurvivingseedsTHP, 1)


total_selective_deaths_survival_MLO <- sum(survivalselectivedeathsMLO$Selective_deaths_survival)
total_selective_deaths_survival_MLP <- sum(survivalselectivedeathsMLP$Selective_deaths_survival)
total_selective_deaths_survival_MHO <- sum(survivalselectivedeathsMHO$Selective_deaths_survival)
total_selective_deaths_survival_MHP <- sum(survivalselectivedeathsMHP$Selective_deaths_survival)
total_selective_deaths_survival_TLO <- sum(survivalselectivedeathsTLO$Selective_deaths_survival)
total_selective_deaths_survival_TLP <- sum(survivalselectivedeathsTLP$Selective_deaths_survival)
total_selective_deaths_survival_THO <- sum(survivalselectivedeathsTHO$Selective_deaths_survival)
total_selective_deaths_survival_THP <- sum(survivalselectivedeathsTHP$Selective_deaths_survival)
total_selective_deaths_survival_MLO
total_selective_deaths_survival_MLP

total_selective_deaths_births_MLO <- sum(birthsselectivedeathsMLO$Selective_deaths_births)
total_selective_deaths_births_MLP <- sum(birthsselectivedeathsMLP$Selective_deaths_births)
total_selective_deaths_births_MHO <- sum(birthsselectivedeathsMHO$Selective_deaths_births)
total_selective_deaths_births_MHP <- sum(birthsselectivedeathsMHP$Selective_deaths_births)
total_selective_deaths_births_TLO <- sum(birthsselectivedeathsTLO$Selective_deaths_births)
total_selective_deaths_births_TLP <- sum(birthsselectivedeathsTLP$Selective_deaths_births)
total_selective_deaths_births_THO <- sum(birthsselectivedeathsTHO$Selective_deaths_births)
total_selective_deaths_births_THP <- sum(birthsselectivedeathsTHP$Selective_deaths_births)
total_selective_deaths_births_MHP

total_seeds_for_next_generation_MLO <- sum(birthsselectivedeathsMLO$Seeds_for_next_generation)
total_seeds_for_next_generation_MLP <- sum(birthsselectivedeathsMLP$Seeds_for_next_generation)
total_seeds_for_next_generation_MHO <- sum(birthsselectivedeathsMHO$Seeds_for_next_generation)
total_seeds_for_next_generation_MHP <- sum(birthsselectivedeathsMHP$Seeds_for_next_generation)
total_seeds_for_next_generation_TLO <- sum(birthsselectivedeathsTLO$Seeds_for_next_generation)
total_seeds_for_next_generation_TLP <- sum(birthsselectivedeathsTLP$Seeds_for_next_generation)
total_seeds_for_next_generation_THO <- sum(birthsselectivedeathsTHO$Seeds_for_next_generation)
total_seeds_for_next_generation_THP <- sum(birthsselectivedeathsTHP$Seeds_for_next_generation)
total_seeds_for_next_generation_MHP

proportion_deaths_selective_births_MLO <- total_selective_deaths_births_MLO / (sum(birthsselectivedeathsMLO$Surviving_individuals_pot) * max(birthsselectivedeathsMLO$Seeds_by_ind) - total_seeds_for_next_generation_MLO)
proportion_deaths_selective_births_MLP <- total_selective_deaths_births_MLP / (sum(birthsselectivedeathsMLP$Surviving_individuals_pot) * max(birthsselectivedeathsMLP$Seeds_by_ind) - total_seeds_for_next_generation_MLP)
proportion_deaths_selective_births_MHO <- total_selective_deaths_births_MHO / (sum(birthsselectivedeathsMHO$Surviving_individuals_pot) * max(birthsselectivedeathsMHO$Seeds_by_ind) - total_seeds_for_next_generation_MHO)
proportion_deaths_selective_births_MHP <- total_selective_deaths_births_MHP / (sum(birthsselectivedeathsMHP$Surviving_individuals_pot) * max(birthsselectivedeathsMHP$Seeds_by_ind) - total_seeds_for_next_generation_MHP)
proportion_deaths_selective_births_TLO <- total_selective_deaths_births_TLO / (sum(birthsselectivedeathsTLO$Surviving_individuals_pot) * max(birthsselectivedeathsTLO$Seeds_by_ind) - total_seeds_for_next_generation_TLO)
proportion_deaths_selective_births_TLP <- total_selective_deaths_births_TLP / (sum(birthsselectivedeathsTLP$Surviving_individuals_pot) * max(birthsselectivedeathsTLP$Seeds_by_ind) - total_seeds_for_next_generation_TLP)
proportion_deaths_selective_births_THO <- total_selective_deaths_births_THO / (sum(birthsselectivedeathsTHO$Surviving_individuals_pot) * max(birthsselectivedeathsTHO$Seeds_by_ind) - total_seeds_for_next_generation_THO)
proportion_deaths_selective_births_THP <- total_selective_deaths_births_THP / (sum(birthsselectivedeathsTHP$Surviving_individuals_pot) * max(birthsselectivedeathsTHP$Seeds_by_ind) - total_seeds_for_next_generation_THP)

proportion_deaths_selective_births_MLO
proportion_deaths_selective_births_MLP
proportion_deaths_selective_births_MHO

proportion_deaths_selective_deaths_MLO <- total_selective_deaths_survival_MLO / (sum(survivalselectivedeathsMLO$Replicates))
proportion_deaths_selective_deaths_MLP <- total_selective_deaths_survival_MLP / (sum(survivalselectivedeathsMLP$Replicates)*30)
proportion_deaths_selective_deaths_MHO <- total_selective_deaths_survival_MHO / (sum(survivalselectivedeathsMHO$Replicates))
proportion_deaths_selective_deaths_MHP <- total_selective_deaths_survival_MHP / (sum(survivalselectivedeathsMHP$Replicates)*30)
proportion_deaths_selective_deaths_TLO <- total_selective_deaths_survival_TLO / (sum(survivalselectivedeathsTLO$Replicates))
proportion_deaths_selective_deaths_TLP <- total_selective_deaths_survival_TLP / (sum(survivalselectivedeathsTLP$Replicates)*30)
proportion_deaths_selective_deaths_THO <- total_selective_deaths_survival_THO / (sum(survivalselectivedeathsTHO$Replicates))
proportion_deaths_selective_deaths_THP <- total_selective_deaths_survival_THP / (sum(survivalselectivedeathsTHP$Replicates)*30)
proportion_deaths_selective_deaths_MLO
proportion_deaths_selective_deaths_MLP

total_surviving_seeds_MLO <- sum(survivalselectivedeathsMLO$Total_surviving_seeds)
total_surviving_seeds_MLP <- sum(survivalselectivedeathsMLP$Total_surviving_seeds)
total_surviving_seeds_MHO <- sum(survivalselectivedeathsMHO$Total_surviving_seeds)
total_surviving_seeds_MHP <- sum(survivalselectivedeathsMHP$Total_surviving_seeds)
total_surviving_seeds_TLO <- sum(survivalselectivedeathsTLO$Total_surviving_seeds)
total_surviving_seeds_TLP <- sum(survivalselectivedeathsTLP$Total_surviving_seeds)
total_surviving_seeds_THO <- sum(survivalselectivedeathsTHO$Total_surviving_seeds)
total_surviving_seeds_THP <- sum(survivalselectivedeathsTHP$Total_surviving_seeds)

total_seeds_produced_MLO <- sum(birthsselectivedeathsMLO$Seeds_total_pot)
total_seeds_produced_MLP <- sum(birthsselectivedeathsMLP$Seeds_total_pot)
total_seeds_produced_MHO <- sum(birthsselectivedeathsMHO$Seeds_total_pot)
total_seeds_produced_MHP <- sum(birthsselectivedeathsMHP$Seeds_total_pot)
total_seeds_produced_TLO <- sum(birthsselectivedeathsTLO$Seeds_total_pot)
total_seeds_produced_TLP <- sum(birthsselectivedeathsTLP$Seeds_total_pot)
total_seeds_produced_THO <- sum(birthsselectivedeathsTHO$Seeds_total_pot)
total_seeds_produced_THP <- sum(birthsselectivedeathsTHP$Seeds_total_pot)

kfecundityMLO <- (total_seeds_produced_MLO / total_surviving_seeds_MLO)
kfecundityMLP <- (total_seeds_produced_MLP / total_surviving_seeds_MLP)
kfecundityMHO <- (total_seeds_produced_MHO / total_surviving_seeds_MHO)
kfecundityMHP <- (total_seeds_produced_MHP / total_surviving_seeds_MHP)
kfecundityTLO <- (total_seeds_produced_TLO / total_surviving_seeds_TLO)
kfecundityTLP <- (total_seeds_produced_TLP / total_surviving_seeds_TLP)
kfecundityTHO <- (total_seeds_produced_THO / total_surviving_seeds_THO)
kfecundityTHP <- (total_seeds_produced_THP / total_surviving_seeds_THP)


kfecundityMLO
kfecundityMLP
kfecundityMHO
kfecundityMHP
kfecundityTLO
kfecundityTLP
kfecundityTHO
kfecundityTHP

kbestsurvivalMLO <- 1 - (min(survivalselectivedeathsMLO$Genotype_death_rate))
kbestsurvivalMHO <- 1 - (min(survivalselectivedeathsMHO$Genotype_death_rate))
kbestsurvivalTLO <- 1 - (min(survivalselectivedeathsTLO$Genotype_death_rate))
kbestsurvivalTHO <- 1 - (min(survivalselectivedeathsTHO$Genotype_death_rate))
kbestsurvivalMLP <- 1 - (min(survivalselectivedeathsMLP$Genotype_death_rate))
kbestsurvivalMHP <- 1 - (min(survivalselectivedeathsMHP$Genotype_death_rate))
kbestsurvivalTLP <- 1 - (min(survivalselectivedeathsTLP$Genotype_death_rate))
kbestsurvivalTHP <- 1 - (min(survivalselectivedeathsTHP$Genotype_death_rate))

kmeansurvivalMLO <- total_surviving_seeds_MLO / (sum(survivalselectivedeathsMLO$Replicates))
kmeansurvivalMHO <- total_surviving_seeds_MHO / (sum(survivalselectivedeathsMHO$Replicates))
kmeansurvivalTLO <- total_surviving_seeds_TLO / (sum(survivalselectivedeathsTLO$Replicates))
kmeansurvivalTHO <- total_surviving_seeds_THO / (sum(survivalselectivedeathsTHO$Replicates))
kmeansurvivalMLP <- total_surviving_seeds_MLP / (30*sum(survivalselectivedeathsMLP$Replicates))
kmeansurvivalMHP <- total_surviving_seeds_MHP / (30*sum(survivalselectivedeathsMHP$Replicates))
kmeansurvivalTLP <- total_surviving_seeds_TLP / (30*sum(survivalselectivedeathsTLP$Replicates))
kmeansurvivalTHP <- total_surviving_seeds_THP / (30*sum(survivalselectivedeathsTHP$Replicates))

reproductive_excess_fecundity_lowdispersal_MLO <- total_surviving_seeds_MLO * (kfecundityMLO - (1 / (0.001 * kbestsurvivalMLO)))
reproductive_excess_fecundity_lowdispersal_MHO <- total_surviving_seeds_MHO * (kfecundityMHO - (1 / (0.001 * kbestsurvivalMHO)))
reproductive_excess_fecundity_lowdispersal_TLO <- total_surviving_seeds_TLO * (kfecundityTLO - (1 / (0.001 * kbestsurvivalTLO)))
reproductive_excess_fecundity_lowdispersal_THO <- total_surviving_seeds_THO * (kfecundityTHO - (1 / (0.001 * kbestsurvivalTHO)))
reproductive_excess_fecundity_lowdispersal_MLP <- total_surviving_seeds_MLP * (kfecundityMLP - (1 / (0.01 * kbestsurvivalMLP)))
reproductive_excess_fecundity_lowdispersal_MHP <- total_surviving_seeds_MHP * (kfecundityMHP - (1 / (0.01 * kbestsurvivalMHP)))
reproductive_excess_fecundity_lowdispersal_TLP <- total_surviving_seeds_TLP * (kfecundityTLP - (1 / (0.01 * kbestsurvivalTLP)))
reproductive_excess_fecundity_lowdispersal_THP <- total_surviving_seeds_THP * (kfecundityTHP - (1 / (0.01 * kbestsurvivalTHP)))

reproductive_excess_fecundity_highdispersal_MLO <- total_surviving_seeds_MLO * (kfecundityMLO - (1 / (0.01 * kbestsurvivalMLO)))
reproductive_excess_fecundity_highdispersal_MHO <- total_surviving_seeds_MHO * (kfecundityMHO - (1 / (0.01 * kbestsurvivalMHO)))
reproductive_excess_fecundity_highdispersal_TLO <- total_surviving_seeds_TLO * (kfecundityTLO - (1 / (0.01 * kbestsurvivalTLO)))
reproductive_excess_fecundity_highdispersal_THO <- total_surviving_seeds_THO * (kfecundityTHO - (1 / (0.01 * kbestsurvivalTHO)))
reproductive_excess_fecundity_highdispersal_MLP <- total_surviving_seeds_MLP * (kfecundityMLP - (1 / (0.1 * kbestsurvivalMLP)))
reproductive_excess_fecundity_highdispersal_MHP <- total_surviving_seeds_MHP * (kfecundityMHP - (1 / (0.1 * kbestsurvivalMHP)))
reproductive_excess_fecundity_highdispersal_TLP <- total_surviving_seeds_TLP * (kfecundityTLP - (1 / (0.1 * kbestsurvivalTLP)))
reproductive_excess_fecundity_highdispersal_THP <- total_surviving_seeds_THP * (kfecundityTHP - (1 / (0.1 * kbestsurvivalTHP)))

reproductive_excess_survival_lowdispersal_MLO <- sum(survivalselectivedeathsMLO$Replicates) * (kbestsurvivalMLO - (1 / (0.001 * kfecundityMLO)))
reproductive_excess_survival_lowdispersal_MHO <- sum(survivalselectivedeathsMHO$Replicates) * (kbestsurvivalMHO - (1 / (0.001 * kfecundityMHO)))
reproductive_excess_survival_lowdispersal_TLO <- sum(survivalselectivedeathsTLO$Replicates) * (kbestsurvivalTLO - (1 / (0.001 * kfecundityTLO)))
reproductive_excess_survival_lowdispersal_THO <- sum(survivalselectivedeathsTHO$Replicates) * (kbestsurvivalTHO - (1 / (0.001 * kfecundityTHO)))
reproductive_excess_survival_lowdispersal_MLP <- (30*sum(survivalselectivedeathsMLP$Replicates)) * (kbestsurvivalMLP - (1 / (0.01 * kfecundityMLP)))
reproductive_excess_survival_lowdispersal_MHP <- (30*sum(survivalselectivedeathsMHP$Replicates)) * (kbestsurvivalMHP - (1 / (0.01 * kfecundityMHP)))
reproductive_excess_survival_lowdispersal_TLP <- (30*sum(survivalselectivedeathsTLP$Replicates)) * (kbestsurvivalTLP - (1 / (0.01 * kfecundityTLP)))
reproductive_excess_survival_lowdispersal_THP <- (30*sum(survivalselectivedeathsTHP$Replicates)) * (kbestsurvivalTHP - (1 / (0.01 * kfecundityTHP)))

reproductive_excess_survival_highdispersal_MLO <- sum(survivalselectivedeathsMLO$Replicates) * (kbestsurvivalMLO - (1 / (0.01 * kfecundityMLO)))
reproductive_excess_survival_highdispersal_MHO <- sum(survivalselectivedeathsMHO$Replicates) * (kbestsurvivalMHO - (1 / (0.01 * kfecundityMHO)))
reproductive_excess_survival_highdispersal_TLO <- sum(survivalselectivedeathsTLO$Replicates) * (kbestsurvivalTLO - (1 / (0.01 * kfecundityTLO)))
reproductive_excess_survival_highdispersal_THO <- sum(survivalselectivedeathsTHO$Replicates) * (kbestsurvivalTHO - (1 / (0.01 * kfecundityTHO)))
reproductive_excess_survival_highdispersal_MLP <- (30*sum(survivalselectivedeathsMLP$Replicates)) * (kbestsurvivalMLP - (1 / (0.1 * kfecundityMLP)))
reproductive_excess_survival_highdispersal_MHP <- (30*sum(survivalselectivedeathsMHP$Replicates)) * (kbestsurvivalMHP - (1 / (0.1 * kfecundityMHP)))
reproductive_excess_survival_highdispersal_TLP <- (30*sum(survivalselectivedeathsTLP$Replicates)) * (kbestsurvivalTLP - (1 / (0.1 * kfecundityTLP)))
reproductive_excess_survival_highdispersal_THP <- (30*sum(survivalselectivedeathsTHP$Replicates)) * (kbestsurvivalTHP - (1 / (0.1 * kfecundityTHP)))

reproductive_excess_fecundity_lowdispersal_MLO / sum(survivalselectivedeathsMLO$Replicates)
reproductive_excess_fecundity_lowdispersal_MHO / sum(survivalselectivedeathsMHO$Replicates)
reproductive_excess_fecundity_lowdispersal_TLO / sum(survivalselectivedeathsTLO$Replicates)
reproductive_excess_fecundity_lowdispersal_THO / sum(survivalselectivedeathsTHO$Replicates)
reproductive_excess_fecundity_lowdispersal_MLP / sum(survivalselectivedeathsMLP$Replicates)
reproductive_excess_fecundity_lowdispersal_MHP / sum(survivalselectivedeathsMHP$Replicates)
reproductive_excess_fecundity_lowdispersal_TLP / sum(survivalselectivedeathsTLP$Replicates)
reproductive_excess_fecundity_lowdispersal_THP / sum(survivalselectivedeathsTHP$Replicates)

reproductive_excess_survival_lowdispersal_MLO
reproductive_excess_survival_lowdispersal_MHO
reproductive_excess_survival_lowdispersal_TLO
reproductive_excess_survival_lowdispersal_THO
reproductive_excess_survival_lowdispersal_MLP
reproductive_excess_survival_lowdispersal_MHP
reproductive_excess_survival_lowdispersal_TLP
reproductive_excess_survival_lowdispersal_THP

reproductive_excess_survival_lowdispersal_MLO / sum(survivalselectivedeathsMLO$Replicates)
reproductive_excess_survival_lowdispersal_MHO / sum(survivalselectivedeathsMHO$Replicates)
reproductive_excess_survival_lowdispersal_TLO / sum(survivalselectivedeathsTLO$Replicates)
reproductive_excess_survival_lowdispersal_THO / sum(survivalselectivedeathsTHO$Replicates)
reproductive_excess_survival_lowdispersal_MLP / sum(survivalselectivedeathsMLP$Replicates)
reproductive_excess_survival_lowdispersal_MHP / sum(survivalselectivedeathsMHP$Replicates)
reproductive_excess_survival_lowdispersal_TLP / sum(survivalselectivedeathsTLP$Replicates)
reproductive_excess_survival_lowdispersal_THP / sum(survivalselectivedeathsTHP$Replicates)

reproductive_excess_fecundity_highdispersal_MLO / sum(survivalselectivedeathsMLO$Replicates)
reproductive_excess_fecundity_highdispersal_MHO / sum(survivalselectivedeathsMHO$Replicates)
reproductive_excess_fecundity_highdispersal_TLO / sum(survivalselectivedeathsTLO$Replicates)
reproductive_excess_fecundity_highdispersal_THO / sum(survivalselectivedeathsTHO$Replicates)
reproductive_excess_fecundity_highdispersal_MLP / sum(survivalselectivedeathsMLP$Replicates)
reproductive_excess_fecundity_highdispersal_MHP / sum(survivalselectivedeathsMHP$Replicates)
reproductive_excess_fecundity_highdispersal_TLP / sum(survivalselectivedeathsTLP$Replicates)
reproductive_excess_fecundity_highdispersal_THP / sum(survivalselectivedeathsTHP$Replicates)

reproductive_excess_survival_highdispersal_MLO
reproductive_excess_survival_highdispersal_MHO
reproductive_excess_survival_highdispersal_TLO
reproductive_excess_survival_highdispersal_THO
reproductive_excess_survival_highdispersal_MLP
reproductive_excess_survival_highdispersal_MHP
reproductive_excess_survival_highdispersal_TLP
reproductive_excess_survival_highdispersal_THP

reproductive_excess_survival_highdispersal_MLO / sum(survivalselectivedeathsMLO$Replicates)
reproductive_excess_survival_highdispersal_MHO / sum(survivalselectivedeathsMHO$Replicates)
reproductive_excess_survival_highdispersal_TLO / sum(survivalselectivedeathsTLO$Replicates)
reproductive_excess_survival_highdispersal_THO / sum(survivalselectivedeathsTHO$Replicates)
reproductive_excess_survival_highdispersal_MLP / sum(survivalselectivedeathsMLP$Replicates)
reproductive_excess_survival_highdispersal_MHP / sum(survivalselectivedeathsMHP$Replicates)
reproductive_excess_survival_highdispersal_TLP / sum(survivalselectivedeathsTLP$Replicates)
reproductive_excess_survival_highdispersal_THP / sum(survivalselectivedeathsTHP$Replicates)



diffwithandwithoutsinglepotgenotypesMLO <- sum(birthsselectivedeathsMLO$Seeds_total_pot)/sum(birthsselectivedeathsMLO$Surviving_individuals_pot) - mean(aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Mean)
diffwithandwithoutsinglepotgenotypesMLP <- sum(birthsselectivedeathsMLP$Seeds_total_pot)/sum(birthsselectivedeathsMLP$Surviving_individuals_pot) - mean(aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Mean)
diffwithandwithoutsinglepotgenotypesMHO <- sum(birthsselectivedeathsMHO$Seeds_total_pot)/sum(birthsselectivedeathsMHO$Surviving_individuals_pot) - mean(aggregatedonlysurvivingseedsnozeroesinvarianceMHO$Mean)
diffwithandwithoutsinglepotgenotypesMHP <- sum(birthsselectivedeathsMHP$Seeds_total_pot)/sum(birthsselectivedeathsMHP$Surviving_individuals_pot) - mean(aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Mean)
diffwithandwithoutsinglepotgenotypesTLO <- sum(birthsselectivedeathsTLO$Seeds_total_pot)/sum(birthsselectivedeathsTLO$Surviving_individuals_pot) - mean(aggregatedonlysurvivingseedsnozeroesinvarianceTLO$Mean)
diffwithandwithoutsinglepotgenotypesTLP <- sum(birthsselectivedeathsTLP$Seeds_total_pot)/sum(birthsselectivedeathsTLP$Surviving_individuals_pot) - mean(aggregatedonlysurvivingseedsnozeroesinvarianceTLP$Mean)
diffwithandwithoutsinglepotgenotypesTHO <- sum(birthsselectivedeathsTHO$Seeds_total_pot)/sum(birthsselectivedeathsTHO$Surviving_individuals_pot) - mean(aggregatedonlysurvivingseedsnozeroesinvarianceTHO$Mean)
diffwithandwithoutsinglepotgenotypesTHP <- sum(birthsselectivedeathsTHP$Seeds_total_pot)/sum(birthsselectivedeathsTHP$Surviving_individuals_pot) - mean(aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Mean)


diffwithandwithoutsinglepotgenotypesMLO
diffwithandwithoutsinglepotgenotypesMLP
diffwithandwithoutsinglepotgenotypesMHO
diffwithandwithoutsinglepotgenotypesMHP
diffwithandwithoutsinglepotgenotypesTLO
diffwithandwithoutsinglepotgenotypesTLP
diffwithandwithoutsinglepotgenotypesTHO
diffwithandwithoutsinglepotgenotypesTHP

diffwithandwithoutsinglepotgenotypesMLO/(sum(birthsselectivedeathsMLO$Seeds_total_pot)/sum(birthsselectivedeathsMLO$Surviving_individuals_pot))
diffwithandwithoutsinglepotgenotypesMLP/(sum(birthsselectivedeathsMLP$Seeds_total_pot)/sum(birthsselectivedeathsMLP$Surviving_individuals_pot))
diffwithandwithoutsinglepotgenotypesMHO/(sum(birthsselectivedeathsMHO$Seeds_total_pot)/sum(birthsselectivedeathsMHO$Surviving_individuals_pot))
diffwithandwithoutsinglepotgenotypesMHP/(sum(birthsselectivedeathsMHP$Seeds_total_pot)/sum(birthsselectivedeathsMHP$Surviving_individuals_pot))
diffwithandwithoutsinglepotgenotypesTLO/(sum(birthsselectivedeathsTLO$Seeds_total_pot)/sum(birthsselectivedeathsTLO$Surviving_individuals_pot))
diffwithandwithoutsinglepotgenotypesTLP/(sum(birthsselectivedeathsTLP$Seeds_total_pot)/sum(birthsselectivedeathsTLP$Surviving_individuals_pot))
diffwithandwithoutsinglepotgenotypesTHO/(sum(birthsselectivedeathsTHO$Seeds_total_pot)/sum(birthsselectivedeathsTHO$Surviving_individuals_pot))
diffwithandwithoutsinglepotgenotypesTHP/(sum(birthsselectivedeathsTHP$Seeds_total_pot)/sum(birthsselectivedeathsTHP$Surviving_individuals_pot))

#The following analysis performs a Box-cox transformation on the data and then does an ANOVA on the resulting transform
#The intent is to investigate whether there is any evidence for differences in mean fecundity between genotypes
#or if the variance between replicates is so large that it swamps any signal.
bcMLO <- boxcox(lm(onlypotswithsurvivingseedsMLO$Seeds_by_ind ~ 1))
bcMHO <- boxcox(lm(onlypotswithsurvivingseedsMHO$Seeds_by_ind ~ 1))
bcMLP <- boxcox(lm(onlypotswithsurvivingseedsMLP$Seeds_by_ind ~ 1))
bcMHP <- boxcox(lm(onlypotswithsurvivingseedsMHP$Seeds_by_ind ~ 1))
bcTLO <- boxcox(lm(onlypotswithsurvivingseedsTLO$Seeds_by_ind ~ 1))
bcTHO <- boxcox(lm(onlypotswithsurvivingseedsTHO$Seeds_by_ind ~ 1))
bcTLP <- boxcox(lm(onlypotswithsurvivingseedsTLP$Seeds_by_ind ~ 1))
bcTHP <- boxcox(lm(onlypotswithsurvivingseedsTHP$Seeds_by_ind ~ 1))

bcMLO

boxcox(lm(onlypotswithsurvivingseedsMLO$Seeds_by_ind ~ 1))
boxcox(lm(onlypotswithsurvivingseedsMHO$Seeds_by_ind ~ 1))
boxcox(lm(onlypotswithsurvivingseedsMLP$Seeds_by_ind ~ 1))
boxcox(lm(onlypotswithsurvivingseedsMHP$Seeds_by_ind ~ 1))
boxcox(lm(onlypotswithsurvivingseedsTLO$Seeds_by_ind ~ 1))
boxcox(lm(onlypotswithsurvivingseedsTHO$Seeds_by_ind ~ 1))
boxcox(lm(onlypotswithsurvivingseedsTLP$Seeds_by_ind ~ 1))
boxcox(lm(onlypotswithsurvivingseedsTHP$Seeds_by_ind ~ 1))

lambdaMLO <- bcMLO$x[which.max(bcMLO$y)]
lambdaMHO <- bcMLO$x[which.max(bcMHO$y)]
lambdaMLP <- bcMLO$x[which.max(bcMLP$y)]
lambdaMHP <- bcMLO$x[which.max(bcMHP$y)]
lambdaTLO <- bcMLO$x[which.max(bcTLO$y)]
lambdaTHO <- bcMLO$x[which.max(bcTHO$y)]
lambdaTLP <- bcMLO$x[which.max(bcTLP$y)]
lambdaTHP <- bcMLO$x[which.max(bcTHP$y)]

boxcoxtransformonlypotswithsurvivingseedsMLO <- Boxcoxtransform(onlypotswithsurvivingseedsMLO, lambdaMLO)
boxcoxtransformonlypotswithsurvivingseedsMHO <- Boxcoxtransform(onlypotswithsurvivingseedsMHO, lambdaMHO)
boxcoxtransformonlypotswithsurvivingseedsMLP <- Boxcoxtransform(onlypotswithsurvivingseedsMLP, lambdaMLP)
boxcoxtransformonlypotswithsurvivingseedsMHP <- Boxcoxtransform(onlypotswithsurvivingseedsMHP, lambdaMHP)
boxcoxtransformonlypotswithsurvivingseedsTLO <- Boxcoxtransform(onlypotswithsurvivingseedsTLO, lambdaTLO)
boxcoxtransformonlypotswithsurvivingseedsTHO <- Boxcoxtransform(onlypotswithsurvivingseedsTHO, lambdaTHO)
boxcoxtransformonlypotswithsurvivingseedsTLP <- Boxcoxtransform(onlypotswithsurvivingseedsTLP, lambdaTLP)
boxcoxtransformonlypotswithsurvivingseedsTHP <- Boxcoxtransform(onlypotswithsurvivingseedsTHP, lambdaTHP)

anovaMLO <- aov(boxcox_transform_fecundity ~ Genotype_id, data = boxcoxtransformonlypotswithsurvivingseedsMLO)
anovaMHO <- aov(boxcox_transform_fecundity ~ Genotype_id, data = boxcoxtransformonlypotswithsurvivingseedsMHO)
anovaMLP <- aov(boxcox_transform_fecundity ~ Genotype_id, data = boxcoxtransformonlypotswithsurvivingseedsMLP)
anovaMHP <- aov(boxcox_transform_fecundity ~ Genotype_id, data = boxcoxtransformonlypotswithsurvivingseedsMHP)
anovaTLO <- aov(boxcox_transform_fecundity ~ Genotype_id, data = boxcoxtransformonlypotswithsurvivingseedsTLO)
anovaTHO <- aov(boxcox_transform_fecundity ~ Genotype_id, data = boxcoxtransformonlypotswithsurvivingseedsTHO)
anovaTLP <- aov(boxcox_transform_fecundity ~ Genotype_id, data = boxcoxtransformonlypotswithsurvivingseedsTLP)
anovaTHP <- aov(boxcox_transform_fecundity ~ Genotype_id, data = boxcoxtransformonlypotswithsurvivingseedsTHP)

summary(anovaMLO)
summary(anovaMHO)
summary(anovaMLP)
summary(anovaMHP)
summary(anovaTLO)
summary(anovaTHO)
summary(anovaTLP)
summary(anovaTHP)

hist(boxcoxtransformonlypotswithsurvivingseedsMLO$boxcox_transform_fecundity)

qqnorm(boxcoxtransformonlypotswithsurvivingseedsMLO$boxcox_transform_fecundity)
qqnorm(boxcoxtransformonlypotswithsurvivingseedsMHO$boxcox_transform_fecundity)
qqnorm(boxcoxtransformonlypotswithsurvivingseedsTLO$boxcox_transform_fecundity)
qqnorm(boxcoxtransformonlypotswithsurvivingseedsTHO$boxcox_transform_fecundity)
qqnorm(boxcoxtransformonlypotswithsurvivingseedsMLP$boxcox_transform_fecundity)
qqnorm(boxcoxtransformonlypotswithsurvivingseedsMHP$boxcox_transform_fecundity)
qqnorm(boxcoxtransformonlypotswithsurvivingseedsTLP$boxcox_transform_fecundity)
qqnorm(boxcoxtransformonlypotswithsurvivingseedsTHP$boxcox_transform_fecundity)

ggplot(data = boxcoxtransformonlypotswithsurvivingseedsMLO, mapping = aes(x = boxcox_transform_fecundity)) +
  xlab("Transformed fecundity") +
  geom_histogram() +
  theme_classic(base_size = 18)

ggplot(data = boxcoxtransformonlypotswithsurvivingseedsMHO, mapping = aes(x = boxcox_transform_fecundity)) +
  xlab("Transformed fecundity") +
  geom_histogram() +
  theme_classic(base_size = 18)

ggplot(data = boxcoxtransformonlypotswithsurvivingseedsTLO, mapping = aes(x = boxcox_transform_fecundity)) +
  xlab("Transformed fecundity") +
  geom_histogram() +
  theme_classic(base_size = 18)

ggplot(data = boxcoxtransformonlypotswithsurvivingseedsTHO, mapping = aes(x = boxcox_transform_fecundity)) +
  xlab("Transformed fecundity") +
  geom_histogram() +
  theme_classic(base_size = 18)

ggplot(data = boxcoxtransformonlypotswithsurvivingseedsMLP, mapping = aes(x = boxcox_transform_fecundity)) +
  xlab("Transformed fecundity") +
  geom_histogram() +
  theme_classic(base_size = 18)

ggplot(data = boxcoxtransformonlypotswithsurvivingseedsMHP, mapping = aes(x = boxcox_transform_fecundity)) +
  xlab("Transformed fecundity") +
  geom_histogram() +
  theme_classic(base_size = 18)

ggplot(data = boxcoxtransformonlypotswithsurvivingseedsTLP, mapping = aes(x = boxcox_transform_fecundity)) +
  xlab("Transformed fecundity") +
  geom_histogram() +
  theme_classic(base_size = 18)

ggplot(data = boxcoxtransformonlypotswithsurvivingseedsTHP, mapping = aes(x = boxcox_transform_fecundity)) +
  xlab("Transformed fecundity") +
  geom_histogram() +
  theme_classic(base_size = 18)

#Following data frames are to perform a statistical test for significance
#on the relationship between reproductive excess and proportion of deaths selective
#across environmental conditions.
#It will also be used to create a scatter plot of this relationship.
scatterplotlowdensitydeathsdataframe <- data.frame(
  environmental_condition = c("MLI", "MHI", "TLI", "THI"),
  reproductive_excess_fecundity_low = c(reproductive_excess_fecundity_lowdispersal_MLO, reproductive_excess_fecundity_lowdispersal_MHO, reproductive_excess_fecundity_lowdispersal_TLO, reproductive_excess_fecundity_lowdispersal_THO),
  reproductive_excess_fecundity_high = c(reproductive_excess_fecundity_highdispersal_MLO, reproductive_excess_fecundity_highdispersal_MHO, reproductive_excess_fecundity_highdispersal_TLO, reproductive_excess_fecundity_highdispersal_THO),
  reproductive_excess_survival_low = c(reproductive_excess_survival_lowdispersal_MLO, reproductive_excess_survival_lowdispersal_MHO, reproductive_excess_survival_lowdispersal_TLO, reproductive_excess_survival_lowdispersal_THO),
  reproductive_excess_survival_high = c(reproductive_excess_survival_highdispersal_MLO, reproductive_excess_survival_highdispersal_MHO, reproductive_excess_survival_highdispersal_TLO, reproductive_excess_survival_highdispersal_THO),
  proportion_of_deaths_selective = c(proportion_deaths_selective_deaths_MLO, proportion_deaths_selective_deaths_MHO, proportion_deaths_selective_deaths_TLO, proportion_deaths_selective_deaths_THO)
)

scatterplothighdensitydeathsdataframe <- data.frame(
  environmental_condition = c("MLP", "MHP", "TLP", "THP"),
  reproductive_excess_fecundity_low = c(reproductive_excess_fecundity_lowdispersal_MLP, reproductive_excess_fecundity_lowdispersal_MHP, reproductive_excess_fecundity_lowdispersal_TLP, reproductive_excess_fecundity_lowdispersal_THP),
  reproductive_excess_fecundity_high = c(reproductive_excess_fecundity_highdispersal_MLP, reproductive_excess_fecundity_highdispersal_MHP, reproductive_excess_fecundity_highdispersal_TLP, reproductive_excess_fecundity_highdispersal_THP),
  reproductive_excess_survival_low = c(reproductive_excess_survival_lowdispersal_MLP, reproductive_excess_survival_lowdispersal_MHP, reproductive_excess_survival_lowdispersal_TLP, reproductive_excess_survival_lowdispersal_THP),
  reproductive_excess_survival_high = c(reproductive_excess_survival_highdispersal_MLP, reproductive_excess_survival_highdispersal_MHP, reproductive_excess_survival_highdispersal_TLP, reproductive_excess_survival_highdispersal_THP),
  proportion_of_deaths_selective = c(proportion_deaths_selective_deaths_MLP, proportion_deaths_selective_deaths_MHP, proportion_deaths_selective_deaths_TLP, proportion_deaths_selective_deaths_THP)
)

cor.test(scatterplothighdensitydeathsdataframe$reproductive_excess, scatterplothighdensitydeathsdataframe$proportion_of_deaths_selective, method = "spearman")
cor.test(scatterplotlowdensitydeathsdataframe$reproductive_excess, scatterplotlowdensitydeathsdataframe$proportion_of_deaths_selective, method = "spearman")
pvaluesexcessvspropdeathsselective <- c(0.4167, 0.08333)
pchisq(-2*sum(log(pvaluesexcessvspropdeathsselective)), 6, lower.tail = FALSE)


#PLOTS

#Following plots show standard deviation and mean in seed production
#OMITTED POTS WITH NO SURVIVING SEEDS AND GENOTYPES WITH NO VARIANCE
stdevvsmeanplotMLO <- ggplot(data = onlysurvivingseedsnozeroesinvarianceMLO, mapping = aes(x = Mean, y = StDev)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_classic(base_size = 15)

stdevvsmeanplotMLP <- ggplot(data = onlysurvivingseedsnozeroesinvarianceMLP, mapping = aes(x = Mean, y = StDev)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_classic(base_size = 15)

stdevvsmeanplotMHO <- ggplot(data = onlysurvivingseedsnozeroesinvarianceMHO, mapping = aes(x = Mean, y = StDev)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_classic(base_size = 15)

stdevvsmeanplotMHP <- ggplot(data = onlysurvivingseedsnozeroesinvarianceMHP, mapping = aes(x = Mean, y = StDev)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_classic(base_size = 15)

stdevvsmeanplotTLO <- ggplot(data = onlysurvivingseedsnozeroesinvarianceTLO, mapping = aes(x = Mean, y = StDev)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_classic(base_size = 15)

stdevvsmeanplotTLP <- ggplot(data = onlysurvivingseedsnozeroesinvarianceTLP, mapping = aes(x = Mean, y = StDev)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_classic(base_size = 15)

stdevvsmeanplotTHO <- ggplot(data = onlysurvivingseedsnozeroesinvarianceTHO, mapping = aes(x = Mean, y = StDev)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_classic(base_size = 15)

stdevvsmeanplotTHP <- ggplot(data = onlysurvivingseedsnozeroesinvarianceTHP, mapping = aes(x = Mean, y = StDev)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_classic(base_size = 15)

stdevvsmeanplotMLO
stdevvsmeanplotMLP
stdevvsmeanplotMHO
stdevvsmeanplotMHP
stdevvsmeanplotTLO
stdevvsmeanplotTLP
stdevvsmeanplotTHO
stdevvsmeanplotTHP

genotypemeanseedproductionhistMLO <- ggplot(data = aggregatedonlysurvivingseedsnozeroesinvarianceMLO, aes(x = Seeds_by_ind)) +
  geom_histogram() +
  xlab("Genotype mean seed production") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 12000, y = 40, label = "MLI", size = 8)

genotypemeanseedproductionhistMLP <- ggplot(data = aggregatedonlysurvivingseedsnozeroesinvarianceMLP, aes(x = Seeds_by_ind)) +
  geom_histogram() +
  xlab("Genotype mean seed production") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 4000, y = 4, label = "MLP", size = 8)

genotypemeanseedproductionhistMHO <- ggplot(data = aggregatedonlysurvivingseedsnozeroesinvarianceMHO, aes(x = Seeds_by_ind)) +
  geom_histogram() +
  xlab("Genotype mean seed production") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 25000, y = 40, label = "MHI", size = 8)

genotypemeanseedproductionhistMHP <- ggplot(data = aggregatedonlysurvivingseedsnozeroesinvarianceMHP, aes(x = Seeds_by_ind)) +
  geom_histogram() +
  xlab("Genotype mean seed production") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 20000, y = 40, label = "MHP", size = 8)

genotypemeanseedproductionhistTLO <- ggplot(data = aggregatedonlysurvivingseedsnozeroesinvarianceTLO, aes(x = Seeds_by_ind)) +
  geom_histogram() +
  xlab("Genotype mean seed production") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 9000, y = 40, label = "TLI", size = 8)

genotypemeanseedproductionhistTLP <- ggplot(data = aggregatedonlysurvivingseedsnozeroesinvarianceTLP, aes(x = Seeds_by_ind)) +
  geom_histogram() +
  xlab("Genotype mean seed production") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 2500, y = 40, label = "TLP", size = 8)

genotypemeanseedproductionhistTHO <- ggplot(data = aggregatedonlysurvivingseedsnozeroesinvarianceTHO, aes(x = Seeds_by_ind)) +
  geom_histogram() +
  xlab("Genotype mean seed production") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 16000, y = 40, label = "THI", size = 8)

genotypemeanseedproductionhistTHP <- ggplot(data = aggregatedonlysurvivingseedsnozeroesinvarianceTHP, aes(x = Seeds_by_ind)) +
  geom_histogram() +
  xlab("Genotype mean seed production") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 6000, y = 150, label = "THP", size = 8)

genotypemeanseedproductionhistMLO
genotypemeanseedproductionhistMLP
genotypemeanseedproductionhistMHO
genotypemeanseedproductionhistMHP
genotypemeanseedproductionhistTLO
genotypemeanseedproductionhistTLP
genotypemeanseedproductionhistTHO
genotypemeanseedproductionhistTHP

genotypemeanseedsurvivalhistMLO <- ggplot(data = survivalselectivedeathsMLO, aes(x = (1 - Genotype_death_rate))) +
  geom_histogram() +
  xlab("Genotype mean seed survival") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 0.8, y = 80, label = "MLI", size = 8)

genotypemeanseedsurvivalhistMLP <- ggplot(data = survivalselectivedeathsMLP, aes(x = (1 - Genotype_death_rate))) +
  geom_histogram() +
  xlab("Genotype mean seed survival") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 0.8, y = 300, label = "MLP", size = 8)

genotypemeanseedsurvivalhistMHO <- ggplot(data = survivalselectivedeathsMHO, aes(x = (1 - Genotype_death_rate))) +
  geom_histogram() +
  xlab("Genotype mean seed survival") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 0.8, y = 280, label = "MHI", size = 8)

genotypemeanseedsurvivalhistMHP <- ggplot(data = survivalselectivedeathsMHP, aes(x = (1 - Genotype_death_rate))) +
  geom_histogram() +
  xlab("Genotype mean seed survival") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 0.8, y = 40, label = "MHP", size = 8)

genotypemeanseedsurvivalhistTLO <- ggplot(data = survivalselectivedeathsTLO, aes(x = (1 - Genotype_death_rate))) +
  geom_histogram() +
  xlab("Genotype mean seed survival") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 0.8, y = 120, label = "TLI", size = 8)

genotypemeanseedsurvivalhistTLP <- ggplot(data = survivalselectivedeathsTLP, aes(x = (1 - Genotype_death_rate))) +
  geom_histogram() +
  xlab("Genotype mean seed survival") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 0.7, y = 100, label = "TLP", size = 8)

genotypemeanseedsurvivalhistTHO <- ggplot(data = survivalselectivedeathsTHO, aes(x = (1 - Genotype_death_rate))) +
  geom_histogram() +
  xlab("Genotype mean seed survival") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 0.7, y = 180, label = "THI", size = 8)

genotypemeanseedsurvivalhistTHP <- ggplot(data = survivalselectivedeathsTHP, aes(x = (1 - Genotype_death_rate))) +
  geom_histogram() +
  xlab("Genotype mean seed survival") +
  ylab("Count") +
  theme_classic(base_size = 15) +
  annotate("text", x = 0.8, y = 40, label = "THP", size = 8)

genotypemeanseedsurvivalhistMLO
genotypemeanseedsurvivalhistMLP
genotypemeanseedsurvivalhistMHO
genotypemeanseedsurvivalhistMHP
genotypemeanseedsurvivalhistTLO
genotypemeanseedsurvivalhistTLP
genotypemeanseedsurvivalhistTHO
genotypemeanseedsurvivalhistTHP

#Following lines are only to visualize proportion of deaths selective
#and minimum reproductive excess required in a model similar to 
#Desai and Fisher 2007
#the following is the transcendental equation for q,
#expressed as a function where the equation is true at 0
#so I can use uniroot to solve it later
qfunction = function(q, s, N, U) {q - (2*log(N*q*s)) / (log((s/U)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
qapproximatefunction = function(s, N, U) {2*log(N*s) / log(s/U)}

xvaluesfors <- c(0.0035, 0.005, 0.0075, 0.01, 0.02, 0.035, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25)
xvaluesforN <- c(100000, 200000, 350000, 500000, 750000, 1000000, 2000000, 3500000, 5000000, 7500000, 10000000)
xvaluesforU <- c(0.0001, 0.0002, 0.00035, 0.0005, 0.00075, 0.001, 0.002, 0.0025)

approximateqvaluesfors <- c(1:13)
approximateqvaluesforN <- c(1:11)
approximateqvaluesforU <- c(1:8)


for (i in 1:length(approximateqvaluesfors)) {
  approximateqvaluesfors[i] <- qapproximatefunction(xvaluesfors[i], 10^6, 0.001)
}

for (i in 1:length(approximateqvaluesforN)) {
  approximateqvaluesforN[i] <- qapproximatefunction(0.01, xvaluesforN[i], 0.001)
}

for (i in 1:length(approximateqvaluesforU)) {
  approximateqvaluesforU[i] <- qapproximatefunction(0.01, 10^6, xvaluesforU[i])
}

dataframefordesaifisherbysapproximate <- data.frame(
  xvaluesfors,
  approximateqvaluesfors
  
)

dataframefordesaifisherbyNapproximate <- data.frame(
  xvaluesforN,
  approximateqvaluesforN
  
)

dataframefordesaifisherbyUapproximate <- data.frame(
  xvaluesforU,
  approximateqvaluesforU
  
)

#It's appallingly tedious that I have to do this for every data point,
#But I can't figure out a way to pass a function with multiple parameters to uniroot.
#Here's all the roots for ten values of s, using N = 10^6 and U = 0.001
qfunctions002 = function(q) {q - (2*log(10^6*q*0.002)) / (log((0.002/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions002, interval = c(4,50))

qfunctions0035 = function(q) {q - (2*log(10^6*q*0.0035)) / (log((0.0035/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions0035, interval = c(4,50))

qfunctions005 = function(q) {q - (2*log(10^6*q*0.005)) / (log((0.005/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions005, interval = c(4,50))

qfunctions0075 = function(q) {q - (2*log(10^6*q*0.0075)) / (log((0.0075/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions0075, interval = c(4,50))

qfunctions01 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions01, interval = c(4,50))

qfunctions02 = function(q) {q - (2*log(10^6*q*0.02)) / (log((0.02/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions02, interval = c(4,50))

qfunctions035 = function(q) {q - (2*log(10^6*q*0.035)) / (log((0.035/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions035, interval = c(4,50))

qfunctions05 = function(q) {q - (2*log(10^6*q*0.05)) / (log((0.05/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions05, interval = c(4,50))

qfunctions075 = function(q) {q - (2*log(10^6*q*0.075)) / (log((0.075/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions075, interval = c(4,50))

qfunctions1 = function(q) {q - (2*log(10^6*q*0.1)) / (log((0.1/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions1, interval = c(4,50))

qfunctions125 = function(q) {q - (2*log(10^6*q*0.125)) / (log((0.125/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions125, interval = c(4,50))

qfunctions15 = function(q) {q - (2*log(10^6*q*0.15)) / (log((0.15/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions15, interval = c(4,50))

qfunctions2 = function(q) {q - (2*log(10^6*q*0.2)) / (log((0.2/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions2, interval = c(4,50))

qfunctions25 = function(q) {q - (2*log(10^6*q*0.25)) / (log((0.25/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctions25, interval = c(4,50))

#Here's all the roots for eleven different values of N, using s = 0.01 and U = 0.001
qfunctionN100000 = function(q) {q - (2*log(100000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN100000, interval = c(4,50))

qfunctionN200000 = function(q) {q - (2*log(200000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN200000, interval = c(4,50))

qfunctionN350000 = function(q) {q - (2*log(350000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN350000, interval = c(4,50))

qfunctionN500000 = function(q) {q - (2*log(500000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN500000, interval = c(4,50))

qfunctionN750000 = function(q) {q - (2*log(750000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN750000, interval = c(4,50))

qfunctionN1000000 = function(q) {q - (2*log(1000000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN1000000, interval = c(4,50))

qfunctionN2000000 = function(q) {q - (2*log(2000000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN2000000, interval = c(4,50))

qfunctionN3500000 = function(q) {q - (2*log(3500000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN3500000, interval = c(4,50))

qfunctionN5000000 = function(q) {q - (2*log(5000000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN5000000, interval = c(4,50))

qfunctionN7500000 = function(q) {q - (2*log(7500000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN7500000, interval = c(4,50))

qfunctionN10000000 = function(q) {q - (2*log(10000000*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionN10000000, interval = c(4,50))

#Here's the roots for nine different values of U, using s = 0.01 and N = 10^6 (less values because of avoiding s/u > 2)
qfunctionU0001 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.0001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionU0001, interval = c(4,50))

qfunctionU0002 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.0002)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionU0002, interval = c(4,50))

qfunctionU00035 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.00035)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionU00035, interval = c(4,50))

qfunctionU0005 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.0005)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionU0005, interval = c(4,50))

qfunctionU00075 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.00075)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionU00075, interval = c(4,50))

qfunctionU001 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.001)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionU001, interval = c(4,50))

qfunctionU002 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.002)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionU002, interval = c(4,50))

qfunctionU0025 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.0025)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionU0025, interval = c(4,50))

qfunctionU005 = function(q) {q - (2*log(10^6*q*0.01)) / (log((0.01/0.005)*((q-1)*sin(pi/q)/pi*exp(pi/q))))}
uniroot(qfunctionU005, interval = c(4,50))


preciseqvaluesfors <- c(15.83797, 12.50026, 10.20584, 9.091909, 7.326943, 6.427112, 5.994478, 5.593782, 5.354374, 5.189291, 5.06586, 4.889208, 4.765465)
preciseqvaluesforN <- c(6.885207, 7.555128, 8.092186, 8.432907, 8.818877, 9.091909, 9.747214, 10.27389, 10.60856, 10.98813, 11.2569)
preciseqvaluesforU <- c(4.281707, 5.10721, 6.038526, 6.824436, 7.996114, 9.091909, 13.44438, 15.82598)

dataframefordesaifisherbysprecise <- data.frame(
  xvaluesfors,
  preciseqvaluesfors
  
)

dataframefordesaifisherbyNprecise <- data.frame(
  xvaluesforN,
  preciseqvaluesforN
  
)

dataframefordesaifisherbyUprecise <- data.frame(
  xvaluesforU,
  preciseqvaluesforU
  
)

approximatedesaifisherbysgraph <- ggplot(data = dataframefordesaifisherbysapproximate) +
  geom_line(aes(x = xvaluesfors, y = (xvaluesfors*approximateqvaluesfors)), color = "blue") +
  geom_line(aes(x = xvaluesfors, y = (1 - 1/(1 + xvaluesfors*approximateqvaluesfors))), color = "red") +
  scale_y_continuous(limits = c(0, 1.0), sec.axis = sec_axis(~. , name = "Minimum reproductive excess required"))

approximatedesaifisherbysgraph +
  xlab("Selection coefficient") +
  ylab("Proportion of deaths selective") +
  theme_classic(base_size = 15)

approximatedesaifisherbyNgraph <- ggplot(data = dataframefordesaifisherbyNapproximate) +
  geom_line(aes(x = xvaluesforN, y = (0.01*approximateqvaluesforN)), color = "blue") +
  geom_line(aes(x = xvaluesforN, y = (1 - 1/(1 + 0.01*approximateqvaluesforN))), color = "red") +
  scale_y_continuous(limits = c(0, 0.2), sec.axis = sec_axis(~. , name = "Minimum reproductive excess required"))

approximatedesaifisherbyNgraph +
  xlab("Population size") +
  ylab("Proportion of deaths selective") +
  theme_classic(base_size = 15)

approximatedesaifisherbyUgraph <- ggplot(data = dataframefordesaifisherbyUapproximate) +
  geom_line(aes(x = xvaluesforU, y = (0.01*approximateqvaluesforU)), color = "blue") +
  geom_line(aes(x = xvaluesforU, y = (1 - 1/(1 + 0.01*approximateqvaluesforU))), color = "red") +
  scale_y_continuous(limits = c(0, 0.5), sec.axis = sec_axis(~. , name = "Minimum reproductive excess"))

approximatedesaifisherbyUgraph +
  xlab("Beneficial mutation rate") +
  ylab("Proportion of deaths selective") +
  theme_classic(base_size = 15)

precisedesaifisherbysgraph <- ggplot(data = dataframefordesaifisherbysprecise) +
  geom_line(aes(x = xvaluesfors, y = (xvaluesfors*preciseqvaluesfors)), color = "blue", size = 1) +
  geom_line(aes(x = xvaluesfors, y = (1 - 1/(1 + xvaluesfors*preciseqvaluesfors))), color = "red", size = 1) +
  scale_y_continuous(limits = c(0, 1.5), sec.axis = sec_axis(~. , name = "Minimum reproductive excess"))

precisedesaifisherbysgraph +
  xlab("Selection coefficient") +
  ylab("Proportion of deaths selective") +
  theme_classic(base_size = 18)

precisedesaifisherbyNgraph <- ggplot(data = dataframefordesaifisherbyNprecise) +
  geom_line(aes(x = xvaluesforN, y = (0.01*preciseqvaluesforN)), color = "blue", size = 1) +
  geom_line(aes(x = xvaluesforN, y = (1 - 1/(1 + 0.01*preciseqvaluesforN))), color = "red", size = 1) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 0.2), sec.axis = sec_axis(~. , name = "Minimum reproductive excess"))

precisedesaifisherbyNgraph +
  xlab("Population size") +
  ylab("Proportion of deaths selective") +
  theme_classic(base_size = 18)

precisedesaifisherbyUgraph <- ggplot(data = dataframefordesaifisherbyUprecise) +
  geom_line(aes(x = xvaluesforU, y = (0.01*preciseqvaluesforU)), color = "blue", size = 1) +
  geom_line(aes(x = xvaluesforU, y = (1 - 1/(1 + 0.01*preciseqvaluesforU))), color = "red", size = 1) +
  scale_x_continuous(breaks = c(0.0, 0.0005, 0.001, 0.0015, 0.002, 0.0025), labels = c("0.0", "0.0005", "0.001", "0.0015", "0.002", "0.0025")) +
  scale_y_continuous(limits = c(0, 0.2), sec.axis = sec_axis(~. , name = "Minimum reproductive excess"))

precisedesaifisherbyUgraph +
  xlab("Beneficial mutation rate") +
  ylab("Proportion of deaths selective") +
  theme_classic(base_size = 18)

simplifiedprecisedesaifisherbysgraph <- ggplot(data = dataframefordesaifisherbysprecise) +
  geom_line(aes(x = xvaluesfors, y = (xvaluesfors*preciseqvaluesfors)), color = "blue", size = 1) +
  geom_line(aes(x = xvaluesfors, y = (1 - 1/(1 + xvaluesfors*preciseqvaluesfors))), color = "red", size = 1) +
  scale_y_continuous(limits = c(0, 1.5))

simplifiedprecisedesaifisherbysgraph +
  xlab("Selection coefficient") +
  ylab("RE or % deaths") +
  theme_classic(base_size = 19)

simplifiedprecisedesaifisherbyNgraph <- ggplot(data = dataframefordesaifisherbyNprecise) +
  geom_line(aes(x = xvaluesforN, y = (0.01*preciseqvaluesforN)), color = "blue", size = 1) +
  geom_line(aes(x = xvaluesforN, y = (1 - 1/(1 + 0.01*preciseqvaluesforN))), color = "red", size = 1) +
  scale_x_log10(labels = fancy_scientific) +
  scale_y_continuous(limits = c(0, 0.2))

simplifiedprecisedesaifisherbyNgraph +
  xlab("Population size") +
  ylab("RE or % deaths") +
  theme_classic(base_size = 19)

simplifiedprecisedesaifisherbyUgraph <- ggplot(data = dataframefordesaifisherbyUprecise) +
  geom_line(aes(x = xvaluesforU, y = (0.01*preciseqvaluesforU)), color = "blue", size = 1) +
  geom_line(aes(x = xvaluesforU, y = (1 - 1/(1 + 0.01*preciseqvaluesforU))), color = "red", size = 1) +
  scale_x_continuous(breaks = c(0.0, 0.0005, 0.001, 0.0015, 0.002, 0.0025), labels = c("0.0", "0.0005", "0.001", "0.0015", "0.002", "0.0025")) +
  scale_y_continuous(limits = c(0, 0.2))

simplifiedprecisedesaifisherbyUgraph +
  xlab("Beneficial mutation rate") +
  ylab("RE or % deaths") +
  theme_classic(base_size = 19)

dataframeforselectivedeathsvsreproductiveexcessgraph <- data.frame(
  xvalues <- c(kMLO, kMLP, kMHO, kMHP, kTLO, kTLP, kTHO, kTHP),
  y1values <- c(0.5544, 0.5335, 0.0, 0.0338, 0.0796, 0.4952, 0.0003, 0.0),
  y2values <- c(0.8731, 0.9591, 0.3834, 0.9474, 0.5877, 0.9591, 0.4215, 0.5853),
  birthscolorwater <- c("red", "red", "blue", "blue", "red", "red", "blue", "blue"),
  birthslinetypepop <- c("dashed", "solid", "dashed", "solid", "dashed", "solid", "dashed", "solid"),
  birthsfillpop <- c("white", "red", "white", "blue", "white", "red", "white", "blue"),
  birthsshapecity <- c("square", "square", "square", "square", "triangle", "triangle", "triangle", "triangle"),
  birthssizeouter <- c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5),
  birthssizeinner <- c(2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2)
)

scatterplotlowdensitydeathsdataframe <- data.frame(
  environmental_condition <- c("MLI", "MHI", "TLI", "THI"),
  reproductive_excess <- c(kMLO, kMHO, kTLO, kTHO),
  proportion_of_deaths_selective <- c(0.664, 0.0847, 0.545, 0.109),
  colorwater <- c("red", "blue", "red", "blue"),
  fillpop <- c("white", "white", "white", "white"),
  shapecity <- c("square", "square", "triangle", "triangle"),
  sizeouter <- c(3.5, 3.5, 3.5, 3.5),
  sizeinner <- c(2.2, 2.2, 2.2, 2.2)
)

scatterplothighdensitydeathsdataframe <- data.frame(
  environmental_condition <- c("MLP", "MHP", "TLP", "THP"),
  reproductive_excess <- c(kMLP, kMHP, kTLP, kTHP),
  proportion_of_deaths_selective <- c(0.953, 0.556, 0.628, 0.372),
  colorwater <- c("red", "blue", "red", "blue"),
  shapecity <- c("circle", "circle", "triangle", "triangle"),
  sizeouter <- c(3, 3, 3, 3)
)

scatterplotlowdensitydeathsdataframe <- data.frame(
  environmental_condition = c("MLI", "MHI", "TLI", "THI"),
  reproductive_excess_fecundity_low = c(reproductive_excess_fecundity_lowdispersal_MLO / sum(survivalselectivedeathsMLO$Replicates), reproductive_excess_fecundity_lowdispersal_MHO / sum(survivalselectivedeathsMHO$Replicates), reproductive_excess_fecundity_lowdispersal_TLO / sum(survivalselectivedeathsTLO$Replicates), reproductive_excess_fecundity_lowdispersal_THO / sum(survivalselectivedeathsTHO$Replicates)),
  reproductive_excess_fecundity_high = c(reproductive_excess_fecundity_highdispersal_MLO / sum(survivalselectivedeathsMLO$Replicates), reproductive_excess_fecundity_highdispersal_MHO / sum(survivalselectivedeathsMHO$Replicates), reproductive_excess_fecundity_highdispersal_TLO / sum(survivalselectivedeathsTLO$Replicates), reproductive_excess_fecundity_highdispersal_THO / sum(survivalselectivedeathsTHO$Replicates)),
  reproductive_excess_survival_low = c(reproductive_excess_survival_lowdispersal_MLO / sum(survivalselectivedeathsMLO$Replicates), reproductive_excess_survival_lowdispersal_MHO / sum(survivalselectivedeathsMHO$Replicates), reproductive_excess_survival_lowdispersal_TLO / sum(survivalselectivedeathsTLO$Replicates), reproductive_excess_survival_lowdispersal_THO / sum(survivalselectivedeathsTHO$Replicates)),
  reproductive_excess_survival_high = c(reproductive_excess_survival_highdispersal_MLO / sum(survivalselectivedeathsMLO$Replicates), reproductive_excess_survival_highdispersal_MHO / sum(survivalselectivedeathsMHO$Replicates), reproductive_excess_survival_highdispersal_TLO / sum(survivalselectivedeathsTLO$Replicates), reproductive_excess_survival_highdispersal_THO / sum(survivalselectivedeathsTHO$Replicates)),
  proportion_of_deaths_selective = c(proportion_deaths_selective_deaths_MLO, proportion_deaths_selective_deaths_MHO, proportion_deaths_selective_deaths_TLO, proportion_deaths_selective_deaths_THO),
  colorwater <- c("red", "blue", "red", "blue"),
  fillpop <- c("white", "white", "white", "white"),
  shapecity <- c("circle", "circle", "triangle", "triangle"),
  sizeouter <- c(6.5, 6.5, 6.5, 6.5),
  sizeinner <- c(3.2, 3.2, 3.2, 3.2)
)

scatterplothighdensitydeathsdataframe <- data.frame(
  environmental_condition = c("MLP", "MHP", "TLP", "THP"),
  reproductive_excess_fecundity_low = c(reproductive_excess_fecundity_lowdispersal_MLP / sum(survivalselectivedeathsMLP$Replicates), reproductive_excess_fecundity_lowdispersal_MHP / sum(survivalselectivedeathsMHP$Replicates), reproductive_excess_fecundity_lowdispersal_TLP / sum(survivalselectivedeathsTLP$Replicates), reproductive_excess_fecundity_lowdispersal_THP / sum(survivalselectivedeathsTHP$Replicates)),
  reproductive_excess_fecundity_high = c(reproductive_excess_fecundity_highdispersal_MLP / sum(survivalselectivedeathsMLP$Replicates), reproductive_excess_fecundity_highdispersal_MHP / sum(survivalselectivedeathsMHP$Replicates), reproductive_excess_fecundity_highdispersal_TLP / sum(survivalselectivedeathsTLP$Replicates), reproductive_excess_fecundity_highdispersal_THP / sum(survivalselectivedeathsTHP$Replicates)),
  reproductive_excess_survival_low = c(reproductive_excess_survival_lowdispersal_MLP / sum(survivalselectivedeathsMLP$Replicates), reproductive_excess_survival_lowdispersal_MHP / sum(survivalselectivedeathsMHP$Replicates), reproductive_excess_survival_lowdispersal_TLP / sum(survivalselectivedeathsTLP$Replicates), reproductive_excess_survival_lowdispersal_THP / sum(survivalselectivedeathsTHP$Replicates)),
  reproductive_excess_survival_high = c(reproductive_excess_survival_highdispersal_MLP / sum(survivalselectivedeathsMLP$Replicates), reproductive_excess_survival_highdispersal_MHP / sum(survivalselectivedeathsMHP$Replicates), reproductive_excess_survival_highdispersal_TLP / sum(survivalselectivedeathsTLP$Replicates), reproductive_excess_survival_highdispersal_THP / sum(survivalselectivedeathsTHP$Replicates)),
  proportion_of_deaths_selective = c(proportion_deaths_selective_deaths_MLP, proportion_deaths_selective_deaths_MHP, proportion_deaths_selective_deaths_TLP, proportion_deaths_selective_deaths_THP),
  colorwater <- c("red", "blue", "red", "blue"),
  fillpop <- c("white", "white", "white", "white"),
  shapecity <- c("circle", "circle", "triangle", "triangle"),
  sizeouter <- c(6.5, 6.5, 6.5, 6.5),
  sizeinner <- c(3.2, 3.2, 3.2, 3.2)
)

propdeathsselectivevsreproductiveexcessfecundityandsurvivallowdensitygraph <- ggplot(data = scatterplotlowdensitydeathsdataframe) +
  geom_segment(aes(y = reproductive_excess_fecundity_low, x = proportion_of_deaths_selective, yend = reproductive_excess_fecundity_high, xend = proportion_of_deaths_selective, color = colorwater)) +
  geom_segment(aes(y = reproductive_excess_survival_low*20000, x = proportion_of_deaths_selective, yend = reproductive_excess_survival_high*20000, xend = proportion_of_deaths_selective, color = colorwater)) +
  geom_point(aes(y = reproductive_excess_fecundity_low, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(y = reproductive_excess_fecundity_high, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(y = reproductive_excess_fecundity_low, x = proportion_of_deaths_selective, color = fillpop, shape = shapecity, size = sizeinner)) +
  geom_point(aes(y = reproductive_excess_fecundity_high, x = proportion_of_deaths_selective, color = fillpop, shape = shapecity, size = sizeinner)) +
  geom_point(aes(y = reproductive_excess_survival_low*20000, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(y = reproductive_excess_survival_high*20000, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_hline(yintercept = 0.0, color = 'black', linetype = 'dashed', size = 1.5) +
  geom_vline(xintercept = 0.1, color = 'black', linetype = 'dotted', size = 1.5) +
  scale_color_identity() +
  scale_linetype_identity() +
  scale_shape_identity() +
  scale_fill_identity() +
  scale_size_identity() +
  scale_y_continuous(limits = c(0, 20000), sec.axis = sec_axis(~. /20000, name = "Excess seedlings per pot"))
propdeathsselectivevsreproductiveexcessfecundityandsurvivallowdensitygraph +
  xlim(0.0, 1.0) +
  ylab("Excess seeds produced per pot") +
  xlab("Proportion of seedling deaths selective") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_line(colour = "light gray"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20)) +
  annotate("text", y = 2500, x = 0.35, size = 12, label = "Low density")

propdeathsselectivevsreproductiveexcessfecundityandsurvivalhighdensitygraph <- ggplot(data = scatterplothighdensitydeathsdataframe) +
  geom_segment(aes(y = reproductive_excess_fecundity_low, x = proportion_of_deaths_selective, yend = reproductive_excess_fecundity_high, xend = proportion_of_deaths_selective, color = colorwater)) +
  geom_segment(aes(y = (reproductive_excess_survival_low)*1000, x = proportion_of_deaths_selective, yend = (reproductive_excess_survival_high)*1000, xend = proportion_of_deaths_selective, color = colorwater)) +
  geom_point(aes(y = reproductive_excess_fecundity_low, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(y = reproductive_excess_fecundity_high, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(y = reproductive_excess_fecundity_low, x = proportion_of_deaths_selective, color = fillpop, shape = shapecity, size = sizeinner)) +
  geom_point(aes(y = reproductive_excess_fecundity_high, x = proportion_of_deaths_selective, color = fillpop, shape = shapecity, size = sizeinner)) +
  geom_point(aes(y = (reproductive_excess_survival_low)*1000, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(y = (reproductive_excess_survival_high)*1000, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_hline(yintercept = 0.0, color = 'black', linetype = 'dashed', size = 1.5) +
  geom_vline(xintercept = 0.1, color = 'black', linetype = 'dotted', size = 1.5) +
  scale_color_identity() +
  scale_linetype_identity() +
  scale_shape_identity() +
  scale_fill_identity() +
  scale_size_identity() +
  scale_y_continuous(limits = c(0, 30000), sec.axis = sec_axis(~. /1000, name = "Excess seedlings per pot"))
propdeathsselectivevsreproductiveexcessfecundityandsurvivalhighdensitygraph +
  xlim(0.0, 1.0) +
  ylab("Excess seeds produced per pot") +
  xlab("Proportion of seed and seedling deaths selective") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_line(colour = "light gray"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20)) +
  annotate("text", y = 5000, x = 0.35, size = 12, label = "High density")


propdeathsselectivevsreproductiveexcessfecunditylowdensitygraph <- ggplot(data = scatterplotlowdensitydeathsdataframe) +
  geom_segment(aes(y = reproductive_excess_fecundity_low, x = proportion_of_deaths_selective, yend = reproductive_excess_fecundity_high, xend = proportion_of_deaths_selective, color = colorwater)) +
  geom_point(aes(y = reproductive_excess_fecundity_low, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(y = reproductive_excess_fecundity_high, x = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(y = reproductive_excess_fecundity_low, x = proportion_of_deaths_selective, color = fillpop, shape = shapecity, size = sizeinner)) +
  geom_point(aes(y = reproductive_excess_fecundity_high, x = proportion_of_deaths_selective, color = fillpop, shape = shapecity, size = sizeinner)) +
  scale_color_identity() +
  scale_linetype_identity() +
  scale_shape_identity() +
  scale_fill_identity() +
  scale_size_identity() 
propdeathsselectivevsreproductiveexcessfecunditylowdensitygraph +
  ylim(0, 20000) +
  xlim(0.0, 1.0) +
  ylab("Reproductive excess") +
  xlab("Proportion of deaths selective") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_line(colour = "light gray"),
        panel.grid.minor = element_blank()) +
  annotate("text", y = 14000, x = 0.9, size = 6.5, label = "Fecundity, low density")

propdeathsselectivevsreproductiveexcesssurvivallowdensitygraph <- ggplot(data = scatterplotlowdensitydeathsdataframe) +
  geom_segment(aes(x = reproductive_excess_survival_low, y = proportion_of_deaths_selective, xend = reproductive_excess_survival_high, yend = proportion_of_deaths_selective, color = colorwater)) +
  geom_point(aes(x = reproductive_excess_survival_low, y = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(x = reproductive_excess_survival_high, y = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(x = reproductive_excess_survival_low, y = proportion_of_deaths_selective, color = fillpop, shape = shapecity, size = sizeinner)) +
  geom_point(aes(x = reproductive_excess_survival_high, y = proportion_of_deaths_selective, color = fillpop, shape = shapecity, size = sizeinner)) +
  scale_color_identity() +
  scale_linetype_identity() +
  scale_shape_identity() +
  scale_fill_identity() +
  scale_size_identity()
propdeathsselectivevsreproductiveexcesssurvivallowdensitygraph +
  xlim(0, 1.0) +
  ylim(0.0, 1.0) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_line(colour = "light gray"),
        panel.grid.minor = element_blank()) +
  annotate("text", x = 0.7, y = 0.9, size = 6.5, label = "Seedling survival, low density")

propdeathsselectivevsreproductiveexcessfecundityhighdensitygraph <- ggplot(data = scatterplothighdensitydeathsdataframe) +
  geom_segment(aes(x = reproductive_excess_fecundity_low, y = proportion_of_deaths_selective, xend = reproductive_excess_fecundity_high, yend = proportion_of_deaths_selective, color = colorwater)) +
  geom_point(aes(x = reproductive_excess_fecundity_low, y = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(x = reproductive_excess_fecundity_high, y = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  scale_color_identity() +
  scale_linetype_identity() +
  scale_shape_identity() +
  scale_fill_identity() +
  scale_size_identity()
propdeathsselectivevsreproductiveexcessfecundityhighdensitygraph +
  xlim(0, 25000) +
  ylim(0.0, 1.0) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_line(colour = "light gray"),
        panel.grid.minor = element_blank()) +
  annotate("text", x = 18000, y = 0.9, size = 6.5, label = "Fecundity, high density")


propdeathsselectivevsreproductiveexcesssurvivalhighdensitygraph <- ggplot(data = scatterplothighdensitydeathsdataframe) +
  geom_segment(aes(x = reproductive_excess_survival_low, y = proportion_of_deaths_selective, xend = reproductive_excess_survival_high, yend = proportion_of_deaths_selective, color = colorwater)) +
  geom_point(aes(x = reproductive_excess_survival_low, y = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(x = reproductive_excess_survival_high, y = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  scale_color_identity() +
  scale_linetype_identity() +
  scale_shape_identity() +
  scale_fill_identity() +
  scale_size_identity()
propdeathsselectivevsreproductiveexcesssurvivalhighdensitygraph +
  xlim(-10.0, 20.0) +
  ylim(0.0, 1.0) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_line(colour = "light gray"),
        panel.grid.minor = element_blank()) +
  annotate("text", x = 10, y = 0.9, size = 6.5, label = "Seed survival, high density")


propdeathsselectivevsreproductiveexcessgraph <- ggplot(data = dataframeforselectivedeathsvsreproductiveexcessgraph) +
  geom_segment(aes(x = xvalues, y = y1values, xend = xvalues, yend = y2values, color = birthscolorwater, linetype = birthslinetypepop)) +
  geom_point(aes(x = xvalues, y = y1values, color = birthscolorwater, shape = birthsshapecity, size = birthssizeouter)) +
  geom_point(aes(x = xvalues, y = y2values, color = birthscolorwater, shape = birthsshapecity, size = birthssizeouter)) +
  geom_point(aes(x = xvalues, y = y1values, color = birthsfillpop, shape = birthsshapecity, size = birthssizeinner)) +
  geom_point(aes(x = xvalues, y = y2values, color = birthsfillpop, shape = birthsshapecity, size = birthssizeinner)) +
  scale_color_identity() +
  scale_linetype_identity() +
  scale_shape_identity() +
  scale_fill_identity() +
  scale_size_identity()
propdeathsselectivevsreproductiveexcessgraph +
  xlim(0, 20000) +
  ylim(0.0, 1.0) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective") +
  theme_classic() +
  annotate("text", x = 14000, y = 0.9, size = 4.5, label = "Seed production")

resultsscatterplotlowdensitydeaths <- ggplot(data = scatterplotlowdensitydeathsdataframe) +
  geom_point(aes(x = reproductive_excess, y = proportion_of_deaths_selective, color = colorwater, shape = shapecity, size = sizeouter)) +
  geom_point(aes(x = reproductive_excess, y = proportion_of_deaths_selective, color = fillpop, shape = shapecity, size = sizeinner)) +
  scale_color_identity() +
  scale_shape_identity() +
  scale_fill_identity() +
  scale_size_identity()
resultsscatterplotlowdensitydeaths +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_line(colour = "light gray"),
        panel.grid.minor = element_blank()) +
  xlim(0, 20000) +
  ylim(0.0, 1.0) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective") +
  annotate("text", x = 14000, y = 0.9, size = 7.5, label = "Low-density deaths")


resultsscatterplothighdensitydeaths <- ggplot(data = scatterplothighdensitydeathsdataframe, aes(x = reproductive_excess, y = proportion_of_deaths_selective, label = environmental_condition)) +
  geom_point(aes(color = colorwater, shape = shapecity, size = sizeouter)) +
  scale_color_identity() +
  scale_shape_identity() +
  scale_size_identity()
resultsscatterplothighdensitydeaths +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_line(colour = "light gray"),
        panel.grid.minor = element_blank()) +
  xlim(0, 1800) +
  ylim(0.0, 1.0) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective") +
  annotate("text", x = 1300, y = 0.9, size = 7.5, label = "High-density deaths") 

ggarrange(propdeathsselectivevsreproductiveexcessgraph +
            theme(axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank(), 
                  axis.title.x = element_blank(), 
                  axis.title.y = element_text(size = 15),
                  axis.text.y = element_text(size = 12),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(colour = "light gray"),
                  panel.border = element_rect(fill = NA)) +
            xlim(0.0, 20000.0) +
            ylim(0.0, 1.0) +
            xlab("Reproductive excess") +
            ylab("Proportion of deaths selective") +
            annotate("text", x = 14000, y = 0.9, size = 7.5, label = "Seed production"), 
          resultsscatterplotlowdensitydeaths +
            theme_bw(base_size = 15) +
            theme(panel.grid = element_line(colour = "light gray"),
                  panel.grid.minor = element_blank()) +
            xlim(0.0, 20000.0) +
            ylim(0.0, 1.0) +
            xlab("Reproductive excess") +
            ylab("Proportion of deaths selective") +
            annotate("text", x = 14000, y = 0.9, size = 7.5, label = "Low-density deaths"), 
          nrow = 2, widths = c(1, 1))



#EXTREME-VALUE BIAS ANALYSIS
#The first set of analyses is intended to determine the variance between observed genotype averages (for survival and fecundity separately)
#and compare it to the variance between replicates within genotypes. The between-replicate variance is a proxy (probably a conservative proxy)
#for within-genotype variance.
#Since the observed genotype averages in fact include the true genotypic value as well as within-genotype variance,
#the analysis then attempts to determine the between-genotype variance that would recover the observed variance in genotype averages
#after adding within-genotype variance.
#This is accomplished by taking the observed variance, adding within-genotype variance, and noting the resulting between-genotype variance.
#This is repeated for several dozen smaller values for the 'observed' variance, to obtain several dozen resulting between-genotype variances.
#The between-genotype variances are regressed against the starting variances, and the regression equation is used to solve for
#the starting variance required to result in the variance in observed genotype averages from the experiment.


#I want to estimate the bias in the fittest genotype's death rate
#caused by choosing the lowest observed death rate to be our fittest
#So I'm going to simulate drawing deaths for each genotype from a binomial distribution
#with the observed means and number of plants
#I can then calculate a new 'observed' death rate for every genotype, choose the lowest
#and then compare to its 'real' death rate
#This will give me an idea of how much the process of choosing the lowest observed death rate
#biases towards genotypes who 'got lucky' and had a fitter year than they normally would

meanbiasMLO <- ExtremeValueBiasSurvival(survivalselectivedeathsMLO, 0)
meanbiasMLP <- ExtremeValueBiasSurvival(survivalselectivedeathsMLP, 1)
meanbiasMHO <- ExtremeValueBiasSurvival(survivalselectivedeathsMHO, 0)
meanbiasMHP <- ExtremeValueBiasSurvival(survivalselectivedeathsMHP, 1)
meanbiasTLO <- ExtremeValueBiasSurvival(survivalselectivedeathsTLO, 0)
meanbiasTLP <- ExtremeValueBiasSurvival(survivalselectivedeathsTLP, 1)
meanbiasTHO <- ExtremeValueBiasSurvival(survivalselectivedeathsTHO, 0)
meanbiasTHP <- ExtremeValueBiasSurvival(survivalselectivedeathsTHP, 1)

meanbiasMLO
meanbiasMLP
meanbiasMHO
meanbiasMHP
meanbiasTLO
meanbiasTLP
meanbiasTHO
meanbiasTHP

adjustedsurvivalselectivedeathsTLP <- CalculateSelectiveDeathsSurvivalWithBias(aggregatedwithcolumnsnosurvivingseedsTLP, nreplicatestableTLP, meanbiasTLP, 1)
adjustedsurvivalselectivedeathsTHP <- CalculateSelectiveDeathsSurvivalWithBias(aggregatedwithcolumnsnosurvivingseedsTHP, nreplicatestableTHP, meanbiasTHP, 1)

adjusted_total_selective_deaths_survivalTLP <- sum(adjustedsurvivalselectivedeathsTLP$Selective_deaths_survival)
adjusted_total_selective_deaths_survivalTHP <- sum(adjustedsurvivalselectivedeathsTHP$Selective_deaths_survival)

adjusted_total_selective_deaths_survivalTLP
adjusted_total_selective_deaths_survivalTHP

total_selective_deaths_survival_TLP
total_selective_deaths_survival_THP

adjusted_proportion_deaths_selective_survivalTLP <- adjusted_total_selective_deaths_survivalTLP / (sum(survivalselectivedeathsTLP$Replicates)*30)
adjusted_proportion_deaths_selective_survivalTHP <- adjusted_total_selective_deaths_survivalTHP / (sum(survivalselectivedeathsTHP$Replicates)*30)

adjusted_proportion_deaths_selective_survivalTLP
adjusted_proportion_deaths_selective_survivalTHP

























































































#Need two data frames, one with a row for each genotype, one with a row for each simulation run
simulationdataframeMLO <- data.frame(
  genotypesurvivalrates = survivalselectivedeathsMLO$Genotype_death_rate,
  genotypeID = survivalselectivedeathsMLO$Genotype_id,
  resampledgenotypesurvivalrates = survivalselectivedeathsMLO$Genotype_death_rate
)
resultsdataframeMLO <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeMLP <- data.frame(
  genotypesurvivalrates = survivalselectivedeathsMLP$Genotype_death_rate,
  genotypeID = survivalselectivedeathsMLP$Genotype_id,
  resampledgenotypesurvivalrates = survivalselectivedeathsMLP$Genotype_death_rate
)
resultsdataframeMLP <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeMHO <- data.frame(
  genotypesurvivalrates = survivalselectivedeathsMHO$Genotype_death_rate,
  genotypeID = survivalselectivedeathsMHO$Genotype_id,
  resampledgenotypesurvivalrates = survivalselectivedeathsMHO$Genotype_death_rate
)
resultsdataframeMHO <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeMHP <- data.frame(
  genotypesurvivalrates = survivalselectivedeathsMHP$Genotype_death_rate,
  genotypeID = survivalselectivedeathsMHP$Genotype_id,
  resampledgenotypesurvivalrates = survivalselectivedeathsMHP$Genotype_death_rate
)
resultsdataframeMHP <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeTLO <- data.frame(
  genotypesurvivalrates = survivalselectivedeathsTLO$Genotype_death_rate,
  genotypeID = survivalselectivedeathsTLO$Genotype_id,
  resampledgenotypesurvivalrates = survivalselectivedeathsTLO$Genotype_death_rate
)
resultsdataframeTLO <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeTLP <- data.frame(
  genotypesurvivalrates = survivalselectivedeathsTLP$Genotype_death_rate,
  genotypeID = survivalselectivedeathsTLP$Genotype_id,
  resampledgenotypesurvivalrates = survivalselectivedeathsTLP$Genotype_death_rate
)
resultsdataframeTLP <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeTHO <- data.frame(
  genotypesurvivalrates = survivalselectivedeathsTHO$Genotype_death_rate,
  genotypeID = survivalselectivedeathsTHO$Genotype_id,
  resampledgenotypesurvivalrates = survivalselectivedeathsTHO$Genotype_death_rate
)
resultsdataframeTHO <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeTHP <- data.frame(
  genotypesurvivalrates = survivalselectivedeathsTHP$Genotype_death_rate,
  genotypeID = survivalselectivedeathsTHP$Genotype_id,
  resampledgenotypesurvivalrates = survivalselectivedeathsTHP$Genotype_death_rate
)
resultsdataframeTHP <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)


#Following lines are the main simulation loop, which is run 10,000 times.
#Each time, 'true' survival rates are redrawn from the original 517 genotypes (with replacement)
#Then each genotype randomly draws a number of survivors from a binomial with their survival rate.
#Additional lines at the end of the loop are for diagnostic purposes.
for (i in 1:10000) {
  for(j in 1:nrow(survivalselectivedeathsMLO)) {
    simulationdataframeMLO$genotypesurvivalrates[j] <- (1 - simulationdataframeMLO$genotypesurvivalrates[j])
  }
  for(j in 1:517) {
    simulationdataframeMLO$resampledgenotypesurvivalrates[j] <- (rbinom(1, (survivalselectivedeathsMLO$Replicates[j]), simulationdataframeMLO$genotypesurvivalrates[j]) / (survivalselectivedeathsMLO$Replicates[j]))
  }
  resultsdataframeMLO$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeMLO$resampledgenotypesurvivalrates)
  resultsdataframeMLO$vectoroftruesurvivalrates[i] <- max(simulationdataframeMLO$genotypesurvivalrates)
  resultsdataframeMLO$vectorofbestgenotypeID[i] <- simulationdataframeMLO$genotypeID[which.max(simulationdataframeMLO$genotypesurvivalrates)]
  resultsdataframeMLO$vectorofwithingenotypevariances[i] <- resultsdataframeMLO$vectoroftruesurvivalrates[i]*(1 - resultsdataframeMLO$vectoroftruesurvivalrates[i]) / (survivalselectivedeathsMLO$Replicates[which(survivalselectivedeathsMLO$Genotype_id == resultsdataframeMLO$vectorofbestgenotypeID[i], arr.ind = FALSE)])
}
resultsdataframeMLO$vectorofmaximumsurvivalrates
resultsdataframeMLO$vectoroftruesurvivalrates
resultsdataframeMLO$vectorofbestgenotypeID


for (i in 1:10000) {
  for(j in 1:nrow(survivalselectivedeathsMLP)) {
    simulationdataframeMLP$genotypesurvivalrates[j] <- (1 - simulationdataframeMLP$genotypesurvivalrates[j])
  }
  for(j in 1:nrow(survivalselectivedeathsMLP)) {
    simulationdataframeMLP$resampledgenotypesurvivalrates[j] <- (rbinom(1, (survivalselectivedeathsMLP$Replicates[j])*30, simulationdataframeMLP$genotypesurvivalrates[j]) / (survivalselectivedeathsMLP$Replicates[j] * 30))
  }
  resultsdataframeMLP$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeMLP$resampledgenotypesurvivalrates)
  resultsdataframeMLP$vectoroftruesurvivalrates[i] <- max(simulationdataframeMLP$genotypesurvivalrates)
  resultsdataframeMLP$vectorofbestgenotypeID[i] <- simulationdataframeMLP$genotypeID[which.max(simulationdataframeMLP$genotypesurvivalrates)]
  resultsdataframeMLP$vectorofwithingenotypevariances[i] <- resultsdataframeMLP$vectoroftruesurvivalrates[i]*(1 - resultsdataframeMLP$vectoroftruesurvivalrates[i]) / (survivalselectivedeathsMLP$Replicates[which(survivalselectivedeathsMLP$Genotype_id == resultsdataframeMLP$vectorofbestgenotypeID[i], arr.ind = FALSE)]*30)
}
resultsdataframeMLP$vectorofmaximumsurvivalrates
resultsdataframeMLP$vectoroftruesurvivalrates
resultsdataframeMLP$vectorofbestgenotypeID


for (i in 1:10000) {
  sampleddeathrates <- floor(runif(nrow(survivalselectivedeathsMHO), min = 1, max = nrow(survivalselectivedeathsMHO)))
  for(j in 1:nrow(survivalselectivedeathsMHO)) {
    simulationdataframeMHO$resampledgenotypesurvivalrates[j] <- (1 - survivalselectivedeathsMHO$Genotype_death_rate[sampleddeathrates[j]])
  }
  for(j in 1:517) {
    simulationdataframeMHO$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (survivalselectivedeathsMHO$Replicates[j]), simulationdataframeMHO$resampledgenotypesurvivalrates[j]) / (survivalselectivedeathsMHO$Replicates[j]))
  }
  resultsdataframeMHO$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeMHO$simulatedgenotypesurvivalrates)
  resultsdataframeMHO$vectoroftruesurvivalrates[i] <- max(simulationdataframeMHO$resampledgenotypesurvivalrates)
  resultsdataframeMHO$vectorofbestgenotypeID[i] <- simulationdataframeMHO$genotypeID[which.max(simulationdataframeMHO$simulatedgenotypesurvivalrates)]
  resultsdataframeMHO$vectorofwithingenotypevariances[i] <- resultsdataframeMHO$vectoroftruesurvivalrates[i]*(1 - resultsdataframeMHO$vectoroftruesurvivalrates[i]) / (survivalselectivedeathsMHO$Replicates[which(survivalselectivedeathsMHO$Genotype_id == resultsdataframeMHO$vectorofbestgenotypeID[i], arr.ind = FALSE)])
}
resultsdataframeMHO$vectorofmaximumsurvivalrates
resultsdataframeMHO$vectoroftruesurvivalrates
resultsdataframeMHO$vectorofbestgenotypeID


for (i in 1:10000) {
  sampleddeathrates <- floor(runif(nrow(survivalselectivedeathsMHP), min = 1, max = nrow(survivalselectivedeathsMHP)))
  for(j in 1:nrow(survivalselectivedeathsMHP)) {
    simulationdataframeMHP$resampledgenotypesurvivalrates[j] <- (1 - survivalselectivedeathsMHP$Genotype_death_rate[sampleddeathrates[j]])
  }
  for(j in 1:nrow(survivalselectivedeathsMHP)) {
    simulationdataframeMHP$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (survivalselectivedeathsMHP$Replicates[j])*30, simulationdataframeMHP$resampledgenotypesurvivalrates[j]) / (survivalselectivedeathsMHP$Replicates[j] * 30))
  }
  resultsdataframeMHP$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeMHP$simulatedgenotypesurvivalrates)
  resultsdataframeMHP$vectoroftruesurvivalrates[i] <- max(simulationdataframeMHP$resampledgenotypesurvivalrates)
  resultsdataframeMHP$vectorofbestgenotypeID[i] <- simulationdataframeMHP$genotypeID[which.max(simulationdataframeMHP$simulatedgenotypesurvivalrates)]
  resultsdataframeMHP$vectorofwithingenotypevariances[i] <- resultsdataframeMHP$vectoroftruesurvivalrates[i]*(1 - resultsdataframeMHP$vectoroftruesurvivalrates[i]) / (survivalselectivedeathsMHP$Replicates[which(survivalselectivedeathsMHP$Genotype_id == resultsdataframeMHP$vectorofbestgenotypeID[i], arr.ind = FALSE)]*30)
}
resultsdataframeMHP$vectorofmaximumsurvivalrates
resultsdataframeMHP$vectoroftruesurvivalrates
resultsdataframeMHP$vectorofbestgenotypeID


for (i in 1:10000) {
  sampleddeathrates <- floor(runif(nrow(survivalselectivedeathsTLO), min = 1, max = nrow(survivalselectivedeathsTLO)))
  for(j in 1:nrow(survivalselectivedeathsTLO)) {
    simulationdataframeTLO$resampledgenotypesurvivalrates[j] <- (1 - survivalselectivedeathsTLO$Genotype_death_rate[sampleddeathrates[j]])
  }
  for(j in 1:517) {
    simulationdataframeTLO$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (survivalselectivedeathsTLO$Replicates[j]), simulationdataframeTLO$resampledgenotypesurvivalrates[j]) / (survivalselectivedeathsTLO$Replicates[j]))
  }
  resultsdataframeTLO$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeTLO$simulatedgenotypesurvivalrates)
  resultsdataframeTLO$vectoroftruesurvivalrates[i] <- max(simulationdataframeTLO$resampledgenotypesurvivalrates)
  resultsdataframeTLO$vectorofbestgenotypeID[i] <- simulationdataframeTLO$genotypeID[which.max(simulationdataframeTLO$simulatedgenotypesurvivalrates)]
  resultsdataframeTLO$vectorofwithingenotypevariances[i] <- resultsdataframeTLO$vectoroftruesurvivalrates[i]*(1 - resultsdataframeTLO$vectoroftruesurvivalrates[i]) / (survivalselectivedeathsTLO$Replicates[which(survivalselectivedeathsTLO$Genotype_id == resultsdataframeTLO$vectorofbestgenotypeID[i], arr.ind = FALSE)])
}
resultsdataframeTLO$vectorofmaximumsurvivalrates
resultsdataframeTLO$vectoroftruesurvivalrates
resultsdataframeTLO$vectorofbestgenotypeID


for (i in 1:10000) {
  sampleddeathrates <- floor(runif(nrow(survivalselectivedeathsTLP), min = 1, max = nrow(survivalselectivedeathsTLP)))
  for(j in 1:nrow(survivalselectivedeathsTLP)) {
    simulationdataframeTLP$resampledgenotypesurvivalrates[j] <- (1 - survivalselectivedeathsTLP$Genotype_death_rate[sampleddeathrates[j]])
  }
  for(j in 1:nrow(survivalselectivedeathsTLP)) {
    simulationdataframeTLP$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (survivalselectivedeathsTLP$Replicates[j])*30, simulationdataframeTLP$resampledgenotypesurvivalrates[j]) / (survivalselectivedeathsTLP$Replicates[j] * 30))
  }
  resultsdataframeTLP$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeTLP$simulatedgenotypesurvivalrates)
  resultsdataframeTLP$vectoroftruesurvivalrates[i] <- max(simulationdataframeTLP$resampledgenotypesurvivalrates)
  resultsdataframeTLP$vectorofbestgenotypeID[i] <- simulationdataframeTLP$genotypeID[which.max(simulationdataframeTLP$simulatedgenotypesurvivalrates)]
  resultsdataframeTLP$vectorofwithingenotypevariances[i] <- resultsdataframeTLP$vectoroftruesurvivalrates[i]*(1 - resultsdataframeTLP$vectoroftruesurvivalrates[i]) / (survivalselectivedeathsTLP$Replicates[which(survivalselectivedeathsTLP$Genotype_id == resultsdataframeTLP$vectorofbestgenotypeID[i], arr.ind = FALSE)]*30)
}
resultsdataframeTLP$vectorofmaximumsurvivalrates
resultsdataframeTLP$vectoroftruesurvivalrates
resultsdataframeTLP$vectorofbestgenotypeID


for (i in 1:10000) {
  sampleddeathrates <- floor(runif(nrow(survivalselectivedeathsTHO), min = 1, max = nrow(survivalselectivedeathsTHO)))
  for(j in 1:nrow(survivalselectivedeathsTHO)) {
    simulationdataframeTHO$resampledgenotypesurvivalrates[j] <- (1 - survivalselectivedeathsTHO$Genotype_death_rate[sampleddeathrates[j]])
  }
  for(j in 1:517) {
    simulationdataframeTHO$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (survivalselectivedeathsTHO$Replicates[j]), simulationdataframeTHO$resampledgenotypesurvivalrates[j]) / (survivalselectivedeathsTHO$Replicates[j]))
  }
  resultsdataframeTHO$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeTHO$simulatedgenotypesurvivalrates)
  resultsdataframeTHO$vectoroftruesurvivalrates[i] <- max(simulationdataframeTHO$resampledgenotypesurvivalrates)
  resultsdataframeTHO$vectorofbestgenotypeID[i] <- simulationdataframeTHO$genotypeID[which.max(simulationdataframeTHO$simulatedgenotypesurvivalrates)]
  resultsdataframeTHO$vectorofwithingenotypevariances[i] <- resultsdataframeTHO$vectoroftruesurvivalrates[i]*(1 - resultsdataframeTHO$vectoroftruesurvivalrates[i]) / (survivalselectivedeathsTHO$Replicates[which(survivalselectivedeathsTHO$Genotype_id == resultsdataframeTHO$vectorofbestgenotypeID[i], arr.ind = FALSE)])
}
resultsdataframeTHO$vectorofmaximumsurvivalrates
resultsdataframeTHO$vectoroftruesurvivalrates
resultsdataframeTHO$vectorofbestgenotypeID


for (i in 1:10000) {
  sampleddeathrates <- floor(runif(nrow(survivalselectivedeathsTHP), min = 1, max = nrow(survivalselectivedeathsTHP)))
  for(j in 1:nrow(survivalselectivedeathsTHP)) {
    simulationdataframeTHP$resampledgenotypesurvivalrates[j] <- (1 - survivalselectivedeathsTHP$Genotype_death_rate[sampleddeathrates[j]])
  }
  for(j in 1:nrow(survivalselectivedeathsTHP)) {
    simulationdataframeTHP$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (survivalselectivedeathsTHP$Replicates[j])*30, simulationdataframeTHP$resampledgenotypesurvivalrates[j]) / (survivalselectivedeathsTHP$Replicates[j] * 30))
  }
  resultsdataframeTHP$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeTHP$simulatedgenotypesurvivalrates)
  resultsdataframeTHP$vectoroftruesurvivalrates[i] <- max(simulationdataframeTHP$resampledgenotypesurvivalrates)
  resultsdataframeTHP$vectorofbestgenotypeID[i] <- simulationdataframeTHP$genotypeID[which.max(simulationdataframeTHP$simulatedgenotypesurvivalrates)]
  resultsdataframeTHP$vectorofwithingenotypevariances[i] <- resultsdataframeTHP$vectoroftruesurvivalrates[i]*(1 - resultsdataframeTHP$vectoroftruesurvivalrates[i]) / (survivalselectivedeathsTHP$Replicates[which(survivalselectivedeathsTHP$Genotype_id == resultsdataframeTHP$vectorofbestgenotypeID[i], arr.ind = FALSE)]*30)
}
resultsdataframeTHP$vectorofmaximumsurvivalrates
resultsdataframeTHP$vectoroftruesurvivalrates
resultsdataframeTHP$vectorofbestgenotypeID





biasvectorMLO <- c(1:10000)
for (i in 1:10000) {
  biasvectorMLO[i] = (resultsdataframeMLO$vectorofmaximumsurvivalrates[i] - resultsdataframeMLO$vectoroftruesurvivalrates[i])
}
mean(biasvectorMLO)

biasvectorMLP <- c(1:10000)
for (i in 1:10000) {
  biasvectorMLP[i] = (resultsdataframeMLP$vectorofmaximumsurvivalrates[i] - resultsdataframeMLP$vectoroftruesurvivalrates[i])
}
mean(biasvectorMLP)

biasvectorMHO <- c(1:10000)
for (i in 1:10000) {
  biasvectorMHO[i] = (resultsdataframeMHO$vectorofmaximumsurvivalrates[i] - resultsdataframeMHO$vectoroftruesurvivalrates[i])
}
mean(biasvectorMHO)

biasvectorMHP <- c(1:10000)
for (i in 1:10000) {
  biasvectorMHP[i] = (resultsdataframeMHP$vectorofmaximumsurvivalrates[i] - resultsdataframeMHP$vectoroftruesurvivalrates[i])
}
mean(biasvectorMHP)

biasvectorTLO <- c(1:10000)
for (i in 1:10000) {
  biasvectorTLO[i] = (resultsdataframeTLO$vectorofmaximumsurvivalrates[i] - resultsdataframeTLO$vectoroftruesurvivalrates[i])
}
mean(biasvectorTLO)

biasvectorTLP <- c(1:10000)
for (i in 1:10000) {
  biasvectorTLP[i] = (resultsdataframeTLP$vectorofmaximumsurvivalrates[i] - resultsdataframeTLP$vectoroftruesurvivalrates[i])
}
mean(biasvectorTLP)

biasvectorTHO <- c(1:10000)
for (i in 1:10000) {
  biasvectorTHO[i] = (resultsdataframeTHO$vectorofmaximumsurvivalrates[i] - resultsdataframeTHO$vectoroftruesurvivalrates[i])
}
mean(biasvectorTHO)

biasvectorTHP <- c(1:10000)
for (i in 1:10000) {
  biasvectorTHP[i] = (resultsdataframeTHP$vectorofmaximumsurvivalrates[i] - resultsdataframeTHP$vectoroftruesurvivalrates[i])
}
mean(biasvectorTHP)






#The following analysis for births may end up in fact not being used.
hist(aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Seeds_by_ind)
withingenotypevarianceinseedproductionMLO <- sum(aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Variance) / nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO)
withingenotypevarianceinseedproductionMLO
amonggenotypevarianceinseedproductionMLO <- var(aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Seeds_by_ind)
amonggenotypevarianceinseedproductionMLO
amonggenotypelogvarianceinseedproductionMLO <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionMLO
amonggenotypestdevinseedproductionMLO <- sqrt(amonggenotypevarianceinseedproductionMLO)
amonggenotypestdevinseedproductionMLO
amonggenotypelogstdevMLO <- sqrt(amonggenotypelogvarianceinseedproductionMLO)
amonggenotypelogstdevMLO
amonggenotypelogmeanMLO <- mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Seeds_by_ind))
amonggenotypelogmeanMLO
MLOstdevtotry <- seq(0.06, 0.45, 0.005)
MLOresultingvariances <- c(1.0:79.0)
MLOaddingwithingenotypevariance <- c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO))
for (i in 1:79) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO), amonggenotypelogmeanMLO, MLOstdevtotry[i])
    for (j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO)) {
    sizeparam <- resampledamonggenotypetruebirths[n] / ((resampledamonggenotypetruebirths[n]*(0.6874)^2.0) - 1)
    MLOaddingwithingenotypevariance[n] <- (rnbinom(1, mu = resampledamonggenotypetruebirths[n], size = sizeparam))
    if(is.na(MLOaddingwithingenotypevariance[n]) == TRUE) {
      MLOaddingwithingenotypevariance[n] <- 0.0
    }
  }
  MLOresultingvariances[i] <- var(MLOaddingwithingenotypevariance)
}
MLOresultingvariances

plot(MLOstdevtotry, MLOresultingvariances)
MLOregressionvariancebystdev <- lm(MLOresultingvariances ~ MLOstdevtotry)
MLOregressionvariancebystdev

#4785617 = (12774194 * correctstdev) + 3757333
MLOinferredlogbirthstdev <- (4785617 - 3757333) / 12774194
MLOinferredlogbirthstdev
#This analysis predicts that fitting a log-normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) would strongly overestimate the resulting variance. It infers that a 
#log-normal distribution with standard deviation = 0.08049698 for the among-genotype variation
#would recover the observed distribution after adding in within-genotype variance.


hist(aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Seeds_by_ind)
withingenotypevarianceinseedproductionMLP <- sum(aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Variance) / nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLP)
amonggenotypevarianceinseedproductionMLP <- var(aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Seeds_by_ind)
withingenotypevarianceinseedproductionMLP
amonggenotypevarianceinseedproductionMLP
amonggenotypelogvarianceinseedproductionMLP <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionMLP
amonggenotypestdevinseedproductionMLP <- sqrt(amonggenotypevarianceinseedproductionMLP)
amonggenotypestdevinseedproductionMLP
amonggenotypelogstdevMLP <- sqrt(amonggenotypelogvarianceinseedproductionMLP)
amonggenotypelogstdevMLP
amonggenotypelogmeanMLP <- mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Seeds_by_ind))
amonggenotypelogmeanMLP
MLPstdevtotry <- seq(0.36, 0.75, 0.005)
MLPresultingvariances <- c(1.0:79.0)
MLPaddingwithingenotypevariance <- c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLP))
for (i in 1:79) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLP), amonggenotypelogmeanMLP, MLPstdevtotry[i])
    for (j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLP)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLP)) {
    sizeparam <- resampledamonggenotypetruebirths[n] / ((resampledamonggenotypetruebirths[n]*(0.9006)^2.0) - 1)
    MLPaddingwithingenotypevariance[n] <- (rnbinom(1, mu = resampledamonggenotypetruebirths[n], size = sizeparam))
    if(is.na(MLPaddingwithingenotypevariance[n]) == TRUE) {
      MLPaddingwithingenotypevariance[n] <- 0.0
    }
  }
  MLPresultingvariances[i] <- var(MLPaddingwithingenotypevariance)
}
MLPresultingvariances

plot(MLPstdevtotry, MLPresultingvariances)
MLPregressionvariancebystdev <- lm(MLPresultingvariances ~ MLPstdevtotry)
MLPregressionvariancebystdev

#968725.5 = (16859744 * correctstdev) - 5947954
MLPinferredlogbirthstdev <- (968725.5 + 5947954) / 16859744
MLPinferredlogbirthstdev
#This analysis predicts that fitting a log-normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) would strongly overestimate the resulting variance. It infers that a 
#log-normal distribution with standard deviation = 0.4102482 for the among-genotype variation
#would recover the observed distribution after adding in within-genotype variance.


hist(aggregatedonlysurvivingseedsnozeroesinvarianceMHO$Seeds_by_ind)
withingenotypevarianceinseedproductionMHO <- sum(aggregatedonlysurvivingseedsnozeroesinvarianceMHO$Variance) / nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHO)
amonggenotypevarianceinseedproductionMHO <- var(aggregatedonlysurvivingseedsnozeroesinvarianceMHO$Seeds_by_ind)
withingenotypevarianceinseedproductionMHO
amonggenotypevarianceinseedproductionMHO
amonggenotypestdevinseedproductionMHO <- sqrt(amonggenotypevarianceinseedproductionMHO)
amonggenotypestdevinseedproductionMHO
amonggenotypemeanMHO <- mean(aggregatedonlysurvivingseedsnozeroesinvarianceMHO$Seeds_by_ind)
amonggenotypemeanMHO
MHOstdevtotry <- seq(3451.0, 3550.0, 1.0)
MHOresultingvariances <- c(1.0:100.0)
MHOaddingwithingenotypevariance <- c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHO))
for (i in MHOstdevtotry) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHO), amonggenotypemeanMHO, i)
    for (j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHO)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHO)) {
    sizeparam <- resampledamonggenotypetruebirths[n] / ((resampledamonggenotypetruebirths[n]*(0.3828)^2.0) - 1)
    MHOaddingwithingenotypevariance[n] <- (rnbinom(1, mu = resampledamonggenotypetruebirths[n], size = sizeparam))
    if(is.na(MHOaddingwithingenotypevariance[n]) == TRUE) {
      MHOaddingwithingenotypevariance[n] <- 0.0
    }
  }
  MHOresultingvariances[i-3450] <- var(MHOaddingwithingenotypevariance)
}
MHOresultingvariances

testmean <- 18318.39
testsizeparam <- testmean / ((testmean*(0.3828)^2.0) - 1)
testMHOwithnoamonggenotypevariance <- rnbinom(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHO), mu = testmean, size = testsizeparam)
testMHOwithnoamonggenotypevariance
var(testMHOwithnoamonggenotypevariance)
hist(testMHOwithnoamonggenotypevariance)

plot(MHOstdevtotry, MHOresultingvariances)
MHOregressionvariancebystdev <- lm(MHOresultingvariances ~ MHOstdevtotry)
MHOregressionvariancebystdev

#12215481 = (25504 * correctstdev) - 25539054
MHOinferredbirthstdev <- (12215481 + 25539054) / 25504
MHOinferredbirthstdev
#This analysis predicts that fitting a normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) would strongly overestimate the resulting variance. It infers that a 
#normal distribution with standard deviation = 1480.338 for the among-genotype variation
#would recover the observed distribution after adding in within-genotype variance.


hist(aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Seeds_by_ind)
withingenotypevarianceinseedproductionMHP <- sum(aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Variance) / nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHP)
withingenotypevarianceinseedproductionMHP
amonggenotypevarianceinseedproductionMHP <- var(aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Seeds_by_ind)
amonggenotypevarianceinseedproductionMHP
amonggenotypelogvarianceinseedproductionMHP <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionMHP
amonggenotypelogstdevMHP <- sqrt(amonggenotypelogvarianceinseedproductionMHP)
amonggenotypelogstdevMHP
amonggenotypestdevinseedproductionMHP <- sqrt(amonggenotypevarianceinseedproductionMHP)
amonggenotypestdevinseedproductionMHP
amonggenotypelogmeanMHP <-mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Seeds_by_ind))
amonggenotypelogmeanMHP
MHPstdevtotry <- seq(0.36, 0.75, 0.005)
MHPresultingvariances <- c(1.0:79.0)
MHPaddingwithingenotypevariance <- c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHP))
for (i in 1:79) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHP), amonggenotypelogmeanMHP, MHPstdevtotry[i])
    for (j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHP)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHP)) {
    sizeparam <- resampledamonggenotypetruebirths[n] / ((resampledamonggenotypetruebirths[n]*(0.9413)^2.0) - 1)
    MHPaddingwithingenotypevariance[n] <- (rnbinom(1, mu = resampledamonggenotypetruebirths[n], size = sizeparam))
    if(is.na(MHPaddingwithingenotypevariance[n]) == TRUE) {
      MHPaddingwithingenotypevariance[n] <- 0.0
    }
  }
  MHPresultingvariances[i] <- var(MHPaddingwithingenotypevariance)
}
MHPresultingvariances

plot(MHPstdevtotry, MHPresultingvariances)
MHPregressionvariancebystdev <- lm(MHPresultingvariances ~ MHPstdevtotry)
MHPregressionvariancebystdev

#23575776 = (161019519 * correctstdev) - 34748183
MHPinferredlogbirthstdev <- (23575776 + 34748183) / 161019519
MHPinferredlogbirthstdev
#This analysis predicts that fitting a log-normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) would strongly overestimate the resulting variance. It infers that a 
#log-normal distribution with standard deviation = 0.3622167 for the among-genotype variation
#would recover the observed distribution after adding in within-genotype variance.


hist(aggregatedonlysurvivingseedsnozeroesinvarianceTLO$Seeds_by_ind)
withingenotypevarianceinseedproductionTLO <- sum(aggregatedonlysurvivingseedsnozeroesinvarianceTLO$Variance) / nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLO)
amonggenotypevarianceinseedproductionTLO <- var(aggregatedonlysurvivingseedsnozeroesinvarianceTLO$Seeds_by_ind)
withingenotypevarianceinseedproductionTLO
amonggenotypevarianceinseedproductionTLO
amonggenotypestdevinseedproductionTLO <- sqrt(amonggenotypevarianceinseedproductionTLO)
amonggenotypestdevinseedproductionTLO
amonggenotypemeanTLO <- mean(aggregatedonlysurvivingseedsnozeroesinvarianceTLO$Seeds_by_ind)
amonggenotypemeanTLO
TLOstdevtotry <- seq(1201.0, 1285.0, 1.0)
TLOresultingvariances <- c(1.0:85.0)
TLOaddingwithingenotypevariance <- c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLO))
for (i in TLOstdevtotry) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLO), amonggenotypemeanTLO, i)
    for (j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLO)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLO)) {
    sizeparam <- resampledamonggenotypetruebirths[n] / ((resampledamonggenotypetruebirths[n]*(0.3651)^2.0) - 1)
    TLOaddingwithingenotypevariance[n] <- (rnbinom(1, mu = resampledamonggenotypetruebirths[n], size = sizeparam))
    if(is.na(TLOaddingwithingenotypevariance[n]) == TRUE) {
      TLOaddingwithingenotypevariance[n] <- 0.0
    }
  }
  TLOresultingvariances[i-1200] <- var(TLOaddingwithingenotypevariance)
}
TLOresultingvariances

plot(TLOstdevtotry, TLOresultingvariances)
TLOregressionvariancebystdev <- lm(TLOresultingvariances ~ TLOstdevtotry)
TLOregressionvariancebystdev

#1654510 = (4383 * correctstdev) - 583572
TLOinferredbirthstdev <- (1654510 + 583572) / 4383
TLOinferredbirthstdev
#This analysis predicts that fitting a normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) would strongly overestimate the resulting variance. It infers that a 
#normal distribution with standard deviation = 510.6279 for the among-genotype variation
#would recover the observed distribution after adding in within-genotype variance.





withingenotypevarianceinseedproductionTLP <- sum(aggregatedonlysurvivingseedsnozeroesinvarianceTLP$Variance) / nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLP)
amonggenotypevarianceinseedproductionTLP <- var(aggregatedonlysurvivingseedsnozeroesinvarianceTLP$Seeds_by_ind)
withingenotypevarianceinseedproductionTLP
amonggenotypevarianceinseedproductionTLP
amonggenotypelogvarianceinseedproductionTLP <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceTLP$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionTLP
amonggenotypestdevinseedproductionTLP <- sqrt(amonggenotypevarianceinseedproductionTLP)
amonggenotypestdevinseedproductionTLP
amonggenotypelogstdevTLP <- sqrt(amonggenotypelogvarianceinseedproductionTLP)
amonggenotypelogstdevTLP
amonggenotypelogmeanTLP <- mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceTLP$Seeds_by_ind))
amonggenotypelogmeanTLP
TLPstdevtotry <- seq(0.46, 0.85, 0.005)
TLPresultingvariances <- c(1.0:79.0)
TLPaddingwithingenotypevariance <- c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLP))
for (i in 1:79) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLP), amonggenotypelogmeanTLP, TLPstdevtotry[i])
    for (j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLP)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLP)) {
    sizeparam <- resampledamonggenotypetruebirths[n] / ((resampledamonggenotypetruebirths[n]*(0.9543)^2.0) - 1)
    TLPaddingwithingenotypevariance[n] <- (rnbinom(1, mu = resampledamonggenotypetruebirths[n], size = sizeparam))
    if(is.na(TLPaddingwithingenotypevariance[n]) == TRUE) {
      TLPaddingwithingenotypevariance[n] <- 0.0
    }
  }
  TLPresultingvariances[i] <- var(TLPaddingwithingenotypevariance)
}
TLPresultingvariances

plot(TLPstdevtotry, TLPresultingvariances)
TLPregressionvariancebystdev <- lm(TLPresultingvariances ~ TLPstdevtotry)
TLPregressionvariancebystdev

#641572.4 = (2356297 * correctstdev) - 701603
TLPinferredlogbirthstdev <- (641572.4 + 701603) / 2356297
TLPinferredlogbirthstdev
#This analysis predicts that fitting a log-normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) would strongly overestimate the resulting variance. It infers that a 
#log-normal distribution with standard deviation = 0.5700365 for the among-genotype variation
#would recover the observed distribution after adding in within-genotype variance.




withingenotypevarianceinseedproductionTHO <- sum(aggregatedonlysurvivingseedsnozeroesinvarianceTHO$Variance) / nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHO)
amonggenotypevarianceinseedproductionTHO <- var(aggregatedonlysurvivingseedsnozeroesinvarianceTHO$Seeds_by_ind)
withingenotypevarianceinseedproductionTHO
amonggenotypevarianceinseedproductionTHO
withingenotypevarianceinseedproductionTHO/amonggenotypevarianceinseedproductionTHO
amonggenotypestdevinseedproductionTHO <- sqrt(amonggenotypevarianceinseedproductionTHO)
amonggenotypestdevinseedproductionTHO
amonggenotypemeanTHO <- mean(aggregatedonlysurvivingseedsnozeroesinvarianceTHO$Seeds_by_ind)
amonggenotypemeanTHO
THOstdevtotry <- seq(2001.0, 2100.0, 1.0)
THOresultingvariances <- c(1.0:100.0)
THOaddingwithingenotypevariance <- c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHO))
for (i in THOstdevtotry) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHO), amonggenotypemeanTHO, i)
    for (j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHO)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHO)) {
    sizeparam <- resampledamonggenotypetruebirths[n] / ((resampledamonggenotypetruebirths[n]*(0.3197)^2.0) - 1)
    THOaddingwithingenotypevariance[n] <- (rnbinom(1, mu = resampledamonggenotypetruebirths[n], size = sizeparam))
    if(is.na(THOaddingwithingenotypevariance[n]) == TRUE) {
      THOaddingwithingenotypevariance[n] <- 0.0
    }
  }
  THOresultingvariances[i-2000] <- var(THOaddingwithingenotypevariance)
}
THOresultingvariances

plot(THOstdevtotry, THOresultingvariances)
THOregressionvariancebystdev <- lm(THOresultingvariances ~ THOstdevtotry)
THOregressionvariancebystdev

#4492312 = (8442 * correctstdev) - 469226
THOinferredbirthstdev <- (4492312 + 469226) / 8442
THOinferredbirthstdev
#This analysis predicts that fitting a normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) would strongly overestimate the resulting variance. It infers that a 
#normal distribution with standard deviation = 587.7207 for the among-genotype variation
#would recover the observed distribution after adding in within-genotype variance.





hist(aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Seeds_by_ind)
hist(log(aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Seeds_by_ind))
withingenotypevarianceinseedproductionTHP <- sum(aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Variance) / nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHP)
amonggenotypevarianceinseedproductionTHP <- var(aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Seeds_by_ind)
withingenotypevarianceinseedproductionTHP
amonggenotypevarianceinseedproductionTHP
withingenotypevarianceinseedproductionTHP/amonggenotypevarianceinseedproductionTHP
amonggenotypelogvarianceinseedproductionTHP <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionTHP
amonggenotypestdevinseedproductionTHP <- sqrt(amonggenotypevarianceinseedproductionTHP)
amonggenotypestdevinseedproductionTHP
amonggenotypelogstdevTHP <- sqrt(amonggenotypelogvarianceinseedproductionTHP)
amonggenotypelogstdevTHP
amonggenotypelogmeanTHP <- mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Seeds_by_ind))
amonggenotypelogmeanTHP
THPstdevtotry <- seq(0.291, 0.39, 0.001)
THPresultingvariances <- c(1.0:100.0)
THPaddingwithingenotypevariance <- c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHP))
for (i in 1:100) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHP), amonggenotypelogmeanTHP, THPstdevtotry[i])
    for (j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHP)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHP)) {
    sizeparam <- resampledamonggenotypetruebirths[n] / ((resampledamonggenotypetruebirths[n]*(0.862)^2.0) - 1)
    THPaddingwithingenotypevariance[n] <- (rnbinom(1, mu = resampledamonggenotypetruebirths[n], size = sizeparam))
    if(is.na(THPaddingwithingenotypevariance[n]) == TRUE) {
      THPaddingwithingenotypevariance[n] <- 0.0
    }
  }
  THPresultingvariances[i] <- var(THPaddingwithingenotypevariance)
}
THPresultingvariances

plot(THPstdevtotry, THPresultingvariances)
THPregressionvariancebystdev <- lm(THPresultingvariances ~ THPstdevtotry)
THPregressionvariancebystdev

#568233.8 = (2638069 * correctstdev) + 192106
THPinferredlogbirthstdev <- (568233.8 - 192106) / 2638069
THPinferredlogbirthstdev
#This analysis predicts that fitting a log-normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) would strongly overestimate the resulting variance. It infers that a 
#log-normal distribution with standard deviation = 0.1425769 for the among-genotype variation
#would recover the observed distribution after adding in within-genotype variance.




simulationdataframebirthsMLO <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO)),
  genotypeID = aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Genotype_id,
  resampledgenotypebirthrates = aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Seeds_by_ind
)
resultsdataframebirthsMLO <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsMHO <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHO)),
  genotypeID = aggregatedonlysurvivingseedsnozeroesinvarianceMHO$Genotype_id,
  resampledgenotypebirthrates = aggregatedonlysurvivingseedsnozeroesinvarianceMHO$Seeds_by_ind
)
resultsdataframebirthsMHO <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsMLP <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLP)),
  genotypeID = aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Genotype_id,
  resampledgenotypebirthrates = aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Seeds_by_ind
)
resultsdataframebirthsMLP <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsMHP <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHP)),
  genotypeID = aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Genotype_id,
  resampledgenotypebirthrates = aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Seeds_by_ind
)
resultsdataframebirthsMHP <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsTLO <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLO)),
  genotypeID = aggregatedonlysurvivingseedsnozeroesinvarianceTLO$Genotype_id,
  resampledgenotypebirthrates = aggregatedonlysurvivingseedsnozeroesinvarianceTLO$Seeds_by_ind
)
resultsdataframebirthsTLO <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsTHO <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHO)),
  genotypeID = aggregatedonlysurvivingseedsnozeroesinvarianceTHO$Genotype_id,
  resampledgenotypebirthrates = aggregatedonlysurvivingseedsnozeroesinvarianceTHO$Seeds_by_ind
)
resultsdataframebirthsTHO <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsTLP <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLP)),
  genotypeID = aggregatedonlysurvivingseedsnozeroesinvarianceTLP$Genotype_id,
  resampledgenotypebirthrates = aggregatedonlysurvivingseedsnozeroesinvarianceTLP$Seeds_by_ind
)
resultsdataframebirthsTLP <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsTHP <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHP)),
  genotypeID = aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Genotype_id,
  resampledgenotypebirthrates = aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Seeds_by_ind
)
resultsdataframebirthsTHP <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

#The following analyses simulate resampling
#1) 517 true genotypic values from a distribution inferred in the analysis above, and
#2) within-genotype noise from a negative binomial distribution with mean equal to the genotypic value and
#variance calculated from the observed relationship between mean and within-genotype variance from the experiment.


#Note that MLO draws from a log-normal.

amonggenotypelogvarianceinseedproductionMLO <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionMLO
amonggenotypestdevinseedproductionMLO <- sqrt(amonggenotypevarianceinseedproductionMLO)
amonggenotypestdevinseedproductionMLO
amonggenotypelogmeanMLO <-mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceMLO$Seeds_by_ind))
amonggenotypelogmeanMLO
for (i in 1:10000) {
  simulationdataframebirthsMLO$resampledgenotypebirthrates <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO), amonggenotypelogmeanMLO, MLOinferredlogbirthstdev)
  for(j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO)) {
    sizeparam <- simulationdataframebirthsMLO$resampledgenotypebirthrates[j] / ((simulationdataframebirthsMLO$resampledgenotypebirthrates[j]*(0.6874)^2.0) - 1)
    if (simulationdataframebirthsMLO$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsMLO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = 1.0, size = (1.0 / ((0.6874^2.0) - 1.0))))
    } else {
      simulationdataframebirthsMLO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = simulationdataframebirthsMLO$resampledgenotypebirthrates[j], size = sizeparam))
      if(is.na(simulationdataframebirthsMLO$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsMLO$simulatedgenotypebirthrates[j] <- 0.0
      }
    }
  }
  resultsdataframebirthsMLO$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsMLO$simulatedgenotypebirthrates)
  resultsdataframebirthsMLO$vectoroftruebirthrates[i] <- max(simulationdataframebirthsMLO$resampledgenotypebirthrates)
  resultsdataframebirthsMLO$vectorofbestgenotypeID[i] <- simulationdataframebirthsMLO$genotypeID[which.max(simulationdataframebirthsMLO$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsMLO$vectorofmaximumbirthrates
resultsdataframebirthsMLO$vectoroftruebirthrates
resultsdataframebirthsMLO$vectorofbestgenotypeID
hist(resultsdataframebirthsMLO$vectorofmaximumbirthrates)
hist(resultsdataframebirthsMLO$vectoroftruebirthrates)


#Note that MHO draws from a normal.

amonggenotypevarianceinseedproductionMHO <- var(aggregatedonlysurvivingseedsnozeroesinvarianceMHO$Seeds_by_ind)
amonggenotypevarianceinseedproductionMHO
amonggenotypestdevinseedproductionMHO <- sqrt(amonggenotypevarianceinseedproductionMHO)
amonggenotypestdevinseedproductionMHO
amonggenotypemeanMHO <- mean(aggregatedonlysurvivingseedsnozeroesinvarianceMHO$Seeds_by_ind)
amonggenotypemeanMHO
for (i in 1:10000) {
  simulationdataframebirthsMHO$resampledgenotypebirthrates <- rnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHO), amonggenotypemeanMHO, MHOinferredbirthstdev)
  for(j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHO)) {
    sizeparam <- simulationdataframebirthsMHO$resampledgenotypebirthrates[j] / ((simulationdataframebirthsMHO$resampledgenotypebirthrates[j]*(0.3828)^2.0) - 1)
    if (simulationdataframebirthsMHO$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsMHO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = 1.0, size = (1.0 / ((0.3828^2.0) - 1.0))))
      if(is.na(simulationdataframebirthsMHO$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsMHO$simulatedgenotypebirthrates[j] <- 0.0
      }
    } else {
      simulationdataframebirthsMHO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = simulationdataframebirthsMHO$resampledgenotypebirthrates[j], size = sizeparam))
      if(is.na(simulationdataframebirthsMHO$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsMHO$simulatedgenotypebirthrates[j] <- 0.0
      }
    }
  }
  resultsdataframebirthsMHO$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsMHO$simulatedgenotypebirthrates, na.rm = TRUE)
  resultsdataframebirthsMHO$vectoroftruebirthrates[i] <- max(simulationdataframebirthsMHO$resampledgenotypebirthrates, na.rm = TRUE)
  resultsdataframebirthsMHO$vectorofbestgenotypeID[i] <- simulationdataframebirthsMHO$genotypeID[which.max(simulationdataframebirthsMHO$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsMHO$vectorofmaximumbirthrates
resultsdataframebirthsMHO$vectoroftruebirthrates
resultsdataframebirthsMHO$vectorofbestgenotypeID
hist(resultsdataframebirthsMHO$vectorofmaximumbirthrates)
hist(resultsdataframebirthsMHO$vectoroftruebirthrates)


#Note that MLP draws from a log-normal.

amonggenotypelogvarianceinseedproductionMLP <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionMLP
amonggenotypestdevinseedproductionMLP <- sqrt(amonggenotypevarianceinseedproductionMLP)
amonggenotypestdevinseedproductionMLP
amonggenotypelogmeanMLP <-mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceMLP$Seeds_by_ind))
amonggenotypelogmeanMLP
for (i in 1:10000) {
  simulationdataframebirthsMLP$resampledgenotypebirthrates <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLP), amonggenotypelogmeanMLP, MLPinferredlogbirthstdev)
  for(j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLP)) {
    sizeparam <- simulationdataframebirthsMLP$resampledgenotypebirthrates[j] / ((simulationdataframebirthsMLP$resampledgenotypebirthrates[j]*(0.9006)^2.0) - 1)
    if (simulationdataframebirthsMLP$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsMLP$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = 1.0, size = (1.0 / ((0.9006^2.0) - 1.0))))
    } else {
      simulationdataframebirthsMLP$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = simulationdataframebirthsMLP$resampledgenotypebirthrates[j], size = sizeparam))
      if(is.na(simulationdataframebirthsMLP$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsMLP$simulatedgenotypebirthrates[j] <- 0.0
      }
    }
  }
  resultsdataframebirthsMLP$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsMLP$simulatedgenotypebirthrates)
  resultsdataframebirthsMLP$vectoroftruebirthrates[i] <- max(simulationdataframebirthsMLP$resampledgenotypebirthrates)
  resultsdataframebirthsMLP$vectorofbestgenotypeID[i] <- simulationdataframebirthsMLP$genotypeID[which.max(simulationdataframebirthsMLP$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsMLP$vectorofmaximumbirthrates
resultsdataframebirthsMLP$vectoroftruebirthrates
resultsdataframebirthsMLP$vectorofbestgenotypeID
hist(resultsdataframebirthsMLP$vectorofmaximumbirthrates)
hist(resultsdataframebirthsMLP$vectoroftruebirthrates)


#Note that MHT draws from a log-normal.

amonggenotypelogvarianceinseedproductionMHP <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionMHP
amonggenotypestdevinseedproductionMHP <- sqrt(amonggenotypevarianceinseedproductionMHP)
amonggenotypestdevinseedproductionMHP
amonggenotypelogmeanMHP <-mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceMHP$Seeds_by_ind))
amonggenotypelogmeanMHP
for (i in 1:10000) {
  simulationdataframebirthsMHP$resampledgenotypebirthrates <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHP), amonggenotypelogmeanMHP, MHPinferredlogbirthstdev)
  for(j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMHP)) {
    sizeparam <- simulationdataframebirthsMHP$resampledgenotypebirthrates[j] / ((simulationdataframebirthsMHP$resampledgenotypebirthrates[j]*(0.9413)^2.0) - 1)
    if (simulationdataframebirthsMHP$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsMHP$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = 1.0, size = (1.0 / ((0.9413^2.0) - 1.0))))
    } else {
      simulationdataframebirthsMHP$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = simulationdataframebirthsMHP$resampledgenotypebirthrates[j], size = sizeparam))
    }
  }
  resultsdataframebirthsMHP$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsMHP$simulatedgenotypebirthrates)
  resultsdataframebirthsMHP$vectoroftruebirthrates[i] <- max(simulationdataframebirthsMHP$resampledgenotypebirthrates)
  resultsdataframebirthsMHP$vectorofbestgenotypeID[i] <- simulationdataframebirthsMHP$genotypeID[which.max(simulationdataframebirthsMHP$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsMHP$vectorofmaximumbirthrates
resultsdataframebirthsMHP$vectoroftruebirthrates
resultsdataframebirthsMHP$vectorofbestgenotypeID
hist(resultsdataframebirthsMHP$vectorofbestgenotypeID)
hist(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_id)
ave(resultsdataframebirthsMHP$vectorofmaximumbirthrates)
hist(resultsdataframebirthsMHP$vectorofmaximumbirthrates)
hist(resultsdataframebirthsMHP$vectoroftruebirthrates)


#Note that TLO draws from a normal.

amonggenotypevarianceinseedproductionTLO <- var(aggregatedonlysurvivingseedsnozeroesinvarianceTLO$Seeds_by_ind)
amonggenotypevarianceinseedproductionTLO
amonggenotypestdevinseedproductionTLO <- sqrt(amonggenotypevarianceinseedproductionTLO)
amonggenotypestdevinseedproductionTLO
amonggenotypemeanTLO <- mean(aggregatedonlysurvivingseedsnozeroesinvarianceTLO$Seeds_by_ind)
amonggenotypemeanTLO
for (i in 1:10000) {
  simulationdataframebirthsTLO$resampledgenotypebirthrates <- rnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLO), amonggenotypemeanTLO, TLOinferredbirthstdev)
  for(j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLO)) {
    sizeparam <- simulationdataframebirthsTLO$resampledgenotypebirthrates[j] / ((simulationdataframebirthsTLO$resampledgenotypebirthrates[j]*(0.3651)^2.0) - 1)
    if (simulationdataframebirthsTLO$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsTLO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = 1.0, size = (1.0 / ((0.3651^2.0) - 1.0))))
      if(is.na(simulationdataframebirthsTLO$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsTLO$simulatedgenotypebirthrates[j] <- 0.0
      }
    } else {
      simulationdataframebirthsTLO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = simulationdataframebirthsTLO$resampledgenotypebirthrates[j], size = sizeparam))
      if(is.na(simulationdataframebirthsTLO$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsTLO$simulatedgenotypebirthrates[j] <- 0.0
      }
    }
  }
  resultsdataframebirthsTLO$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsTLO$simulatedgenotypebirthrates, na.rm = TRUE)
  resultsdataframebirthsTLO$vectoroftruebirthrates[i] <- max(simulationdataframebirthsTLO$resampledgenotypebirthrates, na.rm = TRUE)
  resultsdataframebirthsTLO$vectorofbestgenotypeID[i] <- simulationdataframebirthsTLO$genotypeID[which.max(simulationdataframebirthsTLO$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsTLO$vectorofmaximumbirthrates
resultsdataframebirthsTLO$vectoroftruebirthrates
resultsdataframebirthsTLO$vectorofbestgenotypeID
hist(resultsdataframebirthsTLO$vectorofmaximumbirthrates)
hist(resultsdataframebirthsTLO$vectoroftruebirthrates)

#Note that TLP draws from a log-normal.

amonggenotypelogvarianceinseedproductionTLP <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceTLP$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionTLP
amonggenotypestdevinseedproductionTLP <- sqrt(amonggenotypevarianceinseedproductionTLP)
amonggenotypestdevinseedproductionTLP
amonggenotypelogmeanTLP <-mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceTLP$Seeds_by_ind))
amonggenotypelogmeanTLP
for (i in 1:10000) {
  simulationdataframebirthsTLP$resampledgenotypebirthrates <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLP), amonggenotypelogmeanTLP, TLPinferredlogbirthstdev)
  for(j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTLP)) {
    sizeparam <- simulationdataframebirthsTLP$resampledgenotypebirthrates[j] / ((simulationdataframebirthsTLP$resampledgenotypebirthrates[j]*(0.9543)^2.0) - 1)
    if (simulationdataframebirthsTLP$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsTLP$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = 1.0, size = (1.0 / ((0.9543^2.0) - 1.0))))
    } else {
      simulationdataframebirthsTLP$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = simulationdataframebirthsTLP$resampledgenotypebirthrates[j], size = sizeparam))
      if(is.na(simulationdataframebirthsTLP$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsTLP$simulatedgenotypebirthrates[j] <- 0.0
      }
    }
  }
  resultsdataframebirthsTLP$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsTLP$simulatedgenotypebirthrates)
  resultsdataframebirthsTLP$vectoroftruebirthrates[i] <- max(simulationdataframebirthsTLP$resampledgenotypebirthrates)
  resultsdataframebirthsTLP$vectorofbestgenotypeID[i] <- simulationdataframebirthsTLP$genotypeID[which.max(simulationdataframebirthsTLP$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsTLP$vectorofmaximumbirthrates
resultsdataframebirthsTLP$vectoroftruebirthrates
resultsdataframebirthsTLP$vectorofbestgenotypeID
hist(resultsdataframebirthsTLP$vectorofmaximumbirthrates)
hist(resultsdataframebirthsTLP$vectoroftruebirthrates)


#Note that THO draws from a normal.

amonggenotypevarianceinseedproductionTHO <- var(aggregatedonlysurvivingseedsnozeroesinvarianceTHO$Seeds_by_ind)
amonggenotypevarianceinseedproductionTHO
amonggenotypestdevinseedproductionTHO <- sqrt(amonggenotypevarianceinseedproductionTHO)
amonggenotypestdevinseedproductionTHO
amonggenotypemeanTHO <- mean(aggregatedonlysurvivingseedsnozeroesinvarianceTHO$Seeds_by_ind)
amonggenotypemeanTHO
for (i in 1:10000) {
  simulationdataframebirthsTHO$resampledgenotypebirthrates <- rnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHO), amonggenotypemeanTHO, THOinferredbirthstdev)
  for(j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHO)) {
    sizeparam <- simulationdataframebirthsTHO$resampledgenotypebirthrates[j] / ((simulationdataframebirthsTHO$resampledgenotypebirthrates[j]*(0.3197)^2.0) - 1)
    if (simulationdataframebirthsTHO$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsTHO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = 1.0, size = (1.0 / ((0.3197^2.0) - 1.0))))
      if(is.na(simulationdataframebirthsTHO$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsTHO$simulatedgenotypebirthrates[j] <- 0.0
      }
    } else {
      simulationdataframebirthsTHO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = simulationdataframebirthsTHO$resampledgenotypebirthrates[j], size = sizeparam))
      if(is.na(simulationdataframebirthsTHO$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsTHO$simulatedgenotypebirthrates[j] <- 0.0
      }
    }
  }
  resultsdataframebirthsTHO$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsTHO$simulatedgenotypebirthrates, na.rm = TRUE)
  resultsdataframebirthsTHO$vectoroftruebirthrates[i] <- max(simulationdataframebirthsTHO$resampledgenotypebirthrates, na.rm = TRUE)
  resultsdataframebirthsTHO$vectorofbestgenotypeID[i] <- simulationdataframebirthsTHO$genotypeID[which.max(simulationdataframebirthsTHO$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsTHO$vectorofmaximumbirthrates
resultsdataframebirthsTHO$vectoroftruebirthrates
resultsdataframebirthsTHO$vectorofbestgenotypeID
hist(resultsdataframebirthsTHO$vectorofmaximumbirthrates)
hist(resultsdataframebirthsTHO$vectoroftruebirthrates)


#Note that THP draws from a log-normal.

amonggenotypelogvarianceinseedproductionTHP <- var(log(aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Seeds_by_ind))
amonggenotypelogvarianceinseedproductionTHP
amonggenotypestdevinseedproductionTHP <- sqrt(amonggenotypevarianceinseedproductionTHP)
amonggenotypestdevinseedproductionTHP
amonggenotypelogmeanTHP <- mean(log(aggregatedonlysurvivingseedsnozeroesinvarianceTHP$Seeds_by_ind))
amonggenotypelogmeanTHP
for (i in 1:10000) {
  simulationdataframebirthsTHP$resampledgenotypebirthrates <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHP), amonggenotypelogmeanTHP, THPinferredlogbirthstdev)
  for(j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceTHP)) {
    sizeparam <- simulationdataframebirthsTHP$resampledgenotypebirthrates[j] / ((simulationdataframebirthsTHP$resampledgenotypebirthrates[j]*(0.862)^2.0) - 1)
    if (simulationdataframebirthsTHP$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsTHP$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = 1.0, size = (1.0 / ((0.862^2.0) - 1.0))))
    } else {
      simulationdataframebirthsTHP$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = simulationdataframebirthsTHP$resampledgenotypebirthrates[j], size = sizeparam))
      if(is.na(simulationdataframebirthsTHP$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsTHP$simulatedgenotypebirthrates[j] <- 0.0
      }
    }
  }
  resultsdataframebirthsTHP$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsTHP$simulatedgenotypebirthrates)
  resultsdataframebirthsTHP$vectoroftruebirthrates[i] <- max(simulationdataframebirthsTHP$resampledgenotypebirthrates)
  resultsdataframebirthsTHP$vectorofbestgenotypeID[i] <- simulationdataframebirthsTHP$genotypeID[which.max(simulationdataframebirthsTHP$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsTHP$vectorofmaximumbirthrates
resultsdataframebirthsTHP$vectoroftruebirthrates
resultsdataframebirthsTHP$vectorofbestgenotypeID
hist(resultsdataframebirthsTHP$vectorofmaximumbirthrates)
hist(resultsdataframebirthsTHP$vectoroftruebirthrates)


#The following analysis is similar to the above but instead samples true genotypic values
#directly from the observed genotype averages from the experiment.
#To correct for the fact that the observed averages should already include within-genotype variance,
#the analysis adds a set amount of regression to the (observed total) mean for each true genotypic value
#before simulating the addition of within-genotype variance using the same negative binomial distribution as above.


amonggenotypestdevinseedproductionMLO <- sqrt(amonggenotypevarianceinseedproductionMLO)
amonggenotypestdevinseedproductionMLO
for (i in 1:10000) {
  simulationdataframebirthsMLO$resampledgenotypebirthrates <- rlnorm(nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO), amonggenotypelogmeanMLO, MLOinferredlogbirthstdev)
  for(j in 1:nrow(aggregatedonlysurvivingseedsnozeroesinvarianceMLO)) {
    sizeparam <- simulationdataframebirthsMLO$resampledgenotypebirthrates[j] / ((simulationdataframebirthsMLO$resampledgenotypebirthrates[j]*(0.6874)^2.0) - 1)
    if (simulationdataframebirthsMLO$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsMLO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = 1.0, size = (1.0 / ((0.6874^2.0) - 1.0))))
    } else {
      simulationdataframebirthsMLO$simulatedgenotypebirthrates[j] <- (rnbinom(1, mu = simulationdataframebirthsMLO$resampledgenotypebirthrates[j], size = sizeparam))
      if(is.na(simulationdataframebirthsMLO$simulatedgenotypebirthrates[j]) == TRUE) {
        simulationdataframebirthsMLO$simulatedgenotypebirthrates[j] <- 0.0
      }
    }
  }
  resultsdataframebirthsMLO$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsMLO$simulatedgenotypebirthrates)
  resultsdataframebirthsMLO$vectoroftruebirthrates[i] <- max(simulationdataframebirthsMLO$resampledgenotypebirthrates)
  resultsdataframebirthsMLO$vectorofbestgenotypeID[i] <- simulationdataframebirthsMLO$genotypeID[which.max(simulationdataframebirthsMLO$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsMLO$vectorofmaximumbirthrates
resultsdataframebirthsMLO$vectoroftruebirthrates
resultsdataframebirthsMLO$vectorofbestgenotypeID
hist(resultsdataframebirthsMLO$vectorofmaximumbirthrates)
hist(resultsdataframebirthsMLO$vectoroftruebirthrates)



biasvectorbirthsMLO <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsMLO[i] = (resultsdataframebirthsMLO$vectorofmaximumbirthrates[i] - resultsdataframebirthsMLO$vectoroftruebirthrates[i])
}
meanbiasinmaxbirthrateMLO <- mean(biasvectorbirthsMLO)
meanbiasinmaxbirthrateMLO
hist(biasvectorbirthsMLO)

biasvectorbirthsMLP <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsMLP[i] = (resultsdataframebirthsMLP$vectorofmaximumbirthrates[i] - resultsdataframebirthsMLP$vectoroftruebirthrates[i])
}
meanbiasinmaxbirthrateMLP <- mean(biasvectorbirthsMLP)
meanbiasinmaxbirthrateMLP
hist(biasvectorbirthsMLP)

biasvectorbirthsMHO <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsMHO[i] = (resultsdataframebirthsMHO$vectorofmaximumbirthrates[i] - resultsdataframebirthsMHO$vectoroftruebirthrates[i])
}
meanbiasinmaxbirthrateMHO <- mean(biasvectorbirthsMHO)
meanbiasinmaxbirthrateMHO
hist(biasvectorbirthsMHO)

biasvectorbirthsMHP <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsMHP[i] = (resultsdataframebirthsMHP$vectorofmaximumbirthrates[i] - resultsdataframebirthsMHP$vectoroftruebirthrates[i])
}
meanbiasinmaxbirthrateMHP <- mean(biasvectorbirthsMHP)
meanbiasinmaxbirthrateMHP
hist(biasvectorbirthsMHP)

biasvectorbirthsTLO <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsTLO[i] = (resultsdataframebirthsTLO$vectorofmaximumbirthrates[i] - resultsdataframebirthsTLO$vectoroftruebirthrates[i])
}
meanbiasinmaxbirthrateTLO <- mean(biasvectorbirthsTLO)
meanbiasinmaxbirthrateTLO
hist(biasvectorbirthsTLO)

biasvectorbirthsTLP <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsTLP[i] = (resultsdataframebirthsTLP$vectorofmaximumbirthrates[i] - resultsdataframebirthsTLP$vectoroftruebirthrates[i])
}
meanbiasinmaxbirthrateTLP <- mean(biasvectorbirthsTLP)
meanbiasinmaxbirthrateTLP
hist(biasvectorbirthsTLP)

biasvectorbirthsTHO <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsTHO[i] = (resultsdataframebirthsTHO$vectorofmaximumbirthrates[i] - resultsdataframebirthsTHO$vectoroftruebirthrates[i])
}
meanbiasinmaxbirthrateTHO <- mean(biasvectorbirthsTHO)
meanbiasinmaxbirthrateTHO
hist(biasvectorbirthsTHO)

biasvectorbirthsTHP <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsTHP[i] = (resultsdataframebirthsTHP$vectorofmaximumbirthrates[i] - resultsdataframebirthsTHP$vectoroftruebirthrates[i])
}
meanbiasinmaxbirthrateTHP <- mean(biasvectorbirthsTHP)
meanbiasinmaxbirthrateTHP
hist(biasvectorbirthsTHP)




birthsselectivedeathswithbiasMLO <- CalculateSelectiveDeathsBirthsWithBias(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMLO, nreplicatestableomittingreplicateswithnosurvivingseedsMLO, meanbiasinmaxbirthrateMLO, 0)
birthsselectivedeathswithbiasMLP <- CalculateSelectiveDeathsBirthsWithBias(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMLP, nreplicatestableomittingreplicateswithnosurvivingseedsMLP, meanbiasinmaxbirthrateMLP, 1)
birthsselectivedeathswithbiasMHO <- CalculateSelectiveDeathsBirthsWithBias(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMHO, nreplicatestableomittingreplicateswithnosurvivingseedsMHO, meanbiasinmaxbirthrateMHO, 0)
birthsselectivedeathswithbiasMHP <- CalculateSelectiveDeathsBirthsWithBias(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsMHP, nreplicatestableomittingreplicateswithnosurvivingseedsMHP, meanbiasinmaxbirthrateMHP, 1)
birthsselectivedeathswithbiasTLO <- CalculateSelectiveDeathsBirthsWithBias(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTLO, nreplicatestableomittingreplicateswithnosurvivingseedsTLO, meanbiasinmaxbirthrateTLO, 0)
birthsselectivedeathswithbiasTLP <- CalculateSelectiveDeathsBirthsWithBias(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTLP, nreplicatestableomittingreplicateswithnosurvivingseedsTLP, meanbiasinmaxbirthrateTLP, 1)
birthsselectivedeathswithbiasTHO <- CalculateSelectiveDeathsBirthsWithBias(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTHO, nreplicatestableomittingreplicateswithnosurvivingseedsTHO, meanbiasinmaxbirthrateTHO, 0)
birthsselectivedeathswithbiasTHP <- CalculateSelectiveDeathsBirthsWithBias(aggregatedwithcolumnsonlyreplicateswithsurvivingseedsTHP, nreplicatestableomittingreplicateswithnosurvivingseedsTHP, meanbiasinmaxbirthrateTHP, 1)



total_selective_deaths_births_with_bias_MLO <- sum(birthsselectivedeathswithbiasMLO$Selective_deaths_births)
total_selective_deaths_births_with_bias_MLO
total_selective_deaths_births_with_bias_MLP <- sum(birthsselectivedeathswithbiasMLP$Selective_deaths_births)
total_selective_deaths_births_with_bias_MLP
total_selective_deaths_births_with_bias_MHO <- sum(birthsselectivedeathswithbiasMHO$Selective_deaths_births)
total_selective_deaths_births_with_bias_MHO
total_selective_deaths_births_with_bias_MHP <- sum(birthsselectivedeathswithbiasMHP$Selective_deaths_births)
total_selective_deaths_births_with_bias_MHP
total_selective_deaths_births_with_bias_TLO <- sum(birthsselectivedeathswithbiasTLO$Selective_deaths_births)
total_selective_deaths_births_with_bias_TLO
total_selective_deaths_births_with_bias_TLP <- sum(birthsselectivedeathswithbiasTLP$Selective_deaths_births)
total_selective_deaths_births_with_bias_TLP
total_selective_deaths_births_with_bias_THO <- sum(birthsselectivedeathswithbiasTHO$Selective_deaths_births)
total_selective_deaths_births_with_bias_THO
total_selective_deaths_births_with_bias_THP <- sum(birthsselectivedeathswithbiasTHP$Selective_deaths_births)
total_selective_deaths_births_with_bias_THP


max(birthsselectivedeathsMLO$Seeds_by_ind)
max(birthsselectivedeathsMLO$Seeds_by_ind) - meanbiasinmaxbirthrateMLO
max(birthsselectivedeathsMLP$Seeds_by_ind)
max(birthsselectivedeathsMLP$Seeds_by_ind) - meanbiasinmaxbirthrateMLP
max(birthsselectivedeathsMHO$Seeds_by_ind)
max(birthsselectivedeathsMHO$Seeds_by_ind) - meanbiasinmaxbirthrateMHO
max(birthsselectivedeathsMHP$Seeds_by_ind)
max(birthsselectivedeathsMHP$Seeds_by_ind) - meanbiasinmaxbirthrateMHP
max(birthsselectivedeathsTLO$Seeds_by_ind)
max(birthsselectivedeathsTLO$Seeds_by_ind) - meanbiasinmaxbirthrateTLO
max(birthsselectivedeathsTLP$Seeds_by_ind)
max(birthsselectivedeathsTLP$Seeds_by_ind) - meanbiasinmaxbirthrateTLP
max(birthsselectivedeathsTHO$Seeds_by_ind)
max(birthsselectivedeathsTHO$Seeds_by_ind) - meanbiasinmaxbirthrateTHO
max(birthsselectivedeathsTHP$Seeds_by_ind)
max(birthsselectivedeathsTHP$Seeds_by_ind) - meanbiasinmaxbirthrateTHP

mean(birthsselectivedeathsMLO$Seeds_by_ind)
mean(birthsselectivedeathsMLP$Seeds_by_ind)
mean(birthsselectivedeathsMHO$Seeds_by_ind)
mean(birthsselectivedeathsMHP$Seeds_by_ind)
mean(birthsselectivedeathsTLO$Seeds_by_ind)
mean(birthsselectivedeathsTLP$Seeds_by_ind)
mean(birthsselectivedeathsTHO$Seeds_by_ind)
mean(birthsselectivedeathsTHP$Seeds_by_ind)

adjusted_proportion_deaths_selective_births_MLO <- total_selective_deaths_births_with_bias_MLO / (sum(birthsselectivedeathsMLO$Surviving_individuals_pot) * max(birthsselectivedeathswithbiasMLO$Seeds_by_ind) - total_seeds_for_next_generation_MLO)
adjusted_proportion_deaths_selective_births_MLO
adjusted_proportion_deaths_selective_births_MLP <- total_selective_deaths_births_with_bias_MLP / (sum(birthsselectivedeathsMLP$Surviving_individuals_pot) * max(birthsselectivedeathswithbiasMLP$Seeds_by_ind) - total_seeds_for_next_generation_MLP)
adjusted_proportion_deaths_selective_births_MLP
adjusted_proportion_deaths_selective_births_MHO <- total_selective_deaths_births_with_bias_MHO / (sum(birthsselectivedeathsMHO$Surviving_individuals_pot) * max(birthsselectivedeathswithbiasMHO$Seeds_by_ind) - total_seeds_for_next_generation_MHO)
adjusted_proportion_deaths_selective_births_MHO
adjusted_proportion_deaths_selective_births_MHP <- total_selective_deaths_births_with_bias_MHP / (sum(birthsselectivedeathsMHP$Surviving_individuals_pot) * max(birthsselectivedeathswithbiasMHP$Seeds_by_ind) - total_seeds_for_next_generation_MHP)
adjusted_proportion_deaths_selective_births_MHP
adjusted_proportion_deaths_selective_births_TLO <- total_selective_deaths_births_with_bias_TLO / (sum(birthsselectivedeathsTLO$Surviving_individuals_pot) * max(birthsselectivedeathswithbiasTLO$Seeds_by_ind) - total_seeds_for_next_generation_TLO)
adjusted_proportion_deaths_selective_births_TLO
adjusted_proportion_deaths_selective_births_TLP <- total_selective_deaths_births_with_bias_TLP / (sum(birthsselectivedeathsTLP$Surviving_individuals_pot) * max(birthsselectivedeathswithbiasTLP$Seeds_by_ind) - total_seeds_for_next_generation_TLP)
adjusted_proportion_deaths_selective_births_TLP
adjusted_proportion_deaths_selective_births_THO <- total_selective_deaths_births_with_bias_THO / (sum(birthsselectivedeathsTHO$Surviving_individuals_pot) * max(birthsselectivedeathswithbiasTHO$Seeds_by_ind) - total_seeds_for_next_generation_THO)
adjusted_proportion_deaths_selective_births_THO
adjusted_proportion_deaths_selective_births_THP <- total_selective_deaths_births_with_bias_THP / (sum(birthsselectivedeathsTHP$Surviving_individuals_pot) * max(birthsselectivedeathswithbiasTHP$Seeds_by_ind) - total_seeds_for_next_generation_THP)
adjusted_proportion_deaths_selective_births_THP


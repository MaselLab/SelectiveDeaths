library('ggplot2')
library('scales')
library('grid')
library('MASS')
library('ggrepel')
library('egg')
library('writexl')
install.packages('writexl')
install.packages('ggrepel')
install.packages('egg')

getwd()

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

#Note that the datasets above have replicates for each genotype.
#In order to calculate selective deaths, need to know the optimal genotype at each LH stage
#This requires aggregating data from replicates of each genotype
aggregatedmadridlowwateroneseedselectivedeathsdata <- aggregate(madridlowwateroneseedselectivedeathsdata, by = list(madridlowwateroneseedselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedmadridhighwateroneseedselectivedeathsdata <- aggregate(madridhighwateroneseedselectivedeathsdata, by = list(madridhighwateroneseedselectivedeathsdata$Genotype_id), FUN = mean)

aggregatedmadridlowwatertwentyseedsselectivedeathsdata <- aggregate(madridlowwatertwentyseedsselectivedeathsdata, by = list(madridlowwatertwentyseedsselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedmadridhighwatertwentyseedsselectivedeathsdata <- aggregate(madridhighwatertwentyseedsselectivedeathsdata, by = list(madridhighwatertwentyseedsselectivedeathsdata$Genotype_id), FUN = mean)

aggregatedtuebingenlowwateroneseedselectivedeathsdata <- aggregate(tuebingenlowwateroneseedselectivedeathsdata, by = list(tuebingenlowwateroneseedselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedtuebingenhighwateroneseedselectivedeathsdata <- aggregate(tuebingenhighwateroneseedselectivedeathsdata, by = list(tuebingenhighwateroneseedselectivedeathsdata$Genotype_id), FUN = mean)

aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata <- aggregate(tuebingenlowwatertwentyseedsselectivedeathsdata, by = list(tuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_id), FUN = mean)
aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata <- aggregate(tuebingenhighwatertwentyseedsselectivedeathsdata, by = list(tuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_id), FUN = mean)

#After the data has been aggregated into averages by replicate,
#I still need to know the number of replicates to use in the equation for selective deaths.
#This is easily accomplished in a frequency table called on the non-aggregated data
nreplicatestableMLI <- table(madridlowwateroneseedselectivedeathsdata$Genotype_id)
nreplicatestableMHI <- table(madridhighwateroneseedselectivedeathsdata$Genotype_id)
nreplicatestableMLP <- table(madridlowwatertwentyseedsselectivedeathsdata$Genotype_id)
nreplicatestableMHP <- table(madridhighwatertwentyseedsselectivedeathsdata$Genotype_id)
nreplicatestableTLI <- table(tuebingenlowwateroneseedselectivedeathsdata$Genotype_id)
nreplicatestableTHI <- table(tuebingenhighwateroneseedselectivedeathsdata$Genotype_id)
nreplicatestableTLP <- table(tuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_id)
nreplicatestableTHP <- table(tuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_id)

#I sum the number of replicates for each table
#to count up the total number of pots (replicates) in each environmental condition.
#This isn't necessary for later calculations
#It's to use as a denominator/comparison in the final data table.
sum(nreplicatestableMLO)
sum(nreplicatestableMHO)
sum(nreplicatestableMLT)
sum(nreplicatestableMHT)
sum(nreplicatestableTLO)
sum(nreplicatestableTHO)
sum(nreplicatestableTLT)
sum(nreplicatestableTHT)

#I need to add some new columns to the dataset for the main calculation loop.
AddColumns <- function(environmentalconditiondataframe) {
  environmentalconditiondataframe$Replicates <- c(1:nrow(environmentalconditiondataframe))
  environmentalconditiondataframe$Selective_deaths_survival <- c(1:nrow(environmentalconditiondataframe))
  environmentalconditiondataframe$Selective_deaths_births <- c(1:nrow(environmentalconditiondataframe))
  environmentalconditiondataframe$Genotype_death_rate <- c(1:nrow(environmentalconditiondataframe))
  environmentalconditiondataframe$Seeds_for_next_generation <- c(1:nrow(environmentalconditiondataframe))
  environmentalconditiondataframe <- subset(environmentalconditiondataframe, select = -c(Field_site, Watering_treatment, Group.1))
  names(environmentalconditiondataframe)[2] <- "Total_surviving_individuals"
  names(environmentalconditiondataframe)[3] <- "Total_seeds_produced"
  return(environmentalconditiondataframe)
}

fulldataframeMHI <- AddColumns(aggregatedmadridhighwateroneseedselectivedeathsdata)
fulldataframeMHP <- AddColumns(aggregatedmadridhighwatertwentyseedsselectivedeathsdata)
fulldataframeMLI <- AddColumns(aggregatedmadridlowwateroneseedselectivedeathsdata)
fulldataframeMLP <- AddColumns(aggregatedmadridlowwatertwentyseedsselectivedeathsdata)
fulldataframeTHI <- AddColumns(aggregatedtuebingenhighwateroneseedselectivedeathsdata)
fulldataframeTHP <- AddColumns(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata)
fulldataframeTLI <- AddColumns(aggregatedtuebingenlowwateroneseedselectivedeathsdata)
fulldataframeTLP <- AddColumns(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata)

CalculateFullSelectiveDeaths <- function(environmentalconditiondataframe, nreplicatestable, isPopulation) {
  for (i in 1:nrow(environmentalconditiondataframe)) {
    environmentalconditiondataframe$Replicates[i] <- nreplicatestable[i]
    environmentalconditiondataframe$Total_surviving_individuals[i] <- (environmentalconditiondataframe$Total_surviving_individuals[i] * environmentalconditiondataframe$Replicates[i])
    environmentalconditiondataframe$Total_seeds_produced[i] <- (environmentalconditiondataframe$Total_seeds_produced[i] * environmentalconditiondataframe$Replicates[i])
    if (isPopulation == 1) {
      environmentalconditiondataframe$Genotype_death_rate[i] <- (30*environmentalconditiondataframe$Replicates[i] - environmentalconditiondataframe$Total_surviving_individuals[i]) / (30*environmentalconditiondataframe$Replicates[i])
    } else {
      environmentalconditiondataframe$Genotype_death_rate[i] <- (environmentalconditiondataframe$Replicates[i] - environmentalconditiondataframe$Total_surviving_individuals[i]) / (environmentalconditiondataframe$Replicates[i])
    }
    if (environmentalconditiondataframe$Total_surviving_individuals[i] > 0.0) {
      environmentalconditiondataframe$Seeds_by_ind[i] <- environmentalconditiondataframe$Total_seeds_produced[i] / environmentalconditiondataframe$Total_surviving_individuals[i]
    } else {
      environmentalconditiondataframe$Seeds_by_ind[i] <- 0.0
    }
    if (isPopulation == 1) {
      if (environmentalconditiondataframe$Total_seeds_produced[i] > 30.0*environmentalconditiondataframe$Replicates[i]) {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- 30.0*environmentalconditiondataframe$Replicates[i]
      } else {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- environmentalconditiondataframe$Total_seeds_produced[i]
      }
    } else {
      if (environmentalconditiondataframe$Total_seeds_produced[i] > 10.0*environmentalconditiondataframe$Replicates[i]) {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- 10.0*environmentalconditiondataframe$Replicates[i]
      } else {
        environmentalconditiondataframe$Seeds_for_next_generation[i] <- environmentalconditiondataframe$Total_seeds_produced[i]
      }
    }
  }
  max_birth_rate <- max(environmentalconditiondataframe$Seeds_by_ind)
  min_death_rate <- min(environmentalconditiondataframe$Genotype_death_rate)
  for (i in 1:nrow(environmentalconditiondataframe)) {
    if (isPopulation == 1) {
      environmentalconditiondataframe$Selective_deaths_survival[i] <- (30*environmentalconditiondataframe$Replicates[i]) * (environmentalconditiondataframe$Genotype_death_rate[i] - min_death_rate)
    } else {
      environmentalconditiondataframe$Selective_deaths_survival[i] <- (environmentalconditiondataframe$Replicates[i]) * (environmentalconditiondataframe$Genotype_death_rate[i] - min_death_rate)
    }
    environmentalconditiondataframe$Selective_deaths_births[i] <- (max_birth_rate - environmentalconditiondataframe$Seeds_by_ind[i]) * environmentalconditiondataframe$Total_surviving_individuals[i]
  }
  return(environmentalconditiondataframe)
}

fulldataframeMHI <- CalculateFullSelectiveDeaths(fulldataframeMHI, nreplicatestableMHI, 0)
fulldataframeMHP <- CalculateFullSelectiveDeaths(fulldataframeMHP, nreplicatestableMHP, 1)
fulldataframeMLI <- CalculateFullSelectiveDeaths(fulldataframeMLI, nreplicatestableMLI, 0)
fulldataframeMLP <- CalculateFullSelectiveDeaths(fulldataframeMLP, nreplicatestableMLP, 1)
fulldataframeTHI <- CalculateFullSelectiveDeaths(fulldataframeTHI, nreplicatestableTHI, 0)
fulldataframeTHP <- CalculateFullSelectiveDeaths(fulldataframeTHP, nreplicatestableTHP, 1)
fulldataframeTLI <- CalculateFullSelectiveDeaths(fulldataframeTLI, nreplicatestableTLI, 0)
fulldataframeTLP <- CalculateFullSelectiveDeaths(fulldataframeTLP, nreplicatestableTLP, 1)

write_xlsx(fulldataframeMHI, 'Selective deaths/fulldataframeMHI.xlsx')
write_xlsx(fulldataframeMHP, 'Selective deaths/fulldataframeMHP.xlsx')
write_xlsx(fulldataframeMLI, 'Selective deaths/fulldataframeMLI.xlsx')
write_xlsx(fulldataframeMLP, 'Selective deaths/fulldataframeMLP.xlsx')
write_xlsx(fulldataframeTHI, 'Selective deaths/fulldataframeTHI.xlsx')
write_xlsx(fulldataframeTHP, 'Selective deaths/fulldataframeTHP.xlsx')
write_xlsx(fulldataframeTLI, 'Selective deaths/fulldataframeTLI.xlsx')
write_xlsx(fulldataframeTLP, 'Selective deaths/fulldataframeTLP.xlsx')

total_selective_deaths_survival_MHI <- sum(fulldataframeMHI$Selective_deaths_survival)
total_selective_deaths_survival_MHP <- sum(fulldataframeMHP$Selective_deaths_survival)

total_selective_deaths_births_MHI <- sum(fulldataframeMHI$Selective_deaths_births)
total_selective_deaths_births_MHP <- sum(fulldataframeMHP$Selective_deaths_births)
total_selective_deaths_births_MLI <- sum(fulldataframeMLI$Selective_deaths_births)
total_selective_deaths_births_MLP <- sum(fulldataframeMLP$Selective_deaths_births)
total_selective_deaths_births_THI <- sum(fulldataframeTHI$Selective_deaths_births)
total_selective_deaths_births_THP <- sum(fulldataframeTHP$Selective_deaths_births)
total_selective_deaths_births_TLI <- sum(fulldataframeTLI$Selective_deaths_births)
total_selective_deaths_births_TLP <- sum(fulldataframeTLP$Selective_deaths_births)

total_seeds_for_next_generation_MHI <- sum(fulldataframeMHI$Seeds_for_next_generation)
total_seeds_for_next_generation_MHP <- sum(fulldataframeMHP$Seeds_for_next_generation)
total_seeds_for_next_generation_MLI <- sum(fulldataframeMLI$Seeds_for_next_generation)
total_seeds_for_next_generation_MLP <- sum(fulldataframeMLP$Seeds_for_next_generation)
total_seeds_for_next_generation_THI <- sum(fulldataframeTHI$Seeds_for_next_generation)
total_seeds_for_next_generation_THP <- sum(fulldataframeTHP$Seeds_for_next_generation)
total_seeds_for_next_generation_TLI <- sum(fulldataframeTLI$Seeds_for_next_generation)
total_seeds_for_next_generation_TLP <- sum(fulldataframeTLP$Seeds_for_next_generation)

proportion_deaths_selective_births_MHI <- total_selective_deaths_births_MHI / (sum(fulldataframeMHI$Total_surviving_individuals) * max(fulldataframeMHI$Seeds_by_ind) - total_seeds_for_next_generation_MHI)
proportion_deaths_selective_births_MHP <- total_selective_deaths_births_MHP / (sum(fulldataframeMHP$Total_surviving_individuals) * max(fulldataframeMHP$Seeds_by_ind) - total_seeds_for_next_generation_MHP)
proportion_deaths_selective_births_MLI <- total_selective_deaths_births_MLI / (sum(fulldataframeMLI$Total_surviving_individuals) * max(fulldataframeMLI$Seeds_by_ind) - total_seeds_for_next_generation_MLI)
proportion_deaths_selective_births_MLP <- total_selective_deaths_births_MLP / (sum(fulldataframeMLP$Total_surviving_individuals) * max(fulldataframeMLP$Seeds_by_ind) - total_seeds_for_next_generation_MLP)
proportion_deaths_selective_births_THI <- total_selective_deaths_births_THI / (sum(fulldataframeTHI$Total_surviving_individuals) * max(fulldataframeTHI$Seeds_by_ind) - total_seeds_for_next_generation_THI)
proportion_deaths_selective_births_THP <- total_selective_deaths_births_THP / (sum(fulldataframeTHP$Total_surviving_individuals) * max(fulldataframeTHP$Seeds_by_ind) - total_seeds_for_next_generation_THP)
proportion_deaths_selective_births_TLI <- total_selective_deaths_births_TLI / (sum(fulldataframeTLI$Total_surviving_individuals) * max(fulldataframeTLI$Seeds_by_ind) - total_seeds_for_next_generation_TLI)
proportion_deaths_selective_births_TLP <- total_selective_deaths_births_TLP / (sum(fulldataframeTLP$Total_surviving_individuals) * max(fulldataframeTLP$Seeds_by_ind) - total_seeds_for_next_generation_TLP)

proportion_deaths_selective_births_MHI
proportion_deaths_selective_births_MHP
proportion_deaths_selective_births_MLI
proportion_deaths_selective_births_MLP
proportion_deaths_selective_births_THI
proportion_deaths_selective_births_THP
proportion_deaths_selective_births_TLI
proportion_deaths_selective_births_TLP

#The aggregation by replicate averages across all numbers in each replicate.
#The problem is that the seed numbers average across all numbers, including the zeroes.
#Since the zeroes come from pots with no surviving adults, it skews the seed numbers.
#We only want the numbers of seeds produced per surviving adult individual.
#Each of the following for loops fills in the Replicates and both selective deaths columns.
#The equations for selective deaths in each genotype require comparison to the optimal genotype.
#The optimal genotype in the survival case is just 100% survival,
#which I manually confirmed to have occurred in at least one genotype in each dataset.
#The optimal genotype in the birth case is the genotype with the most seeds produced per individual,
#which I calculate directly before the for loop in each condition.
#I then sum the selective deaths columns to find the total number of selective deaths in each LH stage.
for (i in 1:nrow(aggregatedmadridlowwateroneseedselectivedeathsdata)) {
  aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[i] <- nreplicatestableMLO[i]
  aggregatedmadridlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] <- (aggregatedmadridlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] * aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[i])
  aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_death_rate[i] <- (aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[i] - aggregatedmadridlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]) / (aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[i])
  aggregatedmadridlowwateroneseedselectivedeathsdata$Seeds_total_pot[i] <- aggregatedmadridlowwateroneseedselectivedeathsdata$Seeds_total_pot[i] * aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[i]
  if(aggregatedmadridlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] > 0.0) {
    aggregatedmadridlowwateroneseedselectivedeathsdata$Seeds_by_ind[i] <- aggregatedmadridlowwateroneseedselectivedeathsdata$Seeds_total_pot[i] / aggregatedmadridlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]
  } else {
    aggregatedmadridlowwateroneseedselectivedeathsdata$Seeds_by_ind[i] <- 0.0
  }
}
max_birth_rate_MLO <- max(aggregatedmadridlowwateroneseedselectivedeathsdata$Seeds_by_ind)
min_death_rate_MLO <- min(aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_death_rate)
for (i in 1:nrow(aggregatedmadridlowwateroneseedselectivedeathsdata)) {
  aggregatedmadridlowwateroneseedselectivedeathsdata$Selective_deaths_survival[i] <- (aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[i] - aggregatedmadridlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i])
  aggregatedmadridlowwateroneseedselectivedeathsdata$Selective_deaths_births[i] <- (max_birth_rate_MLO - aggregatedmadridlowwateroneseedselectivedeathsdata$Seeds_by_ind[i]) * (aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[i] - aggregatedmadridlowwateroneseedselectivedeathsdata$Selective_deaths_survival[i])
}
total_selective_deaths_survival_MLO <- sum(aggregatedmadridlowwateroneseedselectivedeathsdata$Selective_deaths_survival)
total_selective_deaths_births_MLO <- sum(aggregatedmadridlowwateroneseedselectivedeathsdata$Selective_deaths_births)
#Following line confirms that at least one genotype has the maximum possible survival rate,
#by confirming that at least one genotype has zero selective deaths when calculated this way.
min(aggregatedmadridlowwateroneseedselectivedeathsdata$Selective_deaths_survival)

for (i in 1:nrow(aggregatedmadridhighwateroneseedselectivedeathsdata)) {
  aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[i] <- nreplicatestableMHO[i]
  aggregatedmadridhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] <- (aggregatedmadridhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] * aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[i])
  aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_death_rate[i] <- (aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[i] - aggregatedmadridhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]) / (aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[i])
  aggregatedmadridhighwateroneseedselectivedeathsdata$Seeds_total_pot[i] <- aggregatedmadridhighwateroneseedselectivedeathsdata$Seeds_total_pot[i] * aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[i]
  if(aggregatedmadridhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] > 0.0) {
    aggregatedmadridhighwateroneseedselectivedeathsdata$Seeds_by_ind[i] <- aggregatedmadridhighwateroneseedselectivedeathsdata$Seeds_total_pot[i] / aggregatedmadridhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]
  } else {
    aggregatedmadridhighwateroneseedselectivedeathsdata$Seeds_by_ind[i] <- 0.0
  }
}
max_birth_rate_MHO <- max(aggregatedmadridhighwateroneseedselectivedeathsdata$Seeds_by_ind)
min_death_rate_MHO <- min(aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_death_rate)
for (i in 1:nrow(aggregatedmadridhighwateroneseedselectivedeathsdata)) {
  aggregatedmadridhighwateroneseedselectivedeathsdata$Selective_deaths_survival[i] <- (aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[i] - aggregatedmadridhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i])
  aggregatedmadridhighwateroneseedselectivedeathsdata$Selective_deaths_births[i] <- (max_birth_rate_MHO - aggregatedmadridhighwateroneseedselectivedeathsdata$Seeds_by_ind[i]) * (aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[i] - aggregatedmadridhighwateroneseedselectivedeathsdata$Selective_deaths_survival[i])
}
total_selective_deaths_survival_MHO <- sum(aggregatedmadridhighwateroneseedselectivedeathsdata$Selective_deaths_survival)
total_selective_deaths_births_MHO <- sum(aggregatedmadridhighwateroneseedselectivedeathsdata$Selective_deaths_births)
#Following line confirms that at least one genotype has the maximum possible survival rate,
#by confirming that at least one genotype has zero selective deaths when calculated this way.
min(aggregatedmadridhighwateroneseedselectivedeathsdata$Selective_deaths_survival)

#So turns out there's actually THIRTY seeds in a pot, and I'm too lazy to change my variable names.
#The max survival is therefore 30 per pot, not twenty, hence the 30 in the 2nd line of for loop.
#Manually confirmed that at least one genotype has the max survival rate surviving to adulthood.
for (i in 1:nrow(aggregatedmadridlowwatertwentyseedsselectivedeathsdata)) {
  aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[i] <- nreplicatestableMLT[i]
  aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] <- (aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] * aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[i])
  aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate[i] <- (30*aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[i] - aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) / (30*aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[i])
  aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] <- aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] * aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[i]
  if(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] > 0.0) {
    aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i] <- aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] / aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
  } else {
    aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i] <- 0.0
  }
}
max_birth_rate_MLT <- max(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind)
min_death_rate_MLT <- min(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate)
for (i in 1:nrow(aggregatedmadridlowwatertwentyseedsselectivedeathsdata)) {
  aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Selective_deaths_survival[i] <- ((30*aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[i]) - aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i])
  aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Selective_deaths_births[i] <- (max_birth_rate_MLT - aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i]) * aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
}
total_selective_deaths_survival_MLT <- sum(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Selective_deaths_survival)
total_selective_deaths_births_MLT <- sum(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Selective_deaths_births)
#Following line confirms that at least one genotype has the maximum possible survival rate,
#by confirming that at least one genotype has zero selective deaths when calculated this way.
min(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Selective_deaths_survival)

for (i in 1:nrow(aggregatedmadridhighwatertwentyseedsselectivedeathsdata)) {
  aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[i] <- nreplicatestableMHT[i]
  aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] <- (aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] * aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[i])
  aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate[i] <- (30*aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[i] - aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) / (30*aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[i])
  aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] <- aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] * aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[i]
  if(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] > 0.0) {
    aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i] <- aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] / aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
  } else {
    aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i] <- 0.0
  }
}
max_birth_rate_MHT <- max(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind)
min_death_rate_MHT <- min(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate)
for (i in 1:nrow(aggregatedmadridhighwatertwentyseedsselectivedeathsdata)) {
  aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Selective_deaths_survival[i] <- ((30*aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[i]) - aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i])
  aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Selective_deaths_births[i] <- (max_birth_rate_MHT - aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i]) * aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
}
total_selective_deaths_survival_MHT <- sum(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Selective_deaths_survival)
total_selective_deaths_births_MHT <- sum(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Selective_deaths_births)
#Following line confirms that at least one genotype has the maximum possible survival rate,
#by confirming that at least one genotype has zero selective deaths when calculated this way.
min(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Selective_deaths_survival)

aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate[1] <- 0
for (i in 1:nrow(aggregatedtuebingenlowwateroneseedselectivedeathsdata)) {
  aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[i] <- nreplicatestableTLO[i]
  aggregatedtuebingenlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] <- (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] * aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[i])
  aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate[i] <- (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[i] - aggregatedtuebingenlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]) / (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[i])
  aggregatedtuebingenlowwateroneseedselectivedeathsdata$Seeds_total_pot[i] <- aggregatedtuebingenlowwateroneseedselectivedeathsdata$Seeds_total_pot[i] * aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[i]
  if(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] > 0.0) {
    aggregatedtuebingenlowwateroneseedselectivedeathsdata$Seeds_by_ind[i] <- aggregatedtuebingenlowwateroneseedselectivedeathsdata$Seeds_total_pot[i] / aggregatedtuebingenlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]
  } else {
    aggregatedtuebingenlowwateroneseedselectivedeathsdata$Seeds_by_ind[i] <- 0.0
  }
}
max_birth_rate_TLO <- max(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Seeds_by_ind)
min_death_rate_TLO <- min(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate)
for (i in 1:nrow(aggregatedtuebingenlowwateroneseedselectivedeathsdata)) {
  aggregatedtuebingenlowwateroneseedselectivedeathsdata$Selective_deaths_survival[i] <- (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[i] - aggregatedtuebingenlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i])
  aggregatedtuebingenlowwateroneseedselectivedeathsdata$Selective_deaths_births[i] <- (max_birth_rate_TLO - aggregatedtuebingenlowwateroneseedselectivedeathsdata$Seeds_by_ind[i]) * (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[i] - aggregatedtuebingenlowwateroneseedselectivedeathsdata$Selective_deaths_survival[i])
}
total_selective_deaths_survival_TLO <- sum(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Selective_deaths_survival)
total_selective_deaths_births_TLO <- sum(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Selective_deaths_births)
#Following line confirms that at least one genotype has the maximum possible survival rate,
#by confirming that at least one genotype has zero selective deaths when calculated this way.
min(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Selective_deaths_survival)

aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate[1] <- 0
for (i in 1:nrow(aggregatedtuebingenhighwateroneseedselectivedeathsdata)) {
  aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[i] <- nreplicatestableTHO[i]
  aggregatedtuebingenhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] <- (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] * aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[i])
  aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate[i] <- (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[i] - aggregatedtuebingenhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]) / (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[i])
  aggregatedtuebingenhighwateroneseedselectivedeathsdata$Seeds_total_pot[i] <- aggregatedtuebingenhighwateroneseedselectivedeathsdata$Seeds_total_pot[i] * aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[i]
  if(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i] > 0.0) {
    aggregatedtuebingenhighwateroneseedselectivedeathsdata$Seeds_by_ind[i] <- aggregatedtuebingenhighwateroneseedselectivedeathsdata$Seeds_total_pot[i] / aggregatedtuebingenhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]
  } else {
    aggregatedtuebingenhighwateroneseedselectivedeathsdata$Seeds_by_ind[i] <- 0.0
  }
}
max_birth_rate_THO <- max(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Seeds_by_ind)
min_death_rate_THO <- min(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate)
for (i in 1:nrow(aggregatedtuebingenhighwateroneseedselectivedeathsdata)) {
  aggregatedtuebingenhighwateroneseedselectivedeathsdata$Selective_deaths_survival[i] <- (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[i] - aggregatedtuebingenhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i])
  aggregatedtuebingenhighwateroneseedselectivedeathsdata$Selective_deaths_births[i] <- (max_birth_rate_THO - aggregatedtuebingenhighwateroneseedselectivedeathsdata$Seeds_by_ind[i]) * (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[i] - aggregatedtuebingenhighwateroneseedselectivedeathsdata$Selective_deaths_survival[i])
}
total_selective_deaths_survival_THO <- sum(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Selective_deaths_survival)
total_selective_deaths_births_THO <- sum(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Selective_deaths_births)
#Following line confirms that at least one genotype has the maximum possible survival rate,
#by confirming that at least one genotype has zero selective deaths when calculated this way.
min(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Selective_deaths_survival)

#So turns out there's actually THIRTY seeds in a pot, and I'm too lazy to change my variable names.
#In the following two treatments, no genotype has the maximum survival rate to adulthood,
#so things need to be changed.
aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate[1] <- 0
for (i in 1:nrow(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata)) {
  aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[i] <- nreplicatestableTLT[i]
  aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] <- (aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] * aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[i])
  aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate[i] <- (30*aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[i] - aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) / (30*aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[i])
  aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] <- aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] * aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[i]
  if(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] > 0.0) {
    aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i] <- aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] / aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
  } else {
    aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i] <- 0.0
  }
}
max_birth_rate_TLT <- max(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind)
min_death_rate_TLT <- min(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate)
for (i in 1:nrow(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata)) {
  aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Selective_deaths_survival[i] <- ((30*aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[i]) - aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) - (30*aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[i]*min_death_rate_TLT)
  aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Selective_deaths_births[i] <- (max_birth_rate_TLT - aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i]) * aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
}
total_selective_deaths_survival_TLT <- sum(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Selective_deaths_survival)
total_selective_deaths_births_TLT <- sum(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Selective_deaths_births)
#Following line confirms that at least one genotype has the maximum survival rate,
#by confirming that at least one genotype has zero selective deaths when calculated this way.
min(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Selective_deaths_survival)

aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate[1] <- 0
for (i in 1:nrow(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata)) {
  aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[i] <- nreplicatestableTHT[i]
  aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] <- (aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] * aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[i])
  aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate[i] <- (30*aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[i] - aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) / (30*aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[i])  
  aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] <- aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] * aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[i]
  if(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i] > 0.0) {
    aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i] <- aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_total_pot[i] / aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
  } else {
    aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i] <- 0.0
  }
}
max_birth_rate_THT <- max(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind)
min_death_rate_THT <- min(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate)
for (i in 1:nrow(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata)) {
  aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Selective_deaths_survival[i] <- ((30*aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[i]) - aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) - (30*aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[i]*min_death_rate_THT)
  aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Selective_deaths_births[i] <- (max_birth_rate_THT - aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i]) * aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
}
total_selective_deaths_survival_THT <- sum(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Selective_deaths_survival)
total_selective_deaths_births_THT <- sum(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Selective_deaths_births)
#Following line confirms that at least one genotype has the maximum survival rate,
#by confirming that at least one genotype has zero selective deaths.
min(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Selective_deaths_survival)

#I'm also interested in the overall death rate to see how it's related to the number of selective deaths.
total_starting_seeds_MLO <- nrow(madridlowwateroneseedselectivedeathsdata)
total_surviving_adults_MLO <- sum(madridlowwateroneseedselectivedeathsdata$Surviving_individuals_pot)
average_adult_death_rate_MLO <- 1 - (total_surviving_adults_MLO / total_starting_seeds_MLO)
proportion_deaths_selective_MLO <- total_selective_deaths_survival_MLO / total_starting_seeds_MLO

total_seeds_possible_MLO <- (total_surviving_adults_MLO * max_birth_rate_MLO)
total_seeds_produced_MLO <- sum(madridlowwateroneseedselectivedeathsdata$Seeds_total_pot)
average_seed_death_rate_MLO <- 1 - (total_seeds_produced_MLO / total_seeds_possible_MLO)

total_starting_seeds_MHO <- nrow(madridhighwateroneseedselectivedeathsdata)
total_surviving_adults_MHO <- sum(madridhighwateroneseedselectivedeathsdata$Surviving_individuals_pot)
average_adult_death_rate_MHO <- 1 - (total_surviving_adults_MHO / total_starting_seeds_MHO)
proportion_deaths_selective_MHO <- total_selective_deaths_survival_MHO / total_starting_seeds_MHO

total_seeds_possible_MHO <- (total_surviving_adults_MHO * max_birth_rate_MHO)
total_seeds_produced_MHO <- sum(madridhighwateroneseedselectivedeathsdata$Seeds_total_pot)
average_seed_death_rate_MHO <- 1 - (total_seeds_produced_MHO / total_seeds_possible_MHO)

total_starting_seeds_MLT <- 30*nrow(madridlowwatertwentyseedsselectivedeathsdata)
total_surviving_adults_MLT <- sum(madridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot)
average_adult_death_rate_MLT <- 1 - (total_surviving_adults_MLT / total_starting_seeds_MLT)
proportion_deaths_selective_MLT <- total_selective_deaths_survival_MLT / total_starting_seeds_MLT

total_seeds_possible_MLT <- (total_surviving_adults_MLT * max_birth_rate_MLT)
total_seeds_produced_MLT <- sum(madridlowwatertwentyseedsselectivedeathsdata$Seeds_total_pot)
average_seed_death_rate_MLT <- 1 - (total_seeds_produced_MLT / total_seeds_possible_MLT)

total_starting_seeds_MHT <- 30*nrow(madridhighwatertwentyseedsselectivedeathsdata)
total_surviving_adults_MHT <- sum(madridhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot)
average_adult_death_rate_MHT <- 1 - (total_surviving_adults_MHT / total_starting_seeds_MHT)
proportion_deaths_selective_MHT <- total_selective_deaths_survival_MHT / total_starting_seeds_MHT

total_seeds_possible_MHT <- (total_surviving_adults_MHT * max_birth_rate_MHT)
total_seeds_produced_MHT <- sum(madridhighwatertwentyseedsselectivedeathsdata$Seeds_total_pot)
average_seed_death_rate_MHT <- 1 - (total_seeds_produced_MHT / total_seeds_possible_MHT)

total_starting_seeds_TLO <- nrow(tuebingenlowwateroneseedselectivedeathsdata)
total_surviving_adults_TLO <- sum(tuebingenlowwateroneseedselectivedeathsdata$Surviving_individuals_pot)
average_adult_death_rate_TLO <- 1 - (total_surviving_adults_TLO / total_starting_seeds_TLO)
proportion_deaths_selective_TLO <- total_selective_deaths_survival_TLO / total_starting_seeds_TLO

total_seeds_possible_TLO <- (total_surviving_adults_TLO * max_birth_rate_TLO)
total_seeds_produced_TLO <- sum(tuebingenlowwateroneseedselectivedeathsdata$Seeds_total_pot)
average_seed_death_rate_TLO <- 1 - (total_seeds_produced_TLO / total_seeds_possible_TLO)

total_starting_seeds_THO <- nrow(tuebingenhighwateroneseedselectivedeathsdata)
total_surviving_adults_THO <- sum(tuebingenhighwateroneseedselectivedeathsdata$Surviving_individuals_pot)
average_adult_death_rate_THO <- 1 - (total_surviving_adults_THO / total_starting_seeds_THO)
proportion_deaths_selective_THO <- total_selective_deaths_survival_THO / total_starting_seeds_THO

total_seeds_possible_THO <- (total_surviving_adults_THO * max_birth_rate_THO)
total_seeds_produced_THO <- sum(tuebingenhighwateroneseedselectivedeathsdata$Seeds_total_pot)
average_seed_death_rate_THO <- 1 - (total_seeds_produced_THO / total_seeds_possible_THO)

total_starting_seeds_TLT <- 30*nrow(tuebingenlowwatertwentyseedsselectivedeathsdata)
total_surviving_adults_TLT <- sum(tuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot)
average_adult_death_rate_TLT <- 1 - (total_surviving_adults_TLT / total_starting_seeds_TLT)
proportion_deaths_selective_TLT <- total_selective_deaths_survival_TLT / total_starting_seeds_TLT

total_seeds_possible_TLT <- (total_surviving_adults_TLT * max_birth_rate_TLT)
total_seeds_produced_TLT <- sum(tuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_total_pot)
average_seed_death_rate_TLT <- 1 - (total_seeds_produced_TLT / total_seeds_possible_TLT)

total_starting_seeds_THT <- 30*nrow(tuebingenhighwatertwentyseedsselectivedeathsdata)
total_surviving_adults_THT <- sum(tuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot)
average_adult_death_rate_THT <- 1 - (total_surviving_adults_THT / total_starting_seeds_THT)
proportion_deaths_selective_THT <- total_selective_deaths_survival_THT / total_starting_seeds_THT

total_seeds_possible_THT <- (total_surviving_adults_THT * max_birth_rate_THT)
total_seeds_produced_THT <- sum(tuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_total_pot)
average_seed_death_rate_THT <- 1 - (total_seeds_produced_THT / total_seeds_possible_THT)

#Summing up all of these to compare average seed selective deaths and average juvenile selective deaths
total_starting_seeds <- total_starting_seeds_MHO + total_starting_seeds_MHT + total_starting_seeds_MLO + total_starting_seeds_MLT + total_starting_seeds_THO + total_starting_seeds_THT + total_starting_seeds_TLO + total_starting_seeds_TLT
total_surviving_adults <- total_surviving_adults_MHO + total_surviving_adults_MHT + total_surviving_adults_MLO + total_surviving_adults_MLT + total_surviving_adults_THO + total_surviving_adults_THT + total_surviving_adults_TLO + total_surviving_adults_TLT
total_seeds_possible <- total_seeds_possible_MHO + total_seeds_possible_MHT + total_seeds_possible_MLO + total_seeds_possible_MLT + total_seeds_possible_THO + total_seeds_possible_THT + total_seeds_possible_TLO + total_seeds_possible_TLT
total_seeds_produced <- total_seeds_produced_MHO + total_seeds_produced_MHT + total_seeds_produced_MLO + total_seeds_produced_MLT + total_seeds_produced_THO + total_seeds_produced_THT + total_seeds_produced_TLO + total_seeds_produced_TLT
proportion_deaths_selective_seeds <- 1 - (total_seeds_produced / total_seeds_possible)
proportion_deaths_selective_juveniles <- 1 - ((total_surviving_adults + 29503.6) / total_starting_seeds) #the 29,503.6 is the number of non-selective deaths during the juvenile phase

#Comparing starting seeds to seeds produced at the end by survivors
#gives an idea of the reproductive excess available for selective deaths
#It's not real reproductive excess, since this estimate includes selective deaths
seedtoseedreproductiveexcessMLO <- (total_seeds_produced_MLO / total_starting_seeds_MLO) - 1
seedtoseedreproductiveexcessMLT <- (total_seeds_produced_MLT / total_starting_seeds_MLT) - 1
seedtoseedreproductiveexcessMHO <- (total_seeds_produced_MHO / total_starting_seeds_MHO) - 1
seedtoseedreproductiveexcessMHT <- (total_seeds_produced_MHT / total_starting_seeds_MHT) - 1
seedtoseedreproductiveexcessTLO <- (total_seeds_produced_TLO / total_starting_seeds_TLO) - 1
seedtoseedreproductiveexcessTLT <- (total_seeds_produced_TLT / total_starting_seeds_TLT) - 1
seedtoseedreproductiveexcessTHO <- (total_seeds_produced_THO / total_starting_seeds_THO) - 1
seedtoseedreproductiveexcessTHT <- (total_seeds_produced_THT / total_starting_seeds_THT) - 1

kMLO <- (total_seeds_produced_MLO / total_surviving_adults_MLO) - 1
kMLT <- (total_seeds_produced_MLT / total_surviving_adults_MLT) - 1
kMHO <- (total_seeds_produced_MHO / total_surviving_adults_MHO) - 1
kMHT <- (total_seeds_produced_MHT / total_surviving_adults_MHT) - 1
kTLO <- (total_seeds_produced_TLO / total_surviving_adults_TLO) - 1
kTLT <- (total_seeds_produced_TLT / total_surviving_adults_TLT) - 1
kTHO <- (total_seeds_produced_THO / total_surviving_adults_THO) - 1
kTHT <- (total_seeds_produced_THT / total_surviving_adults_THT) - 1

reproductiveexcessMLO <- ((total_seeds_produced_MLO / total_surviving_adults_MLO) - 1) * total_surviving_adults_MLO
reproductiveexcessMLT <- ((total_seeds_produced_MLT / total_surviving_adults_MLT) - 1) * total_surviving_adults_MLT

kMLO
kMHO
kTLO
kTHO
kMLT
kMHT
kTLT
kTHT
 seedtoseedreproductiveexcessMLO
reproductiveexcessMHO
reproductiveexcessTLO
reproductiveexcessTHO
reproductiveexcessMLT
reproductiveexcessMHT
reproductiveexcessTLT
reproductiveexcessTHT

total_surviving_adults
total_starting_seeds
proportion_deaths_selective_seeds
proportion_deaths_selective_juveniles

#I'm going to quickly plot a histogram of survival rates in the MLT condition.
#It has a clear pattern of a small handful of successful genotypes and many failed genotypes.
hist(1-aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate)
aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate
sum(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate < 0.2)

#Final death rate numbers
average_adult_death_rate_MLO
average_seed_death_rate_MLO
average_adult_death_rate_MHO
average_seed_death_rate_MHO
average_adult_death_rate_MLT
average_seed_death_rate_MLT
average_adult_death_rate_MHT
average_seed_death_rate_MHT
average_adult_death_rate_TLO
average_seed_death_rate_TLO
average_adult_death_rate_THO
average_seed_death_rate_THO
average_adult_death_rate_TLT
average_seed_death_rate_TLT
average_adult_death_rate_THT
average_seed_death_rate_THT

#Proportion of deaths selective, where different from death rate
proportion_deaths_selective_TLT
proportion_deaths_selective_THT

#Final selective death numbers
total_selective_deaths_survival_MLO
total_selective_deaths_births_MLO
total_selective_deaths_survival_MHO
total_selective_deaths_births_MHO
total_selective_deaths_survival_MLT
total_selective_deaths_births_MLT
total_selective_deaths_survival_MHT
total_selective_deaths_births_MHT
total_selective_deaths_survival_TLO
total_selective_deaths_births_TLO
total_selective_deaths_survival_THO
total_selective_deaths_births_THO
total_selective_deaths_survival_TLT
total_selective_deaths_births_TLT
total_selective_deaths_survival_THT
total_selective_deaths_births_THT

#Comparing total number of selective deaths of juveniles vs seeds
total_selective_deaths_survival <- total_selective_deaths_survival_MHO + total_selective_deaths_survival_MHT + total_selective_deaths_survival_MLO + total_selective_deaths_survival_MLT + total_selective_deaths_survival_THO + total_selective_deaths_survival_TLO + total_selective_deaths_survival_THT + total_selective_deaths_survival_TLT
total_selective_deaths_births <- total_selective_deaths_births_MHO + total_selective_deaths_births_MHT + total_selective_deaths_births_MLO + total_selective_deaths_births_MLT + total_selective_deaths_births_THO + total_selective_deaths_births_THT + total_selective_deaths_births_TLO + total_selective_deaths_births_TLT
total_selective_deaths_survival
total_selective_deaths_births

#I want to estimate the bias in the fittest genotype's death rate
#caused by choosing the lowest observed death rate to be our fittest
#So I'm going to simulate drawing deaths for each genotype from a binomial distribution
#with the observed means and number of plants
#I can then calculate a new 'observed' death rate for every genotype, choose the lowest
#and then compare to its 'real' death rate
#This will give me an idea of how much the process of choosing the lowest observed death rate
#biases towards genotypes who 'got lucky' and had a fitter year than they normally would

#Need two data frames, one with a row for each genotype, one with a row for each simulation run
simulationdataframeTHT <- data.frame(
  simulatedgenotypedeathrates = c(1:nrow(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata)),
  genotypeID = aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_id,
  resampledgenotypesurvivalrates = aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate
)
resultsdataframeTHT <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeTLT <- data.frame(
  simulatedgenotypedeathrates = c(1:nrow(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata)),
  genotypeID = aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_id,
  resampledgenotypesurvivalrates = aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate
)
resultsdataframeTLT <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeTHO <- data.frame(
  simulatedgenotypedeathrates = c(1:nrow(aggregatedtuebingenhighwateroneseedselectivedeathsdata)),
  genotypeID = aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_id,
  resampledgenotypesurvivalrates = aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate
)
resultsdataframeTHO <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeTLO <- data.frame(
  simulatedgenotypedeathrates = c(1:nrow(aggregatedtuebingenlowwateroneseedselectivedeathsdata)),
  genotypeID = aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_id,
  resampledgenotypesurvivalrates = aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate
)
resultsdataframeTLO <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeMHT <- data.frame(
  simulatedgenotypedeathrates = c(1:nrow(aggregatedmadridhighwatertwentyseedsselectivedeathsdata)),
  genotypeID = aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_id,
  resampledgenotypesurvivalrates = aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate
)
resultsdataframeMHT <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeMLT <- data.frame(
  simulatedgenotypedeathrates = c(1:nrow(aggregatedmadridlowwatertwentyseedsselectivedeathsdata)),
  genotypeID = aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_id,
  resampledgenotypesurvivalrates = aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate
)
resultsdataframeMLT <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeMHO <- data.frame(
  simulatedgenotypedeathrates = c(1:nrow(aggregatedmadridhighwateroneseedselectivedeathsdata)),
  genotypeID = aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_id,
  resampledgenotypesurvivalrates = aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_death_rate
)
resultsdataframeMHO <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframeMLO <- data.frame(
  simulatedgenotypedeathrates = c(1:nrow(aggregatedmadridlowwateroneseedselectivedeathsdata)),
  genotypeID = aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_id,
  resampledgenotypesurvivalrates = aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_death_rate
)
resultsdataframeMLO <- data.frame(
  vectorofmaximumsurvivalrates = c(1.0:10000.0),
  vectoroftruesurvivalrates = c(1.0:10000.0),
  vectorofwithingenotypevariances = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)
#These lines set any death rates that are 0 or 1 to a tiny fraction above or below 0 and 1
#This is because the R function to fit a beta distribution crashes on any 0 or 1 values.
aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate[aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate == 1.0] <- 1.0 - 1e-6
aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate[aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate == 0.0] <- 0.0 + 1e-6
aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate[aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate == 1.0] <- 1.0 - 1e-6
aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate[aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate == 0.0] <- 0.0 + 1e-6
aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate[aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate == 1.0] <- 1.0 - 1e-6
aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate[aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate == 0.0] <- 0.0 + 1e-6
aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate[aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate == 1.0] <- 1.0 - 1e-6
aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate[aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate == 0.0] <- 0.0 + 1e-6
aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate[aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate == 1.0] <- 1.0 - 1e-6
aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate[aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate == 0.0] <- 0.0 + 1e-6
aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate[aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate == 1.0] <- 1.0 - 1e-6
aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate[aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate == 0.0] <- 0.0 + 1e-6
aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_death_rate[aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_death_rate == 1.0] <- 1.0 - 1e-6
aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_death_rate[aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_death_rate == 0.0] <- 0.0 + 1e-6
aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_death_rate[aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_death_rate == 1.0] <- 1.0 - 1e-6
aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_death_rate[aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_death_rate == 0.0] <- 0.0 + 1e-6


#We also need to resample the mean death rates of each genotype, rather than use the observed values.
#We assume that the mean death rates of each genotype follow a beta distribution.
#Following lines fit a beta distribution to each environmental condition to be used in the simulation loop 
#Some of these lines give some warnings, but visual inspection of the resulting distributions look fine.
betafitforbetweengenotypevarianceTLT <- fitdistr(1 - aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate, "beta", list(shape1 = 1, shape2 = 1))
betafitforbetweengenotypevarianceTHT <- fitdistr(1 - aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate, "beta", list(shape1 = 1, shape2 = 1))
betafitforbetweengenotypevarianceTLO <- fitdistr(1 - aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate, "beta", list(shape1 = 2, shape2 = 2))
betafitforbetweengenotypevarianceTHO <- fitdistr(1 - aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate, "beta", list(shape1 = 0.1, shape2 = 1))
betafitforbetweengenotypevarianceMLT <- fitdistr(1 - aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate, "beta", list(shape1 = 1, shape2 = 1))
betafitforbetweengenotypevarianceMHT <- fitdistr(1 - aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate, "beta", list(shape1 = 1, shape2 = 1))
betafitforbetweengenotypevarianceMLO <- fitdistr(1 - aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_death_rate, "beta", list(shape1 = 2, shape2 = 2))
betafitforbetweengenotypevarianceMHO <- fitdistr(1 - aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_death_rate, "beta", list(shape1 = 0.1, shape2 = 1))

betafitforbetweengenotypevarianceTLT
betafitforbetweengenotypevarianceTHT
betafitforbetweengenotypevarianceTLO
betafitforbetweengenotypevarianceTHO
betafitforbetweengenotypevarianceMLT
betafitforbetweengenotypevarianceMHT
betafitforbetweengenotypevarianceMLO
betafitforbetweengenotypevarianceMHO

hist(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate)
hist(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate)

#Following lines are the main simulation loop, which is run 10,000 times.
#Each time, 'true' survival rates are redrawn from the beta distribution estimated earlier
#Then each genotype randomly draws a number of survivors from a binomial with their survival rate.
#Additional lines at the end of the loop are for diagnostic purposes.
for (i in 1:10000) {
  simulationdataframeTHT$resampledgenotypesurvivalrates <- rbeta(517, 3.997, 5.363)
  for(j in 1:517) {
    simulationdataframeTHT$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[j])*30, simulationdataframeTHT$resampledgenotypesurvivalrates[j]) / (aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[j] * 30))
  }
  resultsdataframeTHT$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeTHT$simulatedgenotypesurvivalrates)
  resultsdataframeTHT$vectoroftruesurvivalrates[i] <- max(simulationdataframeTHT$resampledgenotypesurvivalrates)
  resultsdataframeTHT$vectorofbestgenotypeID[i] <- simulationdataframeTHT$genotypeID[which.max(simulationdataframeTHT$simulatedgenotypesurvivalrates)]
  resultsdataframeTHT$vectorofwithingenotypevariances[i] <- resultsdataframeTHT$vectoroftruesurvivalrates[i]*(1 - resultsdataframeTHT$vectoroftruesurvivalrates[i]) / (aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[which(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_id == resultsdataframeTHT$vectorofbestgenotypeID[i], arr.ind = FALSE)]*30)
}
resultsdataframeTHT$vectorofmaximumsurvivalrates
resultsdataframeTHT$vectoroftruesurvivalrates
resultsdataframeTHT$vectorofbestgenotypeID
hist(resultsdataframeTHT$vectorofbestgenotypeID)
hist(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_id)
hist(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates)
resultsdataframeTHT$vectorofwithingenotypevariances
mean(resultsdataframeTHT$vectorofwithingenotypevariances)
mean(resultsdataframeTHT$vectorofmaximumsurvivalrates)
mean(resultsdataframeTHT$vectoroftruesurvivalrates)

for (i in 1:10000) {
  simulationdataframeTLT$resampledgenotypesurvivalrates <- rbeta(517, 0.285, 1.72)
  for(j in 1:517) {
    simulationdataframeTLT$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[j])*30, simulationdataframeTLT$resampledgenotypesurvivalrates[j]) / (aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[j] * 30))
  }
  resultsdataframeTLT$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeTLT$simulatedgenotypesurvivalrates)
  resultsdataframeTLT$vectoroftruesurvivalrates[i] <- max(simulationdataframeTLT$resampledgenotypesurvivalrates)
  resultsdataframeTLT$vectorofbestgenotypeID[i] <- simulationdataframeTLT$genotypeID[which.max(simulationdataframeTLT$resampledgenotypesurvivalrates)]
  resultsdataframeTLT$vectorofwithingenotypevariances[i] <- resultsdataframeTLT$vectoroftruesurvivalrates[i]*(1 - resultsdataframeTLT$vectoroftruesurvivalrates[i]) / (aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[which(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_id == resultsdataframeTLT$vectorofbestgenotypeID[i], arr.ind = FALSE)]*30)
}
resultsdataframeTLT$vectorofmaximumsurvivalrates
resultsdataframeTLT$vectoroftruesurvivalrates
resultsdataframeTLT$vectorofbestgenotypeID
hist(resultsdataframeTLT$vectorofbestgenotypeID)
hist(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_id)
hist(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates)
resultsdataframeTLT$vectorofwithingenotypevariances
mean(resultsdataframeTLT$vectorofwithingenotypevariances)
mean(resultsdataframeTLT$vectorofmaximumsurvivalrates)
mean(resultsdataframeTLT$vectoroftruesurvivalrates)

for (i in 1:10000) {
  simulationdataframeTLO$resampledgenotypesurvivalrates <- rbeta(517, 1.13, 1.67)
  for(j in 1:517) {
    simulationdataframeTLO$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[j]), simulationdataframeTLO$resampledgenotypesurvivalrates[j]) / (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[j]))
  }
  resultsdataframeTLO$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeTLO$simulatedgenotypesurvivalrates)
  resultsdataframeTLO$vectoroftruesurvivalrates[i] <- max(simulationdataframeTLO$resampledgenotypesurvivalrates)
  resultsdataframeTLO$vectorofbestgenotypeID[i] <- simulationdataframeTLO$genotypeID[which.max(simulationdataframeTLO$resampledgenotypesurvivalrates)]
  resultsdataframeTLO$vectorofwithingenotypevariances[i] <- resultsdataframeTLO$vectoroftruesurvivalrates[i]*(1 - resultsdataframeTLO$vectoroftruesurvivalrates[i]) / (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[which(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_id == resultsdataframeTLO$vectorofbestgenotypeID[i], arr.ind = FALSE)])
}
resultsdataframeTLO$vectorofmaximumsurvivalrates
resultsdataframeTLO$vectoroftruesurvivalrates
resultsdataframeTLO$vectorofbestgenotypeID
hist(resultsdataframeTLO$vectorofbestgenotypeID)
hist(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_id)
hist(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates)
resultsdataframeTLO$vectorofwithingenotypevariances
mean(resultsdataframeTLO$vectorofwithingenotypevariances)

for (i in 1:10000) {
  simulationdataframeTHO$resampledgenotypesurvivalrates <- rbeta(517, 1.595, 0.155)
  for(j in 1:517) {
    simulationdataframeTHO$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[j]), simulationdataframeTHO$resampledgenotypesurvivalrates[j]) / (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[j]))
  }
  resultsdataframeTHO$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeTHO$simulatedgenotypesurvivalrates)
  resultsdataframeTHO$vectoroftruesurvivalrates[i] <- max(simulationdataframeTHO$resampledgenotypesurvivalrates)
  resultsdataframeTHO$vectorofbestgenotypeID[i] <- simulationdataframeTHO$genotypeID[which.max(simulationdataframeTHO$resampledgenotypesurvivalrates)]
  resultsdataframeTHO$vectorofwithingenotypevariances[i] <- resultsdataframeTHO$vectoroftruesurvivalrates[i]*(1 - resultsdataframeTHO$vectoroftruesurvivalrates[i]) / (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[which(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_id == resultsdataframeTHO$vectorofbestgenotypeID[i], arr.ind = FALSE)])
}
resultsdataframeTHO$vectorofmaximumsurvivalrates
resultsdataframeTHO$vectoroftruesurvivalrates
resultsdataframeTHO$vectorofbestgenotypeID
hist(resultsdataframeTHO$vectorofbestgenotypeID)
hist(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_id)
hist(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates)
resultsdataframeTHO$vectorofwithingenotypevariances
mean(resultsdataframeTHO$vectorofwithingenotypevariances)

for (i in 1:10000) {
  simulationdataframeMHT$resampledgenotypesurvivalrates <- rbeta(517, 1.594, 2.031)
  for(j in 1:517) {
    simulationdataframeMHT$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[j])*30, simulationdataframeMHT$resampledgenotypesurvivalrates[j]) / (aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[j] * 30))
  }
  resultsdataframeMHT$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeMHT$simulatedgenotypesurvivalrates)
  resultsdataframeMHT$vectoroftruesurvivalrates[i] <- max(simulationdataframeMHT$resampledgenotypesurvivalrates)
  resultsdataframeMHT$vectorofbestgenotypeID[i] <- simulationdataframeMHT$genotypeID[which.max(simulationdataframeMHT$resampledgenotypesurvivalrates)]
  resultsdataframeMHT$vectorofwithingenotypevariances[i] <- resultsdataframeMHT$vectoroftruesurvivalrates[i]*(1 - resultsdataframeMHT$vectoroftruesurvivalrates[i]) / (aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[which(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_id == resultsdataframeMHT$vectorofbestgenotypeID[i], arr.ind = FALSE)]*30)
}
resultsdataframeMHT$vectorofmaximumsurvivalrates
resultsdataframeMHT$vectoroftruesurvivalrates
resultsdataframeMHT$vectorofbestgenotypeID
hist(resultsdataframeMHT$vectorofbestgenotypeID)
hist(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_id)
hist(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates)
resultsdataframeMHT$vectorofwithingenotypevariances
mean(resultsdataframeMHT$vectorofwithingenotypevariances)

#MLT is the only environmental condition to have only 516 genotypes, hence some changes.
#Best coding practice would have been to use nrows on each dataset instead, but I'm lazy.
for (i in 1:10000) {
  simulationdataframeMLT$resampledgenotypesurvivalrates <- rbeta(516, 0.111, 1.894)
  for(j in 1:516) {
    simulationdataframeMLT$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[j])*30, simulationdataframeMLT$resampledgenotypesurvivalrates[j]) / (aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[j] * 30))
  }
  resultsdataframeMLT$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeMLT$simulatedgenotypesurvivalrates)
  resultsdataframeMLT$vectoroftruesurvivalrates[i] <- max(simulationdataframeMLT$resampledgenotypesurvivalrates)
  resultsdataframeMLT$vectorofbestgenotypeID[i] <- simulationdataframeMLT$genotypeID[which.max(simulationdataframeMLT$resampledgenotypesurvivalrates)]
  resultsdataframeMLT$vectorofwithingenotypevariances[i] <- resultsdataframeMLT$vectoroftruesurvivalrates[i]*(1 - resultsdataframeMLT$vectoroftruesurvivalrates[i]) / (aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[which(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_id == resultsdataframeMLT$vectorofbestgenotypeID[i], arr.ind = FALSE)]*30)
}
resultsdataframeMLT$vectorofmaximumsurvivalrates
resultsdataframeMLT$vectoroftruesurvivalrates
resultsdataframeMLT$vectorofbestgenotypeID
hist(resultsdataframeMLT$vectorofbestgenotypeID)
hist(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_id)
hist(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates)
resultsdataframeMLT$vectorofwithingenotypevariances
mean(resultsdataframeMLT$vectorofwithingenotypevariances)

for (i in 1:10000) {
  simulationdataframeMLO$resampledgenotypesurvivalrates <- rbeta(517, 0.2697, 0.6757)
  for(j in 1:517) {
    simulationdataframeMLO$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[j]), simulationdataframeMLO$resampledgenotypesurvivalrates[j]) / (aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[j]))
  }
  resultsdataframeMLO$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeMLO$simulatedgenotypesurvivalrates)
  resultsdataframeMLO$vectoroftruesurvivalrates[i] <- max(simulationdataframeMLO$resampledgenotypesurvivalrates)
  resultsdataframeMLO$vectorofbestgenotypeID[i] <- simulationdataframeMLO$genotypeID[which.max(simulationdataframeMLO$resampledgenotypesurvivalrates)]
  resultsdataframeMLO$vectorofwithingenotypevariances[i] <- resultsdataframeMLO$vectoroftruesurvivalrates[i]*(1 - resultsdataframeMLO$vectoroftruesurvivalrates[i]) / (aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[which(aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_id == resultsdataframeMLO$vectorofbestgenotypeID[i], arr.ind = FALSE)])
}
resultsdataframeMLO$vectorofmaximumsurvivalrates
resultsdataframeMLO$vectoroftruesurvivalrates
resultsdataframeMLO$vectorofbestgenotypeID
hist(resultsdataframeMLO$vectorofbestgenotypeID)
hist(aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_id)
hist(aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates)
resultsdataframeMLO$vectorofwithingenotypevariances
mean(resultsdataframeMLO$vectorofwithingenotypevariances)

for (i in 1:10000) {
  simulationdataframeMHO$resampledgenotypesurvivalrates <- rbeta(517, 1.3345, 0.1195)
  for(j in 1:517) {
    simulationdataframeMHO$simulatedgenotypesurvivalrates[j] <- (rbinom(1, (aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[j]), simulationdataframeMHO$resampledgenotypesurvivalrates[j]) / (aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[j]))
  }
  resultsdataframeMHO$vectorofmaximumsurvivalrates[i] <- max(simulationdataframeMHO$simulatedgenotypesurvivalrates)
  resultsdataframeMHO$vectoroftruesurvivalrates[i] <- max(simulationdataframeMHO$resampledgenotypesurvivalrates)
  resultsdataframeMHO$vectorofbestgenotypeID[i] <- simulationdataframeMHO$genotypeID[which.max(simulationdataframeMHO$resampledgenotypesurvivalrates)]
  resultsdataframeMHO$vectorofwithingenotypevariances[i] <- resultsdataframeMHO$vectoroftruesurvivalrates[i]*(1 - resultsdataframeMHO$vectoroftruesurvivalrates[i]) / (aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[which(aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_id == resultsdataframeMHO$vectorofbestgenotypeID[i], arr.ind = FALSE)])
}
resultsdataframeMHO$vectorofmaximumsurvivalrates
resultsdataframeMHO$vectoroftruesurvivalrates
resultsdataframeMHO$vectorofbestgenotypeID
hist(resultsdataframeMHO$vectorofbestgenotypeID)
hist(aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_id)
hist(aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates)
resultsdataframeMHO$vectorofwithingenotypevariances
mean(resultsdataframeMHO$vectorofwithingenotypevariances)

#Records the true death rate minus observed death rate and averages to see the bias.
biasvectorTHT <- c(1:10000)
for (i in 1:10000) {
  biasvectorTHT[i] = (resultsdataframeTHT$vectorofmaximumsurvivalrates[i] - resultsdataframeTHT$vectoroftruesurvivalrates[i])
}
mean(biasvectorTHT)

biasvectorTLT <- c(1:10000)
for (i in 1:10000) {
  biasvectorTLT[i] = (resultsdataframeTLT$vectorofmaximumsurvivalrates[i] - resultsdataframeTLT$vectoroftruesurvivalrates[i])
}
mean(biasvectorTLT)

biasvectorTHO <- c(1:10000)
for (i in 1:10000) {
  biasvectorTHO[i] = (resultsdataframeTHO$vectorofmaximumsurvivalrates[i] - resultsdataframeTHO$vectoroftruesurvivalrates[i])
}
mean(biasvectorTHO)

biasvectorTLO <- c(1:10000)
for (i in 1:10000) {
  biasvectorTLO[i] = (resultsdataframeTLO$vectorofmaximumsurvivalrates[i] - resultsdataframeTLO$vectoroftruesurvivalrates[i])
}
mean(biasvectorTLO)

biasvectorMHT <- c(1:10000)
for (i in 1:10000) {
  biasvectorMHT[i] = (resultsdataframeMHT$vectorofmaximumsurvivalrates[i] - resultsdataframeMHT$vectoroftruesurvivalrates[i])
}
mean(biasvectorMHT)

biasvectorMLT <- c(1:10000)
for (i in 1:10000) {
  biasvectorMLT[i] = (resultsdataframeMLT$vectorofmaximumsurvivalrates[i] - resultsdataframeMLT$vectoroftruesurvivalrates[i])
}
mean(biasvectorMLT)

biasvectorMHO <- c(1:10000)
for (i in 1:10000) {
  biasvectorMHO[i] = (resultsdataframeMHO$vectorofmaximumsurvivalrates[i] - resultsdataframeMHO$vectoroftruesurvivalrates[i])
}
mean(biasvectorMHO)

biasvectorMLO <- c(1:10000)
for (i in 1:10000) {
  biasvectorMLO[i] = (resultsdataframeMLO$vectorofmaximumsurvivalrates[i] - resultsdataframeMLO$vectoroftruesurvivalrates[i])
}
mean(biasvectorMLO)

#Ideally, the variance in fitness from the same genotype being simulated repeatedly should be
#much less than the overall variance between genotypes observed in the original study.
var(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate)
var(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate)
var(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_death_rate)
var(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_death_rate)
var(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_death_rate)
var(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_death_rate)
var(aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_death_rate)
var(aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_death_rate)

#Following lines use the estimate of bias from above to adjust number of selective deaths
#A line at the end shows the unadjusted number of selective deaths for comparison
aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[1] <- 0
for (i in 1:nrow(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata)) {
  aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- ((30*aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[i]) - aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) - (30*aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Replicates[i]*(min_death_rate_THT+0.0131))
  if (aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] < 0.0) {
    aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- 0.0
  }
}
adjusted_total_selective_deaths_survival_THT <- sum(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival)
adjusted_total_selective_deaths_survival_THT
total_selective_deaths_survival_THT

aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[1] <- 0
for (i in 1:nrow(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata)) {
  aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- ((30*aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[i]) - aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) - (30*aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Replicates[i]*(min_death_rate_TLT+0.0054))
  if (aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] < 0.0) {
    aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- 0.0
  }
}
adjusted_total_selective_deaths_survival_TLT <- sum(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival)
adjusted_total_selective_deaths_survival_TLT
total_selective_deaths_survival_TLT

aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[1] <- 0
for (i in 1:nrow(aggregatedtuebingenhighwateroneseedselectivedeathsdata)) {
  aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[i] - aggregatedtuebingenhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]) - (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Replicates[i]*(min_death_rate_THO+2.16e-15))
  if (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] < 0.0) {
    aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- 0.0
  }
}
adjusted_total_selective_deaths_survival_THO <- sum(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival)
adjusted_total_selective_deaths_survival_THO
total_selective_deaths_survival_THO

aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[1] <- 0
for (i in 1:nrow(aggregatedtuebingenlowwateroneseedselectivedeathsdata)) {
  aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[i] - aggregatedtuebingenlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]) - (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Replicates[i]*(min_death_rate_TLO+0.0191))
  if (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] < 0.0) {
    aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- 0.0
  }
}
adjusted_total_selective_deaths_survival_TLO <- sum(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival)
adjusted_total_selective_deaths_survival_TLO
total_selective_deaths_survival_TLO

aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[1] <- 0
for (i in 1:nrow(aggregatedmadridhighwatertwentyseedsselectivedeathsdata)) {
  aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- ((30*aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[i]) - aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) - (30*aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Replicates[i]*(min_death_rate_MHT+0.00739))
  if (aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] < 0.0) {
    aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- 0.0
  }
}
adjusted_total_selective_deaths_survival_MHT <- sum(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival)
adjusted_total_selective_deaths_survival_MHT
total_selective_deaths_survival_MHT

aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[1] <- 0
for (i in 1:nrow(aggregatedmadridlowwatertwentyseedsselectivedeathsdata)) {
  aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- ((30*aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[i]) - aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]) - (30*aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Replicates[i]*(min_death_rate_MLT+0.00758))
  if (aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] < 0.0) {
    aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- 0.0
  }
}
adjusted_total_selective_deaths_survival_MLT <- sum(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_survival)
adjusted_total_selective_deaths_survival_MLT
total_selective_deaths_survival_MLT

aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[1] <- 0
for (i in 1:nrow(aggregatedmadridhighwateroneseedselectivedeathsdata)) {
  aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- (aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[i] - aggregatedmadridhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]) - (aggregatedmadridhighwateroneseedselectivedeathsdata$Replicates[i]*(min_death_rate_MHO+0.0149))
  if (aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] < 0.0) {
    aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- 0.0
  }
}
adjusted_total_selective_deaths_survival_MHO <- sum(aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival)
adjusted_total_selective_deaths_survival_MHO
total_selective_deaths_survival_MHO

aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[1] <- 0
for (i in 1:nrow(aggregatedmadridlowwateroneseedselectivedeathsdata)) {
  aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- (aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[i] - aggregatedmadridlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]) - (aggregatedmadridlowwateroneseedselectivedeathsdata$Replicates[i]*(min_death_rate_MLO+0.000609))
  if (aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] < 0.0) {
    aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival[i] <- 0.0
  }
}
adjusted_total_selective_deaths_survival_MLO <- sum(aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_survival)
adjusted_total_selective_deaths_survival_MLO
total_selective_deaths_survival_MLO

#Now, using the above totals, I calculate an adjusted proportion of deaths which are selective.
#Just for fun, I include a line with the unadjusted proportion (repeated from earlier) for comparison.
adjusted_proportion_deaths_selective_THT <- adjusted_total_selective_deaths_survival_THT / total_starting_seeds_THT
adjusted_proportion_deaths_selective_TLT <- adjusted_total_selective_deaths_survival_TLT / total_starting_seeds_TLT
adjusted_proportion_deaths_selective_THO <- adjusted_total_selective_deaths_survival_THO / total_starting_seeds_THO
adjusted_proportion_deaths_selective_TLO <- adjusted_total_selective_deaths_survival_TLO / total_starting_seeds_TLO
adjusted_proportion_deaths_selective_MHT <- adjusted_total_selective_deaths_survival_MHT / total_starting_seeds_MHT
adjusted_proportion_deaths_selective_MLT <- adjusted_total_selective_deaths_survival_MLT / total_starting_seeds_MLT
adjusted_proportion_deaths_selective_MHO <- adjusted_total_selective_deaths_survival_MHO / total_starting_seeds_MHO
adjusted_proportion_deaths_selective_MLO <- adjusted_total_selective_deaths_survival_MLO / total_starting_seeds_MLO

proportion_deaths_selective_THT
adjusted_proportion_deaths_selective_THT
proportion_deaths_selective_TLT
adjusted_proportion_deaths_selective_TLT
proportion_deaths_selective_THO
adjusted_proportion_deaths_selective_THO
proportion_deaths_selective_TLO
adjusted_proportion_deaths_selective_TLO
proportion_deaths_selective_MHT
adjusted_proportion_deaths_selective_MHT
proportion_deaths_selective_MLT
adjusted_proportion_deaths_selective_MLT
proportion_deaths_selective_MHO
adjusted_proportion_deaths_selective_MHO
proportion_deaths_selective_MLO
adjusted_proportion_deaths_selective_MLO

adjusted_total_proportion_deaths_selective_oneseed <- (adjusted_total_selective_deaths_survival_THO + adjusted_total_selective_deaths_survival_TLO + adjusted_total_selective_deaths_survival_MHO + adjusted_total_selective_deaths_survival_MLO) / (total_starting_seeds_THO + total_starting_seeds_TLO + total_starting_seeds_MHO + total_starting_seeds_MLO)
adjusted_total_proportion_deaths_selective_thirtyseeds <- (adjusted_total_selective_deaths_survival_MHT + adjusted_total_selective_deaths_survival_MLT + adjusted_total_selective_deaths_survival_THT + adjusted_total_selective_deaths_survival_TLT) / (total_starting_seeds_MHT + total_starting_seeds_MLT + total_starting_seeds_THT + total_starting_seeds_TLT)
adjusted_total_proportion_deaths_selective_oneseed
adjusted_total_proportion_deaths_selective_thirtyseeds

#Next section of code is for adjusting the selective deaths that occur during seed production.


#Following code removes genotypes with no survivors from the calculations of selective deaths on births.
#All the following calculations will be done on only the set of genotypes with at least some survivors.
TLTdfnozeroesinbirthrate <- aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata[which(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind != 0.0), names(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata) %in% c("Genotype_id", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Replicates")]
THTdfnozeroesinbirthrate <- aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata[which(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind != 0.0), names(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata) %in% c("Genotype_id", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Replicates")]
TLOdfnozeroesinbirthrate <- aggregatedtuebingenlowwateroneseedselectivedeathsdata[which(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Seeds_by_ind != 0.0), names(aggregatedtuebingenlowwateroneseedselectivedeathsdata) %in% c("Genotype_id", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Replicates")]
THOdfnozeroesinbirthrate <- aggregatedtuebingenhighwateroneseedselectivedeathsdata[which(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Seeds_by_ind != 0.0), names(aggregatedtuebingenhighwateroneseedselectivedeathsdata) %in% c("Genotype_id", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Replicates")]
MLTdfnozeroesinbirthrate <- aggregatedmadridlowwatertwentyseedsselectivedeathsdata[which(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind != 0.0), names(aggregatedmadridlowwatertwentyseedsselectivedeathsdata) %in% c("Genotype_id", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Replicates")]
MHTdfnozeroesinbirthrate <- aggregatedmadridhighwatertwentyseedsselectivedeathsdata[which(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind != 0.0), names(aggregatedmadridhighwatertwentyseedsselectivedeathsdata) %in% c("Genotype_id", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Replicates")]
MLOdfnozeroesinbirthrate <- aggregatedmadridlowwateroneseedselectivedeathsdata[which(aggregatedmadridlowwateroneseedselectivedeathsdata$Seeds_by_ind != 0.0), names(aggregatedmadridlowwateroneseedselectivedeathsdata) %in% c("Genotype_id", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Replicates")]
MHOdfnozeroesinbirthrate <- aggregatedmadridhighwateroneseedselectivedeathsdata[which(aggregatedmadridhighwateroneseedselectivedeathsdata$Seeds_by_ind != 0.0), names(aggregatedmadridhighwateroneseedselectivedeathsdata) %in% c("Genotype_id", "Surviving_individuals_pot", "Seeds_total_pot", "Seeds_by_ind", "Replicates")]


#For births, we'll be resampling genotype 'true' birth rates from a normal or log-normal,
#then resampling seeds per individual from a poisson with mean of the 'true' birth rate.
#However, there's some concern that the within-genotype variance in seed number might be
#on the same order of magnitude as the among-genotype variance in fitness,
#meaning that fitting a normal/log-normal directly to observed genotype birth rates
#would overestimate the real among-genotype variance.
#So instead I resample genotype birth rates from distributions with LESS variance than observed,
#draw from poissons for each genotype, and calculate variance on the resulting 'observations'.
#I'll then plot the variances, and infer that the starting distribution which results
#in a variance close to the observed variance approximates the real among-genotype variance.
THTbirthvariancenozeroes <- var(THTdfnozeroesinbirthrate$Seeds_by_ind)
THTbirthstdevnozeroes <- sqrt(THTbirthvariancenozeroes)
THTbirthmeannozeroes <- mean(THTdfnozeroesinbirthrate$Seeds_by_ind)
THTbirthmeannozeroes
THTbirthvariancenozeroes
THTbirthstdevnozeroes
THTstdevtotry <- c(201.0:270.0)
THTresultingvariances <- c(1.0:70.0)
THTaddingwithingenotypevariance <- c(1:nrow(THTdfnozeroesinbirthrate))
for (i in THTstdevtotry) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rnorm(nrow(THTdfnozeroesinbirthrate), THTbirthmeannozeroes, i)
    for (j in 1:nrow(THTdfnozeroesinbirthrate)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(THTdfnozeroesinbirthrate)) {
    THTaddingwithingenotypevariance[n] <- rpois(1, resampledamonggenotypetruebirths[n])
  }
  THTresultingvariances[i-200] <- var(THTaddingwithingenotypevariance)
}
THTresultingvariances

THTregressionvariancebystdev <- lm(THTresultingvariances ~ THTstdevtotry)
THTregressionvariancebystdev

#70694 = (438.9 * correctstdev) - 46945
(70694 + 46945) / 438.9
THTbirthstdevnozeroes
#In short, this analysis predicts that fitting a normal distribution to the observed
#distribution of genotypes (which includes within-genotype variance) and then adding 
#within-genotype variance in again by drawing from a Poisson distribution has NO effect
#on the resulting variance between genotypes. The standard deviation that produces the
#same resulting variance as was observed in the experiment is slightly LARGER than 
#the observed standard deviation. Going forward, I'll use the observed variance for THT.


TLTbirthvariancenozeroes <- var(TLTdfnozeroesinbirthrate$Seeds_by_ind)
TLTbirthvariancenozeroes
TLTbirthlogvariancenozeroes <- var(log(TLTdfnozeroesinbirthrate$Seeds_by_ind))
TLTbirthlogstdevnozeroes <- sqrt(TLTbirthlogvariancenozeroes)
TLTbirthlogmeannozeroes <- mean(log(TLTdfnozeroesinbirthrate$Seeds_by_ind))
TLTbirthlogmeannozeroes
TLTbirthlogvariancenozeroes
TLTbirthlogstdevnozeroes
TLTstdevtotry <- seq(0.75, 0.99, 0.005)
TLTresultingvariances <- c(1.0:49.0)
TLTaddingwithingenotypevariance <- c(1:nrow(TLTdfnozeroesinbirthrate))
for (i in 1:49) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rlnorm(nrow(TLTdfnozeroesinbirthrate), TLTbirthlogmeannozeroes, TLTstdevtotry[i])
    for (j in 1:nrow(TLTdfnozeroesinbirthrate)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(TLTdfnozeroesinbirthrate)) {
    TLTaddingwithingenotypevariance[n] <- rpois(1, resampledamonggenotypetruebirths[n])
  }
  TLTresultingvariances[i] <- var(TLTaddingwithingenotypevariance)
}
TLTresultingvariances

TLTregressionvariancebystdev <- lm(TLTresultingvariances ~ TLTstdevtotry)
TLTregressionvariancebystdev

#812420.5 = (1790981 * correctstdev) - 1194366
(812420.5 + 1194366) / 1790981

#This analysis predicts that fitting a log-normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) not only doesn't overestimate variance between genotypes, it actually noticeably
#UNDERESTIMATES it (probably because the observed distribution is slightly more spread out than a log-normal)
#Going forward, I'll use the observed variance for TLT.

THObirthvariancenozeroes <- var(THOdfnozeroesinbirthrate$Seeds_by_ind)
THObirthstdevnozeroes <- sqrt(THObirthvariancenozeroes)
THObirthmeannozeroes <- mean(THOdfnozeroesinbirthrate$Seeds_by_ind)
THObirthmeannozeroes
THObirthvariancenozeroes
THObirthstdevnozeroes
THOstdevtotry <- seq(2001.0, 2130.0, 1.0)
THOresultingvariances <- c(1.0:130.0)
THOaddingwithingenotypevariance <- c(1:nrow(THOdfnozeroesinbirthrate))
for (i in THOstdevtotry) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rnorm(nrow(THOdfnozeroesinbirthrate), THObirthmeannozeroes, i)
    for (j in 1:nrow(THOdfnozeroesinbirthrate)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(THOdfnozeroesinbirthrate)) {
    THOaddingwithingenotypevariance[n] <- rpois(1, resampledamonggenotypetruebirths[n])
  }
  THOresultingvariances[i-2000] <- var(THOaddingwithingenotypevariance)
}
THOresultingvariances

THOregressionvariancebystdev <- lm(THOresultingvariances ~ THOstdevtotry)
THOregressionvariancebystdev

#4492312 = (4350 * correctstdev) - 4689168
(4492312 + 4689168) / 4350
THObirthstdevnozeroes

#This analysis predicts that fitting a normal distribution to the observed distribution of
#genotype birth rates overestimates the correct standard deviation in birth rates between genotypes
#by about 9, out of an average standard deviation of over 2,000. This isn't large enough to 
#meaningfully affect results, so going forward I'll use the observed standard deviation for THO.

TLObirthvariancenozeroes <- var(TLOdfnozeroesinbirthrate$Seeds_by_ind)
TLObirthstdevnozeroes <- sqrt(TLObirthvariancenozeroes)
TLObirthmeannozeroes <- mean(TLOdfnozeroesinbirthrate$Seeds_by_ind)
TLObirthmeannozeroes
TLObirthvariancenozeroes
TLObirthstdevnozeroes
TLOstdevtotry <- seq(1201.0, 1375.0, 1.0)
TLOresultingvariances <- c(1.0:175.0)
TLOaddingwithingenotypevariance <- c(1:nrow(TLOdfnozeroesinbirthrate))
for (i in TLOstdevtotry) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rnorm(nrow(TLOdfnozeroesinbirthrate), TLObirthmeannozeroes, i)
    for (j in 1:nrow(TLOdfnozeroesinbirthrate)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(TLOdfnozeroesinbirthrate)) {
    TLOaddingwithingenotypevariance[n] <- rpois(1, resampledamonggenotypetruebirths[n])
  }
  TLOresultingvariances[i-1200] <- var(TLOaddingwithingenotypevariance)
}
TLOresultingvariances

TLOregressionvariancebystdev <- lm(TLOresultingvariances ~ TLOstdevtotry)
TLOregressionvariancebystdev

#1837291 = (2717 * correctstdev) - 1840052
(1837291 + 1840052) / 2717
TLObirthstdevnozeroes

#This analysis predicts that fitting a normal distribution to the observed distribution of 
#birth rates between genotypes and then adding in within-genotype variance has next to no effect
#on the resulting standard deviation. Going forward, I'll use the observed standard deviation for TLO.

MHTbirthvariancenozeroes <- var(MHTdfnozeroesinbirthrate$Seeds_by_ind)
MHTbirthvariancenozeroes
MHTbirthlogvariancenozeroes <- var(log(MHTdfnozeroesinbirthrate$Seeds_by_ind))
MHTbirthlogstdevnozeroes <- sqrt(MHTbirthlogvariancenozeroes)
MHTbirthlogmeannozeroes <- mean(log(MHTdfnozeroesinbirthrate$Seeds_by_ind))
MHTbirthlogmeannozeroes
MHTbirthlogvariancenozeroes
MHTbirthlogstdevnozeroes
MHTstdevtotry <- seq(0.36, 0.75, 0.005)
MHTresultingvariances <- c(1.0:79.0)
MHTaddingwithingenotypevariance <- c(1:nrow(MHTdfnozeroesinbirthrate))
for (i in 1:79) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rlnorm(nrow(MHTdfnozeroesinbirthrate), MHTbirthlogmeannozeroes, MHTstdevtotry[i])
    for (j in 1:nrow(MHTdfnozeroesinbirthrate)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(MHTdfnozeroesinbirthrate)) {
    MHTaddingwithingenotypevariance[n] <- rpois(1, resampledamonggenotypetruebirths[n])
  }
  MHTresultingvariances[i] <- var(MHTaddingwithingenotypevariance)
}
MHTresultingvariances

MHTregressionvariancebystdev <- lm(MHTresultingvariances ~ MHTstdevtotry)
MHTregressionvariancebystdev

#14460181 = (10372387 * correctstdev) - 3633876
(14460181 + 3633876) / 10372387

#This analysis predicts that fitting a log-normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) not only doesn't overestimate variance between genotypes, it seriously
#UNDERESTIMATES it (probably because the distribution of log(birth rate) has a long tail)
#Going forward, I'll use the observed variance for MHT.

MHObirthvariancenozeroes <- var(MHOdfnozeroesinbirthrate$Seeds_by_ind)
MHObirthstdevnozeroes <- sqrt(MHObirthvariancenozeroes)
MHObirthmeannozeroes <- mean(MHOdfnozeroesinbirthrate$Seeds_by_ind)
MHObirthmeannozeroes
MHObirthvariancenozeroes
MHObirthstdevnozeroes
MHOstdevtotry <- seq(3201.0, 3550.0, 1.0)
MHOresultingvariances <- c(1.0:350.0)
MHOaddingwithingenotypevariance <- c(1:nrow(MHOdfnozeroesinbirthrate))
for (i in MHOstdevtotry) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rnorm(nrow(MHOdfnozeroesinbirthrate), MHObirthmeannozeroes, i)
    for (j in 1:nrow(MHOdfnozeroesinbirthrate)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(MHOdfnozeroesinbirthrate)) {
    MHOaddingwithingenotypevariance[n] <- rpois(1, resampledamonggenotypetruebirths[n])
  }
  MHOresultingvariances[i-3200] <- var(MHOaddingwithingenotypevariance)
}
MHOresultingvariances

MHOregressionvariancebystdev <- lm(MHOresultingvariances ~ MHOstdevtotry)
MHOregressionvariancebystdev

#12215481 = (6289 * correctstdev) - 9838179
(12215481 + 9838179) / 6289
MHObirthstdevnozeroes

#This analysis predicts that fitting a normal distribution to the observed distribution of
#genotype birth rates overestimates the correct standard deviation in birth rates between genotypes
#by about 10, out of an average standard deviation of over 3,000. This isn't large enough to 
#meaningfully affect results, so going forward I'll use the observed standard deviation for MHO.

MLTbirthvariancenozeroes <- var(MLTdfnozeroesinbirthrate$Seeds_by_ind)
MLTbirthvariancenozeroes
MLTbirthlogvariancenozeroes <- var(log(MLTdfnozeroesinbirthrate$Seeds_by_ind))
MLTbirthlogstdevnozeroes <- sqrt(MLTbirthlogvariancenozeroes)
MLTbirthlogmeannozeroes <- mean(log(MLTdfnozeroesinbirthrate$Seeds_by_ind))
MLTbirthlogmeannozeroes
MLTbirthlogvariancenozeroes
MLTbirthlogstdevnozeroes
MLTstdevtotry <- seq(1.005, 1.40, 0.005)
MLTresultingvariances <- c(1.0:80.0)
MLTaddingwithingenotypevariance <- c(1:nrow(MLTdfnozeroesinbirthrate))
for (i in 1:80) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rlnorm(nrow(MLTdfnozeroesinbirthrate), MLTbirthlogmeannozeroes, MLTstdevtotry[i])
    for (j in 1:nrow(MLTdfnozeroesinbirthrate)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(MLTdfnozeroesinbirthrate)) {
    MLTaddingwithingenotypevariance[n] <- rpois(1, resampledamonggenotypetruebirths[n])
  }
  MLTresultingvariances[i] <- var(MLTaddingwithingenotypevariance)
}
MLTresultingvariances

MLTregressionvariancebystdev <- lm(MLTresultingvariances ~ MLTstdevtotry)
MLTregressionvariancebystdev

#1986315 = (29369067 * correctstdev) - 29643097
(1986315 + 29643097) / 29369067

#This analysis predicts that fitting a log-normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) does overestimate the true variance between genotypes.
#Going forward, for MLT I'll draw 'true' genotype birth rates from a log-normal with the
#standard deviation predicted from this model: 1.076963.

MLObirthvariancenozeroes <- var(MLOdfnozeroesinbirthrate$Seeds_by_ind)
MLObirthvariancenozeroes
MLObirthlogvariancenozeroes <- var(log(MLOdfnozeroesinbirthrate$Seeds_by_ind))
MLObirthlogstdevnozeroes <- sqrt(MLObirthlogvariancenozeroes)
MLObirthlogmeannozeroes <- mean(log(MLOdfnozeroesinbirthrate$Seeds_by_ind))
MLObirthlogmeannozeroes
MLObirthlogvariancenozeroes
MLObirthlogstdevnozeroes
MLOstdevtotry <- seq(0.4025, 0.70, 0.0025)
MLOresultingvariances <- c(1.0:120.0)
MLOaddingwithingenotypevariance <- c(1:nrow(MLOdfnozeroesinbirthrate))
for (i in 1:120) {
  checkfornegativebirthrates <- TRUE
  while(checkfornegativebirthrates) {
    atleastonenegativebirthrate <- FALSE
    resampledamonggenotypetruebirths <- rlnorm(nrow(MLOdfnozeroesinbirthrate), MLObirthlogmeannozeroes, MLOstdevtotry[i])
    for (j in 1:nrow(MLOdfnozeroesinbirthrate)) {
      if (resampledamonggenotypetruebirths[j] < 1.0) {
        atleastonenegativebirthrate <- TRUE
      }
    }
    if (!atleastonenegativebirthrate) {
      checkfornegativebirthrates <- FALSE
    }
  }
  for (n in 1:nrow(MLOdfnozeroesinbirthrate)) {
    MLOaddingwithingenotypevariance[n] <- rpois(1, resampledamonggenotypetruebirths[n])
  }
  MLOresultingvariances[i] <- var(MLOaddingwithingenotypevariance)
}
MLOresultingvariances

MLOregressionvariancebystdev <- lm(MLOresultingvariances ~ MLOstdevtotry)
MLOregressionvariancebystdev

#12356463 = (27027911 * correctstdev) - 9408809
(12356463 + 9408809) / 27027911

#This analysis predicts that fitting a log-normal distribution to the observed distribution
#of genotypes and then adding in within-genotype variance (which is already included in observed 
#distribution) not only doesn't overestimate variance between genotypes, it UNDERESTIMATES it.
#I'm not totally sure why, the distribution seems pretty log-normal to me.
#Going forward, I'll use the observed variance for MLO.

#Following code is simply to visualize the distribution of birth rates actually observed.
hist(TLTdfnozeroesinbirthrate$Seeds_by_ind)
hist(THTdfnozeroesinbirthrate$Seeds_by_ind)
hist(THOdfnozeroesinbirthrate$Seeds_by_ind)
hist(TLOdfnozeroesinbirthrate$Seeds_by_ind)
hist(MLTdfnozeroesinbirthrate$Seeds_by_ind)
hist(MLOdfnozeroesinbirthrate$Seeds_by_ind)
hist(MHTdfnozeroesinbirthrate$Seeds_by_ind)
hist(MHOdfnozeroesinbirthrate$Seeds_by_ind)
hist(log(MLTdfnozeroesinbirthrate$Seeds_by_ind))
hist(log(MLOdfnozeroesinbirthrate$Seeds_by_ind))
hist(log(MHTdfnozeroesinbirthrate$Seeds_by_ind))
hist(log(TLTdfnozeroesinbirthrate$Seeds_by_ind))

#After all that other stuff is done, now we can do the actual simulations.
#First, create the appropriate data frames.
simulationdataframebirthsTHT <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(THTdfnozeroesinbirthrate)),
  genotypeID = THTdfnozeroesinbirthrate$Genotype_id,
  resampledgenotypebirthrates = THTdfnozeroesinbirthrate$Seeds_by_ind
)
resultsdataframebirthsTHT <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsTLT <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(TLTdfnozeroesinbirthrate)),
  genotypeID = TLTdfnozeroesinbirthrate$Genotype_id,
  resampledgenotypebirthrates = TLTdfnozeroesinbirthrate$Seeds_by_ind
)
resultsdataframebirthsTLT <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsTHO <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(THOdfnozeroesinbirthrate)),
  genotypeID = THOdfnozeroesinbirthrate$Genotype_id,
  resampledgenotypebirthrates = THOdfnozeroesinbirthrate$Seeds_by_ind
)
resultsdataframebirthsTHO <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsTLO <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(TLOdfnozeroesinbirthrate)),
  genotypeID = TLOdfnozeroesinbirthrate$Genotype_id,
  resampledgenotypebirthrates = TLOdfnozeroesinbirthrate$Seeds_by_ind
)
resultsdataframebirthsTLO <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsMHT <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(MHTdfnozeroesinbirthrate)),
  genotypeID = MHTdfnozeroesinbirthrate$Genotype_id,
  resampledgenotypebirthrates = MHTdfnozeroesinbirthrate$Seeds_by_ind
)
resultsdataframebirthsMHT <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsMLT <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(MLTdfnozeroesinbirthrate)),
  genotypeID = MLTdfnozeroesinbirthrate$Genotype_id,
  resampledgenotypebirthrates = MLTdfnozeroesinbirthrate$Seeds_by_ind
)
resultsdataframebirthsMLT <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsMHO <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(MHOdfnozeroesinbirthrate)),
  genotypeID = MHOdfnozeroesinbirthrate$Genotype_id,
  resampledgenotypebirthrates = MHOdfnozeroesinbirthrate$Seeds_by_ind
)
resultsdataframebirthsMHO <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

simulationdataframebirthsMLO <- data.frame(
  simulatedgenotypebirthrates = c(1:nrow(MLOdfnozeroesinbirthrate)),
  genotypeID = MLOdfnozeroesinbirthrate$Genotype_id,
  resampledgenotypebirthrates = MLOdfnozeroesinbirthrate$Seeds_by_ind
)
resultsdataframebirthsMLO <- data.frame(
  vectorofmaximumbirthrates = c(1.0:10000.0),
  vectoroftruebirthrates = c(1.0:10000.0),
  vectorofbestgenotypeID = c(1.0:10000.0)
)

#Now perform the simulations.
#I also include the assignment statements for relevant variables from the previous section,
#since that section isn't otherwise required for this section to run.
THTbirthvariancenozeroes <- var(THTdfnozeroesinbirthrate$Seeds_by_ind)
THTbirthstdevnozeroes <- sqrt(THTbirthvariancenozeroes)
THTbirthmeannozeroes <- mean(THTdfnozeroesinbirthrate$Seeds_by_ind)
for (i in 1:10000) {
  simulationdataframebirthsTHT$resampledgenotypebirthrates <- rnorm(nrow(THTdfnozeroesinbirthrate), THTbirthmeannozeroes, THTbirthstdevnozeroes)
  for(j in 1:nrow(THTdfnozeroesinbirthrate)) {
    if (simulationdataframebirthsTHT$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsTHT$simulatedgenotypebirthrates[j] <- (rpois(1, 1.0))
    } else {
      simulationdataframebirthsTHT$simulatedgenotypebirthrates[j] <- (rpois(1, simulationdataframebirthsTHT$resampledgenotypebirthrates[j]))
    }
  }
  resultsdataframebirthsTHT$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsTHT$simulatedgenotypebirthrates)
  resultsdataframebirthsTHT$vectoroftruebirthrates[i] <- max(simulationdataframebirthsTHT$resampledgenotypebirthrates)
  resultsdataframebirthsTHT$vectorofbestgenotypeID[i] <- simulationdataframebirthsTHT$genotypeID[which.max(simulationdataframebirthsTHT$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsTHT$vectorofmaximumbirthrates
resultsdataframebirthsTHT$vectoroftruebirthrates
resultsdataframebirthsTHT$vectorofbestgenotypeID
hist(resultsdataframebirthsTHT$vectorofbestgenotypeID)
hist(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_id)

#Note that TLT draws from a log-normal. 
TLTbirthlogvariancenozeroes <- var(log(TLTdfnozeroesinbirthrate$Seeds_by_ind))
TLTbirthlogstdevnozeroes <- sqrt(TLTbirthlogvariancenozeroes)
TLTbirthlogmeannozeroes <- mean(log(TLTdfnozeroesinbirthrate$Seeds_by_ind))
for (i in 1:10000) {
  simulationdataframebirthsTLT$resampledgenotypebirthrates <- rlnorm(nrow(TLTdfnozeroesinbirthrate), TLTbirthlogmeannozeroes, TLTbirthlogstdevnozeroes)
  for(j in 1:nrow(TLTdfnozeroesinbirthrate)) {
    if (simulationdataframebirthsTLT$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsTLT$simulatedgenotypebirthrates[j] <- (rpois(1, 1.0))
    } else {
      simulationdataframebirthsTLT$simulatedgenotypebirthrates[j] <- (rpois(1, simulationdataframebirthsTLT$resampledgenotypebirthrates[j]))
    }
  }
  resultsdataframebirthsTLT$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsTLT$simulatedgenotypebirthrates)
  resultsdataframebirthsTLT$vectoroftruebirthrates[i] <- max(simulationdataframebirthsTLT$resampledgenotypebirthrates)
  resultsdataframebirthsTLT$vectorofbestgenotypeID[i] <- simulationdataframebirthsTLT$genotypeID[which.max(simulationdataframebirthsTLT$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsTLT$vectorofmaximumbirthrates
resultsdataframebirthsTLT$vectoroftruebirthrates
resultsdataframebirthsTLT$vectorofbestgenotypeID
hist(resultsdataframebirthsTLT$vectorofbestgenotypeID)
hist(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_id)

THObirthvariancenozeroes <- var(THOdfnozeroesinbirthrate$Seeds_by_ind)
THObirthstdevnozeroes <- sqrt(THObirthvariancenozeroes)
THObirthmeannozeroes <- mean(THOdfnozeroesinbirthrate$Seeds_by_ind)
for (i in 1:10000) {
  simulationdataframebirthsTHO$resampledgenotypebirthrates <- rnorm(nrow(THOdfnozeroesinbirthrate), THObirthmeannozeroes, THObirthstdevnozeroes)
  for(j in 1:nrow(THOdfnozeroesinbirthrate)) {
    if (simulationdataframebirthsTHO$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsTHO$simulatedgenotypebirthrates[j] <- (rpois(1, 1.0))
    } else {
      simulationdataframebirthsTHO$simulatedgenotypebirthrates[j] <- (rpois(1, simulationdataframebirthsTHO$resampledgenotypebirthrates[j]))
    }
  }
  resultsdataframebirthsTHO$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsTHO$simulatedgenotypebirthrates)
  resultsdataframebirthsTHO$vectoroftruebirthrates[i] <- max(simulationdataframebirthsTHO$resampledgenotypebirthrates)
  resultsdataframebirthsTHO$vectorofbestgenotypeID[i] <- simulationdataframebirthsTHO$genotypeID[which.max(simulationdataframebirthsTHO$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsTHO$vectorofmaximumbirthrates
resultsdataframebirthsTHO$vectoroftruebirthrates
resultsdataframebirthsTHO$vectorofbestgenotypeID
hist(resultsdataframebirthsTHO$vectorofbestgenotypeID)
hist(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_id)

TLObirthvariancenozeroes <- var(TLOdfnozeroesinbirthrate$Seeds_by_ind)
TLObirthstdevnozeroes <- sqrt(TLObirthvariancenozeroes)
TLObirthmeannozeroes <- mean(TLOdfnozeroesinbirthrate$Seeds_by_ind)
for (i in 1:10000) {
  simulationdataframebirthsTLO$resampledgenotypebirthrates <- rnorm(nrow(TLOdfnozeroesinbirthrate), TLObirthmeannozeroes, TLObirthstdevnozeroes)
  for(j in 1:nrow(TLOdfnozeroesinbirthrate)) {
    if (simulationdataframebirthsTLO$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsTLO$simulatedgenotypebirthrates[j] <- (rpois(1, 1.0))
    } else {
      simulationdataframebirthsTLO$simulatedgenotypebirthrates[j] <- (rpois(1, simulationdataframebirthsTLO$resampledgenotypebirthrates[j]))
    }
  }
  resultsdataframebirthsTLO$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsTLO$simulatedgenotypebirthrates)
  resultsdataframebirthsTLO$vectoroftruebirthrates[i] <- max(simulationdataframebirthsTLO$resampledgenotypebirthrates)
  resultsdataframebirthsTLO$vectorofbestgenotypeID[i] <- simulationdataframebirthsTLO$genotypeID[which.max(simulationdataframebirthsTLO$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsTLO$vectorofmaximumbirthrates
resultsdataframebirthsTLO$vectoroftruebirthrates
resultsdataframebirthsTLO$vectorofbestgenotypeID
hist(resultsdataframebirthsTLO$vectorofbestgenotypeID)
hist(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_id)

#Note that MHT draws from a log-normal. 
MHTbirthlogvariancenozeroes <- var(log(MHTdfnozeroesinbirthrate$Seeds_by_ind))
MHTbirthlogstdevnozeroes <- sqrt(MHTbirthlogvariancenozeroes)
MHTbirthlogmeannozeroes <- mean(log(MHTdfnozeroesinbirthrate$Seeds_by_ind))
for (i in 1:10000) {
  simulationdataframebirthsMHT$resampledgenotypebirthrates <- rlnorm(nrow(MHTdfnozeroesinbirthrate), MHTbirthlogmeannozeroes, MHTbirthlogstdevnozeroes)
  for(j in 1:nrow(MHTdfnozeroesinbirthrate)) {
    if (simulationdataframebirthsMHT$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsMHT$simulatedgenotypebirthrates[j] <- (rpois(1, 1.0))
    } else {
      simulationdataframebirthsMHT$simulatedgenotypebirthrates[j] <- (rpois(1, simulationdataframebirthsMHT$resampledgenotypebirthrates[j]))
    }
  }
  resultsdataframebirthsMHT$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsMHT$simulatedgenotypebirthrates)
  resultsdataframebirthsMHT$vectoroftruebirthrates[i] <- max(simulationdataframebirthsMHT$resampledgenotypebirthrates)
  resultsdataframebirthsMHT$vectorofbestgenotypeID[i] <- simulationdataframebirthsMHT$genotypeID[which.max(simulationdataframebirthsMHT$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsMHT$vectorofmaximumbirthrates
resultsdataframebirthsMHT$vectoroftruebirthrates
resultsdataframebirthsMHT$vectorofbestgenotypeID
hist(resultsdataframebirthsMHT$vectorofbestgenotypeID)
hist(aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_id)

#Note that MLT draws from a log-normal.
#The simulation protocol overestimates the true variance in the log-normal, so I'm using 1.076963
#for the standard deviation of the log-normal to compensate, as found earlier.
MLTbirthlogvariancenozeroes <- var(log(MLTdfnozeroesinbirthrate$Seeds_by_ind))
MLTbirthlogstdevnozeroes <- sqrt(MLTbirthlogvariancenozeroes)
MLTbirthlogmeannozeroes <- mean(log(MLTdfnozeroesinbirthrate$Seeds_by_ind))
for (i in 1:10000) {
  simulationdataframebirthsMLT$resampledgenotypebirthrates <- rlnorm(nrow(MLTdfnozeroesinbirthrate), MLTbirthlogmeannozeroes, 1.076963)
  for(j in 1:nrow(MLTdfnozeroesinbirthrate)) {
    if (simulationdataframebirthsMLT$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsMLT$simulatedgenotypebirthrates[j] <- (rpois(1, 1.0))
    } else {
      simulationdataframebirthsMLT$simulatedgenotypebirthrates[j] <- (rpois(1, simulationdataframebirthsMLT$resampledgenotypebirthrates[j]))
    }
  }
  resultsdataframebirthsMLT$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsMLT$simulatedgenotypebirthrates)
  resultsdataframebirthsMLT$vectoroftruebirthrates[i] <- max(simulationdataframebirthsMLT$resampledgenotypebirthrates)
  resultsdataframebirthsMLT$vectorofbestgenotypeID[i] <- simulationdataframebirthsMLT$genotypeID[which.max(simulationdataframebirthsMLT$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsMLT$vectorofmaximumbirthrates
resultsdataframebirthsMLT$vectoroftruebirthrates
resultsdataframebirthsMLT$vectorofbestgenotypeID
hist(resultsdataframebirthsMLT$vectorofbestgenotypeID)
hist(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_id)

#Note that MLO draws from a log-normal.
MLObirthlogvariancenozeroes <- var(log(MLOdfnozeroesinbirthrate$Seeds_by_ind))
MLObirthlogstdevnozeroes <- sqrt(MLObirthlogvariancenozeroes)
MLObirthlogmeannozeroes <- mean(log(MLOdfnozeroesinbirthrate$Seeds_by_ind))
for (i in 1:10000) {
  simulationdataframebirthsMLO$resampledgenotypebirthrates <- rlnorm(nrow(MLOdfnozeroesinbirthrate), MLObirthlogmeannozeroes, MLObirthlogstdevnozeroes)
  for(j in 1:nrow(MLOdfnozeroesinbirthrate)) {
    if (simulationdataframebirthsMLO$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsMLO$simulatedgenotypebirthrates[j] <- (rpois(1, 1.0))
    } else {
      simulationdataframebirthsMLO$simulatedgenotypebirthrates[j] <- (rpois(1, simulationdataframebirthsMLO$resampledgenotypebirthrates[j]))
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
hist(resultsdataframebirthsMLO$vectorofbestgenotypeID)
hist(aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_id)


MHObirthvariancenozeroes <- var(MHOdfnozeroesinbirthrate$Seeds_by_ind)
MHObirthstdevnozeroes <- sqrt(MHObirthvariancenozeroes)
MHObirthmeannozeroes <- mean(MHOdfnozeroesinbirthrate$Seeds_by_ind)
for (i in 1:10000) {
  simulationdataframebirthsMHO$resampledgenotypebirthrates <- rnorm(nrow(MHOdfnozeroesinbirthrate), MHObirthmeannozeroes, MHObirthstdevnozeroes)
  for(j in 1:nrow(MHOdfnozeroesinbirthrate)) {
    if (simulationdataframebirthsMHO$resampledgenotypebirthrates[j] < 1.0) {
      simulationdataframebirthsMHO$simulatedgenotypebirthrates[j] <- (rpois(1, 1.0))
    } else {
      simulationdataframebirthsMHO$simulatedgenotypebirthrates[j] <- (rpois(1, simulationdataframebirthsMHO$resampledgenotypebirthrates[j]))
    }
  }
  resultsdataframebirthsMHO$vectorofmaximumbirthrates[i] <- max(simulationdataframebirthsMHO$simulatedgenotypebirthrates)
  resultsdataframebirthsMHO$vectoroftruebirthrates[i] <- max(simulationdataframebirthsMHO$resampledgenotypebirthrates)
  resultsdataframebirthsMHO$vectorofbestgenotypeID[i] <- simulationdataframebirthsMHO$genotypeID[which.max(simulationdataframebirthsMHO$simulatedgenotypebirthrates)]
}
warnings()
resultsdataframebirthsMHO$vectorofmaximumbirthrates
resultsdataframebirthsMHO$vectoroftruebirthrates
resultsdataframebirthsMHO$vectorofbestgenotypeID
hist(resultsdataframebirthsMHO$vectorofbestgenotypeID)
hist(aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_id)


#Calculate expected bias in fittest genotype due to extreme value issue.
biasvectorbirthsTHT <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsTHT[i] = (resultsdataframebirthsTHT$vectorofmaximumbirthrates[i] - resultsdataframebirthsTHT$vectoroftruebirthrates[i])
}
mean(biasvectorbirthsTHT)

biasvectorbirthsTLT <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsTLT[i] = (resultsdataframebirthsTLT$vectorofmaximumbirthrates[i] - resultsdataframebirthsTLT$vectoroftruebirthrates[i])
}
mean(biasvectorbirthsTLT)

biasvectorbirthsTHO <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsTHO[i] = (resultsdataframebirthsTHO$vectorofmaximumbirthrates[i] - resultsdataframebirthsTHO$vectoroftruebirthrates[i])
}
mean(biasvectorbirthsTHO)

biasvectorbirthsTLO <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsTLO[i] = (resultsdataframebirthsTLO$vectorofmaximumbirthrates[i] - resultsdataframebirthsTLO$vectoroftruebirthrates[i])
}
mean(biasvectorbirthsTLO)

biasvectorbirthsMHT <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsMHT[i] = (resultsdataframebirthsMHT$vectorofmaximumbirthrates[i] - resultsdataframebirthsMHT$vectoroftruebirthrates[i])
}
mean(biasvectorbirthsMHT)

biasvectorbirthsMLT <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsMLT[i] = (resultsdataframebirthsMLT$vectorofmaximumbirthrates[i] - resultsdataframebirthsMLT$vectoroftruebirthrates[i])
}
mean(biasvectorbirthsMLT)

biasvectorbirthsMHO <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsMHO[i] = (resultsdataframebirthsMHO$vectorofmaximumbirthrates[i] - resultsdataframebirthsMHO$vectoroftruebirthrates[i])
}
mean(biasvectorbirthsMHO)

biasvectorbirthsMLO <- c(1:10000)
for (i in 1:10000) {
  biasvectorbirthsMLO[i] = (resultsdataframebirthsMLO$vectorofmaximumbirthrates[i] - resultsdataframebirthsMLO$vectoroftruebirthrates[i])
}
mean(biasvectorbirthsMLO)

#Following lines use the estimate of bias from above to adjust number of selective deaths
#A line at the end shows the unadjusted number of selective deaths for comparison


aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[1] <- 0
for (i in 1:nrow(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata)) {
  aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[i] <- ((max_birth_rate_THT - 9.71) - aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i]) * aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
  if (aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[i] < 0.0) {
    aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[i] <- 0.0
  }
}
adjusted_total_selective_deaths_births_THT <- sum(aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births)
adjusted_total_selective_deaths_births_THT
total_selective_deaths_births_THT

aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[1] <- 0
for (i in 1:nrow(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata)) {
  aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[i] <- ((max_birth_rate_TLT - 0.9323) - aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i]) * aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
  if (aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[i] < 0.0) {
    aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[i] <- 0.0
  }
}
adjusted_total_selective_deaths_births_TLT <- sum(aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births)
adjusted_total_selective_deaths_births_TLT
total_selective_deaths_births_TLT

aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[1] <- 0
for (i in 1:nrow(aggregatedtuebingenhighwateroneseedselectivedeathsdata)) {
  aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] <- ((max_birth_rate_THO - 14.08) - aggregatedtuebingenhighwateroneseedselectivedeathsdata$Seeds_by_ind[i]) * aggregatedtuebingenhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]
  if (aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] < 0.0) {
    aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] <- 0.0
  }
}
adjusted_total_selective_deaths_births_THO <- sum(aggregatedtuebingenhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births)
adjusted_total_selective_deaths_births_THO
total_selective_deaths_births_THO

aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[1] <- 0
for (i in 1:nrow(aggregatedtuebingenlowwateroneseedselectivedeathsdata)) {
  aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] <- ((max_birth_rate_TLO - 10.36369) - aggregatedtuebingenlowwateroneseedselectivedeathsdata$Seeds_by_ind[i]) * aggregatedtuebingenlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]
  if (aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] < 0.0) {
    aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] <- 0.0
  }
}
adjusted_total_selective_deaths_births_TLO <- sum(aggregatedtuebingenlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births)
adjusted_total_selective_deaths_births_TLO
total_selective_deaths_births_TLO

#Note that the MHT condition had no positive bias, so I forego this analysis.
total_selective_deaths_births_MHT

aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[1] <- 0
for (i in 1:nrow(aggregatedmadridlowwatertwentyseedsselectivedeathsdata)) {
  aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[i] <- ((max_birth_rate_MLT - 1.630724) - aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Seeds_by_ind[i]) * aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Surviving_individuals_pot[i]
  if (aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[i] < 0.0) {
    aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births[i] <- 0.0
  }
}
adjusted_total_selective_deaths_births_MLT <- sum(aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Adjusted_selective_deaths_births)
adjusted_total_selective_deaths_births_MLT
total_selective_deaths_births_MLT

aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[1] <- 0
for (i in 1:nrow(aggregatedmadridhighwateroneseedselectivedeathsdata)) {
  aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] <- ((max_birth_rate_MHO - 12.831) - aggregatedmadridhighwateroneseedselectivedeathsdata$Seeds_by_ind[i]) * aggregatedmadridhighwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]
  if (aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] < 0.0) {
    aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] <- 0.0
  }
}
adjusted_total_selective_deaths_births_MHO <- sum(aggregatedmadridhighwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births)
adjusted_total_selective_deaths_births_MHO
total_selective_deaths_births_MHO

aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[1] <- 0
for (i in 1:nrow(aggregatedmadridlowwateroneseedselectivedeathsdata)) {
  aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] <- ((max_birth_rate_MLO - 1.77246) - aggregatedmadridlowwateroneseedselectivedeathsdata$Seeds_by_ind[i]) * aggregatedmadridlowwateroneseedselectivedeathsdata$Surviving_individuals_pot[i]
  if (aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] < 0.0) {
    aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births[i] <- 0.0
  }
}
adjusted_total_selective_deaths_births_MLO <- sum(aggregatedmadridlowwateroneseedselectivedeathsdata$Adjusted_selective_deaths_births)
adjusted_total_selective_deaths_births_MLO
total_selective_deaths_births_MLO

#Some things I might like to put in a table of results.
max_birth_rate_THT
max_birth_rate_TLT
max_birth_rate_THO
max_birth_rate_TLO
max_birth_rate_MHT
max_birth_rate_MHO
max_birth_rate_MLT
max_birth_rate_MLO

1 - min_death_rate_MHO
1 - min_death_rate_MHT
1 - min_death_rate_MLO
1 - min_death_rate_MLT
1 - min_death_rate_THO
1 - min_death_rate_THT
1 - min_death_rate_TLO
1 - min_death_rate_TLT

mean(resultsdataframeMHO$vectorofmaximumsurvivalrates)
mean(resultsdataframeMHT$vectorofmaximumsurvivalrates)
mean(resultsdataframeMLO$vectorofmaximumsurvivalrates)
mean(resultsdataframeMLT$vectorofmaximumsurvivalrates)
mean(resultsdataframeTHO$vectorofmaximumsurvivalrates)
mean(resultsdataframeTHT$vectorofmaximumsurvivalrates)
mean(resultsdataframeTLO$vectorofmaximumsurvivalrates)
mean(resultsdataframeTLT$vectorofmaximumsurvivalrates)

mean(resultsdataframeMLO$vectoroftruesurvivalrates)
mean(resultsdataframeMHO$vectoroftruesurvivalrates)
mean(resultsdataframeTLO$vectoroftruesurvivalrates)
mean(resultsdataframeTHO$vectoroftruesurvivalrates)
mean(resultsdataframeMLT$vectoroftruesurvivalrates)
mean(resultsdataframeMHT$vectoroftruesurvivalrates)
mean(resultsdataframeTLT$vectoroftruesurvivalrates)
mean(resultsdataframeTHT$vectoroftruesurvivalrates)

mean(resultsdataframebirthsMLO$vectorofmaximumbirthrates)
mean(resultsdataframebirthsMHO$vectorofmaximumbirthrates)
mean(resultsdataframebirthsTLO$vectorofmaximumbirthrates)
mean(resultsdataframebirthsTHO$vectorofmaximumbirthrates)
mean(resultsdataframebirthsMLT$vectorofmaximumbirthrates)
mean(resultsdataframebirthsMHT$vectorofmaximumbirthrates)
mean(resultsdataframebirthsTLT$vectorofmaximumbirthrates)
mean(resultsdataframebirthsTHT$vectorofmaximumbirthrates)

mean(resultsdataframebirthsMLO$vectoroftruebirthrates)
mean(resultsdataframebirthsMHO$vectoroftruebirthrates)
mean(resultsdataframebirthsTLO$vectoroftruebirthrates)
mean(resultsdataframebirthsTHO$vectoroftruebirthrates)
mean(resultsdataframebirthsMLT$vectoroftruebirthrates)
mean(resultsdataframebirthsMHT$vectoroftruebirthrates)
mean(resultsdataframebirthsTLT$vectoroftruebirthrates)
mean(resultsdataframebirthsTHT$vectoroftruebirthrates)

THTbirthvariancenozeroes
TLTbirthvariancenozeroes
TLObirthvariancenozeroes
THObirthvariancenozeroes
MHTbirthvariancenozeroes
MLTbirthvariancenozeroes
MHObirthvariancenozeroes
MLObirthvariancenozeroes

average_seed_death_rate_THT
average_seed_death_rate_TLT
average_seed_death_rate_THO
average_seed_death_rate_TLO
average_seed_death_rate_MHT
average_seed_death_rate_MLT
average_seed_death_rate_MHO
average_seed_death_rate_MLO

adjusted_proportion_juvenile_deaths_selective_THT <- adjusted_total_selective_deaths_births_THT / (total_seeds_possible_THT - total_seeds_for_next_generation_THP)
adjusted_proportion_juvenile_deaths_selective_TLT <- adjusted_total_selective_deaths_births_TLT / (total_seeds_possible_TLT - total_seeds_for_next_generation_TLP)
adjusted_proportion_juvenile_deaths_selective_THO <- adjusted_total_selective_deaths_births_THO / (total_seeds_possible_THO - total_seeds_for_next_generation_THI)
adjusted_proportion_juvenile_deaths_selective_TLO <- adjusted_total_selective_deaths_births_TLO / (total_seeds_possible_TLO - total_seeds_for_next_generation_TLI)
adjusted_proportion_juvenile_deaths_selective_MHT <- total_selective_deaths_births_MHT / (total_seeds_possible_MHT - total_seeds_for_next_generation_MHP)
adjusted_proportion_juvenile_deaths_selective_MLT <- adjusted_total_selective_deaths_births_MLT / (total_seeds_possible_MLT - total_seeds_for_next_generation_MLP)
adjusted_proportion_juvenile_deaths_selective_MHO <- adjusted_total_selective_deaths_births_MHO / (total_seeds_possible_MHO - total_seeds_for_next_generation_MHI)
adjusted_proportion_juvenile_deaths_selective_MLO <- adjusted_total_selective_deaths_births_MLO / (total_seeds_possible_MLO - total_seeds_for_next_generation_MLI)

adjusted_proportion_juvenile_deaths_selective_THT
adjusted_proportion_juvenile_deaths_selective_TLT
adjusted_proportion_juvenile_deaths_selective_THO
adjusted_proportion_juvenile_deaths_selective_TLO
adjusted_proportion_juvenile_deaths_selective_MHT
adjusted_proportion_juvenile_deaths_selective_MLT
adjusted_proportion_juvenile_deaths_selective_MHO
adjusted_proportion_juvenile_deaths_selective_MLO

adjusted_total_selective_deaths_births <- adjusted_total_selective_deaths_births_MHO + adjusted_total_selective_deaths_births_MLO + total_selective_deaths_births_MHT + adjusted_total_selective_deaths_births_MLT + adjusted_total_selective_deaths_births_THO + adjusted_total_selective_deaths_births_TLO + adjusted_total_selective_deaths_births_THT + adjusted_total_selective_deaths_births_TLT
adjusted_total_selective_deaths_births

total_deaths_births <- (total_seeds_possible_MHO - total_seeds_for_next_generation_MHI) + (total_seeds_possible_MLO - total_seeds_for_next_generation_MLI) + (total_seeds_possible_MHT - total_seeds_for_next_generation_MHP) + (total_seeds_possible_MLT - total_seeds_for_next_generation_MLP) + (total_seeds_possible_THO - total_seeds_for_next_generation_THI) + (total_seeds_possible_TLO - total_seeds_for_next_generation_TLI) + (total_seeds_possible_THT - total_seeds_for_next_generation_THP) + (total_seeds_possible_TLT - total_seeds_for_next_generation_TLP)
total_deaths_births

adjusted_total_selective_deaths_births / total_deaths_births

#Following lines make a data frame out of reproductive excess and proportion of deaths which are selective
#at the life history stage of seed production, to produce a scatterplot of the two
#Ideally, the data from both life history stages would be in the same data frame, but because
#the density conditions need to be treated separately for the other life history stage,
#it's easier for me to just give those different data frames.
scatterplotbirthsdataframe <- data.frame(
  density = c("individual", "low", "individual", "high", "individual", "low", "individual", "high"),
  water = c("low", "low", "high", "high", "low", "low", "high", "high"),
  city = c("Madrid", "Madrid", "Madrid", "Madrid", "Tuebingen", "Tuebingen", "Tuebingen", "Tuebingen"),
  environmental_condition = c("MLI", "MLP", "MHI", "MHP", "TLI", "TLP", "THI", "THP"),
  reproductive_excess = c(kMLO, kMLT, kMHO, kMHT, kTLO, kTLT, kTHO, kTHT),
  proportion_of_deaths_selective = c(adjusted_proportion_juvenile_deaths_selective_MLO, adjusted_proportion_juvenile_deaths_selective_MLT, adjusted_proportion_juvenile_deaths_selective_MHO, adjusted_proportion_juvenile_deaths_selective_MHT, adjusted_proportion_juvenile_deaths_selective_TLO, adjusted_proportion_juvenile_deaths_selective_TLT, adjusted_proportion_juvenile_deaths_selective_THO, adjusted_proportion_juvenile_deaths_selective_THT)
  )

scatterplotlowdensitydeathsdataframe <- data.frame(
  water = c("low", "high", "low", "high"),
  city = c("Madrid", "Madrid", "Tuebingen", "Tuebingen"),
  environmental_condition = c("MLI", "MHI", "TLI", "THI"),
  reproductive_excess = c(kMLO, kMHO, kTLO, kTHO),
  proportion_of_deaths_selective = c(adjusted_proportion_deaths_selective_MLO, adjusted_proportion_deaths_selective_MHO, adjusted_proportion_deaths_selective_TLO, adjusted_proportion_deaths_selective_THO)
)

scatterplothighdensitydeathsdataframe <- data.frame(
  water = c("low", "high", "low", "high"),
  city = c("Madrid", "Madrid", "Tuebingen", "Tuebingen"),
  environmental_condition = c("MLP", "MHP", "TLP", "THP"),
  reproductive_excess = c(kMLT, kMHT, kTLT, kTHT),
  proportion_of_deaths_selective = c(adjusted_proportion_deaths_selective_MLT, adjusted_proportion_deaths_selective_MHT, adjusted_proportion_deaths_selective_TLT, adjusted_proportion_deaths_selective_THT)
)

cor.test(scatterplotbirthsdataframe$reproductive_excess, scatterplotbirthsdataframe$proportion_of_deaths_selective, method = "spearman")
cor.test(scatterplothighdensitydeathsdataframe$reproductive_excess, scatterplothighdensitydeathsdataframe$proportion_of_deaths_selective, method = "spearman")
cor.test(scatterplotlowdensitydeathsdataframe$reproductive_excess, scatterplotlowdensitydeathsdataframe$proportion_of_deaths_selective, method = "spearman")
pvaluesexcessvspropdeathsselective <- c(0.01071, 0.4167, 0.08333)
pchisq(-2*sum(log(pvaluesexcessvspropdeathsselective)), 6, lower.tail = FALSE)

resultsscatterplotbirths <- ggplot(data = scatterplotbirthsdataframe, aes(x = reproductive_excess, y = proportion_of_deaths_selective)) +
  geom_text(aes(label = environmental_condition))
resultsscatterplotbirths

resultsscatterplotbirths2 <- ggplot(data = scatterplotbirthsdataframe, aes(x = reproductive_excess, y = proportion_of_deaths_selective)) +
  geom_point()
resultsscatterplotbirths2

resultsscatterplotbirths3 <- ggplot(data = scatterplotbirthsdataframe, aes(x = reproductive_excess, y = proportion_of_deaths_selective, label = environmental_condition)) +
  geom_point()
resultsscatterplotbirths3 +
  geom_text_repel() +
  theme_classic() +
  ylim(0.0, 1.2) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective")

resultsscatterplotbirths4 <- ggplot(data = scatterplotbirthsdataframe, aes(x = reproductive_excess, y = proportion_of_deaths_selective, shape = city, color = water, fill = density))
resultsscatterplotbirths4 +
  ylim(0.0, 1.2) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective")


resultsscatterplotlowdensitydeaths <- ggplot(data = scatterplotlowdensitydeathsdataframe, aes(x = reproductive_excess, y = proportion_of_deaths_selective, label = environmental_condition)) +
  geom_point()
resultsscatterplotlowdensitydeaths +
  geom_text_repel() +
  theme_classic() +
  ylim(0.0, 1.2) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective")

resultsscatterplotlowdensitydeaths2 <- ggplot(data = scatterplotlowdensitydeathsdataframe, aes(x = reproductive_excess, y = proportion_of_deaths_selective, shape = city, color = water))
resultsscatterplotlowdensitydeaths2 +
  geom_point(size = 4, stroke = 1.15) +
  theme_classic() +
  guides(fill = "none", color = "none", shape = "none") +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_manual(name = "Water condition", values = c("low" = "red", "high" = "blue")) +
  ylim(0.0, 1.2) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective")

resultsscatterplothighdensitydeaths <- ggplot(data = scatterplothighdensitydeathsdataframe, aes(x = reproductive_excess, y = proportion_of_deaths_selective, label = environmental_condition)) +
  geom_point()
resultsscatterplothighdensitydeaths +
  geom_text_repel() +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 13.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "light gray"),
        panel.border = element_rect(fill = NA)) +
  ylim(0.0, 1.2) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective") +
  annotate("text", x = 1200, y = 1.0, size = 6.0, label = "High-density deaths\nSpearmans rho = -0.6")

resultsscatterplothighdensitydeaths2 <- ggplot(data = scatterplothighdensitydeathsdataframe, aes(x = reproductive_excess, y = proportion_of_deaths_selective, shape = city, color = water, fill = water))
resultsscatterplothighdensitydeaths2 +
  geom_point(size = 4, stroke = 1.15) +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 13.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "light gray"),
        panel.border = element_rect(fill = NA)) +
  guides(fill = "none", color = "none", shape = "none") +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_manual(name = "Water condition", values = c("low" = "red", "high" = "blue")) +
  scale_fill_manual(name = "Water condition", values = c("low" = "red", "high" = "blue")) +
  ylim(0.0, 1.2) +
  xlab("Reproductive excess") +
  ylab("Proportion of deaths selective") +
  annotate("text", x = 1200, y = 1.0, size = 6.0, label = "High-density deaths\nSpearmans rho = -0.6")

ggarrange(resultsscatterplotbirths3 +
            geom_text_repel() +
            theme(axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank(), 
                  axis.title.x = element_blank(), 
                  axis.title.y = element_text(size = 17),
                  axis.text.y = element_text(size = 13.5),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(colour = "light gray"),
                  panel.border = element_rect(fill = NA)) +
            xlim(0.0, 20000.0) +
            ylim(0.0, 1.2) +
            xlab("Reproductive excess") +
            ylab("Proportion of deaths selective") +
            annotate("text", x = 13500, y = 1.0, size = 6, label = "Seed production\nSpearmans rho = -0.8571"), 
          resultsscatterplotlowdensitydeaths +
            geom_text_repel() +
            theme(axis.text.x = element_text(size = 13.5), 
                  axis.ticks.x = element_blank(), 
                  axis.title.x = element_text(size = 17), 
                  axis.title.y = element_text(size = 17),
                  axis.text.y = element_text(size = 13.5),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(colour = "light gray"),
                  panel.border = element_rect(fill = NA)) +
            xlim(0.0, 20000.0) +
            ylim(0.0, 1.2) +
            xlab("Reproductive excess") +
            ylab("Proportion of deaths selective") +
            annotate("text", x = 13500, y = 1.0, size = 6, label = "Low-density deaths\nSpearmans rho = -1.0"), 
          nrow = 2)

ggarrange(resultsscatterplotbirths4 +
            theme(axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank(), 
                  axis.title.x = element_blank(), 
                  axis.title.y = element_text(size = 17),
                  axis.text.y = element_text(size = 13.5),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(colour = "light gray"),
                  panel.border = element_rect(fill = NA)) +
            geom_point(size = 4, stroke = 1.15) +
            guides(fill = "none", color = "none", shape = "none") +
            scale_shape_manual(values = c(21, 24)) +
            scale_color_manual(name = "Water condition", values = c("low" = "red", "high" = "blue")) +
            scale_fill_manual(name = "Density", values = c("individual" = "white", "low" = "red", "high" = "blue")) +
            xlim(0.0, 20000.0) +
            ylim(0.0, 1.2) +
            xlab("Reproductive excess") +
            ylab("Proportion of deaths selective") +
            annotate("text", x = 13500, y = 1.0, size = 6, label = "Seed production\nSpearmans rho = -0.8571"), 
          resultsscatterplotlowdensitydeaths2 +
            geom_point(size = 4, stroke = 1.15) +
            theme(axis.text.x = element_text(size = 13.5), 
                  axis.ticks.x = element_blank(), 
                  axis.title.x = element_text(size = 17), 
                  axis.title.y = element_text(size = 17),
                  axis.text.y = element_text(size = 13.5),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(colour = "light gray"),
                  panel.border = element_rect(fill = NA)) +
            guides(fill = "none", color = "none", shape = "none") +
            scale_shape_manual(values = c(21, 24)) +
            scale_color_manual(name = "Water condition", values = c("low" = "red", "high" = "blue")) +
            scale_fill_manual(name = "Water condition", values = c("low" = "red", "high" = "blue")) +
            xlim(0.0, 20000.0) +
            ylim(0.0, 1.2) +
            xlab("Reproductive excess") +
            ylab("Proportion of deaths selective") +
            annotate("text", x = 13500, y = 1.0, size = 6, label = "Low-density deaths\nSpearmans rho = -1.0"), 
          nrow = 2)

zoominonsmallbirthratesTLTplot <- ggplot(data = TLTnozeroesinbirthrate, aes(x = Seeds_by_ind)) +
  geom_histogram()
zoominonsmallbirthratesTLTplot +
  theme_classic() +
  xlim(0.0, 2000.0)

zoominonsmallbirthratesMLTplot <- ggplot(data = MLTnozeroesinbirthrate, aes(x = Seeds_by_ind)) +
  geom_histogram()
zoominonsmallbirthratesMLTplot +
  theme_classic() +
  xlim(0.0, 4000.0)

zoominonsmallbirthratesMLOplot <- ggplot(data = MLOnozeroesinbirthrate, aes(x = Seeds_by_ind)) +
  geom_histogram()
zoominonsmallbirthratesMLOplot +
  theme_classic() +
  xlim(0.0, 10000.0)

zoominonsmallbirthratesMHTplot <- ggplot(data = MHTnozeroesinbirthrate, aes(x = Seeds_by_ind)) +
  geom_histogram()
zoominonsmallbirthratesMHTplot +
  theme_classic() +
  xlim(0.0, 6000.0)

colnamesforgeneticmatrix <- unlist(read.table("Selective deaths/Genetic_distance_matrix.txt", nrow = 1, as.is = TRUE))
geneticdistancematrix <- as.matrix(read.table("Selective deaths/Genetic_distance_matrix.txt", fill = TRUE, skip = 1, row.names = 1))

#Note that genotype ID numbers 9098 and 9992 are in the selective deaths dataset,
#but not the genetic distance dataset.
#Those genotype ID numbers correspond to row numbers 224 and 514 in the selective deaths dataset.
#To simplify things, the following lines create new data frames without those genotypes.
MHIdataformatrix <- aggregatedmadridhighwateroneseedselectivedeathsdata[aggregatedmadridhighwateroneseedselectivedeathsdata$Genotype_id %in% colnamesforgeneticmatrix, ]
MHPdataformatrix <- aggregatedmadridhighwatertwentyseedsselectivedeathsdata[aggregatedmadridhighwatertwentyseedsselectivedeathsdata$Genotype_id %in% colnamesforgeneticmatrix, ]
MLIdataformatrix <- aggregatedmadridlowwateroneseedselectivedeathsdata[aggregatedmadridlowwateroneseedselectivedeathsdata$Genotype_id %in% colnamesforgeneticmatrix, ]
MLPdataformatrix <- aggregatedmadridlowwatertwentyseedsselectivedeathsdata[aggregatedmadridlowwatertwentyseedsselectivedeathsdata$Genotype_id %in% colnamesforgeneticmatrix, ]
THIdataformatrix <- aggregatedtuebingenhighwateroneseedselectivedeathsdata[aggregatedtuebingenhighwateroneseedselectivedeathsdata$Genotype_id %in% colnamesforgeneticmatrix, ]
THPdataformatrix <- aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata[aggregatedtuebingenhighwatertwentyseedsselectivedeathsdata$Genotype_id %in% colnamesforgeneticmatrix, ]
TLIdataformatrix <- aggregatedtuebingenlowwateroneseedselectivedeathsdata[aggregatedtuebingenlowwateroneseedselectivedeathsdata$Genotype_id %in% colnamesforgeneticmatrix, ]
TLPdataformatrix <- aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata[aggregatedtuebingenlowwatertwentyseedsselectivedeathsdata$Genotype_id %in% colnamesforgeneticmatrix, ]

typeof(colnamesforgeneticmatrix)

#This dataframe needs to have one entry for every pair of genotypes
dataframeforgeneticdistanceMHI <- data.frame(
  firstgenotypeID = c(1.0:132355.0),
  secondgenotypeID = c(1.0:132355.0),
  hammingdistance = c(1:132355),
  survivalpairwiseselectivedeaths = c(1.0:132355.0),
  birthspairwiseselectivedeaths = c(1.0:132355.0),
  survivalpercentdeathsselective = c(1.0:132355.0),
  birthspercentdeathsselective = c(1.0:132355.0)
)

dataframeforgeneticdistanceMLI <- data.frame(
  firstgenotypeID = c(1.0:132355.0),
  secondgenotypeID = c(1.0:132355.0),
  hammingdistance = c(1:132355),
  survivalpairwiseselectivedeaths = c(1.0:132355.0),
  birthspairwiseselectivedeaths = c(1.0:132355.0),
  survivalpercentdeathsselective = c(1.0:132355.0),
  birthspercentdeathsselective = c(1.0:132355.0)
)

dataframeforgeneticdistanceTHI <- data.frame(
  firstgenotypeID = c(1.0:132355.0),
  secondgenotypeID = c(1.0:132355.0),
  hammingdistance = c(1:132355),
  survivalpairwiseselectivedeaths = c(1.0:132355.0),
  birthspairwiseselectivedeaths = c(1.0:132355.0),
  survivalpercentdeathsselective = c(1.0:132355.0),
  birthspercentdeathsselective = c(1.0:132355.0)
)

dataframeforgeneticdistanceTLI <- data.frame(
  firstgenotypeID = c(1.0:132355.0),
  secondgenotypeID = c(1.0:132355.0),
  hammingdistance = c(1:132355),
  survivalpairwiseselectivedeaths = c(1.0:132355.0),
  birthspairwiseselectivedeaths = c(1.0:132355.0),
  survivalpercentdeathsselective = c(1.0:132355.0),
  birthspercentdeathsselective = c(1.0:132355.0)
)

dataframeforgeneticdistanceMHP <- data.frame(
  firstgenotypeID = c(1.0:132355.0),
  secondgenotypeID = c(1.0:132355.0),
  hammingdistance = c(1:132355),
  survivalpairwiseselectivedeaths = c(1.0:132355.0),
  birthspairwiseselectivedeaths = c(1.0:132355.0),
  survivalpercentdeathsselective = c(1.0:132355.0),
  birthspercentdeathsselective = c(1.0:132355.0)
)

dataframeforgeneticdistanceMLP <- data.frame(
  firstgenotypeID = c(1.0:131841.0),
  secondgenotypeID = c(1.0:131841.0),
  hammingdistance = c(1:131841),
  survivalpairwiseselectivedeaths = c(1.0:131841.0),
  birthspairwiseselectivedeaths = c(1.0:131841.0),
  survivalpercentdeathsselective = c(1.0:131841.0),
  birthspercentdeathsselective = c(1.0:131841.0)
)

dataframeforgeneticdistanceTHP <- data.frame(
  firstgenotypeID = c(1.0:132355.0),
  secondgenotypeID = c(1.0:132355.0),
  hammingdistance = c(1:132355),
  survivalpairwiseselectivedeaths = c(1.0:132355.0),
  birthspairwiseselectivedeaths = c(1.0:132355.0),
  survivalpercentdeathsselective = c(1.0:132355.0),
  birthspercentdeathsselective = c(1.0:132355.0)
)

dataframeforgeneticdistanceTLP <- data.frame(
  firstgenotypeID = c(1.0:132355.0),
  secondgenotypeID = c(1.0:132355.0),
  hammingdistance = c(1:132355),
  survivalpairwiseselectivedeaths = c(1.0:132355.0),
  birthspairwiseselectivedeaths = c(1.0:132355.0),
  survivalpercentdeathsselective = c(1.0:132355.0),
  birthspercentdeathsselective = c(1.0:132355.0)
)

#Following function calculates selective deaths for survival and births for all combinations of genotypes,
#as if each pair of genotypes was only competing against each other.
#The function takes two dataframes as parameters.
#The first df is one of the dataframes above, with two columns for genotype ID to cover all combinations of genotypes,
#one column for Hamming distance, and columns for absolute number of selective deaths and percent deaths selective for survival and births.
#The second df is one of the dataformatrix dfs above, which is the original dataframe for each environmental condition,
#after having calculated replicates, but with the genotypes not in the Hamming distance matrix removed.
#In addition, the function takes a third parameter simply to declare whether the environmental condition being calculated
#has an individual plant per pot or a population of thirty seeds per pot, since that changes some calculations.
#Note that this function will take minutes to run.
CalculatePairwiseSelectiveDeaths <- function(dataframeforcalculation, environmentalconditiondataframe, isPopulation) {
  counter <- 0
  for (i in 1:(nrow(environmentalconditiondataframe) - 1)) {
    for (j in (i+1):nrow(environmentalconditiondataframe)) {
      counter <- counter + 1
      dataframeforcalculation$firstgenotypeID[counter] <- environmentalconditiondataframe$Genotype_id[i]
      dataframeforcalculation$secondgenotypeID[counter] <- environmentalconditiondataframe$Genotype_id[j]
      dataframeforcalculation$hammingdistance[counter] <- geneticdistancematrix[j, i]
      deathrate1 <- environmentalconditiondataframe$Genotype_death_rate[i]
      deathrate2 <- environmentalconditiondataframe$Genotype_death_rate[j]
      if (isPopulation == 0) {
        if (deathrate1 < deathrate2) {
          dataframeforcalculation$survivalpairwiseselectivedeaths[counter] <- environmentalconditiondataframe$Replicates[j] * (deathrate2 - deathrate1)
          dataframeforcalculation$survivalpercentdeathsselective[counter] <- (dataframeforcalculation$survivalpairwiseselectivedeaths[counter] / (environmentalconditiondataframe$Replicates[j] + environmentalconditiondataframe$Replicates[i]))
        } else {
          dataframeforcalculation$survivalpairwiseselectivedeaths[counter] <- environmentalconditiondataframe$Replicates[i] * (deathrate1 - deathrate2)
          if (dataframeforcalculation$survivalpairwiseselectivedeaths[counter] > 0.0) {
            dataframeforcalculation$survivalpercentdeathsselective[counter] <- (dataframeforcalculation$survivalpairwiseselectivedeaths[counter] / (environmentalconditiondataframe$Replicates[i] + environmentalconditiondataframe$Replicates[j]))
          } else {
            dataframeforcalculation$survivalpercentdeathsselective[counter] <- 0.0
          }
        }
      } else {
        if (deathrate1 < deathrate2) {
          dataframeforcalculation$survivalpairwiseselectivedeaths[counter] <- 30.0 * environmentalconditiondataframe$Replicates[j] * (deathrate2 - deathrate1)
          dataframeforcalculation$survivalpercentdeathsselective[counter] <- (dataframeforcalculation$survivalpairwiseselectivedeaths[counter] / (30.0 * (environmentalconditiondataframe$Replicates[i] + environmentalconditiondataframe$Replicates[j])))
        } else {
          dataframeforcalculation$survivalpairwiseselectivedeaths[counter] <- 30.0 * environmentalconditiondataframe$Replicates[i] * (deathrate1 - deathrate2)
          if (dataframeforcalculation$survivalpairwiseselectivedeaths[counter] > 0.0) {
            dataframeforcalculation$survivalpercentdeathsselective[counter] <- (dataframeforcalculation$survivalpairwiseselectivedeaths[counter] / (30.0 * (environmentalconditiondataframe$Replicates[i] + environmentalconditiondataframe$Replicates[j])))
          } else {
            dataframeforcalculation$survivalpercentdeathsselective[counter] <- 0.0
          }
        }
      }
      birthrate1 <- environmentalconditiondataframe$Seeds_by_ind[i]
      birthrate2 <- environmentalconditiondataframe$Seeds_by_ind[j]
      if (birthrate1 > birthrate2) {
        dataframeforcalculation$birthspairwiseselectivedeaths[counter] <- (birthrate1 - birthrate2) * environmentalconditiondataframe$Surviving_individuals_pot[j]
        dataframeforcalculation$birthspercentdeathsselective[counter] <- dataframeforcalculation$birthspairwiseselectivedeaths[counter] / (birthrate1 * (environmentalconditiondataframe$Surviving_individuals_pot[i] + environmentalconditiondataframe$Surviving_individuals_pot[j]))
      } else {
        if (birthrate2 > 0.0) {
          dataframeforcalculation$birthspairwiseselectivedeaths[counter] <- (birthrate2 - birthrate1) * environmentalconditiondataframe$Surviving_individuals_pot[i]
          dataframeforcalculation$birthspercentdeathsselective[counter] <- dataframeforcalculation$birthspairwiseselectivedeaths[counter] / (birthrate2 * (environmentalconditiondataframe$Surviving_individuals_pot[i] + environmentalconditiondataframe$Surviving_individuals_pot[j]))
        } else {
          dataframeforcalculation$birthspairwiseselectivedeaths[counter] <- 0.0
          dataframeforcalculation$birthspercentdeathsselective[counter] <- 0.0
        }
      }
    }
  }
  return(dataframeforcalculation)
}

dataframeforgeneticdistanceMHI <- CalculatePairwiseSelectiveDeaths(dataframeforgeneticdistanceMHI, MHIdataformatrix, 0)
dataframeforgeneticdistanceMLI <- CalculatePairwiseSelectiveDeaths(dataframeforgeneticdistanceMLI, MLIdataformatrix, 0)
dataframeforgeneticdistanceTHI <- CalculatePairwiseSelectiveDeaths(dataframeforgeneticdistanceTHI, THIdataformatrix, 0)
dataframeforgeneticdistanceTLI <- CalculatePairwiseSelectiveDeaths(dataframeforgeneticdistanceTLI, TLIdataformatrix, 0)
dataframeforgeneticdistanceMHP <- CalculatePairwiseSelectiveDeaths(dataframeforgeneticdistanceMHP, MHPdataformatrix, 1)
dataframeforgeneticdistanceMLP <- CalculatePairwiseSelectiveDeaths(dataframeforgeneticdistanceMLP, MLPdataformatrix, 1)
dataframeforgeneticdistanceTHP <- CalculatePairwiseSelectiveDeaths(dataframeforgeneticdistanceTHP, THPdataformatrix, 1)
dataframeforgeneticdistanceTLP <- CalculatePairwiseSelectiveDeaths(dataframeforgeneticdistanceTLP, TLPdataformatrix, 1)

for (i in 1:515) {
  if (is.na(MLPdataformatrix$Genotype_death_rate[i])) {
    print(i)
  }
}

#This code is backup code in case the function above doesn't work properly.
#Make sure to reset the counter to zero before running the loop.
#I can't think of a more elegant way to iterate through the different combinations of genotypes.
counter <- 0
for (i in 1:514) {
  
  for (j in (i+1):515) {
    counter <- counter + 1
    dataframeforgeneticdistance$firstgenotypeID[counter] <- MHIdataformatrix$Genotype_id[i]
    dataframeforgeneticdistance$secondgenotypeID[counter] <- MHIdataformatrix$Genotype_id[j]
    dataframeforgeneticdistance$hammingdistance[counter] <- geneticdistancematrix[j, i]
    deathrate1 <- MHIdataformatrix$Genotype_death_rate[i]
    deathrate2 <- MHIdataformatrix$Genotype_death_rate[j]
    if (deathrate1 < deathrate2) {
      dataframeforgeneticdistance$survivalpairwiseselectivedeaths[counter] <- (MHIdataformatrix$Replicates[j] - MHIdataformatrix$Surviving_individuals_pot[j]) - (MHIdataformatrix$Replicates[j] * deathrate1)
      dataframeforgeneticdistance$survivalpercentdeathsselective[counter] <- (dataframeforgeneticdistance$survivalpairwiseselectivedeaths[counter] / (MHIdataformatrix$Replicates[j] + MHIdataformatrix$Replicates[i]))
    } else {
      dataframeforgeneticdistance$survivalpairwiseselectivedeaths[counter] <- (MHIdataformatrix$Replicates[i] - MHIdataformatrix$Surviving_individuals_pot[i]) - (MHIdataformatrix$Replicates[i] * deathrate2)
      if (dataframeforgeneticdistance$survivalpairwiseselectivedeaths[counter] > 0.0) {
        dataframeforgeneticdistance$survivalpercentdeathsselective[counter] <- (dataframeforgeneticdistance$survivalpairwiseselectivedeaths[counter] / (MHIdataformatrix$Replicates[i] + MHIdataformatrix$Replicates[j]))
      } else {
        dataframeforgeneticdistance$survivalpercentdeathsselective[counter] <- 0.0
      }
    }
    birthrate1 <- MHIdataformatrix$Seeds_by_ind[i]
    birthrate2 <- MHIdataformatrix$Seeds_by_ind[j]
    if (birthrate1 > birthrate2) {
      dataframeforgeneticdistance$birthspairwiseselectivedeaths[counter] <- (birthrate1 - birthrate2) * MHIdataformatrix$Surviving_individuals_pot[j]
    } else {
      dataframeforgeneticdistance$birthspairwiseselectivedeaths[counter] <- (birthrate2 - birthrate1) * MHIdataformatrix$Surviving_individuals_pot[i]
    }
    dataframeforgeneticdistance$birthspercentdeathsselective[counter] <- dataframeforgeneticdistance$birthspairwiseselectivedeaths[counter] / (MHIdataformatrix$Seeds_total_pot[i] + MHIdataformatrix$Seeds_total_pot[j])
  }
}


#This version of the graph draws a violin plot in bins across the x-axis where the width
#of the violin plot shows how many data points are in that bin.
#It also draws a smoothed line through the data. A loess line would be ideal, but the loess method
#in geom_smooth crashes my R when it tries to run on the 133k data points here.
MHIsurvivalhammingvspropboxplot <- ggplot() +
  geom_violin(data = dataframeforgeneticdistance, aes(group=cut_interval(hammingdistance, n = 18), x = hammingdistance, y = survivalpercentdeathsselective),
              fill= "azure2", colour = "darkslategray4", scale = "count") +
  geom_jitter(data = dataframeforgeneticdistance, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistance, aes(x = hammingdistance, y = survivalpercentdeathsselective),
              se = TRUE, color = 'gold')

#Given the lack of a relationship observed, a single scatterplot should suffice.
MHIsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMHI, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2)

#This version of the graph draws a violin plot in bins across the x-axis where the width
#of the violin plot shows how many data points are in that bin.
#It also draws a smoothed line through the data. A loess line would be ideal, but the loess method
#in geom_smooth crashes my R when it tries to run on the 133k data points here.
MHIbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMHI, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2)


MHIsurvivalfinalgraph <- MHIsurvivalhammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.35, size = 5.5, label = "MHI\nSpearmans rho = 0.0004\np-value = 0.8811")
MHIsurvivalfinalgraph

MHIbirthsfinalgraph <- MHIbirthshammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.4, size = 5.5, label = "MHI\nSpearmans rho = 0.0139\np-value = 3.988e-7")
MHIbirthsfinalgraph

MLIsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMLI, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceMLI, aes(x = hammingdistance, y = survivalpercentdeathsselective),
              se = TRUE, color = 'gold')
MLIsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMLI, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2)
MLIsurvivalfinalgraph <- MLIsurvivalhammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.63, size = 5.5, label = "MLI\nSpearmans rho = 0.096\np-value < 2.2e-16")
MLIsurvivalfinalgraph

MLIbirthshammingvspropboxplot <- ggplot() +
  geom_violin(data = dataframeforgeneticdistanceMLI, aes(group=cut_interval(hammingdistance, n = 18), x = hammingdistance, y = birthspercentdeathsselective),
              fill= "azure2", colour = "darkslategray4", scale = "count") +
  geom_jitter(data = dataframeforgeneticdistanceMLI, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceMLI, aes(x = hammingdistance, y = birthspercentdeathsselective),
              se = TRUE, color = 'gold')
MLIbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMLI, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2)
MLIbirthsfinalgraph <- MLIbirthshammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size = 20)) +
  annotate("text", x = 140000, y = 0.7, size = 5.5, label = "MLI\nSpearmans rho = 0.029\np-value < 2.2e-16")
MLIbirthsfinalgraph

THIsurvivalhammingvspropboxplot <- ggplot() +
  geom_violin(data = dataframeforgeneticdistanceTHI, aes(group=cut_interval(hammingdistance, n = 18), x = hammingdistance, y = survivalpercentdeathsselective),
              fill= "azure2", colour = "darkslategray4", scale = "count") +
  geom_jitter(data = dataframeforgeneticdistanceTHI, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceTHI, aes(x = hammingdistance, y = survivalpercentdeathsselective),
              se = TRUE, color = 'gold')
THIsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTHI, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2)
THIsurvivalfinalgraph <- THIsurvivalhammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.35, size = 5.5, label = "THI\nSpearmans rho = -0.0056\np-value = 0.0416")
THIsurvivalfinalgraph

THIbirthshammingvspropboxplot <- ggplot() +
  geom_violin(data = dataframeforgeneticdistanceTHI, aes(group=cut_interval(hammingdistance, n = 18), x = hammingdistance, y = birthspercentdeathsselective),
              fill= "azure2", colour = "darkslategray4", scale = "count") +
  geom_jitter(data = dataframeforgeneticdistanceTHI, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceTHI, aes(x = hammingdistance, y = birthspercentdeathsselective),
              se = TRUE, color = 'gold')
THIbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTHI, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2)
THIbirthsfinalgraph <- THIbirthshammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.4, size = 5.5, label = "THI\nSpearmans rho = -0.0056\np-value = 0.04149")
THIbirthsfinalgraph

TLIsurvivalhammingvspropboxplot <- ggplot() +
  geom_violin(data = dataframeforgeneticdistanceTLI, aes(group=cut_interval(hammingdistance, n = 18), x = hammingdistance, y = survivalpercentdeathsselective),
              fill= "azure2", colour = "darkslategray4", scale = "count") +
  geom_jitter(data = dataframeforgeneticdistanceTLI, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceTLI, aes(x = hammingdistance, y = survivalpercentdeathsselective),
              se = TRUE, color = 'gold')
TLIsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTLI, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2)
TLIsurvivalfinalgraph <- TLIsurvivalhammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.42, size = 5.5, label = "TLI\nSpearmans rho = -0.00237\np-value = 0.3884")
TLIsurvivalfinalgraph

TLIbirthshammingvspropboxplot <- ggplot() +
  geom_violin(data = dataframeforgeneticdistanceTLI, aes(group=cut_interval(hammingdistance, n = 18), x = hammingdistance, y = birthspercentdeathsselective),
              fill= "azure2", colour = "darkslategray4", scale = "count") +
  geom_jitter(data = dataframeforgeneticdistanceTLI, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceTLI, aes(x = hammingdistance, y = birthspercentdeathsselective),
              se = TRUE, color = 'gold')
TLIbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTLI, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2)
TLIbirthsfinalgraph <- TLIbirthshammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.6, size = 5.5, label = "TLI\nSpearmans rho = 0.0413\np-value < 2.2e-16")
TLIbirthsfinalgraph

MHPsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMHP, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2)
MHPsurvivalfinalgraph <- MHPsurvivalhammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 180000, y = 0.6, size = 5.5, label = "MHP\nSpearmans rho = 0.0111\np-value = 4.909e-5")
MHPsurvivalfinalgraph

MHPbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMHP, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceMHP, aes(x = hammingdistance, y = birthspercentdeathsselective),
              se = TRUE, color = 'gold')
MHPbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMHP, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2)
MHPbirthsfinalgraph <- MHPbirthshammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.77, size = 5.5, label = "MHP\nSpearmans rho = -0.0258\np-value < 2.2e-16")
MHPbirthsfinalgraph

MLPsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMLP, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceMLP, aes(x = hammingdistance, y = survivalpercentdeathsselective),
              se = TRUE, color = 'gold')
MLPsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMLP, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2)
MLPsurvivalfinalgraph <- MLPsurvivalhammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  ylim(-0.01, 1.01) +
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.85, size = 5.5, label = "MLP\nSpearmans rho = 0.0547\np-value < 2.2e-16")
MLPsurvivalfinalgraph

MLPbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMLP, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceMLP, aes(x = hammingdistance, y = birthspercentdeathsselective),
              se = TRUE, color = 'gold')
MLPbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceMLP, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2)
MLPbirthsfinalgraph <- MLPbirthshammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") +
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.77, size = 5.5, label = "MLP\nSpearmans rho = 0.0455\np-value < 2.2e-16")
MLPbirthsfinalgraph

THPsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTHP, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceTHP, aes(x = hammingdistance, y = survivalpercentdeathsselective),
              se = TRUE, color = 'gold')
THPsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTHP, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2)
THPsurvivalfinalgraph <- THPsurvivalhammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.37, size = 5.5, label = "THP\nSpearmans rho = -0.0085\np-value = 0.001995")
THPsurvivalfinalgraph

THPbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTHP, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceTHP, aes(x = hammingdistance, y = birthspercentdeathsselective),
              se = TRUE, color = 'gold')
THPbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTHP, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2)
THPbirthsfinalgraph <- THPbirthshammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.6, size = 5.5, label = "THP\nSpearmans rho = -0.00288\np-value = 0.2942")
THPbirthsfinalgraph

TLPsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTLP, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceTLP, aes(x = hammingdistance, y = survivalpercentdeathsselective),
              se = TRUE, color = 'gold')
TLPsurvivalhammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTLP, aes(x = hammingdistance, y = survivalpercentdeathsselective), color="black", size=0.4, alpha=0.2)
TLPsurvivalfinalgraph <- TLPsurvivalhammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 175000, y = 0.365, size = 5.5, label = "TLP\nSpearmans rho = -0.00813\np-value = 0.003085")
TLPsurvivalfinalgraph

TLPbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTLP, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2) +
  geom_smooth(data = dataframeforgeneticdistanceTLP, aes(x = hammingdistance, y = birthspercentdeathsselective),
              se = TRUE, color = 'gold')
TLPbirthshammingvspropboxplot <- ggplot() +
  geom_jitter(data = dataframeforgeneticdistanceTLP, aes(x = hammingdistance, y = birthspercentdeathsselective), color="black", size=0.4, alpha=0.2)
TLPbirthsfinalgraph <- TLPbirthshammingvspropboxplot +
  scale_x_continuous(name = "Hamming distance", labels = comma) +
  labs(y = "Proportion of deaths selective") + 
  theme_bw() + 
  theme_classic(base_size = 9) + 
  theme(text = element_text(size=20)) +
  annotate("text", x = 165000, y = 0.77, size = 5.5, label = "TLP\nSpearmans rho = -0.0577\np-value < 2.2e-16")
TLPbirthsfinalgraph

for (i in nrow(dataframeforgeneticdistanceMLP)) {
  if (dataframeforgeneticdistanceMLP$birthspercentdeathsselective[i] > 1.0) {
    print(dataframeforgeneticdistanceMLP$firstgenotypeID[i])
    print(dataframeforgeneticdistanceMLP$secondgenotypeID[i])
  }
}

cor.test(dataframeforgeneticdistanceMHI$hammingdistance, dataframeforgeneticdistanceMHI$survivalpercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceMHI$hammingdistance, dataframeforgeneticdistanceMHI$birthspercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceMLI$hammingdistance, dataframeforgeneticdistanceMLI$survivalpercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceMLI$hammingdistance, dataframeforgeneticdistanceMLI$birthspercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceTHI$hammingdistance, dataframeforgeneticdistanceTHI$survivalpercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceTHI$hammingdistance, dataframeforgeneticdistanceTHI$birthspercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceTLI$hammingdistance, dataframeforgeneticdistanceTLI$survivalpercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceTLI$hammingdistance, dataframeforgeneticdistanceTLI$birthspercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceMHP$hammingdistance, dataframeforgeneticdistanceMHP$survivalpercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceMHP$hammingdistance, dataframeforgeneticdistanceMHP$birthspercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceMLP$hammingdistance, dataframeforgeneticdistanceMLP$survivalpercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceMLP$hammingdistance, dataframeforgeneticdistanceMLP$birthspercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceTHP$hammingdistance, dataframeforgeneticdistanceTHP$survivalpercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceTHP$hammingdistance, dataframeforgeneticdistanceTHP$birthspercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceTLP$hammingdistance, dataframeforgeneticdistanceTLP$survivalpercentdeathsselective, method = "spearman")
cor.test(dataframeforgeneticdistanceTLP$hammingdistance, dataframeforgeneticdistanceTLP$birthspercentdeathsselective, method = "spearman")

fulldataframeMLP$Genotype_id[fulldataframeMLP$Genotype_death_rate == 0]

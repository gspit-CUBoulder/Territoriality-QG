# This is a script of MCMCglmm models to test for repeatability and heritability in 4 territorial phenotypes

library(tidyverse)
library(lubridate)
library(ggplot2)
library(viridis)
library(MCMCglmm)
library(pedtricks)
library(pedigree)

load("pedigree.rdata")

load("cache_size.rdata")

cache_sizes <- cache_sizes %>%
  mutate(animal = squirrel_id) %>%
  filter(!is.na(Grid),
         !is.na(Sex),
         !is.na(std_density),
         !is.na(std_local.density),
         !is.na(std_age))

load("IntrusionRate.rdata")

intruders_model <- intruders_model %>%
  filter(!is.na(Grid),
         !is.na(Sex),
         !is.na(std_density),
         !is.na(std_julian),
         !is.na(std_age),
         !is.na(std_local.density)) %>%
  mutate(animal = squirrel_id)

load("RattleProbability.rdata")

rattle_probability_final <- rattle_probability_final %>%
  filter(!is.na(Grid),
         !is.na(Sex),
         !is.na(std_density),
         !is.na(std_julian),
         !is.na(std_age)) %>%
  mutate(animal = squirrel_id)

load("TerritorySize.rdata")

data <- data %>%
  mutate(squirrel_id = as.factor(squirrel_id),
         animal = squirrel_id) %>%
  as.data.frame()

# Set Phenotypic Priors
PhenoPrior <- list(G=list(G1=list(V = diag(1), nu = 1, 
                                  alpha.mu = 0, alpha.V = diag(1)*1000),
                          G2=list(V = diag(1), nu = 1, 
                                  alpha.mu = 0, alpha.V = diag(1)*1000)),
                   R=list(V = diag(1) , nu = 1))

PhenoPrior_bi <- list(G=list(G1=list(V = diag(1), nu = 1, 
                                     alpha.mu = 0, alpha.V = diag(1)*1000),
                             G2=list(V = diag(1), nu = 1, 
                                     alpha.mu = 0, alpha.V = diag(1)*1000)),
                      R=list(V = diag(1) , fix = 1))

# Territory Size
TSModel1 <- MCMCglmm(area_ha_log ~ Grid + Sex + std_density,
                     random = ~squirrel_id + year.f,
                     prior = PhenoPrior,
                     family = "gaussian",
                     data = data,
                     nitt = 250000, burnin = 50000, thin = 200,  verbose=F)

RattleModel1 <- MCMCglmm(rattle_binary ~ Grid + Sex + spl(std_density, k=7) + 
                           spl(std_julian, k=4) + std_age,
                         random = ~squirrel_id + Year.f,
                         prior = PhenoPrior_bi,
                         family = "threshold",
                         data = rattle_probability_final, 
                         nitt = 250000, burnin = 50000, thin = 200, verbose=F)



IntrusionModel1 <- MCMCglmm(intruded ~ Grid + Sex + spl(std_density, k=8) + 
                              spl(std_julian, k=4) + std_age + std_local.density,
                            random = ~squirrel_id + Year.f,
                            prior = PhenoPrior_bi,
                            family = "threshold",
                            data = intruders_model, 
                            nitt = 250000, burnin = 50000, thin = 200, verbose=F)


CacheSizeModel1 <- MCMCglmm(cache_size_lastyear ~ Grid + Sex + std_density + 
                              std_age + I(std_age^2) + std_local.density,
                            random = ~squirrel_id + Year.f,
                            prior = PhenoPrior,
                            family = "gaussian",
                            data = cache_sizes, 
                            nitt = 250000, burnin = 50000, thin = 200, verbose=F)


# Set genetic prior

GenoPrior <- list(G=list(G1=list(V = diag(1), nu = 1, 
                                 alpha.mu = 0, alpha.V = diag(1)*1000),
                         G2=list(V = diag(1), nu = 1, 
                                 alpha.mu = 0, alpha.V = diag(1)*1000),
                         G3=list(V = diag(1), nu = 1, 
                                 alpha.mu = 0, alpha.V = diag(1)*1000)),
                  R=list(V = diag(1) , nu = 1))

GenoPrior_bi <- list(G=list(G1=list(V = diag(1), nu = 1, 
                                    alpha.mu = 0, alpha.V = diag(1)*1000),
                            G2=list(V = diag(1), nu = 1, 
                                    alpha.mu = 0, alpha.V = diag(1)*1000),
                            G3=list(V = diag(1), nu = 1, 
                                    alpha.mu = 0, alpha.V = diag(1)*1000)),
                     R=list(V = diag(1) , fix = 1))


TSModel2 <- MCMCglmm(area_ha_log ~ Grid + Sex + std_density,
                     random = ~squirrel_id + year.f + animal,
                     prior = GenoPrior,
                     pedigree = ped,
                     family = "gaussian",
                     data = data,
                     nitt = 2500000, burnin = 500000, thin = 2000,  verbose=F)


RattleModel2 <- MCMCglmm(rattle_binary ~ Grid + Sex + spl(std_density, k=7) + 
                           spl(std_julian, k=4) + std_age,
                         random = ~squirrel_id + Year.f + animal,
                         prior = GenoPrior_bi,
                         family = "threshold",
                         pedigree = ped,
                         data = rattle_probability_final, 
                         nitt = 250000, burnin = 50000, thin = 200, verbose=F)



IntrusionModel2 <- MCMCglmm(intruded ~ Grid + Sex + spl(std_density,k=7) + 
                              spl(std_julian,k=4) + std_age + std_local.density,
                            random = ~squirrel_id + Year.f + animal,
                            prior = GenoPrior_bi,
                            pedigree = ped,
                            family = "threshold",
                            data = intruders_model, 
                            nitt = 250000, burnin = 50000, thin = 200, verbose=F)



CacheSizeModel2 <- MCMCglmm(cache_size_lastyear ~ Grid + Sex + std_density + std_age +
                              I(std_age^2) + std_local_density,
                            random = ~squirrel_id + Year.f + animal,
                            prior = GenoPrior,
                            pedigree = ped,
                            family = "gaussian",
                            data = cache_sizes, 
                            nitt = 250000, burnin = 50000, thin = 200, verbose=F)


save(TSModel1, TSModel2, RattleModel1, RattleModel2,
     IntrusionModel1, IntrusionModel2, CacheSizeModel1, CacheSizeModel2, file ="Univariate_PEPrior_output.rdata")
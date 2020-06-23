# Make a list of the default parameters from Dietze et al. paper.

timestep <- 1800 # seconds

params <- list()
## hydrology
params$gevap <- 0.005 ## m2/s (Bonan p 201) [was 0.01, tuned]
params$Wthresh <- 1
params$Kroot <- 0.2 # umolH20*ha/(m3*kgRoot*s) -> uptake per kg root per m available water
params$SLA <- 10
params$alpha <- 0.8
params$Vcmax <- 18
params$Jmax <- params$Vcmax * 1.67
params$m <- 4
params$g0 <- 0
## allometry
params$allomB0 <- exp(-2.5355) / 0.8 ## Jenkins 2004, Pine, adj for BGB
params$allomB1 <- 2.4349
params$allomL0 <- exp(-2.907) ## Jenkins full database softwoods, Table 3
params$allomL1 <- 1.674
## plant respiration
params$Rleaf <- 0.04 * params$Vcmax
params$Rroot <- 1.2
params$Rstem <- 0.05
params$Rg <- 0.33
## soil respiration
params$Q10 <- 2
params$Rbasal <- 0.2 / (params$Q10^2.5) # umol/m2/sec per Mg/ha of SOM
## turnover
params$leafLitter <- 0.33 / 365 / 86400 * timestep
params$CWD <- 0.0001 / 365 / 86400 * timestep
params$rootLitter <- 1.0 / 365 / 86400 * timestep
## mortality
params$mort1 <- 1
params$mort2 <- 5
params$NSCthreshold <- 0.01
## NSC Allocation
params$Lmin <- 0.75
params$q <- 1
params$StoreMinDay <- 2
params$Smax <- 1
params$Rfrac <- 0.2
params$SeedlingMort <- 0.99
params$Kleaf <- (1 / 21 / 48) / 2^2.5 ## assumes it takes 21 days to regrow at 25C

default_parameters <- params
usethis::use_data(default_parameters, overwrite = TRUE)



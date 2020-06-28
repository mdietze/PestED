
#' Simple Ecosystem Model
#' @param X array [leaf,wood,root,storage,som,SoilWater,stem density]
#' @param params params
#' @param timestep is in seconds, defaults to 30 min
#' @param inputs PAR, temp, VPD
#' @param pest [phloem, xylem, leaf, root, stem] a vector of pest impacts, each entry in the vector
#' represents the percent change in thepholme, yxlem, leaf, root, and stem based on the type of disruption.
#' @author Michael C, Dietze <dietze@bu.edu>
#' @return X
#' @export
SEM <- function(X, params, inputs, pest = c(0, 0, 0, 1, 0), timestep = 1800) {

  # Check the parameter inputs
  extra_params <- which(!names(params) %in% names(default_parameters))
  msg_params   <- paste(names(params)[extra_params], collapse = ', ')
  assertthat::assert_that(length(extra_params) == 0, msg = paste0('params contains unknown parameters: ', msg_params))
  missing_params <- which(!names(default_parameters) %in% names(params))
  msg_params <- paste(names(params)[missing_params], collapse = ', ')
  assertthat::assert_that(length(missing_params) == 0, msg = paste0('params is missing the following: ', msg_params))

  # Define constants used to conver from unit A -> unit B
  ## convert umol/m2/sec -> Mg/ha/sec
  k <- 1e-6 * 12 * 1e-6 * 10000 # mol/umol*gC/mol*Mg/g*m2/ha
  ktree <- 1e-6 * 12 * 1e-3 # mol/umol*gC/mol*kg/g -> kg/m2/sec
  kH2O <- 1e-6 * 18 * 10^-6 # mol/umol*gH2O/mol*m/g
  P <- 101.325 ## average atm pressure (kPa)

  ### Hydrology ###

  ## evaporation (patch)
  EVAP <- min(X[6] / timestep, 1.15 * params$gevap * (0.622 * inputs$VPD / P) / 1000) ## Bonan Eqn 13.10: rate = m/s
  ## 1.15 kg/m3 = rho = density of air
  ## 1000 at end converts kg/m2/s to m/s
  ##**** CHECK UNITS (want m/s)
  X[6] <- X[6] + inputs$precip / 1000 - EVAP * timestep

  ## plant available water (trapazoidal response) (patch)
  paw <- ifelse(X[6] > params$Wthresh, X[6] - 0.5 * params$Wthresh,
                0.5 * X[6] * X[6] / params$Wthresh
  )
  supply <- (1 - pest[2]) * params$Kroot * X[3] * paw * X[7] ## potential rate of water uptake, umol/m2Ground/s

  ### turnover (tree)  ###
  leafLitter <- X[1] * min(1, params$leafLitter + pest[3])
  CWD <- X[2] * params$CWD
  rootLitter <- X[3] * pest[4] * params$rootLitter
  X[1] <- X[1] - leafLitter
  X[2] <- X[2] - CWD
  X[3] <- X[3] - rootLitter


  ## Leaf Respiration (m2 leaf)
  Rleaf <- arrhenius(params$Rleaf, inputs$temp) # umol/m2/sec

  ## LAI & Canopy Optics (patch)
  LAI <- X[1] * params$SLA * (X[7] / 10000) ## m2/m2
  PARmid <- (1 - 0.5 * exp(-0.5 * LAI)) * inputs$PAR ## estimate mid-canopy PAR,  umol/m2/sec

  ## photosynthesis & Transpiration (plot)
  Ags <- c(0, 0)
  if (inputs$PAR > 1e-20) {
    Ags <- solve.FVcB(
      arrhenius(params$Vcmax, inputs$temp),
      arrhenius(params$Jmax, inputs$temp),
      Rleaf,
      42.75 * exp(37830 * (inputs$temp + 273.15 - 298) / (298 * R * (inputs$temp + 273.15))), ## Gamma star
      params$alpha, ## quantum yield
      params$m, ## stomatal slope
      params$g0, ## cuticular conductance
      inputs$VPD,
      PARmid
    )
    demand <- max(Ags[2] * 0.622 * inputs$VPD / P * LAI * 1e6, 1e-10) ## transpiration without water limitation, umol/m2/s
    fopen <- max(0, min(1, supply / demand))
    GPP <- (Ags[1] + Rleaf) * fopen * LAI # umol/m2/sec
    TRANSP <- demand * fopen # umol/m2/sec
    X[6] <- X[6] - TRANSP * timestep * kH2O
  } else {
    fopen <- GPP <- TRANSP <- 0
  }

  ### allocation rules:
  #  1 - carbon is only available once it's in the phloem
  #  2 - maintainence respiration
  #  3 - minimum storage: maintain enough C for n days of Rm
  #  4 - minimum leaf & root
  #      roots:leaves constant (such that for climatic mean Waterbalance, maintain E[supply/demand]~1)
  #      but if stressed, allowed to loose leaves
  #  5 - storage: maintain a proportion of potential Bleaf
  #  6 - allometric maximum leaf & root
  #  7 - stem growth

  ###### respiration & allocation ########

  ## maintainence respiration (priority #2) (umol/s/tree)
  GPP <- GPP * 10000 / X[7]
  Rleaf <- Rleaf * LAI * 10000 / X[7]
  Rstem <- X[2] * arrhenius(params$Rstem, inputs$temp)
  Rroot <- X[3] * arrhenius(params$Rroot, inputs$temp)
  Rg <- 0 ## growth respiration: kg per plant per timestep

  ## update storage for priorities 1 & 2 (kg/tree)
  X[4] <- X[4] + ((1 - pest[1]) * (GPP - Rleaf) - Rstem - Rroot) * ktree * timestep

  ## calculate allometric potentials
  DBH <- (X[2] / params$allomB0)^(1 / params$allomB1) ## infer DBH from woody biomas
  Lmax <- params$allomL0 * DBH^params$allomL1 ## set maximum leaf biomass from DBH
  Lmin <- params$Lmin * Lmax ## set minimum leaf and root biomass as a fraction of maximum
  Rmin <- Lmin * params$q
  Smin <- (Rleaf * LAI * 10000 / X[7] + Rstem + Rroot) * ktree * 86400 * params$StoreMinDay ## set minimum storage based on the number of days the plant could survive
  Smax <- params$Smax * Lmax ## Set maximum storage biomas as a multiplier to maximum leaf biomass (similar to Fisher et al 2010)

  ## priority 3: only allocate if store above minimum
  if (X[4] > Smin) {
    leafGrowMax <- Lmax * params$Kleaf * 2^(inputs$temp / 10) # thermal limit to growth

    ## priority #4 mimimum leaf and root
    if (X[1] < Lmin) {
      leafAlloc <- max(min(Lmin - X[1], (X[4] - Smin) / (1 + params$Rg), leafGrowMax), 0) # demand,supply
      X[1] <- X[1] + leafAlloc
      X[4] <- X[4] - leafAlloc * (1 + params$Rg)
      Rg <- Rg + leafAlloc * params$Rg
    }
    if (X[3] < Rmin) {
      rootAlloc <- max(min(Lmin - X[1], (X[4] - Smin) / (1 + params$Rg)), 0)
      X[3] <- X[3] + rootAlloc
      X[4] <- X[4] - rootAlloc * (1 + params$Rg)
      Rg <- Rg + rootAlloc * params$Rg
    }

    ## priority #5, maximum storage
    if (X[4] > Smax) {

      ## priority #6: Leaf and stem growth
      if (X[1] < Lmax | X[3] < Lmax * params$q) {
        leafDemand <- Lmax - X[1]
        rootDemand <- Lmax * params$q - X[3]
        storeSupply <- (X[4] - Smax) / (1 + params$Rg)
        falloc <- min(1, max(0, storeSupply / (leafDemand + rootDemand)))
        leafAlloc <- leafDemand * falloc
        rootAlloc <- rootDemand * falloc
        X[1] <- X[1] + leafAlloc
        X[3] <- X[3] + rootAlloc
        X[4] <- X[4] - (leafAlloc + rootAlloc) * (1 + params$Rg)
        Rg <- Rg + (leafAlloc + rootAlloc) * params$Rg
      }

      ## priority #7: Growth & reproduction
      if (X[4] > Smax) {
        growAlloc <- (X[4] - Smax) / (1 + params$Rg)
        reproAlloc <- growAlloc * params$Rfrac
        stemAlloc <- growAlloc - reproAlloc
        X[2] <- X[2] + stemAlloc
        X[4] <- X[4] - growAlloc * (1 + params$Rg)
        X[5] <- X[5] + reproAlloc * params$SeedlingMort ## bulk of reproductive allocation dies
        X[7] <- X[7] + reproAlloc * (1 - params$SeedlingMort) * X[7] / sum(X[1:4]) ## naive reproduction (new trees enter as adults)
        Rg <- Rg + growAlloc * params$Rg
      }
    } ## end Store > Smax
  } ## end Store > Smin

  ## mortality
  if (X[4] <= params$NSCthreshold * Smax) {
    X[7] <- 0
  } else {
    mortRate <- (pest[5] + params$mort1 * exp(-params$mort2 * X[4] / Smax)) * timestep / 86400 / 365
    X[5] <- X[5] + X[7] * mortRate * sum(X[1:4]) / 1000 ## dead trees go to SOM
    X[7] <- X[7] * (1 - mortRate) ## reduce density but not per-tree pools
  }

  ## soil respiration
  Rh <- min(params$Rbasal * X[5] * params$Q10^(inputs$temp / 10), X[5] / (k * timestep)) ## min ensures SOM never goes negative
  X[5] <- X[5] - Rh * k * timestep

  return(c(X, GPP, fopen, Rleaf, Rstem + Rroot, Rg / ktree / timestep))
  #  return(data.frame(X1=Xnew[,1],X2=Xnew[,2],X3=Xnew[,3],LAI,NEP=GPP-alloc[,1]-Rh,GPP,
  #                    Ra=alloc[,1],NPPw=alloc[,2],NPPl=alloc[,3],Rh,litter,CWD))


}

## Farquhar-Ball Berry Optimization Functions
farquhar <- function(Ci, Fparams, I) {

  ## photosynthesis related constants.
  R <- 8.3144621 ## ideal gas constant in J/K/mol
  Tleaf <- 298
  Kc <- 404.9 * exp(79430 * (Tleaf - 298) / (298 * R * Tleaf))
  Ko <- 278.4 * exp(36380 * (Tleaf - 298) / (298 * R * Tleaf))
  Km <- Kc * (1 + 210 / Ko)


  a <- 0.9 ## curvature parameter
  b <- -(Fparams[5] * I + Fparams[2])
  c <- Fparams[5] * I * Fparams[2]
  J <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
  aJ <- J * (Ci - Fparams[4]) / (4 * Ci + 8 * Fparams[4]) ## electron transport limited without covariates
  aC <- Fparams[1] * (Ci - Fparams[4]) / (Ci + Km)
  min(aJ, aC) - Fparams[3]
}

ballberry <- function(input, BBparams, Fparams, obs) {

  ## photosynthesis related constants.
  R <- 8.3144621 ## ideal gas constant in J/K/mol
  Tleaf <- 298
  Kc <- 404.9 * exp(79430 * (Tleaf - 298) / (298 * R * Tleaf))
  Ko <- 278.4 * exp(36380 * (Tleaf - 298) / (298 * R * Tleaf))
  Km <- Kc * (1 + 210 / Ko)

  ## is actually the Medlyn et al 2011 model
  Ci <- obs[1] - 1.6 * input[1] / input[2]
  e1 <- (farquhar(Ci, Fparams, obs[3]) - input[1])
  #  e2 <- (BBparams[1] + (1+BBparams[2]/sqrt(obs[2]))*input[1]/obs[1] - input[2])*100
  e2 <- (BBparams[1] + BBparams[2] * input[1] / ((obs[1] - Fparams[4]) * (1 + obs[2])) - input[2]) * 100
  return(e1^2 + e2^2)
}

solve.FVcB <- function(Vcmax, Jmax, Rleaf, Gstar, alpha, m, g0, VPD, PAR) {

  ## photosynthesis related constants.
  R <- 8.3144621 ## ideal gas constant in J/K/mol
  Tleaf <- 298
  Kc <- 404.9 * exp(79430 * (Tleaf - 298) / (298 * R * Tleaf))
  Ko <- 278.4 * exp(36380 * (Tleaf - 298) / (298 * R * Tleaf))
  Km <- Kc * (1 + 210 / Ko)


  Ca <- 400
  out <- stats::optim(c(15, 0.1), # solve simultaneously for An.pred and gs.pred
               ballberry,
               BBparams = c(g0, m), # Ball-Berry params
               Fparams = c(Vcmax, Jmax, Rleaf, Gstar, alpha), # Farquhar params
               obs = c(Ca, VPD, PAR)
  ) # data
  if (out$par[2] >= 0) {
    return(out$par)
  } else {
    return(c(0, 0))
  }
}

arrhenius <- function(observed.value, new.temp, old.temp = 25) {
  return(observed.value / exp(3000 * (1 / (273.15 + new.temp) - 1 / (273.15 + old.temp))))
}

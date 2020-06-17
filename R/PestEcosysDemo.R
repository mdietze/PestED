
#' Simple Ecosystem Model
#' @param X = [leaf,wood,root,storage,som,SoilWater,stem density]
#' @param params params
#' @param timestep is in seconds, defaults to 30 min
#' @param inputs: PAR, temp, VPD
#' @param pest [phloem, xylem, leaf, root, stem]
#' @author Michael C, Dietze <dietze@bu.edu>
#' @return X
#' @export
SEM <- function(X, params, inputs, pest = c(0, 0, 0, 1, 0), timestep = 1800) {

  # Define constants used to conver from unit A -> unit B
  ## convert umol/m2/sec -> Mg/ha/sec
  k <- 1e-6 * 12 * 1e-6 * 10000 # mol/umol*gC/mol*Mg/g*m2/ha
  ktree <- 1e-6 * 12 * 1e-3 # mol/umol*gC/mol*kg/g -> kg/m2/sec
  kH2O <- 1e-6 * 18 * 10^-6 # mol/umol*gH2O/mol*m/g
  P <- 101.325 ## average atm pressure (kPa)

  ## pest impacts:
  ## phloem feaders: % tax flux of carbon out of (GPP-Rl) and into Bstore
  ## xylem disrupters (bark beetle, canker, wilt, girdling): % decrease water supply
  ## defoliators: % foliage removed
  ## root rot: factor accelleration of root turnover
  ## stem rot: increase in background mortality

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
  a <- 0.9 ## curvature parameter
  b <- -(Fparams[5] * I + Fparams[2])
  c <- Fparams[5] * I * Fparams[2]
  J <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
  aJ <- J * (Ci - Fparams[4]) / (4 * Ci + 8 * Fparams[4]) ## electron transport limited without covariates
  aC <- Fparams[1] * (Ci - Fparams[4]) / (Ci + Km)
  min(aJ, aC) - Fparams[3]
}

ballberry <- function(input, BBparams, Fparams, obs) {
  ## is actually the Medlyn et al 2011 model
  Ci <- obs[1] - 1.6 * input[1] / input[2]
  e1 <- (farquhar(Ci, Fparams, obs[3]) - input[1])
  #  e2 <- (BBparams[1] + (1+BBparams[2]/sqrt(obs[2]))*input[1]/obs[1] - input[2])*100
  e2 <- (BBparams[1] + BBparams[2] * input[1] / ((obs[1] - Fparams[4]) * (1 + obs[2])) - input[2]) * 100
  return(e1^2 + e2^2)
}


solve.FVcB <- function(Vcmax, Jmax, Rleaf, Gstar, alpha, m, g0, VPD, PAR) {
  Ca <- 400
  out <- optim(c(15, 0.1), # solve simultaneously for An.pred and gs.pred
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


#' main
#'
#' @return nothing
#' @export
#' @importFrom ncdf4 nc_open ncvar_get nc_close
main <- function() {

  ### paramters
  timestep <- 1800 # seconds
  params <- list()
  ## hydrology
  params$gevap <- 0.005 ## m2/s (Bonan p 201) [was 0.01, tuned]
  params$Wthresh <- 1
  params$Kroot <- 0.2 # umolH20*ha/(m3*kgRoot*s) -> uptake per kg root per m available water
  ## photosynthesis
  R <- 8.3144621 ## ideal gas constant in J/K/mol
  Tleaf <- 298
  Kc <- 404.9 * exp(79430 * (Tleaf - 298) / (298 * R * Tleaf))
  Ko <- 278.4 * exp(36380 * (Tleaf - 298) / (298 * R * Tleaf))
  Km <- Kc * (1 + 210 / Ko)
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

  ## initialize state variables
  DBH <- 10
  X <- rep(1, 7)
  X[1] <- X[3] <- X[4] <- params$allomL0 * DBH^params$allomL1
  X[2] <- params$allomB0 * DBH^params$allomB1
  X[5] <- 10
  X[7] <- 700


  if (!exists("inputs")) {
    met <- nc_open("AMF_USMe2_2005_L2_GF_V006.nc")
    print(met)
    PAR <- ncvar_get(met, "PAR")
    for (i in which(PAR < -10)) {
      PAR[i] <- PAR[i - 1]
    } ## uber-naive gapfilling
    temp <- ncvar_get(met, "TA")
    VPD <- ncvar_get(met, "VPD")
    precip <- ncvar_get(met, "PREC")
    time <- ncvar_get(met, "DOY")
    nc_close(met)
    plot(PAR, type = "l")
    plot(temp, type = "l")
    plot(VPD, type = "l")
    plot(precip, type = "l")
    inputs <- data.frame(PAR = PAR, temp = temp, VPD = VPD, precip = precip)
  }

  varnames <- c("Bleaf", "Bwood", "Broot", "Bstore", "BSOM", "Water", "density", "GPP", "fopen", "Rleaf", "RstemRroot", "Rgrow")
  units <- c("kg/plant", "kg/plant", "kg/plant", "kg/plant", "Mg/ha", "m", "stems/ha")

  iterate.SEM <- function(pest, t.start = 7000, years = 1) {
    pest.orig <- pest
    pest <- c(0, 0, 0, 1, 0)
    nt <- length(time) * years
    output <- array(NA, c(nt, 12))
    for (t in 1:nt) {
      ## turn pests on/off
      if (t %in% t.start) {
        pest <- pest.orig
      } else {
        pest[3] <- 0 ## defoliation turned off
      }

      ti <- (t - 1) %% nrow(inputs) + 1 ## indexing to allow met to loop
      output[t, ] <- SEM(X, params, inputs[ti, ], pest)
      X <- output[t, 1:7]
      if ((t %% (48 * 7)) == 0) print(t / 48) ## day counter
      if (X[7] == 0) break
    }

    colnames(output) <- varnames
    return(output)
  } # end iterate.SEM

  plot.SEM <- function(output) {
    output <- as.data.frame(output)
    for (i in 1:ncol(output)) {
      if (i <= 7) {
        plot(output[, i], type = "l", main = varnames[i], ylab = units[i])
      } else {
        plot(tapply(output[, i], rep(1:366, each = 48)[1:nrow(output)], mean), main = varnames[i], type = "l")
      }
    }
  }

  if (FALSE) {
    default <- iterate.SEM(c(0, 0, 0, 1, 0))
    plot.SEM(default)

    defol <- iterate.SEM(c(0, 0, 1, 1, 0)) ## assume a one-time 100% defoliation
    plot.SEM(defol)
    plot.SEM(default - defol)

    beetle <- iterate.SEM(c(0, 0.8, 0, 1, 0), years = 4) ## assume a 80% reduction in conductivity
    plot.SEM(beetle)
    plot.SEM(default - beetle)

    1 - apply(default, 2, min) / apply(default, 2, max)
    1 - apply(defol, 2, min) / apply(defol, 2, max)

    L4 <- read.csv("AMF_USMe2_2005_L4_h_V002.txt", header = TRUE, na.strings = "-9999")
    L4[L4 == -9999] <- NA

    default <- as.data.frame(default)

    ## GPP: model and observed
    GPP <- default[, 8] * default[, 7] / 10000 ## convert back to umol/m2/sec to compart to tower
    plot(GPP, type = "l")
    points(L4$GPP_st_MDS, pch = ".", col = 2)
    plot(GPP, L4$GPP_st_MDS, pch = ".")
    mean(GPP) / mean(L4$GPP_st_MDS)

    ## GPP Diurnal
    GPP.mod.diurnal <- tapply(GPP, L4$Hour, mean)
    GPP.obs.diurnal <- tapply(L4$GPP_st_MDS, L4$Hour, mean)
    ylim <- range(c(GPP.mod.diurnal, GPP.obs.diurnal))
    tod <- sort(unique(L4$Hour))
    plot(tod, GPP.mod.diurnal, ylim = ylim, col = 2, xlab = "Time of Day", ylab = "GPP", main = "Diurnal Cycle", type = "l", lwd = 3)
    lines(tod, GPP.obs.diurnal, lwd = 3)
    legend("topleft", legend = c("obs", "mod"), col = 1:2, pch = 20, cex = 0.75)

    ## RA & NPP (umol/sec/tree)
    RA <- default$Rleaf + default$RstemRroot + default$Rgrow
    NPP <- default[, 8] - RA
    mean(NPP) / mean(default[, 8])
    Rplant <- apply(default[, c("Rleaf", "RstemRroot", "Rgrow")], 2, mean)
    Rplant / sum(Rplant)

    ## woody increment
    DBH <- (default$Bwood / params$allomB0)^(1 / params$allomB1) ## infer DBH from woody biomas
    plot(DBH)
    inc <- DBH[length(DBH)] - DBH[1]
    inc
  }
}

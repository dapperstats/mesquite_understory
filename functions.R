run_model <- function(mdata, prior, reps, output = FALSE, save = TRUE){
  out <- list()
  for(i in 1:reps){
    mod <- MCMCglmm(cbind(UN, OP2, TL) ~ trait -1 +
                                       trait:mesq + trait:mesq2 + 
                                       trait:soil +
                                       trait:Burned +
                                       trait:Grazed +
                                       trait:year,
                 family = c("multinomial2", "poisson"), 
                 rcov = ~ us(trait):units, 
                 random = ~ idh(trait):Pasture +
                            idh(trait):Exclosure +
                            idh(trait):Transect + 
                            us(trait):Species + 
                            idh(trait):Functional_group +
                            idh(trait):Lifecycle +
                            idh(trait):Native + 
                            us(mesq:trait):Species + us(mesq2:trait):Species +
                            idh(soil:trait):Species +
                            idh(Burned:trait):Species + 
                            idh(Grazed:trait):Species, 
                 data = mdata, verbose = FALSE, prior = prior, pl = TRUE,
                 pr = TRUE, nitt = 200000, thin = 10, burnin = 100000)
    if(save){
      fname <- paste0("mod_run", i, ".RData")
      save(mod, file = fname)
    }
    if(output){
      out[[i]] <- mod
    }
  }
  if(output){
    out
  }
}


set_prior <- function(){
  list(G = list(G1 = list(V = diag(2), nu = 0.2),
                G2 = list(V = diag(2), nu = 0.2),
                G3 = list(V = diag(2), nu = 0.2),
                G4 = list(V = diag(2), nu = 0.2),
                G5 = list(V = diag(2), nu = 0.2),
                G6 = list(V = diag(2), nu = 0.2),
                G7 = list(V = diag(2), nu = 0.2),
                G8 = list(V = diag(2), nu = 0.2),
                G9 = list(V = diag(2), nu = 0.2),
                G10 = list(V = diag(6), nu = 0.6),
                G11 = list(V = diag(4), nu = 0.4),
                G12 = list(V = diag(4), nu = 0.4)))
}

load_and_list_runs <- function(){
  load("mod_run1.RData")
  mod1 <- mod
  load("mod_run2.RData")
  mod2 <- mod
  load("mod_run3.RData")
  mod3 <- mod
  list(run1= mod1, run2 = mod2, run3 = mod3)
}

combine_op <- function(dens){

  opend <- dens[dens$Species=="OPEN",]
  opfud <- dens[dens$Species=="OPFU",]
  opspd <- dens[dens$Species=="OPSP",]
  opendid <- paste0(opend$Pasture, opend$Exclosure, opend$Transect)
  opfudid <- paste0(opfud$Pasture, opfud$Exclosure, opfud$Transect) 
  opspdid <- paste0(opspd$Pasture, opspd$Exclosure, opspd$Transect) 

  op <- opend
  op$Species <- "OP"
  op$Species.1 <- "spp"
  op$full <- "Opuntia_spp"
  cols_to_comb <- c("YR2011_OP", "YR2011_UN", "YR2011_TL", "per_under2011",
                    "YR2014_OP", "YR2014_UN", "YR2014_TL", "per_under2014",
                    "YR2017_OP", "YR2017_UN", "YR2017_TL", "per_under2017")


  for(i in 1:nrow(op)){
    row_in1 <- which(opendid == opendid[i])
    row_in2 <- which(opendid == opfudid[i])
    row_in3 <- which(opendid == opspdid[i])
    op[row_in1, cols_to_comb] <-   op[row_in1, cols_to_comb] +
                                   op[row_in2, cols_to_comb] + 
                                   op[row_in3, cols_to_comb] 
  
  }
  op$per_under2011 <- op$YR2011_UN / op$YR2011_TL * 100
  op$per_under2014 <- op$YR2014_UN / op$YR2014_TL * 100 
  op$per_under2017 <- op$YR2017_UN / op$YR2017_TL * 100 

  op$per_under2011[is.na(op$per_under2011)] <- NA
  op$per_under2014[is.na(op$per_under2014)] <- NA
  op$per_under2017[is.na(op$per_under2017)] <- NA

  rows_out <- which(dens$Species %in% c("OPEN", "OPFU", "OPSP"))
  dens <- dens[-rows_out,]
  rbind(dens, op)
}

prep_data <- function(dens, mesq){

  dens$Lifecycle[dens$Lifecycle == "Perennail"] <- "Perennial"
  dens$Lifecycle[dens$Lifecycle == "perennial"] <- "Perennial"

  dens$Functional_group[dens$Functional_group == "shrub"] <- "Shrub"
  dens$Functional_group[dens$Functional_group == "Grass "] <- "Grass"

  densm <- dens
  densm$mesq2011 <- NA
  densm$mesq2014 <- NA
  densm$mesq2017 <- NA
  for(i in 1:nrow(densm)){
    spot <- which(mesq$Pasture == dens$Pasture[i] &
                  mesq$Exclosure == dens$Exclosure[i] &
                  mesq$Transect == dens$Transect[i])
    if (length(spot) == 1){
      densm$mesq2011[i] <- mesq$YR2011[spot]/1000
      densm$mesq2014[i] <- mesq$YR2014[spot]/1000
      densm$mesq2017[i] <- mesq$YR2017[spot]/1000
    }
  }

  densm3 <- rbind(densm, densm, densm)
  densm3$year <- rep(c(1, 4, 7), each = nrow(densm))
  densm3$mesq <- c(densm$mesq2011, densm$mesq2014, densm$mesq2017)
  densm3$mesq2 <- densm3$mesq^2
  densm3$per_under <- c(densm$per_under2011, densm$per_under2014, 
                        densm$per_under2017)/100
  densm3$OP <- round(c(densm$YR2011_OP, densm$YR2014_OP, densm$YR2017_OP))
  densm3$UN <- round(c(densm$YR2011_UN, densm$YR2014_UN, densm$YR2017_UN))
  densm3$TL <- densm3$UN + densm3$OP
  densm3$OP2 <- densm3$OP
  densm3$OP2[densm3$mesq==0] <-0

  densm3$Burned <- 0
  burnt <- which(densm3$Burn_2017 == "Yes" & densm3$year == 7)
  densm3$Burned[burnt] <- 1
  densm3$Burned <- as.factor(densm3$Burned)
  densm3$Grazed <- as.factor(as.numeric(as.factor(densm3$Grazed)) - 1)

  densm3$soil <- as.factor(densm3$soil)
  contrasts(densm3$soil) <- matrix(c(-2, 0, 1, 1, 1, -1), 
                                   nrow = 3, ncol = 2, byrow = TRUE)



  #spp_of_int <- c("BORO", "MUPO", "APTE", "SEMA", "HECO", "ERLE", "ARIS",
  #                "OP")
  densm3#[which(densm3$Species %in% spp_of_int), ]
  
}

ilogit <- function(x){
  exp(x) / (exp(x) + 1)
}

fr0 <- function(x){
  length(x[x==0])/length(x)
}

table1 <- function(dat){
  
  spco <- sort(c("BORO", "MUPO", "APTE", "SEMA", "HECO", "ERLE", "ARIS",
               "OPEN", "OPFU", "OPSP"))

  nsps <- length(spco)
  
  matr <- matrix(NA, nrow = nsps, ncol = 6)
  for(i in 1:nsps){
    op <- dat$OP2[dat$Species == spco[i]]
    un <- dat$UN[dat$Species == spco[i]]
    wts <- op + un
    xx <- un / (un + op)
    xx <- xx[wts!=0]
    wts <- wts[wts!=0]
    xxx <- dat$TL[dat$Species == spco[i]]

    matr[i,1] <- spco[i]
    matr[i,5] <- round(mean(xx), 3)
    matr[i,6] <- round(sd(xx), 3)
    matr[i,3] <- round(mean(xxx), 3)
    matr[i,4] <- round(sd(xxx), 3)
    matr[i,2] <- round(1- fr0(xxx), 3)   

  }
  tabl <- data.frame(matr)
  colnames(tabl) <- c("Species", "Presence", "Density (mean)", "Density (sd)",
                      "% under mesquite (mean)", "% under mesquite (sd)")

  write.csv(tabl, "output/table1.csv", row.names = FALSE)
}

table2 <- function(modl){

  mm <- mcmc.list(modl[[1]]$VCV[,c(42, 39, 40)], 
                  modl[[2]]$VCV[,c(42, 39, 40)], 
                  modl[[3]]$VCV[,c(42, 39, 40)])
  smm <- summary(mm)
  sizes <- effectiveSize(mm)
  psrfs <- gelman.diag(mm)[[1]][,1]

  resids <- matrix(NA, nrow = 3, ncol = 6)
  row.names(resids) <- c("Density", "Distribution", "Covariance")
  colnames(resids) <- c("Mean", "Median", "95% HPD: low", "95% HPD: up", 
                        "ESS", "psrf")
  resids[,1] <- smm$statistics[,"Mean"]
  resids[,2] <- smm$quantiles[,"50%"]
  resids[,3] <- smm$quantiles[,"2.5%"]
  resids[,4] <- smm$quantiles[,"97.5%"]
  resids[,5] <- sizes
  resids[,6] <- psrfs
  resids[,1:4] <- round(resids[,1:4], 3)
  resids[,5] <- round(resids[,5], 1)
  resids[,6] <- round(resids[,6], 3)
  write.csv(resids, "output/table2.csv")
}

table4 <- function(modl){

  mm <- mcmc.list(modl[[1]]$Sol[,1:16], modl[[2]]$Sol[,1:16], 
                  modl[[3]]$Sol[,1:16])
  smm <- summary(mm)
  sizes <- effectiveSize(mm)
  psrfs <- gelman.diag(mm)[[1]][,1]
  nF <- modl[[1]]$Fixed$nfl
  ests <- rbind(modl[[1]]$Sol[,1:nF, drop = FALSE],
                modl[[2]]$Sol[,1:nF, drop = FALSE],
                modl[[3]]$Sol[,1:nF, drop = FALSE])
  npests <- colSums(ests > 0)
  nsamps <- dim(ests)[1]
  tails <- pmin(npests/nsamps, 1 - npests/nsamps)
  ps <- 2 * pmax(0.5/nsamps, tails)

  fixeds <- matrix(NA, nrow = 18, ncol = 8)
  colnames(fixeds) <- c("Term", "Mean", "Median", "95% HPD: low", 
                        "95% HPD: up", "ESS", "psrf", "p")
  fnames <- c("Intercept", "Mesquite cover", "Mesquite cover2", 
                "Soil (sandy loam deep)", "Soil (sandy loam upland)", 
                "Burned (yes)", "Grazed (yes)", "Year")
  sp1 <- seq(2, 16, 2)
  sp2 <- seq(1, 15, 2)
  fixeds[,1] <- c("Density", fnames, "Distribution", fnames)
  fixeds[,2] <- c(NA, smm$statistics[sp1,"Mean"], 
                  NA, smm$statistics[sp2,"Mean"])
  fixeds[,3] <- c(NA, smm$quantiles[sp1,"50%"], 
                  NA, smm$quantiles[sp2,"50%"])
  fixeds[,4] <- c(NA, smm$quantiles[sp1,"2.5%"], 
                  NA, smm$quantiles[sp2,"2.5%"])
  fixeds[,5] <- c(NA, smm$quantiles[sp1,"97.5%"], 
                  NA, smm$quantiles[sp2,"97.5%"])
  fixeds[,6] <- c(NA, sizes[sp1], NA, sizes[sp2])
  fixeds[,7] <- c(NA, psrfs[sp1], NA, psrfs[sp2])  
  fixeds[,8] <- c(NA, ps[sp1], NA, ps[sp2])   

  fixeds[,2:5] <- round(as.numeric(fixeds[,2:5]), 3)
  fixeds[,6] <- round(as.numeric(fixeds[,6]), 1)
  fixeds[,7] <- round(as.numeric(fixeds[,7]), 2)
  fixeds[,8] <- round(as.numeric(fixeds[,8]), 4)
  write.csv(fixeds, "output/table4.csv", row.names = FALSE)
}

table3 <- function(mod){

  cin <- c(6, 4, 2, 5, 3, 1, 
           10, 12, 14, 16, 7, 11, 13, 15, 8,
           20, 24, 17, 21, 18, 22,
           28:30, 25:27,
           33, 34, 31, 32,
           37, 38, 35, 36)
           

  mm <- mcmc.list(modl[[1]]$VCV[,cin], 
                  modl[[2]]$VCV[,cin], 
                  modl[[3]]$VCV[,cin])

  smm <- summary(mm)
  sizes <- effectiveSize(mm)
  psrfs <- gelman.diag(mm)[[1]][,1]

  rands <- matrix(NA, nrow = length(cin), ncol = 8)
  colnames(rands) <- c("Response", "Term", "Mean", "Median", "95% HPD: low", 
                        "95% HPD: up", "ESS", "psrf")
  rtraits <- rownames(smm[[1]])
  rtraits[grepl("UN", rownames(smm[[1]]))] <- "Distribution"
  rtraits[grepl("TL", rownames(smm[[1]]))] <- "Density"
  rtraits[grepl("UN", rownames(smm[[1]])) & 
          grepl("TL", rownames(smm[[1]]))] <- "Covariance"

  rterms <- rownames(smm[[1]])
  rterms <- gsub("traitUN.", "", rterms)
  rterms <- gsub("traitTL.", "", rterms)
  rterms <- gsub("mesq:mesq", "Mesquite Cover", rterms)
  rterms <- gsub("mesq2:mesq2", "Mesquite Cover2", rterms)
  rterms <- gsub("0", " (no)", rterms)
  rterms <- gsub("1", " (yes)", rterms)
  rterms <- gsub("Functional_group", "Functional group",rterms)
  rterms <- gsub("soil", "Soil ", rterms)
  rterms <- gsub("loamy_upland", "(LU)", rterms)
  rterms <- gsub("sandy_loam_deep", "(SLD)", rterms)
  rterms <- gsub("sandy_loam_upland", "(SLU)", rterms)
  rterms <- gsub(":Species", " [slope]", rterms)

  rands[,1] <- rtraits
  rands[,2] <- rterms
  rands[,3] <- smm$statistics[,"Mean"]
  rands[,4] <- smm$quantiles[,"50%"]
  rands[,5] <- smm$quantiles[,"2.5%"]
  rands[,6] <- smm$quantiles[,"97.5%"]
  rands[,7] <- sizes
  rands[,8] <- psrfs
  rands[,3:6] <- round(as.numeric(rands[,3:6]), 3)
  rands[,7] <- round(as.numeric(rands[,7]), 1)
  rands[,8] <- round(as.numeric(rands[,8]), 3)

  write.csv(rands, "output/table3.csv", row.names = FALSE)
}


fig1 <- function(dens){

  spco <- sort(c("BORO", "MUPO", "APTE", "SEMA", "HECO", "ERLE", "ARIS",
               "OPEN", "OPFU", "OPSP"))

  sprest <- sort(c("ARIS", "BORO", "HECO", "MUPO", "SEMA"))
  spnot <- spco[!(spco %in% sprest)]

  nsp <- length(spco)
  nsprest <- length(sprest)
  nspnot <- length(spnot)

  cols <- magma(nsp, 0.6, 0.3, 0.85)

  ltys <- rep(1, nsp)
  ltys[spco %in% sprest] <- 3

  tiff("output/fig1.tif", width = 9, height = 6, units = "in", res = 400)

  par(mfrow = c(1, 2), mar = c(4,4,2,1), bty = "L")


  plot(1, 1, type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
       xlim = c(0, 150), ylim = c(0, 1))
  axis(1, at = seq(0, 150, 50))
  axis(2, at = seq(0, 1, 0.1), las = 1)
  mtext("Density of plants", side = 1, line = 2.25, cex = 1.25)
  mtext("Relative probability density", side = 2, line = 2.75, cex = 1.25)

  for(i in 1:nsps){
    xx <- dens$TL[dens$Species %in% spco[i]]
    xy <- density(xx, from = 0, to = 150, n = 100)
    xy$y <- xy$y / max(xy$y)   
    points(xy$x, xy$y, type = "l", col = cols[i], lwd = 2, lty = ltys[i]) 
  }

  ids <- paste0(dens$year, dens$Transect, dens$Exclosure, dens$Pasture)
  uids <- unique(ids)
  nuids <- length(uids)
  tdens <- rep(NA, nuids)
  udens <- rep(NA, nuids)
  for(i in 1:nuids){
    tdens[i] <- sum(dens$TL[ids == uids[i]])
    udens[i] <- sum(dens$UN[ids == uids[i]]) 
  }
  xy <- density(tdens, from = 0, to = 150, n = 100)
  xy$y <- xy$y / max(xy$y)   
  points(xy$x, xy$y, type = "l", col = 1, lwd = 3, lty = 2)

  mtext(side = 3, c(spco, "Total"), col = c(cols, 1), line = 0.75, cex = 0.9, 
        at = seq(20, 310, length.out = 1 + nsps))

  plot(1, 1, type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  axis(1, at = seq(0, 1, 0.1), labels = FALSE)
  axis(1, at = seq(0, 1, 0.5), labels = c(0, 50, 100))
  axis(2, at = seq(0, 1, 0.1), las = 1)
  mtext("Fraction under mesquite (%)", side = 1, line = 2.25, cex = 1.25)
  mtext("Relative probability density", side = 2, line = 2.75, cex = 1.25)

  for(i in 1:nsps){
    op <- dens$OP2[dens$Species == spco[i]]
    un <- dens$UN[dens$Species == spco[i]]
    wts <- op + un
    xx <- un / (un + op)
    xx <- xx[wts!=0]
    wts <- wts[wts!=0]
    xy <- suppressWarnings(
              density(xx, weights = wts, from = 0, to = 1, n = 100)
            )
    xy$y <- xy$y / max(xy$y)   
    points(xy$x, xy$y, type = "l", col = cols[i], lwd = 2, lty = ltys[i]) 
  }
  xy <- suppressWarnings(
               density(udens, weights = tdens, from = 0, to = 1, n = 100)
        )
  xy$y <- xy$y / max(xy$y)   
  points(xy$x, xy$y, type = "l", col = 1, lwd = 3, lty = 2)


  dev.off()

}


fig2 <- function(mod, mdata){


  spco <- sort(c("BORO", "MUPO", "APTE", "SEMA", "HECO", "ERLE", "ARIS",
               "OPEN", "OPFU", "OPSP"))
  nsps <- length(spco)
  tab <- matrix(NA, nrow = 4, ncol = nsps)
  vec <- rep(NA, nsps)
  spots <- list(tab, vec, vec, tab, vec, vec)
  cnames <- colnames(modl[[1]]$Sol)
  for(i in 1:nsps){

    fgr <- unique(mdata$Functional_group[mdata$Species == spco[i]])
    nat <- unique(mdata$Native[mdata$Species == spco[i]])
    lc <- unique(mdata$Lifecycle[mdata$Species == spco[i]])
    spin <- grepl(spco[i], cnames)
    spots[[1]][1,i] <- which(cnames == paste0("traitTL.Species.", spco[i]))
    spots[[1]][2,i] <- which(cnames == 
                                    paste0("traitTL.Functional_group.", fgr))
    spots[[1]][3,i] <- which(cnames == paste0("traitTL.Native.", nat))
    spots[[1]][4,i] <- which(cnames == paste0("traitTL.Lifecycle.", lc))
    spots[[2]][i] <- which(cnames == paste0("mesq:traitTL.Species.", spco[i]))
    spots[[3]][i] <- which(cnames == 
                                  paste0("mesq2:traitTL.Species.", spco[i]))
    spots[[4]][1,i] <- which(cnames == paste0("traitUN.Species.", spco[i]))
    spots[[4]][2,i] <- which(cnames == 
                                    paste0("traitUN.Functional_group.", fgr))
    spots[[4]][3,i] <- which(cnames == paste0("traitUN.Native.", nat))
    spots[[4]][4,i] <- which(cnames == paste0("traitUN.Lifecycle.", lc))
    spots[[5]][i] <- which(cnames == paste0("mesq:traitUN.Species.", spco[i]))
    spots[[6]][i] <- which(cnames == 
                                  paste0("mesq2:traitUN.Species.", spco[i]))
  }

  mm <- rbind(modl[[1]]$Sol, modl[[2]]$Sol, modl[[3]]$Sol)
  cols <- cividis(nsps, 0.8, 0.1, 0.9)

  tiff("output/fig2.tif", width = 9, height = 15, units = "in", res = 400)
  par(mar = c(4, 4, 2, 1.5), mfrow = c(3, 2))

  st <- list()
  for(i in 1:10){
    st[[i]] <- as.numeric(mm[,spots[[1]][1,i]]) +
               as.numeric(mm[,spots[[1]][2,i]]) +
               as.numeric(mm[,spots[[1]][3,i]]) +
               as.numeric(mm[,spots[[1]][4,i]])
  }
  boxplot(st, outline = FALSE, ylim = c(-12, 12), col = cols,
          xlab = "", xaxt = "n", ylab = "", yaxt = "n")

  mtext(side = 3, "Density, Intercept", line = 0.25, cex = 1.1)
  labs <- spco
  axis(side = 1, labels = FALSE, at = 1:10)
  axis(side = 2, at = seq(-10, 10, 5), las = 1, cex.axis = 1.25)
  text(1:10, -14, labs, cex = 1.2, srt = 30, xpd = TRUE, adj = 1)


  st <- list()
  for(i in 1:10){
    st[[i]] <- as.numeric(mm[,spots[[4]][1,i]]) +
               as.numeric(mm[,spots[[4]][2,i]]) +
               as.numeric(mm[,spots[[4]][3,i]]) +
               as.numeric(mm[,spots[[4]][4,i]])
  }
  boxplot(st, outline = FALSE, ylim = c(-10, 10), col = cols,
          xlab = "", xaxt = "n", ylab = "", yaxt = "n")

  labs <- spco
  mtext(side = 3, "Distribution, Intercept", line = 0.25, cex = 1.1)
  axis(side = 1, labels = FALSE, at = 1:10)
  axis(side = 2, at = seq(-10, 10, 5), las = 1, cex.axis = 1.25)
  text(1:10, rep(-11.7, 8), labs, cex = 1.2, srt = 30, xpd = TRUE, adj = 1)


  st <- list()
  for(i in 1:10){
    st[[i]] <- as.numeric(mm[,spots[[2]][i]])
  }
  boxplot(st, outline = FALSE, ylim = c(-25, 25), col = cols,
          xlab = "", xaxt = "n", ylab = "", yaxt = "n")
  mtext(side = 3, "Density, Mesquite Cover", line = 0.25, cex = 1.1)
  axis(side = 1, labels = FALSE, at = 1:10)
  axis(side = 2, at = seq(-20, 20, 10), las = 1, cex.axis = 1.25)
  text(1:10, rep(-29.2, 8), labs, cex = 1.2, srt = 30, xpd = TRUE, adj = 1)


  st <- list()
  for(i in 1:10){
    st[[i]] <- as.numeric(mm[,spots[[5]][i]])
  }
  boxplot(st, outline = FALSE, ylim = c(-8, 8), col = cols,
          xlab = "", xaxt = "n", ylab = "", yaxt = "n")

  mtext(side = 3, "Distribution, Mesquite Cover", line = 0.25, cex = 1.1)
  axis(side = 1, labels = FALSE, at = 1:10)
  axis(side = 2, at = seq(-8, 8, 4), las = 1, cex.axis = 1.25)
  text(1:10, rep(-9.3, 8), labs, cex = 1.2, srt = 30, xpd = TRUE, adj = 1)


  st <- list()
  for(i in 1:10){
    st[[i]] <- as.numeric(mm[,spots[[3]][i]])
  }
  boxplot(st, outline = FALSE, ylim = c(-45, 45), col = cols,
          xlab = "", xaxt = "n", ylab = "", yaxt = "n")

  mtext(side = 3, expression(paste("Density, Mesquite ", Cover^2 )), 
        cex = 1.1)
  axis(side = 1, labels = FALSE, at = 1:10)
  axis(side = 2, at = seq(-40, 40, 20), las = 1, cex.axis = 1.25)
  text(1:10, rep(-52.5, 8), labs, cex = 1.2, srt = 30, xpd = TRUE, adj = 1)


  st <- list()
  for(i in 1:10){
    st[[i]] <- as.numeric(mm[,spots[[6]][i]])
  }
  boxplot(st, outline = FALSE, ylim = c(-20, 20), col = cols,
          xlab = "", xaxt = "n", ylab = "", yaxt = "n")

  mtext(side = 3, expression(paste("Distribution, Mesquite ", Cover^2 )), 
        cex = 1.1)
  axis(side = 1, labels = FALSE, at = 1:10)
  axis(side = 2, at = seq(-20, 20, 10), las = 1, cex.axis = 1.25)
  text(1:10, rep(-23.3, 8), labs, cex = 1.2, srt = 30, xpd = TRUE, adj = 1)


  dev.off()
}

fig3 <- function(mod, mdata){


  spco <- sort(c("BORO", "MUPO", "APTE", "SEMA", "HECO", "ERLE", "ARIS",
               "OPEN", "OPFU", "OPSP"))
  nsps <- length(spco)
  tab <- matrix(NA, nrow = 4, ncol = nsps)
  vec <- rep(NA, nsps)
  spots <- list(tab, vec, vec, tab, vec, vec)
  cnames <- colnames(mod[[1]]$Sol)
  for(i in 1:nsps){

    fgr <- unique(mdata$Functional_group[mdata$Species == spco[i]])
    nat <- unique(mdata$Native[mdata$Species == spco[i]])
    lc <- unique(mdata$Lifecycle[mdata$Species == spco[i]])
    spin <- grepl(spco[i], cnames)
    spots[[1]][1,i] <- which(cnames == paste0("traitTL.Species.", spco[i]))
    spots[[1]][2,i] <- which(cnames == 
                                    paste0("traitTL.Functional_group.", fgr))
    spots[[1]][3,i] <- which(cnames == paste0("traitTL.Native.", nat))
    spots[[1]][4,i] <- which(cnames == paste0("traitTL.Lifecycle.", lc))
    spots[[2]][i] <- which(cnames == paste0("mesq:traitTL.Species.", spco[i]))
    spots[[3]][i] <- which(cnames == 
                                  paste0("mesq2:traitTL.Species.", spco[i]))
    spots[[4]][1,i] <- which(cnames == paste0("traitUN.Species.", spco[i]))
    spots[[4]][2,i] <- which(cnames == 
                                    paste0("traitUN.Functional_group.", fgr))
    spots[[4]][3,i] <- which(cnames == paste0("traitUN.Native.", nat))
    spots[[4]][4,i] <- which(cnames == paste0("traitUN.Lifecycle.", lc))
    spots[[5]][i] <- which(cnames == paste0("mesq:traitUN.Species.", spco[i]))
    spots[[6]][i] <- which(cnames == 
                                  paste0("mesq2:traitUN.Species.", spco[i]))
  }


  spco2 <- unique(mdata$Species)
  nsps2 <- length(spco2)
  tab <- matrix(NA, nrow = 4, ncol = nsps2)
  vec <- rep(NA, nsps2)
  spots2 <- list(tab, vec, vec, tab, vec, vec)
  for(i in 1:nsps2){

    fgr <- unique(mdata$Functional_group[mdata$Species == spco2[i]])
    nat <- unique(mdata$Native[mdata$Species == spco2[i]])
    lc <- unique(mdata$Lifecycle[mdata$Species == spco2[i]])
    spots2[[1]][1,i] <- which(cnames == paste0("traitTL.Species.", spco2[i]))
    spots2[[1]][2,i] <- which(cnames == 
                                    paste0("traitTL.Functional_group.", fgr))
    spots2[[1]][3,i] <- which(cnames == paste0("traitTL.Native.", nat))
    spots2[[1]][4,i] <- which(cnames == paste0("traitTL.Lifecycle.", lc))
    spots2[[2]][i] <- which(cnames == 
                                  paste0("mesq:traitTL.Species.", spco2[i]))
    spots2[[3]][i] <- which(cnames == 
                                  paste0("mesq2:traitTL.Species.", spco2[i]))
    spots2[[4]][1,i] <- which(cnames == paste0("traitUN.Species.", spco2[i]))
    spots2[[4]][2,i] <- which(cnames == 
                                    paste0("traitUN.Functional_group.", fgr))
    spots2[[4]][3,i] <- which(cnames == paste0("traitUN.Native.", nat))
    spots2[[4]][4,i] <- which(cnames == paste0("traitUN.Lifecycle.", lc))
    spots2[[5]][i] <- which(cnames == 
                                 paste0("mesq:traitUN.Species.", spco2[i]))
    spots2[[6]][i] <- which(cnames == 
                                  paste0("mesq2:traitUN.Species.", spco2[i]))
  }


  x <- seq(0, 0.56, 0.001)
  mm <- rbind(mod[[1]]$Sol, mod[[2]]$Sol, mod[[3]]$Sol)
  nsamps <- nrow(mm)
 
  spco <- sort(c("BORO", "MUPO", "APTE", "SEMA", "HECO", "ERLE", "ARIS",
               "OPEN", "OPFU", "OPSP"))

  sprest <- sort(c("ARIS", "BORO", "HECO", "MUPO", "SEMA"))
  spnot <- spco[!(spco %in% sprest)]

  nsp <- length(spco)
  nsprest <- length(sprest)
  nspnot <- length(spnot)

  cols <- magma(nsp, 0.6, 0.3, 0.85)

  ltys <- rep(1, nsp)
  ltys[spco %in% sprest] <- 3

  tiff("output/fig3.tif", width = 9, height = 6, units = "in", res = 400)


  par(mfrow = c(1,2), mar = c(4, 4, 2, 1), bty = "L")

  plot(1, 1, type = "n", ylim = c(0, 5.5), xlim = c(0,0.55), xaxt = "n", 
       yaxt = "n", xlab = "", ylab = "")
  for(i in 1:10){
    ymat <- matrix(NA, nrow = nsamps, ncol = length(x))
    for(j in 1:nsamps){
      y <- (mm[j,2] + mm[j,spots[[1]][1,i]] + mm[j,spots[[1]][2,i]] +
                mm[j,spots[[1]][3,i]] + mm[j,spots[[1]][4,i]]) + 
           (mm[j,4] + mm[j,spots[[2]][i]]) * x + 
           (mm[j,6] + mm[j,spots[[3]][i]]) * x^2 
      ymat[j,] <- exp(y)
    }
    points(x, apply(ymat, 2, mean), type = "l", col = cols[i], lwd = 2,
           lty = ltys[i])
  }
  yymat <- matrix(NA, nrow = nsps2, ncol = length(x))  
  for(i in 1:nsps2){
    ymat <- matrix(NA, nrow = nsamps, ncol = length(x))
    for(j in 1:nsamps){
      y <- (mm[j,2] + mm[j,spots2[[1]][1,i]] + mm[j,spots2[[1]][2,i]] +
                mm[j,spots2[[1]][3,i]] + mm[j,spots2[[1]][4,i]]) + 
           (mm[j,4] + mm[j,spots2[[2]][i]]) * x + 
           (mm[j,6] + mm[j,spots2[[3]][i]]) * x^2 
      ymat[j,] <- exp(y)
    }
    yymat[i,] <- (apply(ymat, 2, mean))
  }
  points(x, apply(yymat, 2, sum), type = "l", col = 1, lwd = 3, lty = 2)

  axis(1, at = seq(0, 0.5, 0.1), labels = seq(0, 50, 10), las = 1)
  axis(2, at = seq(0, 5, 1), las = 1)

  mtext(side = 1, "Mesquite cover (%)", line = 2.5, cex = 1.25)
  mtext(side = 2, "Density (number of plants)", line = 2.5, cex = 1.25)



  plot(1, 1, type = "n", ylim = c(0, 1), xlim = c(0,0.55), xaxt = "n", 
       yaxt = "n", xlab = "", ylab = "")
  for(i in 1:10){
    ymat <- matrix(NA, nrow = nsamps, ncol = length(x))
    for(j in 1:nsamps){
      y <- (mm[j,1] + mm[j,spots[[4]][1,i]] + mm[j,spots[[4]][2,i]] +
                mm[j,spots[[4]][3,i]] + mm[j,spots[[4]][4,i]]) + 
           (mm[j,3] + mm[j,spots[[5]][i]]) * x + 
           (mm[j,5] + mm[j,spots[[6]][i]]) * x^2 
      ymat[j,] <- ilogit(y)
    }
    points(x, apply(ymat, 2, mean), type = "l", col = cols[i], lwd = 2,
           lty = ltys[i])
  }

  yymat1 <- array(NA, dim = c(nsamps, length(x)))
  yymat2 <- array(NA, dim = c(nsamps, length(x)))

  for(j in 1:nsamps){
    tmat1 <- matrix(NA, nrow = nsps2, ncol = length(x))
    tmat2 <- matrix(NA, nrow = nsps2, ncol = length(x))
    for(i in 1:nsps2){
      y <- (mm[j,2] + mm[j,spots2[[1]][1,i]] + mm[j,spots2[[1]][2,i]] +
                mm[j,spots2[[1]][3,i]] + mm[j,spots2[[1]][4,i]]) + 
           (mm[j,4] + mm[j,spots2[[2]][i]]) * x + 
           (mm[j,6] + mm[j,spots2[[3]][i]]) * x^2 
      y2 <- (mm[j,1] + mm[j,spots2[[4]][1,i]] + mm[j,spots2[[4]][2,i]] +
                mm[j,spots2[[4]][3,i]] + mm[j,spots2[[4]][4,i]]) + 
           (mm[j,3] + mm[j,spots2[[5]][i]]) * x + 
           (mm[j,5] + mm[j,spots2[[6]][i]]) * x^2 
      tmat1[i,] <- exp(y)
      tmat2[i,] <- exp(y) * ilogit(y2)

    }
    yymat1[j,] <- apply(tmat1, 2, sum)
    yymat2[j,] <- apply(tmat2, 2, sum)
  }
  yymat3 <- yymat2/yymat1
  yyt <- apply(yymat3, 2, mean)
  points(x, yyt, type = "l", col = 1, lwd = 3, lty = 2)

  axis(1, at = seq(0, 0.5, 0.1), labels = seq(0, 50, 10), las = 1)
  axis(2, at = seq(0, 1, 0.1), labels = FALSE, las = 1)
  axis(2, at = seq(0, 1, 0.5), labels = seq(0, 100, 50), las = 1)

  mtext(side = 1, "Mesquite cover (%)", line = 2.5, cex = 1.25)
  mtext(side = 2, "Percentage under mesquite (%)", line = 2.5, cex = 1.25)

  mtext(side = 3, c(spco, "Total"), col = c(cols, 1), line = 0.75, cex = 0.9, 
        at = seq(-0.8, 0.45, length.out = nsps+1))

  dev.off() 
}

table5 <- function(mod, mdata){

  spco <- sort(c("BORO", "MUPO", "APTE", "SEMA", "HECO", "ERLE", "ARIS",
               "OPEN", "OPFU", "OPSP"))
  nsps <- length(spco)
  tab <- matrix(NA, nrow = 4, ncol = nsps)
  vec <- rep(NA, nsps)
  spots <- list(tab, vec, vec, tab, vec, vec)
  cnames <- colnames(modl[[1]]$Sol)
  for(i in 1:nsps){

    fgr <- unique(mdata$Functional_group[mdata$Species == spco[i]])
    nat <- unique(mdata$Native[mdata$Species == spco[i]])
    lc <- unique(mdata$Lifecycle[mdata$Species == spco[i]])
    spin <- grepl(spco[i], cnames)
    spots[[1]][1,i] <- which(cnames == paste0("traitTL.Species.", spco[i]))
    spots[[1]][2,i] <- which(cnames == 
                                    paste0("traitTL.Functional_group.", fgr))
    spots[[1]][3,i] <- which(cnames == paste0("traitTL.Native.", nat))
    spots[[1]][4,i] <- which(cnames == paste0("traitTL.Lifecycle.", lc))
    spots[[2]][i] <- which(cnames == paste0("mesq:traitTL.Species.", spco[i]))
    spots[[3]][i] <- which(cnames == 
                                  paste0("mesq2:traitTL.Species.", spco[i]))
    spots[[4]][1,i] <- which(cnames == paste0("traitUN.Species.", spco[i]))
    spots[[4]][2,i] <- which(cnames == 
                                    paste0("traitUN.Functional_group.", fgr))
    spots[[4]][3,i] <- which(cnames == paste0("traitUN.Native.", nat))
    spots[[4]][4,i] <- which(cnames == paste0("traitUN.Lifecycle.", lc))
    spots[[5]][i] <- which(cnames == paste0("mesq:traitUN.Species.", spco[i]))
    spots[[6]][i] <- which(cnames == 
                                  paste0("mesq2:traitUN.Species.", spco[i]))
  }


  spco2 <- unique(mdata$Species)
  nsps2 <- length(spco2)
  tab <- matrix(NA, nrow = 4, ncol = nsps2)
  vec <- rep(NA, nsps2)
  spots2 <- list(tab, vec, vec, tab, vec, vec)
  for(i in 1:nsps2){

    fgr <- unique(mdata$Functional_group[mdata$Species == spco2[i]])
    nat <- unique(mdata$Native[mdata$Species == spco2[i]])
    lc <- unique(mdata$Lifecycle[mdata$Species == spco2[i]])
    spots2[[1]][1,i] <- which(cnames == paste0("traitTL.Species.", spco2[i]))
    spots2[[1]][2,i] <- which(cnames == 
                                    paste0("traitTL.Functional_group.", fgr))
    spots2[[1]][3,i] <- which(cnames == paste0("traitTL.Native.", nat))
    spots2[[1]][4,i] <- which(cnames == paste0("traitTL.Lifecycle.", lc))
    spots2[[2]][i] <- which(cnames == 
                                  paste0("mesq:traitTL.Species.", spco2[i]))
    spots2[[3]][i] <- which(cnames == 
                                  paste0("mesq2:traitTL.Species.", spco2[i]))
    spots2[[4]][1,i] <- which(cnames == paste0("traitUN.Species.", spco2[i]))
    spots2[[4]][2,i] <- which(cnames == 
                                    paste0("traitUN.Functional_group.", fgr))
    spots2[[4]][3,i] <- which(cnames == paste0("traitUN.Native.", nat))
    spots2[[4]][4,i] <- which(cnames == paste0("traitUN.Lifecycle.", lc))
    spots2[[5]][i] <- which(cnames == 
                                 paste0("mesq:traitUN.Species.", spco2[i]))
    spots2[[6]][i] <- which(cnames == 
                                  paste0("mesq2:traitUN.Species.", spco2[i]))
  }


  x <- seq(0, 0.56, 0.001)
  mm <- rbind(modl[[1]]$Sol, modl[[2]]$Sol, modl[[3]]$Sol)
  nsamps <- nrow(mm)

  out <- matrix(NA, nrow = 13, ncol = 4)
  colnames(out) <- c("Taxon", "Max slope (mean [95% HPD])",
                     "Max % under (mean [95% HPD])", 
                     "Mesquite cover at max % under (mean [95% HPD])")
  out[,1] <- c(spco, "Native", "Non-Native", "Total")

  for(i in 1:10){
    ymat <- matrix(NA, nrow = nsamps, ncol = length(x))
    ydmat <- matrix(NA, nrow = nsamps, ncol = length(x) - 1)
    for(j in 1:nsamps){
      y <- (mm[j,1] + mm[j,spots[[4]][1,i]] + mm[j,spots[[4]][2,i]] +
                mm[j,spots[[4]][3,i]] + mm[j,spots[[4]][4,i]]) + 
           (mm[j,3] + mm[j,spots[[5]][i]]) * x + 
           (mm[j,5] + mm[j,spots[[6]][i]]) * x^2 
      ymat[j,] <- ilogit(y)
      ydmat[j,] <- diff(ilogit(y))
    }
    ms <- apply(ydmat, 1, max) *1000 
    m1 <- apply(ymat, 1, max)
    wm1 <- x[apply(ymat, 1, which.max)]
    out[i,2] <- paste0(round(mean(ms), 3),
                ", [",
                 paste(round(HPDinterval(as.mcmc(ms)), 4), collapse = ", "),
                 "]")
    out[i,3] <- paste0(round(mean(m1), 3),
                ", [",
                 paste(round(HPDinterval(as.mcmc(m1)), 4), collapse = ", "),
                 "]")
    out[i,4] <- paste0(round(mean(wm1), 3), 
                ", [",
                 paste(round(HPDinterval(as.mcmc(wm1)), 4), collapse = ", "),
                 "]")
  }

  yymat1 <- array(NA, dim = c(nsamps, length(x)))
  yymat2 <- array(NA, dim = c(nsamps, length(x)))
  yymat1n <- array(NA, dim = c(nsamps, length(x)))
  yymat2n <- array(NA, dim = c(nsamps, length(x)))
  yymat1y <- array(NA, dim = c(nsamps, length(x)))
  yymat2y <- array(NA, dim = c(nsamps, length(x)))
  nn <- which(spco2 %in% unique(mdata$Species[mdata$Native == "N"]))
  yn <- which(spco2 %in% unique(mdata$Species[mdata$Native == "Y"])) 
  for(j in 1:nsamps){
    tmat1 <- matrix(NA, nrow = nsps2, ncol = length(x))
    tmat2 <- matrix(NA, nrow = nsps2, ncol = length(x))
    for(i in 1:nsps2){
      y <- (mm[j,2] + mm[j,spots2[[1]][1,i]] + mm[j,spots2[[1]][2,i]] +
                mm[j,spots2[[1]][3,i]] + mm[j,spots2[[1]][4,i]]) + 
           (mm[j,4] + mm[j,spots2[[2]][i]]) * x + 
           (mm[j,6] + mm[j,spots2[[3]][i]]) * x^2 
      y2 <- (mm[j,1] + mm[j,spots2[[4]][1,i]] + mm[j,spots2[[4]][2,i]] +
                mm[j,spots2[[4]][3,i]] + mm[j,spots2[[4]][4,i]]) + 
           (mm[j,3] + mm[j,spots2[[5]][i]]) * x + 
           (mm[j,5] + mm[j,spots2[[6]][i]]) * x^2 
      tmat1[i,] <- exp(y)
      tmat2[i,] <- exp(y) * ilogit(y2)

    }
    yymat1[j,] <- apply(tmat1, 2, sum)
    yymat2[j,] <- apply(tmat2, 2, sum)
    yymat1n[j,] <- apply(tmat1[nn,], 2, sum)
    yymat2n[j,] <- apply(tmat2[nn,], 2, sum)
    yymat1y[j,] <- apply(tmat1[yn,], 2, sum)
    yymat2y[j,] <- apply(tmat2[yn,], 2, sum)
  }
  yymat3 <- yymat2/yymat1
  yymat3y <- yymat2y/yymat1y
  yymat3n <- yymat2n/yymat1n
  yydmat3 <- matrix(NA, nrow = nrow(yymat3), ncol = ncol(yymat3) - 1)
  yydmat3y <- matrix(NA, nrow = nrow(yymat3), ncol = ncol(yymat3) - 1)
  yydmat3n <- matrix(NA, nrow = nrow(yymat3), ncol = ncol(yymat3) - 1)
  for(k in 1:ncol(yydmat3)){
    yydmat3[,k] <- yymat3[,k+1] - yymat3[,k]
    yydmat3y[,k] <- yymat3y[,k+1] - yymat3y[,k] 
    yydmat3n[,k] <- yymat3n[,k+1] - yymat3n[,k]
  }

  m1 <- apply(yymat3y, 1, max)
  wm1 <- x[apply(yymat3y, 1, which.max)]
  ms <- apply(yydmat3y, 1, max) *1000 

  out[11,2] <- paste0(round(mean(ms), 3),
                ", [",
                 paste(round(HPDinterval(as.mcmc(ms)), 4), collapse = ", "),
                 "]")
  out[11,3] <- paste0(round(mean(m1), 3),
                ", [",
                 paste(round(HPDinterval(as.mcmc(m1)), 4), collapse = ", "),
                 "]")
  out[11,4] <- paste0(round(mean(wm1), 3), 
                ", [",
                 paste(round(HPDinterval(as.mcmc(wm1)), 4), collapse = ", "),
                 "]")

  m1 <- apply(yymat3n, 1, max)
  wm1 <- x[apply(yymat3n, 1, which.max)]
  ms <- apply(yydmat3n, 1, max) *1000 

  out[12,2] <- paste0(round(mean(ms), 3),
                ", [",
                 paste(round(HPDinterval(as.mcmc(ms)), 4), collapse = ", "),
                 "]")
  out[12,3] <- paste0(round(mean(m1), 3),
                ", [",
                 paste(round(HPDinterval(as.mcmc(m1)), 4), collapse = ", "),
                 "]")
  out[12,4] <- paste0(round(mean(wm1), 3), 
                ", [",
                 paste(round(HPDinterval(as.mcmc(wm1)), 4), collapse = ", "),
                 "]")


  m1 <- apply(yymat3, 1, max)
  wm1 <- x[apply(yymat3, 1, which.max)]
  ms <- apply(yydmat3, 1, max) *1000 

  out[13,2] <- paste0(round(mean(ms), 3),
                ", [",
                 paste(round(HPDinterval(as.mcmc(ms)), 4), collapse = ", "),
                 "]")
  out[13,3] <- paste0(round(mean(m1), 3),
                ", [",
                 paste(round(HPDinterval(as.mcmc(m1)), 4), collapse = ", "),
                 "]")
  out[13,4] <- paste0(round(mean(wm1), 3), 
                ", [",
                 paste(round(HPDinterval(as.mcmc(wm1)), 4), collapse = ", "),
                 "]")

  write.csv(out, "output/table5.csv", row.names = FALSE)

}


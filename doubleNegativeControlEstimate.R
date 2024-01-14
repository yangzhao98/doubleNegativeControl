library(splines)
library(sandwich)
# - Prepare the data using moving average method
MovingAvg <- function(X, dat, magT) {
  # X: character, variable name of exposure in the data set
  # dat: original data set 
  # maT: integer, time-period for moving averaging 
  nameList <- names(dat)
  lnth <- nrow(dat)
  datNew <- dat[(magT+1):lnth,]
  datNew$X <- unlist(lapply(1:(lnth-magT), function(i) mean(dat[[X]][i:(i+magT)], na.rm = TRUE)))
  names(datNew) <- c(nameList, paste(X, "magT", magT, sep = ""))
  return(datNew)
}

# - Prepare the data using different lag structure
lagStructure <- function(X = "o3", dat = datExample, lagT = 1) {
  # X: character, variable name of exposure in the data set
  # dat: original data set 
  # lagT: integer, time-period for moving averaging 
  nameList <- names(dat)
  lnth <- nrow(dat)
  datNew <- dat[(lagT+1):lnth,]
  datNew$X <- unlist(lapply(1:(lnth-lagT), function(i) dat[[X]][i]))
  names(datNew) <- c(nameList, paste(X, "lagT", lagT, sep = ""))
  return(datNew)
}

# - Construct data for performing double negative control analysis
constructNegativeControlDat <- function(Y, X, C, lagNegCont) {
  
  # Original data
  Y0 <- Y                                # vector, [Outcome]
  X0 <- X                                # vector, [Exposure]
  V0 <- C                                # matrix, [Observed confounders]
  lagNegCont <- lagNegCont               # No. of lag for constructing negative controls
  
  # Reconstruct the index
  lnth  <- length(Y0)                    # No. of sample size
  lnthW <- 1:(lnth - lagNegCont - 1)     # Negative control outcome (lag lagNegCont, backward, Y_i-lagNegCont)
  lnthY <- (lagNegCont + 1):(lnth - 1)   # Outcome
  lnthZ <- (lagNegCont + 2):lnth         # Negative control exposure (lag lagNegCont, forward, X_i+lagNegCont)

  # Reconstruct the data
  Y  <- Y0[lnthY]                        # [Outcome]
  yX <- X0[lnthY]                        # [Exposure]
  W  <- Y0[lnthW]                        # [Negative control outcome] with lagged-k (backward, Y_i-k) 
  wX <- X0[lnthW]                        # [Exposure] for negative control outcome
  Z  <- X0[lnthZ]                        # [Negative control exposure] with lagged-k (forward, X_i+k)
  if (!is.matrix(V0)) {
    yV <- V0[lnthY]                      # [Observed confounders] for outcome
    wV <- V0[lnthW]                      # [Observed confounders] for negative control outcome
  } else {
    yV <- V0[lnthY,]
    wV <- V0[lnthW,]
    if(!is.null(colnames(C))) {
      colnames(yV) <- paste("yV", colnames(C), sep = "")
      colnames(wV) <- paste("wV", colnames(C), sep = "")
    } else {
      colnames(yV) <- paste("yV", 1:ncol(C), sep = "")
      colnames(wV) <- paste("wV", 1:ncol(C), sep = "")
    }
  }
  return(as.data.frame(cbind(Y, yX, yV, W, wX, wV, Z)))
}

# - Perform double negative control analysis using two-stage least estimator
doubleNegativeControlEstimate <- function(formulaDNC, formulaCT, dat) {
  
  # Double negative control estimator using two-stage least square estimator via AER::ivreg()
  fitModelDNC <- as.formula(formulaDNC)  # Model for double negative control estimator
  aerACE <- AER::ivreg(fitModelDNC, data = dat)
  ACE <- as.data.frame(cbind(lmtest::coeftest(aerACE, vcov = vcovHAC, df = Inf),
                             lmtest::coefci(aerACE, vcov = vcovHAC, df = Inf)))
  vcovACE <- sandwich::vcovHAC(aerACE)
  ## Confounding test using least square estimator
  fitModelCT <- as.formula(formulaCT)   # Model for confounding test
  confoundingtest <- stats::lm(fitModelCT, data = dat)
  CT <- as.data.frame(cbind(lmtest::coeftest(confoundingtest, vcov = vcovHAC),
                            lmtest::coefci(confoundingtest, vcov = vcovHAC)))
  CT <- CT[rownames(CT) == "Z",]
  return(list(doubleNC = ACE, vcovACE = vcovACE, confoundingEffect = CT))
}

# ---------------------------------------------------------------------------------------------
# - Corss-checking the previous functions using the simulated data in the original paper
# ---------------------------------------------------------------------------------------------
# Constructing negative control outcome and exposure with a pre-specified lagged time period
# datDNC <- constructNegativeControlDat(Y = Y0, X = X0, C = V0, lagNegCont = 1)
# 
# Obtain the double negative control estimator
# uACE <- doubleNegativeControlEstimate(Y ~ W + yX + yV + wX + wV | yX + yV + wX + wV + Z,
#                                       W ~ wX + wV + Z,
#                                       dat = datDNC)
# ---------------------------------------------------------------------------------------------
# 
# # - Reanalyze the trial data in the 'dlnm' package
# # create the crossbasis objects and summarize their contents
# cb1.pm <- crossbasis(chicagoNMMAPS$pm10, lag=15, argvar=list(fun="lin"),
#                      arglag=list(fun="poly",degree=4))
# cb1.temp <- crossbasis(chicagoNMMAPS$temp, lag=3, argvar=list(df=5),
#                        arglag=list(fun="strata",breaks=1))
# summary(cb1.pm)
# summary(cb1.temp)
# library(splines)
# model1 <- glm(death ~ cb1.pm + cb1.temp + ns(time, 7*14) + dow,
#               family=quasipoisson(), chicagoNMMAPS)
# pred1.pm <- crosspred(cb1.pm, model1, at=0:20, bylag=0.2, cumul=TRUE)
# 
# # plot the lag-response curves for specific and incremental cumulative effects
# plot(pred1.pm, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
#      main="Lag-response curve for a 10-unit increase in PM10")
# plot(pred1.pm, "slices", var=10, col=2, cumul=TRUE, ylab="Cumulative RR",
#      main="Lag-response curve of incremental cumulative effects")
# 
# 
# library(tidyverse)
# library(data.table)
# library(lubridate)
# datExample <- chicagoNMMAPS %>%
#   mutate(week = week(ymd(date))) %>%
#   mutate(weekIdx = paste(year,"-Week ", week, sep = "")) %>%
#   group_by(weekIdx) %>%
#   summarise(date = first(date),
#             year = first(year),
#             month = first(month),
#             death = sum(death),
#             cvd = sum(cvd),
#             resp = sum(resp),
#             temp = mean(temp, na.rm = TRUE),
#             dptp = mean(dptp, na.rm = TRUE),
#             rhum = mean(rhum, na.rm = TRUE),
#             pm10 = mean(pm10, na.rm = TRUE),
#             o3 = mean(o3, na.rm = TRUE)) %>%
#   ungroup(weekIdx) %>%
#   mutate(time = 1:n())
# 
# 
# 
# Y <- sqrt(datExample$death)
# X <- datExample$o3
# C <- as.matrix(datExample[,c("time","temp", "pm10")])
# 
# datDNC_eg <- constructNegativeControlDat(Y = Y, X = X, C = C, lagNegCont = 1)
# datDNC_eg <- datDNC_eg[complete.cases(datDNC_eg),]
# 
# # Main single-pollutant model
# spm_uACE <- doubleNegativeControlEstimate(Y ~ W + yX + yVtemp + ns(yVtime,7) + wX + wVtemp + ns(wVtime,7) | 
#                                             yX + yVtemp + ns(yVtime,7) + wX + wVtemp + ns(wVtime,7) + Z,
#                                           W ~ wX + wVtemp + ns(wVtime,7) + Z,
#                                           dat = datDNC_eg)
# 
# # Main two-pollutant model
# tpm_uACE <- doubleNegativeControlEstimate(Y ~ W + yX + yVtemp + ns(yVtime,7) + yVpm10 + wX + wVtemp + ns(wVtime,7) + wVpm10 | 
#                                             yX + yVtemp + ns(yVtime,7) + yVpm10 + wX + wVtemp + ns(wVtime,7) + wVpm10 + Z,
#                                           W ~ wX + wVtemp + ns(wVtime,7) + wVpm10 + Z,
#                                           dat = datDNC_eg)
# 
# 
# # Non-linear effect, concentration-response curve
# yXknots <- quantile(datDNC_eg$yX, probs = c(0.25, 0.75))
# wXknots <- quantile(datDNC_eg$wX, probs = c(0.25, 0.75))
# 
# nl_uACE <- doubleNegativeControlEstimate(Y ~ W + bs(yX, df = 3, knots = yXknots) + yVtemp + ns(yVtime,7) + yVpm10 + bs(wX, df = 3, knots = wXknots) + wVtemp + ns(wVtime,7) + wVpm10 | 
#                                            bs(yX, df = 3, knots = yXknots) + yVtemp + ns(yVtime,7) + yVpm10 + bs(wX, df = 3, knots = wXknots) + wVtemp + ns(wVtime,7) + wVpm10 + Z,
#                                          W ~ bs(wX, df = 3, knots = wXknots) + wVtemp + ns(wVtime,7) + wVpm10 + Z,
#                                          dat = datDNC_eg)
# 
# covMat <- nl_uACE$vcovACE
# covMat <- covMat[rownames(covMat) %in% paste("bs(yX, df = 3, knots = yXknots)", 1:5, sep = ""),
#                  colnames(covMat) %in% paste("bs(yX, df = 3, knots = yXknots)", 1:5, sep = "")]
# coef <- nl_uACE$doubleNC
# coef <- coef[rownames(coef) %in% paste("bs(yX, df = 3, knots = yXknots)", 1:5, sep = ""), "Estimate"]
# 
# yXGrid <- seq(min(datDNC_eg$yX), max(datDNC_eg$yX), length = 100)
# yXpred <- onebasis(datDNC_eg$yX, "bs", degree = 3, knots = yXknots)
# # byX <- bs(datDNC$yX, df = 3, knots = yXknots) 
# yXPredE <- crosspred(yXpred, coef = coef,  vcov = covMat, at = yXGrid, cen = 10)
# dat_crcurve <- as.data.frame(cbind(o3 = yXGrid, 
#                                    ACE = yXPredE$allfit, 
#                                    ACELCI = yXPredE$alllow, 
#                                    ACEUCI = yXPredE$allhigh))
# rownames(dat_crcurve) <- NULL
# # plot(yXPredE)
# library(ggplot2)
# pCRCurve <- ggplot(data = dat_crcurve, aes(x = o3)) +
#   geom_line(aes(y = ACE), linetype = 1, size = 1) +
#   geom_line(aes(y = ACELCI), linetype = 2) +
#   geom_line(aes(y = ACEUCI), linetype = 2) +
#   theme_classic() +
#   scale_x_continuous(
#     name = expression(paste(O^{3}, " Concentration (", mu ,"/", m^{3},")", sep = "")),
#     limits = c(0,55),
#     breaks = c(0,10,20,30,40,50,60),
#     labels = c("0", "10", "20", "30", "40", "50", "60"),
#     expand = c(0,0)) +
#   scale_y_continuous(
#     name = "Absolute difference in morbidity",
#     limits = c(-22,32),
#     breaks = c(-20,-10,0,10,20,30),
#     expand = c(0,0)) +
#   geom_hline(aes(yintercept = 0), color = "gray")
# pCRCurve
# 
# 
# # Delayed effects
# # S1. Moving average 
# datExamplemaT1 <- MovingAvg(X = "o3", dat = datExample, maT = 1)  
# Y <- sqrt(datExamplemaT1$death)
# X <- datExamplemaT1$o3maT1
# C <- as.matrix(datExamplemaT1[,c("time","temp")])
# datDNC_maT1 <- constructNegativeControlDat(Y = Y, X = X, C = C, lagNegCont = 1)
# head(datDNC_maT1)
# 
# # S2. Lag-structure
# datExamplelagT1 <- lagStructure(X = "o3", dat = datExample, lagT = 1)
# 
#
# **************************************************************************** #
# the pair z-test
# Ref. https://mgimond.github.io/Stats-in-R/z_t_tests.html#3_references
# **************************************************************************** #
pairZTest <- function(X1, X2, seX1, seX2) {
  z <- (X1-X2)/sqrt(seX1^2 + seX2^2)
  pval <- 2*pnorm(abs(z), lower.tail = FALSE)
  return(pval)
}

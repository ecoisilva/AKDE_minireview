
library(ctmm)

# Define and simulate from a movement model -------------------------------

# The position autocorrelation timescale (in sec)
tau_p <- 1 %#% 'day'

# The velocity autocorrelation timescale (in sec)
tau_v <- 1 %#% 'hour'

# The spatial variance in m^2
sig <- 100000

# True 95% range area
true.area <- -2*log(0.05)*pi*sig

# Generate the model
MODEL <- ctmm(tau = c(tau_p, tau_v),
              isotropic = TRUE,
              sigma = sig,
              mu = c(0,0))

# Define the sampling schedule (e.g., 10 days sampled hourly)
st <- seq(1,10 %#% 'days', # The duration
          1 %#% 'hours') # The sampling interval

sim <- simulate(MODEL, t = st)
plot(sim)


# Estimate a KDE HR from the data -----------------------------------------

START <- Sys.time() # Start time

# Fit an IID model
IID.mod <- ctmm.fit(sim, 
                    method = "ML")
summary(IID.mod)

# KDE range estimate
kdeRngEst <- akde(data = sim, CTMM = IID.mod)

END <- Sys.time() # End time

# KDE 95% range area
kdeRngArea <- sum(kdeRngEst$CDF <= 0.95) * prod(kdeRngEst$dr)

# KDE 95% range area error
KDE.error <- ( kdeRngArea - true.area ) / true.area

# The computation time
KDE.time <- END - START
KDE.time

# Reset the clock
rm(START); rm(END)


# Estimate an AKDE HR from the data ---------------------------------------

START <- Sys.time() # Start time

# Guesstimate the parameters of an autocorrelated model
GUESS <- ctmm.guess(data = sim, interactive = FALSE)

# Fit and select autocorrelated models
OUF.mod <- ctmm.select(data = sim, 
                       CTMM = GUESS,
                       method = "ML")
summary(OUF.mod)

# AKDE range estimate
akdeRngEst <- akde(data = sim, CTMM = OUF.mod, debias = FALSE)
summary(akdeRngEst)

END <- Sys.time() # End time

# AKDE 95% range area
akdeRngArea <- sum(akdeRngEst$CDF <= 0.95) * prod(akdeRngEst$dr)

# AKDE 95% range area error
AKDE.error <- ( akdeRngArea - true.area ) / true.area

# The computation time
AKDE.time <- END - START

#Reset the clock
rm(START); rm(END)


# Estimate an AKDEc HR from the data --------------------------------------

START <- Sys.time() # Start time

#Guesstimate the parameters of an autocorrelated model
GUESS <- ctmm.guess(data = sim,
                    interactive = FALSE)

#Fit and select autocorrelated models
OUF.mod <- ctmm.select(data = sim,
                       CTMM = GUESS,
                       method = "ML")

summary(OUF.mod)

#AKDEc range estimate
akdecRngEst <- akde(data = sim,
                    CTMM = OUF.mod,
                    debias = TRUE)

END <- Sys.time() # End time

# AKDEc 95% range area
akdecRngArea <- sum(akdecRngEst$CDF <= 0.95) * prod(akdecRngEst$dr)

# AKDEc 95% range area error
AKDEc.error <- ( akdecRngArea - true.area ) / true.area

# The computation time
AKDEc.time <- END - START


#Reset the clock
rm(START); rm(END)


# Estimate a weighted AKDEc HR from the data ------------------------------

START <- Sys.time() # Start time

# Guesstimate the parameters of an autocorrelated model
GUESS <- ctmm.guess(data = sim,
                    interactive = FALSE)

# Fit and select autocorrelated models
OUF.mod <- ctmm.select(data = sim,
                       CTMM = GUESS,
                       method = "ML")
summary(OUF.mod)

# Weighted AKDEc range estimate
wakdecRngEst <- akde(data = sim,
                     CTMM = OUF.mod,
                     debias = TRUE,
                     weights = TRUE)

END <- Sys.time() # End time

# Weighted AKDEc 95% range area
wakdecRngArea <- sum(wakdecRngEst$CDF <= 0.95) * prod(wakdecRngEst$dr)

# Weighted AKDEc 95% range area error
wAKDEc.error <- ( wakdecRngArea - true.area ) / true.area

# AKDE 95% error
AKDE.error <- ( akdeRngArea - true.area ) / true.area

# The computation time
wAKDEc.time <- END - START

#Reset the clock
rm(START); rm(END)


# Estimate a pHREML weighted AKDEc HR from the data -----------------------

START <- Sys.time() # Start time

# Guesstimate the parameters of an autocorrelated model
GUESS <- ctmm.guess(data = sim,
                    interactive = FALSE)

# Fit and select autocorrelated models
OUF.mod.phreml <- ctmm.select(data = sim,
                              CTMM = GUESS,
                              method = "pHREML")

summary(OUF.mod.phreml)

# pHREML weighted AKDEc range estimate
pHREML_wakdecRngEst <- akde(data = sim,
                            CTMM = OUF.mod.phreml,
                            debias = TRUE,
                            weights = TRUE)

END <- Sys.time() # End time

# pHREML weighted AKDEc 95% range area
pHREML_wakdecRngArea <- sum(pHREML_wakdecRngEst$CDF <= 0.95) * 
  prod(pHREML_wakdecRngEst$dr)

# pHREML weighted AKDEc 95% range area error
pHREML_wAKDEc.error <- ( pHREML_wakdecRngArea - true.area ) / true.area

# The computation time
pHREML_wAKDEc.time <- END - START

# Reset the clock
rm(START); rm(END)


# Estimate a bootstrapped pHREML weighted AKDEc HR from the data ----------

START <- Sys.time() # Start time

# Guesstimate the parameters of an autocorrelated model
GUESS <- ctmm.guess(data = sim,
                    interactive = FALSE)

# Fit and select autocorrelated models
OUF.mod.phreml <- ctmm.select(data = sim,
                              CTMM = GUESS,
                              method = "pHREML")

summary(OUF.mod.phreml)

# Bootstrap the fit
boot.fit <- ctmm.boot(data = sim,
                      CTMM = OUF.mod.phreml,
                      method = "pHREML",
                      iterate = TRUE)

summary(boot.fit)

# Bootstrapped weighted AKDEc range estimate
boot_pHREML_wakdecRngEst <- akde(data = sim,
                                 CTMM = boot.fit,
                                 debias = TRUE,
                                 weights = TRUE)

END <- Sys.time() # End time

# Bootstrapped pHREML weighted AKDEc 95% range area
boot_pHREML_wakdecRngArea <- sum(boot_pHREML_wakdecRngEst$CDF <= 0.95) *
  prod(boot_pHREML_wakdecRngEst$dr)

# Bootstrapped pHREML weighted AKDEc 95% range area error
boot_pHREML_wAKDEc.error <- ( boot_pHREML_wakdecRngArea - true.area ) / 
  true.area

# The computation time
boot_pHREML_wAKDEc.time <- END - START

# Reset the clock
rm(START); rm(END)


# Method comparisons ------------------------------------------------------

# Compare all of the relative area error side by side
ERROR <- c(KDE.error,
          AKDE.error,
          AKDEc.error,
          wAKDEc.error,
          pHREML_wAKDEc.error,
          boot_pHREML_wAKDEc.error)

ERROR

# Compare all of computation times side by side
COST <- c(KDE.time,
          AKDE.time,
          AKDEc.time,
          wAKDEc.time,
          pHREML_wAKDEc.time,
          boot_pHREML_wAKDEc.time)

COST

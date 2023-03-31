
library(ctmm)

# Define and simulate from a movement model -------------------------------

# The position autocorrelation timescale (in sec)
tau_p <- 1 %#% 'day'

# The velocity autocorrelation timescale (in sec)
tau_v <- 1 %#% 'hour'

# The spatial variance in m^2
sig <- 100000

# True 95% home range area
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
kde <- akde(data = sim, CTMM = IID.mod)

END <- Sys.time() # End time

# KDE 95% range area
KDE.area <- sum(kde$CDF <= 0.95) * prod(kde$dr)

# Metrics
KDE.bias <- KDE.area - true.area
KDE.error <- KDE.bias / true.area

# Computational cost
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
AKDE <- akde(data = sim, CTMM = OUF.mod, debias = FALSE)
summary(AKDE)

END <- Sys.time() # End time

# AKDE 95% range area
AKDE.area <- sum(AKDE$CDF <= 0.95) * prod(AKDE$dr)

# Metrics
AKDE.bias <- AKDE.area - true.area
AKDE.error <- KDE.bias / true.area

# Computational cost
AKDE.time <- END - START
AKDE.time

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
AKDEc <- akde(data = sim,
                    CTMM = OUF.mod,
                    debias = TRUE)

END <- Sys.time() # End time

# AKDEc 95% range area
AKDEc.area <- sum(AKDEc$CDF <= 0.95) * prod(AKDEc$dr)

# Metrics
AKDEc.bias <- AKDEc.area - true.area
AKDEc.error <- KDE.bias / true.area

# Computational cost
AKDEc.time <- END - START
AKDEc.time

# Reset the clock
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
wAKDEc <- akde(data = sim,
                     CTMM = OUF.mod,
                     debias = TRUE,
                     weights = TRUE)

END <- Sys.time() # End time

# Weighted AKDEc 95% range area
wAKDEc.area <- sum(wAKDEc$CDF <= 0.95) * prod(wAKDEc$dr)

# Metrics
wAKDEc.bias <- wAKDEc.area - true.area
wAKDEc.error <- KDE.bias / true.area

# Computational cost
wAKDEc.time <- END - START
wAKDEc.time

# Reset the clock
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
pHREML_wAKDEc <- akde(data = sim,
                            CTMM = OUF.mod.phreml,
                            debias = TRUE,
                            weights = TRUE)

END <- Sys.time() # End time

# pHREML weighted AKDEc 95% range area
pHREML_wAKDEc.area <- sum(pHREML_wAKDEc$CDF <= 0.95) * 
  prod(pHREML_wAKDEc$dr)

# Metrics
pHREML_wAKDEc.bias <- pHREML_wAKDEc.area - true.area
pHREML_wAKDEc.error <- KDE.bias / true.area

# Computational cost
pHREML_wAKDEc.time <- END - START
pHREML_wAKDEc.time

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
boot_pHREML_wAKDEc <- akde(data = sim,
                                 CTMM = boot.fit,
                                 debias = TRUE,
                                 weights = TRUE)

END <- Sys.time() # End time

# Bootstrapped pHREML weighted AKDEc 95% range area
boot_pHREML_wAKDEc.area <- sum(boot_pHREML_wAKDEc$CDF <= 0.95) *
  prod(boot_pHREML_wAKDEc$dr)

# Metrics
boot_pHREML_wAKDEc.bias <- boot_pHREML_wAKDEc.area - true.area
boot_pHREML_wAKDEc.error <- KDE.bias / true.area

# Computational cost
boot_pHREML_wAKDEc.time <- END - START
boot_pHREML_wAKDEc.time

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

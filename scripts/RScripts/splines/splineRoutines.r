#library(ggformula)
library(mgcv)
source("dateFormatRoutines.r")

# t, time points for which to infer the derivative
# gamToDerive, the spline/gam for which the derivatives are inferred
#gamDerivative <- function(tt.df, gamToDerive, outputpath) {
computeSplineDerivativeTable <- function(t, gamToDerive) {
  tt.df = data.frame(t=t)

  ## finite difference interval, delta t
  eps <- 1e-7

  #When type="lpmatrix" then a matrix is returned which yields the values of the linear predictor (minus any offset)
  #when postmultiplied by the parameter vector (in this case se.fit is ignored). The latter option is most useful for
  #getting variance estimates for quantities derived from the model: for example integrated quantities, or derivatives of smooths.
  X0 <- predict(gamToDerive, tt.df)
  #print(X0)
  X0_e <- exp(X0)
  #print(X0_e)
  # super mall interval for derivation
  tt_eps.df <- tt.df + eps

  X1 <- predict(gamToDerive, tt_eps.df)
  X1_e <- exp(X1)
  # finite difference approximation of first derivative
  # the design matrix
  Xp <- (X1_e - X0_e) / eps

  # first derivative
  d1_gam <- Xp#%*% coef(gamToDerive)

  d1_gam.table <- data.frame(t=t, value=exp(d1_gam))

  return(d1_gam.table)
}

# Using gratia package - compute confidence interval on the derivatives
# t, time points for which to infer the derivative
# gamToDerive, the spline/gam for which the derivatives are inferred
computeSplineDerivativeTableGratia <- function(t, gamToDerive, infDur) {
  tt.df = data.frame(t=t)
  #print(t)
  deriv <- fderiv(gamToDerive,tt.df)
  #print(deriv)
  #print(smooths(gamToDerive))
  # Normal CI
  ci <- confint(deriv, type = "confidence")
  #print(ci)
  # Simultaneous CI
  ci.sim <- confint(deriv, type = "simultaneous", nsim = 2500)
  # Size of interval
  ci.size <- ci.sim$est-ci.sim$lower
  #print(tt.df$t)
  #print(exp(as.vector(infDur*ci$lower)))
  #print(exp(as.vector(infDur*ci$upper)))
  d1_gam.table <- data.frame(t=tt.df$t, value=exp(infDur*ci$est), lower=exp(as.vector(infDur*ci$lower)), upper=exp(as.vector(infDur*ci$upper)),
                             si.lower=exp(as.vector(infDur*ci.sim$lower)),si.upper=exp(as.vector(infDur*ci.sim$upper)))
  return(d1_gam.table)
}


#generate random values from a multivariate normal distribution
random_mvn <- function(n, mu, sig) {
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}

# infer the confidence interval by taking the simultaneous interval for a function f(x) at a set of M locations
# simulation-based approach by Rupert et. al (2003)
calculateCI <- function(gam_mod_cs, xx) {
  #  Bayesian covariance matrix of the model coefficients Vb, with unknown (to be estimated) smoothing parameter
  Vb <- vcov(gam_mod_cs, unconditional = TRUE)
  
  # generate fitted values and CI in link function...
  splinePred_gam_cs <- predict(gam_mod_cs, data.frame(t=xx), se.fit = TRUE, type = "link")
  
  N <- 10000
  
  # N draws from approximately multivariate normaldistribution with mean vector 0 and covariance matrix Vb
  BUdiff <- random_mvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  
  # basis function at x
  Cx <- predict(gam_mod_cs, data.frame(t=xx), type = "lpmatrix")
  #deviations between the fitted and true parameters
  simDev <- Cx %*% t(BUdiff)
  #find the absolute values of the standardized deviations from the true model
  absDev <- abs(sweep(simDev, 1, splinePred_gam_cs$se.fit, FUN = "/"))
  #maximum of the absolute standardized deviations at the grid of x values for each simulation
  max_abs_sd <- apply(absDev, 2L, max)
  #calculate the critical value for a 95% simultaneous confidence interval/band
  crit <- quantile(max_abs_sd, prob = 0.95, type = 8)
  
  # ...and backtransform onto response scale
  return(data.frame(t=xx, value=gam_mod_cs$family$linkinv(splinePred_gam_cs$fit), se=splinePred_gam_cs$se.fit,
                    lower = gam_mod_cs$family$linkinv((splinePred_gam_cs$fit-2*splinePred_gam_cs$se.fit)),
                    upper = gam_mod_cs$family$linkinv((splinePred_gam_cs$fit+2*splinePred_gam_cs$se.fit)),
                    lowerSim = gam_mod_cs$family$linkinv((splinePred_gam_cs$fit-crit*splinePred_gam_cs$se.fit)),
                    upperSim = gam_mod_cs$family$linkinv((splinePred_gam_cs$fit+crit*splinePred_gam_cs$se.fit))
  ))
}

computeSpline <- function(input.table) {
  #input table must include a column with t, value and variance
  
  #number of data points
  N<-nrow(input.table)
  #weight: choose 1/variance, maybe better standard deviation?
  weights <- 1/(input.table$variance)
  #normalize splines by sum of weights
  weights = weights/sum(weights)
  
  # optimze k! (default k=10)
  # increase k=degree of freedom (number of base functions) 
  # until edf doesn't change substantially
  k=5
  edf_act=1
  repeat {
    edf = edf_act
    k<-k+5
    #compute splines, s=smooth function
    gam_mod_cs <- gam(value ~ s(t,bs="cs",k=k),
                      weights = weights,
                      data=input.table, method="REML", family=gaussian(link="log"))
    edf_act<- summary(gam_mod_cs)$edf
    if((edf_act-edf)<0.5 || k>= nrow(input.table))
      break
  }
  
  #return spline model
  return(gam_mod_cs)
}

computeSplineTable <- function(input.table) {
  gam_mod_cs <- computeSpline(input.table)
  # vector for which the outcomes should be predicted
  xx<- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)

  gam.table<-calculateCI(gam_mod_cs, xx)
  return(gam.table)
}

addSplineValuesForTrueN <- function(input.table, gam.table) {
  # optimize k (degree of freedom)
  k=5
  edf_act=1
  repeat {
    edf = edf_act
    k<-k+5
    #trueN comes from poisson distribution, hence no negative values (family=poisson -> link=logs)
    trueN_gam_cs <- gam(round(trueN) ~ s(t,bs="cs", k=k),
                        data=input.table, method="REML", family = poisson())
    edf_act<- summary(trueN_gam_cs)$edf
    if((edf_act-edf)<0.5 || k>= min(30,nrow(input.table)))
      break
  }
  
  xx<- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)
  
  gam.table_true<-calculateCI(trueN_gam_cs, xx)
  gam.table["value_trueN"] <-gam.table_true$value
  gam.table["value_trueN_lower"] <- gam.table_true$lowerSe
  gam.table["value_trueN_upper"] <- gam.table_true$upperSe
  gam.table["value_trueN_lowerSim"] <- gam.table_true$lowerSim
  gam.table["value_trueN_upperSim"] <- gam.table_true$upperSim
  
  return(gam.table)
}

computeSplineNewCasesTable <- function(input.table) {
  # optimize k (degree of freedom)
  k=5
  edf_act=1
  repeat {
    edf = edf_act
    k<-k+5
    #trueN comes from poisson distribution, hence no negative values (family=poisson -> link=logs)
    trueN_gam_cs <- gam(round(new_cases) ~ s(t,bs="cs", k=k),
                        data=input.table, method="REML", family = poisson())
    edf_act<- summary(trueN_gam_cs)$edf
    if((edf_act-edf)<0.5 || k>= min(30,nrow(input.table)))
      break
  }
  
  xx<- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)
  
  gam.table<-calculateCI(trueN_gam_cs, xx)
  
  return(gam.table)
}


computeInterpolation <- function (input.table, ts) {
  N_sampling <- 500
  p_sample <- 0.5
  width <- 10

  
  samplings <- matrix(ncol=length(ts), nrow=0)
  
  maxV =0
  for(i in seq(N_sampling)) {
    # sample for each row if it is in the sampling set
    sub_input.table <- input.table[runif(nrow(input.table))<p_sample, ]
    maxV = max(maxV, sub_input.table$value, na.rm = T)
    
    #### first interpolate, collect and smooth median afterwards
    interpol <- approx(sub_input.table$t, sub_input.table$value, xout=ts)
    samplings <-rbind(samplings, interpol$y)
    #samplings.df <- rbind(samplings.df, data.frame(t=interpol$x, value=interpol$y, simu=i))
  }
  
  # taking the quantiles of all samplings
  quantiles <- apply(samplings, 2, quantile, c(0.05, 0.5, 0.95), na.rm=T)
  # ... and smooth them
  smoothed_quantiles <- t(apply(quantiles, 1, filter, filter=rep(1,width)/width, sides=2))
  
  return(smoothed_quantiles)
}



# computeSplineNewCasesTable <- function(input.table) {
#   repCases_gam_cs <- gam(round(new_cases) ~ s(t,bs="cs"),
#                       data=input.table, method="REML")
#   xx<- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)
#   repCases_spline_gam_cs <- predict.gam(repCases_gam_cs, data.frame(t=xx), type = "response")
#   gam.table <- data.frame(t=xx, value=repCases_spline_gam_cs)
# 
#   return(gam.table)
# 
# }

computeRatio <- function(values, values_trueN) {
  return(values/values_trueN)
}

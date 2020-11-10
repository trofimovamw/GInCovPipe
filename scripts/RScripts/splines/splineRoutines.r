#library(ggformula)
library(mgcv)
source("dateFormatRoutines.r")

# t, time points for which to infer the derivative
# gamToDerive, the spline/gam for which the derivatives are inferred
#gamDerivative <- function(tt.df, gamToDerive, outputpath) {
computeSplineDerivativeTable <- function(t, gamToDerive, infDur) {
  tt.df = data.frame(t=t)

  ## finite difference interval, delta t
  eps <- 1e-7

  #When type="lpmatrix" then a matrix is returned which yields the values of the linear predictor (minus any offset)
  #when postmultiplied by the parameter vector (in this case se.fit is ignored). The latter option is most useful for
  #getting variance estimates for quantities derived from the model: for example integrated quantities, or derivatives of smooths.
  X0 <- predict(gamToDerive, tt.df, type = "lpmatrix")

  # super small interval for derivation
  tt_eps.df <- tt.df + eps
  X1 <- predict(gamToDerive, tt_eps.df, type = "lpmatrix")

  # finite difference approximation of first derivative
  # the design matrix
  Xp <- (X1 - X0) / eps
  # first derivative
  d1_gam <- Xp %*% coef(gamToDerive)
  # Standard error of prediction
  # Covariance matrix
  Vb <- vcov(gamToDerive)
  # SE
  d1_se <- sqrt(diag(Xp %*% Vb %*% t(Xp)))
  # Taking the exponential of the derivative times infectiousness period
  d1_gam.table <- data.frame(t=t, value=exp(infDur*d1_gam), lower=exp(infDur*d1_gam-2*d1_se), upper=exp(infDur*d1_gam+2*d1_se))

  return(d1_gam.table)
}

computeSpline <- function(input.table) {
    #input table must include a column with meanBinDate, value and variance
  # add day of year if not present yet
  #input.table <- addDayOfYearToTable(input.table)

  #number of data points
  N<-nrow(input.table)


  #weight: choose 1/variance, maybe better standard deviation?
  weights <- 1/input.table$variance
  #normalize splines by sum of weights
  weights = weights/sum(weights)

  #compute splines
  #s=smooth function
  gam_mod_cs <- gam(value ~ s(t,k = -1,bs="cs"),
                    weights = weights,
                    data=input.table, method="REML")

  #return spline model
  return(gam_mod_cs)
}

computeSplineTable <- function(input.table) {
  pseudoNumber=10^-10
  #log transform, to avoid negative values
  input.table$value <- log(input.table$value + pseudoNumber)
  gam_mod_cs <- computeSpline(input.table)
  # vector for which the outcomes should be predicted
  xx<- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)

  splinePred_gam_cs <- predict.gam(gam_mod_cs, data.frame(t=xx), se=T)

  splinePred_gam_cs$lower <- exp(splinePred_gam_cs$fit)-1.96*splinePred_gam_cs$se.fit
  splinePred_gam_cs$upper <- exp(splinePred_gam_cs$fit)+1.96*splinePred_gam_cs$se.fit
  #retransfrom
  gam.table <- data.frame(t=xx, value=exp(splinePred_gam_cs$fit), se=splinePred_gam_cs$se.fit, lower=splinePred_gam_cs$lower, upper=splinePred_gam_cs$upper)

  return(gam.table)
}

addSplineValuesForTrueN <- function(input.table, gam.table) {
  #trueN comes from poisson distribution, hence no negative values (family=poisson -> link=logs)
  trueN_gam_cs <- gam(round(trueN) ~ s(t,bs="cs"),
                      data=input.table, method="REML")#, family = poisson())
  xx <- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)
  # predict on same scale as response variables (type=response)
  trueN_spline_gam_cs <- predict.gam(object = trueN_gam_cs, newdata = data.frame(t=xx), type = "response")
  gam.table["value_trueN"] <- trueN_spline_gam_cs

  return(gam.table)
}

computeSplineNewCasesTable <- function(input.table) {
  repCases_gam_cs <- gam(round(new_cases) ~ s(t,bs="cs"),
                      data=input.table, method="REML")
  xx<- seq(min(input.table$t), max(input.table$t), len = max(input.table$t) - min(input.table$t)+1)
  repCases_spline_gam_cs <- predict.gam(repCases_gam_cs, data.frame(t=xx), type = "response")
  gam.table <- data.frame(t=xx, value=repCases_spline_gam_cs)

  return(gam.table)

}

computeRatio <- function(values, values_trueN) {
  return(values/values_trueN)
}

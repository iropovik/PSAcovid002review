# Load libraries
if (!require(lme4)) {install.packages('lme4')}
if (!require(psych)) {install.packages('psych')}
if (!require(reshape2)) {install.packages('reshape2')}

# Custom robust multivariate RE meta-analytic model
rmaCustom <- function(data = NA){
  if(nrow(data) > 20){
    rmaObject <- rma.mv(yi = yi, V = vi, data = data, method = "REML", random = ~ 1|study)
    rmaObject
  } else {
    rmaObject <- rma(yi = yi, vi = vi, data = data, method = "REML")
    rmaObject
  }
}

# 95% prediction interval
pi95 <- function(rmaObject = NA){
  pi95Out <- c("95% PI LB" = round(predict.rma(rmaObject)$cr.lb, 3), "95% PI UB" = round(predict.rma(rmaObject)$cr.ub, 3))
  pi95Out
}

# Proportion of significant results ---------------------------------------

propSig <- function(p.values = NA){
  as.integer(table(p.values < .05)[2])/length(p.values < .05)
}


# Heterogeneity -----------------------------------------------------------

heterogeneity <- function(rmaObject = NA){

  # Total heterogeneity - tau
  tau <- sqrt(sum(rmaObject$sigma2))


  # I^2
  W <- diag(1/rmaObject$vi)
  X <- model.matrix(rmaObject)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2<- 100 * sum(rmaObject$sigma2) / (sum(rmaObject$sigma2) + (rmaObject$k-rmaObject$p)/sum(diag(P)))

  c("Tau" = tau,
    "I^2" = I2)
}


# Permutation p-curve -----------------------------------------------------

# Evidential value; Permutation p-curve
# Subseting only the effects that are focal for the published study

# Choose effects from a single study by random (for the purpose of permutation p-curve)
duplicated.random = function(x, incomparables = FALSE, ...)
{
  if ( is.vector(x) )
  {
    permutation = sample(length(x))
    x.perm      = x[permutation]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
  else ( is.matrix(x) )
  {
    permutation = sample(nrow(x))
    x.perm      = x[permutation,]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
}

pcurve <- function(data = NA){
  #p-curve data export
  set.seed(123)
  # pcurveData <- na.omit(gsub("^.*?: ",": ", x = data$zPcurve[!duplicated.random(as.numeric(data$study))], replacement = ""))
  #write(pcurveData, "pcurve.data.txt")

  #Permutation p-curve
  cols <- 15
  pCurve <- matrix(ncol=cols, nrow=nsim)
  for(i in 1:nsim){
    pCurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = data$zPcurve[!duplicated.random(as.numeric(data$study))], replacement = "")))
  }
  pCurveOut <- data.frame(pCurve)
  pCurveOutCut <- as.data.frame(describe(pCurveOut[c(4, 6, 8, 10)]))[,3, drop = FALSE]
  rownames(pCurveOutCut) <- c("fullp", "fullp33", "halfp", "halfp33")
  pCurveOutCut
}


######################
# Code adapted from Carter, E. C., Schönbrodt, F. D., Hilgard, J., & Gervais, W. (2018). Correcting for bias in psychology: A comparison of meta-analytic methods. Retrieved from https://osf.io/rf3ys/.
# https://github.com/nicebread/meta-showdown/blob/master/MA-methods/7-Selection%20Models.R
# Return a result data frame either in wide or long format (for the 3PSM output)
library(dplyr)
returnRes <- function(res, long=TRUE, reduce=TRUE) {
  if (is.null(res)) return(NULL)

  # convert all factor columns to characters
  res %>% mutate_if(is.factor, as.character) -> res

  if (long==FALSE) {
    # return wide format
    return(res)
  } else {
    # transform to long format
    longRes <- melt(res, id.vars=c("method", "term"))
    if (reduce==TRUE & nrow(res) > 1) {longRes <- longRes %>% filter(!is.na(value)) %>% arrange(method, term, variable)}
    return(longRes)
  }
}

# 3-parameter selection model (3PSM)
# p-value intervals may be re-specified if they contain too few values
if (!require(weightr)) {
  install.packages('weightr')
}

threePSM.est <- function(d, v, min.pvalues=1, long=FALSE) {

  w1 <- tryCatch(
    weightfunct(d, v, steps = c(0.025, 1), mods = NULL, weights = NULL, fe = FALSE, table = TRUE),
    error = function(e) NULL
  )

  res.NA <- data.frame(
    method = "3PSM",
    term = c("tau2", "b0", "pr.nonsig"),
    estimate = NA,
    std.error = NA,
    statistic = NA,
    p.value = NA,
    conf.low = NA,
    conf.high = NA
  )

  if (is.null(w1)) return(returnRes(res.NA))

  # If <= 3 p-values in an interval: return NA
  p.table <- table(cut(w1$p, breaks=c(0, .025, 1)))
  if (any(p.table < min.pvalues)) {
    return(returnRes(res.NA))
  } else {
    est <- w1[[2]]$par

    # Compute standard errors from hessian
    std.err <- sqrt(abs(diag(solve(w1[[2]]$hessian))))


    res.wide <- data.frame(
      method = "3PSM",
      term = c("tau2", "b0", "pr.nonsig"),
      estimate = round(est, 4),
      std.error = round(std.err, 4),
      statistic = round(est/std.err, 4),
      p.value = round((pnorm(abs(est/std.err), lower.tail=FALSE)*2), 4),
      conf.low = round((est + qnorm(.025)*std.err), 4),
      conf.high = round((est + qnorm(1-.025)*std.err), 4)
    )
  }

  return(returnRes(res.wide))
}

fourPSM.est <- function(d, v, min.pvalues=0, long=TRUE, fallback = FALSE) {
  w1 <- tryCatch(
    weightfunct(d, v, steps = c(0.025, 0.5, 1), mods = NULL, weights = NULL, fe = FALSE, table = TRUE),
    error = function(e) NULL
  )

  res.NA <- data.frame(
    method = "4PSM",
    term = c("tau2", "b0", "pr.nonsig", "pr.opposite"),
    estimate = NA,
    std.error = NA,
    statistic = NA,
    p.value = NA,
    conf.low = NA,
    conf.high = NA
  )

  if (is.null(w1)) return(returnRes(res.NA))

  # if <= min.pvalues p-values in an interval: return NA
  p.table <- table(cut(w1$p, breaks=c(0, .025, 0.5, 1)))
  if (any(p.table < min.pvalues)) {
    if (fallback==TRUE) {
      return(threePSM.est(d, v, min.pvalues=min.pvalues, long=long))
    } else {
      return(returnRes(res.NA))
    }
  } else {
    est <- w1[[2]]$par

    # compute standard errors from hessian
    std.err <- sqrt(abs(diag(solve(w1[[2]]$hessian))))

    res.wide <- data.frame(
      method = "4PSM",
      term = c("tau2", "b0", "pr.nonsig", "pr.opposite"),
      estimate = round(est, 4),
      std.error = round(std.err, 4),
      statistic = round(est/std.err, 4),
      p.value = round(pnorm(est/std.err, lower.tail=FALSE)*2, 4),
      conf.low = round(est + qnorm(.025)*std.err, 4),
      conf.high = round(est + qnorm(1-.025)*std.err, 4)
    )
  }

  return(returnRes(res.wide))
}

#################
#PET-PEESE with 3PSM as the conditional estimator instead of PET

pet.peese <- function(d, v, study){
  if(length(d) < 30){
    pet <<- rma(yi = d ~ sqrt(v), vi = v, method="REML")
    pet.out <- round(c(pet$b[1], pet$se[1], pet$zval[1], pet$pval[1], pet$ci.lb[1], pet$ci.ub[1]), 3)
    names(pet.out) <- c("PET estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
    pet.out

    pet <- rma(yi = yi~ sqrt(vi), vi = vi,data = data, method="REML")
    pet$pval

    peese <<- rma(yi = d ~ v, vi = v, method="REML")
    peese.out <- round(c(peese$b[1], peese$se[1], peese$zval[1], peese$pval[1], peese$ci.lb[1], peese$ci.ub[1]), 3)
    names(peese.out) <- c("PEESE estimate", "se", "zval", "pval", "ci.lb", "ci.ub")

    ifelse(pet$pval < .05 & pet$b[1] > 0,
           return(peese.out),  return(pet.out))

  } else {
    pet <<- rma.mv(yi = d ~ sqrt(v), V = v, random = ~ 1|study, method="REML")
    pet.out <- round(c(pet$b[1], pet$se[1], pet$zval[1], pet$pval[1], pet$ci.lb[1], pet$ci.ub[1]), 3)
    names(pet.out) <- c("PET estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
    pet.out

    peese <<- rma.mv(yi = d ~ v, V = v, random = ~ 1|study, method="REML")
    peese.out <- round(c(peese$b[1], peese$se[1], peese$zval[1], peese$pval[1], peese$ci.lb[1], peese$ci.ub[1]), 3)
    names(peese.out) <- c("PEESE estimate", "se", "zval", "pval", "ci.lb", "ci.ub")

    fourPSM <- fourPSM.est(d, v)

    ifelse(fourPSM$value[4] < .05 & fourPSM$value[1] > 0,
           return(peese.out),  return(pet.out))
  }
}

# Publication bias --------------------------------------------------------

bias <- function(rmaObject = NA, alpha = .05, briefBias = TRUE){
  # Correlation between the ES and precision (SE)
  esPrec <- cor(rmaObject$yi, sqrt(rmaObject$vi), method = "kendall")
  # Small-study effects correction
  # 3-parameter selection model
  fourPSM <- fourPSM.est(rmaObject$yi, rmaObject$vi)

  # PET-PEESE
  pp <- pet.peese(rmaObject$yi, rmaObject$vi, rmaObject$mf.r[[1]]$study)

  if(briefBias == TRUE){
    return(list("4PSM ES estimate" = fourPSM[1, 4],
                "4PSM confidence interval" = c(fourPSM[5, 4], fourPSM[6, 4]),
                "4PSM p-value" = fourPSM[4, 4],
                "Whether PET or PEESE was used" = ifelse(fourPSM$value[4] < .05 & fourPSM$value[1] > 0, "PEESE", "PET"),
                "PET-PEESE ES estimate" = as.numeric(pp[1]),
                "PET-PEESE confidence interval" = as.numeric(c(pp[5], pp[6])),
                "PET-PEESE p-value" = pp[4]))}
  else{
    return(list("4PSM" = fourPSM,
                "PET-PEESE" = pp))
  }
}

# Summary results ---------------------------------------------------------

maResults <- function(rmaObject = NA, data = NA, alpha = .05, briefBias = F){
  list(
    "RMA results" = rmaObject,
    "Prediction interval" = pi95(rmaObject),
    "Heterogeneity" = heterogeneity(rmaObject),
    "p-curve" = pcurve(data),
    "Proportion of significant results" = propSig(data$p),
    "Publication bias" = bias(rmaObject, briefBias = briefBias))
}

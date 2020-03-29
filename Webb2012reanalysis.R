library(tidyr)
library(metafor)
library(compute.es)
library(psych)

rm(list = ls())
# No of simulations for the permutation p-curve
nsim <- 5 # Set to 5 just to make code checking/running fast.
# Assuming pre-post correlation of r = .5
corr <- .5
# For a more elaborate output from the pub bias tests/correction methods, set briefBias to FALSE
briefBias <- TRUE
# Subset only the designs with a control group?
RCTonly <- 1 # If no, set to 0

source("pcurve.R")
source("functions.R")

# read in the data
dataWide <- read.csv("Webb 2012 meta re-analysis.csv")

# define factors
dataWide$study <- as.factor(dataWide$study)
levels(dataWide$effect) <- c("reappraise emotional response", "reappraise emotional stimulus",
                         "reappraise via perspective taking", "reappraisal–mixed")

# convert from wide to long and order
data <- gather(dataWide, domain, yiD, experiential:physiological, na.rm = T)
data <- data[order(data$study),]

# compute ES, variances, and p-values for two-group designs
es2g <- with(data[!is.na(data$nCtrl),], des(yiD, nExp, nCtrl, id = study, verbose = F, dig = 1e10))
data$yi[!is.na(data$nCtrl)] <- es2g$g
data$vi[!is.na(data$nCtrl)] <- es2g$var.g
data$p[!is.na(data$nCtrl)] <- es2g$pval.g

# compute ES, variances, and p-values for one-group designs
corr <- .5 # Assuming pre-post correlation of r = .5
data$yi[is.na(data$nCtrl)] <- with(data[is.na(data$nCtrl),], (1 - (3/(4*nExp - 3))) * yiD)
data$vi[is.na(data$nCtrl)] <- with(data[is.na(data$nCtrl),], (1 - (3/(4*nExp - 3)))^2 * ((1 / nExp) + ((yiD^2) / (2 * nExp))) * 2 * (1 - corr))
data$p[is.na(data$nCtrl)] <- with(data[is.na(data$nCtrl),], 2*pt(abs(yi/sqrt(vi)), nExp - 1, lower.tail = FALSE))

# Create a results label based on p-values, converted to z-score.
data$zPcurve <- paste(data$study, ": ", "Z=", qnorm(1-(data$p)/2), sep = "")

# Subset only data for designs with a control group
if(RCTonly == 1){
  data <- data[!is.na(data$nCtrl),]
} else {
  paste("Result including designs lacking a control group")
}

# Subset data objects
dataResp <- data[data$effect == "reappraise emotional response",]
dataStim <- data[data$effect == "reappraise emotional stimulus",]
dataPers <- data[data$effect == "reappraise via perspective taking",]
dataMix <- data[data$effect == "reappraisal–mixed",]

# Meta-analysis -----------------------------------------------------------

#'## Meta-analysis
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
namesObjects <- levels(data$effect)[-4]
dataObjects <- list("reappraise emotional response" = dataResp, "reappraise emotional stimulus" = dataStim, "reappraise via perspective taking" = dataPers)
rmaObjects <- setNames(lapply(dataObjects, function(x){rmaCustom(x)}), nm = namesObjects)

# RMA for mixed reappraisal
rmaCustom(dataMix)

results <- list(NA)
for(i in 1:length(rmaObjects)){
  results[[i]] <- maResults(rmaObject = rmaObjects[[i]], data = dataObjects[[i]], alpha = .05, briefBias = briefBias)
}

results <- setNames(results, nm = namesObjects)
results

rmaObjects <- setNames(lapply(rmaObjects, function(x){metafor::funnel(x, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei")}), nm = namesObjects)

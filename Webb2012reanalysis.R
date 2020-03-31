#' ---
#' title: "Re-analysis of Webb et al. 2012 meta-analysis"
#' author: "Ivan Ropovik"
#' output:
#'    html_document:
#'       code_folding: hide
#'       toc: true
#'       toc_float: true
#'       fig_retina: 2
#'       theme: paper
#' always_allow_html: yes
#' ---
#+ setup, include=FALSE
knitr::opts_chunk$set(echo=FALSE)

# Libraries
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
# Subset only the experiential designs?
Experiential <- 1

# Source funcions
source("pcurve.R")
source("functions.R")

#'#### Re-analysis settings
paste("Subsetting only the designs with a control group set to ", RCTonly, ".", sep = "")
paste("Subsetting only experiential designs set to ", Experiential, ".", sep = "")

# Read in the data
dataWide <- read.csv("Webb 2012 meta re-analysis.csv")

# Define factors
dataWide$study <- as.factor(dataWide$study)
levels(dataWide$effect) <- c("reappraise emotional response", "reappraise emotional stimulus",
                         "reappraise via perspective taking", "reappraisal–mixed")

# Convert from wide to long and order
data <- gather(dataWide, domain, yiD, experiential:physiological, na.rm = T)
data <- data[order(data$study),]

# Compute ES, variances, and p-values for two-group designs
es2g <- with(data[!is.na(data$nCtrl),], des(yiD, nExp, nCtrl, id = study, verbose = F, dig = 1e10))
data$yi[!is.na(data$nCtrl)] <- es2g$g
data$vi[!is.na(data$nCtrl)] <- es2g$var.g
data$p[!is.na(data$nCtrl)] <- es2g$pval.g

# Compute ES, variances, and p-values for one-group designs
corr <- .5 # Assuming pre-post correlation of r = .5
data$yi[is.na(data$nCtrl)] <- with(data[is.na(data$nCtrl),], (1 - (3/(4*nExp - 3))) * yiD)
data$vi[is.na(data$nCtrl)] <- with(data[is.na(data$nCtrl),], (1 - (3/(4*nExp - 3)))^2 * ((1 / nExp) + ((yiD^2) / (2 * nExp))) * 2 * (1 - corr))
data$p[is.na(data$nCtrl)] <- with(data[is.na(data$nCtrl),], 2*pt(abs(yi/sqrt(vi)), nExp - 1, lower.tail = FALSE))

# Create a results label based on p-values, converted to z-score.
data$zPcurve <- paste(data$study, ": ", "Z=", qnorm(1-(data$p)/2), sep = "")

#'# RMA for mixed reappraisal
#'
#' Reproducing the RMA for mixed reappraisal
#+ include=TRUE
dataMix <- data[data$effect == "reappraisal–mixed",]
rmaCustom(dataMix)

# Subset only data for designs with a control group
if(RCTonly == 1){
  data <- data[!is.na(data$nCtrl),]
} else {
  paste("Result including designs lacking a control group")
}

# Subset only data for experiential domain
if(Experiential == 1){
  data <- data[data$domain == "experiential",]
} else {
  paste("Result including physiological and behavioral effects")
}

# Subset data objects
dataResp <- data[data$effect == "reappraise emotional response",]
dataStim <- data[data$effect == "reappraise emotional stimulus",]
dataPers <- data[data$effect == "reappraise via perspective taking",]

# Meta-analysis -----------------------------------------------------------
#'# Meta-analysis and bias correction
#'
#'For reappraise emotional response, reappraise emotional stimulus, and reappraise via perspective taking, respectively
namesObjects <- levels(data$effect)[-4]
dataObjects <- list("reappraise emotional response" = dataResp, "reappraise emotional stimulus" = dataStim, "reappraise via perspective taking" = dataPers)
rmaObjects <- setNames(lapply(dataObjects, function(x){rmaCustom(x)}), nm = namesObjects)

results <- list(NA)
for(i in 1:length(rmaObjects)){
  results[[i]] <- maResults(rmaObject = rmaObjects[[i]], data = dataObjects[[i]], alpha = .05, briefBias = briefBias)
}

results <- setNames(results, nm = namesObjects)
results

#'## Funnel plots
#'
#'Funnel plots for reappraise emotional response, reappraise emotional stimulus, and reappraise via perspective taking, respectively
rmaObjects <- lapply(setNames(rmaObjects, nm = namesObjects), function(x){metafor::funnel(x, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei")})

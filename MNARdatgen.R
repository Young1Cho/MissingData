####################################################
# Author: Young Won Cho
# Date: 2022-09-18
# Purpose: Generating Amputed Data using mice::ampute()
#          for simulation of univariate OSC model
####################################################

# Generate missing values with ampute()
# Reference: https://rianneschouten.github.io/mice_ampute/vignette/ampute.html

# packages ----
library("mice")

# How ampute() works ----
set.seed(2016)
testdata <- as.data.frame(MASS::mvrnorm(n = 10000,
                                        mu = c(10, 5, 0),
                                        Sigma = matrix(data = c(1.0, 0.2, 0.2,
                                                                0.2, 1.0, 0.2,
                                                                0.2, 0.2, 1.0),
                                                       nrow = 3,
                                                       byrow = TRUE)))
summary(testdata)
str(testdata)
var(testdata)

result <- ampute(data = testdata)
    result
    head(result$amp)

## Patterns: ----
    # "matrix" to define missing pattern; rows: patterns ; colums: variables.
    # 0: missingness; 1: complete variables
mypatterns <- result$patterns
mypatterns

mypatterns[2, 1] <- 0
mypatterns <- rbind(mypatterns, c(0, 1, 0))
mypatterns
result <- ampute(testdata, patterns = mypatterns)
md.pattern(result$amp)

## Freq ----
# The number of subsets is determined by the number of patterns.
# The size of the subsets is determined with a "frequency vector".
    # The # of elements in the vector must be equal to the number of patterns.
    # A vector value must be between 0 and 1, summing to 1 (100%).
result$freq
myfreq <- c(0.7, 0.1, 0.1, 0.1)
result <- ampute(testdata, freq = myfreq, patterns = mypatterns)
md.pattern(result$amp)

## Mech ----
# What kind of missingness mechanism we will implement?
# "string" with MCAR, MAR or MNAR.
    # MCAR: all data rows have the same probability of being amputed.
    # MAR : the info of missing is in the observed data.
    # MNAR: the info of missing is in the missing itself (dependent variable).

result$mech
# For MAR and MNAR, we calculate a weighted sum score for all data rows.
# It depends on "pre-specified weights"
    # Manipulated by a "k X m matrix" (k: # of patterns; m: # of variables).
    # Weights matrix: row = the patterns; colums = weights for variables.
    # The influence of the weights is relative within a pattern;
        # weight values of 0.4 and 0.2 = weight values of 4 and 2.

result$weights
# By default, the weights of 1: complete variable; 0 to variables with missing.
myweights <- result$weights
myweights[1, ] <- c(0, 2, 1)
myweights[3, ] <- c(3, 1, 0)
myweights
# Specifying the weights matrix
    # non-zero value: which variables affect the missingness?
    # zero value: which variables do not affect?

# mech: "MNAR"
result <- ampute(testdata, freq = myfreq,
                 patterns = mypatterns, mech = "MNAR")
result$patterns
result$weights
# > result$patterns
#   V1 V2 V3
# 1  0  1  1
# 2  0  0  1
# 3  1  1  0
# 4  0  1  0
# MCAR/MAR has non-zero weights for variables without missingness.
# Unlike MCAR/MAR, MNAR shows non-zero weights for variables with missingness.

## Linying code: ----
# The logic for the code is:
# missing data model
logity[i, t] <- phi0 + nmarphi * y[i, t] + marphi * xfull[i, t]
pry[i, t] <- exp(logity[i, t]) / (1 + exp(logity[i, t]))  # prob of missingness
ry[i, t] <- rbinom(1, 1, pry[i, t])                    # missing data indicator
#y is the DV that has missing value. x is fully observed.

# Set up the missing data model as needed:
# - if you want MCAR, just keep phi0 as non-zero.
# - if you want NMAR, keep both phi0 and nmarphi (and may be marphi)
# - if you want MAR, keep phi0 and marphi as non-zero.
# Also, we usually don't want initial observations to be missing in timeseries.

## Type ----
    # decide logistic prob. distribution (applied to the weighted sum scores)
    # This argument is used when cont == TRUE.
    # type: LEFT, MID, RIGHT or TAIL.
        # RIGHT (default): high weighted sum scores -> larger prob of missing
        # LEFT/MID/TAIL(both-tailed)
        # -> larger prob of missing for low/average/extreme weighted sum scores
# when cont == FALSE,
    # We can define the probability values manually.
    # using the "odds" argument.

tdf<- data.frame(Time=1:100, x=rnorm(100))

MyPattern <- matrix(data = c(1,0), nrow = 1, ncol = 2, byrow = FALSE)
MyOdds <- matrix(c(1, 2, 10, 0), nrow = 1, byrow = F)
result<-ampute(data = tdf,
   prop = 0.5,
   patterns = MyPattern,
   freq = NULL,
   mech = "MAR",
   weights = NULL,
   std = TRUE,
   cont = TRUE,
   type = "RIGHT",
   odds = MyOdds,
   bycases = TRUE,
   run = TRUE)

result
md.pattern(result$amp)
# compare the complete vs. amputed data
plot(tdf$x, type='b', col='grey')
lines(result$amp, col=2, type='o')

df<-result$amp
df$missing <- ifelse(is.na(df$x), 1, 0)

model <- glm(missing ~ Time, data = df, family = binomial(link = "logit")); summary(model)
model <- glm(missing ~ Time + I(Time^2), data = df, family = binomial(link = "logit")); summary(model)
plot(df$Time, model$fitted.values)
plot(df$Time, exp(model$fitted.values) / (1 + exp(model$fitted.values)))

## Creat Amputed Data ----

# 'long' data is from uniOSCdatagen.R file
long<-as.data.frame(long)
long$ID <- as.factor(long$ID)
str(long)

mypattern = matrix(c(1,1,0), nrow = 1)
myweight = matrix(c(0,0,1), nrow = 1)
result<-ampute(data = long,
   prop = 0.3,
   patterns = mypattern,
   freq = NULL,
   mech = "MNAR",
   weights = myweight,
   std = TRUE,
   cont = TRUE,
   type = "RIGHT",
   odds = NULL,
   bycases = TRUE,
   run = TRUE)
result

plot(long$Time[1:100], long$x1[1:100], type='b', col='grey')
lines(result$amp[1:100,2], result$amp[1:100,3], col=2, type='o')

## Run Mplus ----
# library(MplusAutomation)
# MplusDir<-getwd()
# dir.exists(MplusDir)
# MplusAutomation::runModels(target = getwd())

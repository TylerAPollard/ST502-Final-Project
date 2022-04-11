#### ST 502 Final Project
## Tyler Pollard and Apostolos Stamenos
## 30 March 2022

## Load required libraries
library(tidyverse)
library(ggplot2)
library(data.table)

## =========================================== PART 1 ===============================================
## Read in data
data <- fread("framingham_data.csv")

## Filter into smoker and nonsmoker data frames
# This dataset represents two independent samples of systolic blood pressure (sysBP) for 
# smokers (currentSmoker=1) and nonsmokers (currentSmoker=0) in the Framingham heart study
smoker_df <- data %>% filter(currentSmoker == 1)
nonsmoker_df <- data %>% filter(currentSmoker == 0)

## Set alpha level for all tests
alpha <- 0.05

## Pooled Two Sample t-test p-value method

## Pooled Two Sample t-test confidence interval method

## Satterthwaite Approximation Two Sample t-test p-value method
# Calcuate number of data points for both groups
n_nonsmoker <- nrow(nonsmoker_df)
n_smoker <- nrow(smoker_df)

# Calculate sample variances of both groups
samp_mean_nonsmoker <- mean(nonsmoker_df$sysBP)
samp_mean_smoker <- mean(smoker_df$sysBP)
samp_var_nonsmoker <- sum((nonsmoker_df$sysBP - samp_mean_nonsmoker)^2)/(n_nonsmoker - 1)
samp_var_smoker <- sum((smoker_df$sysBP - samp_mean_smoker)^2)/(n_smoker - 1)
var(nonsmoker_df$sysBP) # check nonsmoker sample variance
var(smoker_df$sysBP) # check nonsmoker sample variance

# Calculate degrees of freedom for t test
nu <- (samp_var_nonsmoker/n_nonsmoker + samp_var_smoker/n_smoker)^2/(
  (samp_var_nonsmoker/n_nonsmoker)^2/(n_nonsmoker - 1) + (samp_var_smoker/n_smoker)^2/(n_smoker - 1)
)
nu_num <- (samp_var_nonsmoker/n_nonsmoker + samp_var_smoker/n_smoker)^2
nu_den <- (samp_var_nonsmoker/n_nonsmoker)^2/(n_nonsmoker - 1) + (samp_var_smoker/n_smoker)^2/(n_smoker - 1)

# Calculate observed t test statistic
t_obs_satterthwaite <- (samp_mean_nonsmoker - samp_mean_smoker)/sqrt(samp_var_nonsmoker/n_nonsmoker + samp_var_smoker/n_smoker)

# Calculate p-value
p_value_satterthwaite <- 2*pt(t_obs_satterthwaite, df = nu, lower.tail = FALSE)

t_pdf <- function(t) {(gamma((nu+1)/2)/(sqrt(nu*pi)*gamma(nu/2)))*(1+t^2/nu)^(-(nu+1)/2)}
t_pdf2 <- function(t) {(factorial((nu+1)/2-1)/(sqrt(nu*pi)*factorial(nu/2-1)))*(1+t^2/nu)^(-(nu+1)/2)}
p_value_satterthwaite2 <- 2*integrate(t_pdf2, lower = t_obs_satterthwaite, upper = Inf)$value

t.test(nonsmoker_df$sysBP, smoker_df$sysBP, var.equal = FALSE, conf.level = .95)
t.test(nonsmoker_df$sysBP, smoker_df$sysBP, var.equal = TRUE, conf.level = .95)

## Satterthwaite Approxiamtion Two Sample t-test confidence interval method
se_satterthwaite <- sqrt(samp_var_nonsmoker/(n_nonsmoker) + samp_var_smoker/(n_smoker))
ci_satterthwaite <- c((samp_mean_nonsmoker - samp_mean_smoker) - qt(1-(alpha/2), df = nu)*se_satterthwaite,
                      (samp_mean_nonsmoker - samp_mean_smoker) + qt(1-(alpha/2), df = nu)*se_satterthwaite)
ci_satterthwaite

## Checking normal assumption
# Nonsmoker data
hist(nonsmoker_df$sysBP, breaks = seq(min(nonsmoker_df$sysBP), max(nonsmoker_df$sysBP), length.out = 20))
qqnorm(nonsmoker_df$sysBP)
qqline(nonsmoker_df$sysBP)

# Smoker data
hist(smoker_df$sysBP, breaks = seq(min(smoker_df$sysBP), max(smoker_df$sysBP), length.out = 20))
qqnorm(smoker_df$sysBP)
qqline(smoker_df$sysBP)

## =========================================== PART 2 ===============================================
## Varying parameters
# True variances
true_var1 <- c(1, 4, 9)
true_var2 <- 1

# Sample Sizes
n1 <- c(10, 30, 70)
n2 <- c(10, 30, 70)

# True mean difference mu1 - mu2
true_mean_diff <- c(-5, -1, 1, 5)
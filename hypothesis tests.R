#### ST 502 Final Project
## Tyler Pollard and Apostolos Stamenos
## 26 March 2022

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

## Satterthwaite Approxiamtion Two Sample t-test confidence interval method


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
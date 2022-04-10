#### ST 502 Final Project
## Tyler Pollard and Apostolos Stamenos
## 30 March 2022

## Load required libraries
library(tidyverse)
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
n1 <- length(nonsmoker_df$sysBP)
n2 <- length(smoker_df$sysBP)
df <- n1 + n2 - 2

ybar_1 <- sum(nonsmoker_df$sysBP)/n1
ybar_2 <- sum(smoker_df$sysBP)/n2
diff <- ybar_1 - ybar_2
D_0 <- 0

s_sq1 <- sum((nonsmoker_df$sysBP-ybar_1)^2)/(n1-1)
s_sq2 <- sum((smoker_df$sysBP-ybar_2)^2)/(n2-1)
s_sq_p <- ((n1-1)*s_sq1+(n2-1)*s_sq2)/df

se <- (sqrt(s_sq_p)*sqrt(1/n1+1/n2))
T <- (diff - D_0)/se
p_val <- 2*pt(T, df)

## Pooled Two Sample t-test confidence interval method
t_quants <- qt(c(alpha/2, 1-alpha/2), df)
CI <- diff+se*t_quants

## Satterthwaite Approximation Two Sample t-test p-value method

## Satterthwaite Approximation Two Sample t-test confidence interval method


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


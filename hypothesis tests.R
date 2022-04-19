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
data <- data %>% 
  mutate(currentSmoker = if_else(currentSmoker == 1, 'Smoker', 'Nonsmoker'))
smoker_df <- data %>% filter(currentSmoker == "Smoker")
nonsmoker_df <- data %>% filter(currentSmoker == "Nonsmoker")

## Set alpha level for all tests
alpha <- 0.05

## Shared statistics
# Calculate number of data points in each group
n_nonsmoker <- length(nonsmoker_df$sysBP)
n_smoker <- length(smoker_df$sysBP)

# Calculate sample means
samp_mean_nonsmoker <- sum(nonsmoker_df$sysBP)/n_nonsmoker
samp_mean_smoker <- sum(smoker_df$sysBP)/n_smoker

# Calculate sample variances
samp_var_nonsmoker <- sum((nonsmoker_df$sysBP - samp_mean_nonsmoker)^2)/(n_nonsmoker - 1)
samp_var_smoker <- sum((smoker_df$sysBP - samp_mean_smoker)^2)/(n_smoker - 1)

## Pooled Two Sample t-test p-value method
df <- n_nonsmoker + n_smoker - 2

diff <- samp_mean_nonsmoker - samp_mean_smoker
D_0 <- 0

samp_var_pooled <- ((n_nonsmoker-1)*samp_var_nonsmoker+(n_smoker-1)*samp_var_smoker)/df

se_pooled <- (sqrt(samp_var_pooled)*sqrt(1/n_nonsmoker+1/n_smoker))
T_pooled <- (diff - D_0)/se_pooled
p_val_pooled <- 2*pt(T_pooled, df, lower.tail = FALSE)

## Pooled Two Sample t-test confidence interval method
t_quants <- qt(c(alpha/2, 1-alpha/2), df)
CI_pooled <- diff+se_pooled*t_quants

## Satterthwaite Approximation Two Sample t-test p-value method
# Calculate degrees of freedom for t test
nu <- (samp_var_nonsmoker/n_nonsmoker + samp_var_smoker/n_smoker)^2/(
  (samp_var_nonsmoker/n_nonsmoker)^2/(n_nonsmoker - 1) + (samp_var_smoker/n_smoker)^2/(n_smoker - 1)
)
nu <- floor(nu)

# Calculate observed t test statistic
t_obs_satterthwaite <- (samp_mean_nonsmoker - samp_mean_smoker)/sqrt(samp_var_nonsmoker/n_nonsmoker + samp_var_smoker/n_smoker)

# Calculate p-value
p_value_satterthwaite <- 2*pt(abs(t_obs_satterthwaite), df = nu, lower.tail = FALSE)

t.test(nonsmoker_df$sysBP, smoker_df$sysBP, var.equal = FALSE, conf.level = .95)
t.test(nonsmoker_df$sysBP, smoker_df$sysBP, var.equal = TRUE, conf.level = .95)

## Satterthwaite Approximation Two Sample t-test confidence interval method
se_satterthwaite <- sqrt(samp_var_nonsmoker/(n_nonsmoker) + samp_var_smoker/(n_smoker))
ci_satterthwaite <- c((samp_mean_nonsmoker - samp_mean_smoker) - qt(1-(alpha/2), df = nu)*se_satterthwaite,
                      (samp_mean_nonsmoker - samp_mean_smoker) + qt(1-(alpha/2), df = nu)*se_satterthwaite)
ci_satterthwaite

## Checking normal assumption
# Nonsmoker data
# Histogram
hist(nonsmoker_df$sysBP)

# Q-Q Plot
ggplot(data = nonsmoker_df, aes(sample = sysBP)) + 
  geom_qq() +
  geom_qq_line()

# ECDF vs CDF
nonsmoker_probs <- pnorm(sort(nonsmoker_df$sysBP), mean = ybar_1, sd = sqrt(s_var1))
nonsmoker_norm_df <- data.frame(x = sort(nonsmoker_df$sysBP), y = nonsmoker_probs)
ggplot(data = nonsmoker_df, aes(sysBP)) +
  stat_ecdf() +
  geom_line(data = nonsmoker_norm_df, aes(x = x, y = y)) 
  
# Smoker data
# Histogram
hist(smoker_df$sysBP)

# Q-Q Plot
ggplot(data = smoker_df, aes(sample = sysBP)) + 
  geom_qq() +
  geom_qq_line()

# ECDF vs CDF
smoker_probs <- pnorm(sort(smoker_df$sysBP), mean = ybar_2, sd = sqrt(s_var2))
smoker_norm_df <- data.frame(x = sort(smoker_df$sysBP), y = smoker_probs)
ggplot(data = smoker_df, aes(sysBP)) +
  stat_ecdf() +
  geom_line(data = smoker_norm_df, aes(x = x, y = y)) 

# Boxplot
ggplot(data = data, aes(x = currentSmoker, y = sysBP)) + 
  geom_boxplot()

summary(nonsmoker_df$sysBP)
summary(smoker_df$sysBP)

## 4
f_stat <- s_var1/s_var2
pf(f_stat, df1 = 225-1, df2 = 75-1, lower.tail = FALSE)


## =========================================== PART 2 ===============================================
## Varying parameters
# True variances
true_var1 <- c(1, 4, 9)
true_var2 <- 1

# Sample Sizes
n1 <- c(10, 30, 70)
n2 <- c(10, 30, 70)

# True mean difference mu1 - mu2
true_mean_diff <- c(-5, -1, 0, 1, 5)

## Function to generate datasets for simulation study
generate_data <- function(replicates = 100, n1, n2, var1, var2, mean_difference, mu1 = 137) {
  nonsmokers <- matrix(NA, nrow = n1, ncol = replicates)
  smokers <- matrix(NA, nrow = n2, ncol = replicates)
  
  # Generate multiple datasets
  for(i in 1:replicates) {
    nonsmokers[, i] <- rnorm(n1, mu1, sqrt(var1))
    smokers[, i] <- rnorm(n2, mu1 - mean_difference, sqrt(var2))
  }
  return(list(nonsmokers = nonsmokers, smokers = smokers))
}

## Generate all possible data sets
generated_data_sets <- list()
for(i in true_mean_diff){
  for(j in n1){
    for(k in n2){
      for(l in true_var1){
        data_name <- paste0( "true_mean_diff=", i, ".", "n1=", j, ".", "n2=", k, ".", "var1=", l, ".", "var2=1")
        generated_data_sets[[data_name]] <- generate_data(replicates = 100, n1 = j, n2 = k, var1 = l, var2 = 1, mean_difference = i, mu1 = 137)
      }
    }
  }
}

## Satterthwaite Approximation Two Sample t-test p-value method function
# Create function for calculating Satterthwaite p-value
satterthwaite_p_value <- function(data1, data2){
  # Determine number of data points in each data set
  n1 <- length(data1)
  n2 <- length(data2)
  
  # Calculate sample means
  samp_mean1 <- sum(data1)/n1
  samp_mean2 <- sum(data2)/n2
  
  # Calculate sample variances
  samp_var1 <- sum((data1 - samp_mean1)^2)/(n1 - 1)
  samp_var2 <- sum((data2 - samp_mean2)^2)/(n2 - 1)
  
  # Calculate degrees of freedom for t test
  nu <- (samp_var1/n1 + samp_var2/n2)^2/(
    (samp_var1/n1)^2/(n1 - 1) + (samp_var2/n2)^2/(n2 - 1)
  )
  nu <- floor(nu)
  
  # Calculate observed t test statistic under null hypothesis mu1 - mu2 = 0
  null_t_obs_satterthwaite <- (samp_mean1 - samp_mean2)/sqrt(samp_var1/n1 + samp_var2/n2)
  
  # Calculate p-value under null hypothesis mu1 - mu2 = 0
  null_p_value_satterthwaite <- 2*pt(abs(null_t_obs_satterthwaite), df = nu, lower.tail = FALSE)
  
  # Calculate observed t test statistic under alternative hypothesis mu1 - mu2 = -5
  alt_t_obs_satterthwaite_neg5 <- ((samp_mean1 - samp_mean2) - (-5))/sqrt(samp_var1/n1 + samp_var2/n2)
  alt_t_obs_satterthwaite_neg1 <- ((samp_mean1 - samp_mean2) - (-1))/sqrt(samp_var1/n1 + samp_var2/n2)
  alt_t_obs_satterthwaite_1 <- ((samp_mean1 - samp_mean2) - (1))/sqrt(samp_var1/n1 + samp_var2/n2)
  alt_t_obs_satterthwaite_5 <- ((samp_mean1 - samp_mean2) - (5))/sqrt(samp_var1/n1 + samp_var2/n2)
  
  # Calculate p-value under null hypothesis mu1 - mu2 = -5
  alt_p_value_satterthwaite_neg5 <- 2*pt(abs(alt_t_obs_satterthwaite_neg5), df = nu, lower.tail = FALSE)
  alt_p_value_satterthwaite_neg1 <- 2*pt(abs(alt_t_obs_satterthwaite_neg1), df = nu, lower.tail = FALSE)
  alt_p_value_satterthwaite_1 <- 2*pt(abs(alt_t_obs_satterthwaite_1), df = nu, lower.tail = FALSE)
  alt_p_value_satterthwaite_5 <- 2*pt(abs(alt_t_obs_satterthwaite_5), df = nu, lower.tail = FALSE)
  # alt_p_value_satterthwaite2 <- 1-((1-pt(abs(alt_t_obs_satterthwaite), df = nu, lower.tail = FALSE)) -
  #   pt(abs(alt_t_obs_satterthwaite), df = nu, lower.tail = FALSE))
  
  data.frame(alpha = null_p_value_satterthwaite,
             power_neg5 = alt_p_value_satterthwaite_neg5,
             power_neg1 = alt_p_value_satterthwaite_neg1,
             power_1 = alt_p_value_satterthwaite_1,
             power_5 = alt_p_value_satterthwaite_5)
}

## Calculate Satterthwaite alpha values
satterthwaite_alpha_values <- vector("numeric", length = length(generated_data_sets))
satterthwaite_power_neg5_values <- vector("numeric", length = length(generated_data_sets))
satterthwaite_power_neg1_values <- vector("numeric", length = length(generated_data_sets))
satterthwaite_power_1_values <- vector("numeric", length = length(generated_data_sets))
satterthwaite_power_5_values <- vector("numeric", length = length(generated_data_sets))
for(m in 1:length(generated_data_sets)){
  simulated_alpha_values_satterthwaite <- vector("numeric", length = 100)
  simulated_power_neg5_values_satterthwaite <- vector("numeric", length = 100)
  simulated_power_neg1_values_satterthwaite <- vector("numeric", length = 100)
  simulated_power_1_values_satterthwaite <- vector("numeric", length = 100)
  simulated_power_5_values_satterthwaite <- vector("numeric", length = 100)
  for(n in 1:100){
    simulated_alpha_values_satterthwaite[n] <- satterthwaite_p_value(generated_data_sets[[m]]$nonsmokers[,n], generated_data_sets[[m]]$smokers[,n])$alpha
    simulated_power_neg5_values_satterthwaite[n] <- satterthwaite_p_value(generated_data_sets[[m]]$nonsmokers[,n], generated_data_sets[[m]]$smokers[,n])$power_neg5
    simulated_power_neg1_values_satterthwaite[n] <- satterthwaite_p_value(generated_data_sets[[m]]$nonsmokers[,n], generated_data_sets[[m]]$smokers[,n])$power_neg1
    simulated_power_1_values_satterthwaite[n] <- satterthwaite_p_value(generated_data_sets[[m]]$nonsmokers[,n], generated_data_sets[[m]]$smokers[,n])$power_1
    simulated_power_5_values_satterthwaite[n] <- satterthwaite_p_value(generated_data_sets[[m]]$nonsmokers[,n], generated_data_sets[[m]]$smokers[,n])$power_5
  }
  simulated_alpha <- sum(simulated_alpha_values_satterthwaite < 0.05)/length(simulated_alpha_values_satterthwaite)
  simulated_power_neg5 <- sum(simulated_power_neg5_values_satterthwaite < 0.05)/length(simulated_power_neg5_values_satterthwaite)
  simulated_power_neg1 <- sum(simulated_power_neg1_values_satterthwaite < 0.05)/length(simulated_power_neg1_values_satterthwaite )
  simulated_power_1 <- sum(simulated_power_1_values_satterthwaite < 0.05)/length(simulated_power_1_values_satterthwaite)
  simulated_power_5 <- sum(simulated_power_5_values_satterthwaite < 0.05)/length(simulated_power_5_values_satterthwaite)
  satterthwaite_alpha_values[m] <- simulated_alpha
  satterthwaite_power_neg5_values[m] <- simulated_power_neg5
  satterthwaite_power_neg1_values[m] <- simulated_power_neg1
  satterthwaite_power_1_values[m] <- simulated_power_1
  satterthwaite_power_5_values[m] <- simulated_power_5
}

simulated_alpha_df <- data.frame(test = names(generated_data_sets), 
                                 alpha = satterthwaite_alpha_values,
                                 power_neg5 = satterthwaite_power_neg5_values,
                                 power_neg1 = satterthwaite_power_neg1_values,
                                 power_1 = satterthwaite_power_1_values,
                                 power_5 = satterthwaite_power_5_values
                                 )

rm(simulated_alpha_df)
##### Testing satterthwaite function
t1 <- generated_data_sets[[55]]$nonsmokers[,1]
t2 <- generated_data_sets[[55]]$smokers[,1]

n1 <- length(t1)
n2 <- length(t2)

# Calculate sample means
samp_mean1 <- sum(t1)/n1
samp_mean2 <- sum(t2)/n2

# Calculate sample variances
samp_var1 <- sum((t1 - samp_mean1)^2)/(n1 - 1)
samp_var2 <- sum((t2 - samp_mean2)^2)/(n2 - 1)

# Calculate degrees of freedom for t test
nu <- (samp_var1/n1 + samp_var2/n2)^2/(
  (samp_var1/n1)^2/(n1 - 1) + (samp_var2/n2)^2/(n2 - 1)
)
#nu <- floor(nu)

# Calculate observed t test statistic under null hypothesis mu1 - mu2 = 0
null_t_obs_satterthwaite <- (samp_mean1 - samp_mean2)/sqrt(samp_var1/n1 + samp_var2/n2)

# Calculate p-value under null hypothesis mu1 - mu2 = 0
null_p_value_satterthwaite <- 2*pt(abs(null_t_obs_satterthwaite), df = nu, lower.tail = FALSE)

# Calculate observed t test statistic under alternative hypothesis mu1 - mu2 /= 0
alt_t_obs_satterthwaite <- ((samp_mean1 - samp_mean2) - (-5))/sqrt(samp_var1/n1 + samp_var2/n2)

# Calculate p-value under null hypothesis mu1 - mu2 = 0
alt_p_value_satterthwaite <- 2*pt(abs(alt_t_obs_satterthwaite), df = nu, lower.tail = FALSE)
alt_p_value_satterthwaite2 <- 1-((1-pt(abs(alt_t_obs_satterthwaite), df = nu, lower.tail = FALSE)) -
                                   pt(abs(alt_t_obs_satterthwaite), df = nu, lower.tail = FALSE))
t.test(t1, t2, var.equal = FALSE)
##### end testing function

## Satterthwaite Approximation Two Sample t-test p-value method
## Satterthwaite Approximation Two Sample t-test confidence interval method
#### ST 502 Final Project
## Tyler Pollard and Apostolos Stamenos
## 30 March 2022

## Load required libraries
library(tidyverse)
library(data.table)
library(MESS)

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
p_val_pooled <- 2*pt(abs(T_pooled), df = df, lower.tail = FALSE)

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
col_nonsmoker <- 'springgreen4'
col_smoker <- 'darkred'

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

# Combined visualizations
ggplot(data = data) + 
  geom_histogram(aes(x = sysBP, fill = currentSmoker), alpha = 0.5, bins = 20) +
  scale_fill_manual(name = '', values = c('Smoker' = col_smoker, 'Nonsmoker' = col_nonsmoker))

ggplot(data = data, aes(sample = sysBP, color = currentSmoker)) + 
  geom_qq(alpha = 0.5) +
  geom_qq_line() +
  scale_colour_manual(name = '', values = c('Smoker' = col_smoker, 'Nonsmoker' = col_nonsmoker))


# ECDF vs CDF
smoker_probs <- pnorm(sort(smoker_df$sysBP), mean = ybar_2, sd = sqrt(s_var2))
smoker_norm_df <- data.frame(x = sort(smoker_df$sysBP), y = smoker_probs)
ggplot(data = smoker_df, aes(sysBP)) +
  stat_ecdf() +
  geom_line(data = smoker_norm_df, aes(x = x, y = y)) 

# Boxplot
ggplot(data = data, aes(x = currentSmoker, y = sysBP, fill = currentSmoker)) + 
  scale_fill_manual(name = '', values = c('Smoker' = col_smoker, 'Nonsmoker' = col_nonsmoker)) +
  geom_boxplot(alpha = 0.7)

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

## Pooled and Satterthwaite Approximation Two Sample t-test p-value method function
# Create function for calculating p-value
p_value <- function(data1, data2, equal.variance = FALSE){
  # Determine number of data points in each data set
  n1 <- length(data1)
  n2 <- length(data2)
  
  # Calculate sample means
  samp_mean1 <- sum(data1)/n1
  samp_mean2 <- sum(data2)/n2
  
  # Calculate sample variances
  samp_var1 <- sum((data1 - samp_mean1)^2)/(n1 - 1)
  samp_var2 <- sum((data2 - samp_mean2)^2)/(n2 - 1)
  
  # Calculate degrees of freedom and observed t statistic for Pooled t test
  if(equal.variance == TRUE){
    ## Pooled Two Sample t-test p-value method
    df <- n1 + n2 - 2
    
    diff <- samp_mean1 - samp_mean2
    D_0 <- 0
    
    samp_var_pooled <- ((n1-1)*samp_var1+(n2-1)*samp_var2)/df
    
    se_pooled <- (sqrt(samp_var_pooled)*sqrt(1/n1+1/n2))
    T_pooled <- (diff - D_0)/se_pooled
    p_value <- 2*pt(abs(T_pooled), df = df, lower.tail = FALSE)
  }
  # Calculate degrees of freedom and observed t statistic for Satterthwaite t test
  else{
    # Calculate degrees of freedom for t test
    nu <- (samp_var1/n1 + samp_var2/n2)^2/(
      (samp_var1/n1)^2/(n1 - 1) + (samp_var2/n2)^2/(n2 - 1)
    )
    nu <- floor(nu)
    
    # Calculate observed t test statistic
    t_obs_satterthwaite <- (samp_mean1 - samp_mean2)/sqrt(samp_var1/n1 + samp_var2/n2)
    
    # Calculate p-value 
    p_value <- 2*pt(abs(t_obs_satterthwaite), df = nu, lower.tail = FALSE)
  }
  
  return(p_value)
}

## Calculate Simulated alpha and power values
simulated_values <- function(data, equal.variance = FALSE){
  # Create empty vectors to store simulated alpha and power values for each unique combination of parameters 
  # in generated data
  alpha_values <- c(rep(NA, length(data)))
  power_values <- c(rep(NA, length(data)))
  # Create for loop to run through each of the data sets in generated data
  for(m in 1:length(data)){
    # Create empty vectors to store simulated alpha and power values for each of the 100 generated data sets 
    # within a unique combination of parameter values in generated data
    simulated_alpha_values <- c(rep(NA, 100))
    simulated_power_values <- c(rep(NA, 100))
    # Create for loop to run through all 100 generated data sets within a unique combination of parameter 
    # values in generated data
    for(n in 1:100){
      # If statement to determine if generated data is from null hypothesis. If it is calculate p-value and 
      # store in simulated_alpha_values vector
      if(grepl("true_mean_diff=0", names(data)[m])){
        # Calculate simulated alpha values for Pooled T Test
        if(equal.variance == TRUE){
          simulated_alpha_values[n] <- p_value(data[[m]]$nonsmokers[,n], data[[m]]$smokers[,n], equal.variance = TRUE)
        }
        # Calculate simulated alpha values for Satterthwaite T Test
        else{
          simulated_alpha_values[n] <- p_value(data[[m]]$nonsmokers[,n], data[[m]]$smokers[,n], equal.variance = FALSE)
        }
        # If statement to determine if generated data is from alternative hypothesis. If it is calculate p-value and 
        # store in simulated_power_values vector
      }else{
        # Calculate simulated power values for Pooled T Test
        if(equal.variance == TRUE){
          simulated_power_values[n] <- p_value(data[[m]]$nonsmokers[,n], data[[m]]$smokers[,n], equal.variance = TRUE)
        }
        # Calculate simulated power values for Satterthwaite T Test
        else{
          simulated_power_values[n] <- p_value(data[[m]]$nonsmokers[,n], data[[m]]$smokers[,n], equal.variance = FALSE)
        }
      }
    }
    # Calculate proportion of p-values that would result in rejecting null hypothesis
    simulated_alpha <- sum(simulated_alpha_values < 0.05)/length(simulated_alpha_values)
    simulated_power <- sum(simulated_power_values < 0.05)/length(simulated_alpha_values)
    # Store simulated proportion into alpha and power vectors. If the true mean difference is 0 then the simulated 
    # proportion will correspond to alpha and stored, otherwise it is NA. If the true mean difference is not 0 then the
    # simulated proportion will correspond to power and stored, otherwise it is NA.
    alpha_values[m] <- simulated_alpha
    power_values[m] <- simulated_power
  }
  # Combine results into data frame to output
  simulated_df <- data.frame(test = names(data), 
                             alpha = alpha_values,
                             sim_power = power_values
  )
  # Output simulated alpha and power values data frame
  return(simulated_df)
}

## Calculate Exact alpha and power values
exact_power_values <- function(n1, n2, true_var1, true_var2, true_mean_diff, equal.variance = FALSE){
  # Create empty vector for generated data set name with parameter values identified and empty vector for exact power
  test = c() 
  calc_power = c()  
  # For loop heirarchy to run through all possible combinations of parameters used to generate data
  for(i in true_mean_diff){
    for(j in n1){
      for(k in n2){
        for(l in true_var1){
          for(p in true_var2){
            test <- c(test, paste0( "true_mean_diff=", i, ".", "n1=", j, ".", "n2=", k, ".", "var1=", l, ".", "var2=", p))
            # If statement to determine lower sample size of the two groups
            # n1 >= n2
            if(j >= k){
              # Calculate power for Pooled T Test
              if(equal.variance == TRUE){
                calc_power <- c(calc_power, power_t_test(
                  n = k, 
                  ratio = j/k, 
                  delta = i, 
                  sd = sqrt(((j-1)*l+(k-1)*p)/(j+k-2)), 
                  #sd.ratio = sqrt(l)/sqrt(p),
                  sig.level = 0.05,
                  alternative = "two.sided")$power)
              }
              # Calculate power for Satterthwaite T Test
              else{
                calc_power <- c(calc_power, power_t_test(
                  n = k, 
                  ratio = j/k, 
                  delta = i, 
                  sd = sqrt(p), 
                  sd.ratio = sqrt(l)/sqrt(p),
                  sig.level = 0.05,
                  alternative = "two.sided")$power)
              }
            }
            # n1 < n2
            else{
              # Calculate power for Pooled T Test
              if(equal.variance == TRUE){
                calc_power <- c(calc_power, power_t_test(
                  n = j, 
                  ratio = k/j, 
                  delta = i, 
                  sd = sqrt(((j-1)*l+(k-1)*p)/(j+k-2)), 
                  #sd.ratio = sqrt(p)/sqrt(l),
                  sig.level = 0.05,
                  alternative = "two.sided")$power)
              }
              # Calculate power for Satterthwaite T Test
              else{
                calc_power <- c(calc_power, power_t_test(
                  n = j, 
                  ratio = k/j, 
                  delta = i, 
                  sd = sqrt(l), 
                  sd.ratio = sqrt(p)/sqrt(l),
                  sig.level = 0.05,
                  alternative = "two.sided")$power)
              }
            }
          }
        }
      }
    }
  }
  exact_power_df <- data.frame(test = test, 
                               calc_power = calc_power
  )
  return(exact_power_df)
}

# Calculate simulated Satterthwaite alpha and power values data frame from generated data 
satterthwaite_simulated_power_df <- simulated_values(generated_data_sets, equal.variance = FALSE)

# Calculate exact power for Satterthwaite T test
satterthwaite_exact_power_df <- exact_power_values(n1 = n1, n2 = n2, true_var1 = true_var1, true_var2 = true_var2, true_mean_diff = true_mean_diff, equal.variance = FALSE)

# Create power table for Satterthwaite T Test
satterthwaite_power_df <- full_join(satterthwaite_simulated_power_df, satterthwaite_exact_power_df)

# Calculate simulated Pooled alpha and power values data frame from generated data 
pooled_simulated_power_df <- simulated_values(generated_data_sets, equal.variance = TRUE)

# Calculate exact power for Pooled T test
pooled_exact_power_df <- exact_power_values(n1 = n1, n2 = n2, true_var1 = true_var1, true_var2 = true_var2, true_mean_diff = true_mean_diff, equal.variance = TRUE)

# Create power table for Pooled T Test
pooled_power_df <- full_join(pooled_simulated_power_df, pooled_exact_power_df)



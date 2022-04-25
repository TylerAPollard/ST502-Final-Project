#### ST 502 Final Project
## Tyler Pollard and Apostolos Stamenos
## 30 March 2022

## Load required libraries
library(tidyverse)
library(data.table)
library(MESS)
library(car)

## =========================================== PART 1 ===============================================
## Read in data
# This dataset represents two independent samples of systolic blood pressure (sysBP) for 
# smokers (currentSmoker=1) and nonsmokers (currentSmoker=0) in the Framingham heart study
data <- fread("framingham_data.csv")

## Filter into smoker and nonsmoker data frames
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

## -------------------- Pooled Two Sample t-test p-value method -----------------------------
# Calculate the degrees of freedom
df <- n_nonsmoker + n_smoker - 2

# Calculate the point estimate for difference of sample means
diff <- samp_mean_nonsmoker - samp_mean_smoker
# Set the value of true difference of population means to 0 under the null hypothesis
D_0 <- 0

# Calculate pooled sample variance
samp_var_pooled <- ((n_nonsmoker-1)*samp_var_nonsmoker+(n_smoker-1)*samp_var_smoker)/df

# Calculate the standard error 
se_pooled <- (sqrt(samp_var_pooled)*sqrt(1/n_nonsmoker+1/n_smoker))
# Calculate the observed t statistic
T_pooled <- (diff - D_0)/se_pooled
# Calculate the p-value using observed t statistic and degrees of freedom
p_val_pooled <- 2*pt(abs(T_pooled), df = df, lower.tail = FALSE)

## -------------------- Pooled Two Sample t-test confidence interval method -----------------
# Calculate the 0.025 and 0.975 t quantiles for pooled degrees of freedom
t_quants <- qt(c(alpha/2, 1-alpha/2), df)
# Calculate confidence interval using difference of sample means, t quantiles, and standard error
CI_pooled <- diff+se_pooled*t_quants

## -------------- Satterthwaite Approximation Two Sample t-test p-value method --------------
# Calculate degrees of freedom for t test using Satterthwaite approximation
nu <- (samp_var_nonsmoker/n_nonsmoker + samp_var_smoker/n_smoker)^2/(
  (samp_var_nonsmoker/n_nonsmoker)^2/(n_nonsmoker - 1) + (samp_var_smoker/n_smoker)^2/(n_smoker - 1)
)
# Round down the degrees of freedom to nearest integer
nu <- floor(nu)

# Calculate observed t test statistic
t_obs_satterthwaite <- (samp_mean_nonsmoker - samp_mean_smoker)/sqrt(samp_var_nonsmoker/n_nonsmoker + samp_var_smoker/n_smoker)

# Calculate the p-value using observed t statistic and degrees of freedom
p_value_satterthwaite <- 2*pt(abs(t_obs_satterthwaite), df = nu, lower.tail = FALSE)

## ---------- Satterthwaite Approximation Two Sample t-test confidence interval method ------
# Calculate the standard error
se_satterthwaite <- sqrt(samp_var_nonsmoker/(n_nonsmoker) + samp_var_smoker/(n_smoker))
# Calculate confidence interval using difference of sample means, 0.025 and 0.975 t quantiles, and standard error
ci_satterthwaite <- c((samp_mean_nonsmoker - samp_mean_smoker) - qt(1-(alpha/2), df = nu)*se_satterthwaite,
                      (samp_mean_nonsmoker - samp_mean_smoker) + qt(1-(alpha/2), df = nu)*se_satterthwaite)

## ------------------------- Checking normal assumption ------------------------------------
## Combined visualizations
# Set colors for each data set
col_nonsmoker <- 'springgreen4'
col_smoker <- 'darkred'

# Histogram
ggplot(data = data) + 
  geom_histogram(aes(x = sysBP, fill = currentSmoker), alpha = 0.5, bins = 20) +
  scale_fill_manual(name = '', values = c('Smoker' = col_smoker, 'Nonsmoker' = col_nonsmoker)) +
  labs(x = "Systolic Blood Pressure", y = "Number of Participants") +
  theme_bw()

# QQ Plots
ggplot(data = data, aes(sample = sysBP, color = currentSmoker)) + 
  geom_qq(alpha = 0.5) +
  geom_qq_line() +
  facet_grid(cols = vars(currentSmoker)) +
  scale_colour_manual(name = '', values = c('Smoker' = col_smoker, 'Nonsmoker' = col_nonsmoker)) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_bw() +
  theme(legend.position = "none") 

# Boxplots
ggplot(data = data, aes(x = currentSmoker, y = sysBP, fill = currentSmoker)) + 
  scale_fill_manual(name = '', values = c('Smoker' = col_smoker, 'Nonsmoker' = col_nonsmoker)) +
  labs(x = "Current Smoking Status", y = "Systolic Blood Pressure") +
  geom_boxplot(alpha = 0.7) +
  theme_bw()

## Check equal variance assumption using Levene Test for medians
leveneTest(sysBP ~ as.factor(currentSmoker), data = data, center = "median")


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

# Setting the seed
set.seed(4242022)

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
    # Calculate degrees of freedom
    df <- n1 + n2 - 2
    
    # Calculate sample mean difference
    diff <- samp_mean1 - samp_mean2
    # Assign hypothesized difference equal to 0 to represent Null hypothesis
    D_0 <- 0
    
    # Calculate pooled sample variance
    samp_var_pooled <- ((n1-1)*samp_var1+(n2-1)*samp_var2)/df
    
    # Calculate pooled standard error
    se_pooled <- (sqrt(samp_var_pooled)*sqrt(1/n1+1/n2))
    
    # Calculate pooled t statistic
    T_pooled <- (diff - D_0)/se_pooled
    
    # Calculate pooled probability
    p_value <- 2*pt(abs(T_pooled), df = df, lower.tail = FALSE)
  }
  # Calculate degrees of freedom and observed t statistic for Satterthwaite t test
  else{
    # Calculate degrees of freedom for t test
    nu <- (samp_var1/n1 + samp_var2/n2)^2/(
      (samp_var1/n1)^2/(n1 - 1) + (samp_var2/n2)^2/(n2 - 1)
    )
    # Round down nu value to whole number
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
  prob_values <- c(rep(NA, length(data)))
  metric <- c(rep(NA, length(data)))
  # Create for loop to run through each of the data sets in generated data
  for(m in 1:length(data)){
    # Create empty vectors to store simulated alpha and power values for each of the 100 generated data sets 
    # within a unique combination of parameter values in generated data
    simulated_prob_values <- c(rep(NA, 100))
    # Create for loop to run through all 100 generated data sets within a unique combination of parameter 
    # values in generated data
    for(n in 1:100){
      # If statement to determine if generated data is from null hypothesis. If it is calculate p-value and 
      # store in simulated_alpha_values vector
      if(grepl("true_mean_diff=0", names(data)[m])){
        # Calculate simulated alpha values for Pooled T Test
        if(equal.variance == TRUE){
          simulated_prob_values[n] <- p_value(data[[m]]$nonsmokers[,n], data[[m]]$smokers[,n], equal.variance = TRUE)
          simulated_metric <- "alpha"
        }
        # Calculate simulated alpha values for Satterthwaite T Test
        else{
          simulated_prob_values[n] <- p_value(data[[m]]$nonsmokers[,n], data[[m]]$smokers[,n], equal.variance = FALSE)
          simulated_metric <- "alpha"
        }
        # If statement to determine if generated data is from alternative hypothesis. If it is calculate p-value and 
        # store in simulated_power_values vector
      }else{
        # Calculate simulated power values for Pooled T Test
        if(equal.variance == TRUE){
          simulated_prob_values[n] <- p_value(data[[m]]$nonsmokers[,n], data[[m]]$smokers[,n], equal.variance = TRUE)
          simulated_metric <- "power"
        }
        # Calculate simulated power values for Satterthwaite T Test
        else{
          simulated_prob_values[n] <- p_value(data[[m]]$nonsmokers[,n], data[[m]]$smokers[,n], equal.variance = FALSE)
          simulated_metric <- "power"
        }
      }
    }
    # Calculate proportion of p-values that would result in rejecting null hypothesis
    simulated_prob <- sum(simulated_prob_values < 0.05)/length(simulated_prob_values)
    # Store simulated proportion into alpha and power vectors. If the true mean difference is 0 then the simulated 
    # proportion will correspond to alpha and stored, otherwise it is NA. If the true mean difference is not 0 then the
    # simulated proportion will correspond to power and stored, otherwise it is NA.
    prob_values[m] <- simulated_prob
    metric[m] <- simulated_metric
  }
  # Combine results into data frame to output
  simulated_df <- data.frame(test = names(data), 
                             probs = prob_values,
                             metric = metric
  )
  # Output simulated alpha and power values data frame
  return(simulated_df)
}

# Calculate simulated Satterthwaite alpha and power values data frame from generated data 
satterthwaite_simulated_power_df <- simulated_values(generated_data_sets, equal.variance = FALSE)
# Add column for test type identifier
satterthwaite_simulated_power_df$type <- "Satterthwaite"

# Calculate simulated Pooled alpha and power values data frame from generated data 
pooled_simulated_power_df <- simulated_values(generated_data_sets, equal.variance = TRUE)
# Add column for test type identifier
pooled_simulated_power_df$type <- "Pooled"

## Data Visuals for Part 2
# Combine the Pooled and Satterthwaite power data frames into one 
combined_df <- rbind(pooled_simulated_power_df, satterthwaite_simulated_power_df)
# Create column identifiers for the varying parameters that will be used to parse data in the visuals
combined_df2 <- combined_df %>% 
  mutate(
    delta = case_when(
      grepl("true_mean_diff=-5", combined_df$test) ~ -5,
      grepl("true_mean_diff=-1", combined_df$test) ~ -1,
      grepl("true_mean_diff=0", combined_df$test) ~ 0,
      grepl("true_mean_diff=1", combined_df$test) ~ 1,
      grepl("true_mean_diff=5", combined_df$test) ~ 5
    ),
    n1 = case_when(
      grepl("n1=10", combined_df$test) ~ 10,
      grepl("n1=30", combined_df$test) ~ 30,
      grepl("n1=70", combined_df$test) ~ 70,
    ),
    n2 = case_when(
      grepl("n2=10", combined_df$test) ~ 10,
      grepl("n2=30", combined_df$test) ~ 30,
      grepl("n2=70", combined_df$test) ~ 70,
    ),
    var1 = case_when(
      grepl("var1=1", combined_df$test) ~ 1,
      grepl("var1=4", combined_df$test) ~ 4,
      grepl("var1=9", combined_df$test) ~ 9,
    ),
    var2 = 1
  )

# Filter combined data frame to create separate power and alpha data frames 
combined_power_df <- combined_df2 %>% 
  filter(metric != "alpha") %>%
  mutate(type = as.factor(type),
         var1 = as.factor(var1))
combined_alpha_df <- combined_df2 %>% filter(metric == "alpha")

# Create labels names to facet the alpha and power plots by. This will allow us to produce a 3x3 grid
# of plots faceted by the varying sample sizes
n1.labs <- c('n1 = 10', 'n1 = 30', 'n1 = 70')
n2.labs <- c('n2 = 10', 'n2 = 30', 'n2 = 70')
names(n1.labs) <- c(10, 30, 70)
names(n2.labs) <- c(10, 30, 70)

# For alpha create 3x3 figure of varying n's and on x axis plot var1 and y axis alpha and parse by test 
# Create using points with line 
ggplot(data = combined_alpha_df) +
  geom_point(aes(x = var1, y = probs, color = type), size = 1.5) +
  geom_line(aes(x = var1, y = probs, color = type), size = 1) +
  facet_grid(n2 ~ n1, labeller = labeller(
    n1 = n1.labs,
    n2 = n2.labs
  )) +
  scale_x_continuous(breaks = c(1,4,9)) +
  scale_colour_manual(name = "Test", values = c('Pooled' = 'darkorange1', 'Satterthwaite' = 'mediumpurple1')) +
  labs(x = expression(sigma[1]^2), y = expression(paste("Alpha, ", alpha))) +
  theme_bw()

# For power create 3x3 figure of varying n's and on x axis plot true mean diff and y axis power and parse by test and var
# Create 6 line plots for each of the 9 graphs. Different line types for each test and different color for each
# variance. 
ggplot(data = combined_power_df) + 
  geom_line(aes(x = delta, y = probs, color = var1, linetype = type), size = 1) + 
  facet_grid(n2 ~ n1, labeller = labeller(
    n1 = n1.labs,
    n2 = n2.labs
  )) +
  scale_x_continuous(breaks = c(-5,-1,1,5)) +
  scale_linetype_discrete(name = "Test") +
  scale_colour_manual(name = expression(sigma[1]^2), values = c('1' = 'springgreen4', '4' = 'darkred', '9' = 'steelblue')) +
  labs(x = expression(paste("True Mean Difference, ", Delta)), y = expression(paste('Power = ', 1-beta))) +
  theme_bw()




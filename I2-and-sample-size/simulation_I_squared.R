library(metafor)
library(dplyr)
library(data.table)
library(ggplot2)
library(faux)
library(scales)

set.seed(123)

# Setting parameters
N <- 200000
mu_d <- 1
tau <- 0.1
tau2 <- tau^2  

# Function to calculate I² 
I2 <- function(v, tau2) {
  wi <- 1 / v
  m <- length(v)  
  v_tilde <- ((m-1) * sum(wi)) / (sum(wi)^2 - sum(wi^2))  # Calculate estimated sampling variance
  return(100 * (tau2 / (tau2 + v_tilde)))  
}

# Initialize data frame for storing I² values
I2_data <- data.frame(n = integer(), m = integer(), I2 = numeric())

# Calculate I² values for varying sample sizes
for (n in seq(5, 10000, by = 5)) {
  m <- floor(N / n)  
  d_values <- rnorm(m, mean = mu_d, sd = tau)  
  var_d_values <- 1 / n + d_values^2 / (2 * n)  # Calculate variance for each study
  I2_value <- I2(var_d_values, tau2)  
  I2_data <- rbind(I2_data, data.frame(n = n, m = m, I2 = I2_value)) 
}

# Set up up the simulation.We will use Treatment vs. Control as manipulated within subject, thus the will be correlated.
# Correlation of 0.5 is set up to obtain a Cohen's d of 1 when the standard deviation for both condition is of 1.

sample_size <- c(25, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 5000, 10000) 
effect_sd <- 1
correlation <- 0.5  

all_sim_results <- data.table()
study_id_counter <- 1

for (n in sample_size) {
  m <- floor(N / n)
  
  for (study_id in 1:m) {
    effect <- rnorm(1, mu_d, tau)
    data <- rnorm_multi(
      n = n,
      mu = c(0, effect),
      sd = c(effect_sd, effect_sd),
      r = correlation,
      varnames = c("Control", "Treatment")
    )
    
    data <- as.data.table(data)
    data[, `:=` (Study_ID = study_id_counter, Participant_ID = seq_len(n), Sample_Size = n)]
    all_sim_results <- rbindlist(list(all_sim_results, data), use.names = TRUE, fill = TRUE)
    study_id_counter <- study_id_counter + 1
  }
}

# Function to compute effect sizes
calculate_summary <- function(x, y) {
  n <- length(x)
  mean_diff <- mean(x - y)
  sd_diff <- sd(x - y)
  d <- mean_diff / sd_diff
  var_d <- (1/n + d^2/(2*n)) * 2 * (1 - cor(x, y))
  sd_d <- sqrt(var_d)
  return(c(d = d, var_d = var_d, sd_d = sd_d, sd_diff = sd_diff))
}

summary_data <- all_sim_results %>%
  group_by(Study_ID, Sample_Size) %>%
  summarise(
    yi = calculate_summary(Treatment, Control)[1],
    vi = calculate_summary(Treatment, Control)[2],
    .groups = 'drop'
  )

models_fitting <- data.table()

for (n in sample_size) {
  current_data <- subset(summary_data, Sample_Size == n)
  
  model <- rma(yi, vi , method="REML", data = current_data)
  I2 <- model$I2
  tau <- sqrt(model$tau2)
  d <- coef(model)[1]
  
  models_fitting <- rbindlist(list(models_fitting, data.table(
    Sample_Size = n,
    I2 = I2,
    d = d,
    tau = tau
  )), use.names = TRUE, fill = TRUE)
}




p <- ggplot() +
  geom_line(data = I2_data, aes(x = n, y = I2), color = "#1f77b4", linewidth = 1, alpha = 0.8) +
  geom_point(data = models_fitting, aes(x = Sample_Size, y = I2), color = "#ff7f0e", size = 3, shape = 19, alpha = 0.9) +
  scale_x_log10(labels = label_number(), breaks = trans_breaks("log10", function(x) 10^x)) +
  labs(title = "I² values across different sample sizes",
       x = "sample size (n)",
       y = "I² (%)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  scale_colour_manual("Legend", 
                      labels = c("Theory", "Observed"),
                      values = c("#1f77b4", "#ff7f0e"))
print(p)

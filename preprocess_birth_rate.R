setwd('/home/rstudio/workspace/interest_rate_birth_rate')

# Data preprocess
# 조출생률 data (annual data)
library(KFAS)
temp <- read.csv('출생아수__합계출산율__자연증가_등_20250328222127.csv')
temp <- t(temp)[-1, 3]
annual_rates <- as.numeric(temp)

# Create monthly data by uniformly distributing annual values across 12 months
n_years <- length(annual_rates)
monthly_rates <- rep(annual_rates, each = 12) / 12

# Create a date sequence for these monthly values.
# (Assume a starting year, e.g., 2000. Adjust as needed based on your data.)
start_year <- 1970
dates <- seq(as.Date(paste0(start_year, "-01-01")), 
             by = "month", 
             length.out = length(monthly_rates))

# Combine the dates and rates into a data frame
monthly_data <- data.frame(date = dates, monthly_rate = monthly_rates)

# Plot the uniformly distributed monthly 조출생률
library(ggplot2)
ggplot(monthly_data, aes(x = date, y = monthly_rate)) +
  geom_line(color = "blue") +
  labs(title = "Uniformly Distributed Monthly 조출생률",
       x = "Date", y = "Monthly Rate") +
  theme_minimal()

# -------------------------------
# Gaussian Filtering Implementation
# -------------------------------

# Define the Gaussian kernel parameters
sigma <- 3
kernel_half <- 6  # Half-length of kernel (you can adjust this value)
x <- -kernel_half:kernel_half
gaussian_kernel <- dnorm(x, mean = 0, sd = sigma)
# Normalize the kernel so that its sum is 1
gaussian_kernel <- gaussian_kernel / sum(gaussian_kernel)

# Apply the Gaussian filter using convolution.
# Using 'sides = 2' gives a centered moving average.
monthly_data$gaussian_smoothed <- stats::filter(monthly_data$monthly_rate,
                                                filter = gaussian_kernel,
                                                sides = 2)

# Save the monthly_data data frame to a CSV file without row names
write.csv(monthly_data, file = "조출생률_smoothed_monthly_data.csv", row.names = FALSE)

# Plot the original uniform values and the Gaussian smoothed values
ggplot(monthly_data, aes(x = date)) +
  geom_line(aes(y = monthly_rate), color = "blue", linetype = "dashed", 
            size = 1, alpha = 0.7) +
  geom_line(aes(y = gaussian_smoothed), color = "green", size = 1) +
  labs(title = "Monthly 조출생률: Original (Dashed Blue) vs Gaussian Smoothed (Green)",
       x = "Date", y = "Monthly Rate") +
  theme_minimal()

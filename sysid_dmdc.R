###############################################
#
# Rage: 2003.11 - 2024-06
#
###############################################



######################
#
# Preprocessing data
#
######################


# 조출생률
birth_rate <- read.csv('조출생률_smoothed_monthly_data.csv')[407:654,]
row.names(birth_rate) <-NULL


# 기준금리
library(dplyr)
library(tidyr)
library(lubridate)
interest_rate <- read.csv('기준금리.csv', stringsAsFactors = FALSE)

# Convert the date column to Date objects
interest_rate$변경일자 <- as.Date(interest_rate$변경일자)

# Create a complete sequence of months from the minimum to maximum date
interest_rate <- interest_rate %>%
  complete(변경일자 = seq(from = min(변경일자), to = max(변경일자), by = "month")) %>%
  arrange(변경일자) %>%
  # Fill missing interest rate values with the last observation carried forward
  fill(기준금리, .direction = "down")

interest_rate <- interest_rate %>%
  # Create a new column that represents the month (e.g., 2025-02-01 for February 2025)
  mutate(month = floor_date(변경일자, "month")) %>%
  # Group by the month
  group_by(month) %>%
  # Select the row with the latest modification date in that month
  slice_max(변경일자, with_ties = FALSE) %>%
  ungroup()

interest_rate <- interest_rate[55:302,]


# 실업률
unemploy_rate <- data.frame(t(read.csv('성_연령별_실업률_20250328224250.csv'))[-(1:2),])
colnames(unemploy_rate) <- c('실업률')
head(unemploy_rate)
# Extract the row names from your data frame
date_str <- rownames(unemploy_rate)

# Remove the leading 'X' and add the day (".01")
date_str_clean <- paste0(sub("^X", "", date_str), ".01")

# Convert the string to Date objects; format is year.month.day
dates <- as.Date(date_str_clean, format = "%Y.%m.%d")

# Optionally, add these dates as a new column in your data frame
unemploy_rate$date <- dates

rownames(unemploy_rate) <- NULL

unemploy_rate <- unemploy_rate[54:301,]


# 소비자물가지수
price <- data.frame(t(read.csv('소비자물가지수_2020100__20250328223530.csv'))[-(1:2),])
colnames(price) <- c('소비자물가지수')
# Extract the row names from your data frame
date_str <- rownames(price)

# Remove the leading 'X' and add the day (".01")
date_str_clean <- paste0(sub("^X", "", date_str), ".01")

# Convert the string to Date objects; format is year.month.day
dates <- as.Date(date_str_clean, format = "%Y.%m.%d")

# Optionally, add these dates as a new column in your data frame
price$date <- dates

rownames(price) <- NULL

price <- price[466:713,]


# 주택매매가격
house <- data.frame(t(read.csv('유형별_매매가격지수_20250328223335.csv'))[-(1:2),])
colnames(house) <- c('매매가격')
# Extract the row names from your data frame
date_str <- rownames(house)

# Remove the leading 'X' and add the day (".01")
date_str_clean <- paste0(sub("^X", "", date_str), ".01")

# Convert the string to Date objects; format is year.month.day
dates <- as.Date(date_str_clean, format = "%Y.%m.%d")

# Optionally, add these dates as a new column in your data frame
house$date <- dates

rownames(house) <- NULL

house <- house[1:248,]


######################
#
# Construct data matrix X
#
######################
X_whole <- t(data.matrix(data.frame(birth=birth_rate$gaussian_smoothed,
                unemploy=unemploy_rate$실업률,
                price=price$소비자물가지수,
                house=house$매매가격)))
U <- t(data.matrix(data.frame(interest=interest_rate$기준금리)))[,1:247]

X <- X_whole[,1:247]
Xnext <- X_whole[,2:248]

O <- rbind(X,U)

Odecomposed <- svd(O)
A <- Xnext %*% Odecomposed$v %*% diag(Odecomposed$d^(-1)) %*% t(Odecomposed$u[1:4,])
B <- Xnext %*% Odecomposed$v %*% diag(Odecomposed$d^(-1)) %*% Odecomposed$u[5,]



######################
#
# Simulate
#
######################

sys1 = function(x,u) {
  return(A%*%x+B%*%u)
}

T <- ncol(X)

# Initialize the state matrix (4 variables x T time steps)
x_sim <- matrix(0, nrow = 4, ncol = T)

# Set the initial condition to be the first column of X (ground truth)
x_sim[, 1] <- X[, 1]

# Simulate the system forward using sys1
for (t in 1:(T - 1)) {
  x_sim[, t + 1] <- sys1(x_sim[, t], U[t])
}

# Create a date sequence assuming the simulation starts in November 2003.
time_index <- seq.Date(from = as.Date("2003-11-01"), length.out = T, by = "month")

# Construct a data frame comparing the ground truth and simulated birth rates.
# Note: we assume the birth rate is the first row of the state vector.
sim_birth_rate <- data.frame(
  date = time_index,
  simulated = x_sim[1, ],
  ground_truth = X[1, ]
)

# Load ggplot2 and plot the results
library(ggplot2)
ggplot(sim_birth_rate, aes(x = date)) +
  geom_line(aes(y = ground_truth, color = "Ground Truth")) +
  geom_line(aes(y = simulated, color = "Simulated")) +
  labs(
    title = "Comparison of Simulated vs. Ground Truth Birth Rate",
    x = "Date",
    y = "Birth Rate"
  ) +
  scale_color_manual(name = "Legend", values = c("Ground Truth" = "blue", "Simulated" = "red")) +
  theme_minimal()

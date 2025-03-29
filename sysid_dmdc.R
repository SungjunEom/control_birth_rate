###############################################
#
# Refactored Code: Aggregated Switched System Simulation
# over all intervals (2003.11 - 2024.06)
#
###############################################

######################
#
# Helper functions
#
######################

# Function to compute the rank of the controllability matrix
controllabilityRank <- function(A, B) {
  # Check that A is a square matrix
  if (nrow(A) != ncol(A)) {
    stop("Matrix A must be square.")
  }
  
  n <- nrow(A)
  # Initialize the controllability matrix with B
  ctrb <- B
  temp <- B
  
  # Build the controllability matrix [B, AB, A^2B, ..., A^(n-1)B]
  for (i in 1:(n - 1)) {
    temp <- A %*% temp
    ctrb <- cbind(ctrb, temp)
  }
  
  # Compute the rank using the QR decomposition
  rank <- qr(ctrb)$rank
  return(rank)
}

######################
#
# Preprocessing Data
#
######################

# Birth Rate (조출생률)
birth_rate <- read.csv('조출생률_smoothed_monthly_data.csv')[407:654,]
row.names(birth_rate) <- NULL

# Interest Rate (기준금리)
library(dplyr)
library(tidyr)
library(lubridate)
interest_rate <- read.csv('기준금리.csv', stringsAsFactors = FALSE)
interest_rate$변경일자 <- as.Date(interest_rate$변경일자)
interest_rate <- interest_rate %>%
  complete(변경일자 = seq(from = min(변경일자), to = max(변경일자), by = "month")) %>%
  arrange(변경일자) %>%
  fill(기준금리, .direction = "down")
interest_rate <- interest_rate %>%
  mutate(month = floor_date(변경일자, "month")) %>%
  group_by(month) %>%
  slice_max(변경일자, with_ties = FALSE) %>%
  ungroup()
interest_rate <- interest_rate[55:302,]

# Unemployment Rate (실업률)
unemploy_rate <- data.frame(t(read.csv('성_연령별_실업률_20250328224250.csv'))[-(1:2),])
colnames(unemploy_rate) <- c('실업률')
date_str <- rownames(unemploy_rate)
date_str_clean <- paste0(sub("^X", "", date_str), ".01")
dates <- as.Date(date_str_clean, format = "%Y.%m.%d")
unemploy_rate$date <- dates
rownames(unemploy_rate) <- NULL
unemploy_rate <- unemploy_rate[54:301,]

# Consumer Price Index (소비자물가지수)
price <- data.frame(t(read.csv('소비자물가지수_2020100__20250328223530.csv'))[-(1:2),])
colnames(price) <- c('소비자물가지수')
date_str <- rownames(price)
date_str_clean <- paste0(sub("^X", "", date_str), ".01")
dates <- as.Date(date_str_clean, format = "%Y.%m.%d")
price$date <- dates
rownames(price) <- NULL
price <- price[466:713,]

# Housing Price (주택매매가격)
house <- data.frame(t(read.csv('유형별_매매가격지수_20250328223335.csv'))[-(1:2),])
colnames(house) <- c('매매가격')
date_str <- rownames(house)
date_str_clean <- paste0(sub("^X", "", date_str), ".01")
dates <- as.Date(date_str_clean, format = "%Y.%m.%d")
house$date <- dates
rownames(house) <- NULL
house <- house[1:248,]

######################
#
# Construct Data Matrices
#
######################
# X_whole: State variables [birth, unemployment, price, house] (4 x 248)
X_whole <- t(data.matrix(data.frame(
  birth    = birth_rate$gaussian_smoothed,
  unemploy = unemploy_rate$실업률,
  price    = price$소비자물가지수,
  house    = house$매매가격
)))

# U: Input variable (interest rate) - note U has 247 columns for state transitions.
U <- t(data.matrix(data.frame(interest = interest_rate$기준금리)))[, 1:247]

######################
#
# Define Reusable System Identification Function
#
######################
system_identification <- function(start_idx, end_idx, X_whole, U) {
  m <- end_idx - start_idx + 1
  if(m < 2) {
    stop("Interval too short for system identification.")
  }
  
  X_interval <- X_whole[, start_idx:end_idx]
  U_interval <- U[start_idx:(end_idx - 1)]
  
  # Define state transition: X (time t) and Xnext (time t+1)
  X <- X_interval[, 1:(m - 1)]
  Xnext <- X_interval[, 2:m]
  
  # Stack state and input data
  O <- rbind(X, U_interval)
  
  Odecomposed <- svd(O)
  
  A <- Xnext %*% Odecomposed$v %*% diag(Odecomposed$d^(-1)) %*% t(Odecomposed$u[1:4, ])
  B <- Xnext %*% Odecomposed$v %*% diag(Odecomposed$d^(-1)) %*% Odecomposed$u[5, ]
  
  return(list(A = A, B = B))
}

# Define the five time intervals (by column indices in X_whole)
intervals <- list(
  "2003.11-2008.06" = c(1, 56),
  "2008.07-2013.06" = c(57, 116),
  "2013.07-2018.01" = c(117, 171),
  "2018.02-2023.09" = c(172, 239),
  "2023.10-2024.06" = c(240, 248)
)

######################
#
# Run System Identification for Each Interval
#
######################
# Store each interval's A and B matrices in a list
results <- lapply(names(intervals), function(interval_name) {
  idx <- intervals[[interval_name]]
  cat("Processing interval:", interval_name, "\n")
  res <- system_identification(start_idx = idx[1], end_idx = idx[2], X_whole = X_whole, U = U)
  return(res)
})
names(results) <- names(intervals)

######################
#
# Aggregate Simulation of the Switched System
#
######################
# Create a global date sequence (Nov 2003 to Jun 2024)
global_time_index <- seq.Date(from = as.Date("2003-11-01"), length.out = ncol(X_whole), by = "month")

# Initialize aggregated simulation as empty
aggregated_sim <- NULL
# Record switching time indices to add vertical lines later.
switch_indices <- c()

# Loop over intervals in chronological order.
for(interval_name in names(intervals)) {
  idx <- intervals[[interval_name]]
  # Number of time steps in the current interval
  T_interval <- idx[2] - idx[1] + 1
  
  # Extract the input for this interval (U indices align with state transitions)
  U_interval <- U[idx[1]:(idx[2]-1)]
  
  # Retrieve system matrices A and B for this interval
  A <- results[[interval_name]]$A
  B <- results[[interval_name]]$B
  
  # Allocate simulation matrix for this interval
  x_sim_interval <- matrix(0, nrow = nrow(X_whole), ncol = T_interval)
  # Use the ground truth initial condition for the current interval
  x_sim_interval[, 1] <- X_whole[, idx[1]]
  
  # Simulate the system for the interval: x(t+1) = A %*% x(t) + B * u(t)
  for (t in 1:(T_interval - 1)) {
    x_sim_interval[, t + 1] <- A %*% x_sim_interval[, t] + B * U_interval[t]
  }
  
  # Append the simulation for this interval (keep all columns)
  aggregated_sim <- cbind(aggregated_sim, x_sim_interval)
  
  # Record the global index of the interval's ending (switch point)
  switch_indices <- c(switch_indices, idx[2])
}



# Check that the aggregated simulation has the same number of time steps as the original data.
cat("Aggregated simulation steps:", ncol(aggregated_sim), "\n")
cat("Expected total steps:", ncol(X_whole), "\n")


######################
#
# Plot Aggregated Simulation vs. Ground Truth
#
######################
library(ggplot2)
# Create a data frame for plotting.
# Create a data frame for plotting.
sim_data <- data.frame(
  date = global_time_index,
  simulated = aggregated_sim[1, ],         # first row: birth rate
  ground_truth = X_whole[1, ],
  input = interest_rate$기준금리            # input series over the full period
)
# Define parameters for input shifting and expansion.
offset_value <- -2
input_scale_factor <- max(sim_data$ground_truth, na.rm = TRUE) / max(sim_data$input - offset_value, na.rm = TRUE)
expansion_factor <- 1.5  # Set >1 to expand both series

# Compute expanded series for ground truth and simulation.
expanded_ground_truth <- mean(sim_data$ground_truth, na.rm = TRUE) +
  expansion_factor * (sim_data$ground_truth - mean(sim_data$ground_truth, na.rm = TRUE))
expanded_simulated <- mean(sim_data$simulated, na.rm = TRUE) +
  expansion_factor * (sim_data$simulated - mean(sim_data$simulated, na.rm = TRUE))

# Create an interval label vector for each time index
sim_data$sim_interval <- factor(c(
  rep("2003.11-2008.06", 56), 
  rep("2008.07-2013.06", 60), 
  rep("2013.07-2018.01", 55), 
  rep("2018.02-2023.09", 68), 
  rep("2023.10-2024.06", 9)
))

# Build the plot:
p_total <- ggplot(sim_data, aes(x = date)) +
  # Ground truth with fixed color (will appear in the legend)
  geom_line(aes(y = expanded_ground_truth, color = "Ground Truth"), 
            size = 1, linetype = "solid") +
  # Simulated series colored by interval
  geom_line(aes(y = expanded_simulated, color = sim_interval), 
            size = 1, linetype = "dashed") +
  # Input series with fixed color (appears in legend)
  geom_line(aes(y = (input - offset_value) * input_scale_factor, color = "Input"), 
            size = 1, linetype = "dotted") +
  scale_y_continuous(
    name = "Expanded Birth Rate",
    sec.axis = sec_axis(~ . / input_scale_factor + offset_value, name = "Interest Rate")
  ) +
  labs(
    title = "Aggregated Switched System Simulation: Expanded Series, Shifted Input, and Intervals",
    x = "Date"
  ) +
  scale_color_manual(name = "Legend",
                     # Define colors for the fixed series and each interval:
                     values = c("Ground Truth"      = "blue",
                                "Input"             = "green",
                                "2003.11-2008.06"   = "red",
                                "2008.07-2013.06"   = "orange",
                                "2013.07-2018.01"   = "purple",
                                "2018.02-2023.09"   = "brown",
                                "2023.10-2024.06"   = "pink")) +
  theme_minimal()

# Add vertical lines at switching points (except the last one)
for(switch_idx in switch_indices[-length(switch_indices)]) {
  p_total <- p_total + geom_vline(xintercept = as.numeric(global_time_index[switch_idx]), 
                                  linetype = "dotdash", color = "gray")
}

print(p_total)


######################
#
# Check controllability
#
######################

rank_value <- controllabilityRank(results$`2023.10-2024.06`$A,results$`2023.10-2024.06`$B)
print(paste("Rank of the controllability matrix:", rank_value))


######################
#
# Simulate the future
#
######################
# ----- Future Simulation Using the "2023.10-2024.06" System -----

# Extract system matrices from the last interval ("2023.10-2024.06")
A_future <- results[["2023.10-2024.06"]][["A"]]
B_future <- results[["2023.10-2024.06"]][["B"]]

# Define the future simulation horizon (e.g., 24 months)
T_future <- 24

# Use the last observed state from the overall data as the initial condition.
# X_whole is 4 x 248, so we take the last column.
x0 <- X_whole[, ncol(X_whole)]  

# Define a vector of virtual inputs for the future.
# For example, use a constant virtual input equal to the mean interest rate.
# virtual_inputs <- rep(mean(interest_rate$기준금리, na.rm = TRUE), T_future)
# Alternatively, you could define a different pattern, e.g.:
# virtual_inputs <- seq(1.5, 2.5, length.out = T_future)
virtual_inputs <- rep(2.5, length.out = T_future)

# Initialize the simulation matrix for the future.
# We simulate T_future steps, so we need T_future+1 columns (including the initial state).
x_future <- matrix(0, nrow = length(x0), ncol = T_future + 1)
x_future[, 1] <- x0

# Simulate the future states using the model: x(t+1) = A_future %*% x(t) + B_future * u(t)
for(t in 1:T_future) {
  x_future[, t+1] <- A_future %*% x_future[, t] + B_future * virtual_inputs[t]
}

# Create a time vector for the future simulation.
# Use the last date in global_time_index (which corresponds to the last observed month)
# and add months to generate future dates.
future_start_date <- global_time_index[ncol(X_whole)] %m+% months(1)
future_dates <- seq.Date(from = future_start_date, by = "month", length.out = T_future + 1)

# Build a data frame for plotting.
# We only visualize the birth rate variable (first row of the state vector).
future_sim_data <- data.frame(
  date = future_dates,
  birth_rate_sim = x_future[1, ]
)

# Plot the simulated future birth rate.
library(ggplot2)
p_future <- ggplot(future_sim_data, aes(x = date, y = birth_rate_sim)) +
  geom_line(color = "red", size = 1) +
  labs(
    title = "Future Simulation: Birth Rate",
    x = "Date",
    y = "Simulated Birth Rate"
  ) +
  theme_minimal()

print(p_future)


###############################################
#
# Future Closed-Loop Simulation with Output Tracking,
# Input Saturation, and Recording of Inputs
# (Using the acker() function from the control package)
#
###############################################


# Install and load the control package (which provides the acker() function)
if (!require(control)) install.packages("control")
library(control)

# Choose desired closed-loop poles.
p_desired <- c(0.5, 0.4, 0.3, 0.2)

# Compute the state feedback gain K using Ackermann's formula.
K <- acker(A_future, B_future, p_desired)
cat("Designed state feedback gain K:\n")
print(K)

# Define the output matrix. We want to track the birth rate (first state).
C <- matrix(c(1, 0, 0, 0), nrow = 1)

# Compute the reference tracking gain Nbar.
# The closed-loop system is: x(t+1) = (A - B*K)x(t) + B * (Nbar * r_target)
# For steady state, we need r_target = C * (I - (A - B*K))^{-1} * B * (Nbar * r_target).
# Hence, Nbar = 1 / [C*(I - (A - B*K))^{-1}*B].
Nbar <- as.numeric(1 / (C %*% solve(diag(4) - (A_future - B_future %*% t(K))) %*% B_future))
cat("Reference tracking gain Nbar:\n")
print(Nbar)

# --- Closed-Loop Simulation Setup ---
# Assume X_whole is already defined from previous processing (a 4 x 248 matrix).
# Use the last observed state as the initial condition.
x0 <- X_whole[, ncol(X_whole)]  # last observed state

# Define the desired reference for the output (birth rate).
r_target <- 0.5  # Change this value to the desired birth rate

# Set the future simulation horizon (e.g., 24 months)
T_future <- 48

# Preallocate a matrix for the closed-loop simulation.
x_closed <- matrix(0, nrow = nrow(A_future), ncol = T_future + 1)
x_closed[, 1] <- x0

# Preallocate a vector to record the saturated control inputs.
u_record <- numeric(T_future)

# Simulation loop with input saturation (u constrained between 0 and 100)
for (t in 1:T_future) {
  # Compute control input: u(t) = -K*x(t) + Nbar * r_target.
  # Note: use -t(K) since K from acker() is returned as a row vector.
  u <- -t(K) %*% x_closed[, t] + Nbar * r_target
  # Saturate the input: ensure u is between 0 and 100.
  u_sat <- pmin(pmax(u, 0), 100)
  # Record the saturated input (as a scalar).
  u_record[t] <- u_sat[1,1]
  # Update the state using the saturated control input.
  x_closed[, t + 1] <- A_future %*% x_closed[, t] + B_future * u_sat[1,1]
}

# Create a future time vector.
future_start_date <- global_time_index[ncol(X_whole)] %m+% months(1)
future_dates <- seq.Date(from = future_start_date, by = "month", length.out = T_future + 1)

# Build a data frame to hold the simulation results.
# The birth_rate is the first state variable.
# For inputs, we prepend an NA for the initial time (since control starts at t=1).
future_closed_loop <- data.frame(
  date = future_dates,
  birth_rate = x_closed[1, ],
  input = c(NA, u_record)
)

# Determine a scaling factor to overlay the input on the same plot.
scale_factor_input <- max(future_closed_loop$birth_rate, na.rm = TRUE) / max(u_record, na.rm = TRUE)

# Plot the simulated birth rate and the saturated input.
library(ggplot2)
p_closed_loop <- ggplot(future_closed_loop, aes(x = date)) +
  geom_line(aes(y = birth_rate, color = "Birth Rate"), size = 1) +
  geom_line(aes(y = input * scale_factor_input, color = "Input"), 
            size = 1, linetype = "dotted") +
  scale_y_continuous(
    name = "Simulated Birth Rate",
    sec.axis = sec_axis(~ . / scale_factor_input, name = "Input")
  ) +
  labs(title = "Future Closed-Loop Simulation: Output Tracking & Input Saturation",
       x = "Date") +
  scale_color_manual(name = "Legend", 
                     values = c("Birth Rate" = "red", "Input" = "blue")) +
  theme_minimal()

print(p_closed_loop)



###############################################
#
# Optimal Control for Output Tracking via CVXR
#
###############################################

# Install and load the CVXR package if not already installed
if (!require(CVXR)) install.packages("CVXR")
library(CVXR)

# Define the output matrix (we track the birth rate, i.e., the first state)
C <- matrix(c(1, 0, 0, 0), nrow = 1)

# Define weight parameters (you can adjust these)
Q    <- 1      # weight on tracking error (stage cost)
R    <- 0.01   # weight on control effort (stage cost)
Q_f  <- 10     # terminal weight on tracking error

# Define the desired target for the birth rate
r_target <- 1

# Set the simulation (control) horizon (e.g., 24 steps/months)
T_horizon <- 24

# System dimension (n = 4)
n <- 4

# x0 is the initial state at time 0 (for instance, the last column of X_whole)
# (Assume X_whole is defined; here we use its last column as x0)
x0 <- X_whole[, ncol(X_whole)]

# Define optimization variables:
# x: state trajectory (n x (T_horizon+1))
# u: control input trajectory (length T_horizon)
x <- Variable(n, T_horizon + 1)
u <- Variable(T_horizon)

# Build the list of constraints:
constraints <- list()
# Initial condition:
constraints <- c(constraints, x[, 1] == x0)
# Dynamics and input constraints for t = 0,...,T_horizon-1:
for(t in 1:T_horizon) {
  constraints <- c(constraints,
                   x[, t + 1] == A_future %*% x[, t] + B_future * u[t],
                   u[t] >= 0,
                   u[t] <= 100)
}

# Define the cost function
# Note: x[,1] is x(0), x[,t+1] is x(t) for t=0,...,T_horizon.
cost <- 0
# Sum stage costs for t = 0,...,T_horizon-1:
for(t in 1:T_horizon) {
  # Tracking error: C*x[,t] - r_target; we use squared error multiplied by Q.
  cost <- cost + Q * (C %*% x[, t] - r_target)^2 + R * (u[t])^2
}
# Terminal cost at t = T_horizon:
cost <- cost + Q_f * (C %*% x[, T_horizon + 1] - r_target)^2

# Formulate the optimization problem:
problem <- Problem(Minimize(cost), constraints)

# Solve the problem
result <- solve(problem)
cat("Optimal cost:", result$value, "\n")

# Extract the optimal trajectories
x_opt <- result$getValue(x)
u_opt <- result$getValue(u)

# Create a future time vector.
# future_start_date: one month after the last observed date in global_time_index.
future_start_date <- global_time_index[ncol(X_whole)] %m+% months(1)
future_dates <- seq.Date(from = future_start_date, by = "month", length.out = T_horizon + 1)

# Build a data frame for plotting the birth rate (first state variable) and the input.
opt_data <- data.frame(
  date       = future_dates,
  birth_rate = as.numeric(x_opt[1, ]),      # first state = birth rate
  input      = c(NA, as.numeric(u_opt))       # input: no input at time 0
)

# For visualization, we overlay the birth rate and the input.
# We use a scaling factor to put the input on a similar scale.
scale_factor_input <- max(opt_data$birth_rate, na.rm = TRUE) / max(opt_data$input, na.rm = TRUE)

library(ggplot2)
p_opt <- ggplot(opt_data, aes(x = date)) +
  geom_line(aes(y = birth_rate, color = "Birth Rate"), size = 1) +
  geom_line(aes(y = input * scale_factor_input, color = "Input"), 
            size = 1, linetype = "dotted") +
  scale_y_continuous(
    name = "Birth Rate",
    sec.axis = sec_axis(~ . / scale_factor_input, name = "Input")
  ) +
  labs(title = "Optimal Output Tracking: Birth Rate and Control Input",
       x = "Date") +
  scale_color_manual(name = "Legend", 
                     values = c("Birth Rate" = "red", "Input" = "blue")) +
  theme_minimal()

print(p_opt)


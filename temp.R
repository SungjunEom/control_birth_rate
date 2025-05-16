###########################
# Install and load the CVXR package if not already installed
if (!require(CVXR)) install.packages("CVXR")
library(CVXR)

# Define the output matrix (we track the birth rate, i.e., the first state)
C <- matrix(c(1, 0, 0, 0), nrow = 1)

# Define weight parameters (you can adjust these)
Q    <- 10      # weight on tracking error (stage cost)
R    <- 0.01   # weight on control effort (stage cost)
Q_f  <- 10     # terminal weight on tracking error

# Set the simulation (control) horizon (e.g., 48 steps/months)
T_horizon <- 48

# System dimension (n = 4)
n <- 4

# x0 is the initial state at time 0 (for instance, the last column of X_whole)
x0 <- X_whole[, ncol(X_whole)]

# Define a time-varying reference trajectory as an increasing function.
# For example, we can let it increase linearly from 0.5 to 1.5.
r_target_vec <- seq(0.4, 0.5, length.out = T_horizon + 1)

# Define optimization variables:
# x: state trajectory (n x (T_horizon+1))
# u: control input trajectory (length T_horizon)
x <- Variable(n, T_horizon + 1)
u <- Variable(T_horizon)

# Build the list of constraints:
constraints <- list()
# Initial condition:
constraints <- c(constraints, x[, 1] == x0)
# Dynamics and input constraints for t = 1,...,T_horizon:
for(t in 1:T_horizon) {
  constraints <- c(constraints,
                   x[, t + 1] == A_future %*% x[, t] + B_future * u[t],
                   u[t] >= 0,
                   u[t] <= 100)
}

# Define the cost function using the time-varying reference.
# Here the stage cost at time t compares the output C*x[, t] with r_target_vec[t]
cost <- 0
for(t in 1:T_horizon) {
  cost <- cost + Q * (C %*% x[, t] - r_target_vec[t])^2 + R * (u[t])^2
}
# Terminal cost at t = T_horizon+1 uses the last element of the reference trajectory.
cost <- cost + Q_f * (C %*% x[, T_horizon + 1] - r_target_vec[T_horizon + 1])^2

# Formulate and solve the optimization problem:
problem <- Problem(Minimize(cost), constraints)
result <- solve(problem)
cat("Optimal cost:", result$value, "\n")

# Extract the optimal trajectories:
x_opt <- result$getValue(x)
u_opt <- result$getValue(u)


# Create a future time vector.
# future_start_date: one month after the last observed date in global_time_index.
future_start_date <- global_time_index[ncol(X_whole)] %m+% months(1)
future_dates <- seq.Date(from = future_start_date, by = "month", length.out = T_horizon + 1)

#########
# Plot only birth-rate vs. inputs
###########
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
  geom_line(aes(y = birth_rate, color = "월별 가상 조출생률"), size = 1) +
  geom_line(aes(y = input * scale_factor_input, color = "기준금리"), 
            size = 1, linetype = "dotted") +
  scale_y_continuous(
    name = "월별 가상 조출생률",
    sec.axis = sec_axis(~ . / scale_factor_input, name = "기준금리")
  ) +
  labs(title = "출력 추종 최적 제어: 점진적 감소",
       x = "날짜") +
  scale_color_manual(name = "", 
                     values = c("월별 가상 조출생률" = "red", "기준금리" = "blue")) +
  theme_minimal()

print(p_opt)

library(lubridate)
library(ggplot2)

# --- Part 1: Prepare the Historical Ground Truth Data ---

# Ground truth birth rate from the original data (first row of X_whole)
historical_ground_truth <- X_whole[1, ]
# Historical dates (assumed to cover the entire period, e.g., 2003-11 to 2024-06)
historical_dates <- global_time_index  # length should match the number of columns in X_whole

# Historical control input (actual interest rate over the historical period)
historical_input <- interest_rate$기준금리[1:length(historical_ground_truth)]

# --- Part 2: Prepare the Future Optimal Control Simulation Data ---

# Optimal control simulation using CVXR produced x_opt and u_opt.
# Note: x_opt includes the initial state (same as last historical value),
# so for the future portion, we remove the duplicate initial state.
future_birth_rate <- as.numeric(x_opt[1, -1])  # optimal birth rate simulation for future time steps
future_input      <- as.numeric(u_opt)         # optimal control inputs over the horizon

# Future dates: start one month after the last historical date
future_dates <- seq.Date(from = historical_dates[length(historical_dates)] %m+% months(1),
                         by = "month", length.out = length(future_birth_rate))

# --- Part 3: Combine Historical and Future Data ---

# Combined birth rate: ground truth for historical, then future optimal control simulation
combined_birth_rate <- c(historical_ground_truth, future_birth_rate)
# Combined input: historical interest rate for the historical period and optimal control input for future
combined_input <- c(historical_input, future_input)
# Combined date vector:
combined_dates <- c(historical_dates, future_dates)

# Create a phase label to distinguish historical from future simulation
phase <- c(rep("조출생률 실제값", length(historical_ground_truth)),
           rep("최적 제어 시뮬레이션", length(future_birth_rate)))

# Build a combined data frame for plotting.
combined_df <- data.frame(
  date = combined_dates,
  birth_rate = combined_birth_rate,
  input = combined_input,
  phase = phase
)

# --- Part 4: Plot the Combined Data ---

# Compute a scaling factor to overlay the interest rate on the same plot.
scale_factor_input <- max(combined_birth_rate, na.rm = TRUE) / max(combined_input, na.rm = TRUE)

p_combined_gt <- ggplot(combined_df, aes(x = date)) +
  # Plot the birth rate trajectory, colored by phase.
  geom_line(aes(y = birth_rate, color = phase), size = 1) +
  # Overlay the control input (scaled) with a fixed color.
  geom_line(aes(y = input * scale_factor_input), color = "darkgreen", linetype = "dotted", size = 1) +
  scale_y_continuous(
    name = "월별 가상 조출생률",
    sec.axis = sec_axis(~ . / scale_factor_input, name = "금리 (%)")
  ) +
  labs(title = "출력 추종 최적 제어: 점진적 증가",
       x = "날짜", color = "") +
  theme_minimal()

print(p_combined_gt)

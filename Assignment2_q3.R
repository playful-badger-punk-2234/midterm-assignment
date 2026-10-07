## Assignment 2 Q3 ##
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(123)

# -------------------------------------------------------------------------
# Define Design Configurations
# -------------------------------------------------------------------------
mu_true <- 5.0
sigma_tot_sq <- 1.0
rho_vals <- c(0.00, 0.45, 0.95)

# Unbalanced group setups
g_sizes_small <- c(4, 5, 6, 7, 8, 10) # G = 6, N = 40
set.seed(42)
g_sizes_large <- c(14, 15, 17, 18, 19, 20, 21, 22, 24, 25,
                   15, 16, 18, 19, 21, 22, 23, 25, 26, 20) # G = 20, N = 400

scenarios <- expand.grid(
  Size = c("Small", "Large"),
  Rho_True = rho_vals,
  stringsAsFactors = FALSE
)

# Simulator function
simulate_lmm_data <- function(group_sizes, mu, rho, sig_tot = 1.0) {
  sig_a_sq   <- rho * sig_tot
  sig_eps_sq <- (1 - rho) * sig_tot
  
  G <- length(group_sizes)
  a_i <- rnorm(G, mean = 0, sd = sqrt(sig_a_sq))
  
  group_vec <- rep(1:G, times = group_sizes)
  a_expanded <- rep(a_i, times = group_sizes)
  eps_ij <- rnorm(sum(group_sizes), mean = 0, sd = sqrt(sig_eps_sq))
  
  data.frame(
    group = factor(group_vec),
    y = mu + a_expanded + eps_ij
  )
}

# -------------------------------------------------------------------------
# Part (i): Single Realization per Scenario
# -------------------------------------------------------------------------
single_run_results <- list()
sim_datasets <- list()

for (k in 1:nrow(scenarios)) {
  sz  <- scenarios$Size[k]
  rho <- scenarios$Rho_True[k]
  g_sizes <- if (sz == "Small") g_sizes_small else g_sizes_large
  
  df_sim <- simulate_lmm_data(g_sizes, mu_true, rho)
  sim_datasets[[paste(sz, rho, sep = "_")]] <- df_sim %>%
    mutate(Scenario = paste0(sz, " (rho = ", rho, ")"))
  
  fit <- lmer(y ~ 1 + (1 | group), data = df_sim)
  
  # Point estimates
  mu_est <- fixef(fit)["(Intercept)"]
  se_mu  <- summary(fit)$coefficients["(Intercept)", "Std. Error"]
  ci_mu  <- mu_est + c(-1, 1) * 1.96 * se_mu
  
  vc <- as.data.frame(VarCorr(fit))
  sig_a_est_sq   <- ifelse(any(vc$grp == "group"), vc$vcov[vc$grp == "group"], 0)
  sig_eps_est_sq <- vc$vcov[vc$grp == "Residual"]
  rho_est        <- sig_a_est_sq / (sig_a_est_sq + sig_eps_est_sq)
  
  single_run_results[[k]] <- data.frame(
    Size          = sz,
    N             = sum(g_sizes),
    G             = length(g_sizes),
    True_Rho      = rho,
    True_Mu       = mu_true,
    Est_Mu        = round(mu_est, 3),
    CI_Lower      = round(ci_mu[1], 3),
    CI_Upper      = round(ci_mu[2], 3),
    True_Sig_a2   = round(rho * sigma_tot_sq, 3),
    Est_Sig_a2    = round(sig_a_est_sq, 3),
    True_Sig_eps2 = round((1 - rho) * sigma_tot_sq, 3),
    Est_Sig_eps2  = round(sig_eps_est_sq, 3),
    Est_Rho       = round(rho_est, 3)
  )
}
table_part_i <- do.call(rbind, single_run_results)

# -------------------------------------------------------------------------
# Part (ii): B = 100 Repetitions
# -------------------------------------------------------------------------
B <- 100
mc_results <- list()
all_mc_draws <- list()

for (k in 1:nrow(scenarios)) {
  sz  <- scenarios$Size[k]
  rho <- scenarios$Rho_True[k]
  g_sizes <- if (sz == "Small") g_sizes_small else g_sizes_large
  
  mu_hats <- numeric(B)
  
  for (b in 1:B) {
    df_b <- simulate_lmm_data(g_sizes, mu_true, rho)
    fit_b <- suppressMessages(lmer(y ~ 1 + (1 | group), data = df_b))
    mu_hats[b] <- fixef(fit_b)["(Intercept)"]
  }
  
  avg_mu <- mean(mu_hats)
  var_mu <- var(mu_hats)
  se_avg <- sqrt(var_mu / B)
  ci_avg_mu <- avg_mu + c(-1, 1) * 1.96 * se_avg
  
  mc_results[[k]] <- data.frame(
    Size             = sz,
    True_Rho         = rho,
    Avg_Est_Mu_B100  = round(avg_mu, 3),
    Var_Est_Mu_B100  = round(var_mu, 4),
    CI_Avg_Lower     = round(ci_avg_mu[1], 3),
    CI_Avg_Upper     = round(ci_avg_mu[2], 3)
  )
  
  all_mc_draws[[k]] <- data.frame(
    Estimate = mu_hats,
    Size     = sz,
    Rho      = factor(paste0("rho == ", rho))
  )
}
table_part_ii <- do.call(rbind, mc_results)

# Assemble Master Table 1
Table_1 <- left_join(table_part_i, table_part_ii, by = c("Size", "True_Rho"))

combined_data <- bind_rows(sim_datasets)

ggplot(combined_data, aes(x = group, y = y)) +
  geom_hline(yintercept = mu_true, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.3, color = "slategray") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.8, color = "midnightblue") +
  facet_wrap(~ Scenario, scales = "free_x", ncol = 3) +
  labs(
    title = "Figure 1: Simulated Hierarchical Data Across Six Scenarios",
    subtitle = "Diamonds = Group means; Gray dots = Individual observations; Dashed red line = True mean (5.0)",
    x = "Group ID",
    y = expression(y[ij])
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

mc_plot_data <- bind_rows(all_mc_draws)

ggplot(mc_plot_data, aes(x = Estimate, fill = Size)) +
  geom_density(alpha = 0.45) +
  geom_vline(xintercept = mu_true, linetype = "dashed", color = "black") +
  facet_wrap(~ Rho, labeller = label_parsed, ncol = 3) +
  labs(
    title = "Figure 2: Empirical Sampling Distributions of Fixed Effect Mean (B = 100)",
    subtitle = "Sampling variance expands drastically as rho increases, especially for small G",
    x = expression(hat(mu)),
    y = "Density",
    fill = "Sample Size"
  ) +
  scale_fill_manual(values = c("Small" = "coral2", "Large" = "royalblue3")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

##########################################################################

# Ensure Size and Rho have clean factor ordering for display
mc_plot_data$Size <- factor(mc_plot_data$Size, levels = c("Small", "Large"))
mc_plot_data$Rho_Label <- factor(
  mc_plot_data$Rho,
  levels = c("rho == 0", "rho == 0.45", "rho == 0.95"),
  labels = c("rho == 0.00", "rho == 0.45", "rho == 0.95")
)

ggplot(mc_plot_data, aes(x = Estimate)) +
  # Histogram with density scaling
  geom_histogram(aes(y = after_stat(density), fill = Size), 
                 bins = 18, color = "white", alpha = 0.55) +
  # Smoothed empirical density
  geom_density(color = "black", linewidth = 0.8) +
  # Reference line at true mean mu = 5.0
  geom_vline(xintercept = 5.0, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  # 2 rows (Size) by 3 columns (Rho)
  facet_grid(Size ~ Rho_Label, labeller = label_parsed, scales = "free_y") +
  scale_fill_manual(values = c("Small" = "#E66101", "Large" = "#5E3C99")) +
  labs(
    title = "Figure 2: Empirical Sampling Distributions of Mean Estimator (B = 100)",
    subtitle = "Dashed red line indicates true mean (mu = 5.0)",
    x = expression(hat(mu)),
    y = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 11)
  )


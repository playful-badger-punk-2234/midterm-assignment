## Assignment 2 ##
# Load Libaries
library(lme4)
library(multcomp)
library(lmerTest)
library(emmeans)
library(patchwork)
library(ggplot2)

# setwd("C:/Users/lachl/Documents/Uni - DTU/02429 - Analysis of Correlated Data MLP/Data")
# remove(list=ls())

# Q1


# Create data
df <- data.frame(
  plant = rep(1:4, each = 3),
  leaf  = rep(1:3, times = 4),
  det_1 = c(3.28, 3.52, 2.88, 2.46, 1.87, 2.19, 2.77, 3.74, 2.55, 3.78, 4.07, 3.31),
  det_2 = c(3.09, 3.48, 2.80, 2.44, 1.92, 2.19, 2.66, 3.44, 2.55, 3.87, 4.12, 3.31)
)
# Add row totals matching X_ij
df$total <- df$det_1 + df$det_2

# Alternative expression of df
df <- data.frame(
  plant = factor(rep(1:4, each = 6)),
  leaf  = factor(rep(rep(1:3, each = 2), times = 4)),
  calcium = c(
    3.28, 3.09, 3.52, 3.48, 2.88, 2.80, # Plant 1
    2.46, 2.44, 1.87, 1.92, 2.19, 2.19, # Plant 2
    2.77, 2.66, 3.74, 3.44, 2.55, 2.55, # Plant 3
    3.78, 3.87, 4.07, 4.12, 3.31, 3.31  # Plant 4
  )
)

# Data Visulisation
ggplot(df, aes(x = leaf, y = calcium)) +
  geom_line(aes(group = leaf), color = "gray60", linewidth = 0.8) +
  geom_point(aes(color = plant), size = 3) +
  facet_wrap(~ plant, labeller = label_both, nrow = 1) +
  labs(
    title = "Calcium Determinations by Leaf and Plant",
    subtitle = "Connected points represent duplicate determinations on the same leaf",
    x = "Leaf ID (within plant)",
    y = "Calcium Concentration (%)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

library(tidyr)

# Reshape wide to compare det_1 vs det_2
diff_df <- df %>%
  mutate(determination = rep(c("d1", "d2"), 12)) %>%
  pivot_wider(names_from = determination, values_from = calcium) %>%
  mutate(
    diff = d1 - d2,
    mean_val = (d1 + d2) / 2
  )

ggplot(diff_df, aes(x = mean_val, y = diff, color = plant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size = 3) +
  labs(
    title = "Measurement Error: Difference between Duplicate Determinations",
    subtitle = "Points scattered symmetrically around 0 indicate unbiased analytical error",
    x = "Mean Determination Value (%)",
    y = "Determination 1 - Determination 2 (%)",
    color = "Plant"
  ) +
  theme_minimal()

# part i)
# Fit the Linear Mixed Model
# In lme4, nesting is specified as (1 | plant / leaf)
fit <- lmer(calcium ~ 1 + (1 | plant/leaf), data = df)
summary(fit)

# Extract variance components
var_comp <- as.data.frame(VarCorr(fit))
var_comp[, c("grp", "vcov", "sdcor")]

# Using aov with nested error terms
aov_fit <- aov(calcium ~ 1 + Error(plant/leaf), data = df)
summary(aov_fit)

# part ii)


# part iii)
# a) Direct formula using coef summary
se_mean <- summary(fit)$coefficients["(Intercept)", "Std. Error"]
mu_hat  <- fixef(fit)["(Intercept)"]

# Normal Wald CI (default in large samples)
ci_wald_z <- mu_hat + c(-1, 1) * qnorm(0.975) * se_mean
ci_wald_z

# Or using confint with method = "Wald"
confint(fit, parm = "(Intercept)", method = "Wald")

# Note: With only 4 plants (3 degrees of freedom for the top level),
# a small-sample t-distribution (df = 3) is often preferred:
ci_wald_t <- mu_hat + c(-1, 1) * qt(0.975, df = 3) * se_mean
ci_wald_t

# b) Profile Confidence Interval
ci_profile <- confint(fit, parm = "(Intercept)", method = "profile")
ci_profile

# c) Bootstrap Confidence Interval (set seed for reproducibility)
set.seed(42)
ci_boot <- confint(fit, parm = "(Intercept)", method = "boot", nsim = 1000)
ci_boot




















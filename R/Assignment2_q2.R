## Assignment 2 Question 2 ##
# Load necessary libraries
library(lme4)
library(lmerTest) # Provides p-values for lmer models
library(ggplot2)
library(dplyr)
library(performance) # For diagnostics

setwd("C:/Users/lachl/Documents/Uni - DTU/02429 - Analysis of Correlated Data MLP/Assignment 2")
remove(list=ls())

load("Pups_Dempster.RData")
pups <- Pups_Dempster

# Assuming the data is loaded into a dataframe named 'pups'
pups <- pups %>%
  mutate(
    Dose = as.factor(Dose),
    Sex = as.factor(Sex),
    Dam = as.factor(Dam),
    littersize = as.numeric(littersize),
    weight = as.numeric(weight)
  )

# Exploratory Visualizations
# 1. Treatment effect on weight
ggplot(pups, aes(x = Dose, y = weight, fill = Dose)) +
  geom_boxplot() +
  labs(title = "Pup Weight by Dose", x = "Dose", y = "Weight") +
  theme_minimal()

# 2. Littersize vs Weight 
ggplot(pups, aes(x = littersize, y = weight, color = Dose)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Pup Weight vs Litter Size", x = "Litter Size", y = "Weight") +
  theme_minimal()



# Fit the model
model_lmm <- lmer(weight ~ Dose + Sex + littersize + (1 | Dam), data = pups)

# View the summary
summary(model_lmm)


# Extract residuals and random effects
res <- residuals(model_lmm)
ranef_dam <- ranef(model_lmm)$Dam[[1]]

# 1. Homoscedasticity: Residuals vs Fitted values
plot(model_lmm, main = "Residuals vs Fitted")

# 2. Normality of Residuals
qqnorm(res, main = "Q-Q Plot of Residuals")
qqline(res, col = "red")

# 3. Normality of Random Effects
qqnorm(ranef_dam, main = "Q-Q Plot of Dam Random Effects")
qqline(ranef_dam, col = "blue")

# Alternative: Use the performance package for an automated comprehensive check
# check_model(model_lmm)

# ANOVA table for fixed effects (Type III Wald F-tests)
anova(model_lmm, type = 3, ddf = "Satterthwaite")

# Calculate 95% Confidence Intervals for the fixed effects
confint(model_lmm, method = "Wald")



# Extract variance components
var_comps <- as.data.frame(VarCorr(model_lmm))
sigma_dam_sq <- var_comps$vcov[1]
sigma_eps_sq <- var_comps$vcov[2]

# Calculate Intraclass Correlation Coefficient (ICC)
icc <- sigma_dam_sq / (sigma_dam_sq + sigma_eps_sq)
print(paste("ICC:", round(icc, 3)))






























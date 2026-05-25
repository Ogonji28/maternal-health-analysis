# =============================================================================
# Multilevel Modelling of Underutilization using DHS 2014 & 2022 Data
# Method: lmer() from lme4 package | REML & ML estimation
# Based on Module 5: Introduction to Multilevel Modelling (CMM, Bristol)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. INSTALL & LOAD REQUIRED PACKAGES
# -----------------------------------------------------------------------------
# Install once if not already installed:
#install.packages("lme4")
# install.packages("lmerTest")   # Adds p-values to lmer output
# install.packages("dplyr")
# install.packages("haven")      # For reading Stata .dta / .sav DHS files
# install.packages("tidyverse")
#install.packages("reformulas")

#library(Matrix)
library(lme4)
#library(lmerTest)   # Extends lmer with Satterthwaite p-values
library(haven)
library(dplyr)
library(tidyverse)
library(lattice)

# -----------------------------------------------------------------------------
# 2. LOAD DATA
# -----------------------------------------------------------------------------
# Option A: If your DHS data is in Stata format (.dta)
 #mydata <- read_dta("C:\Users\norah.ogonji_evidenc\Desktop\Norah_Learning\Strathmore_Semester Study\MY MAsters Project\DHS Data\2-Clean Data\KDHS_pooled_2014_2022_final.dta")

# Option B: If already saved as CSV
# mydata <- read.csv("your_dhs_file.csv")

# Option C: Load and combine 2014 and 2022 separately then merge
 dhs2014 <- read_dta("C:/Users/norah.ogonji_evidenc/Desktop/Norah_Learning/Strathmore_Semester Study/MY MAsters Project/DHS Data/2-Clean Data/KDHS2014_analysis_final.dta")
 
 dhs2022 <- read_dta("C:/Users/norah.ogonji_evidenc/Desktop/Norah_Learning/Strathmore_Semester Study/MY MAsters Project/DHS Data/2-Clean Data/KDHS2022_analysis_final.dta")
 dhs2014$survey_year <- 2014
 dhs2022$survey_year <- 2022
 mydata <- bind_rows(dhs2014, dhs2022)
 
 


 dhs2014 <- dhs2014 %>%
   mutate(v024 = stringr::str_to_title(as.character(as_factor(v024))))
 
 dhs2022 <- dhs2022 %>%
   mutate(v024 = stringr::str_to_title(as.character(as_factor(v024))))
 
 # Now merge
 mydata <- bind_rows(dhs2014, dhs2022)
 
 # Verify — should now show 8 regions with counts for BOTH years
 table(mydata$v024, mydata$survey_year)
 
 

# Quick check of the data structure:
#str(mydata)
#head(mydata)

# -----------------------------------------------------------------------------
# 3. DATA PREPARATION
# -----------------------------------------------------------------------------
# Ensure key variables are correctly coded:

# Underutilization: binary outcome (1 = Yes, 0 = No)
mydata$underutilized <- as.numeric(mydata$underutilized)  # Adjust var name

# Survey year: create indicator variable
mydata$year2022 <- ifelse(mydata$survey_year == 2022, 1, 0)

# Cluster/PSU identifier (Level 2 grouping variable in DHS)
# In DHS, 'v001' is the cluster number (PSU) — adjust to your variable name
#mydata$cluster_id <- mydata$cluster_id   # Replace with your actual cluster var

# Check missingness
colSums(is.na(mydata[, c("underutilized", "cluster_id", "survey_year")]))

# Drop missing on outcome or cluster
mydata_clean <- mydata %>%
  filter(!is.na(underutilized), !is.na(cluster_id))

# -----------------------------------------------------------------------------
# 4. NULL MODEL (No Predictors) — Partitioning Variance
# Model: underutilized_ij = β0 + u0j + eij
# This tests whether clustering (PSU) explains variance in underutilization
# REML = TRUE is the default and preferred for variance estimation
# -----------------------------------------------------------------------------
null_model_REML <- lmer(
  underutilized ~ 1 + (1 | cluster_id),
  data   = mydata_clean,
  REML   = TRUE   # <-- REML estimation for unbiased variance components
)

summary(null_model_REML)

# --- Variance Partition Coefficient (VPC / ICC) ---
# Tells you: what % of variation in underutilization is between clusters
var_components <- as.data.frame(VarCorr(null_model_REML))
var_cluster    <- var_components$vcov[1]   # Level 2 (cluster) variance
var_residual   <- var_components$vcov[2]   # Level 1 (individual) variance
total_var      <- var_cluster + var_residual

VPC <- var_cluster / total_var
cat("\n--- Variance Partition Coefficient (VPC/ICC) ---\n")
cat(sprintf("Between-cluster variance : %.4f\n", var_cluster))
cat(sprintf("Within-cluster variance  : %.4f\n", var_residual))
cat(sprintf("VPC (ICC)                : %.4f (%.1f%%)\n", VPC, VPC * 100))

# -----------------------------------------------------------------------------
# 5. TESTING REML vs ML — Why it Matters
# -----------------------------------------------------------------------------
# REML  = TRUE  → Best for estimating variance components (default)
# REML  = FALSE → Required when comparing fixed-effects models via LRT

null_model_ML <- lmer(
  underutilized ~ 1 + (1 | cluster_id),
  data = mydata_clean,
  REML = FALSE   # <-- ML for model comparison
)

# Compare fit statistics
cat("\n--- Model Fit Comparison: REML vs ML ---\n")
cat("REML AIC:", AIC(null_model_REML), "\n")
cat("ML   AIC:", AIC(null_model_ML),   "\n")

# -----------------------------------------------------------------------------
# 6. RANDOM INTERCEPT MODEL — Adding Survey Year
# Model: underutilized_ij = β0 + β1(year2022) + u0j + eij
# Clusters have different baseline underutilization, but year effect is fixed
# Use ML (REML=FALSE) for comparing fixed-effect structures
# -----------------------------------------------------------------------------
model_year_ML <- lmer(
  underutilized ~ year2022 + (1 | cluster_id),
  data = mydata_clean,
  REML = FALSE
)

summary(model_year_ML)

# Likelihood Ratio Test: Does adding survey year improve fit?
anova(null_model_ML, model_year_ML)

# Refit with REML for final coefficient interpretation
model_year_REML <- lmer(
  underutilized ~ year2022 + (1 | cluster_id),
  data = mydata_clean,
  REML = TRUE
)

summary(model_year_REML)

# -----------------------------------------------------------------------------
# 7. RANDOM INTERCEPT MODEL — Adding Individual-Level Covariates
# Common DHS variables (adjust names to match your dataset):
#   v106 = education level
#   v190 = wealth index
#   v013 = age group
#   v025 = urban/rural
# -----------------------------------------------------------------------------
model_covariates_ML <- lmer(
  underutilized ~ year2022 +
    factor(v106) +    # Education level
    factor(v190) +    # Wealth quintile
    factor(v013) +    # Age group
    v025 +            # Urban/Rural (1=urban, 2=rural)
    (1 | cluster_id),
  data = mydata_clean,
  REML = FALSE
)

summary(model_covariates_ML)

# LRT: Does adding covariates improve over year-only model?
anova(model_year_ML, model_covariates_ML)

# Final REML refit for reporting
model_covariates_REML <- lmer(
  underutilized ~ year2022 +
    factor(v106) +
    factor(v190) +
    factor(v013) +
    v025 +
    (1 | cluster_id),
  data = mydata_clean,
  REML = TRUE
)

summary(model_covariates_REML)

# -----------------------------------------------------------------------------
# 8. RANDOM SLOPE MODEL — Does Year Effect Vary Across Clusters?
# Model: underutilized_ij = β0 + β1(year2022) + u0j + u1j(year2022) + eij
# Allows the 2014→2022 change to differ by cluster
# -----------------------------------------------------------------------------
model_randomslope_ML <- lmer(
  underutilized ~ year2022 + (1 + year2022 | cluster_id),
  data    = mydata_clean,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")   # Helps with convergence
)

summary(model_randomslope_ML)

# LRT: Is the random slope significant? (Compare to random intercept only)
anova(model_year_ML, model_randomslope_ML)

# Final REML refit
model_randomslope_REML <- lmer(
  underutilized ~ year2022 + (1 + year2022 | cluster_id),
  data    = mydata_clean,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(model_randomslope_REML)

# Correlation between intercept and slope (from random effects output):
# Negative correlation → clusters with high baseline underutilization
#   improved more (converging toward the mean over time)
# Positive correlation → high-baseline clusters stayed high

# -----------------------------------------------------------------------------
# 9. EXTRACT & DISPLAY CLUSTER-LEVEL RESIDUALS (School effects equivalent)
# -----------------------------------------------------------------------------
# Level 2 residuals: each cluster's deviation from overall mean
cluster_residuals <- ranef(model_covariates_REML)$cluster_id
cluster_residuals$cluster_id <- rownames(cluster_residuals)
colnames(cluster_residuals)[1] <- "u0j"

# Sort and view top and bottom clusters
cluster_residuals <- cluster_residuals %>% arrange(desc(u0j))
cat("\nTop 10 clusters (highest underutilization above average):\n")
print(head(cluster_residuals, 10))
cat("\nBottom 10 clusters (most below average underutilization):\n")
print(tail(cluster_residuals, 10))

# Caterpillar plot of cluster residuals
dotplot(ranef(model_covariates_REML, condVar = TRUE))

# -----------------------------------------------------------------------------
# 10. MODEL SUMMARY TABLE (for reporting)
# -----------------------------------------------------------------------------
# Use lmerTest for p-values in summary
summary(model_covariates_REML)   # Gives Satterthwaite df + p-values via lmerTest

# Optionally use broom.mixed for tidy output:
# install.packages("broom.mixed")
# library(broom.mixed)
# tidy(model_covariates_REML, effects = "fixed", conf.int = TRUE)

# -----------------------------------------------------------------------------
# 11. DECISION GUIDE: REML vs ML
# -----------------------------------------------------------------------------
# +---------------------------+-------------+----------------------------------+
# | Purpose                   | Use REML?   | Reason                           |
# +---------------------------+-------------+----------------------------------+
# | Estimate variance comps   | YES (TRUE)  | REML gives unbiased estimates    |
# | Compare fixed effects     | NO (FALSE)  | LRT only valid under ML          |
# | Compare random effects    | YES (TRUE)  | REML LRT valid for random parts  |
# | Final reported model      | YES (TRUE)  | Convention; better var estimates |
# +---------------------------+-------------+----------------------------------+

cat("\nScript complete. Adjust variable names to match your DHS dataset.\n")

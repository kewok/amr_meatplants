library(mgcv)
library(dplyr)
library(tidyr)

# Preliminaries, copied from multiMethod_comparison
data <- read.csv('../Datasets2024/April11_2025/seasonal_amr_noGuam.csv')

# Concatenate the columns into a single column for each factor
# https://stackoverflow.com/a/71590918/1790399

df  <- data |>
  mutate(ID = row_number()) |>
  pivot_longer(cols = c("Large","Small","Very_Small"), 
  names_to = "Size") |>
  filter(value == 1)

# Get rid of the "value" column 
df <- df %>% select(-value)

df  <- df |>
  mutate(ID = row_number()) |>
  pivot_longer(cols = c("Animal_Chicken_Broiler___Young_Chicken_Carcass_Rinse","Product_Raw_Ground__Comminuted_or_Otherwise_Nonintact_Chicken","Product_Raw_Intact_Chicken","Product_Raw_Intact_Beef","Product_Raw_Ground__Comminuted_or_Otherwise_Nonintact_Beef","Product_Raw_Ground__Comminuted_or_Otherwise_Nonintact_Pork","Product_Raw_Intact_Pork","Product_Raw_Ground__Comminuted_or_Otherwise_Nonintact_Turkey","Animal_Turkey_Turkey_Carcass_Sponge"), 
  names_to = "Product") |>
  filter(value == 1)

# Get rid of the "value" column 
df <- df %>% select(-value)

library(GWmodel3)
library(sp)
spdf <- df
coordinates(spdf) <- ~longitude + latitude
provoliko = "ESRI:54009"
proj4string(spdf) <- provoliko
sp_df <- sf::st_as_sf(spdf)

ols <- readRDS("OLS_model.rds")
basic_GAM <- readRDS("basic_GAM_Poisson_model.rds")
spatiotemporal_GAM <- readRDS("GAM_model.rds")
classic_gwr <- readRDS("basicGWR.rds")
classic_gtwr <- readRDS("basicGTWR.rds")

# R^2
summary(ols)$r.squared
1-sum(residuals(basic_GAM)^2)/sum((basic_GAM$y - mean(basic_GAM$y))^2)
1-sum(residuals(spatiotemporal_GAM)^2)/sum((spatiotemporal_GAM$y - mean(spatiotemporal_GAM$y))^2)
classic_gwr$diagnostic$RSquare
classic_gtwr$diagnostic$RSquare

# deviance explained for the two GAM models
summary(basic_GAM)$dev.expl
summary(spatiotemporal_GAM)$dev.expl

# RMSE
sqrt(mean(residuals(ols)^2))
sqrt(mean(residuals(basic_GAM)^2))
sqrt(mean(residuals(spatiotemporal_GAM)^2))
sqrt(classic_gwr$diagnostic$RSS / nrow(sp_df))
sqrt(classic_gtwr$diagnostic$RSS / nrow(sp_df))

# wMAPE
sum(abs(residuals(ols))) / sum(ols$model[[1]]) # A hack to get the y values without reloading the dataset
sum(abs(residuals(basic_GAM))) / sum(basic_GAM$y)
sum(abs(residuals(spatiotemporal_GAM))) / sum(spatiotemporal_GAM$y)
sum(abs(classic_gwr$diagnostic$RSS)) / sum(abs(sp_df$Sum_resistants))
sum(abs(classic_gtwr$diagnostic$RSS)) / sum(abs(sp_df$Sum_resistants))

# MAE
mean(abs(residuals(ols)))
mean(abs(residuals(basic_GAM)))
mean(abs(residuals(spatiotemporal_GAM)))
betas_gwr <- as.matrix(sf::st_drop_geometry(classic_gwr$SDF[classic_gwr$indep_vars]))
yhat_gwr  <- rowSums(classic_gwr$args$x * betas_gwr)
mean(abs(yhat_gwr - classic_gwr$args$y))
betas_gtwr <- as.matrix(sf::st_drop_geometry(classic_gtwr$SDF[classic_gtwr$indep_vars]))
yhat_gtwr  <- rowSums(classic_gtwr$args$x * betas_gtwr)
mean(abs(yhat_gtwr - classic_gtwr$args$y))

# AICc
AIC(ols) + 2*length(coef(ols))*(length(coef(ols)) + 1) / (nrow(sp_df) - length(coef(ols)) - 1)
AIC(basic_GAM) + 2*length(coef(basic_GAM))*(length(coef(basic_GAM)) + 1) / (nrow(sp_df) - length(coef(basic_GAM)) - 1)
AIC(spatiotemporal_GAM) + 2*length(coef(spatiotemporal_GAM))*(length(coef(spatiotemporal_GAM)) + 1) / (nrow(sp_df) - length(coef(spatiotemporal_GAM)) - 1)
classic_gwr$diagnostic$AICc
classic_gtwr$diagnostic$AICc

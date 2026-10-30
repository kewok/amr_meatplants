library(mgcv)
library(dplyr)
library(tidyr)

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

# OLS regression with no account for spatiotemporal processes
ols_model <- lm(Sum_resistants ~ factor(Size) + factor(Product), data = df)
saveRDS(ols_model, file = "OLS_model.rds")

# Poisson regression with no account for spatiotemporal processes
poisson_model <- gam(Sum_resistants ~ factor(Size) + factor(Product), data = df,  family = poisson(link = "log"))
saveRDS(poisson_model, file = "basic_GAM_Poisson_model.rds")

# GAM poisson regression with spatiotemporal smoothing and interactions
gam_model <- gam(Sum_resistants ~ s(longitude, latitude, bs='tp') + s(season) + ti(longitude, latitude, season , d = c(2, 1)) + factor(Size) + factor(Product), data = df, family = poisson(link = "log"), method='REML')
saveRDS(gam_model, file = "GAM_model.rds")

library(GWmodel3)
library(sp)
spdf <- df
coordinates(spdf) <- ~longitude + latitude
provoliko = "ESRI:54009"
proj4string(spdf) <- provoliko
sp_df <- sf::st_as_sf(spdf)

# Classical GWR, using the default Gaussian kernel and an adaptive bandwidth
classic_gwr <- gwr_basic(Sum_resistants ~ factor(Size) + factor(Product),data=sp_df, parallel_method='omp',parallel_arg = 8)
saveRDS(classic_gwr, file = "basicGWR.rds")

# Classical GTWR, using the default Gaussian kernel and an adaptive bandwidth
classic_gtwr <- gtwr(Sum_resistants ~ factor(Size) + factor(Product), timestamps=as.numeric(sp_df$season), data=sp_df, parallel_method='omp',parallel_arg = 8)
saveRDS(classic_gtwr, file = "basicGTWR.rds")
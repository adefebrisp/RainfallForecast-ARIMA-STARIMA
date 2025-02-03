# Set the working directory 
# Replace the path below with the directory where your shapefile data is stored
setwd("C:/Users/LENOVO/Downloads/Code_and_Data_for_FORECASTING_RAINFALL_TO_ANTICIPATE_CROP_FAILURE_USING_ARIMA_AND_STARIMA__A_CASE_STUDY_OF_JAWA_TIMUR__INDONE_442840142/Data")

# Install necessary R packages 
install.packages("sp")        # For spatial data structures and methods
install.packages("rgdal")     # For reading and writing spatial data formats
install.packages("raster")    # For working with raster data
install.packages("readGDAL")  #For reading GDAL data formats
install.packages("tmap")      # For thematic mapping
install.packages("ggplot2")   # For data visualization
install.packages("plyr")      # For data manipulation
install.packages("gstat")     # For geostatistical analysis
install.packages("spDataLarge") # For large spatial datasets
install.packages("knitr")     # For dynamic report generation
install.packages("forecast")  # For time series forecasting

# Import the 'starima' package from the specified file location
source("path_file")

# Load the required libraries for analysis
library(lattice)      # For lattice plots
library(spdep)       # For spatial dependence analysis
library(sf)          # For spatial data handling
library(sftime)      # For handling time in sf objects
library(tmap)        # For thematic mapping
library(ggplot2)     # For data visualization
library(gridExtra)   # For arranging plots
library(gstat)       # For geostatistical analysis
library(OpenStreetMap)  # For accessing OpenStreetMap data
library(spacetime)   # For space-time analysis
library(stars)       # For handling raster data
library(reshape)     # For reshaping data
library(knitr)       # For dynamic report generation
library(forecast)    # For time series forecasting
library(urca)        # For unit root and cointegration analysis

# Read the spatial dataset for Jawa Timur (East Java)
jatim <- st_read(dsn = "Jatim.shp", layer = "Jatim")

#--------------------------------------------------------------------------
# Purpose: Generate column names for a spatial dataset spanning from 
#          2002 to 2023 with month numbers prefixed, and replace existing 
#          column names in the dataset accordingly.

# Generate column names for the range 2002 to 2023 with month numbers prefixed
start_year <- 2002
end_year <- 2023
num_years <- end_year - start_year + 1

# Initialize an empty vector to store column names
column_names <- c()

# Generate column names with month numbers prefixed
for (year in start_year:end_year) {
  for (month in 1:12) {
    column_names <- c(column_names, sprintf("%02d_%d", month, year))
  }
}

# Determine the number of columns to replace
num_cols_to_replace <- min(length(column_names), ncol(jatim) - 3)

# Assign the correct column names starting from the fourth column
colnames(jatim)[4:(3 + num_cols_to_replace)] <- column_names[1:num_cols_to_replace]

#-------------------------------------------------------------------------
# Create a matrix from the data without geometry information
jatim_rainfall_matrix <- data.matrix(st_drop_geometry(jatim[, -c(1:3)]))

# Set row names of the matrix to administrative division names
rownames(jatim_rainfall_matrix) <- st_drop_geometry(jatim[, "NOZONA_LAM"])[[1]]

#=========================================================================
#Exploratory Data Analysis
#=========================================================================
# JATIM_18
# Calculate summary statistic 
summary(jatim_rainfall_matrix[1, ])
# Calculate the mean rainfall
mu_18 <- mean(jatim_rainfall_matrix[1,])
# Calculate the standard deviation of rainfall
sdev_18 <- sd(jatim_rainfall_matrix[1,])

# Explore the distribution of the data using a histogram
hist(jatim_rainfall_matrix[1,], main = "Histogram of rainfall data in JATIM_18 from 2002 to 2023", xlab = "Rainfall (mm)")
# Optionally, add a vertical line at the mean (mu) on the histogram
abline(v = mu_18, col = "red")

# Check the deviation of the distribution from normality using a QQ plot
qqnorm(jatim_rainfall_matrix[1,])
qqline(jatim_rainfall_matrix[1,], col = "red")

# JATIM_21---------------------------------------------------------------------
# Calculate summary statistic 
summary(jatim_rainfall_matrix[2, ])
# Calculate the mean rainfall
mu_21 <- mean(jatim_rainfall_matrix[2,])
# Calculate the standard deviation of rainfall
sdev_21 <- sd(jatim_rainfall_matrix[2,])

# Explore the distribution of the data using a histogram
hist(jatim_rainfall_matrix[2,], main = "Histogram of rainfall data in JATIM_21 from 2002 to 2023", xlab = "Rainfall (mm)")
# Optionally, add a vertical line at the mean (mu) on the histogram
abline(v = mu_21, col = "red")

# Check the deviation of the distribution from normality using a QQ plot
qqnorm(jatim_rainfall_matrix[2,])
qqline(jatim_rainfall_matrix[2,], col = "red")

# JATIM_22----------------------------------------------------------------------
# Calculate summary statistic 
summary(jatim_rainfall_matrix[3, ])
# Calculate the mean rainfall
mu_22 <- mean(jatim_rainfall_matrix[3,])
# Calculate the standard deviation of rainfall
sdev_22 <- sd(jatim_rainfall_matrix[3,])

# Explore the distribution of the data using a histogram
hist(jatim_rainfall_matrix[3,], main = "Histogram of rainfall data in JATIM_22 from 2002 to 2023", xlab = "Rainfall (mm)")
# Optionally, add a vertical line at the mean (mu) on the histogram
abline(v = mu_22, col = "red")

# Check the deviation of the distribution from normality using a QQ plot
qqnorm(jatim_rainfall_matrix[3,])
qqline(jatim_rainfall_matrix[3,], col = "red")

#-------------------------------------------------------------------------
# Create time series plots for selected seasonal zone in Jawa Timur 
# JATIM_18
ts_data <- ts(jatim_rainfall_matrix[1,])
p1 <- autoplot(ts_data) +
  labs(x = "Time (in months)", y = "Average rainfall (mm)") +
  ggtitle("JATIM_18")

# JATIM_21
ts_data <- ts(jatim_rainfall_matrix[2,])
p2 <- autoplot(ts_data) +
  labs(x = "Time (in months)", y = "Average rainfall (mm)") +
  ggtitle("JATIM_21")

# JATIM_22
ts_data <- ts(jatim_rainfall_matrix[3,])
p3 <- autoplot(ts_data) +
  labs(x = "Time (in months)", y = "Average rainfall (mm)") +
  ggtitle("JATIM_22")

# Arrange the plots in a grid layout
grid.arrange(p1, p2, p3)

#-------------------------------------------------------------------------
# Temporal Autocorrelation Analysis
# Analyze and visualize lagged relationships in the mean rainfall data 

# Calculate the column means for the rainfall data starting from column 2
MeanRainfallJatim <- colMeans(jatim_rainfall_matrix[, 2:ncol(jatim_rainfall_matrix)])

# Extract the month_year information from the column names
month_year <- gsub("\\..*", "", colnames(jatim_rainfall_matrix)[-1])

# Create a lagged dataset based on the mean values
Lagged <- data.frame(month_year = month_year[-1],  # Remove the first element as there's no lag for it
                     t = MeanRainfallJatim[2:length(MeanRainfallJatim)], 
                     t_minus_1 = MeanRainfallJatim[1:(length(MeanRainfallJatim)-1)])

# Convert month_year to a factor with appropriate levels
Lagged$month_year <- factor(Lagged$month_year, levels = unique(Lagged$month_year))

# Plotting
p1 <- ggplot(Lagged, aes(x = month_year, y = t, group = 1)) + 
  geom_line() +
  labs(x = "Month_Year", y = "t") +
  scale_x_discrete(breaks = Lagged$month_year[c(1, which(as.integer(sub("_.*", "", Lagged$month_year)) %% 12 == 0), length(Lagged$month_year))],
                   labels = Lagged$month_year[c(1, which(as.integer(sub("_.*", "", Lagged$month_year)) %% 12 == 0), length(Lagged$month_year))]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(Lagged, aes(x = t, y = t_minus_1)) + 
  geom_point() + 
  labs(x = "t", y = "t-1") +
  geom_smooth(method = "lm") +
  annotate("text", x = 100, y = 400, label = paste("r =", round(cor(Lagged$t, Lagged$t_minus_1), 3))) 

# Combine plots side by side
grid.arrange(p1, p2, nrow = 1)

#-------------------------------------------------------------------------
# Spatial Autocorrelation Analysis
# Convert the polygon neighborhood list to a spatial weights matrix
W <- nb2listw(poly2nb(jatim))

# Display the spatial weights matrix
W

# Convert the spatial weights matrix to a matrix format
kable(listw2mat(W))

# Calculate global Moran's I correlation and Monte Carlo Simulation
jatim_rainfall_avg <- rowMeans(jatim_rainfall_matrix)
moran.test(x = jatim_rainfall_avg, listw = W)
moran.mc(x = jatim_rainfall_avg, listw = W, nsim = 20)

# Calculate local Moran's I statistics
lm <- localmoran(x = rowMeans(jatim_rainfall_matrix), listw = W)
lm

#-------------------------------------------------------------------------
# Test the stationarity of time series data for selected seasonal zone in Jatim
# Using Augmented Dickey-Fuller (ADF) and Kwiatkowski-Phillips-Schmidt-Shin (KPSS) tests

# Time series data for JATIM_18
ts_jatim_18 <- ts(jatim_rainfall_matrix[1,])
summary(ur.df(ts_jatim_18))  # ADF
summary(ur.kpss(ts_jatim_18))  # KPSS

# Time series data for JATIM_21
ts_jatim_21 <- ts(jatim_rainfall_matrix[2,])
summary(ur.df(ts_jatim_21))  # ADF
summary(ur.kpss(ts_jatim_21))  # KPSS

# Time series data for JATIM_22
ts_jatim_22 <- ts(jatim_rainfall_matrix[3,])
summary(ur.df(ts_jatim_22))  # ADF
summary(ur.kpss(ts_jatim_22))  # KPSS

#=========================================================================
#ARIMA
#=========================================================================
#JATIM_18
# Purpose: Analyze and model the time series data for rainfall in JATIM_18
# Training set -> jatim_rainfall_matrix [1:258] -> 01_2002 - 06_2023
# Test set -> jatim_rainfall_matrix [259:264] -> 07_2023 - 12_2023

# Plot the autocorrelation function (ACF) 
acf(jatim_rainfall_matrix[1,], lag.max=50, main="ACF JATIM_18")

# Plot the partial autocorrelation function (PACF)
pacf(jatim_rainfall_matrix[1,], lag.max=50, main="PACF ACF JATIM_18")

# AutoARIMA model fitting
fit.auto.ar.18 <- auto.arima(jatim_rainfall_matrix["JATIM_18",1:258])
fit.auto.ar.18

# Model fitting with different orders
fit.ar.18.1 <- Arima(jatim_rainfall_matrix["JATIM_18", 1:258], order = c(1, 0, 1))
fit.ar.18.2 <- Arima(jatim_rainfall_matrix["JATIM_18", 1:258], order = c(1, 0, 2))

# Display model summaries
fit.ar.18.1
fit.ar.18.2
fit.auto.ar.18

# Calculate Normalized Root Mean Squared Error (NRMSE) for each model
NRMSE_fit.ar.18.1 <- NRMSE(res=fit.ar.18.1$residuals, obs=jatim_rainfall_matrix["JATIM_18",1:258])
NRMSE_fit.ar.18.2 <- NRMSE(res=fit.ar.18.2$residuals, obs=jatim_rainfall_matrix["JATIM_18",1:258])
NRMSE_fit.auto.ar.18 <- NRMSE(res=fit.auto.ar.18$residuals, obs=jatim_rainfall_matrix["JATIM_18",1:258])

# Display NRMSE values for each model
NRMSE_fit.ar.18.1
NRMSE_fit.ar.18.2
NRMSE_fit.auto.ar.18

# Purpose: Check the normality of residuals for diagnostic purposes.
# Time series diagnostics for the best ARIMA model -> ARIMA(1,0,2)
tsdiag(fit.ar.18.2)

# Box-Ljung test for residual autocorrelation
# Setting the lag parameter to twice the order of the autoregressive part of the model
# This allows the test to capture any residual autocorrelation in the data beyond what is captured by the ARIMA model
Box.test(fit.ar.18.2$residuals, lag = 3)

# Purpose: Perform time series forecasting for rainfall in JATIM_18
# Fit the best ARIMA model to the test set of rainfall data
pre.Ar.18 <- Arima(jatim_rainfall_matrix["JATIM_18", 259:(ncol(jatim_rainfall_matrix))], model=fit.ar.18.2)

# Plot the observed and fitted values from the ARIMA model
matplot(cbind(pre.Ar.18$fitted, pre.Ar.18$x), type="l", main="Prediction Rainfall of JATIM_18")

# Calculate the Normalized Root Mean Squared Error (NRMSE) for the forecasted values
NRMSE_fit.18 <- NRMSE(res=pre.Ar.18$residuals, obs=jatim_rainfall_matrix["JATIM_18", 259:264])
NRMSE_fit.18

#JATIM_21------------------------------------------------------------
# Purpose: Analyze and model the time series data for rainfall in JATIM_21
# Training set -> jatim_rainfall_matrix [1:258] -> 01_2002 - 06_2023
# Test set -> jatim_rainfall_matrix [259:264] -> 07_2023 - 12_2023

# Plot the autocorrelation function (ACF) 
acf(jatim_rainfall_matrix[2,], lag.max=50, main="ACF JATIM_21")

# Plot the partial autocorrelation function (PACF)
pacf(jatim_rainfall_matrix[2,], lag.max=50, main="PACF ACF JATIM_21")

# AutoARIMA model fitting
fit.auto.ar.21 <- auto.arima(jatim_rainfall_matrix["JATIM_21",1:258])
fit.auto.ar.21

# Model fitting with different orders
fit.ar.21.1 <- Arima(jatim_rainfall_matrix["JATIM_21", 1:258], order = c(1, 0, 2))
fit.ar.21.2 <- Arima(jatim_rainfall_matrix["JATIM_21", 1:258], order = c(2, 0, 2))
fit.ar.21.3 <- Arima(jatim_rainfall_matrix["JATIM_21", 1:258], order = c(5, 0, 2))
fit.ar.21.4 <- Arima(jatim_rainfall_matrix["JATIM_21", 1:258], order = c(6, 0, 2))
fit.ar.21.5 <- Arima(jatim_rainfall_matrix["JATIM_21", 1:258], order = c(1, 0, 1))
fit.ar.21.6 <- Arima(jatim_rainfall_matrix["JATIM_21", 1:258], order = c(2, 0, 1))
fit.ar.21.7 <- Arima(jatim_rainfall_matrix["JATIM_21", 1:258], order = c(3, 0, 1))
fit.ar.21.8 <- Arima(jatim_rainfall_matrix["JATIM_21", 1:258], order = c(4, 0, 1))
fit.ar.21.9 <- Arima(jatim_rainfall_matrix["JATIM_21", 1:258], order = c(6, 0, 1))

# Display model summaries
fit.ar.21.1
fit.ar.21.2
fit.ar.21.3
fit.ar.21.4
fit.ar.21.5
fit.ar.21.6
fit.ar.21.7
fit.ar.21.8
fit.ar.21.9
fit.auto.ar.21

# Calculate Normalized Root Mean Squared Error (NRMSE) for each model
NRMSE_fit.ar.21.1 <- NRMSE(res=fit.ar.21.1$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])
NRMSE_fit.ar.21.2 <- NRMSE(res=fit.ar.21.2$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])
NRMSE_fit.ar.21.3 <- NRMSE(res=fit.ar.21.3$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])
NRMSE_fit.ar.21.4 <- NRMSE(res=fit.ar.21.4$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])
NRMSE_fit.ar.21.5 <- NRMSE(res=fit.ar.21.5$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])
NRMSE_fit.ar.21.6 <- NRMSE(res=fit.ar.21.6$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])
NRMSE_fit.ar.21.7 <- NRMSE(res=fit.ar.21.7$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])
NRMSE_fit.ar.21.8 <- NRMSE(res=fit.ar.21.8$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])
NRMSE_fit.ar.21.9 <- NRMSE(res=fit.ar.21.9$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])
NRMSE_fit.auto.ar.21 <- NRMSE(res=fit.auto.ar.21$residuals, obs=jatim_rainfall_matrix["JATIM_21",1:258])

# Display NRMSE values for each model
NRMSE_fit.ar.21.1
NRMSE_fit.ar.21.2
NRMSE_fit.ar.21.3
NRMSE_fit.ar.21.4
NRMSE_fit.ar.21.5
NRMSE_fit.ar.21.6
NRMSE_fit.ar.21.7
NRMSE_fit.ar.21.8
NRMSE_fit.ar.21.9
NRMSE_fit.auto.ar.21

# Purpose: Check the normality of residuals for diagnostic purposes.
# Time series diagnostics for the best ARIMA model -> ARIMA(5,0,2)
tsdiag(fit.ar.21.3)

# Box-Ljung test for residual autocorrelation
# Setting the lag parameter to twice the order of the autoregressive part of the model
# This allows the test to capture any residual autocorrelation in the data beyond what is captured by the ARIMA model
Box.test(fit.ar.21.3$residuals, lag = 10)

# Purpose: Perform time series forecasting for rainfall in JATIM_21
# Fit the best ARIMA model to the test set of rainfall data
pre.Ar.21 <- Arima(jatim_rainfall_matrix["JATIM_21", 259:(ncol(jatim_rainfall_matrix))], model=fit.ar.21.3)

# Plot the observed and fitted values from the ARIMA model
matplot(cbind(pre.Ar.21$fitted, pre.Ar.21$x), type="l", main="Prediction Rainfall of JATIM_21")

# Calculate the Normalized Root Mean Squared Error (NRMSE) for the forecasted values
NRMSE_fit.21 <- NRMSE(res=pre.Ar.21$residuals, obs=jatim_rainfall_matrix["JATIM_21", 259:264])
NRMSE_fit.21

#JATIM_22------------------------------------------------------------
# Purpose: Analyze and model the time series data for rainfall in JATIM_22
# Training set -> jatim_rainfall_matrix [1:258] -> 01_2002 - 06_2023
# Test set -> jatim_rainfall_matrix [259:264] -> 07_2023 - 12_2023

# Plot the autocorrelation function (ACF) 
acf(jatim_rainfall_matrix[3,], lag.max=50, main="ACF JATIM_22")

# Plot the partial autocorrelation function (PACF)
pacf(jatim_rainfall_matrix[3,], lag.max=50, main="PACF ACF JATIM_22")

# AutoARIMA model fitting
fit.auto.ar.22 <- auto.arima(jatim_rainfall_matrix["JATIM_22",1:258])
fit.auto.ar.22

# Model fitting with different orders
fit.ar.22.1 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(1, 0, 2))
fit.ar.22.2 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(2, 0, 2))
fit.ar.22.3 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(4, 0, 2))
fit.ar.22.4 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(5, 0, 2))
fit.ar.22.5 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(6, 0, 2))
fit.ar.22.6 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(1, 0, 1))
fit.ar.22.7 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(2, 0, 1))
fit.ar.22.8 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(3, 0, 1))
fit.ar.22.9 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(4, 0, 1))
fit.ar.22.10 <- Arima(jatim_rainfall_matrix["JATIM_22", 1:258], order = c(6, 0, 1))

# Display model summaries
fit.ar.22.1
fit.ar.22.2
fit.ar.22.3
fit.ar.22.4
fit.ar.22.5
fit.ar.22.6
fit.ar.22.7
fit.ar.22.8
fit.ar.22.9
fit.ar.22.10
fit.auto.ar.22

# Calculate Normalized Root Mean Squared Error (NRMSE) for each model
NRMSE_fit.ar.22.1 <- NRMSE(res=fit.ar.22.1$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.ar.22.2 <- NRMSE(res=fit.ar.22.2$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.ar.22.3 <- NRMSE(res=fit.ar.22.3$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.ar.22.4 <- NRMSE(res=fit.ar.22.4$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.ar.22.5 <- NRMSE(res=fit.ar.22.5$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.ar.22.6 <- NRMSE(res=fit.ar.22.6$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.ar.22.7 <- NRMSE(res=fit.ar.22.7$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.ar.22.8 <- NRMSE(res=fit.ar.22.8$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.ar.22.9 <- NRMSE(res=fit.ar.22.9$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.ar.22.10 <- NRMSE(res=fit.ar.22.10$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])
NRMSE_fit.auto.ar.22 <- NRMSE(res=fit.auto.ar.22$residuals, obs=jatim_rainfall_matrix["JATIM_22",1:258])

# Display NRMSE values for each model
NRMSE_fit.ar.22.1
NRMSE_fit.ar.22.2
NRMSE_fit.ar.22.3
NRMSE_fit.ar.22.4
NRMSE_fit.ar.22.5
NRMSE_fit.ar.22.6
NRMSE_fit.ar.22.7
NRMSE_fit.ar.22.8
NRMSE_fit.ar.22.9
NRMSE_fit.ar.22.10
NRMSE_fit.auto.ar.22

# Purpose: Check the normality of residuals for diagnostic purposes.
# Time series diagnostics for the best ARIMA model -> ARIMA(4,0,2)
tsdiag(fit.ar.22.3)

# Box-Ljung test for residual autocorrelation
# Setting the lag parameter to twice the order of the autoregressive part of the model
# This allows the test to capture any residual autocorrelation in the data beyond what is captured by the ARIMA model
Box.test(fit.ar.22.3$residuals, lag = 8)

# Purpose: Perform time series forecasting for rainfall in JATIM_22
# Fit the best ARIMA model to the test set of rainfall data
pre.Ar.22 <- Arima(jatim_rainfall_matrix["JATIM_22", 259:(ncol(jatim_rainfall_matrix))], model=fit.ar.22.3)

# Plot the observed and fitted values from the ARIMA model
matplot(cbind(pre.Ar.22$fitted, pre.Ar.22$x), type="l", main="Prediction Rainfall of JATIM_22")

# Calculate the Normalized Root Mean Squared Error (NRMSE) for the forecasted values
NRMSE_fit.22 <- NRMSE(res=pre.Ar.22$residuals, obs=jatim_rainfall_matrix["JATIM_22", 259:264])
NRMSE_fit.22

#=========================================================================
#STARIMA
#=========================================================================
# Purpose: Convert polygon neighbors to a binary spatial weights matrix.
# Create neighborhood structure from polygon neighbors
nbrs <- poly2nb(jatim)

# Convert the neighborhood structure to a binary spatial weights matrix
W1 <- nb2mat(nbrs)

# Create an identity matrix with dimensions matching the spatial weights matrix
W0 <- diag(x=1, nrow(W1), ncol(W1)) # Spatial order zero

# Purpose: Calculate spatial-temporal autocorrelation functions (STACF and STPACF) for rainfall data in Jawa Timur.
# Transpose the rainfall matrix to have time series in columns (STACF and STPACF require this format)
jatim_rainfall_mat2 <- t(jatim_rainfall_matrix)

#Model Identification--------------------------------------------------------
# Calculate the spatial-temporal autocorrelation function (STACF) for the rainfall data using a spatial weights matrix
stacf(jatim_rainfall_mat2, W0, 50) # Spatial order zero

# Calculate the spatial-temporal partial autocorrelation function (STPACF) for the rainfall data using a spatial weights matrix
stpacf(jatim_rainfall_mat2, W1, 50) # Spatial order one

#Parameter Estimation and Fitting--------------------------------------------
# Create a list of spatial weight matrices, where zero elements are not needed
W_fit <- list(w1 = W1)

# Fit space-time autoregressive integrated moving average (STARIMA) models with different orders
fit.star1 <- starima_fit(Z = jatim_rainfall_mat2[1:258, ], W = W_fit, p = 1, d = 0, q = 1)
fit.star2 <- starima_fit(Z = jatim_rainfall_mat2[1:258, ], W = W_fit, p = 1, d = 0, q = 2)

fit.star1
fit.star2

# Extract NRMSE values from STARIMA model fits
# NRMSE: Normalized Root Mean Squared Error, indicating the goodness-of-fit of the models
fit.star1$NRMSE
fit.star2$NRMSE # The best ST-ARIMA(1,0,2) model (The lowest NRMSE)

# Diagnostic Checking----------------------------------------------------------
stacf(fit.star2$RES,W1,50) 

# Histogram of residual
hist(fit.star2$RES[,1])
hist(fit.star2$RES[,2])
hist(fit.star2$RES[,3])

# Perform the Ljung-Box test on the residuals of the STARIMA model
#JATIM_18 -> STARIMA(1,0,2) -> Best model
Box.test(fit.star2$RES[,1],lag=1, type="Ljung")
#JATIM_21 -> STARIMA(1,0,2) -> Best model
Box.test(fit.star2$RES[,2],lag=1, type="Ljung")
#JATIM_22 -> STARIMA(1,0,2) -> Best model
Box.test(fit.star2$RES[,3],lag=1, type="Ljung")

#PREDICTION--------------------------------------------------------------------
pre.star <- starima_pre(jatim_rainfall_mat2[(258-0-3+1):264,],
                        model=fit.star2)
matplot(1:6,cbind(jatim_rainfall_mat2[259:264,1],pre.star$PRE[,1]),type="l", main="Prediction Rainfall of JATIM_18 (ST-ARIMA)")
matplot(1:6,cbind(jatim_rainfall_mat2[259:264,2],pre.star$PRE[,2]),type="l", main="Prediction Rainfall of JATIM_21 (ST-ARIMA)")
matplot(1:6,cbind(jatim_rainfall_mat2[259:264,3],pre.star$PRE[,3]),type="l", main="Prediction Rainfall of JATIM_22 (ST-ARIMA)")

#Check the NRMSE of the ST-ARIMA Prediction
pre.star$NRMSE

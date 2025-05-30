######################################################
################# Linear Time Series #################
###### BENABDESADOK Nayel, CIANFARANI Ana-Sofia #######
######################################################

####### - Part 0 : Installations & Data Download - ######

# Load required packages
library(fUnitRoots)
library(zoo)
library(tseries)
library(forecast)
library(ggplot2)
library(scales)      
library(RColorBrewer)
library(gridExtra)
library(astsa)
library(tsoutliers)
library(ellipse)

# # Load the data on personal computer (Sofia)
# path <- "C:/Users/cianf/OneDrive/Documents/STL"
# path <- "C:/Users/cianf/OneDrive/Documents/STL/data.csv"
# # Load the data on personnal computer (Nayel)
path_donnees <- "/Users/nayelbenabdesadok/GitProjects/Time_Series_ENSAE/reprise/data.csv"
path <- "/Users/nayelbenabdesadok/GitProjects/Time_Series_ENSAE/reprise"
setwd(path)
getwd()


# Load the data 
data <- read.csv(file = path_donnees, sep = ";")

# Clean the table
data_all <- data[4:425, ]
data_all$Codes <- NULL
colnames(data_all) <- c("Date_MY", "IPI")


# Converting date to the correct format
data_all$Date <- as.Date(paste(data_all$Date_MY, 1, sep = "-"), format = "%Y-%m-%d")
data_all$Date_MY <- NULL

# Extracting month and year from date
Year <- as.numeric(format(data_all$Date, format = "%Y"))
Month <- format(data_all$Date, format = "%m")
data_all <- cbind(data_all, Year, Month)
data_all$Year <- as.numeric(data_all$Year)
data_all$IPI <- as.numeric(data_all$IPI)
sort(data_all$IPI)
data_all <- data_all[order(data_all$Date), ]
rownames(data_all) <- seq(length = nrow(data_all))


# Create a time series for easier analysis, and plot
Xt.ts <- ts(data_all$IPI, start = c(1990, 1), end = c(2025, 2), frequency = 12)
par(mfrow = c(1, 1))
plot.ts(Xt.ts, xlab = "Years", ylab = "IPI Construction (CVS-CJO, base 2021)")

png("Xt.png", width = 1000, height = 600)
plot.ts(Xt.ts,
        main = "Main Series",
        ylab = "IPI Construction (CVS-CJO, base 2021)",
        xlab = "Time",
        lwd = 2)
dev.off()
#==============================================================================#
# The dates range from January 1990 (1990m1) to February 2025 (2025m2) as shown above.
# The series represents the Index of Industrial Production (IPI) in the construction sector, 
# seasonally and working-day adjusted (CVS-CJO), with a base of 100 in 2021.
# At first glance, the series appears non-stationary, with multiple structural breaks.
# A sharp drop is clearly visible around April 2020, likely linked to the COVID-19 crisis.
# Before applying any transformations, we will examine the overall trend and variability.
#==============================================================================#



####### - Part I : The Data - ######

####### - Part 1.1 : Trend/Cycle/Seasonal Decomposition (What does the chosen series represent ?)- ######
# Decomposition of the time series to analyze the trend component.
decompose <- decompose(Xt.ts,type="multiplicative")
png("decompose.png", width = 1000, height = 600)
autoplot(decompose, main = "Time Series Decomposition", xlab ="Year")
dev.off()


#==============================================================================#
# Trend Analysis: the series displays three distinct phases –
# a decline from 1990 to 2000, a growth period from 2000 to 2008,
# and a significant downward trend from 2008 to 2025.
# A sharp outlier is observed in early 2020, likely due to the COVID-19 crisis.
#==============================================================================#

#===========================#
#  Is there Seasonality ?   #
#===========================#
n_years <- length(unique(floor(time(Xt.ts))))
colors <- scales::hue_pal()(n_years)

# Standard season plot (monthly profile)
png("season_plot.png", width = 1000, height = 600)
ggseasonplot(Xt.ts,
             year.labels = TRUE,
             year.labels.left = TRUE) +
  ylab("Production Index - Construction") +
  xlab("Month") +
  ggtitle("Monthly Profiles (Season Plot)") +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14)
dev.off()

# Polar season plot
png("season_polar_plot.png", width = 1000, height = 600)
ggseasonplot(Xt.ts,
             polar = TRUE) +
  ylab("Production Index - Construction") +
  xlab("Month") +
  ggtitle("Seasonality in Polar Coordinates") +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14)
dev.off()

# Monthly boxplots
png("monthly_boxplots.png", width = 1000, height = 600)
cols <- RColorBrewer::brewer.pal(12, "Set3")
boxplot(Xt.ts ~ cycle(Xt.ts),
        col = cols,
        pch = 20,
        cex = 0.5,
        main = "Monthly Boxplots",
        ylab = "Production Index - Construction",
        xlab = "Month",
        names = month.abb)
dev.off()
#==============================================================================#
# Regardless of the plot displayed, we observe no influence of the month
# on the behavior of the series.
# Absence of seasonality.
#==============================================================================#



####### - Part 1.2 : Stationary ######
#===========================#
#  Is Xt.ts Stationary ?    #
#===========================#
# ACF
acf_plot_orig <- ggAcf(Xt.ts, lag.max = 100, plot = TRUE) +
  ggtitle("ACF of Original Series")
#PACF
pacf_plot_orig <- ggPacf(Xt.ts, lag.max = 100, plot = TRUE) +
  ggtitle("PACF of Original Series")

png("ACF_PACF_Xt.png", width = 1000, height = 600)
gridExtra::grid.arrange(acf_plot_orig, pacf_plot_orig, ncol = 2)
dev.off()

#==============================================================================#
# ACF decays slowly + PACF behaves like white noise.
# Suggests non-stationarity → Let's run some tests: ADF, PP, KPSS.
#==============================================================================#

# Center the series
Xt_centered <- Xt.ts - mean(Xt.ts)
ggAcf(Xt_centered,lag.max=100, plot=T)

#==============================================================================#
# We decided to center the series but this is optional : 
# the series X_t and X_t_centered are the same in our analysis ! 
#==============================================================================#


#==============================================================================#
# To determine the lag to use for our ADF test, we refer to Schwert (1989).
# The code below implements their lag selection algorithm. 
# For more details see our report, Appendix : Stationarity tests.
# We set type = ct as there is a trend.
#==============================================================================#

# Initialize the lag value at pmax = 17
lag <- 17

# Store the ADF test result
adf_result <- NULL

# Loop until the t-statistic is above 1.6 in absolute value
while (lag > 0) {
  # Perform the ADF test
  adf_result <- adfTest(Xt_centered, lags = lag, type = "ct")
  
  # Check the t-statistic
  if (abs(adf_result@test$statistic) > 1.6) {
    break
  }
  
  # Decrease the lag value
  lag <- lag - 1
}
cat("Final lag value:", lag, "\n")
print(adf_result)

# Phillips-Perron Test (PP)
library(tseries)
cat("\nPhillips-Perron Test Results:\n")
print(pp.test(Xt_centered))

# KPSS Test
cat("\nKPSS Test Results:\n")
kpss_result <- kpss.test(Xt_centered)
print(kpss_result)

#==============================================================================#
# Stationarity Tests Summary:
# - ADF Test (lag = 8, type = "ct") → p-value = 0.7114 → FAIL to reject H0
# - PP Test → p-value < 0.01 → REJECT H0 (suggests stationarity but less reliable for finite samples)
# - KPSS Test → p-value = 0.01 → REJECT H0 of stationarity
#
# Conclusion: Two out of three tests (ADF and KPSS) indicate non-stationarity.
# Although the PP test suggests stationarity, the overall evidence leans toward
# non-stationarity of the centered series, likely due to trend or structural changes.
#==============================================================================#

#==============================================================================#
# Stationarity Tests on the Differenced Centered Series
#==============================================================================#

diff_Xt_centered <- diff(Xt_centered)

png("diff_Xt_centered.png", width = 1000, height = 600)
plot.ts(diff_Xt_centered,
        main = "Differenced Centered Series",
        ylab = "Differenced IPI",
        xlab = "Time",
        lwd = 2)
dev.off()

#==============================================================================#
# This time, our ADF test works with the max lag immediately.
# We set type = nc as there is the series is centered and stationary.
#==============================================================================#

# Augmented Dickey-Fuller Test (ADF)
adf_result_diff <- adfTest(diff_Xt_centered, lags = 17, type = "nc")
cat("ADF Test Results (Differenced Series):\n")
print(adf_result_diff)

# Phillips-Perron Test (PP)
cat("\nPhillips-Perron Test Results (Differenced Series):\n")
print(pp.test(diff_Xt_centered))

# KPSS Test
cat("\nKPSS Test Results (Differenced Series):\n")
kpss_result_diff <- kpss.test(diff_Xt_centered)
print(kpss_result_diff)

#==============================================================================#
## Stationarity Tests Summary:
# - ADF Test (lag = 17, type = "nc") → p-value = 0.01 → REJECT H0 - Stationarity
# - PP Test → p-value < 0.01 → REJECT H0 - Stationarity
# - KPSS Test → p-value > 0.1 → DO NOT REJECT H0 - Stationarity
#
# All three stationarity tests (ADF, PP, KPSS) applied to the differenced
# centered series confirm stationarity (p-values all in favor).
# Conclusion: The original series is integrated of order 1 (I(1)).
#==============================================================================#
####### - Part 1.3 : Before/after Stationarity ######

png("Xt_vs_diffXt.png", width = 1000, height = 300)
par(mfrow = c(1, 2))  

plot.ts(Xt.ts,
        main = expression("Original Series " ~ X[t]),
        ylab = "IPI Construction",
        xlab = "Time",
        lwd = 2,
        col = "steelblue")

plot.ts(diff(Xt.ts),
        main = expression("Differenced Series " ~ Y[t] == Delta * X[t]),
        ylab = "Differenced IPI",
        xlab = "Time",
        lwd = 2,
        col = "darkorange")

dev.off()
####### - Part II : ARMA Models ######

####### - Part 2.1 : Pick p and q ######


acf_diff_plot <- ggAcf(diff_Xt_centered, lag.max = 40, plot = TRUE) +
  ggtitle("ACF of Differenced Series")
pacf_diff_plot <- ggPacf(diff_Xt_centered, lag.max = 40, plot = TRUE) +
  ggtitle("PACF of Differenced Series")

png("ACF_PACF_diff_Xt_centered.png", width = 1000, height = 300)
gridExtra::grid.arrange(acf_diff_plot, pacf_diff_plot, ncol = 2)
dev.off()

#============================================================================#
# Property of MA(2) models: if the autocorrelation function (ACF) becomes
# zero for all lags h > 2, the process can be identified as MA(2).
# Similarly, for AR(5) models: the partial autocorrelation function (PACF)
# cuts off after lag 5.
#
# p_max=5 and q_max=2. Let's test all ARMA(p,q) such that p<=pmax and q<qmax
#============================================================================#

# we call the differents models arma_p_q
#==============================================================================#
# We test all ARIMA(p,1,q) models for the construction series Xt.ts,
# where p ≤ 5 and q ≤ 2. We compute AIC and BIC for each combination.
#==============================================================================#
arma_0_1 <- sarima(Xt.ts, 0, 1, 1)
arma_0_2 <- sarima(Xt.ts, 0, 1, 2)

arma_1_0 <- sarima(Xt.ts, 1, 1, 0)
arma_1_1 <- sarima(Xt.ts, 1, 1, 1)
arma_1_2 <- sarima(Xt.ts, 1, 1, 2)

arma_2_0 <- sarima(Xt.ts, 2, 1, 0)
arma_2_1 <- sarima(Xt.ts, 2, 1, 1)
arma_2_2 <- sarima(Xt.ts, 2, 1, 2)

arma_3_0 <- sarima(Xt.ts, 3, 1, 0)
arma_3_1 <- sarima(Xt.ts, 3, 1, 1)
arma_3_2 <- sarima(Xt.ts, 3, 1, 2)

arma_4_0 <- sarima(Xt.ts, 4, 1, 0)
arma_4_1 <- sarima(Xt.ts, 4, 1, 1)
arma_4_2 <- sarima(Xt.ts, 4, 1, 2)

arma_5_0 <- sarima(Xt.ts, 5, 1, 0)
arma_5_1 <- sarima(Xt.ts, 5, 1, 1)
arma_5_2 <- sarima(Xt.ts, 5, 1, 2)


#==============================================================================#
# To eliminate unsuitable models, we first exclude models with at least 
# one non-significant coefficient (p-value > 5%), based on hypothesis testing (H0 vs H1).
#==============================================================================#


arma_0_1$ttable 
arma_0_2$ttable

arma_1_0$ttable
arma_1_1$ttable
arma_1_2$ttable # out

arma_2_0$ttable
arma_2_1$ttable # out
arma_2_2$ttable # out

arma_3_0$ttable
arma_3_1$ttable # out
arma_3_2$ttable # out

arma_4_0$ttable
arma_4_1$ttable # out
arma_4_2$ttable # out

arma_5_0$ttable
arma_5_1$ttable # out
arma_5_2$ttable # out

#==============================================================================#
# We still keep the excluded models ("out models") to compute the Ljung-Box 
# statistics on their residuals. 
# We will later exclude models whose Ljung-Box p-value is less than 5%.
#==============================================================================#
models <- c("arma_0_1", "arma_0_2",
            "arma_1_0", "arma_1_1", "arma_1_2",
            "arma_2_0", "arma_2_1", "arma_2_2",
            "arma_3_0", "arma_3_1", "arma_3_2",
            "arma_4_0", "arma_4_1", "arma_4_2",
            "arma_5_0", "arma_5_1", "arma_5_2")

for (m in models) {
  cat("\n--- Ljung-Box Test for", m, "---\n")
  print(Box.test(residuals(get(m)$fit), lag = 10, type = "Ljung-Box"))
}
# The remaining models are ARMA(1,1), ARMA(0,2), and ARMA(5,0).
# To select the final model, we choose the one that minimizes both the AIC and BIC criteria.


info_criteria <- data.frame(
  Model = character(),
  AIC = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)

for (m in models) {
  fit <- get(m)$fit
  model_aic <- AIC(fit)
  model_bic <- BIC(fit)
  
  info_criteria <- rbind(info_criteria, data.frame(
    Model = m,
    AIC = model_aic,
    BIC = model_bic
  ))
}
print(info_criteria)

# ARMA (1,1) is choosen ! 

####### - Part 2.2 : Taking outliers into account  #######

# Detect the outliers
arima111 <- arima(Xt.ts, order = c(1, 1, 1))
outlier_result <- tso(Xt.ts, types = c("AO", "LS", "TC", "IO"),
                      tsmethod = "arima", args.tsmethod = list(order = c(1, 1, 1)))
all_outliers <- outlier_result$outliers
print(all_outliers)

#==============================================================================#
# A total of 10 outliers are detected: 8 additive outliers (AO) and 2 transient changes (TC).
# The last outlier, however, does not correspond to any visually identifiable shock in the series.
# Therefore, we retain only the first 9 outliers for the adjustment.
# We then correct the series accordingly using the corresponding estimated effects.
#==============================================================================#

selected_indices <- head(all_outliers$ind, 9)
adjusted_effects <- outlier_result$effects
keep_positions <- rep(FALSE, length(adjusted_effects))
keep_positions[selected_indices] <- TRUE
adjusted_effects[!keep_positions] <- 0

yadj_9out <- Xt.ts - adjusted_effects
# Adjusted coefficients : phi1 = 0.1781 , theta1 = -0.7285


final_model_9out <- arima(yadj_9out, order = c(1, 1, 1))
summary(final_model_9out)
# AIC = 2395.14 < 2432.691  (AIC of the initial model) --> improvement confirmed

png("outliers_plot_full.png", width = 800, height = 600)
plot(outlier_result)
dev.off()



####### - Part III — Forecasting ######

####### - Part 3.1 : X_T+1 and X_T+2  #######
forecast_2 <- forecast(final_model_9out, h = 2)
print(forecast_2)

#==============================================================================#
# Point Forecast     80% CI                  95% CI
#------------------------------------------------------------------------------#
# Mar 2025:  93.86069   [88.56977 ; 99.15160]   [85.76894 ; 101.95240]
# Apr 2025:  93.86615   [88.06500 ; 99.66730]   [84.99406 ; 102.73820]
#==============================================================================#


#===========================#
# Zoom Forecast Plot (Last 36 Months)
#===========================#

start_zoom <- time(yadj_9out)[length(yadj_9out) - 35]
end_zoom <- time(forecast_2$mean)[2]

recent_values <- window(yadj_9out, start = start_zoom)
forecast_values <- forecast_2$mean
forecast_lower <- forecast_2$lower[,2]  # 95% CI
forecast_upper <- forecast_2$upper[,2]  # 95% CI

y_min <- min(recent_values, forecast_lower)
y_max <- max(recent_values, forecast_upper)

# Get the last observed point and first forecast point 
# Then create a line between them, as autoplot does not do this automatically.
last_obs_time <- time(tail(recent_values, 1))
last_obs_value <- tail(recent_values, 1)
first_fc_time <- time(forecast_values)[1]
first_fc_value <- forecast_values[1]
connect_df <- data.frame(
  Time = c(last_obs_time, first_fc_time),
  Value = c(last_obs_value, first_fc_value)
)

# Plot, including the  connecting line
png("forecast_zoom_adjusted.png", width = 1000, height = 600)
autoplot(forecast_2, series = "Forecast") +
  autolayer(recent_values, series = "Observed", color = "black") +
  geom_line(data = connect_df, aes(x = Time, y = Value), linetype = "solid", color = "black") +
  ggtitle("Two-Step Forecast: Last 36 Months (Zoom)") +
  ylab("IPI Construction (Corrected Series)") +
  xlab("Time") +
  coord_cartesian(xlim = c(start_zoom, end_zoom), ylim = c(y_min, y_max)) +
  scale_color_manual(name = "Series",
                     values = c("Observed" = "black", "Forecast" = "blue")) +
  theme_minimal(base_size = 14)
dev.off()

####### - Part 3.2 : Confidence Region  #######

phi <- coef(final_model_9out)["ar1"]
theta <- coef(final_model_9out)["ma1"]
sigma2 <- final_model_9out$sigma2

sigma_g1 <- sqrt(sigma2)
sigma_g2 <- sqrt(sigma2 * (1 + (1 + phi - theta)^2))
rho <- sigma2 * (1 + phi - theta)

Sigma <- matrix(c(sigma_g1^2, rho,
                  rho, sigma_g2^2), nrow = 2)

centre <- c(forecast_2$mean[1], forecast_2$mean[2])

ell <- ellipse(Sigma, centre = centre, level = 0.95, npoints = 1000)

png("confidence_ellipse_forecast.png", width = 800, height = 600)
plot(ell,
     type = 'l',
     xlab = expression(hat(X)[T+1]),
     ylab = expression(hat(X)[T+2]),
     main = "95% Confidence Ellipse for Two-Step Forecast")
points(x = centre[1], y = centre[2], pch = 19, col = "red", cex = 1.5)
grid()
dev.off()

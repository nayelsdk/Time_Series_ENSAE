######################################################
################# Linear Time Series #################
###### BENABDESADOK Nayel, CIANFARANI Ana-Sofia ######
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

# # Load the data on personal computer (Sofia)
# path <- "C:/Users/cianf/OneDrive/Documents/STL"
path <- "/Users/nayelbenabdesadok/GitProjects/Time_Series_ENSAE/reprise"
# # Load the data on personnal computer (Nayel)
path_donnees <- "/Users/nayelbenabdesadok/GitProjects/Time_Series_ENSAE/reprise/data.csv"
setwd(path)
getwd()
# list.files()

# path_donnees <- "C:/Users/cianf/OneDrive/Documents/STL/CVSCJO.csv"
# data <- read_delim(path_donnees,";", escape_double = FALSE)


# Load the data online
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



####### - Part 1 : The Data - ######

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

# Center the series (optional for visual clarity)
Xt_centered <- Xt.ts - mean(Xt.ts)
ggAcf(Xt_centered,lag.max=100, plot=T)

# Augmented Dickey-Fuller Test (ADF)
adf_result <- adfTest(Xt_centered, lags = 10, type = "ct")
cat("ADF Test Results:\n")
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
# - ADF Test (lag = 10, type = "ct") → p-value = 0.827 → FAIL to reject H0
# - PP Test → p-value < 0.01 → REJECT H0 (suggests stationarity)
# - KPSS Test → p-value < 0.01 → REJECT H0 of stationarity
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



# Augmented Dickey-Fuller Test (ADF)
adf_result_diff <- adfTest(diff_Xt_centered, lags = 10, type = "ct")
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
# All three stationarity tests (ADF, PP, KPSS) applied to the differenced
# centered series confirm stationarity (p-values all in favor).
# Conclusion: The original series is integrated of order 1 (I(1)).
#==============================================================================#



####### - Part 2 : ARMA Models ######

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
# We will later exclude models whose Ljung-Box p-value is less than 15%.
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

####### - Part 2.2 : Taking outliers into account  ######

outlier_result <- tso(Xt.ts, types = c("AO", "LS", "TC","IO"))
summary(outlier_result)
outlier_result$outliers
png("outliers_plot.png", width = 800, height = 600)
plot(outlier_result)
dev.off()



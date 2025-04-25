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
library(readr)

# # Load the data on personal computer (Sofia)
# path <- "C:/Users/cianf/OneDrive/Documents/STL"
# setwd(path)
# getwd() 
# list.files() 

# path_donnees <- "C:/Users/cianf/OneDrive/Documents/STL/CVSCJO.csv"
# data <- read_delim(path_donnees,";", escape_double = FALSE)


# Load the data online
path_donnees <- "/home/onyxia/work/Time_Series_ENSAE/Partie_1_clean/CVSCJO.csv"
data <- read_delim(file = path_donnees, delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Clean the table
data_all=data[4:425,]
data_all$Codes<-NULL
colnames(data_all)=c("Date_MY","IPI")


# Converting date to the correct format
data_all$Date<-as.Date(paste(data_all$Date_MY,1,sep="-"), format = "%Y-%m-%d")
data_all$Date_MY<-NULL

# Extracting month and year from date 
Year=as.numeric(format(data_all$Date, format = "%Y"))
Month=format(data_all$Date, format = "%m")
data_all=cbind(data_all,Year,Month)
data_all$Year<-as.numeric(data_all$Year)
data_all$IPI<-as.numeric(data_all$IPI)
sort(data_all$IPI)
data_all<-data_all[order(data_all$Date),] 
rownames(data_all) <- seq(length=nrow(data_all)) 


# Create a time series for easier analysis, and plot
Xt.ts<-ts(data_all$IPI,start=c(1990,1), end=c(2025,2), frequency=12)
par(mfrow=c(1,1))
plot.ts(Xt.ts, xlab="Years", ylab="IPI (CVS-CJO)")


####### - Part 1 : Preliminary Analysis - ######

####### - Part 1.1 : Trend/Cycle/Seasonal Decomposition - ######

# Decomposing data into seasonal + trend + random
decompose_Xt = decompose(Xt.ts,"additive")

plot(as.ts(decompose_Xt$seasonal))
plot(as.ts(decompose_Xt$trend))
plot(as.ts(decompose_Xt$random))
plot(decompose_Xt)

####### - Part 1.2 : Tranforming the whole series - ######

# Set up a 2-row, 1-column layout
par(mfrow = c(2, 1))

# 1. Plot original series
plot(Xt.ts, type = "l")

# 2. First difference to remove trend
Xt_diff1 <- diff(Xt.ts)

# 3. Plot differenced series
plot(Xt_diff1, type = "l")

# 4. Check stationarity 
# With Augmented Dickey-Fuller Test (ADF)
adf <- adfTest(Xt.ts, lag=1, type="ct") #
adf
# With Philipps-Perron Test (PP)
pp.test(Xt_diff1)
# With Kwiatkowski-Phillips-Schmidt-Shin Test (KPSS)
kpss_test <- kpss.test(Xt_diff1)
print("KPSS Test Results:")
print(kpss_test)


## We then save ACF plot
png("acf_plot_Xt.png")  
acf(Xt_diff1, main = "")
dev.off()

# Save PACF plot
png("pacf_plot_Xt.png")
pacf(Xt_diff1, main = "")
dev.off()


# As we find a noticeable change in the data before and after January 2009,
# We split the data into two subsets of the time series
Xt1 <- window(Xt.ts, end=c(2008, 12))
Xt2 <- window(Xt.ts, start=c(2009, 1))

# And plot them
plot.ts(Xt1, xlab="Years", ylab="IPI",main="1990-2008")
plot.ts(Xt2, xlab="Years", ylab="IPI",main="2009-2025")


####### - Part 1.3 : Tranforming Xt1 - ######

decompose_Xt1 = decompose(Xt1,"additive")

plot(as.ts(decompose_Xt1$seasonal))
plot(as.ts(decompose_Xt1$trend))
plot(as.ts(decompose_Xt1$random))
plot(decompose_Xt1)


# Set up a 2-row, 1-column layout
par(mfrow = c(2, 1))

# 1. Plot original series
plot(Xt1)

# 2. First difference to remove trend
Xt1_diff1 <- Xt.ts

# 3. Plot differenced series
plot(Xt1_diff1)

# 4. Check stationarity 
# With Augmented Dickey-Fuller Test (ADF)
adf <- adfTest(Xt1_diff1, lag=14, type="ct") #
adf
# With Philipps-Perron Test (PP)
pp.test(Xt1_diff1)
# With Kwiatkowski-Phillips-Schmidt-Shin Test (KPSS)
kpss_test <- kpss.test(Xt1_diff1)
print("KPSS Test Results:")
print(kpss_test)



## We then print ACF 
png("acf_plot_Xt1.png")
acf(Xt1_diff1, main = "")
dev.off()

# Save PACF plot
png("pacf_plot_Xt1.png")
pacf(Xt1_diff1, main = "")
dev.off()

####### - Part 1.4 : Tranforming Xt2 - ######

decompose_Xt2 = decompose(Xt2,"additive")

plot(as.ts(decompose_Xt2$seasonal))
plot(as.ts(decompose_Xt2$trend))
plot(as.ts(decompose_Xt2$random))
plot(decompose_Xt2)

# Set up a 2-row, 1-column layout
par(mfrow = c(2, 1))

# 1. Plot original series
plot(Xt2, main = "Original Series", type = "l")

# 2. First difference to remove trend
Xt2_diff1 <- diff(Xt2)

# 3. Plot differenced series
plot(Xt2_diff1, main = "Differenced Series", type = "l")

# 4. Check stationarity 
# With Augmented Dickey-Fuller Test (ADF)
adf <- adfTest(Xt2_diff1, lag = 14, type = "ct")
adf

# With Phillips-Perron Test (PP)
pp.test(Xt2_diff1)

# With Kwiatkowski-Phillips-Schmidt-Shin Test (KPSS)
kpss_test <- kpss.test(Xt2_diff1)
print("KPSS Test Results:")
print(kpss_test)


## We then print ACF 
png("acf_plot_Xt2.png") 
acf(Xt2_diff1, main = "")
dev.off()

# Save PACF plot
png("pacf_plot_Xt2.png")
pacf(Xt2_diff1, main = "")
dev.off()

####### - Part 1.5 : Simultaneous visualisation - ######

# Save a grid of both series and their differentiated versions
png("timeseries_grid.png", width = 800, height = 800)

# Set up 2x2 plotting area
par(mfrow = c(2, 2))

# Plot the series
plot(Xt1, main = "Xt1")
plot(Xt2, main = "Xt2")
plot(Xt1_diff1, main = "Xt1_diff1")
plot(Xt2_diff1, main = "Xt2_diff1")

# Close the graphics device
dev.off()
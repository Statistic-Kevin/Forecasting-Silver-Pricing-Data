# Import libraries
library(dplyr)
library(zoo)
library(tidyverse)
library(lubridate)
library(forecast)
library(MLmetrics)
library(tseries)#adf test
library(TSstudio)
library(padr)
library(imputeTS)
library(DataExplorer)
library(Rbeast)#Bayesian time series decomposition for changepoint, trend, and periodicity or seasonality
library(tsoutliers)#automatic detection of outliers in time series
library(autostsm)#Automatic model selection for structural time series decomposition into trend, cycle, and seasonal components using the Kalman filter.

# Load and prepare data
Silver <- read.csv("C:/Users/user/Desktop/Msc Statistics/Time Series Analysis/Project/LBMA-SILVER.csv", header = TRUE)

# Function to check data structure and summary
check_data <- function(data) {
  #introduce(data)
  str(data)
  summary(data)
  dim(data)
}
check_data(Silver)#4 variables (3 numerical and character)
# Changing Date from character to date (daily date)
Silver <- Silver %>% 
  mutate(Date = ymd(Date)) %>%
  arrange(Date)
check_data(Silver)
# Reinspect data
summarize_missing_values <- function(df) {
  df %>%
    summarise_all(~sum(is.na(.))) %>%
    gather(key = "Variable", value = "MissingValues") %>%
    arrange(desc(MissingValues))
}
summarize_missing_values(Silver)#No missing values in GBP and date unlike EUros

#Check for date intervals
count_date_intervals <- function(df, date_column) {
  df %>%
    arrange(!!sym(date_column)) %>%
    mutate(Interval = c(NA, diff(ymd(!!sym(date_column))))) %>%
    filter(Interval > 1) %>% 
    nrow()
}

count_date_intervals(Silver,"Date")#2871 intervals that are greater than 1


# Pad and clean data
Silver <- Silver %>% pad(interval = "day")## Pad the data to fill in missing dates
summarize_missing_values(Silver)#6098 values owing to date padding to be imputed based on last observation carried forward

Silver_clean<-Silver %>% 
  dplyr::select(Date, GBP) %>% #Select argument was not being used
  na_locf()

summarize_missing_values(Silver_clean)#No missing values for price in GBP and Dates


# Time series creation
Silver_ts <- ts(Silver_clean$GBP, start = c(1968, 1, 2), frequency = 365)
autoplot(Silver_ts)

#<<<<<<<<<<<<<<<<<<<<<<Data Inspection>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#The data seems to have trend that is not linear.Seasonality is however not clear but seems to be present at certain time periods.

p <- ggplot(Silver_clean, aes(x = Date, y = GBP)) + 
  geom_line(color = "#00AFBB") + ylab("Silver Price in GBP") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle=90)) +
  theme_classic()
#Smooth Line for trend
p + stat_smooth(
  color = "#FC4E07", fill = "#FC4E07",
  method = "loess" 
) 
#Detecting Seasonality
f<-decompose(Silver_ts);f#suggests additive seasonality

# Decompose and inspect seasonality
Silver_decompose <- decompose(Silver_ts, type = "additive")
adjust_silver <- Silver_ts / Silver_decompose$seasonal

par(mfrow = c(2, 1))
plot(Silver_ts, main = "Original Time Series")
plot(adjust_silver, main = "Seasonally Adjusted Time Series")
#Given the increasing amplitude of the seasonal effects with the trend level in the original time series, multiplicative seasonality is likely more appropriate
# Decompose and inspect seasonality
Silver_decompose <- decompose(Silver_ts, type = "multiplicative")
adjust_silver <- Silver_ts / Silver_decompose$seasonal
par(mfrow=c(2,1))
plot(Silver_ts, main = "Original Time Series")
plot(adjust_silver, main = "Seasonally Adjusted Time Series")

par(mfrow=c(2,1))
acf(Silver_ts,main="ACF of Original Time Series")
pacf(Silver_ts,main="PACF of Original Time Series")

# Bayesian time series decomposition
out1 <- beast(Silver_ts, start = c(1968, 1, 2), season = 'harmonic')#Season is harmonic owing to the time periodicity of 365 days
plot(out1, vars = c('s', 't'), main = "Bayesian Time Series Decomposition")#Plotting Trend and Seasonality
#The plot shows significant abrupt disruptions in the late 1970's  to the mid 1980's  as well as 2010 to 2011 relating to sharp silver prices around period of recession i.e (https://www.bullionbypost.co.uk/index/silver/silver-price-during-recession/)
#variance often change over time in time series, i.e., time-series data suffer from a distribution shift problem. This change in temporal distribution is one of the main challenges that prevent accurate time-series forecasting. 
#To address this issue,(Taesung et al, 2022) proposes a simple yet effective normalization method called reversible instance normalization (RevIN), a generally-applicable normalization-and-denormalization method with learnable affine transformation.(https://openreview.net/forum?id=cGDAkQo1C0p) 
#The proposed method is symmetrically structured to remove and restore the statistical information of a time-series instance, leading to significant performance improvements in time-series forecasting,
#However naively discarding time points where data had disruptions by only considering the data from 2012 onwards. Taking this decision is highly based on simplicity of analysis, as well as that the fact that recession can in no way be considered a sure indicator of an impending price rise for silver.
#hencetherefore modelling done with this data is still valid

# Filter data post-2012
New_silver <- Silver_clean %>% filter(Date >= "2012-01-01")
New_Silver_ts <- ts(New_silver$GBP, start = c(2012, 1, 1), frequency = 365)
autoplot(New_Silver_ts)

# Plot new data
q <- ggplot(New_silver, aes(x = Date, y = GBP)) + 
  geom_line(color = "#00AFBB") + ylab("Silver Price in GBP") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_classic() +
  stat_smooth(color = "#FC4E07", fill = "#FC4E07", method = "loess")
print(q)

#Checking for seasonality
decomposed <- decompose(New_Silver_ts, type="additive")
plot(decomposed)
ggseasonplot(New_Silver_ts,col=rainbow(12), year.labels = TRUE,main="Seasonal Plot Decomposition")#No apparent seasonal component across over the years.

# Time series modeling
#We will split our data into train and test and we will try to forecast last 6 months price.
#6 months=7*4*6=168 which represents 0.86% 
train <- head(New_Silver_ts, -168)
test <- tail(New_Silver_ts, 168)
acf(train,main="ACF Plot of Train Time Series")

New_silver_diff <- diff(train, lag = 1)
autoplot(New_silver_diff, main = "Differenced Train Time Series (lag=1)")#the series is centred
par(mfrow=c(2,1))
acf(New_silver_diff,main="ACF of Differenced Series")#The differenced time series drops to zero after the first lag showing stationarity
pacf(New_silver_diff,main="PACF of Differenced Series")
mean(New_silver_diff)#The mean (-0.0002223922) is close to zero implying a centred time series
#Pararmetric Tests For STationarity
kpss.test(New_silver_diff, null = "Trend")#the p.value is 0.1 hence fail to reject stationarity no trend (https://www.statology.org/kpss-test-in-r/)
adf.test(New_silver_diff)#P value is 0.01 less than alpha=0.05 hence we  reject null hypothesis of non-stationarity
# Model fitting
obj1 <- auto.arima(y = New_silver_diff, d = 1, D = 0, max.p = 2, max.q = 5, seasonal = FALSE, trace = TRUE)
#Fitting Model
obj1 <- Arima(y = train, order = c(1, 1, 0), seasonal = c(0, 0, 0), include.constant = FALSE)
#Model Coefficients and respective var
obj1$coef
obj1$sigma2
obj1$var.coef 
# The tsdiag function is compatible with Arima :
tsdiag(obj1, main = "Model Diagnostic Plot")#P value is 0.991 not auto corrlation
#Testing Auto Correlation
#H0: residual has no-autocorrelation
#   Vs
#H1: residual has autocorrelation

#Normality of residuals
Box.test(obj1$residuals, type = "Ljung-Box")#p- value 
#Testing for Normality of Residuals
#H0: residual normally distributed
#  Vs
#H1: residual not normally distributed
shapiro.test(obj1$residuals)#p value <0.005 hence we fail to reject the null hypothesis of normality

# Forecasting
#Forecasting
#Main goal is to forecast 6 months
h<-7*4*6;h

Silver_forecast <- forecast(object=obj1,h=168)

# Visualization
New_Silver_ts %>% 
  autoplot() +
  autolayer(obj1$fitted,lwd=0.5,
            series = "ARIMA Model") +
  autolayer(Silver_forecast$mean,lwd=0.5,
            series="Forecast 6 Months")

# Calculate MAPE
train_error <- MAPE(obj1$fitted, train) * 100
test_error <- MAPE(Silver_forecast$mean, test) * 100

print(paste("Train MAPE: ", round(train_error, 2), "%"))
print(paste("Test MAPE: ", round(test_error, 2), "%"))

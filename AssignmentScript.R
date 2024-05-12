##### Data ###########

#Clean Global Environment
rm(list=ls())

#Set Working Directory 
setwd("/Users/bryanbastian/Desktop/Monash /ETF3600")

#Setup library
library(readxl)
library(dplyr)
library(tidyverse)
library(tangram)
library(lmtest)
library(sandwich)
library(car)
library(maxLik)
library(jtools)
library(pscl)
library(nnet)
library(fastDummies)
library(mlbench)
library(DescTools)
library(readxl)
library(olsrr)
library(ggpubr)
library(MASS)
library(pROC)
library(pscl)

#Load data needed
Stock_Price <- read.csv("SP500.csv", header = TRUE)
Volatility <- read.csv("VXVCLS.csv", header = TRUE)
Risk_Free <- read.csv("DTB3.csv", header = TRUE)
Extra <- read.csv("DDFUELNYH.csv", header = TRUE)

#Change the header name
colnames(Stock_Price) <- c("Date", "Price")
colnames(Volatility) <- c("Date", "VIX")
colnames(Risk_Free) <-c("Date", "Rf")
colnames(Extra) <- c("Date", "Ex")

#Merge the dataset
Stock <- merge(Stock_Price, Volatility, by.x ="Date", by.y = "Date")
Stock <- merge(Stock, Risk_Free, by.x = "Date", by.y = "Date")
Stock <- merge(Stock, Extra, by.x = "Date", by.y = "Date")

#Mark the missing value to NA
Stock[is.na(Stock) | Stock=="."] = NA
sum(is.na(Stock$Price))

#Convert DATE to Date class:
Stock[, 1] <- as.Date(Stock[, 1], format="%Y-%m-%d")

#Convert SP500, VIX, risk-free, and Diesel Fuel Price to numeric:
Stock[, 2] <- as.numeric(Stock[, 2])
Stock[, 3] <- as.numeric(Stock[, 3])
Stock[, 4] <- as.numeric(Stock[, 4])
Stock[, 5] <- as.numeric(Stock[, 5])

#Create time lag for VIX
Stock$VIX_1 <- lag(Stock$VIX,1)
Stock$VIX_2 <- lag(Stock$VIX,2)

#Create squared variables 
Stock$VIXsq <- (Stock$VIX)^2
Stock$VIX1sq <- (Stock$VIX_1)^2
Stock$VIX2sq <- (Stock$VIX_2)^2

#Defining annualised log returns, annualised log price rate, and excess return
Stock$Price_1 <- lag(Stock$Price,1)
Stock$anlogr <- 1200*(log(Stock$Price)-log(Stock$Price_1))

Stock$excess <- Stock$anlogr - Stock$Rf

Stock$Ex_1 <- lag(Stock$Ex, 1)
Stock$anlogp <- 1200*(log(Stock$Ex)-log(Stock$Ex_1))

#Defining binary models and trichotomous variable
Stock$rawreturn <- ifelse(Stock$anlogr>0,1,0)
Stock$excessreturn <- ifelse(Stock$excess>0,1,0)
Stock$tri <- ifelse(Stock$anlogr<0,1, ifelse(Stock$excess<=0, 2, ifelse(Stock$excess>0, 3, NA)))

#Filter the NA
tidy_Stock <- Stock
for (col in names(tidy_Stock)[-1]) {
  # Filter out rows where the column has NA values
  tidy_Stock <- tidy_Stock |> filter(!is.na(.data[[col]]))
}
remove(col)
summary(tidy_Stock)

#Find the difference in VIX

tidy_Stock$diff_VIX <- tidy_Stock$VIX - tidy_Stock$VIX_1
tidy_Stock$diff_VIX_2 <- tidy_Stock$VIX - tidy_Stock$VIX_2
tidy_Stock$diff_VIXsq <- tidy_Stock$diff_VIX^2
tidy_Stock$diff_VIXsq2 <- tidy_Stock$diff_VIX_2^2

##### Modelling#######
attach(tidy_Stock)

##### Checking correlation between explanatory variables####
tidy_Stock$logVIX <- log(VIX)
tidy_Stock$logVIX1 <- log(VIX_1)
tidy_Stock$logVIX2 <- log(VIX_2)
tidy_Stock$logdiff1 <- log(VIX/VIX_1)
tidy_Stock$logdiff2 <- log(VIX/VIX_2)

covtable <- round(cor(tidy_Stock[-1]), 3)
View(covtable)

##### RawReturnLogit######

#Create empty dataframe to store goodness of fit
RawLogit <- data.frame(
  AIC = numeric(),
  McFadden = numeric(),
  HitRate = numeric())

RawLModels <- list(
  Rlogit_v1 = glm(rawreturn~VIX, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v2 = glm(rawreturn~VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v3 = glm(rawreturn~VIX+VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v4 = glm(rawreturn~log(VIX), data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v5 = glm(rawreturn~VIX_1, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v6 = glm(rawreturn~VIX1sq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v7 = glm(rawreturn~VIX_1+VIX1sq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v8 = glm(rawreturn~log(VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v9 = glm(rawreturn~VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v10 = glm(rawreturn~VIX2sq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v11 = glm(rawreturn~VIX_2+VIX2sq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v12 = glm(rawreturn~log(VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v13 = glm(rawreturn~diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v14 = glm(rawreturn~diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v15 = glm(rawreturn~diff_VIX+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v16 = glm(rawreturn~log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v17 = glm(rawreturn~diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v18 = glm(rawreturn~diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v19 = glm(rawreturn~diff_VIX_2+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v20 = glm(rawreturn~log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v21 = glm(rawreturn~VIX+diff_VIX, data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v22 = glm(rawreturn~VIXsq+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v23 = glm(rawreturn~log(VIX)+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v24 = glm(rawreturn~VIX_1+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v25 = glm(rawreturn~VIX1sq+diff_VIX, data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v26 = glm(rawreturn~log(VIX_1)+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v27 = glm(rawreturn~VIX_2+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v28 = glm(rawreturn~VIX2sq+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v29 = glm(rawreturn~log(VIX_2)+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v30 = glm(rawreturn~VIX+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v31 = glm(rawreturn~VIXsq+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v32 = glm(rawreturn~log(VIX)+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v33 = glm(rawreturn~VIX_1+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v34 = glm(rawreturn~VIX1sq+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v35 = glm(rawreturn~log(VIX_1)+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v36 = glm(rawreturn~VIX_2+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v37 = glm(rawreturn~VIX2sq+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v38 = glm(rawreturn~log(VIX_2)+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v39 = glm(rawreturn~VIX+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v40 = glm(rawreturn~VIXsq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v41 = glm(rawreturn~log(VIX)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v42 = glm(rawreturn~VIX_1+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v43 = glm(rawreturn~VIX1sq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v44 = glm(rawreturn~log(VIX_1)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v45 = glm(rawreturn~VIX_2+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v46 = glm(rawreturn~VIX2sq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v47 = glm(rawreturn~log(VIX_2)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v48 = glm(rawreturn~VIX+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v49 = glm(rawreturn~VIXsq+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v50 = glm(rawreturn~log(VIX)+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v51 = glm(rawreturn~VIX_1+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v52 = glm(rawreturn~VIX1sq+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v53 = glm(rawreturn~log(VIX_1)+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v54 = glm(rawreturn~VIX_2+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v55 = glm(rawreturn~VIX2sq+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v56 = glm(rawreturn~log(VIX_2)+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v57 = glm(rawreturn~VIX+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v58 = glm(rawreturn~VIXsq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v59 = glm(rawreturn~log(VIX)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v60 = glm(rawreturn~VIX_1+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v61 = glm(rawreturn~VIX1sq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v62 = glm(rawreturn~log(VIX_1)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v63 = glm(rawreturn~VIX_2+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v64 = glm(rawreturn~VIX2sq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v65 = glm(rawreturn~log(VIX_2)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  
  Rlogit_v66 = glm(rawreturn~VIX+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v67 = glm(rawreturn~VIXsq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v68 = glm(rawreturn~log(VIX)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v69 = glm(rawreturn~VIX_1+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v70 = glm(rawreturn~VIX1sq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")), 
  Rlogit_v71 = glm(rawreturn~log(VIX_1)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v72 = glm(rawreturn~VIX_2+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v73 = glm(rawreturn~VIX2sq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Rlogit_v74 = glm(rawreturn~log(VIX_2)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")))

for (i in 1:length(RawLModels)){
  model = RawLModels[[i]]
  model_AIS = summary(model)
  model_McFadden = pR2(model)
  model_hitmiss = hitmiss(model)
  
  
  RawLogit <- rbind(RawLogit, data.frame(
    AIC = model_AIS$aic,
    McFadden = model_McFadden[4],
    HitRate = model_hitmiss[1]
  ), make.row.names = FALSE)
}

####Check the table and find the summary of the best goodness of fit###
View(RawLogit)
best_goodness_fit = glm(rawreturn~log(VIX)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit"))
summary(best_goodness_fit)

##### RawReturnProbit#####
#Create Empty Data Frame
RawProbit <- data.frame(
  AIC = numeric(),
  McFadden = numeric(),
  HitRate = numeric())

RawPModels <- list(
  RProbit_v1 = glm(rawreturn~VIX, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v2 = glm(rawreturn~VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v3 = glm(rawreturn~VIX+VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v4 = glm(rawreturn~log(VIX), data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v5 = glm(rawreturn~VIX_1, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v6 = glm(rawreturn~VIX1sq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v7 = glm(rawreturn~VIX_1+VIX1sq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v8 = glm(rawreturn~log(VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v9 = glm(rawreturn~VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v10 = glm(rawreturn~VIX2sq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v11 = glm(rawreturn~VIX_2+VIX2sq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v12 = glm(rawreturn~log(VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v13 = glm(rawreturn~diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v14 = glm(rawreturn~diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v15 = glm(rawreturn~diff_VIX+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v16 = glm(rawreturn~log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v17 = glm(rawreturn~diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v18 = glm(rawreturn~diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v19 = glm(rawreturn~diff_VIX_2+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v20 = glm(rawreturn~log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v21 = glm(rawreturn~VIX+diff_VIX, data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v22 = glm(rawreturn~VIXsq+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v23 = glm(rawreturn~log(VIX)+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v24 = glm(rawreturn~VIX_1+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v25 = glm(rawreturn~VIX1sq+diff_VIX, data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v26 = glm(rawreturn~log(VIX_1)+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v27 = glm(rawreturn~VIX_2+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v28 = glm(rawreturn~VIX2sq+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v29 = glm(rawreturn~log(VIX_2)+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v30 = glm(rawreturn~VIX+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v31 = glm(rawreturn~VIXsq+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v32 = glm(rawreturn~log(VIX)+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v33 = glm(rawreturn~VIX_1+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v34 = glm(rawreturn~VIX1sq+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v35 = glm(rawreturn~log(VIX_1)+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v36 = glm(rawreturn~VIX_2+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v37 = glm(rawreturn~VIX2sq+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v38 = glm(rawreturn~log(VIX_2)+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v39 = glm(rawreturn~VIX+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v40 = glm(rawreturn~VIXsq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v41 = glm(rawreturn~log(VIX)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v42 = glm(rawreturn~VIX_1+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v43 = glm(rawreturn~VIX1sq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v44 = glm(rawreturn~log(VIX_1)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v45 = glm(rawreturn~VIX_2+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v46 = glm(rawreturn~VIX2sq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v47 = glm(rawreturn~log(VIX_2)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v48 = glm(rawreturn~VIX+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v49 = glm(rawreturn~VIXsq+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v50 = glm(rawreturn~log(VIX)+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v51 = glm(rawreturn~VIX_1+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v52 = glm(rawreturn~VIX1sq+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v53 = glm(rawreturn~log(VIX_1)+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v54 = glm(rawreturn~VIX_2+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v55 = glm(rawreturn~VIX2sq+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v56 = glm(rawreturn~log(VIX_2)+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v57 = glm(rawreturn~VIX+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v58 = glm(rawreturn~VIXsq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v59 = glm(rawreturn~log(VIX)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v60 = glm(rawreturn~VIX_1+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v61 = glm(rawreturn~VIX1sq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v62 = glm(rawreturn~log(VIX_1)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v63 = glm(rawreturn~VIX_2+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v64 = glm(rawreturn~VIX2sq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v65 = glm(rawreturn~log(VIX_2)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  
  RProbit_v66 = glm(rawreturn~VIX+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v67 = glm(rawreturn~VIXsq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v68 = glm(rawreturn~log(VIX)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v69 = glm(rawreturn~VIX_1+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v70 = glm(rawreturn~VIX1sq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")), 
  RProbit_v71 = glm(rawreturn~log(VIX_1)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v72 = glm(rawreturn~VIX_2+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v73 = glm(rawreturn~VIX2sq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  RProbit_v74 = glm(rawreturn~log(VIX_2)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")))

for (i in 1:length(RawPModels)){
  model = RawPModels[[i]]
  model_AIS = summary(model)
  model_McFadden = pR2(model)
  model_hitmiss = hitmiss(model)
  
  
  RawProbit <- rbind(RawProbit, data.frame(
    AIC = round(model_AIS$aic,3),
    McFadden = round(model_McFadden[4]*100,3),
    HitRate = round(model_hitmiss[1],3)
  ), make.row.names = FALSE)
}

####Check the table and find the summary of the best goodness of fit###
View(RawProbit)
best_fit_raw_probit = glm(rawreturn~log(VIX)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit"))

##### ExcessReturnLogit#####
#For Excess Logit 
ExLogit <- data.frame(
  AIC = numeric(),
  McFadden = numeric(),
  HitRate = numeric())

ExLModels <- list(
  Elogit_v1 = glm(excessreturn~VIX, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v2 = glm(excessreturn~VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v3 = glm(excessreturn~VIX+VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v4 = glm(excessreturn~log(VIX), data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v5 = glm(excessreturn~VIX_1, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v6 = glm(excessreturn~VIX1sq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v7 = glm(excessreturn~VIX_1+VIX1sq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v8 = glm(excessreturn~log(VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v9 = glm(excessreturn~VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v10 = glm(excessreturn~VIX2sq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v11 = glm(excessreturn~VIX_2+VIX2sq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v12 = glm(excessreturn~log(VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v13 = glm(excessreturn~diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v14 = glm(excessreturn~diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v15 = glm(excessreturn~diff_VIX+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v16 = glm(excessreturn~log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v17 = glm(excessreturn~diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v18 = glm(excessreturn~diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v19 = glm(excessreturn~diff_VIX_2+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v20 = glm(excessreturn~log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v21 = glm(excessreturn~VIX+diff_VIX, data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v22 = glm(excessreturn~VIXsq+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v23 = glm(excessreturn~log(VIX)+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v24 = glm(excessreturn~VIX_1+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v25 = glm(excessreturn~VIX1sq+diff_VIX, data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v26 = glm(excessreturn~log(VIX_1)+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v27 = glm(excessreturn~VIX_2+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v28 = glm(excessreturn~VIX2sq+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v29 = glm(excessreturn~log(VIX_2)+diff_VIX, data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v30 = glm(excessreturn~VIX+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v31 = glm(excessreturn~VIXsq+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v32 = glm(excessreturn~log(VIX)+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v33 = glm(excessreturn~VIX_1+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v34 = glm(excessreturn~VIX1sq+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v35 = glm(excessreturn~log(VIX_1)+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v36 = glm(excessreturn~VIX_2+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v37 = glm(excessreturn~VIX2sq+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v38 = glm(excessreturn~log(VIX_2)+diff_VIXsq, data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v39 = glm(excessreturn~VIX+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v40 = glm(excessreturn~VIXsq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v41 = glm(excessreturn~log(VIX)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v42 = glm(excessreturn~VIX_1+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v43 = glm(excessreturn~VIX1sq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v44 = glm(excessreturn~log(VIX_1)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v45 = glm(excessreturn~VIX_2+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v46 = glm(excessreturn~VIX2sq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v47 = glm(excessreturn~log(VIX_2)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v48 = glm(excessreturn~VIX+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v49 = glm(excessreturn~VIXsq+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v50 = glm(excessreturn~log(VIX)+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v51 = glm(excessreturn~VIX_1+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v52 = glm(excessreturn~VIX1sq+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v53 = glm(excessreturn~log(VIX_1)+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v54 = glm(excessreturn~VIX_2+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v55 = glm(excessreturn~VIX2sq+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v56 = glm(excessreturn~log(VIX_2)+diff_VIX_2, data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v57 = glm(excessreturn~VIX+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v58 = glm(excessreturn~VIXsq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v59 = glm(excessreturn~log(VIX)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v60 = glm(excessreturn~VIX_1+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v61 = glm(excessreturn~VIX1sq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v62 = glm(excessreturn~log(VIX_1)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v63 = glm(excessreturn~VIX_2+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v64 = glm(excessreturn~VIX2sq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v65 = glm(excessreturn~log(VIX_2)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="logit")),
  
  Elogit_v66 = glm(excessreturn~VIX+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v67 = glm(excessreturn~VIXsq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v68 = glm(excessreturn~log(VIX)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v69 = glm(excessreturn~VIX_1+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v70 = glm(excessreturn~VIX1sq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")), 
  Elogit_v71 = glm(excessreturn~log(VIX_1)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v72 = glm(excessreturn~VIX_2+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v73 = glm(excessreturn~VIX2sq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")),
  Elogit_v74 = glm(excessreturn~log(VIX_2)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="logit")))



for (i in 1:length(ExLModels)){
  model = ExLModels[[i]]
  model_AIS = summary(model)
  model_McFadden = pR2(model)
  model_hitmiss = hitmiss(model)
  
  
  ExLogit <- rbind(ExLogit, data.frame(
    AIC = round(model_AIS$aic,3),
    McFadden = round(model_McFadden[4]*100,3),
    HitRate = round(model_hitmiss[1],3)
  ), make.row.names = FALSE)
}

View(ExLogit)

##### ExcessReturnProbit######
#For Excess Probit
ExProbit <- data.frame(
  AIC = numeric(),
  McFadden = numeric(),
  HitRate = numeric())

ExPModels <- list(
  ExProbit_v1 = glm(excessreturn~VIX, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v2 = glm(excessreturn~VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v3 = glm(excessreturn~VIX+VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v4 = glm(excessreturn~log(VIX), data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v5 = glm(excessreturn~VIX_1, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v6 = glm(excessreturn~VIX1sq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v7 = glm(excessreturn~VIX_1+VIX1sq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v8 = glm(excessreturn~log(VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v9 = glm(excessreturn~VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v10 = glm(excessreturn~VIX2sq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v11 = glm(excessreturn~VIX_2+VIX2sq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v12 = glm(excessreturn~log(VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v13 = glm(excessreturn~diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v14 = glm(excessreturn~diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v15 = glm(excessreturn~diff_VIX+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v16 = glm(excessreturn~log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v17 = glm(excessreturn~diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v18 = glm(excessreturn~diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v19 = glm(excessreturn~diff_VIX_2+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v20 = glm(excessreturn~log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v21 = glm(excessreturn~VIX+diff_VIX, data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v22 = glm(excessreturn~VIXsq+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v23 = glm(excessreturn~log(VIX)+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v24 = glm(excessreturn~VIX_1+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v25 = glm(excessreturn~VIX1sq+diff_VIX, data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v26 = glm(excessreturn~log(VIX_1)+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v27 = glm(excessreturn~VIX_2+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v28 = glm(excessreturn~VIX2sq+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v29 = glm(excessreturn~log(VIX_2)+diff_VIX, data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v30 = glm(excessreturn~VIX+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v31 = glm(excessreturn~VIXsq+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v32 = glm(excessreturn~log(VIX)+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v33 = glm(excessreturn~VIX_1+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v34 = glm(excessreturn~VIX1sq+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v35 = glm(excessreturn~log(VIX_1)+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v36 = glm(excessreturn~VIX_2+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v37 = glm(excessreturn~VIX2sq+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v38 = glm(excessreturn~log(VIX_2)+diff_VIXsq, data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v39 = glm(excessreturn~VIX+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v40 = glm(excessreturn~VIXsq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v41 = glm(excessreturn~log(VIX)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v42 = glm(excessreturn~VIX_1+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v43 = glm(excessreturn~VIX1sq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v44 = glm(excessreturn~log(VIX_1)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v45 = glm(excessreturn~VIX_2+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v46 = glm(excessreturn~VIX2sq+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v47 = glm(excessreturn~log(VIX_2)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v48 = glm(excessreturn~VIX+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v49 = glm(excessreturn~VIXsq+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v50 = glm(excessreturn~log(VIX)+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v51 = glm(excessreturn~VIX_1+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v52 = glm(excessreturn~VIX1sq+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v53 = glm(excessreturn~log(VIX_1)+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v54 = glm(excessreturn~VIX_2+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v55 = glm(excessreturn~VIX2sq+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v56 = glm(excessreturn~log(VIX_2)+diff_VIX_2, data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v57 = glm(excessreturn~VIX+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v58 = glm(excessreturn~VIXsq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v59 = glm(excessreturn~log(VIX)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v60 = glm(excessreturn~VIX_1+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v61 = glm(excessreturn~VIX1sq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v62 = glm(excessreturn~log(VIX_1)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v63 = glm(excessreturn~VIX_2+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v64 = glm(excessreturn~VIX2sq+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v65 = glm(excessreturn~log(VIX_2)+diff_VIXsq2, data=tidy_Stock, family=binomial(link="probit")),
  
  ExProbit_v66 = glm(excessreturn~VIX+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v67 = glm(excessreturn~VIXsq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v68 = glm(excessreturn~log(VIX)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v69 = glm(excessreturn~VIX_1+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v70 = glm(excessreturn~VIX1sq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")), 
  ExProbit_v71 = glm(excessreturn~log(VIX_1)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v72 = glm(excessreturn~VIX_2+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v73 = glm(excessreturn~VIX2sq+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")),
  ExProbit_v74 = glm(excessreturn~log(VIX_2)+log(VIX/VIX_2), data=tidy_Stock, family=binomial(link="probit")))

for (i in 1:length(ExPModels)){
  model = ExPModels[[i]]
  model_AIS = summary(model)
  model_McFadden = pR2(model)
  model_hitmiss = hitmiss(model)
  
  
  ExProbit <- rbind(ExProbit, data.frame(
    AIC = round(model_AIS$aic,3),
    McFadden = round(model_McFadden[4]*100,3),
    HitRate = round(model_hitmiss[1],3)
  ), make.row.names = FALSE)
}

View(ExProbit)

best_fit_ex_probit =glm(excessreturn~log(VIX)+log(VIX/VIX_1), data=tidy_Stock, family=binomial(link="probit"))

##### TriLogit#####
####Ordered Logit

tri = as.factor(tri)
levels(tri)

TriLogit <- data.frame(
  AIC = numeric(),
  McFadden = numeric(),
  HitMiss = numeric())

TriLModels <- list(
  TriLogit_v1 = polr(tri~VIX, method ="logistic"),
  TriLogit_v2 = polr(tri~VIXsq, method ="logistic"),
  TriLogit_v3 = polr(tri~VIX+VIXsq, method ="logistic"),
  TriLogit_v4 = polr(tri~log(VIX), method ="logistic"),
  
  TriLogit_v5 = polr(tri~VIX_1, method ="logistic"),
  TriLogit_v6 = polr(tri~VIX1sq, method ="logistic"),
  TriLogit_v7 = polr(tri~VIX_1+VIX1sq, method ="logistic"),
  TriLogit_v8 = polr(tri~log(VIX_1), method ="logistic"),
  
  TriLogit_v9 = polr(tri~VIX_2, method ="logistic"),
  TriLogit_v10 = polr(tri~VIX2sq, method ="logistic"),
  TriLogit_v11 = polr(tri~VIX_2+VIX2sq, method ="logistic"),
  TriLogit_v12 = polr(tri~log(VIX_2), method ="logistic"),
  
  TriLogit_v13 = polr(tri~diff_VIX, method ="logistic"),
  TriLogit_v14 = polr(tri~diff_VIXsq, method ="logistic"),
  TriLogit_v15 = polr(tri~diff_VIX+diff_VIXsq, method ="logistic"),
  TriLogit_v16 = polr(tri~log(VIX/VIX_1), method ="logistic"),
  
  TriLogit_v17 = polr(tri~diff_VIX_2, method ="logistic"),
  TriLogit_v18 = polr(tri~diff_VIXsq2, method ="logistic"),
  TriLogit_v19 = polr(tri~diff_VIX_2+diff_VIXsq2, method ="logistic"),
  TriLogit_v20 = polr(tri~log(VIX/VIX_2), method ="logistic"),
  
  TriLogit_v21 = polr(tri~VIX+diff_VIX, method ="logistic"), 
  TriLogit_v22 = polr(tri~VIXsq+diff_VIX, method ="logistic"),
  TriLogit_v23 = polr(tri~log(VIX)+diff_VIX, method ="logistic"),
  TriLogit_v24 = polr(tri~VIX_1+diff_VIX, method ="logistic"),
  TriLogit_v25 = polr(tri~VIX1sq+diff_VIX, method ="logistic"), 
  TriLogit_v26 = polr(tri~log(VIX_1)+diff_VIX, method ="logistic"),
  TriLogit_v27 = polr(tri~VIX_2+diff_VIX, method ="logistic"),
  TriLogit_v28 = polr(tri~VIX2sq+diff_VIX, method ="logistic"),
  TriLogit_v29 = polr(tri~log(VIX_2)+diff_VIX, method ="logistic"),
  
  TriLogit_v30 = polr(tri~VIX+diff_VIXsq, method ="logistic"), 
  TriLogit_v31 = polr(tri~VIXsq+diff_VIXsq, method ="logistic"),
  TriLogit_v32 = polr(tri~log(VIX)+diff_VIXsq, method ="logistic"),
  TriLogit_v33 = polr(tri~VIX_1+diff_VIXsq, method ="logistic"),
  TriLogit_v34 = polr(tri~VIX1sq+diff_VIXsq, method ="logistic"), 
  TriLogit_v35 = polr(tri~log(VIX_1)+diff_VIXsq, method ="logistic"),
  TriLogit_v36 = polr(tri~VIX_2+diff_VIXsq, method ="logistic"),
  TriLogit_v37 = polr(tri~VIX2sq+diff_VIXsq, method ="logistic"),
  TriLogit_v38 = polr(tri~log(VIX_2)+diff_VIXsq, method ="logistic"),
  
  TriLogit_v39 = polr(tri~VIX+log(VIX/VIX_1), method ="logistic"), 
  TriLogit_v40 = polr(tri~VIXsq+log(VIX/VIX_1), method ="logistic"),
  TriLogit_v41 = polr(tri~log(VIX)+log(VIX/VIX_1), method ="logistic"),
  TriLogit_v42 = polr(tri~VIX_1+log(VIX/VIX_1), method ="logistic"),
  TriLogit_v43 = polr(tri~VIX1sq+log(VIX/VIX_1), method ="logistic"), 
  TriLogit_v44 = polr(tri~log(VIX_1)+log(VIX/VIX_1), method ="logistic"),
  TriLogit_v45 = polr(tri~VIX_2+log(VIX/VIX_1), method ="logistic"),
  TriLogit_v46 = polr(tri~VIX2sq+log(VIX/VIX_1), method ="logistic"),
  TriLogit_v47 = polr(tri~log(VIX_2)+log(VIX/VIX_1), method ="logistic"),
  
  TriLogit_v48 = polr(tri~VIX+diff_VIX_2, method ="logistic"), 
  TriLogit_v49 = polr(tri~VIXsq+diff_VIX_2, method ="logistic"),
  TriLogit_v50 = polr(tri~log(VIX)+diff_VIX_2, method ="logistic"),
  TriLogit_v51 = polr(tri~VIX_1+diff_VIX_2, method ="logistic"),
  TriLogit_v52 = polr(tri~VIX1sq+diff_VIX_2, method ="logistic"), 
  TriLogit_v53 = polr(tri~log(VIX_1)+diff_VIX_2, method ="logistic"),
  TriLogit_v54 = polr(tri~VIX_2+diff_VIX_2, method ="logistic"),
  TriLogit_v55 = polr(tri~VIX2sq+diff_VIX_2, method ="logistic"),
  TriLogit_v56 = polr(tri~log(VIX_2)+diff_VIX_2, method ="logistic"),
  
  TriLogit_v57 = polr(tri~VIX+diff_VIXsq2, method ="logistic"), 
  TriLogit_v58 = polr(tri~VIXsq+diff_VIXsq2, method ="logistic"),
  TriLogit_v59 = polr(tri~log(VIX)+diff_VIXsq2, method ="logistic"),
  TriLogit_v60 = polr(tri~VIX_1+diff_VIXsq2, method ="logistic"),
  TriLogit_v61 = polr(tri~VIX1sq+diff_VIXsq2, method ="logistic"), 
  TriLogit_v62 = polr(tri~log(VIX_1)+diff_VIXsq2, method ="logistic"),
  TriLogit_v63 = polr(tri~VIX_2+diff_VIXsq2, method ="logistic"),
  TriLogit_v64 = polr(tri~VIX2sq+diff_VIXsq2, method ="logistic"),
  TriLogit_v65 = polr(tri~log(VIX_2)+diff_VIXsq2, method ="logistic"),
  
  TriLogit_v66 = polr(tri~VIX+log(VIX/VIX_2), method ="logistic"), 
  TriLogit_v67 = polr(tri~VIXsq+log(VIX/VIX_2), method ="logistic"),
  TriLogit_v68 = polr(tri~log(VIX)+log(VIX/VIX_2), method ="logistic"),
  TriLogit_v69 = polr(tri~VIX_1+log(VIX/VIX_2), method ="logistic"),
  TriLogit_v70 = polr(tri~VIX1sq+log(VIX/VIX_2), method ="logistic"), 
  TriLogit_v71 = polr(tri~log(VIX_1)+log(VIX/VIX_2), method ="logistic"),
  TriLogit_v72 = polr(tri~VIX_2+log(VIX/VIX_2), method ="logistic"),
  TriLogit_v73 = polr(tri~VIX2sq+log(VIX/VIX_2), method ="logistic"),
  TriLogit_v74 = polr(tri~log(VIX_2)+log(VIX/VIX_2), method ="logistic"))  

for (i in 1:length(TriLModels)){
  model = TriLModels[[i]]
  model_AIC = AIC(model)
  model_McFadden = pR2(model)[4]
  predicted_classes <- predict(model, tidy_Stock, type = "class")
  actual_classes <- tidy_Stock$tri  
  correct_predictions <- predicted_classes == actual_classes
  accuracy <- sum(correct_predictions) / length(correct_predictions) * 100
  
  TriLogit <- rbind(TriLogit, data.frame(
    AIC = model_AIC,
    McFadden = model_McFadden,
    HitMiss = accuracy
  ), make.row.names = FALSE)
}

##### TriProbit######
#Ordered Probit
TriProbit <- data.frame(
  AIC = numeric(),
  McFadden = numeric(),
  HitMiss = numeric())

TriPModels <- list(
  TriProbit_v1 = polr(tri~VIX, method ="probit"),
  TriProbit_v2 = polr(tri~VIXsq, method ="probit"),
  TriProbit_v3 = polr(tri~VIX+VIXsq, method ="probit"),
  TriProbit_v4 = polr(tri~log(VIX), method ="probit"),
  
  TriProbit_v5 = polr(tri~VIX_1, method ="probit"),
  TriProbit_v6 = polr(tri~VIX1sq, method ="probit"),
  TriProbit_v7 = polr(tri~VIX_1+VIX1sq, method ="probit"),
  TriProbit_v8 = polr(tri~log(VIX_1), method ="probit"),
  
  TriProbit_v9 = polr(tri~VIX_2, method ="probit"),
  TriProbit_v10 = polr(tri~VIX2sq, method ="probit"),
  TriProbit_v11 = polr(tri~VIX_2+VIX2sq, method ="probit"),
  TriProbit_v12 = polr(tri~log(VIX_2), method ="probit"),
  
  TriProbit_v13 = polr(tri~diff_VIX, method ="probit"),
  TriProbit_v14 = polr(tri~diff_VIXsq, method ="probit"),
  TriProbit_v15 = polr(tri~diff_VIX+diff_VIXsq, method ="probit"),
  TriProbit_v16 = polr(tri~log(VIX/VIX_1), method ="probit"),
  
  TriProbit_v17 = polr(tri~diff_VIX_2, method ="probit"),
  TriProbit_v18 = polr(tri~diff_VIXsq2, method ="probit"),
  TriProbit_v19 = polr(tri~log(VIX/VIX_2), method ="probit"),
  
  TriProbit_v20 = polr(tri~VIX+diff_VIX, method ="probit"), 
  TriProbit_v21 = polr(tri~VIXsq+diff_VIX, method ="probit"),
  TriProbit_v22 = polr(tri~log(VIX)+diff_VIX, method ="probit"),
  TriProbit_v23 = polr(tri~VIX_1+diff_VIX, method ="probit"),
  TriProbit_v24 = polr(tri~VIX1sq+diff_VIX, method ="probit"), 
  TriProbit_v25 = polr(tri~log(VIX_1)+diff_VIX, method ="probit"),
  TriProbit_v26 = polr(tri~VIX_2+diff_VIX, method ="probit"),
  TriProbit_v27 = polr(tri~VIX2sq+diff_VIX, method ="probit"),
  TriProbit_v28 = polr(tri~log(VIX_2)+diff_VIX, method ="probit"),
  
  TriProbit_v29 = polr(tri~VIX+diff_VIXsq, method ="probit"), 
  TriProbit_v30 = polr(tri~VIXsq+diff_VIXsq, method ="probit"),
  TriProbit_v31 = polr(tri~log(VIX)+diff_VIXsq, method ="probit"),
  TriProbit_v32 = polr(tri~VIX_1+diff_VIXsq, method ="probit"),
  TriProbit_v33 = polr(tri~VIX1sq+diff_VIXsq, method ="probit"), 
  TriProbit_v34= polr(tri~log(VIX_1)+diff_VIXsq, method ="probit"),
  TriProbit_v35 = polr(tri~VIX_2+diff_VIXsq, method ="probit"),
  TriProbit_v36 = polr(tri~VIX2sq+diff_VIXsq, method ="probit"),
  TriProbit_v37 = polr(tri~log(VIX_2)+diff_VIXsq, method ="probit"),
  
  TriProbit_v38 = polr(tri~VIX+log(VIX/VIX_1), method ="probit"), 
  TriProbit_v39 = polr(tri~VIXsq+log(VIX/VIX_1), method ="probit"),
  TriProbit_v40 = polr(tri~log(VIX)+log(VIX/VIX_1), method ="probit"),
  TriProbit_v41 = polr(tri~VIX_1+log(VIX/VIX_1), method ="probit"),
  TriProbit_v42 = polr(tri~VIX1sq+log(VIX/VIX_1), method ="probit"), 
  TriProbit_v43 = polr(tri~log(VIX_1)+log(VIX/VIX_1), method ="probit"),
  TriProbit_v44 = polr(tri~VIX_2+log(VIX/VIX_1), method ="probit"),
  TriProbit_v45 = polr(tri~VIX2sq+log(VIX/VIX_1), method ="probit"),
  TriProbit_v46 = polr(tri~log(VIX_2)+log(VIX/VIX_1), method ="probit"),
  
  TriProbit_v47 = polr(tri~VIX+diff_VIX_2, method ="probit"), 
  TriProbit_v48 = polr(tri~VIXsq+diff_VIX_2, method ="probit"),
  TriProbit_v49 = polr(tri~log(VIX)+diff_VIX_2, method ="probit"),
  TriProbit_v50 = polr(tri~VIX_1+diff_VIX_2, method ="probit"),
  TriProbit_v51 = polr(tri~VIX1sq+diff_VIX_2, method ="probit"), 
  TriProbit_v52 = polr(tri~log(VIX_1)+diff_VIX_2, method ="probit"),
  TriProbit_v53 = polr(tri~VIX_2+diff_VIX_2, method ="probit"),
  TriProbit_v54 = polr(tri~VIX2sq+diff_VIX_2, method ="probit"),
  TriProbit_v55 = polr(tri~log(VIX_2)+diff_VIX_2, method ="probit"),
  
  TriProbit_v56 = polr(tri~VIX+diff_VIXsq2, method ="probit"), 
  TriProbit_v57 = polr(tri~VIXsq+diff_VIXsq2, method ="probit"),
  TriProbit_v58 = polr(tri~log(VIX)+diff_VIXsq2, method ="probit"),
  TriProbit_v59 = polr(tri~VIX_1+diff_VIXsq2, method ="probit"),
  TriProbit_v60 = polr(tri~VIX1sq+diff_VIXsq2, method ="probit"), 
  TriProbit_v61 = polr(tri~log(VIX_1)+diff_VIXsq2, method ="probit"),
  TriProbit_v62 = polr(tri~VIX_2+diff_VIXsq2, method ="probit"),
  TriProbit_v63 = polr(tri~VIX2sq+diff_VIXsq2, method ="probit"),
  TriProbit_v64 = polr(tri~log(VIX_2)+diff_VIXsq2, method ="probit"),
  
  TriProbit_v65 = polr(tri~VIX+log(VIX/VIX_2), method ="probit"), 
  TriProbit_v66 = polr(tri~VIXsq+log(VIX/VIX_2), method ="probit"),
  TriProbit_v67 = polr(tri~log(VIX)+log(VIX/VIX_2), method ="probit"),
  TriProbit_v68 = polr(tri~VIX_1+log(VIX/VIX_2), method ="probit"),
  TriProbit_v69 = polr(tri~VIX1sq+log(VIX/VIX_2), method ="probit"), 
  TriProbit_v70 = polr(tri~log(VIX_1)+log(VIX/VIX_2), method ="probit"),
  TriProbit_v71 = polr(tri~VIX_2+log(VIX/VIX_2), method ="probit"),
  TriProbit_v72 = polr(tri~VIX2sq+log(VIX/VIX_2), method ="probit"),
  TriProbit_v73 = polr(tri~log(VIX_2)+log(VIX/VIX_2), method ="probit")) 

for (i in 1:length(TriPModels)){
  model = TriPModels[[i]]
  model_AIC = AIC(model)
  model_McFadden = pR2(model)[4]
  predicted_classes <- predict(model, tidy_Stock, type = "class")
  actual_classes <- tidy_Stock$tri  
  correct_predictions <- predicted_classes == actual_classes
  accuracy <- sum(correct_predictions) / length(correct_predictions) * 100
  
  TriProbit <- rbind(TriProbit, data.frame(
    AIC = model_AIC,
    McFadden = model_McFadden,
    HitMiss = accuracy
  ), make.row.names = FALSE)
}

##### Considering Explanatory Variable#####
RProbit_v41 = glm(rawreturn~log(VIX)+log(VIX/VIX_1)+anlogp, data=tidy_Stock, family=binomial(link="probit"))
summary(RProbit_v41)
pR2(RProbit_v41)
hitmiss(RProbit_v41)

##### Plotting####
#ROC Curve

library("pROC")
tidy_Stock$predicted_prob <- predict(best_fit_raw_probit, type = "response")
test_roc = roc(tidy_Stock$rawreturn ~ tidy_Stock$predicted_prob, plot = TRUE, 
               main="Raw Return Probit Model",
              xlab="False Positive Rate", ylab="True Positive Rate", print.auc = TRUE)
as.numeric(test_roc$auc)

tidy_Stock$predicted_again <- predict(best_fit_ex_probit, type = "response")
test_roc_2 = roc(tidy_Stock$excessreturn ~ tidy_Stock$predicted_again, plot = TRUE, 
               main="Excess Return Probit Model",
               xlab="False Positive Rate", ylab="True Positive Rate", print.auc = TRUE)
as.numeric(test_roc_2$auc)


tidy_Stock$predicted_again_2 <- predict(RProbit_v41, type = "response")
test_roc_3 = roc(tidy_Stock$rawreturn ~ tidy_Stock$predicted_again_2, plot = TRUE, 
                 main="Proposed Explanatory Model",
                 xlab="False Positive Rate", ylab="True Positive Rate", print.auc = TRUE)
as.numeric(test_roc_3$auc)

barplot(table(tri), main="Frequency of Each Category",
        ylab = "Frequency",
        xlab = "Category")

#### Export table
install.packages("writexl")
library("writexl")
write_xlsx(ExProbit, "ExProbit.xlsx")

write.csv(RawProbit, "RawReturnProbit.csv")

################################################################################
### (2) Model analysis R code for:

### "The relative influence of cross-seasonal and local weather effects on the 
### breeding success of a migratory songbird"

###  DR de Zwaan, A Drake, AF Camfield, EC MacDonald, and K Martin (2022)

###  Journal of Animal Ecology	

### Run with R version 4.1.1

################################################################################


### Set working directory
setwd("")

### Load required packages

library(plyr)
library(dplyr)
library(lme4)
library(doBy)
library(corrplot)
library(Hmisc)
library(RColorBrewer)
library(MuMIn)
library(lmerTest)
library(ggplot2)

### Read in best climate predictors for clutch initiation and dev time (code 01)

ci_predictors <- read.csv("ci_weather_summarized.csv")
dev_predictors <- read.csv("dev_weather_summarized.csv")

### Read in nest data

nest_data <- read.csv('Nest data.csv')

### Subset nest data

nest_data <- nest_data %>% select(nest_id, year, clutch_initiation, age_at_fledge:brood_size)


### Merge nest data with predictor datasets

ci_nest_data <- merge(nest_data, ci_predictors, by=c("nest_id","year","clutch_initiation","age_at_fledge",
                                                     "exposure_days","brood_size"), all.y=TRUE)

dev_nest_data <- merge(nest_data, dev_predictors, by=c("nest_id","year","clutch_initiation","age_at_fledge",
                                                       "exposure_days","brood_size"), all.y=TRUE)


#################
### Summarize response variables among and within years
#################

### Clutch initiation date

ci_summary <- summaryBy(clutch_initiation ~ year, data = ci_nest_data, FUN=c(min, max, mean, sd,length))
ci_summary$se <- ci_summary$clutch_initiation.sd/sqrt(ci_summary$clutch_initiation.length)

### Minimum and maximum clutch initiation date across years

min(ci_summary$clutch_initiation.mean)
max(ci_summary$clutch_initiation.mean)

### Annual range
ci_summary$range <- ci_summary$clutch_initiation.max - ci_summary$clutch_initiation.min
ci_summary

median(ci_summary$range)
min(ci_summary$range)
max(ci_summary$range)


### Development time

### Remove 2010, as only has 1 nest

dev_nest_data <- subset(dev_nest_data, year!=2010)

### Summarize

dev_summary <- summaryBy(age_at_fledge ~ year, data=dev_nest_data, FUN=c(min, max, mean, sd,length))
dev_summary$se <- dev_summary$age_at_fledge.sd/sqrt(dev_summary$age_at_fledge.length)

### Minimum and maximum clutch initiation date across years

min(dev_summary$age_at_fledge.mean)
max(dev_summary$age_at_fledge.mean)



##########################
### Add in snow depth data
##########################

### Read in annual snow depth measured on April 15th.

snow_data <- read.csv("Snow depth data.csv")

head(snow_data)

### Rename column headers

colnames(snow_data) <- c("year","snow_depth_m")

### Merge with CI and dev data

ci_nest_data <- merge(ci_nest_data, snow_data, by=c("year"), all.x=TRUE)

dev_nest_data <- merge(dev_nest_data, snow_data, by=c("year"), all.x=TRUE)


###################
### Scale variables
###################

ci_nest_data_sc <- scale(ci_nest_data[,c(9:20)], center=TRUE, scale=TRUE)

dev_nest_data_sc <- scale(dev_nest_data[,c(9:12)], center=TRUE, scale=TRUE)


### Make year into a factor

ci_nest_data$year_f <- factor(ci_nest_data$year)
dev_nest_data$year_f <- factor(dev_nest_data$year)

# Recombine with unscaled values

ci_nest_data_sc <- as.data.frame(ci_nest_data_sc)
dev_nest_data_sc <- as.data.frame(dev_nest_data_sc)

non_scaled_var <- ci_nest_data %>% select(year,nest_id,clutch_initiation,year_f)
non_scaled_var_dev <- dev_nest_data %>% select(year,nest_id,age_at_fledge,clutch_initiation,brood_size,year_f)

ci_nest_data_use <- bind_cols(non_scaled_var,ci_nest_data_sc)
dev_nest_data_use <- bind_cols(non_scaled_var_dev,dev_nest_data_sc)


#############
### Assess correlations among predictors
#############

ci_cor_use <- ci_nest_data_use %>% select(min_avg_temp_w:min_avg_temp_s,
                                          max_avg_temp_s,avg_temp_b,
                                          min_avg_temp_b,snow_depth_m)

dev_cor_use <- dev_nest_data_use %>% select(freeze_days_s,precip_days_b,snow_depth_m)


### Rename columns for both datasets

colnames(ci_cor_use) <- c("Winter min avg temp", "Stopover avg min temp",
                          "Stopover freeze days", "Stopover precip", 
                          "Stopover min avg temp", "Stopover max avg temp", 
                          "Breeding avg temp", "Breeding min avg temp", 
                          "Snow depth")

colnames(dev_cor_use) <- c("Stopover freeze days",
                          "Breeding precip", "Snow depth")


### Clutch initiation correlation matrix

ci_res2 <- rcorr(as.matrix(ci_cor_use), type="pearson")
ci_res2$P <- ifelse(is.na(ci_res2$P) == TRUE, 1, ci_res2$P)

# Create plot

windowsFonts(Times=windowsFont("Times New Roman"))

png(file="correlogram_CI.png",width=3750,height=2750, res=600)
corrplot(ci_res2$r, method = "circle", type="lower", order="original", 
         col=brewer.pal(n=6, name="RdBu"), p.mat = ci_res2$P, sig.level = 0.001, 
         insig = "n", tl.col="black", addgrid.col=NA, rect.col = "black", 
         cl.pos="b", cl.cex = 0.9, tl.srt = 0.5, 
         tl.offset = 0.7, tl.cex = 0.4, addCoef.col = "black", 
         number.digits = 2, number.cex = 0.90)
dev.off()



### Development correlation matrix

dev_res2 <- rcorr(as.matrix(dev_cor_use), type="pearson")
dev_res2$P <- ifelse(is.na(dev_res2$P) == TRUE, 1, dev_res2$P)

# Create plot

windowsFonts(Times=windowsFont("Times New Roman"))

png(file="correlogram_dev.png",width=3750,height=2750, res=600)
corrplot(dev_res2$r, method = "circle", type="lower", order="original", 
         col=brewer.pal(n=6, name="RdBu"), p.mat = dev_res2$P, sig.level = 0.001, 
         insig = "n", tl.col="black", addgrid.col=NA, rect.col = "black", 
         cl.pos="b", cl.cex = 0.9, tl.srt = 0.5, 
         tl.offset = 0.7, tl.cex = 0.4, addCoef.col = "black", 
         number.digits = 2, number.cex = 0.85)
dev.off()


################################################################################
### Model selection and model averaging
################################################################################


#####################
### Clutch initiation
#####################


### Fit global model

ci_global <- lmer(clutch_initiation ~ min_avg_temp_w + avg_min_temp_s + 
                        freeze_days_s +
                        precip_s + min_avg_temp_s + 
                        max_avg_temp_s + avg_temp_b + 
                        min_avg_temp_b +  
                        snow_depth_m + year + 
                       (1|year_f), data=ci_nest_data_use)


summary(ci_global)


### Set global optiosn to na.fail
options(na.action = "na.fail")

### Set cor matrix that says TRUE for all combinations except highly correlated one
ci_cor <- ci_nest_data_use %>% select(min_avg_temp_w:min_avg_temp_s,
                                      max_avg_temp_s, avg_temp_b,
                                      min_avg_temp_b, snow_depth_m)

ci_cor_use <- abs(cor(ci_cor)) <= .6 
ci_cor_use[!lower.tri(ci_cor_use)] <- NA

### Fit all possible sub-models
ci_dredged <- dredge(ci_global, beta = "sd", evaluate = TRUE,
                     rank = "AIC", fixed = "year", subset = ci_cor_use) # subset using the matrix created above


length(ci_dredged$delta) # 60 submodels

### Summarize list of top models
ci_top_models <- subset(ci_dredged, delta <= 2)
ci_top_models

### Model average
ci_avg_model <- model.avg(ci_dredged, delta <= 2)
summary(ci_avg_model)
confint(ci_avg_model, level = 0.85)


######################
### Development time
######################

### Resolve missing values with average

dev_nest_data_use$brood_size <- ifelse(is.na(dev_nest_data_use$brood_size) == TRUE,
                                       4, dev_nest_data_use$brood_size)

### Scale clutch initiation

dev_nest_data_use$clutch_initiation_sc <- scale(dev_nest_data_use$clutch_initiation, center=T, scale=T)


### Fit global model

dev_global <- lmer(age_at_fledge ~ clutch_initiation_sc + brood_size + freeze_days_s + 
                      precip_days_b + snow_depth_m +  year + (1|year_f), data=dev_nest_data_use)


summary(dev_global)

### Set global option to na.fail
options(na.action = "na.fail")

### Create cor matrix
dev_cor <- dev_nest_data_use %>% dplyr::select(freeze_days_s, precip_days_b, snow_depth_m)
dev_cor_use <- abs(cor(dev_cor)) <= .6 
dev_cor_use[!lower.tri(dev_cor_use)] <- NA

### Fit all sub models
dev_dredged <- dredge(dev_global, beta = "sd", evaluate = TRUE,
                     rank = "AIC", fixed = "year", 
                     subset = dev_cor_use) 


length(dev_dredged$delta) # 32 submodels

### Summarize top models
dev_top_models <- subset(dev_dredged, delta <= 2)
dev_top_models

### Model averaeg
dev_avg_model <- model.avg(dev_dredged, delta <= 2)
summary(dev_avg_model)
confint(dev_avg_model, level = 0.85)


################################################################################
### Predation risk
################################################################################

head(nest_data)

### Create fledging success variable

nest_data$success <- ifelse(is.na(nest_data$age_at_fledge),0,1)

### Standardize clutch initiation across years

nest_data$clutch_initiation_overall_std <- nest_data$clutch_initiation - mean(nest_data$clutch_initiation)

### Standardize clutch initiation within each year

ci_2003 <- subset(nest_data, year==2003)
ci_2004 <- subset(nest_data, year==2004)
ci_2005 <- subset(nest_data, year==2005)
ci_2006 <- subset(nest_data, year==2006)
ci_2007 <- subset(nest_data, year==2007)
ci_2010 <- subset(nest_data, year==2010)
ci_2011 <- subset(nest_data, year==2011)
ci_2015 <- subset(nest_data, year==2015)
ci_2016 <- subset(nest_data, year==2016)
ci_2017 <- subset(nest_data, year==2017)
ci_2018 <- subset(nest_data, year==2018)
ci_2019 <- subset(nest_data, year==2019)


ci_2003$clutch_initiation_year_std <- ci_2003$clutch_initiation - mean(ci_2003$clutch_initiation)
ci_2004$clutch_initiation_year_std <- ci_2004$clutch_initiation - mean(ci_2004$clutch_initiation)
ci_2005$clutch_initiation_year_std <- ci_2005$clutch_initiation - mean(ci_2005$clutch_initiation)
ci_2006$clutch_initiation_year_std <- ci_2006$clutch_initiation - mean(ci_2006$clutch_initiation)
ci_2007$clutch_initiation_year_std <- ci_2007$clutch_initiation - mean(ci_2007$clutch_initiation)
ci_2010$clutch_initiation_year_std <- ci_2010$clutch_initiation - mean(ci_2010$clutch_initiation)
ci_2011$clutch_initiation_year_std <- ci_2011$clutch_initiation - mean(ci_2011$clutch_initiation)
ci_2015$clutch_initiation_year_std <- ci_2015$clutch_initiation - mean(ci_2015$clutch_initiation)
ci_2016$clutch_initiation_year_std <- ci_2016$clutch_initiation - mean(ci_2016$clutch_initiation)
ci_2017$clutch_initiation_year_std <- ci_2017$clutch_initiation - mean(ci_2017$clutch_initiation)
ci_2018$clutch_initiation_year_std <- ci_2018$clutch_initiation - mean(ci_2018$clutch_initiation)
ci_2019$clutch_initiation_year_std <- ci_2019$clutch_initiation - mean(ci_2019$clutch_initiation)


### Combine together

ci_int <- rbind(ci_2003, ci_2004)
ci_int <- rbind(ci_int, ci_2005)
ci_int <- rbind(ci_int, ci_2006)
ci_int <- rbind(ci_int, ci_2007)
ci_int <- rbind(ci_int, ci_2010)
ci_int <- rbind(ci_int, ci_2011)
ci_int <- rbind(ci_int, ci_2015)
ci_int <- rbind(ci_int, ci_2016)
ci_int <- rbind(ci_int, ci_2017)
ci_int <- rbind(ci_int, ci_2018)
ci_int <- rbind(ci_int, ci_2019)


ci_use <- ci_int


##################################################
### Daily nest survival analysis
##################################################

### Reset global action

options(na.action = "na.omit")

### Convert year into a factor

ci_use$year_f <- factor(ci_use$year)

### Create alternate complementary log-log exposure models

linear_model <- glmer(success ~ clutch_initiation_year_std + offset(log(exposure_days)) + (1|year_f), 
                   family = binomial(link="cloglog"), data = ci_use) 

quadratic_model <- glmer(success ~ poly(clutch_initiation_year_std,2) + offset(log(exposure_days)) + (1|year_f), 
                   family = binomial(link="cloglog"), data = ci_use) 

### Model comparison
AICc(linear_model, quadratic_model)

summary(linear_model)


### Create back-calculation function to estimate DNS and overall survival

cloglogit<-function( l ) 1-exp(-exp(l))


### Average nest success and DNS across years 

summary(linear_model)

cloglogit(-3.71844 + -0.05798*0)

### Survival
1-(1-0.02397961)^24  # 44.2%

### DNS
0.4415126^(1/24) #0.967



#################
### Annual DNS estimates
#################

annual_DNS <- glm(success ~ clutch_initiation_year_std + year_f +  offset(log(exposure_days)), 
                       family = binomial(link="cloglog"), data = ci_use) 
summary(annual_DNS)


### Estimate overall survival percentage and DNS for each year using a loop 

DNS_vec <- c(0, 0.74857, 0.34597, -0.38619, -17.26569, -2.79573, -1.32106,
              -1.26104, -1.13012, -0.75955, -1.40932, 0.23081)

years <- c(2003,2004,2005,2006,2007,2010,2011,2015,2016,2017,2018,2019)

overall_survival <- c()
DNS_year <- c()


for (i in 1:length(DNS_vec)) {
  
      # Annual estimates
      DNS_int <- cloglogit(-2.91328 + DNS_vec[i])
      overall_survival[i] <- 1-(1-DNS_int)^24
      DNS_year[i] <- overall_survival[i]^(1/24)
      
   }



### Create a dataframe
yearly_DNS_list <- list(years,DNS_year,overall_survival)
yearly_DNS_df <- bind_cols(yearly_DNS_list)
colnames(yearly_DNS_df) <- c("year", "DNS", "survival")

yearly_DNS_df



#########
### Bootstrap estimates for survival SE
#########

### Remove NAs from exposure days first

ci_use_new <- subset(ci_use, exposure_days > 0)

prob <- c(11/347, 15/347, 37/347, 45/347, 13/347, 20/347, 20/347, 56/347, 53/347,
          35/347, 32/347, 10/347)

generated_pop <- data.frame(clutch_initiation_year_std = rnorm(1000, mean(ci_use_new$clutch_initiation_year_std), sd(ci_use_new$clutch_initiation_year_std)),
                            exposure_days = rnorm(1000, mean(ci_use_new$exposure_days), sd(ci_use_new$exposure_days)),
                            year =  sample(years, 1000, replace=TRUE, prob=prob)) 


### Remove the few negative exposure days

generated_pop_new <- subset(generated_pop, exposure_days > 0)

### Convert year to a factor

generated_pop_new$year_f <- factor(generated_pop_new$year)

summary(annual_DNS)

### Standard error
predicted <- predict(annual_DNS, newdata = generated_pop_new,
                     type = "response", na.action=na.pass)

generated_pop_new$predicted_dns <- predicted

### Table
predicted_table_y <- summaryBy(predicted_dns ~ year_f, data=generated_pop_new, FUN=c(mean,sd,length))
predicted_table_y$survival_se <- predicted_table_y$predicted_dns.sd/sqrt(predicted_table_y$predicted_dns.length)
predicted_table_y


###########
### Estimates for DNS se (based on estimated confidence interval)
###########

DNS_vec <- yearly_DNS_df$DNS
survival_vec <- yearly_DNS_df$survival*100
survival_se <- predicted_table_y$survival_se*100

### Calculate the upper and lower CI for prob of survival
survival_lower_ci <- survival_vec - (survival_se*1.96)  
survival_upper_ci <- survival_vec + (survival_se*1.96)

### Convert to DNS CI
DNS_lower_ci <- survival_lower_ci^(1/24)
DNS_upper_ci <- survival_upper_ci^(1/24)

# For 2007 because creates an NA
9.8^(1/24) #1.099768
0^(1/24) #0.00

### Average the interval

DNS_ci_interval <- (DNS_upper_ci - DNS_lower_ci)/2

#For 2007
1.099768/2 # 0.549884

### Convert from CI to se

DNS_se <- DNS_ci_interval/1.96

# For 2007
DNS_se[5] <- 0.549884/1.96

### Table

predicted_table_y$DNS_se <- DNS_se
predicted_table_y


#############################
### Predation risk simulation
#############################

### Generate bootstrapped population

ci_use_new <- subset(ci_use, exposure_days>0)

generated_pop_new <- data.frame(clutch_initiation_year_std = rnorm(1000, median(ci_use_new$clutch_initiation_year_std), sd(ci_use_new$clutch_initiation_year_std)),
                            exposure_days = rnorm(1000, median(ci_use_new$exposure_days), sd(ci_use_new$exposure_days)))

head(generated_pop_new)

### Round values because these are integers

generated_pop_new$clutch_initiation_year_std <- round(generated_pop_new$clutch_initiation_year_std)
generated_pop_new$exposure_days <- round(generated_pop_new$exposure_days)

### Remove sub zero exposure days which are impossible

generated_pop_new <- subset(generated_pop_new, exposure_days > 0)


### Fit same best model from above

linear_model_predict <- glm(success ~ clutch_initiation_year_std +  offset(log(exposure_days)), 
                       family = binomial(link="cloglog"), data = ci_use) 
summary(linear_model_predict)


### Now predict values

predicted_new <- predict(linear_model_predict, newdata = generated_pop_new,
                     type = "response",
                     se.fit = TRUE, na.action=na.pass)

### Save predictions and SE
generated_pop_new$predicted_dns <- predicted_new$fit
generated_pop_new$predicted_dns_se <- predicted_new$se.fit

head(generated_pop_new)

max(generated_pop_new$predicted.dns)

### Back calculate
generated_pop_new$back_calc <- generated_pop_new$predicted_dns^(1/24) 

### Calculate for average (9), fast (7), and prolong (13) time in the nest

generated_pop_new$days_9 <- generated_pop_new$back_calc^9
generated_pop_new$days_7 <- generated_pop_new$back_calc^7
generated_pop_new$days_13 <- generated_pop_new$back_calc^13


head(generated_pop_new)


### Now subset to remove all DFE that are greater than 14 or less than -14 
### (outside observed range)

generated_pop_use <- subset(generated_pop_new, clutch_initiation_year_std < 15 & 
                              clutch_initiation_year_std > -15)

### Calculate survival percentage (used to plot)
generated_pop_use$days_9_prop <- generated_pop_use$days_9 * 100
generated_pop_use$days_13_prop <- generated_pop_use$days_13 * 100
generated_pop_use$days_7_prop <- generated_pop_use$days_7 * 100

head(generated_pop_use)


################################ End of code ###################################
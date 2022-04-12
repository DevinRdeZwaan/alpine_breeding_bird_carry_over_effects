################################################################################
### (1) Sliding window analysis R code for:

### "The relative influence of cross-seasonal and local weather effects on the 
### breeding success of a migratory songbird"

###  DR de Zwaan, A Drake, AF Camfield, EC MacDonald, and K Martin (2022)

###  Journal of Animal Ecology	

### Run with R version 4.1.1

################################################################################

### Set working directory

setwd("")

### Load required packages

library(lubridate)
library(plyr)
library(dplyr)
library(lme4)
library(climwin)
library(doBy)
library(beepr)

### Load climate data

non_breeding_weather <- read.csv("Regional non-breeding weather data.csv")
breeding_weather <- read.csv("Local breeding weather data.csv")

### Split weather data into different non-breeding regions

winter_weather <- subset(non_breeding_weather, region=="winter")
stopover_weather <- subset(non_breeding_weather, region=="stopover")


### Convert date-time to proper format (d/m/y)

winter_weather$date <- ymd(winter_weather$datetime)
winter_weather$date_new <- format(winter_weather$date, "%d/%m/%Y")

stopover_weather$date <- ymd(stopover_weather$datetime)
stopover_weather$date_new <- format(stopover_weather$date, "%d/%m/%Y")

breeding_weather$date <- strptime(paste(breeding_weather$year, breeding_weather$jdate), format="%Y %j") 
breeding_weather$date_new <- format(breeding_weather$date, format= "%d/%m/%Y")


### Create year, month, and day columns

winter_weather$year <- year(winter_weather$date)
winter_weather$month <- month(winter_weather$date)
winter_weather$day <- day(winter_weather$date)

stopover_weather$year <- year(stopover_weather$date)
stopover_weather$month <- month(stopover_weather$date)
stopover_weather$day <- day(stopover_weather$date)

### Create a 'days below freezing' variable

winter_weather$sub_zero <- ifelse(winter_weather$avg_daily_temp <= 0, 1, 0)
stopover_weather$sub_zero <- ifelse(stopover_weather$avg_daily_temp <= 0, 1, 0)

breeding_weather$sub_zero <- ifelse(breeding_weather$temp_avg <=0, 1,0)


### Read in nest data

nest_data <- read.csv('Nest data.csv')

### Convert ordinal date to proper format date-time formate

nest_data$clutch_initiation_int <- strptime(paste(nest_data$year, nest_data$clutch_initiation), format="%Y %j") 
nest_data$first_egg <- format(nest_data$clutch_initiation_int, format= "%d/%m/%Y")

### For age at fledge, select successful nests only

fledge_data <- subset(nest_data, fledge > 0)

### Convert to proper date-time

fledge_data$fledge_int <- strptime(paste(fledge_data$year, fledge_data$fledge), format="%Y %j") 
fledge_data$fledge_date <- format(fledge_data$fledge_int, format= "%d/%m/%Y")

### Convert to proper date-time for hatch date as well

fledge_data$hatch_int <- strptime(paste(fledge_data$year, fledge_data$hatch), format="%Y %j") 
fledge_data$hatch_date <- format(fledge_data$hatch_int, format= "%d/%m/%Y")

head(fledge_data)

### Convert year into a factor

nest_data$year_f <- factor(nest_data$year)
fledge_data$year_f <- factor(fledge_data$year)

### Sample sizes - total

length(nest_data$nest_id) # Number of nests for clutch initiation = 382
length(fledge_data$nest_id) # Number of nests for fledging data = 118

### Sample sizes - annual

summaryBy(nest_id ~ year, data=nest_data, FUN=length)
summaryBy(nest_id ~ year, data=fledge_data, FUN=length)


################################################################################
#### Sliding window analysis
################################################################################

##########################
### Clutch initiation date
##########################

#############
### A) Winter
#############

### Average values (all weather variables)

winter_dfe_sw_avg <- slidingwin(baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data),  
                                xvar = list(TempMean=winter_weather$avg_daily_temp, 
                                            Frozen = winter_weather$sub_zero,
                                            TempMin = winter_weather$min_daily_temp, 
                                            TempMax = winter_weather$max_daily_temp, 
                                            RainTotal = winter_weather$precip_daily_sum), #precip_daily_sum is the hourly precip summed within days
                                type = "absolute", 
                                refday = c(01,06),
                                range = c(182,31), # December 1 to May 1
                                stat = c("mean"), 
                                cmissing = "method1",
                                func = "lin",
                                exclude = c(21,0), # Excludes windows <21 days at any time.
                                cinterval = "day",
                                cdate = winter_weather$date_new, bdate = nest_data$first_egg) 
beep(2)


### Minimum temperature values 
### (temperature variables only, as does not make sense for precipitation and freeze days)

winter_dfe_sw_min <- slidingwin(baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data),  
                                xvar = list(TempMean=winter_weather$avg_daily_temp,
                                            TempMin = winter_weather$min_daily_temp),
                                type = "absolute", 
                                refday = c(01,06),
                                range = c(182,31), # December 1 to May 1
                                stat = c("min"), 
                                cmissing = "method1",
                                func = "lin",
                                exclude = c(21,0), 
                                cinterval = "day",
                                cdate = winter_weather$date_new, bdate = nest_data$first_egg) 
beep(2)


### Maximum temperature values

winter_dfe_sw_max <- slidingwin(baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data),  
                                xvar = list(TempMean=winter_weather$avg_daily_temp,
                                            TempMax = winter_weather$max_daily_temp),
                                type = "absolute", 
                                refday = c(01,06),
                                range = c(182,31), # December 1 to May 1
                                stat = c("max"), 
                                cmissing = "method1",
                                func = "lin",
                                exclude = c(21,0), 
                                cinterval = "day",
                                cdate = winter_weather$date_new, bdate = nest_data$first_egg) 
beep(2)


### Winter results - average

# Temperature

head(winter_dfe_sw_avg[[1]]$Dataset, 40) ## Mean temp (52-31)  (B = -1.44)
head(winter_dfe_sw_avg[[2]]$Dataset, 40) ## Freeze days (none better than null)
head(winter_dfe_sw_avg[[3]]$Dataset, 40) ## Min temp (none better than the null)
head(winter_dfe_sw_avg[[4]]$Dataset, 40) ## Max temp (52-31)  (B = -1.31)

# Precipitation
head(winter_dfe_sw_avg[[5]]$Dataset, 40) ## Precip (102-80; 103-52)  (B = 10.21; 19.23)


### Winter results - min

# Temperature
head(winter_dfe_sw_min[[1]]$Dataset, 40) ## Avg temp (180-159)  (B = -0.54)
head(winter_dfe_sw_min[[2]]$Dataset, 40) ## Min temp (180-159)  (B = -0.48)


### Winter results - max

# Temperature
head(winter_dfe_sw_max[[1]]$Dataset, 40) ## Avg temp (119-84)  (B = -1.30)
head(winter_dfe_sw_max[[2]]$Dataset, 40) ## Max temp (119-84)  (B = -0.96)


##################
### Randomizations
##################

### Winter - average 

winter_dfe_ran_avg <- randwin(repeats = 100, 
                              xvar = list(TempMean=winter_weather$avg_daily_temp, 
                                          #Frozen = winter_weather$sub_zero,
                                          #TempMin = winter_weather$min_daily_temp, 
                                          #TempMax = winter_weather$max_daily_temp, 
                                          RainTotal = winter_weather$precip_daily_sum),
                         cdate = winter_weather$date_new, 
                         bdate = nest_data$first_egg, 
                         baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data), 
                         type = "absolute",  refday = c(01,06),
                         range = c(182,31), stat = c("mean"), 
                         exclude = c(21,0),cmissing = "method1",
                         func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = winter_dfe_sw_avg[[1]]$Dataset, datasetrand = winter_dfe_ran_avg[[1]], 
       metric="AIC", sample.size=382) # < 0.12

# Total precip
pvalue(dataset = winter_dfe_sw_avg[[5]]$Dataset, datasetrand = winter_dfe_ran_avg[[2]], 
       metric="AIC", sample.size=382) # <0.2


### Winter - min

winter_dfe_ran_min <- randwin(repeats = 100, 
                              xvar = list(TempMean=winter_weather$avg_daily_temp,
                                          TempMin = winter_weather$min_daily_temp),
                              cdate = winter_weather$date_new, 
                              bdate = nest_data$first_egg, 
                              baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data), 
                              type = "absolute",  refday = c(01,06),
                              range = c(182,31), stat = c("min"), 
                              exclude = c(21,0),cmissing = "method1",
                              func = "lin", cinterval = "day")

#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = winter_dfe_sw_min[[1]]$Dataset, datasetrand = winter_dfe_ran_min[[1]], 
       metric="AIC", sample.size=382) # 0.04

# Min temp
pvalue(dataset = winter_dfe_sw_min[[2]]$Dataset, datasetrand = winter_dfe_ran_min[[2]], 
       metric="AIC", sample.size=382) # 0.08


### Winter - max

winter_dfe_ran_max <- randwin(repeats = 100, 
                              xvar = list(TempMean=winter_weather$avg_daily_temp,
                                          TempMin = winter_weather$max_daily_temp),
                              cdate = winter_weather$date_new, 
                              bdate = nest_data$first_egg, 
                              baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data), 
                              type = "absolute",  refday = c(01,06),
                              range = c(182,31), stat = c("max"), 
                              exclude = c(21,0),cmissing = "method1",
                              func = "lin", cinterval = "day")

#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = winter_dfe_sw_max[[1]]$Dataset, datasetrand = winter_dfe_ran_max[[1]], 
       metric="AIC", sample.size=382) # 0.12

# Max temp
pvalue(dataset = winter_dfe_sw_max[[2]]$Dataset, datasetrand = winter_dfe_ran_max[[2]], 
       metric="AIC", sample.size=382) # 0.12



#################
### B) Stopover
#################


### Average values (all weather variables)

stopover_dfe_sw_avg <- slidingwin(baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data),  
                                  xvar = list(TempMean= stopover_weather$avg_daily_temp, 
                                              Frozen = stopover_weather$sub_zero,
                                              TempMin = stopover_weather$min_daily_temp, 
                                              TempMax = stopover_weather$max_daily_temp, 
                                              RainTotal = stopover_weather$precip_daily_sum), #precip_daily_sum is the hourly precip summed within days
                                type = "absolute", 
                                refday = c(01,06),
                                range = c(92,31), # March 1 to May 1
                                stat = c("mean"), 
                                cmissing = "method1",
                                func = "lin",
                                exclude = c(21,0),
                                cinterval = "day",
                                cdate = stopover_weather$date_new, bdate = nest_data$first_egg) 
beep(2)


### Minimum temp 
stopover_dfe_sw_min <- slidingwin(baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data),  
                                xvar = list(TempMean=stopover_weather$avg_daily_temp,
                                            TempMin = stopover_weather$min_daily_temp),
                                type = "absolute", 
                                refday = c(01,06),
                                range = c(92,31), # December 1 to May 1
                                stat = c("min"), 
                                cmissing = "method1",
                                func = "lin",
                                exclude = c(21,0),
                                cinterval = "day",
                                cdate = stopover_weather$date_new, bdate = nest_data$first_egg) 
beep(2)


### Maximum temp 
stopover_dfe_sw_max <- slidingwin(baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data),  
                                  xvar = list(TempMean=stopover_weather$avg_daily_temp,
                                              TempMax = stopover_weather$max_daily_temp),
                                  type = "absolute", 
                                  refday = c(01,06),
                                  range = c(92,31), # December 1 to May 1
                                  stat = c("max"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  exclude = c(21,0),
                                  cinterval = "day",
                                  cdate = stopover_weather$date_new, bdate = nest_data$first_egg) 

### Stopover results - average

# Temperature
head(stopover_dfe_sw_avg[[1]]$Dataset, 40) ## Mean temp (67-32)  (B = -1.47)
head(stopover_dfe_sw_avg[[2]]$Dataset, 50) ## Freeze days (66-45) (B = 34.32)
head(stopover_dfe_sw_avg[[3]]$Dataset, 20) ## Min temp (72-51) (B = -1.65)
head(stopover_dfe_sw_avg[[4]]$Dataset, 20) ## Max temp (67-32)  (B = -1.29)

# Precipitation
head(stopover_dfe_sw_avg[[5]]$Dataset, 20) ## Precip (88-52)  (B = 7.94)


### Stopover results - min

# Temperature
head(stopover_dfe_sw_min[[1]]$Dataset, 40) ## Avg temp (60-39)  (B = -1.61)
head(stopover_dfe_sw_min[[2]]$Dataset, 40) ## Min temp (79-47)  (B = -1.78)

### Stopover results - max

# Temperature
head(stopover_dfe_sw_max[[1]]$Dataset, 40) ## Avg temp (52-31)  (B = -0.95)
head(stopover_dfe_sw_max[[2]]$Dataset, 40) ## Max temp (52-31)  (B = -0.82)


##################
### Randomizations
##################

### Stopover - average 

stopover_dfe_ran_avg <- randwin(repeats = 100, 
                              xvar = list(TempMean=stopover_weather$avg_daily_temp, 
                                          Frozen = stopover_weather$sub_zero,
                                          TempMin = stopover_weather$min_daily_temp, 
                                          TempMax = stopover_weather$max_daily_temp, 
                                          RainTotal = stopover_weather$precip_daily_sum),
                              cdate = stopover_weather$date_new, 
                              bdate = nest_data$first_egg, 
                              baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data), 
                              type = "absolute",  refday = c(01,06),
                              range = c(92,31), stat = c("mean"), 
                              exclude = c(21,0),
                              func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = stopover_dfe_sw_avg[[1]]$Dataset, datasetrand = stopover_dfe_ran_avg[[1]], 
       metric="AIC", sample.size=382) # 0.12

# Freeze days
pvalue(dataset = stopover_dfe_sw_avg[[2]]$Dataset, datasetrand = stopover_dfe_ran_avg[[2]], 
       metric="AIC", sample.size=382) # 0.04

# Min temp
pvalue(dataset = stopover_dfe_sw_avg[[3]]$Dataset, datasetrand = stopover_dfe_ran_avg[[3]], 
       metric="AIC", sample.size=382) # 0.08

# Max temp
pvalue(dataset = stopover_dfe_sw_avg[[4]]$Dataset, datasetrand = stopover_dfe_ran_avg[[4]], 
       metric="AIC", sample.size=382) # 0.28

# Total precip
pvalue(dataset = stopover_dfe_sw_avg[[5]]$Dataset, datasetrand = stopover_dfe_ran_avg[[5]], 
       metric="AIC", sample.size=382) # 0.08


### Stopover - min 

stopover_dfe_ran_min <- randwin(repeats = 100, 
                              xvar = list(TempMean= stopover_weather$avg_daily_temp, 
                                          TempMin = stopover_weather$min_daily_temp),
                              cdate = stopover_weather$date_new, 
                              bdate = nest_data$first_egg, 
                              baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data), 
                              type = "absolute",  refday = c(01,06),
                              range = c(92,31), stat = c("min"), 
                              exclude = c(21,0),
                              func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = stopover_dfe_sw_min[[1]]$Dataset, datasetrand = stopover_dfe_ran_min[[1]], 
       metric="AIC", sample.size=382) # 0.08

# Min temp
pvalue(dataset = stopover_dfe_sw_min[[2]]$Dataset, datasetrand = stopover_dfe_ran_min[[2]], 
       metric="AIC", sample.size=382) # 0.08


### Stopover - max 

stopover_dfe_ran_max <- randwin(repeats = 100, 
                                xvar = list(TempMean= stopover_weather$avg_daily_temp, 
                                            TempMax = stopover_weather$max_daily_temp),
                                cdate = stopover_weather$date_new, 
                                bdate = nest_data$first_egg, 
                                baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data), 
                                type = "absolute",  refday = c(01,06),
                                range = c(92,31), stat = c("max"), 
                                exclude = c(21,0),
                                func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = stopover_dfe_sw_max[[1]]$Dataset, datasetrand = stopover_dfe_ran_max[[1]], 
       metric="AIC", sample.size=382) # < 0.001

# Max temp
pvalue(dataset = stopover_dfe_sw_max[[2]]$Dataset, datasetrand = stopover_dfe_ran_max[[2]], 
       metric="AIC", sample.size=382) # 0.04



###########################
### C) Breeding site
###########################


### Average values (all weather variables)

breeding_dfe_sw_avg <- slidingwin(baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data),  
                              xvar = list(TempMean=breeding_weather$temp_avg, 
                                          Frozen = breeding_weather$sub_zero,
                                          RainTotal = breeding_weather$precip_days),
                              type = "absolute", 
                              refday = c(01,06),
                              range= c(47,0),
                              stat = c("mean"), 
                              cmissing = "method1",
                              func = "lin",
                              cinterval = "day",
                              exclude = c(7,0), 
                              cdate = breeding_weather$date_new, bdate = nest_data$first_egg) 
beep(2)


### Minimum values (only temperature)

breeding_dfe_sw_min <- slidingwin(baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data),  
                                  xvar = list(TempMean=breeding_weather$temp_avg),
                                  type = "absolute", 
                                  refday = c(01,06),
                                  range= c(47,0),
                                  stat = c("min"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  cinterval = "day",
                                  exclude = c(7,0),
                                  cdate = breeding_weather$date_new, bdate = nest_data$first_egg) 
beep(2)


### Max values (only temperature)

breeding_dfe_sw_max <- slidingwin(baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data),  
                                  xvar = list(TempMean=breeding_weather$temp_avg),
                                  type = "absolute", 
                                  refday = c(01,06),
                                  range= c(47,0),
                                  stat = c("max"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  cinterval = "day",
                                  exclude = c(7,0),
                                  cdate = breeding_weather$date_new, bdate = nest_data$first_egg) 


### Breeding results - average

# Temperature
head(breeding_dfe_sw_avg[[1]]$Dataset, 40) ## Mean temp (39-4) (B = -1.84)
head(breeding_dfe_sw_avg[[2]]$Dataset, 40) ## Freeze days (14-8) (B = 24.4)

# Precipitation
head(breeding_dfe_sw_avg[[3]]$Dataset, 40) ## Precip (none better than the null)


### Breeding results - min

# Temperature
head(breeding_dfe_sw_min[[1]]$Dataset, 40) ## Min temp (14-8) (B = -1.20)


### Breeding results - max

# Temperature
head(breeding_dfe_sw_max[[1]]$Dataset, 40) ## Mean temp (20-12) (B = -0.89)



##################
### Randomizations
##################


### Breeding - average 

breeding_dfe_ran_avg <- randwin(repeats = 100, 
                                xvar = list(TempMean=breeding_weather$temp_avg, 
                                            Frozen = breeding_weather$sub_zero),
                                            #RainTotal = breeding_weather$precip.days),
                                cdate = breeding_weather$date_new, 
                                bdate = nest_data$first_egg, 
                                baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data), 
                                type = "absolute",  refday = c(01,06),
                                range = c(47,0), stat = c("mean"), 
                                exclude = c(7,0), cmissing = "method1",
                                func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = breeding_dfe_sw_avg[[1]]$Dataset, datasetrand = breeding_dfe_ran_avg[[1]], 
       metric="AIC", sample.size=382) # 0.04

# Freeze days
pvalue(dataset = breeding_dfe_sw_avg[[2]]$Dataset, datasetrand = breeding_dfe_ran_avg[[2]], 
       metric="AIC", sample.size=382) # 0.04



### Breeding - min 

breeding_dfe_ran_min <- randwin(repeats = 100, 
                                xvar = list(TempMean=breeding_weather$temp_avg),
                                cdate = breeding_weather$date_new, 
                                bdate = nest_data$first_egg, 
                                baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data), 
                                type = "absolute",  refday = c(01,06),
                                range = c(47,0), stat = c("min"), 
                                exclude = c(7,0), cmissing = "method1",
                                func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = breeding_dfe_sw_min[[1]]$Dataset, datasetrand = breeding_dfe_ran_min[[1]], 
       metric="AIC", sample.size=382) # <0.001


### Breeding - max 

breeding_dfe_ran_max <- randwin(repeats = 100, 
                                xvar = list(TempMean=breeding_weather$temp_avg),
                                cdate = breeding_weather$date_new, 
                                bdate = nest_data$first_egg, 
                                baseline = lmer(clutch_initiation ~ 1 + (1|year_f), data = nest_data), 
                                type = "absolute",  refday = c(01,06),
                                range = c(47,0), stat = c("max"), 
                                exclude = c(7,0), cmissing = "method1",
                                func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = breeding_dfe_sw_max[[1]]$Dataset, datasetrand = breeding_dfe_ran_max[[1]], 
       metric="AIC", sample.size=382)  


###########################################################################################
### Extract weather variables from time periods identified by sliding window ##############
###########################################################################################


###############
### Clutch initiation date
###############


############
### 1) Make time periods in ordinal date
############


##########
### Winter 
##########

### Min average temperature (180-159)

nest_data$min_avg_temp_start_w <- (152-180) + 365 # the 365 is necessary to convert to a valid ordinal date
nest_data$min_avg_temp_end_w <- (152-159) + 365


############
### Stopover
############

### Avg minimum temp (72-51)

nest_data$avg_min_temp_start_s <- (152-71)
nest_data$avg_min_temp_end_s <- (152-51)


### Freeze days (66-45)

nest_data$freeze_days_start_s <- (152-66)
nest_data$freeze_days_end_s <- (152-45)


### Precip (88-52)

nest_data$precip_start_s <- (152-88)
nest_data$precip_end_s <- (152-52)


### Minimum average temp (60-39)

nest_data$min_avg_temp_start_s <- (152-60)
nest_data$min_avg_temp_end_s <- (152-39)


### Minimum min temp (79-47)

nest_data$min_min_temp_start_s <- (152-79)
nest_data$min_min_temp_end_s <- (152-47)

### Maximum avg temp (52-31)

nest_data$max_avg_temp_start_s <- (152-52)
nest_data$max_avg_temp_end_s <- (152-31)


################
### Breeding
################

### Mean avg temperature (39-4)

nest_data$avg_temp_start_b <- (152-39)
nest_data$avg_temp_end_b <- (152-4)

### Freeze days (14-8)

nest_data$freeze_days_start_b <- (152-14)
nest_data$freeze_days_end_b <- (152-8)

### Min avg temp

nest_data$min_avg_start_b <- (152-14)
nest_data$min_avg_end_b <- (152-8)

### Max avg temp (20-12)

nest_data$max_avg_start_b <- (152-20)
nest_data$max_avg_end_b <- (152-12)


############
### 2) Convert ordinal date time periods into start and end date-time
############


##########
### Winter 
##########

### Minimum average daily temperature

nest_data$min_avg_temp_start_w_int <- strptime(paste(nest_data$year-1, nest_data$min_avg_temp_start_w), format="%Y %j") 
nest_data$min_avg_temp_start_w_date <- as.POSIXct(nest_data$min_avg_temp_start_w_int)

nest_data$min_avg_temp_end_w_int <- strptime(paste(nest_data$year-1, nest_data$min_avg_temp_end_w), format="%Y %j") 
nest_data$min_avg_temp_end_w_date <- as.POSIXct(nest_data$min_avg_temp_end_w_int)


############
### Stopover 
############

### Avg min temp

nest_data$avg_min_temp_start_s_int <- strptime(paste(nest_data$year, nest_data$avg_min_temp_start_s), format="%Y %j") 
nest_data$avg_min_temp_start_s_date <- as.POSIXct(nest_data$avg_min_temp_start_s_int)

nest_data$avg_min_temp_end_s_int <- strptime(paste(nest_data$year, nest_data$avg_min_temp_end_s), format="%Y %j") 
nest_data$avg_min_temp_end_s_date <- as.POSIXct(nest_data$avg_min_temp_end_s_int)


### Freeze days

nest_data$freeze_days_start_s_int <- strptime(paste(nest_data$year, nest_data$freeze_days_start_s), format="%Y %j") 
nest_data$freeze_days_start_s_date <- as.POSIXct(nest_data$freeze_days_start_s_int)

nest_data$freeze_days_end_s_int <- strptime(paste(nest_data$year, nest_data$freeze_days_end_s), format="%Y %j") 
nest_data$freeze_days_end_s_date <- as.POSIXct(nest_data$freeze_days_end_s_int)


### Precip

nest_data$precip_start_s_int <- strptime(paste(nest_data$year, nest_data$precip_start_s), format="%Y %j") 
nest_data$precip_start_s_date <- as.POSIXct(nest_data$precip_start_s_int)

nest_data$precip_end_s_int <- strptime(paste(nest_data$year, nest_data$precip_end_s), format="%Y %j") 
nest_data$precip_end_s_date <- as.POSIXct(nest_data$precip_end_s_int)


### Minimum daily average

nest_data$min_avg_temp_start_s_int <- strptime(paste(nest_data$year, nest_data$min_avg_temp_start_s), format="%Y %j") 
nest_data$min_avg_temp_start_s_date <- as.POSIXct(nest_data$min_avg_temp_start_s_int)

nest_data$min_avg_temp_end_s_int <- strptime(paste(nest_data$year, nest_data$min_avg_temp_end_s), format="%Y %j") 
nest_data$min_avg_temp_end_s_date <- as.POSIXct(nest_data$min_avg_temp_end_s_int)


### Minimum min temp

nest_data$min_min_temp_start_s_int <- strptime(paste(nest_data$year, nest_data$min_min_temp_start_s), format="%Y %j") 
nest_data$min_min_temp_start_s_date <- as.POSIXct(nest_data$min_min_temp_start_s_int)

nest_data$min_min_temp_end_s_int <- strptime(paste(nest_data$year, nest_data$min_min_temp_end_s), format="%Y %j") 
nest_data$min_min_temp_end_s_date <- as.POSIXct(nest_data$min_min_temp_end_s_int)


### Maximum avg temp

nest_data$max_avg_temp_start_s_int <- strptime(paste(nest_data$year, nest_data$max_avg_temp_start_s), format="%Y %j") 
nest_data$max_avg_temp_start_s_date <- as.POSIXct(nest_data$max_avg_temp_start_s_int)

nest_data$max_avg_temp_end_s_int <- strptime(paste(nest_data$year, nest_data$max_avg_temp_end_s), format="%Y %j") 
nest_data$max_avg_temp_end_s_date <- as.POSIXct(nest_data$max_avg_temp_end_s_int)


################
### Breeding
################

### Average temperature

nest_data$avg_temp_start_b_int <- strptime(paste(nest_data$year, nest_data$avg_temp_start_b), format="%Y %j") 
nest_data$avg_temp_start_b_date <- as.POSIXct(nest_data$avg_temp_start_b_int)

nest_data$avg_temp_end_b_int <- strptime(paste(nest_data$year, nest_data$avg_temp_end_b), format="%Y %j") 
nest_data$avg_temp_end_b_date <- as.POSIXct(nest_data$avg_temp_end_b_int)


### Freeze days

nest_data$freeze_days_start_b_int <- strptime(paste(nest_data$year, nest_data$freeze_days_start_b), format="%Y %j") 
nest_data$freeze_days_start_b_date <- as.POSIXct(nest_data$freeze_days_start_b_int)

nest_data$freeze_days_end_b_int <- strptime(paste(nest_data$year, nest_data$freeze_days_end_b), format="%Y %j") 
nest_data$freeze_days_end_b_date <- as.POSIXct(nest_data$freeze_days_end_b_int)


### Min avg temp

nest_data$min_avg_start_b_int <- strptime(paste(nest_data$year, nest_data$min_avg_start_b), format="%Y %j") 
nest_data$min_avg_start_b_date <- as.POSIXct(nest_data$min_avg_start_b_int)

nest_data$min_avg_end_b_int <- strptime(paste(nest_data$year, nest_data$min_avg_end_b), format="%Y %j") 
nest_data$min_avg_end_b_date <- as.POSIXct(nest_data$min_avg_end_b_int)


### Max avg temp

nest_data$max_avg_start_b_int <- strptime(paste(nest_data$year, nest_data$max_avg_start_b), format="%Y %j") 
nest_data$max_avg_start_b_date <- as.POSIXct(nest_data$max_avg_start_b_int)

nest_data$max_avg_end_b_int <- strptime(paste(nest_data$year, nest_data$max_avg_end_b), format="%Y %j") 
nest_data$max_avg_end_b_date <- as.POSIXct(nest_data$max_avg_end_b_int)



###########################################
### Convert ordinal dates to proper posix format
###########################################

### Convert weather dates to POSIXct or else the following code will not work 
### (i.e., it won't work with same date format as sliding window)

# Winter

winter_weather$j_date <- yday(winter_weather$date)
winter_weather$date.int <- strptime(paste(winter_weather$year, winter_weather$j_date), format="%Y %j") 
winter_weather$date_use <- as.POSIXct(winter_weather$date.int)


# Stopover

stopover_weather$j_date <- yday(stopover_weather$date)
stopover_weather$date.int <- strptime(paste(stopover_weather$year, stopover_weather$j_date), format="%Y %j") 
stopover_weather$date_use <- as.POSIXct(stopover_weather$date.int)


# Breeding

breeding_weather$date.int <- strptime(paste(breeding_weather$year, breeding_weather$jdate), format="%Y %j") 
breeding_weather$date_use <- as.POSIXct(breeding_weather$date.int)


############################################################################
### Create loop to extract weather variables for each year and save as a csv
############################################################################

### Empty lists

min_avg_temp_w_list <- list()

avg_min_temp_s_list <- list()
freeze_days_s_list <- list()
precip_s_list <- list()
min_avg_temp_s_list <- list()
min_min_temp_s_list <- list()
max_avg_temp_s_list <- list()

avg_temp_b_list <- list()
freeze_days_b_list <- list()
min_avg_temp_b_list <- list()
max_avg_temp_b_list <- list()


### Run loop across years

for (i in 1:length(nest_data$nest_id)) {
  
  # Winter
  min_avg_temp_w_list[[i]] <- winter_weather %>% 
                              filter(date_use >= nest_data$min_avg_temp_start_w_date[i] 
                                     & date_use <= nest_data$min_avg_temp_end_w_date[i]) %>% 
    
                              summarise(min_avg_temp_w = min(avg_daily_temp))
  
  
  # Stopover
  avg_min_temp_s_list[[i]] <- stopover_weather %>% 
                              filter(date_use >= nest_data$avg_min_temp_start_s_date[i] 
                                     & date_use <= nest_data$avg_min_temp_end_s_date[i]) %>% 
    
                              summarise(avg_min_temp_s = mean(min_daily_temp))
  
  
  freeze_days_s_list[[i]] <- stopover_weather %>% 
                             filter(date_use >= nest_data$freeze_days_start_s_date[i] 
                            & date_use <= nest_data$freeze_days_end_s_date[i]) %>% 
    
                            summarise(freeze_days_s = mean(sub_zero))
  
  
  precip_s_list[[i]] <- stopover_weather %>% 
                        filter(date_use >= nest_data$precip_start_s_date[i] 
                               & date_use <= nest_data$precip_end_s_date[i]) %>% 
    
                        summarise(precip_s = mean(precip_daily_sum))
  
  min_avg_temp_s_list[[i]] <- stopover_weather %>% 
                              filter(date_use >= nest_data$min_avg_temp_start_s_date[i] 
                                     & date_use <= nest_data$min_avg_temp_end_s_date[i]) %>% 
    
                              summarise(min_avg_temp_s = min(avg_daily_temp))
  
  min_min_temp_s_list[[i]] <- stopover_weather %>% 
                              filter(date_use >= nest_data$min_min_temp_start_s_date[i] 
                                     & date_use <= nest_data$min_min_temp_end_s_date[i]) %>% 
    
                              summarise(min_min_temp_s = min(min_daily_temp))
  
  max_avg_temp_s_list[[i]] <- stopover_weather %>% 
                              filter(date_use >= nest_data$max_avg_temp_start_s_date[i] 
                                     & date_use <= nest_data$max_avg_temp_end_s_date[i]) %>% 
    
                              summarise(max_avg_temp_s = max(avg_daily_temp))
  
  
  # Breeding
  
  avg_temp_b_list[[i]] <- breeding_weather %>% 
                          filter(date_use >= nest_data$avg_temp_start_b_date[i] 
                                 & date_use <= nest_data$avg_temp_end_b_date[i]) %>% 
    
                          summarise(avg_temp_b = mean(temp_avg))
  
  freeze_days_b_list[[i]] <- breeding_weather %>% 
                             filter(date_use >= nest_data$freeze_days_start_b_date[i] 
                                    & date_use <= nest_data$freeze_days_end_b_date[i]) %>% 
    
                             summarise(freeze_days_b = mean(sub_zero))
  
  min_avg_temp_b_list[[i]] <- breeding_weather %>% 
                              filter(date_use >= nest_data$min_avg_start_b_date[i] 
                                     & date_use <= nest_data$min_avg_end_b_date[i]) %>% 
    
                              summarise(min_avg_temp_b = min(temp_avg))
  
  max_avg_temp_b_list[[i]] <- breeding_weather %>% 
                              filter(date_use >= nest_data$max_avg_start_b_date[i] 
                                     & date_use <= nest_data$max_avg_end_b_date[i]) %>% 
    
                              summarise(max_avg_temp_b = max(temp_avg))
  
  
  ### Convert lists to rows
  
  min_avg_temp_w_data <- bind_rows(min_avg_temp_w_list)
  
  avg_min_temp_s_data <- bind_rows(avg_min_temp_s_list)
  freeze_days_s_data <- bind_rows(freeze_days_s_list)
  precip_s_data <- bind_rows(precip_s_list)
  min_avg_temp_s_data <- bind_rows(min_avg_temp_s_list)
  min_min_temp_s_data <- bind_rows(min_min_temp_s_list)
  max_avg_temp_s_data <- bind_rows(max_avg_temp_s_list)
  
  avg_temp_b_data <- bind_rows(avg_temp_b_list)
  freeze_days_b_data <- bind_rows(freeze_days_b_list)
  min_avg_temp_b_data <- bind_rows(min_avg_temp_b_list)
  max_avg_temp_b_data <- bind_rows(max_avg_temp_b_list)
  
  
  ### Convert each variable to a column in a larger dataset
  climate_windows <- bind_cols(min_avg_temp_w_data, avg_min_temp_s_data,
                               freeze_days_s_data, precip_s_data, 
                               min_avg_temp_s_data, min_min_temp_s_data, 
                               max_avg_temp_s_data, avg_temp_b_data, 
                               freeze_days_b_data,  min_avg_temp_b_data,
                               max_avg_temp_b_data)
  
}


### Combine climate windows with nest data

nest_data_new <- nest_data[, c(1:8)] # Extract only columns of interest from original data

nest_data_final <- bind_cols(nest_data_new, climate_windows)


### Save data

#write.csv(nest_data_final, "ci_weather_summarized.csv")



################################################################################
### Nestling development time
################################################################################


#############
### A) Winter
#############

### Average values

winter_dev_sw_avg <- slidingwin(baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data),  
                                xvar = list(TempMean=winter_weather$avg_daily_temp, 
                                            Frozen = winter_weather$sub_zero,
                                            TempMin = winter_weather$min_daily_temp, 
                                            TempMax = winter_weather$max_daily_temp, 
                                            RainTotal = winter_weather$precip_daily_sum),
                                type = "absolute", 
                                refday = c(01,06),
                                range = c(182,31), # December 1 to May 1
                                stat = c("mean"), 
                                cmissing = "method1",
                                func = "lin",
                                exclude = c(21,0), # Excludes windows <7 days at any time.
                                cinterval = "day",
                                cdate = winter_weather$date_new, bdate = fledge_data$first_egg) 


### Minimum values 

winter_dev_sw_min <- slidingwin(baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data),  
                                xvar = list(TempMean=winter_weather$avg_daily_temp, 
                                            TempMin = winter_weather$min_daily_temp), 
                                type = "absolute", 
                                refday = c(01,06),
                                range = c(182,31), # December 1 to May 1
                                stat = c("min"), 
                                cmissing = "method1",
                                func = "lin",
                                exclude = c(21,0), 
                                cinterval = "day",
                                cdate = winter_weather$date_new, bdate = fledge_data$first_egg) 
beep(2)



### Max values 

winter_dev_sw_max <- slidingwin(baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data),  
                                xvar = list(TempMean=winter_weather$avg_daily_temp, 
                                            TempMin = winter_weather$max_daily_temp), 
                                type = "absolute", 
                                refday = c(01,06),
                                range = c(182,31), # December 1 to May 1
                                stat = c("max"), 
                                cmissing = "method1",
                                func = "lin",
                                exclude = c(21,0),
                                cinterval = "day",
                                cdate = winter_weather$date_new, bdate = fledge_data$first_egg) 


### Winter results - average

# Temperature
head(winter_dev_sw_avg[[1]]$Dataset, 40) ## Mean temp (none better than null)
head(winter_dev_sw_avg[[2]]$Dataset, 40) ## Freeze days (65-31) (B = 14.4)
head(winter_dev_sw_avg[[3]]$Dataset, 40) ## Min temp (none better than the null)
head(winter_dev_sw_avg[[4]]$Dataset, 40) ## Max temp (none better than the null)

# Precipitation
head(winter_dev_sw_avg[[5]]$Dataset, 40) ## Precip (none better than the null)


### Winter results - min

# Temperature
head(winter_dev_sw_min[[1]]$Dataset, 40) ## Mean temp (None better than null)
head(winter_dev_sw_min[[2]]$Dataset, 40) ## Min temp (None better than null)


### Winter results - max

# Temperature
head(winter_dev_sw_max[[1]]$Dataset, 40) ## Mean temp (None better than null)
head(winter_dev_sw_max[[2]]$Dataset, 40) ## Max temp (None better than null)



##################
### Randomizations
##################

### Winter - average 

winter_dev_ran_avg <- randwin(repeats = 100, 
                              xvar = list(#TempMean=winter_climate$avg_daily_temp, 
                                          Frozen = winter_weather$sub_zero),
                                          #TempMin = winter_climate$min_daily_temp, 
                                          #TempMax = winter_climate$max_daily_temp, 
                                          #RainTotal = winter_climate$precip_daily_sum),
                              cdate = winter_weather$date_new, 
                              bdate = fledge_data$first_egg, 
                              baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data), 
                              type = "absolute",  refday = c(01,06),
                              range = c(182,31), stat = c("mean"), 
                              exclude = c(21,0), cmissing = "method1",
                              func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########


# Freeze days
pvalue(dataset = winter_dev_sw_avg[[2]]$Dataset, datasetrand = winter_dev_ran_avg[[1]], 
       metric="AIC", sample.size=118) # 0.2



### Winter - min 

winter_dev_ran_min <- randwin(repeats = 100, 
                              xvar = list(TempMean=winter_weather$avg_daily_temp), 
                                          #TempMin = winter_climate$min_daily_temp),
                              cdate = winter_weather$date_new, 
                              bdate = fledge_data$first_egg, 
                              baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data), 
                              type = "absolute",  refday = c(01,06),
                              range = c(182,31), stat = c("min"), 
                              exclude = c(21,0),
                              func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = winter_dev_sw_min[[1]]$Dataset, datasetrand = winter_dev_ran_min[[1]], 
       metric="AIC", sample.size=118) # 0.04



#################
### B) Stopover
#################


### Average values (all weather variables)

stopover_dev_sw_avg <- slidingwin(baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data),  
                                  xvar = list(TempMean=stopover_weather$avg_daily_temp, 
                                              Frozen = stopover_weather$sub_zero,
                                              TempMin = stopover_weather$min_daily_temp, 
                                              TempMax = stopover_weather$max_daily_temp, 
                                              RainTotal = stopover_weather$precip_daily_sum),
                                  type = "absolute", 
                                  refday = c(01,06),
                                  range = c(92,31), # March 1 to May 1
                                  stat = c("mean"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  exclude = c(21,0), 
                                  cinterval = "day",
                                  cdate = stopover_weather$date_new, bdate = fledge_data$first_egg) 
beep(2)


### Minimum values 
### (temperature variables only, as does not make sense for precipitation and freeze days)

stopover_dev_sw_min <- slidingwin(baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data),  
                                  xvar = list(TempMean=stopover_weather$avg_daily_temp, 
                                              TempMin = stopover_weather$min_daily_temp), 
                                  type = "absolute", 
                                  refday = c(01,06),
                                  range = c(92,31), # December 1 to May 1
                                  stat = c("min"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  exclude = c(21,0),
                                  cinterval = "day",
                                  cdate = stopover_weather$date_new, bdate = fledge_data$first_egg) 
beep(2)


### Max values 

stopover_dev_sw_max <- slidingwin(baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data),  
                                  xvar = list(TempMean=stopover_weather$avg_daily_temp, 
                                              TempMax = stopover_weather$max_daily_temp), 
                                  type = "absolute", 
                                  refday = c(01,06),
                                  range = c(92,31), # December 1 to May 1
                                  stat = c("max"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  exclude = c(21,0),
                                  cinterval = "day",
                                  cdate = stopover_weather$date_new, bdate = fledge_data$first_egg) 
beep(2)



### Stopover results - average

# Temperature
head(stopover_dev_sw_avg[[1]]$Dataset, 20) ## Mean temp (none better than null)
head(stopover_dev_sw_avg[[2]]$Dataset, 20) ## Freeze days (66-40) (B = 4.90)
head(stopover_dev_sw_avg[[3]]$Dataset, 20) ## Min temp (none better than the null)
head(stopover_dev_sw_avg[[4]]$Dataset, 20) ## Max temp (none better than the null)

# Precipitation
head(stopover_dev_sw_avg[[5]]$Dataset, 20) ## Precip (none better than the null)


### Winter results - min

# Temperature
head(stopover_dev_sw_min[[1]]$Dataset, 40) ## Mean temp (none better than the null)
head(stopover_dev_sw_min[[2]]$Dataset, 40) ## Min temp (65-44) (B = -0.24)


### Winter results - max

# Temperature
head(stopover_dev_sw_max[[1]]$Dataset, 40) ## Mean temp (none better than the null)
head(stopover_dev_sw_max[[2]]$Dataset, 40) ## Max temp (none better than null)



##################
### Randomizations
##################

### Stopover - average 

stopover_dev_ran_avg <- randwin(repeats = 100, 
                                xvar = list(#TempMean=stopover_climate$avg_daily_temp, 
                                            Frozen = stopover_weather$sub_zero),
                                            #TempMin = stopover_climate$min_daily_temp, 
                                            #TempMax = stopover_climate$max_daily_temp, 
                                            #RainTotal = stopover_climate$precip_daily_sum),
                                cdate = stopover_weather$date_new, 
                                bdate = fledge_data$first_egg, 
                                baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data), 
                                type = "absolute",  refday = c(01,06),
                                range = c(92,31), stat = c("mean"), 
                                exclude = c(21,0),
                                func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Freeze days
pvalue(dataset = stopover_dev_sw_avg[[2]]$Dataset, datasetrand = stopover_dev_ran_avg[[1]], 
       metric="AIC", sample.size=118) # 0.04


### Stopover - min 

stopover_dev_ran_min <- randwin(repeats = 100, 
                                xvar = list(#TempMean= stopover_climate$avg_daily_temp, 
                                            TempMin = stopover_weather$min_daily_temp), 
                                            #TempMax = stopover_climate$max_daily_temp),
                                cdate = stopover_weather$date_new, 
                                bdate = fledge_data$first_egg, 
                                baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data), 
                                type = "absolute",  refday = c(01,06),
                                range = c(92,31), stat = c("min"), 
                                exclude = c(21,0),
                                func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Min temp
pvalue(dataset = stopover_dev_sw_min[[2]]$Dataset, datasetrand = stopover_dev_ran_min[[1]], 
       metric="AIC", sample.size=118) # 0.02


###########################
### C) Breeding site
###########################


### Average values (all weather variables)

breeding_dev_sw_avg <- slidingwin(baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data),  
                                  xvar = list(TempMean=breeding_weather$temp_avg, 
                                              Frozen = breeding_weather$sub_zero,
                                              RainTotal = breeding_weather$precip_days),
                                  type = "absolute", 
                                  refday = c(01,06),
                                  range= c(47,0), 
                                  stat = c("mean"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  cinterval = "day",
                                  exclude = c(7,0),
                                  cdate = breeding_weather$date_new, bdate = fledge_data$first_egg) 
beep(2)



### Minimum values (only temperature)

breeding_dev_sw_min <- slidingwin(baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data),  
                                  xvar = list(TempMean=breeding_weather$temp_avg),
                                  type = "absolute", 
                                  refday = c(01,06),
                                  range= c(47,0),
                                  stat = c("min"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  cinterval = "day",
                                  exclude = c(7,0), 
                                  cdate = breeding_weather$date_new, bdate = fledge_data$first_egg) 
beep(2)


### Max values (only temperature)

breeding_dev_sw_max <- slidingwin(baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data),  
                                  xvar = list(TempMean=breeding_weather$temp_avg),
                                  type = "absolute", 
                                  refday = c(01,06),
                                  range= c(47,0),
                                  stat = c("max"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  cinterval = "day",
                                  exclude = c(7,0),
                                  cdate = breeding_weather$date_new, bdate = fledge_data$first_egg)


### Breeding results - average

# Temperature
head(breeding_dev_sw_avg[[1]]$Dataset, 40) ## Mean temp (29-21) (-0.15)
head(breeding_dev_sw_avg[[2]]$Dataset, 40) ## Freeze days (none better than null)

# Precipitation
head(breeding_dev_sw_avg[[3]]$Dataset, 40) ## Precip (25-12) (B = -1.65)


### Breeding results - min

# Temperature
head(breeding_dev_sw_min[[1]]$Dataset, 40) ## Min temp (none better than the null)


### Breeding results - max

# Temperature
head(breeding_dev_sw_max[[1]]$Dataset, 40) ## Max temp (11-1) (-0.10)


##################
### Randomizations
##################


### Breeding - average 

breeding_dev_ran_avg <- randwin(repeats = 100, 
                                xvar = list(TempMean=breeding_weather$temp_avg), 
                                            #Frozen = breeding_weather$sub_zero,
                                            #RainTotal = breeding_weather$precip_days),
                                cdate = breeding_weather$date_new, 
                                bdate = fledge_data$first_egg, 
                                baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data), 
                                type = "absolute",  refday = c(01,06),
                                range = c(47,0), stat = c("mean"), 
                                cmissing = "method1",
                                exclude = c(7,0),
                                func = "lin", cinterval = "day")
beep(8)


#########
### Compare randomizations to original sliding window
#########

# Mean temp
pvalue(dataset = breeding_dev_sw_avg[[1]]$Dataset, datasetrand = breeding_dev_ran_avg[[1]], 
       metric="AIC", sample.size=118) # 0.24

# Precip
pvalue(dataset = breeding_dev_sw_avg[[3]]$Dataset, datasetrand = breeding_dev_ran_avg[[2]], 
       metric="AIC", sample.size=118) # 0.08



### Breeding - min 

# Note - not necessary because already selected out with null window results.

### Breeding - max

breeding_dev_ran_max <- randwin(repeats = 100, 
                                xvar = list(TempMean=breeding_weather$temp_avg), 
                                cdate = breeding_weather$date_new, 
                                bdate = fledge_data$first_egg, 
                                baseline = lmer(age_at_fledge ~ 1 + (1|year_f), data = fledge_data), 
                                type = "absolute",  refday = c(01,06),
                                range = c(47,0), stat = c("max"), 
                                cmissing = "method1",
                                exclude = c(7,0),
                                func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Max temp
pvalue(dataset = breeding_dev_sw_max[[1]]$Dataset, datasetrand = breeding_dev_ran_max[[1]], 
       metric="AIC", sample.size=118) # 0.36



###########################################################################################
### Extract weather variables from time periods identified by sliding window ##############
###########################################################################################


############
### 1) Make time periods in ordinal date
############


##########
### Winter 
##########

### Minimum average daily temperature (182-161)

fledge_data$min_avg_temp_start_w <- (152-182) + 365 
fledge_data$min_avg_temp_end_w <- (152-161) + 365


############
### Stopover
############


### Freeze days (66-40)

fledge_data$freeze_days_start_s <- (152-66)
fledge_data$freeze_days_end_s <- (152-40)


################
### Breeding
################

### Precipitation (25-12)

fledge_data$precip_start_b <- (152-25)
fledge_data$precip_end_b <- (152-12)


############
### 2) Convert ordinal date time periods into start and end date-time
############


##########
### Winter 
##########

### Minimum average daily temperature

fledge_data$min_avg_temp_start_w_int <- strptime(paste(fledge_data$year-1, fledge_data$min_avg_temp_start_w), format="%Y %j") 
fledge_data$min_avg_temp_start_w_date <- as.POSIXct(fledge_data$min_avg_temp_start_w_int)

fledge_data$min_avg_temp_end_w_int <- strptime(paste(fledge_data$year-1, fledge_data$min_avg_temp_end_w), format="%Y %j") 
fledge_data$min_avg_temp_end_w_date <- as.POSIXct(fledge_data$min_avg_temp_end_w_int)


############
### Stopover 
############

### Freeze days

fledge_data$freeze_days_start_s_int <- strptime(paste(fledge_data$year, fledge_data$freeze_days_start_s), format="%Y %j") 
fledge_data$freeze_days_start_s_date <- as.POSIXct(fledge_data$freeze_days_start_s_int)

fledge_data$freeze_days_end_s_int <- strptime(paste(fledge_data$year, fledge_data$freeze_days_end_s), format="%Y %j") 
fledge_data$freeze_days_end_s_date <- as.POSIXct(fledge_data$freeze_days_end_s_int)


################
### Breeding
################


### Precipitation

fledge_data$precip_start_b_int <- strptime(paste(fledge_data$year, fledge_data$precip_start_b), format="%Y %j") 
fledge_data$precip_start_b_date <- as.POSIXct(fledge_data$precip_start_b_int)

fledge_data$precip_end_b_int <- strptime(paste(fledge_data$year, fledge_data$precip_end_b), format="%Y %j") 
fledge_data$precip_end_b_date <- as.POSIXct(fledge_data$precip_end_b_int)


############################################################################
### Create loop to extract weather variables for each year and save as a csv
############################################################################

### Empty lists

min_avg_temp_w_list <- list()
freeze_days_s_list <- list()
precip_b_list <- list()


### Run loop across years

for (i in 1:length(fledge_data$nest_id)) {
  
  # Winter
  min_avg_temp_w_list[[i]] <- winter_weather %>% 
                              filter(date_use >= fledge_data$min_avg_temp_start_w_date[i] 
                                     & date_use <= fledge_data$min_avg_temp_end_w_date[i]) %>% 
    
                              summarise(min_avg_temp_w = min(avg_daily_temp))
  
  
  # Stopover
  freeze_days_s_list[[i]] <- stopover_weather %>% 
                             filter(date_use >= fledge_data$freeze_days_start_s_date[i] 
                                    & date_use <= fledge_data$freeze_days_end_s_date[i]) %>% 
    
                             summarise(freeze_days_s = mean(sub_zero))
  
  # Breeding
  
  precip_b_list[[i]] <- breeding_weather %>% 
                        filter(date_use >= fledge_data$precip_start_b_date[i] 
                               & date_use <= fledge_data$precip_end_b_date[i]) %>% 
    
                        summarise(precip_days_b = sum(precip_days))
  
  ### Convert lists to rows
  
  min_avg_temp_w_data <- bind_rows(min_avg_temp_w_list)
  
  freeze_days_s_data <- bind_rows(freeze_days_s_list)
  
  precip_b_data <- bind_rows(precip_b_list)
  
 
  ### Convert each variable to a column in a larger dataset
  dev_climate_windows <- bind_cols(min_avg_temp_w_data, freeze_days_s_data,
                                   precip_b_data)
  
}


### Combine climate windows with nest data

head(fledge_data)

fledge_data_new <- fledge_data[, c(1:8)] 

dev_final <- bind_cols(fledge_data_new, dev_climate_windows)


### Save data

#write.csv(dev_final, "dev_weather_summarized.csv")


############################ End of script #####################################

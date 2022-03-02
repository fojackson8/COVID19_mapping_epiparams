#Deconvolution and rt calculation
#set WD
setwd("~/rhys/pystan")

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Main functions to run script
files.sources = list.files(path = "./main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}

#Load packages

library(tidyverse)
library(safejoin)
library(data.table)
library(lubridate)
library(dplyr)
library(zoo) 
library(countrycode) 
library(purrr) 
library(readr) 
library(stringr)   
library(tidyr) 
library(usdata)
library(seqinr)
library(ggplot2)
library(xts)
library(MASS)
library(ggfortify)
library(EpiEstim)
library(grid)

#Load in argentina

argentina = as.data.frame(read_csv("data/argentina_linelist_geocoded_adm1_20220131.csv"))

#Calculate number of cases per day

argentina_variables <-  argentina %>% 
  group_by(events.confirmed.date) %>%
  summarise(cases = n())

#fill in missing gaps in dates
 
argentina_variables <- as.data.frame(argentina_variables) %>%
  mutate(events.confirmed.date = as.Date(events.confirmed.date)) %>%
  complete(events.confirmed.date = seq.Date(min(events.confirmed.date), max(events.confirmed.date), by="day"))

argentina_variables[is.na(argentina_variables)] <- 0

argentina_variables <- filter(argentina_variables, events.confirmed.date >= "2020-03-03" & events.confirmed.date <= "2021-03-31")

#filer for epoch

argentina_epoch1 <- filter (argentina_variables,events.confirmed.date >= "2020-03-03" & events.confirmed.date <= "2020-06-30")
argentina_epoch2 <- filter(argentina_variables, events.confirmed.date >= "2020-07-01" & events.confirmed.date <= "2020-11-30")
argentina_epoch3 <- filter(argentina_variables, events.confirmed.date >= "2020-12-01" & events.confirmed.date <= "2021-03-31")

#number of days in each epoch

argentina_epoch1$length <- 1:120
argentina_epoch2$length <- 120:272
argentina_epoch3$length <- 272:392

#load distribution

onset_to_diagnosis_epoch_1 = fread('fitting_outputs/onset-to-diagnosis_epoch_1_ar-samples-lognormal.csv')
onset_to_diagnosis_epoch_2 = fread('fitting_outputs/onset-to-diagnosis_epoch_2_ar-samples-lognormal.csv')
onset_to_diagnosis_epoch_3 = fread('fitting_outputs/onset-to-diagnosis_epoch_3_ar-samples-lognormal.csv')

#calculate mean delay

onset_to_diagnosis_epoch_1 <- exp((onset_to_diagnosis_epoch_1[,1:20]+0.5)*(onset_to_diagnosis_epoch_1[,21:40]^2)) 
onset_to_diagnosis_epoch_2 <- exp((onset_to_diagnosis_epoch_2[,1:20]+0.5)*(onset_to_diagnosis_epoch_2[,21:40]^2)) 
onset_to_diagnosis_epoch_3 <- exp((onset_to_diagnosis_epoch_3[,1:20]+0.5)*(onset_to_diagnosis_epoch_3[,21:40]^2)) 

#unlist tibble

onset_to_diagnosis_epoch_1 <- as.numeric(unlist(onset_to_diagnosis_epoch_1))
onset_to_diagnosis_epoch_2 <- as.numeric(unlist(onset_to_diagnosis_epoch_2))
onset_to_diagnosis_epoch_3 <- as.numeric(unlist(onset_to_diagnosis_epoch_3))


#function for converting delay 

get_delay_multinom <- function(r_delay){
    tibble(delay =  r_delay %>% round) %>%
      group_by(delay) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      complete(delay = 0:max(delay), fill = list(n = 0)) %>%
      mutate(n = n/sum(n)) %>%
      pull(n)
  }

## calculate deconvolution of observed cases

ar_cases_epoch1 <- get_RL(observed = argentina_epoch1$cases, 
               times = argentina_epoch1$length, 
       p_delay = get_delay_multinom(onset_to_diagnosis_epoch_1),
               out_col_name = 'case_deconvolved',
       right_censor = FALSE)

ar_cases_epoch1 <- head(ar_cases_epoch1,-1)

ar_cases_epoch2 <- get_RL(observed = argentina_epoch2$cases, 
                          times = argentina_epoch2$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_2),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

ar_cases_epoch2 <- head(ar_cases_epoch2,-1)


ar_cases_epoch3 <- get_RL(observed = argentina_epoch3$cases, 
                          times = argentina_epoch3$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_3),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

ar_cases_epoch3 <- head(ar_cases_epoch3,-1)

ar_cases_deconvolved <- na.omit(round(rbind(ar_cases_epoch1,ar_cases_epoch2,ar_cases_epoch3)))


#plot the two incidences against each other 
#to avoid daily fluctuations will use a 5 day average 

ar_cases_deconvolved_average <- as.data.frame(as.numeric(rollapply(ar_cases_deconvolved[,2],7,fill=NA,partial = TRUE,FUN=mean)))

ar_cases_deconvolved_average$date <- seq(as.Date('2020-03-03'), as.Date('2021-03-28'), by = "1 days")
ar_cases_deconvolved_average$group <- 'Deconvolved case counts'
colnames(ar_cases_deconvolved_average) <- c('case','date','group')

ar_cases <- as.data.frame(as.numeric(rollapply(argentina_variables[,2],7,fill=NA,partial = TRUE,FUN=mean)))
ar_cases$date <- seq(as.Date('2020-03-03'), as.Date('2021-03-31'), by = "1 days")
ar_cases$group <- 'Raw case counts'
colnames(ar_cases) <- c('case','date','group')


combined <- rbind(ar_cases_deconvolved_average,ar_cases)

#plot

country <- grobTree(textGrob("Argentina", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))

ar_incidence <- ggplot(combined, aes(x = date, y=case, color = factor(group), fill = factor(group))) + 
  geom_line() +
  labs(x= "Time",y= 'Incidence', fill = "Method") +
  theme_bw() +
  scale_x_date (date_breaks = "3 month", date_labels = "%Y %b %d") +
  geom_vline(xintercept = as.Date("2020-06-30"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-11-30"), color = "black", linetype = "dashed") +
  ylim(0,27000) + 
  scale_color_manual(values = c("royalblue", "deeppink"),name = "Methods") + 
  theme(text = element_text(size=12)) +
  annotation_custom(country) 
  

ar_incidence

#Calculate rt

#format cases

#raw case counts

argentina_cases_rt <-  argentina %>% 
  group_by(date =  as.Date(events.confirmed.date)) %>% 
  summarise(number_of_cases = n())

argentina_cases_rt <- pad(argentina_cases_rt)

argentina_cases_rt[is.na(argentina_cases_rt)] <- 0

argentina_cases_rt[, 2] <- cumsum(argentina_cases_rt[, 2])

#Deconvolved case counts 

ar_cases_deconvolved_rt <- ar_cases_deconvolved

ar_cases_deconvolved_rt$date <- seq(as.Date('2020-03-03'), as.Date('2021-03-28'), by = "1 days")

ar_cases_deconvolved_rt <- pad(ar_cases_deconvolved_rt)

ar_cases_deconvolved_rt[is.na(ar_cases_deconvolved_rt)] <- 0

ar_cases_deconvolved_rt[, 2] <- cumsum(ar_cases_deconvolved_rt[, 2])

#Merge together
case_data <- left_join(argentina_cases_rt,ar_cases_deconvolved_rt, by = c("date"))

case_data[is.na(case_data)] <- 0

case_data <- filter(case_data, date >= "2020-03-03" & date <= "2021-03-27")

colnames(case_data) <- c("Date", "raw","Day","deconvoluted")

case_data <- case_data[,c(3,1,2,4)]
case_data$Day <- as.integer(case_data$Day)
summary(ar_cases_deconvolved_rt$case_deconvolved)

##LA countries:
start <- 17 # start = 17 is 19 March on both cumulative and daily incidence: this is the data published by WHO on 1 Mar. 
# To compute the incidence on day 52 (29 Feb), the file is read from day 51 (the "start-1" appearing in lines 127 and 151)
# note that start cannot be one as needs to compute the difference in culumative incidence 
Lag <- 0 # Use data till end of dataset

## Open data files
datD <- as.data.frame(case_data)
class(datD)

head(datD)
ncountries = ncol(datD)-2

# Default methods as used in the paper (GCV.Cp with day-of-the-week effect, WD) 
# Further possibilities available for FitMethod, FixedEff
FitMethod <-c('GCV.Cp') # c('GCV.Cp','ML') # Explore fit with different options for the GAM
FixedEff <-c('WD') # c('WE', 'WD', 'None') # Explore alternative fit with weekend effect (WE), day-of-the-week effect (WD), or nothing
colval <-c('black', 'red', 'green')
clist <- c(3,4) # List of the columns for the 5 countries in Figure 1(a)
clistIT <- c(1:3,5) # Relevant columns from Italian data for Figure 1(b)
thin <- 2 # Width of thin line
thick <- 3 # Width of thick line

print('Generating figures for the Electronic Supplementary Material:')

rt <- list()

#change number of interations with double_time function

for(k in 1:length(FitMethod)){
  for(j in 1:length(FixedEff)){
    pdf(file=paste('SuppWHO_',FitMethod[k],'_', FixedEff[j],'_incid.pdf', sep=''))
    for(i in 1:ncountries) {
      a <- DoubleTime(diff(datD[(start-1):(dim(datD)[1]-Lag),2+i]),datD$Day[start:(dim(datD)[1]-Lag)], meth=FitMethod[k], FE=FixedEff[j], plt=T, figtitle=colnames(datD)[2+i])
      print(a)
      rt[[i]] = a 
    }
    dev.off()
  }
}

case_rt <- rt[[1]]
case_rt$dates <- seq(as.Date('2020-03-19'), as.Date('2021-03-28'), by = "1 days")
case_rt$group <- 'Raw case counts'
deconvoluted_rt <- rt[[2]]
deconvoluted_rt$dates <- seq(as.Date('2020-03-19'), as.Date('2021-03-28'), by = "1 days")
deconvoluted_rt$group <- 'Deconvolved case counts'

combined_rt <- rbind(case_rt,deconvoluted_rt)

dev.new()

ar_rt <- ggplot(combined_rt, aes(x = dates, y=sdt, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=sdtlow, ymax=sdtup, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Method") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(r[t]), fill = "Method") +
  geom_vline(xintercept = as.Date("2020-06-30"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-11-30"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "3 month", date_labels = "%Y %b %d") +
  ylim(-0.05,0.2) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  theme(text = element_text(size=12)) +
  annotation_custom(country)

ar_rt

#Load in Brazil

brazil = as.data.frame(fread("data/brazil_linelist_geocoded_adm1_20220131.csv", drop = c("V1", "_id", 
                                                                                         "caseReference.additionalSources", "caseReference.sourceEntryId", "caseReference.sourceId", "caseReference.sourceUrl","caseReference.uploadIds", "events.firstClinicalConsultation.date", 
                                                                                         "genomeSequences","genomeSequences.repositoryUrl","genomeSequences.sampleCollectionDate","location.place","pathogens",
                                                                                         "preexistingConditions.hasPreexistingConditions","preexistingConditions.values","revisionMetadata.creationMetadata.date","revisionMetadata.creationMetadata.notes","revisionMetadata.editMetadata.date","revisionMetadata.editMetadata.notes","revisionMetadata.revisionNumber","transmission.linkedCaseIds",
                                                                                         "transmission.places", "transmission.routes","travelHistory.travel.dateRange.end","travelHistory.travel.dateRange.start","travelHistory.travel.location.administrativeAreaLevel1","travelHistory.travel.location.administrativeAreaLevel2","travelHistory.travel.location.administrativeAreaLevel3","travelHistory.travel.location.country","travelHistory.travel.location.geometry.coordinates",
                                                                                         "travelHistory.travel.location.geoResolution", "travelHistory.travel.location.name","travelHistory.travel.location.place","travelHistory.travel.methods","travelHistory.travel.purpose","travelHistory.traveledPrior30Days","vaccines.0.batch","vaccines.0.date","vaccines.0.name",
                                                                                         "vaccines.0.sideEffects", "vaccines.1.batch","vaccines.1.date","vaccines.1.name","vaccines.1.sideEffects","vaccines.2.batch","vaccines.2.date","vaccines.2.name","vaccines.2.sideEffects",
                                                                                         "vaccines.3.batch", "vaccines.3.date", "vaccines.3.name","vaccines.3.sideEffects",
                                                                                         "caseReference.verificationStatus", "demographics.ageRange.end","demographics.ageRange.start","demographics.ethnicity","demographics.gender","demographics.nationalities","demographics.occupation","genomeSequences.sequenceId","genomeSequences.sequenceLength",
                                                                                         "genomeSequences.sequenceName")))

#Calculate number of cases per day

brazil_variables <-  brazil %>% 
  group_by(events.confirmed.date) %>%
  summarise(cases = n())

#fill in missing gaps in dates

brazil_variables <- as.data.frame(brazil_variables) %>%
  mutate(events.confirmed.date = as.Date(events.confirmed.date)) %>%
  complete(events.confirmed.date = seq.Date(min(events.confirmed.date), max(events.confirmed.date), by="day"))

brazil_variables[is.na(brazil_variables)] <- 0

brazil_variables <- filter(brazil_variables, events.confirmed.date >= "2020-02-25" & events.confirmed.date <= "2021-03-31")

#filer for epoch

brazil_epoch1 <- filter (brazil_variables,events.confirmed.date >= "2020-02-25" & events.confirmed.date <= "2020-06-30")
brazil_epoch2 <- filter(brazil_variables, events.confirmed.date >= "2020-07-01" & events.confirmed.date <= "2020-11-30")
brazil_epoch3 <- filter(brazil_variables, events.confirmed.date >= "2020-12-01" & events.confirmed.date <= "2021-03-31")

#number of days in each epoch

brazil_epoch1$length <- 1:127
brazil_epoch2$length <- 127:279
brazil_epoch3$length <- 279:399

#load distribution

onset_to_diagnosis_epoch_1 = fread('fitting_outputs/onset-to-diagnosis_epoch_1_br-samples-lognormal.csv')
onset_to_diagnosis_epoch_2 = fread('fitting_outputs/onset-to-diagnosis_epoch_2_br-samples-lognormal.csv')
onset_to_diagnosis_epoch_3 = fread('fitting_outputs/onset-to-diagnosis_epoch_3_br-samples-lognormal.csv')

#calculate mean delay

onset_to_diagnosis_epoch_1 <- exp((onset_to_diagnosis_epoch_1[,1:27]+0.5)*(onset_to_diagnosis_epoch_1[,28:54]^2)) 
onset_to_diagnosis_epoch_2 <- exp((onset_to_diagnosis_epoch_2[,1:27]+0.5)*(onset_to_diagnosis_epoch_2[,28:54]^2)) 
onset_to_diagnosis_epoch_3 <- exp((onset_to_diagnosis_epoch_3[,1:27]+0.5)*(onset_to_diagnosis_epoch_3[,28:54]^2)) 

#unlist tibble

onset_to_diagnosis_epoch_1 <- as.numeric(unlist(onset_to_diagnosis_epoch_1))
onset_to_diagnosis_epoch_2 <- as.numeric(unlist(onset_to_diagnosis_epoch_2))
onset_to_diagnosis_epoch_3 <- as.numeric(unlist(onset_to_diagnosis_epoch_3))

## calculate deconvolution of observed cases

br_cases_epoch1 <- get_RL(observed = brazil_epoch1$cases, 
                          times = brazil_epoch1$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_1),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

br_cases_epoch1 <- head(br_cases_epoch1,-1)

br_cases_epoch2 <- get_RL(observed = brazil_epoch2$cases, 
                          times = brazil_epoch2$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_2),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

br_cases_epoch2 <- head(br_cases_epoch2,-1)


br_cases_epoch3 <- get_RL(observed = brazil_epoch3$cases, 
                          times = brazil_epoch3$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_3),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

br_cases_epoch3 <- head(br_cases_epoch3,-1)

br_cases_deconvolved <- na.omit(round(rbind(br_cases_epoch1,br_cases_epoch2,br_cases_epoch3)))

#plot the two incidences against each other 
#to avoid daily fluctuations will use a 5 day average 

br_cases_deconvolved_average <- as.data.frame(as.numeric(rollapply(br_cases_deconvolved[,2],7,fill=NA,partial = TRUE,FUN=mean)))

br_cases_deconvolved_average$date <- seq(as.Date('2020-02-25'), as.Date('2021-03-27'), by = "1 days")
br_cases_deconvolved_average$group <- 'Deconvolved case counts'
colnames(br_cases_deconvolved_average) <- c('case','date','group')

br_cases <- as.data.frame(as.numeric(rollapply(brazil_variables[,2],7,fill=NA,partial = TRUE,FUN=mean)))
br_cases$date <- seq(as.Date('2020-02-25'), as.Date('2021-03-31'), by = "1 days")
br_cases$group <- 'Raw case counts'
colnames(br_cases) <- c('case','date','group')


combined <- rbind(br_cases_deconvolved_average,br_cases)

#plot

country <- grobTree(textGrob("Brazil", x=0.05,  y=0.97, hjust=0,
                             gp=gpar(col="Black", fontsize=12, fontface="bold")))

br_incidence <- ggplot(combined, aes(x = date, y=case, color = factor(group), fill = factor(group))) + 
  geom_line() +
  labs(x= "Time",y= 'Incidence', fill = "Method") +
  theme_bw() +
  scale_x_date (date_breaks = "3 month", date_labels = "%Y %b %d") +
  geom_vline(xintercept = as.Date("2020-06-30"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-11-30"), color = "black", linetype = "dashed") +
  ylim(0,27000) + 
  scale_color_manual(values = c("royalblue", "deeppink"),name = "Methods") + 
  theme(text = element_text(size=12)) +
  annotation_custom(country) 
  

br_incidence

#Calculate rt

#format cases

#raw case counts

brazil_cases_rt <-  brazil %>% 
  group_by(date =  as.Date(events.confirmed.date)) %>% 
  summarise(number_of_cases = n())

brazil_cases_rt <- pad(brazil_cases_rt)

brazil_cases_rt[is.na(brazil_cases_rt)] <- 0

brazil_cases_rt[, 2] <- cumsum(brazil_cases_rt[, 2])

#Deconvolved case counts 

br_cases_deconvolved_rt <- br_cases_deconvolved

br_cases_deconvolved_rt$date <- seq(as.Date('2020-02-25'), as.Date('2021-03-27'), by = "1 days")

br_cases_deconvolved_rt <- pad(br_cases_deconvolved_rt)

br_cases_deconvolved_rt[is.na(br_cases_deconvolved_rt)] <- 0

br_cases_deconvolved_rt[, 2] <- cumsum(br_cases_deconvolved_rt[, 2])

#Merge together
case_data <- left_join(brazil_cases_rt,br_cases_deconvolved_rt, by = c("date"))

case_data[is.na(case_data)] <- 0

case_data <- filter(case_data, date >= "2020-02-25" & date <= "2021-03-27")

colnames(case_data) <- c("Date", "raw","Day","deconvoluted")

case_data <- case_data[,c(3,1,2,4)]
case_data$Day <- as.integer(case_data$Day)

##LA countries:
start <- 25 # start = 25 is 19 March on both cumulative and daily incidence: this is the data published by WHO on 1 Mar. 
# To compute the incidence on day 52 (29 Feb), the file is read from day 51 (the "start-1" appearing in lines 127 and 151)
# note that start cannot be one as needs to compute the difference in culumative incidence 
Lag <- 0 # Use data till end of dataset

## Open data files
datD <- as.data.frame(case_data)
class(datD)

head(datD)
ncountries = ncol(datD)-2

# Default methods as used in the paper (GCV.Cp with day-of-the-week effect, WD) 
# Further possibilities available for FitMethod, FixedEff
FitMethod <-c('GCV.Cp') # c('GCV.Cp','ML') # Explore fit with different options for the GAM
FixedEff <-c('WD') # c('WE', 'WD', 'None') # Explore alternative fit with weekend effect (WE), day-of-the-week effect (WD), or nothing
colval <-c('black', 'red', 'green')
clist <- c(3,4) # List of the columns for the 5 countries in Figure 1(a)
clistIT <- c(1:3,5) # Relevant columns from Italian data for Figure 1(b)
thin <- 2 # Width of thin line
thick <- 3 # Width of thick line

print('Generating figures for the Electronic Supplementary Material:')

rt <- list()

#change number of interations with double_time function

for(k in 1:length(FitMethod)){
  for(j in 1:length(FixedEff)){
    pdf(file=paste('SuppWHO_',FitMethod[k],'_', FixedEff[j],'_incid.pdf', sep=''))
    for(i in 1:ncountries) {
      a <- DoubleTime(diff(datD[(start-1):(dim(datD)[1]-Lag),2+i]),datD$Day[start:(dim(datD)[1]-Lag)], meth=FitMethod[k], FE=FixedEff[j], plt=T, figtitle=colnames(datD)[2+i])
      print(a)
      rt[[i]] = a 
    }
    dev.off()
  }
}

case_rt <- rt[[1]]
case_rt$dates <- seq(as.Date('2020-03-19'), as.Date('2021-03-28'), by = "1 days")
case_rt$group <- 'Raw case counts'
deconvoluted_rt <- rt[[2]]
deconvoluted_rt$dates <- seq(as.Date('2020-03-19'), as.Date('2021-03-28'), by = "1 days")
deconvoluted_rt$group <- 'Deconvolved case counts'

combined_rt <- rbind(case_rt,deconvoluted_rt)

dev.new()

br_rt <- ggplot(combined_rt, aes(x = dates, y=sdt, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=sdtlow, ymax=sdtup, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Method") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(r[t]), fill = "Method") +
  theme_bw() +
  scale_x_date (date_breaks = "3 month", date_labels = "%Y %b %d") +
  ylim(-0.05,0.2) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  theme(text = element_text(size=12)) +
  annotation_custom(country) 

br_rt

#Load in Mexico

mexico = as.data.frame(fread("data/mexico_linelist_geocoded_adm1_20220131.csv", drop = c("V1", "_id", 
                                                                                         "caseReference.additionalSources", "caseReference.sourceEntryId", "caseReference.sourceId", "caseReference.sourceUrl","caseReference.uploadIds", "events.firstClinicalConsultation.date", 
                                                                                         "genomeSequences","genomeSequences.repositoryUrl","genomeSequences.sampleCollectionDate","location.place","pathogens",
                                                                                         "preexistingConditions.hasPreexistingConditions","preexistingConditions.values","revisionMetadata.creationMetadata.date","revisionMetadata.creationMetadata.notes","revisionMetadata.editMetadata.date","revisionMetadata.editMetadata.notes","revisionMetadata.revisionNumber","transmission.linkedCaseIds",
                                                                                         "transmission.places", "transmission.routes","travelHistory.travel.dateRange.end","travelHistory.travel.dateRange.start","travelHistory.travel.location.administrativeAreaLevel1","travelHistory.travel.location.administrativeAreaLevel2","travelHistory.travel.location.administrativeAreaLevel3","travelHistory.travel.location.country","travelHistory.travel.location.geometry.coordinates",
                                                                                         "travelHistory.travel.location.geoResolution", "travelHistory.travel.location.name","travelHistory.travel.location.place","travelHistory.travel.methods","travelHistory.travel.purpose","travelHistory.traveledPrior30Days","vaccines.0.batch","vaccines.0.date","vaccines.0.name",
                                                                                         "vaccines.0.sideEffects", "vaccines.1.batch","vaccines.1.date","vaccines.1.name","vaccines.1.sideEffects","vaccines.2.batch","vaccines.2.date","vaccines.2.name","vaccines.2.sideEffects",
                                                                                         "vaccines.3.batch", "vaccines.3.date", "vaccines.3.name","vaccines.3.sideEffects",
                                                                                         "caseReference.verificationStatus", "demographics.ageRange.end","demographics.ageRange.start","demographics.ethnicity","demographics.gender","demographics.nationalities","demographics.occupation","genomeSequences.sequenceId","genomeSequences.sequenceLength",
                                                                                         "genomeSequences.sequenceName")))

#Calculate number of cases per day

mexico_variables <-  mexico %>% 
  group_by(events.confirmed.date) %>%
  summarise(cases = n())

#fill in missing gaps in dates

mexico_variables <- as.data.frame(mexico_variables) %>%
  mutate(events.confirmed.date = as.Date(events.confirmed.date)) %>%
  complete(events.confirmed.date = seq.Date(min(events.confirmed.date), max(events.confirmed.date), by="day"))

mexico_variables[is.na(mexico_variables)] <- 0

mexico_variables <- filter(mexico_variables, events.confirmed.date >= "2020-02-28" & events.confirmed.date <= "2021-03-31")

#filer for epoch

mexico_epoch1 <- filter (mexico_variables,events.confirmed.date >= "2020-02-28" & events.confirmed.date <= "2020-06-30")
mexico_epoch2 <- filter(mexico_variables, events.confirmed.date >= "2020-07-01" & events.confirmed.date <= "2020-11-30")
mexico_epoch3 <- filter(mexico_variables, events.confirmed.date >= "2020-12-01" & events.confirmed.date <= "2021-03-31")

#number of days in each epoch

mexico_epoch1$length <- 1:124
mexico_epoch2$length <- 124:276
mexico_epoch3$length <- 276:396

#load distribution

onset_to_diagnosis_epoch_1 = fread('fitting_outputs/onset-to-diagnosis_epoch_1_mex-samples-lognormal.csv')
onset_to_diagnosis_epoch_2 = fread('fitting_outputs/onset-to-diagnosis_epoch_2_mex-samples-lognormal.csv')
onset_to_diagnosis_epoch_3 = fread('fitting_outputs/onset-to-diagnosis_epoch_3_mex-samples-lognormal.csv')

#calculate mean delay

onset_to_diagnosis_epoch_1 <- exp((onset_to_diagnosis_epoch_1[,1:32]+0.5)*(onset_to_diagnosis_epoch_1[,33:64]^2)) 
onset_to_diagnosis_epoch_2 <- exp((onset_to_diagnosis_epoch_2[,1:32]+0.5)*(onset_to_diagnosis_epoch_2[,33:64]^2)) 
onset_to_diagnosis_epoch_3 <- exp((onset_to_diagnosis_epoch_3[,1:32]+0.5)*(onset_to_diagnosis_epoch_3[,33:64]^2)) 

#unlist tibble

onset_to_diagnosis_epoch_1 <- as.numeric(unlist(onset_to_diagnosis_epoch_1))
onset_to_diagnosis_epoch_2 <- as.numeric(unlist(onset_to_diagnosis_epoch_2))
onset_to_diagnosis_epoch_3 <- as.numeric(unlist(onset_to_diagnosis_epoch_3))

## calculate deconvolution of observed cases

mex_cases_epoch1 <- get_RL(observed = mexico_epoch1$cases, 
                          times = mexico_epoch1$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_1),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

mex_cases_epoch1 <- head(mex_cases_epoch1,-1)

mex_cases_epoch2 <- get_RL(observed = mexico_epoch2$cases, 
                          times = mexico_epoch2$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_2),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

mex_cases_epoch2 <- head(mex_cases_epoch2,-1)


mex_cases_epoch3 <- get_RL(observed = mexico_epoch3$cases, 
                          times = mexico_epoch3$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_3),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

mex_cases_epoch3 <- head(mex_cases_epoch3,-1)

mex_cases_deconvolved <- na.omit(round(rbind(mex_cases_epoch1,mex_cases_epoch2,mex_cases_epoch3)))

#plot the two incidences against each other 
#to avoid daily fluctuations will use a 5 day average 

mex_cases_deconvolved_average <- as.data.frame(as.numeric(rollapply(mex_cases_deconvolved[,2],7,fill=NA,partial = TRUE,FUN=mean)))

mex_cases_deconvolved_average$date <- seq(as.Date('2020-02-28'), as.Date('2021-03-28'), by = "1 days")
mex_cases_deconvolved_average$group <- 'Deconvolved case counts'
colnames(mex_cases_deconvolved_average) <- c('case','date','group')

mex_cases <- as.data.frame(as.numeric(rollapply(mexico_variables[,2],7,fill=NA,partial = TRUE,FUN=mean)))
mex_cases$date <- seq(as.Date('2020-02-28'), as.Date('2021-03-31'), by = "1 days")
mex_cases$group <- 'Raw case counts'
colnames(mex_cases) <- c('case','date','group')


combined <- rbind(mex_cases_deconvolved_average,mex_cases)

#plot

country <- grobTree(textGrob("Mexico", x=0.05,  y=0.97, hjust=0,
                             gp=gpar(col="Black", fontsize=12, fontface="bold")))


mex_incidence <- ggplot(combined, aes(x = date, y=case, color = factor(group), fill = factor(group))) + 
  geom_line() +
  labs(x= "Time",y= 'Incidence', fill = "Method") +
  theme_bw() +
  scale_x_date (date_breaks = "3 month", date_labels = "%Y %b %d") +
  geom_vline(xintercept = as.Date("2020-06-30"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-11-30"), color = "black", linetype = "dashed") +
  ylim(0,27000) + 
  scale_color_manual(values = c("royalblue", "deeppink"),name = "Methods") + 
  theme(text = element_text(size=12)) +
  annotation_custom(country) 
  

mex_incidence

#Calculate rt

#format cases

#raw case counts

mexico_cases_rt <-  mexico %>% 
  group_by(date =  as.Date(events.confirmed.date)) %>% 
  summarise(number_of_cases = n())

mexico_cases_rt <- pad(mexico_cases_rt)

mexico_cases_rt[is.na(mexico_cases_rt)] <- 0

mexico_cases_rt[, 2] <- cumsum(mexico_cases_rt[, 2])

#Deconvolved case counts 

mex_cases_deconvolved_rt <- mex_cases_deconvolved

mex_cases_deconvolved_rt$date <- seq(as.Date('2020-02-28'), as.Date('2021-03-28'), by = "1 days")

mex_cases_deconvolved_rt <- pad(mex_cases_deconvolved_rt)

mex_cases_deconvolved_rt[is.na(mex_cases_deconvolved_rt)] <- 0

mex_cases_deconvolved_rt[, 2] <- cumsum(mex_cases_deconvolved_rt[, 2])

#Merge together
case_data <- left_join(mexico_cases_rt,mex_cases_deconvolved_rt, by = c("date"))

case_data[is.na(case_data)] <- 0

case_data <- filter(case_data, date >= "2020-02-28" & date <= "2021-03-28")

colnames(case_data) <- c("Date", "raw","Day","deconvoluted")

case_data <- case_data[,c(3,1,2,4)]
case_data$Day <- as.integer(case_data$Day)

##LA countries:
start <- 21 # start = 52 is 29 Feb on both cumulative and daily incidence: this is the data published by WHO on 1 Mar. 
# To compute the incidence on day 52 (29 Feb), the file is read from day 51 (the "start-1" appearing in lines 127 and 151)
# note that start cannot be one as needs to compute the difference in culumative incidence 
Lag <- 0 # Use data till end of dataset

## Open data files
datD <- as.data.frame(case_data)
class(datD)

head(datD)
ncountries = ncol(datD)-2

# Default methods as used in the paper (GCV.Cp with day-of-the-week effect, WD) 
# Further possibilities available for FitMethod, FixedEff
FitMethod <-c('GCV.Cp') # c('GCV.Cp','ML') # Explore fit with different options for the GAM
FixedEff <-c('WD') # c('WE', 'WD', 'None') # Explore alternative fit with weekend effect (WE), day-of-the-week effect (WD), or nothing
colval <-c('black', 'red', 'green')
clist <- c(3,4) # List of the columns for the 5 countries in Figure 1(a)
clistIT <- c(1:3,5) # Relevant columns from Italian data for Figure 1(b)
thin <- 2 # Width of thin line
thick <- 3 # Width of thick line

print('Generating figures for the Electronic Supplementary Material:')

rt <- list()

#change number of interations with double_time function

for(k in 1:length(FitMethod)){
  for(j in 1:length(FixedEff)){
    pdf(file=paste('SuppWHO_',FitMethod[k],'_', FixedEff[j],'_incid.pdf', sep=''))
    for(i in 1:ncountries) {
      a <- DoubleTime(diff(datD[(start-1):(dim(datD)[1]-Lag),2+i]),datD$Day[start:(dim(datD)[1]-Lag)], meth=FitMethod[k], FE=FixedEff[j], plt=T, figtitle=colnames(datD)[2+i])
      print(a)
      rt[[i]] = a 
    }
    dev.off()
  }
}

case_rt <- rt[[1]]
case_rt$dates <- seq(as.Date('2020-03-19'), as.Date('2021-03-28'), by = "1 days")
case_rt$group <- 'Raw case counts'
deconvoluted_rt <- rt[[2]]
deconvoluted_rt$dates <- seq(as.Date('2020-03-19'), as.Date('2021-03-28'), by = "1 days")
deconvoluted_rt$group <- 'Deconvolved case counts'

combined_rt <- rbind(case_rt,deconvoluted_rt)

dev.new()

mex_rt <- ggplot(combined_rt, aes(x = dates, y=sdt, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=sdtlow, ymax=sdtup, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Method") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(r[t]), fill = "Method") +
  geom_vline(xintercept = as.Date("2020-06-30"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-11-30"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "3 month", date_labels = "%Y %b %d") +
  ylim(-0.05,0.2) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  theme(text = element_text(size=12)) +
  annotation_custom(country) 

mex_rt

#Load in Colombia

colombia = as.data.frame(fread("data/colombia_linelist_geocoded_adm1_20220131.csv", drop = c("V1", "_id", 
                                                                                             "caseReference.additionalSources", "caseReference.sourceEntryId", "caseReference.sourceId", "caseReference.sourceUrl","caseReference.uploadIds", "events.firstClinicalConsultation.date", 
                                                                                             "genomeSequences","genomeSequences.repositoryUrl","genomeSequences.sampleCollectionDate","location.place","pathogens",
                                                                                             "preexistingConditions.hasPreexistingConditions","preexistingConditions.values","revisionMetadata.creationMetadata.date","revisionMetadata.creationMetadata.notes","revisionMetadata.editMetadata.date","revisionMetadata.editMetadata.notes","revisionMetadata.revisionNumber","transmission.linkedCaseIds",
                                                                                             "transmission.places", "transmission.routes","travelHistory.travel.dateRange.end","travelHistory.travel.dateRange.start","travelHistory.travel.location.administrativeAreaLevel1","travelHistory.travel.location.administrativeAreaLevel2","travelHistory.travel.location.administrativeAreaLevel3","travelHistory.travel.location.country","travelHistory.travel.location.geometry.coordinates",
                                                                                             "travelHistory.travel.location.geoResolution", "travelHistory.travel.location.name","travelHistory.travel.location.place","travelHistory.travel.methods","travelHistory.travel.purpose","travelHistory.traveledPrior30Days","vaccines.0.batch","vaccines.0.date","vaccines.0.name",
                                                                                             "vaccines.0.sideEffects", "vaccines.1.batch","vaccines.1.date","vaccines.1.name","vaccines.1.sideEffects","vaccines.2.batch","vaccines.2.date","vaccines.2.name","vaccines.2.sideEffects",
                                                                                             "vaccines.3.batch", "vaccines.3.date", "vaccines.3.name","vaccines.3.sideEffects",
                                                                                             "caseReference.verificationStatus", "demographics.ageRange.end","demographics.ageRange.start","demographics.ethnicity","demographics.gender","demographics.nationalities","demographics.occupation","genomeSequences.sequenceId","genomeSequences.sequenceLength",
                                                                                             "genomeSequences.sequenceName")))

#Calculate number of cases per day

colombia_variables <-  colombia %>% 
  group_by(events.confirmed.date) %>%
  summarise(cases = n())

#fill in missing gaps in dates

colombia_variables <- as.data.frame(colombia_variables) %>%
  mutate(events.confirmed.date = as.Date(events.confirmed.date)) %>%
  complete(events.confirmed.date = seq.Date(min(events.confirmed.date), max(events.confirmed.date), by="day"))

colombia_variables[is.na(colombia_variables)] <- 0

colombia_variables <- filter(colombia_variables, events.confirmed.date >= "2020-03-06" & events.confirmed.date <= "2021-03-31")

#filer for epoch

colombia_epoch1 <- filter (colombia_variables,events.confirmed.date >= "2020-03-06" & events.confirmed.date <= "2020-06-30")
colombia_epoch2 <- filter(colombia_variables, events.confirmed.date >= "2020-07-01" & events.confirmed.date <= "2020-11-30")
colombia_epoch3 <- filter(colombia_variables, events.confirmed.date >= "2020-12-01" & events.confirmed.date <= "2021-03-31")

#number of days in each epoch

colombia_epoch1$length <- 1:117
colombia_epoch2$length <- 117:269
colombia_epoch3$length <- 269:389

#load distribution

onset_to_diagnosis_epoch_1 = fread('fitting_outputs/onset-to-diagnosis_epoch_1_co-samples-lognormal.csv')
onset_to_diagnosis_epoch_2 = fread('fitting_outputs/onset-to-diagnosis_epoch_2_co-samples-lognormal.csv')
onset_to_diagnosis_epoch_3 = fread('fitting_outputs/onset-to-diagnosis_epoch_3_co-samples-lognormal.csv')

#calculate mean delay

onset_to_diagnosis_epoch_1 <- exp((onset_to_diagnosis_epoch_1[,1:30]+0.5)*(onset_to_diagnosis_epoch_1[,31:60]^2)) 
onset_to_diagnosis_epoch_2 <- exp((onset_to_diagnosis_epoch_2[,1:30]+0.5)*(onset_to_diagnosis_epoch_2[,31:60]^2)) 
onset_to_diagnosis_epoch_3 <- exp((onset_to_diagnosis_epoch_3[,1:30]+0.5)*(onset_to_diagnosis_epoch_3[,31:60]^2)) 

#unlist tibble

onset_to_diagnosis_epoch_1 <- as.numeric(unlist(onset_to_diagnosis_epoch_1))
onset_to_diagnosis_epoch_2 <- as.numeric(unlist(onset_to_diagnosis_epoch_2))
onset_to_diagnosis_epoch_3 <- as.numeric(unlist(onset_to_diagnosis_epoch_3))

## calculate deconvolution of observed cases

co_cases_epoch1 <- get_RL(observed = colombia_epoch1$cases, 
                          times = colombia_epoch1$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_1),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

co_cases_epoch1 <- tail(co_cases_epoch1, -8)

co_cases_epoch1 <- head(co_cases_epoch1,-1)

co_cases_epoch2 <- get_RL(observed = colombia_epoch2$cases, 
                          times = colombia_epoch2$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_2),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

co_cases_epoch2 <- head(co_cases_epoch2,-1)

onset_to_diagnosis_epoch_3 <- as.data.frame(onset_to_diagnosis_epoch_3)

onset_to_diagnosis_epoch_3 <- filter(onset_to_diagnosis_epoch_3, onset_to_diagnosis_epoch_3 <= 30)

onset_to_diagnosis_epoch_3 <- as.numeric(onset_to_diagnosis_epoch_3$onset_to_diagnosis_epoch_3)

co_cases_epoch3 <- get_RL(observed = colombia_epoch3$cases, 
                          times = colombia_epoch3$length, 
                          p_delay = get_delay_multinom(onset_to_diagnosis_epoch_3),
                          out_col_name = 'case_deconvolved',
                          right_censor = FALSE)

co_cases_epoch3 <- head(co_cases_epoch3,-1)

co_cases_deconvolved <- na.omit(round(rbind(co_cases_epoch1,co_cases_epoch2,co_cases_epoch3)))

#plot the two incidences against each other 
#to avoid daily fluctuations will use a 7 day average 


co_cases_deconvolved_average <- as.data.frame(as.numeric(rollapply(co_cases_deconvolved[,2],7,fill=NA,partial = TRUE,FUN=mean)))

co_cases_deconvolved_average$date <- seq(as.Date('2020-03-05'), as.Date('2021-03-29'), by = "1 days")
co_cases_deconvolved_average$group <- 'Deconvolved case counts'
colnames(co_cases_deconvolved_average) <- c('case','date','group')

co_cases <- as.data.frame(as.numeric(rollapply(colombia_variables[,2],7,fill=NA,partial = TRUE,FUN=mean)))
co_cases$date <- seq(as.Date('2020-03-06'), as.Date('2021-03-31'), by = "1 days")
co_cases$group <- 'Raw case counts'
colnames(co_cases) <- c('case','date','group')


combined <- rbind(co_cases_deconvolved_average,co_cases)

#plot

country <- grobTree(textGrob("Colombia", x=0.05,  y=0.97, hjust=0,
                             gp=gpar(col="Black", fontsize=12, fontface="bold")))

co_incidence <- ggplot(combined, aes(x = date, y=case, color = factor(group), fill = factor(group))) + 
  geom_line() +
  labs(x= "Time",y= 'Incidence', fill = "Method") +
  theme_bw() +
  scale_x_date (date_breaks = "3 month", date_labels = "%Y %b %d") +
  geom_vline(xintercept = as.Date("2020-06-30"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-11-30"), color = "black", linetype = "dashed") +
  ylim(0,27000) + 
  scale_color_manual(values = c("royalblue", "deeppink"),name = "Methods") + 
  theme(text = element_text(size=12)) +
  annotation_custom(country) 

co_incidence

#Calculate rt

#format cases

#raw case counts

colombia_cases_rt <-  colombia %>% 
  group_by(date =  as.Date(events.confirmed.date)) %>% 
  summarise(number_of_cases = n())

colombia_cases_rt <- pad(colombia_cases_rt)

colombia_cases_rt[is.na(colombia_cases_rt)] <- 0

colombia_cases_rt[, 2] <- cumsum(colombia_cases_rt[, 2])

#Deconvolved case counts 

co_cases_deconvolved_rt <- co_cases_deconvolved

co_cases_deconvolved_rt$date <- seq(as.Date('2020-03-05'), as.Date('2021-03-29'), by = "1 days")

co_cases_deconvolved_rt <- pad(co_cases_deconvolved_rt)

co_cases_deconvolved_rt[is.na(co_cases_deconvolved_rt)] <- 0

co_cases_deconvolved_rt[, 2] <- cumsum(co_cases_deconvolved_rt[, 2])

#Merge together
case_data <- left_join(colombia_cases_rt,co_cases_deconvolved_rt, by = c("date"))

case_data[is.na(case_data)] <- 0

case_data <- filter(case_data, date >= "2020-03-05" & date <= "2021-03-29")

colnames(case_data) <- c("Date", "raw","Day","deconvoluted")

case_data <- case_data[,c(3,1,2,4)]
case_data$Day <- as.integer(case_data$Day)

summary(ar_cases_deconvolved_rt$case_deconvolved)

##LA countries:
start <- 13 # start = 52 is 29 Feb on both cumulative and daily incidence: this is the data published by WHO on 1 Mar. 
# To compute the incidence on day 52 (29 Feb), the file is read from day 51 (the "start-1" appearing in lines 127 and 151)
# note that start cannot be one as needs to compute the difference in culumative incidence 
Lag <- 0 # Use data till end of dataset

## Open data files
datD <- as.data.frame(case_data)
class(datD)

head(datD)
ncountries = ncol(datD)-2

# Default methods as used in the paper (GCV.Cp with day-of-the-week effect, WD) 
# Further possibilities available for FitMethod, FixedEff
FitMethod <-c('GCV.Cp') # c('GCV.Cp','ML') # Explore fit with different options for the GAM
FixedEff <-c('WD') # c('WE', 'WD', 'None') # Explore alternative fit with weekend effect (WE), day-of-the-week effect (WD), or nothing
colval <-c('black', 'red', 'green')
clist <- c(3,4) # List of the columns for the 5 countries in Figure 1(a)
clistIT <- c(1:3,5) # Relevant columns from Italian data for Figure 1(b)
thin <- 2 # Width of thin line
thick <- 3 # Width of thick line

print('Generating figures for the Electronic Supplementary Material:')

rt <- list()

#change number of interations with double_time function - to 390

for(k in 1:length(FitMethod)){
  for(j in 1:length(FixedEff)){
    pdf(file=paste('SuppWHO_',FitMethod[k],'_', FixedEff[j],'_incid.pdf', sep=''))
    for(i in 1:ncountries) {
      a <- DoubleTime(diff(datD[(start-1):(dim(datD)[1]-Lag),2+i]),datD$Day[start:(dim(datD)[1]-Lag)], meth=FitMethod[k], FE=FixedEff[j], plt=T, figtitle=colnames(datD)[2+i])
      print(a)
      rt[[i]] = a 
    }
    dev.off()
  }
}

case_rt <- rt[[1]]
case_rt$dates <- seq(as.Date('2020-03-19'), as.Date('2021-03-28'), by = "1 days")
case_rt$group <- 'Raw case counts'
deconvoluted_rt <- rt[[2]]
deconvoluted_rt$dates <- seq(as.Date('2020-03-19'), as.Date('2021-03-28'), by = "1 days")
deconvoluted_rt$group <- 'Deconvolved case counts'

combined_rt <- rbind(case_rt,deconvoluted_rt)

dev.new()

co_rt <- ggplot(combined_rt, aes(x = dates, y=sdt, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=sdtlow, ymax=sdtup, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Method") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(r[t]), fill = "Method") +
  theme_bw() +
  scale_x_date (date_breaks = "3 month", date_labels = "%Y %b %d") +
  geom_vline(xintercept = as.Date("2020-06-30"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-11-30"), color = "black", linetype = "dashed") +
  ylim(-0.05,0.2) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  theme(text = element_text(size=12)) +
  annotation_custom(country) 

co_rt

#merge all plots together 

ggarrange(
  ar_rt, br_rt, co_rt, mex_rt, labels = c("A", "B", 'C','D'),
  common.legend = TRUE
)

ggarrange(
  ar_incidence, br_incidence, co_incidence, mex_incidence, labels = c("A", "B", 'C','D'),
  common.legend = TRUE
)

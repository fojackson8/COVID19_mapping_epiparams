#Process Gh data

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

#set WD
setwd("~/rhys/pystan")

# Main functions to run code

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
require(ggdistribute)

###Brazil###

#Load data

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

###Onset to Diagnosis###

#need to get State, date, Onset to diagnosis, Onset to Hospitilsation and Onset to Death

brazil_variables <-  brazil %>% 
  summarise('onset-to-diagnosis' = delta_symptom_onset_to_confirmed,
            'onset-to-death' = delta_symptom_onset_to_death,
            'onset-to-hospitalisation' = delta_symptom_onset_to_hospitalisation,
            'hospitalisation_to_death' = delta_hospitalisation_to_death,
            State = ADM1_PT,
            Date =  events.confirmed.date )

#Split into three epochs start - 30th June 2020, 1st July - 30th November 2020, 1 December - 31st March 2021

brazil_variables$Date <- as.Date (brazil_variables$Date)

brazil_epoch1 <- filter (brazil_variables, Date <= "2020-06-30")

brazil_epoch2 <- filter(brazil_variables, Date >= "2020-07-01" & Date <= "2020-11-30")

brazil_epoch3 <- filter(brazil_variables, Date >= "2020-12-01" & Date <= "2021-03-31")

#split into each variable and select only complete cases and remove outliars 

#epoch1

brazil_epoch1_onset_to_diagnosis <- brazil_epoch1 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_1_br' = `onset-to-diagnosis`,)

brazil_epoch1_onset_to_diagnosis <- data_process(brazil_epoch1_onset_to_diagnosis)

brazil_epoch1_onset_to_death <- brazil_epoch1 %>% 
  summarise(State = State,
            'onset-to-death_epoch_1_br' = `onset-to-death`)

brazil_epoch1_onset_to_death <- data_process_death(brazil_epoch1_onset_to_death)

brazil_epoch1_onset_to_hospitalisation <- brazil_epoch1 %>% 
  summarise(State = State,
            'onset-to-hospitalisation_epoch_1_br' = `onset-to-hospitalisation`)

brazil_epoch1_onset_to_hospitalisation <- data_process(brazil_epoch1_onset_to_hospitalisation)

brazil_epoch1_hospitalisation_to_death <- brazil_epoch1 %>% 
  summarise(State = State,
            'hospitalisation_to_death_epoch_1_br' = `hospitalisation_to_death`)

brazil_epoch1_hospitalisation_to_death <- data_process_death(brazil_epoch1_hospitalisation_to_death)

#epoch2

brazil_epoch2_onset_to_diagnosis <- brazil_epoch2 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_2_br' = `onset-to-diagnosis`,)

brazil_epoch2_onset_to_diagnosis <- data_process(brazil_epoch2_onset_to_diagnosis)

brazil_epoch2_onset_to_death <- brazil_epoch2 %>% 
  summarise(State = State,
    'onset-to-death_epoch_2_br' = `onset-to-death`)

brazil_epoch2_onset_to_death <- data_process(brazil_epoch2_onset_to_death)

brazil_epoch2_onset_to_hospitalisation <- brazil_epoch2 %>% 
  summarise(State = State,
    'onset-to-hospitalisation_epoch_2_br' = `onset-to-hospitalisation`)

brazil_epoch2_onset_to_hospitalisation <- data_process(brazil_epoch2_onset_to_hospitalisation)

brazil_epoch2_hospitalisation_to_death <- brazil_epoch2 %>% 
  summarise(State = State,
            'hospitalisation_to_death_epoch_2_br' = `hospitalisation_to_death`)

brazil_epoch2_hospitalisation_to_death <- data_process_death(brazil_epoch2_hospitalisation_to_death)

#epoch3

brazil_epoch3_onset_to_diagnosis <- brazil_epoch3 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_3_br' = `onset-to-diagnosis`)

brazil_epoch3_onset_to_diagnosis <- data_process(brazil_epoch3_onset_to_diagnosis)

brazil_epoch3_onset_to_death <- brazil_epoch3 %>% 
  summarise(State = State,
    'onset-to-death_epoch_3_br' = `onset-to-death`)

brazil_epoch3_onset_to_death <- data_process(brazil_epoch3_onset_to_death)

brazil_epoch3_onset_to_hospitalisation <- brazil_epoch3 %>% 
  summarise(State = State,
    'onset-to-hospitalisation_epoch_3_br' = `onset-to-hospitalisation`)

brazil_epoch3_onset_to_hospitalisation <- data_process(brazil_epoch3_onset_to_hospitalisation)

brazil_epoch3_hospitalisation_to_death <- brazil_epoch3 %>% 
  summarise(State = State,
            'hospitalisation_to_death_epoch_3_br' = `hospitalisation_to_death`)

brazil_epoch3_hospitalisation_to_death <- data_process_death(brazil_epoch3_hospitalisation_to_death)

#save data

write.csv(brazil_epoch1_onset_to_diagnosis, file = "data_for_pystan/brazil_epoch1_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch2_onset_to_diagnosis, file = "data_for_pystan/brazil_epoch2_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch3_onset_to_diagnosis, file = "data_for_pystan/brazil_epoch3_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch1_onset_to_death, file = "data_for_pystan/brazil_epoch1_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch2_onset_to_death, file = "data_for_pystan/brazil_epoch2_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch3_onset_to_death, file = "data_for_pystan/brazil_epoch3_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch1_onset_to_hospitalisation, file = "data_for_pystan/brazil_epoch1_onset_to_hospitalisation.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch2_onset_to_hospitalisation, file = "data_for_pystan/brazil_epoch2_onset_to_hospitalisation.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch3_onset_to_hospitalisation, file = "data_for_pystan/brazil_epoch3_onset_to_hospitalisation.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch1_hospitalisation_to_death, file = "data_for_pystan/brazil_epoch1_hospitalisation_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch2_hospitalisation_to_death, file = "data_for_pystan/brazil_epoch2_hospitalisation_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(brazil_epoch3_hospitalisation_to_death, file = "data_for_pystan/brazil_epoch3_hospitalisation_to_death.csv",row.names=FALSE, quote = TRUE)

#Colombia 

#Load data

columbia = as.data.frame(fread("data/colombia_linelist_geocoded_adm1_20220131.csv", drop = c("V1", "_id", 
                                                                      "caseReference.additionalSources", "caseReference.sourceEntryId", "caseReference.sourceId", "caseReference.sourceUrl","caseReference.uploadIds", "events.firstClinicalConsultation.date", 
                                                                      "genomeSequences","genomeSequences.repositoryUrl","genomeSequences.sampleCollectionDate","location.place","pathogens",
                                                                      "preexistingConditions.hasPreexistingConditions","preexistingConditions.values","revisionMetadata.creationMetadata.date","revisionMetadata.creationMetadata.notes","revisionMetadata.editMetadata.date","revisionMetadata.editMetadata.notes","revisionMetadata.revisionNumber","transmission.linkedCaseIds",
                                                                      "transmission.places", "transmission.routes","travelHistory.travel.dateRange.end","travelHistory.travel.dateRange.start","travelHistory.travel.location.administrativeAreaLevel1","travelHistory.travel.location.administrativeAreaLevel2","travelHistory.travel.location.administrativeAreaLevel3","travelHistory.travel.location.country","travelHistory.travel.location.geometry.coordinates",
                                                                      "travelHistory.travel.location.geoResolution", "travelHistory.travel.location.name","travelHistory.travel.location.place","travelHistory.travel.methods","travelHistory.travel.purpose","travelHistory.traveledPrior30Days","vaccines.0.batch","vaccines.0.date","vaccines.0.name",
                                                                      "vaccines.0.sideEffects", "vaccines.1.batch","vaccines.1.date","vaccines.1.name","vaccines.1.sideEffects","vaccines.2.batch","vaccines.2.date","vaccines.2.name","vaccines.2.sideEffects",
                                                                      "vaccines.3.batch", "vaccines.3.date", "vaccines.3.name","vaccines.3.sideEffects",
                                                                      "caseReference.verificationStatus", "demographics.ageRange.end","demographics.ageRange.start","demographics.ethnicity","demographics.gender","demographics.nationalities","demographics.occupation","genomeSequences.sequenceId","genomeSequences.sequenceLength",
                                                                      "genomeSequences.sequenceName")))

###Onset to Diagnosis###

#need to get State, date, Onset to diagnosis and Onset to Death
#Note no onset_to_hospitilsation

columbia_variables <-  columbia %>% 
  summarise('onset-to-diagnosis' = delta_symptom_onset_to_confirmed,
            'onset-to-death' = delta_symptom_onset_to_death,
            State = ADM1_ES,
            Date =  events.confirmed.date)


#Split into three epochs start - 30th June 2020, 1st July - 30th November 2020, 1 December - 31st March 2021

columbia_variables$Date <- as.Date (columbia_variables$Date)

#advalability of data for diagnosis and death

col_diagnosis <- dplyr::select(columbia_variables,c(1,3:4))

col_diagnosis <- na.omit(col_diagnosis)

col_date_diagnosis <- col_diagnosis %>%
  group_by(State) %>%
  arrange(Date) %>%
  slice(1L)

col_death <- dplyr::select(columbia_variables,c(2,3:4))

col_death <- na.omit(col_death)

col_date_death <- col_death %>%
  group_by(State) %>%
  arrange(Date) %>%
  slice(1L)

columbia_dates <- left_join(col_date_diagnosis,col_date_death, by =  c("State" = "State"))

columbia_dates <- dplyr::select(columbia_dates,c(2,3,5))

colnames(columbia_dates) <- c("State","Date of first Diagnosis","date of first Death")

#Choose which states can be included

#remove Guaviare from death 

columbia_epoch1 <- filter (columbia_variables, Date <= "2020-06-30")

columbia_epoch2 <- filter(columbia_variables, Date >= "2020-07-01" & Date <= "2020-11-30")

columbia_epoch3 <- filter(columbia_variables, Date >= "2020-12-01" & Date <= "2021-03-31")

#split into each variable and select only complete cases and remove outliars 

#epoch1

columbia_epoch1_onset_to_diagnosis <- columbia_epoch1 %>% 
  dplyr :: summarise(State = State,
            'onset-to-diagnosis_epoch_1_co' = `onset-to-diagnosis`,)

columbia_epoch1_onset_to_diagnosis <- data_process(columbia_epoch1_onset_to_diagnosis)

columbia_epoch1_onset_to_death <- columbia_epoch1 %>% 
  summarise(State = State,
            'onset-to-death_epoch_1_co' = `onset-to-death`)

columbia_epoch1_onset_to_death <- data_process_death(columbia_epoch1_onset_to_death)

columbia_epoch1_onset_to_death <- filter(columbia_epoch1_onset_to_death, State != 'Guaviare')

#epoch2

columbia_epoch2_onset_to_diagnosis <- columbia_epoch2 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_2_co' = `onset-to-diagnosis`,)

columbia_epoch2_onset_to_diagnosis <- data_process(columbia_epoch2_onset_to_diagnosis)

columbia_epoch2_onset_to_diagnosis <- filter(columbia_epoch2_onset_to_diagnosis, State %in% columbia_epoch1_onset_to_diagnosis$State)

columbia_epoch2_onset_to_death <- columbia_epoch2 %>% 
  summarise(State = State,
            'onset-to-death_epoch_2_co' = `onset-to-death`)

columbia_epoch2_onset_to_death <- data_process_death(columbia_epoch2_onset_to_death)

columbia_epoch2_onset_to_death <- filter(columbia_epoch2_onset_to_death, State %in% columbia_epoch1_onset_to_death$State)


#epoch3

columbia_epoch3_onset_to_diagnosis <- columbia_epoch3 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_3_co' = `onset-to-diagnosis`)

columbia_epoch3_onset_to_diagnosis <- data_process(columbia_epoch3_onset_to_diagnosis)

columbia_epoch3_onset_to_diagnosis <- filter(columbia_epoch3_onset_to_diagnosis, State %in% columbia_epoch1_onset_to_diagnosis$State)

columbia_epoch3_onset_to_death <- columbia_epoch3 %>% 
  summarise(State = State,
            'onset-to-death_epoch_3_co' = `onset-to-death`)

columbia_epoch3_onset_to_death <- data_process_death(columbia_epoch3_onset_to_death)

columbia_epoch3_onset_to_death <- filter(columbia_epoch3_onset_to_death, State %in% columbia_epoch1_onset_to_death$State)

#save data

write.csv(columbia_epoch1_onset_to_diagnosis, file = "data_for_pystan/colombia_epoch1_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(columbia_epoch2_onset_to_diagnosis, file = "data_for_pystan/colombia_epoch2_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(columbia_epoch3_onset_to_diagnosis, file = "data_for_pystan/colombia_epoch3_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(columbia_epoch1_onset_to_death, file = "data_for_pystan/colombia_epoch1_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(columbia_epoch2_onset_to_death, file = "data_for_pystan/colombia_epoch2_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(columbia_epoch3_onset_to_death, file = "data_for_pystan/colombia_epoch3_onset_to_death.csv",row.names=FALSE, quote = TRUE)



#Argentina

#Load data

argentina = as.data.frame(read_csv("data/argentina_linelist_geocoded_adm1_20220131.csv"))
                          
                         
###Onset to Diagnosis###

#need to get State, date, Onset to diagnosis, Onset to Hospitilsation and Onset to Death

argentina_variables <-  argentina %>% 
  summarise('onset-to-diagnosis' = delta_symptom_onset_to_confirmed,
            'onset-to-death' = delta_symptom_onset_to_death,
            'onset-to-hospitalisation' = delta_symptom_onset_to_hospitalisation,
            State = ADM1_ES,
            Date =  events.confirmed.date )

#Split into three epochs start - 30th June 2020, 1st July - 30th November 2020, 1 December - 31st March 2021

argentina_variables$Date <- as.Date (argentina_variables$Date)

##determine which states to include in anaylsis

ar_diagnosis <- dplyr::select(argentina_variables,c(1,4:5))

ar_diagnosis <- na.omit(ar_diagnosis)

ar_date_diagnosis <- ar_diagnosis %>%
  group_by(State) %>%
  arrange(Date) %>%
  slice(1L)

ar_hosp <- dplyr::select(argentina_variables,c(3,4:5))

ar_hosp <- na.omit(ar_hosp)

ar_date_hosp <- ar_hosp %>%
  group_by(State) %>%
  arrange(Date) %>%
  slice(1L)

ar_death <- dplyr::select(argentina_variables,c(2,4:5))

ar_death <- na.omit(ar_death)

ar_date_death <- ar_death %>%
  group_by(State) %>%
  arrange(Date) %>%
  slice(1L)

argentina_dates <- left_join(ar_date_diagnosis,ar_date_death, by =  c("State" = "State"))

argentina_dates <- left_join(argentina_dates,ar_date_hosp, by =  c("State" = "State"))

argentina_dates <- dplyr::select(argentina_dates,c(2,3,5,7))

colnames(argentina_dates) <- c("State","Date of first Diagnosis","Date of first Death","Date of first hospitalisation")


##diagnosis include all, deaths exclude - La Pampa, Salta, Santa Cruz, Formosa, Tucumán, San Luis, 
##hospitalisations exclude - Formosa, La Pampa, Santa Cruz, Misiones, San Luis, Salta, Chubut, Tucumán, 

argentina_epoch1 <- filter (argentina_variables, Date <= "2020-06-30")

argentina_epoch2 <- filter(argentina_variables, Date >= "2020-07-01" & Date <= "2020-11-30")

argentina_epoch3 <- filter(argentina_variables, Date >= "2020-12-01" & Date <= "2021-03-31")

#split into each variable and select only complete cases and remove outliars 

#epoch1

argentina_epoch1_onset_to_diagnosis <- argentina_epoch1 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_1_ar' = `onset-to-diagnosis`,)

argentina_epoch1_onset_to_diagnosis <- data_process(argentina_epoch1_onset_to_diagnosis)

argentina_epoch1_onset_to_death <- argentina_epoch1 %>% 
  summarise(State = State,
            'onset-to-death_epoch_1_ar' = `onset-to-death`)

argentina_epoch1_onset_to_death <- data_process_death(argentina_epoch1_onset_to_death)

death_state_exclude <- c('La Pampa', 'Salta', 'Santa Cruz', 'Formosa', 'Tucumán', 'San Luis')

argentina_epoch1_onset_to_death <- filter(argentina_epoch1_onset_to_death, !State %in% death_state_exclude)

argentina_epoch1_onset_to_hospitalisation <- argentina_epoch1 %>% 
  summarise(State = State,
            'onset-to-hospitalisation_epoch_1_ar' = `onset-to-hospitalisation`)

argentina_epoch1_onset_to_hospitalisation <- data_process(argentina_epoch1_onset_to_hospitalisation)

hosp_state_exclude <- c('La Pampa', 'Salta', 'Santa Cruz', 'Formosa', 'Tucumán', 'San Luis','Chubut','Misiones','Chaco')

argentina_epoch1_onset_to_hospitalisation <- filter(argentina_epoch1_onset_to_hospitalisation, !State %in% hosp_state_exclude)

#epoch2

argentina_epoch2_onset_to_diagnosis <- argentina_epoch2 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_2_ar' = `onset-to-diagnosis`,)

argentina_epoch2_onset_to_diagnosis <- data_process(argentina_epoch2_onset_to_diagnosis)

argentina_epoch2_onset_to_death <- argentina_epoch2 %>% 
  summarise(State = State,
            'onset-to-death_epoch_2_ar' = `onset-to-death`)

argentina_epoch2_onset_to_death <- data_process_death(argentina_epoch2_onset_to_death)

argentina_epoch2_onset_to_death <- filter(argentina_epoch2_onset_to_death, !State %in% death_state_exclude)

argentina_epoch2_onset_to_hospitalisation <- argentina_epoch2 %>% 
  summarise(State = State,
            'onset-to-hospitalisation_epoch_2_ar' = `onset-to-hospitalisation`)

argentina_epoch2_onset_to_hospitalisation <- data_process(argentina_epoch2_onset_to_hospitalisation)

argentina_epoch2_onset_to_hospitalisation <- filter(argentina_epoch2_onset_to_hospitalisation, !State %in% hosp_state_exclude)

#epoch3

argentina_epoch3_onset_to_diagnosis <- argentina_epoch3 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_3_ar' = `onset-to-diagnosis`)

argentina_epoch3_onset_to_diagnosis <- data_process(argentina_epoch3_onset_to_diagnosis)

argentina_epoch3_onset_to_death <- argentina_epoch3 %>% 
  summarise(State = State,
            'onset-to-death_epoch_3_ar' = `onset-to-death`)

argentina_epoch3_onset_to_death <- data_process_death(argentina_epoch3_onset_to_death)

argentina_epoch3_onset_to_death <- filter(argentina_epoch3_onset_to_death, !State %in% death_state_exclude)

argentina_epoch3_onset_to_hospitalisation <- argentina_epoch3 %>% 
  summarise(State = State,
            'onset-to-hospitalisation_epoch_3_ar' = `onset-to-hospitalisation`)

argentina_epoch3_onset_to_hospitalisation <- data_process(argentina_epoch3_onset_to_hospitalisation)

argentina_epoch3_onset_to_hospitalisation <- filter(argentina_epoch3_onset_to_hospitalisation, !State %in% hosp_state_exclude)

#save data

write.csv(argentina_epoch1_onset_to_diagnosis, file = "data_for_pystan/argentina_epoch1_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(argentina_epoch2_onset_to_diagnosis, file = "data_for_pystan/argentina_epoch2_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(argentina_epoch3_onset_to_diagnosis, file = "data_for_pystan/argentina_epoch3_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(argentina_epoch1_onset_to_hospitalisation, file = "data_for_pystan/argentina_epoch1_onset_to_hospitalisation.csv",row.names=FALSE, quote = TRUE)
write.csv(argentina_epoch2_onset_to_hospitalisation, file = "data_for_pystan/argentina_epoch2_onset_to_hospitalisation.csv",row.names=FALSE, quote = TRUE)
write.csv(argentina_epoch3_onset_to_hospitalisation, file = "data_for_pystan/argentina_epoch3_onset_to_hospitalisation.csv",row.names=FALSE, quote = TRUE)
write.csv(argentina_epoch1_onset_to_death, file = "data_for_pystan/argentina_epoch1_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(argentina_epoch2_onset_to_death, file = "data_for_pystan/argentina_epoch2_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(argentina_epoch3_onset_to_death, file = "data_for_pystan/argentina_epoch3_onset_to_death.csv",row.names=FALSE, quote = TRUE)

#Mexico

#Load data

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

###Onset to Diagnosis###

#need to get State, date, Onset to diagnosis, Onset to Hospitilsation and Onset to Death

mexico_variables <-  mexico %>% 
  summarise('onset-to-diagnosis' = delta_symptom_onset_to_confirmed,
            'onset-to-death' = delta_symptom_onset_to_death,
            'onset-to-hospitalisation' = delta_symptom_onset_to_hospitalisation,
            'hospitalisation_to_death' = delta_hospitalisation_to_death,
            State = ADM1_ES,
            Date =  events.confirmed.date )

#Split into three epochs start - 30th June 2020, 1st July - 30th November 2020, 1 December - 31st March 2021

mexico_variables$Date <- as.Date (mexico_variables$Date)

mexico_epoch1 <- filter (mexico_variables, Date <= "2020-06-30")

mexico_epoch2 <- filter(mexico_variables, Date >= "2020-07-01" & Date <= "2020-11-30")

mexico_epoch3 <- filter(mexico_variables, Date >= "2020-12-01" & Date <= "2021-03-31")

#split into each variable and select only complete cases and remove outliars 

#epoch1

mexico_epoch1_onset_to_diagnosis <- mexico_epoch1 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_1_mex' = `onset-to-diagnosis`,)

mexico_epoch1_onset_to_diagnosis <- data_process(mexico_epoch1_onset_to_diagnosis)

mexico_epoch1_onset_to_death <- mexico_epoch1 %>% 
  summarise(State = State,
            'onset-to-death_epoch_1_mex' = `onset-to-death`)

mexico_epoch1_onset_to_death <- data_process_death(mexico_epoch1_onset_to_death)

mexico_epoch1_onset_to_hospitalisation <- mexico_epoch1 %>% 
  summarise(State = State,
            'onset-to-hospitalisation_epoch_1_mex' = `onset-to-hospitalisation`)

mexico_epoch1_onset_to_hospitalisation <- data_process(mexico_epoch1_onset_to_hospitalisation)

mexico_epoch1_hospitalisation_to_death <- mexico_epoch1 %>% 
  summarise(State = State,
            'hospitalisation_to_death_epoch_1_mex' = `hospitalisation_to_death`)

mexico_epoch1_hospitalisation_to_death <- data_process_death(mexico_epoch1_hospitalisation_to_death)

#epoch2

mexico_epoch2_onset_to_diagnosis <- mexico_epoch2 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_2_mex' = `onset-to-diagnosis`,)

mexico_epoch2_onset_to_diagnosis <- data_process(mexico_epoch2_onset_to_diagnosis)

mexico_epoch2_onset_to_death <- mexico_epoch2 %>% 
  summarise(State = State,
            'onset-to-death_epoch_2_mex' = `onset-to-death`)

mexico_epoch2_onset_to_death <- data_process(mexico_epoch2_onset_to_death)

mexico_epoch2_onset_to_hospitalisation <- mexico_epoch2 %>% 
  summarise(State = State,
            'onset-to-hospitalisation_epoch_2_mex' = `onset-to-hospitalisation`)

mexico_epoch2_onset_to_hospitalisation <- data_process(mexico_epoch2_onset_to_hospitalisation)

mexico_epoch2_hospitalisation_to_death <- mexico_epoch2 %>% 
  summarise(State = State,
            'hospitalisation_to_death_epoch_2_mex' = `hospitalisation_to_death`)

mexico_epoch2_hospitalisation_to_death <- data_process_death(mexico_epoch2_hospitalisation_to_death)

#epoch3

mexico_epoch3_onset_to_diagnosis <- mexico_epoch3 %>% 
  summarise(State = State,
            'onset-to-diagnosis_epoch_3_mex' = `onset-to-diagnosis`)

mexico_epoch3_onset_to_diagnosis <- data_process(mexico_epoch3_onset_to_diagnosis)

mexico_epoch3_onset_to_death <- mexico_epoch3 %>% 
  summarise(State = State,
            'onset-to-death_epoch_3_mex' = `onset-to-death`)

mexico_epoch3_onset_to_death <- data_process(mexico_epoch3_onset_to_death)

mexico_epoch3_onset_to_hospitalisation <- mexico_epoch3 %>% 
  summarise(State = State,
            'onset-to-hospitalisation_epoch_3_mex' = `onset-to-hospitalisation`)

mexico_epoch3_onset_to_hospitalisation <- data_process(mexico_epoch3_onset_to_hospitalisation)

mexico_epoch3_hospitalisation_to_death <- mexico_epoch3 %>% 
  summarise(State = State,
            'hospitalisation_to_death_epoch_3_mex' = `hospitalisation_to_death`)

mexico_epoch3_hospitalisation_to_death <- data_process_death(mexico_epoch3_hospitalisation_to_death)

#save data

write.csv(mexico_epoch1_onset_to_diagnosis, file = "data_for_pystan/mexico_epoch1_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch2_onset_to_diagnosis, file = "data_for_pystan/mexico_epoch2_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch3_onset_to_diagnosis, file = "data_for_pystan/mexico_epoch3_onset_to_diagnosis.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch1_onset_to_death, file = "data_for_pystan/mexico_epoch1_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch2_onset_to_death, file = "data_for_pystan/mexico_epoch2_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch3_onset_to_death, file = "data_for_pystan/mexico_epoch3_onset_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch1_onset_to_hospitalisation, file = "data_for_pystan/mexico_epoch1_onset_to_hospitalisation.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch2_onset_to_hospitalisation, file = "data_for_pystan/mexico_epoch2_onset_to_hospitalisation.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch3_onset_to_hospitalisation, file = "data_for_pystan/mexico_epoch3_onset_to_hospitalisation.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch1_hospitalisation_to_death, file = "data_for_pystan/mexico_epoch1_hospitalisation_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch2_hospitalisation_to_death, file = "data_for_pystan/mexico_epoch2_hospitalisation_to_death.csv",row.names=FALSE, quote = TRUE)
write.csv(mexico_epoch3_hospitalisation_to_death, file = "data_for_pystan/mexico_epoch3_hospitalisation_to_death.csv",row.names=FALSE, quote = TRUE)


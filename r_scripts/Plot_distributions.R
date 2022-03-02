#Process distributions 

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

#set WD
setwd("~/rhys/pystan")

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
library(plyr)

#Brazil

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

###Onset to death###


brazil_variables <-  brazil %>% 
  summarise('onset-to-confirmed' = delta_symptom_onset_to_confirmed,
            'onset-to-death' = delta_symptom_onset_to_death,
            'onset-to-hospitalisation' = delta_symptom_onset_to_hospitalisation,
            State = ADM1_PT,
            Date =  events.confirmed.date )

states = unique(brazil_variables$State)
states <- sort(states)

#merge all tables together

onset_to_death_epoch_1 = fread('fitting_outputs/onset-to-death_epoch_1_br-samples-gamma.csv')
onset_to_death_epoch_2 = fread('fitting_outputs/onset-to-death_epoch_2_br-samples-gamma.csv')
onset_to_death_epoch_3 = fread('fitting_outputs/onset-to-death_epoch_3_br-samples-gamma.csv')

onset_to_death_epoch_1 <- onset_to_death_epoch_1[,1:27] / onset_to_death_epoch_1[,28:54] 
onset_to_death_epoch_2 <- onset_to_death_epoch_2[,1:27] / onset_to_death_epoch_2[,28:54] 
onset_to_death_epoch_3 <- onset_to_death_epoch_3[,1:27] / onset_to_death_epoch_3[,28:54] 

#get 

onset_states_1 <- onset_to_death_epoch_1[,c(1:27)]

colnames(onset_states_1) <- states

means_br_onset_death_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:27) 

onset_states_2 <- onset_to_death_epoch_2[,c(1:27)]

colnames(onset_states_2) <- states

means_br_onset_death_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:27) 

onset_states_3 <- onset_to_death_epoch_3[,c(1:27)]

colnames(onset_states_3) <- states

means_br_onset_death_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:27) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_br <- rbind(means_br_onset_death_epoch1,means_br_onset_death_epoch2,means_br_onset_death_epoch3)

#plot

br_death_epoch_1 <- ggplot(combined) +
  aes(x=epoch_1, y=reorder(State,-epoch_1)) +
  geom_posterior(
    aes(color="orange"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-death') +
  xlim(5,25)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

br_death_epoch_1 <- br_death_epoch_1 +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
br_death_epoch_1 <- br_death_epoch_1 +   geom_posterior(
  aes(x=epoch_3, y=State, color="blue"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)



br_death_epoch_1 <- br_death_epoch_1 + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                       values=c("orange"="orange","purple"="purple","blue"="blue"),
                       labels = c("Epoch 1", "Epoch 2", "Epoch 3"))
#hospitalization to death

#merge all tables together

hospitalisation_to_death_epoch_1 = fread('fitting_outputs/hospitalisation_to_death_epoch_1_br-samples-gamma.csv')
hospitalisation_to_death_epoch_2 = fread('fitting_outputs/hospitalisation_to_death_epoch_2_br-samples-gamma.csv')
hospitalisation_to_death_epoch_3 = fread('fitting_outputs/hospitalisation_to_death_epoch_3_br-samples-gamma.csv')

hospitalisation_to_death_epoch_1 <- hospitalisation_to_death_epoch_1[,1:27] / hospitalisation_to_death_epoch_1[,28:54] 
hospitalisation_to_death_epoch_2 <- hospitalisation_to_death_epoch_2[,1:27] / hospitalisation_to_death_epoch_2[,28:54] 
hospitalisation_to_death_epoch_3 <- hospitalisation_to_death_epoch_3[,1:27] / hospitalisation_to_death_epoch_3[,28:54] 

#get 
onset_states_1 <- hospitalisation_to_death_epoch_1[,c(1:27)]

colnames(onset_states_1) <- states

means_br_hosp_death_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:27) 

onset_states_2 <- hospitalisation_to_death_epoch_2[,c(1:27)]

colnames(onset_states_2) <- states

means_br_hosp_death_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:27) 

onset_states_3 <- hospitalisation_to_death_epoch_3[,c(1:27)]

colnames(onset_states_3) <- states

means_br_hosp_death_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:27) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_br <- rbind(mean_br,means_br_hosp_death_epoch1,means_br_hosp_death_epoch2,means_br_hosp_death_epoch3)

#plot

br_hosp_death_epoch_1 <- ggplot(combined) +
  aes(x=epoch_1, y=reorder(State,-epoch_1)) +
  geom_posterior(
    aes(color="orange"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from hospitalisation-to-death') +
  xlim(5,25)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

br_hosp_death_epoch_1 <- br_hosp_death_epoch_1 +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
br_hosp_death_epoch_1 <- br_hosp_death_epoch_1 +   geom_posterior(
  aes(x=epoch_3, y=State, color="blue"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

br_hosp_death_epoch_1 <- br_hosp_death_epoch_1 + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                       values=c("orange"="orange","purple"="purple","blue"="blue"),
                       labels = c("Epoch 1", "Epoch 2", "Epoch 3"))


#Onset to Diagnosis 

onset_to_diagnosis_epoch_1 = fread('fitting_outputs/onset-to-diagnosis_epoch_1_br-samples-lognormal.csv')
onset_to_diagnosis_epoch_2 = fread('fitting_outputs/onset-to-diagnosis_epoch_2_br-samples-lognormal.csv')
onset_to_diagnosis_epoch_3 = fread('fitting_outputs/onset-to-diagnosis_epoch_3_br-samples-lognormal.csv')

onset_to_diagnosis_epoch_1 <- exp((onset_to_diagnosis_epoch_1[,1:27]+0.5)*(onset_to_diagnosis_epoch_1[,28:54]^2)) 
onset_to_diagnosis_epoch_2 <- exp((onset_to_diagnosis_epoch_2[,1:27]+0.5)*(onset_to_diagnosis_epoch_2[,28:54]^2)) 
onset_to_diagnosis_epoch_3 <- exp((onset_to_diagnosis_epoch_3[,1:27]+0.5)*(onset_to_diagnosis_epoch_3[,28:54]^2)) 

#get 
onset_states_1 <- onset_to_diagnosis_epoch_1[,c(1:27)]

colnames(onset_states_1) <- states

means_br_onset_diagnosis_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:27) 

onset_states_2 <- onset_to_diagnosis_epoch_2[,c(1:27)]

colnames(onset_states_2) <- states

means_br_onset_diagnosis_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:27) 

onset_states_3 <- onset_to_diagnosis_epoch_3[,c(1:27)]

colnames(onset_states_3) <- states

means_br_onset_diagnosis_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:27) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_br <- rbind(mean_br,means_br_onset_diagnosis_epoch1,means_br_onset_diagnosis_epoch2,means_br_onset_diagnosis_epoch3)

#plot

br_diagnosis <- ggplot(combined) +
  aes(x=epoch_3, y=reorder(State,-epoch_3)) +
  geom_posterior(
    aes(color="blue"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-diagnosis') +
  xlim(0,50)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

br_diagnosis <- br_diagnosis +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
br_diagnosis <- br_diagnosis +   geom_posterior(
  aes(x=epoch_1, y=State, color="orange"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

br_diagnosis <- br_diagnosis + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                       values=c("orange"="orange","purple"="purple","blue"="blue"),
                       labels = c("Epoch 1", "Epoch 2", "Epoch 3"))


#Onset to hospitalisation  

onset_to_hospitalisation_epoch_1 = fread('fitting_outputs/onset-to-hospitalisation_epoch_1_br-samples-lognormal.csv')
onset_to_hospitalisation_epoch_2 = fread('fitting_outputs/onset-to-hospitalisation_epoch_2_br-samples-lognormal.csv')
onset_to_hospitalisation_epoch_3 = fread('fitting_outputs/onset-to-hospitalisation_epoch_3_br-samples-lognormal.csv')

onset_to_hospitalisation_epoch_1 <- exp((onset_to_hospitalisation_epoch_1[,1:27]+0.5)*(onset_to_hospitalisation_epoch_1[,28:54]^2)) 
onset_to_hospitalisation_epoch_2 <- exp((onset_to_hospitalisation_epoch_2[,1:27]+0.5)*(onset_to_hospitalisation_epoch_2[,28:54]^2)) 
onset_to_hospitalisation_epoch_3 <- exp((onset_to_hospitalisation_epoch_3[,1:27]+0.5)*(onset_to_hospitalisation_epoch_3[,28:54]^2)) 

#get 
onset_states_1 <- onset_to_hospitalisation_epoch_1[,c(1:27)]

colnames(onset_states_1) <- states

means_br_onset_hosp_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:27) 

onset_states_2 <- onset_to_hospitalisation_epoch_2[,c(1:27)]

colnames(onset_states_2) <- states

means_br_onset_hosp_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:27) 

onset_states_3 <- onset_to_hospitalisation_epoch_3[,c(1:27)]

colnames(onset_states_3) <- states

means_br_onset_hosp_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:27) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_br <- rbind(mean_br,means_br_onset_hosp_epoch1,means_br_onset_hosp_epoch2,means_br_onset_hosp_epoch3)

#plot

br_hospitalisation <- ggplot(combined) +
  aes(x=epoch_3, y=reorder(State,-epoch_3)) +
  geom_posterior(
    aes(color="blue"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-hospitalisation') +
  xlim(0,25)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

br_hospitalisation <- br_hospitalisation +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
br_hospitalisation <- br_hospitalisation +   geom_posterior(
  aes(x=epoch_1, y=State, color="orange"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

br_hospitalisation <- br_hospitalisation + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                       values=c("orange"="orange","purple"="purple","blue"="blue"),
                       labels = c("Epoch 1", "Epoch 2", "Epoch 3"))

ggarrange(
  br_diagnosis,br_hospitalisation,br_hosp_death_epoch_1,br_death_epoch_1, labels = c("A", "B", 'C','D'),
  common.legend = TRUE
)

write.csv(mean_br, file = "results/mean_br.csv", row.names = TRUE)

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

###Onset to death###


mexico_variables <-  mexico %>% 
  summarise('onset-to-confirmed' = delta_symptom_onset_to_confirmed,
            'onset-to-death' = delta_symptom_onset_to_death,
            'onset-to-hospitalisation' = delta_symptom_onset_to_hospitalisation,
            State = ADM1_ES,
            Date =  events.confirmed.date )

states = unique(mexico_variables$State)
states <- sort(states)

#merge all tables together

onset_to_death_epoch_1 = fread('fitting_outputs/onset-to-death_epoch_1_mex-samples-gamma.csv')
onset_to_death_epoch_2 = fread('fitting_outputs/onset-to-death_epoch_2_mex-samples-gamma.csv')
onset_to_death_epoch_3 = fread('fitting_outputs/onset-to-death_epoch_3_mex-samples-gamma.csv')

onset_to_death_epoch_1 <- onset_to_death_epoch_1[,1:32] / onset_to_death_epoch_1[,33:64] 
onset_to_death_epoch_2 <- onset_to_death_epoch_2[,1:32] / onset_to_death_epoch_2[,33:64] 
onset_to_death_epoch_3 <- onset_to_death_epoch_3[,1:32] / onset_to_death_epoch_3[,33:64] 

#get 

onset_states_1 <- onset_to_death_epoch_1[,c(1:32)]

colnames(onset_states_1) <- states

means_mex_onset_death_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:32) 

onset_states_2 <- onset_to_death_epoch_2[,c(1:32)]

colnames(onset_states_2) <- states

means_mex_onset_death_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:32) 

onset_states_3 <- onset_to_death_epoch_3[,c(1:32)]

colnames(onset_states_3) <- states

means_mex_onset_death_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:32) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_mex <- rbind(means_mex_onset_death_epoch1,means_mex_onset_death_epoch2,means_mex_onset_death_epoch3)

#plot

mex_death_epoch_1 <- ggplot(combined) +
  aes(x=epoch_1, y=reorder(State,-epoch_1)) +
  geom_posterior(
    aes(color="orange"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-death') +
  xlim(5,20)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

mex_death_epoch_1 <- mex_death_epoch_1 +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
mex_death_epoch_1 <- mex_death_epoch_1 +   geom_posterior(
  aes(x=epoch_3, y=State, color="blue"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)



mex_death_epoch_1 <- mex_death_epoch_1 + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                          values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                          labels = c("Epoch 1", "Epoch 2", "Epoch 3"))
#hospitalization to death

#merge all tables together

hospitalisation_to_death_epoch_1 = fread('fitting_outputs/hospitalisation_to_death_epoch_1_mex-samples-gamma.csv')
hospitalisation_to_death_epoch_2 = fread('fitting_outputs/hospitalisation_to_death_epoch_2_mex-samples-gamma.csv')
hospitalisation_to_death_epoch_3 = fread('fitting_outputs/hospitalisation_to_death_epoch_3_mex-samples-gamma.csv')

hospitalisation_to_death_epoch_1 <- hospitalisation_to_death_epoch_1[,1:32] / hospitalisation_to_death_epoch_1[,33:64] 
hospitalisation_to_death_epoch_2 <- hospitalisation_to_death_epoch_2[,1:32] / hospitalisation_to_death_epoch_2[,33:64] 
hospitalisation_to_death_epoch_3 <- hospitalisation_to_death_epoch_3[,1:32] / hospitalisation_to_death_epoch_3[,33:64] 

#get 
onset_states_1 <- hospitalisation_to_death_epoch_1[,c(1:32)]

colnames(onset_states_1) <- states

means_mex_hosp_death_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:32) 

onset_states_2 <- hospitalisation_to_death_epoch_2[,c(1:32)]

colnames(onset_states_2) <- states

means_mex_hosp_death_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:32) 

onset_states_3 <- hospitalisation_to_death_epoch_3[,c(1:32)]

colnames(onset_states_3) <- states

means_mex_hosp_death_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:32) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_mex <- rbind(mean_mex,means_mex_hosp_death_epoch1,means_mex_hosp_death_epoch2,means_mex_hosp_death_epoch3)

#plot

mex_hosp_death_epoch_1 <- ggplot(combined) +
  aes(x=epoch_1, y=reorder(State,-epoch_1)) +
  geom_posterior(
    aes(color="orange"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from hospitalisation-to-death') +
  xlim(5,20)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

mex_hosp_death_epoch_1 <- mex_hosp_death_epoch_1 +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
mex_hosp_death_epoch_1 <- mex_hosp_death_epoch_1 +   geom_posterior(
  aes(x=epoch_3, y=State, color="blue"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

mex_hosp_death_epoch_1 <- mex_hosp_death_epoch_1 + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                                    values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                                    labels = c("Epoch 1", "Epoch 2", "Epoch 3"))


#Onset to Diagnosis 

onset_to_diagnosis_epoch_1 = fread('fitting_outputs/onset-to-diagnosis_epoch_1_mex-samples-lognormal.csv')
onset_to_diagnosis_epoch_2 = fread('fitting_outputs/onset-to-diagnosis_epoch_2_mex-samples-lognormal.csv')
onset_to_diagnosis_epoch_3 = fread('fitting_outputs/onset-to-diagnosis_epoch_3_mex-samples-lognormal.csv')

onset_to_diagnosis_epoch_1 <- exp((onset_to_diagnosis_epoch_1[,1:32]+0.5)*(onset_to_diagnosis_epoch_1[,33:64]^2)) 
onset_to_diagnosis_epoch_2 <- exp((onset_to_diagnosis_epoch_2[,1:32]+0.5)*(onset_to_diagnosis_epoch_2[,33:64]^2)) 
onset_to_diagnosis_epoch_3 <- exp((onset_to_diagnosis_epoch_3[,1:32]+0.5)*(onset_to_diagnosis_epoch_3[,33:64]^2)) 

#get 
onset_states_1 <- onset_to_diagnosis_epoch_1[,c(1:32)]

colnames(onset_states_1) <- states

means_mex_onset_diagnosis_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:32) 

onset_states_2 <- onset_to_diagnosis_epoch_2[,c(1:32)]

colnames(onset_states_2) <- states

means_mex_onset_diagnosis_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:32) 

onset_states_3 <- onset_to_diagnosis_epoch_3[,c(1:32)]

colnames(onset_states_3) <- states

means_mex_onset_diagnosis_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:32) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_mex <- rbind(mean_mex,means_mex_onset_diagnosis_epoch1,means_mex_onset_diagnosis_epoch2,means_mex_onset_diagnosis_epoch3)

#plot

mex_diagnosis <- ggplot(combined) +
  aes(x=epoch_3, y=reorder(State,-epoch_3)) +
  geom_posterior(
    aes(color="blue"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-diagnosis') +
  xlim(0,20)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

mex_diagnosis <- mex_diagnosis +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
mex_diagnosis <- mex_diagnosis +   geom_posterior(
  aes(x=epoch_1, y=State, color="orange"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

mex_diagnosis <- mex_diagnosis + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                  values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                  labels = c("Epoch 1", "Epoch 2", "Epoch 3"))


#Onset to hospitalisation  

onset_to_hospitalisation_epoch_1 = fread('fitting_outputs/onset-to-hospitalisation_epoch_1_mex-samples-lognormal.csv')
onset_to_hospitalisation_epoch_2 = fread('fitting_outputs/onset-to-hospitalisation_epoch_2_mex-samples-lognormal.csv')
onset_to_hospitalisation_epoch_3 = fread('fitting_outputs/onset-to-hospitalisation_epoch_3_mex-samples-lognormal.csv')

onset_to_hospitalisation_epoch_1 <- exp((onset_to_hospitalisation_epoch_1[,1:32]+0.5)*(onset_to_hospitalisation_epoch_1[,33:64]^2)) 
onset_to_hospitalisation_epoch_2 <- exp((onset_to_hospitalisation_epoch_2[,1:32]+0.5)*(onset_to_hospitalisation_epoch_2[,33:64]^2)) 
onset_to_hospitalisation_epoch_3 <- exp((onset_to_hospitalisation_epoch_3[,1:32]+0.5)*(onset_to_hospitalisation_epoch_3[,33:64]^2)) 

#get 
onset_states_1 <- onset_to_hospitalisation_epoch_1[,c(1:32)]

colnames(onset_states_1) <- states

means_mex_onset_hosp_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:32) 

onset_states_2 <- onset_to_hospitalisation_epoch_2[,c(1:32)]

colnames(onset_states_2) <- states

means_mex_onset_hosp_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:32) 

onset_states_3 <- onset_to_hospitalisation_epoch_3[,c(1:32)]

colnames(onset_states_3) <- states

means_mex_onset_hosp_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:32) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_mex <- rbind(mean_mex,means_mex_onset_hosp_epoch1,means_mex_onset_hosp_epoch2,means_mex_onset_hosp_epoch3)

#plot

mex_hospitalisation <- ggplot(combined) +
  aes(x=epoch_3, y=reorder(State,-epoch_3)) +
  geom_posterior(
    aes(color="blue"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-hospitalisation') +
  xlim(0,20)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

mex_hospitalisation <- mex_hospitalisation +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
mex_hospitalisation <- mex_hospitalisation +   geom_posterior(
  aes(x=epoch_1, y=State, color="orange"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

mex_hospitalisation <- mex_hospitalisation + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                              values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                              labels = c("Epoch 1", "Epoch 2", "Epoch 3"))

ggarrange(
  mex_diagnosis,mex_hospitalisation,mex_hosp_death_epoch_1,mex_death_epoch_1, labels = c("A", "B", 'C','D'),
  common.legend = TRUE
)

write.csv(mean_mex, file = "results/mean_mex.csv", row.names = TRUE)

#Argentina

#Load data

argentina = as.data.frame(read_csv("data/argentina_linelist_geocoded_adm1_20220131.csv"))

#Onset to Diagnosis 

argentina_variables <-  argentina %>% 
  summarise('onset-to-confirmed' = delta_symptom_onset_to_confirmed,
            'onset-to-death' = delta_symptom_onset_to_death,
            'onset-to-hospitalisation' = delta_symptom_onset_to_hospitalisation,
            State = ADM1_ES,
            Date =  events.confirmed.date )

#Onset to Diagnosis 

onset_to_diagnosis_epoch_1 = fread('fitting_outputs/onset-to-diagnosis_epoch_1_ar-samples-lognormal.csv')
onset_to_diagnosis_epoch_2 = fread('fitting_outputs/onset-to-diagnosis_epoch_2_ar-samples-lognormal.csv')
onset_to_diagnosis_epoch_3 = fread('fitting_outputs/onset-to-diagnosis_epoch_3_ar-samples-lognormal.csv')

onset_to_diagnosis_epoch_1 <- exp((onset_to_diagnosis_epoch_1[,1:20]+0.5)*(onset_to_diagnosis_epoch_1[,21:40]^2)) 
onset_to_diagnosis_epoch_2 <- exp((onset_to_diagnosis_epoch_2[,1:20]+0.5)*(onset_to_diagnosis_epoch_2[,21:40]^2)) 
onset_to_diagnosis_epoch_3 <- exp((onset_to_diagnosis_epoch_3[,1:20]+0.5)*(onset_to_diagnosis_epoch_3[,21:40]^2)) 

states = unique(argentina_variables$State)
states <- sort(states)

#get 
onset_states_1 <- onset_to_diagnosis_epoch_1[,c(1:20)]

colnames(onset_states_1) <- states

means_ar_onset_diagnosis_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:20) 

onset_states_2 <- onset_to_diagnosis_epoch_2[,c(1:20)]

colnames(onset_states_2) <- states

means_ar_onset_diagnosis_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:20) 

onset_states_3 <- onset_to_diagnosis_epoch_3[,c(1:20)]

colnames(onset_states_3) <- states

means_ar_onset_diagnosis_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:20) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_ar_diagnosis <- rbind(means_ar_onset_diagnosis_epoch1,means_ar_onset_diagnosis_epoch2,means_ar_onset_diagnosis_epoch3)

#plot

ar_diagnosis <- ggplot(combined) +
  aes(x=epoch_3, y=reorder(State,-epoch_3)) +
  geom_posterior(
    aes(color="blue"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-diagnosis') +
  xlim(0,25)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

ar_diagnosis <- ar_diagnosis +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
ar_diagnosis <- ar_diagnosis +   geom_posterior(
  aes(x=epoch_1, y=State, color="orange"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

ar_diagnosis <- ar_diagnosis + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                  values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                  labels = c("Epoch 1", "Epoch 2", "Epoch 3"))

#onset_to_death

death_state_exclude <- c('La Pampa', 'Salta', 'Santa Cruz', 'Formosa', 'Tucumán', 'San Luis')

argentina_variables_death <- filter(argentina_variables, !State %in% death_state_exclude)


states = unique(argentina_variables_death$State)
states <- sort(states)

#merge all tables together

onset_to_death_epoch_1 = fread('fitting_outputs/onset-to-death_epoch_1_ar-samples-gamma.csv')
onset_to_death_epoch_2 = fread('fitting_outputs/onset-to-death_epoch_2_ar-samples-gamma.csv')
onset_to_death_epoch_3 = fread('fitting_outputs/onset-to-death_epoch_3_ar-samples-gamma.csv')

onset_to_death_epoch_1 <- onset_to_death_epoch_1[,1:14] / onset_to_death_epoch_1[,15:28] 
onset_to_death_epoch_2 <- onset_to_death_epoch_2[,1:14] / onset_to_death_epoch_2[,15:28] 
onset_to_death_epoch_3 <- onset_to_death_epoch_3[,1:14] / onset_to_death_epoch_3[,15:28] 

#get 

onset_states_1 <- onset_to_death_epoch_1[,c(1:14)]

colnames(onset_states_1) <- states

means_ar_onset_death_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:14) 

onset_states_2 <- onset_to_death_epoch_2[,c(1:14)]

colnames(onset_states_2) <- states

means_ar_onset_death_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:14) 

onset_states_3 <- onset_to_death_epoch_3[,c(1:14)]

colnames(onset_states_3) <- states

means_ar_onset_death_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:14) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_ar_death <-rbind(means_ar_onset_death_epoch1,means_ar_onset_death_epoch2,means_ar_onset_death_epoch3)

#plot

ar_death_epoch_1 <- ggplot(combined) +
  aes(x=epoch_1, y=reorder(State,-epoch_1)) +
  geom_posterior(
    aes(color="orange"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-death') +
  xlim(5,25)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

ar_death_epoch_1 <- ar_death_epoch_1 +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
ar_death_epoch_1 <- ar_death_epoch_1 +   geom_posterior(
  aes(x=epoch_3, y=State, color="blue"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)



ar_death_epoch_1 <- ar_death_epoch_1 + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                            values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                            labels = c("Epoch 1", "Epoch 2", "Epoch 3"))
#hospitalization to death

#merge all tables together

hospitalisation_to_death_epoch_1 = fread('fitting_outputs/hospitalisation_to_death_epoch_1_ar-samples-gamma.csv')
hospitalisation_to_death_epoch_2 = fread('fitting_outputs/hospitalisation_to_death_epoch_2_ar-samples-gamma.csv')
hospitalisation_to_death_epoch_3 = fread('fitting_outputs/hospitalisation_to_death_epoch_3_ar-samples-gamma.csv')

hospitalisation_to_death_epoch_1 <- hospitalisation_to_death_epoch_1[,1:32] / hospitalisation_to_death_epoch_1[,33:64] 
hospitalisation_to_death_epoch_2 <- hospitalisation_to_death_epoch_2[,1:32] / hospitalisation_to_death_epoch_2[,33:64] 
hospitalisation_to_death_epoch_3 <- hospitalisation_to_death_epoch_3[,1:32] / hospitalisation_to_death_epoch_3[,33:64] 

#get 
onset_states_1 <- hospitalisation_to_death_epoch_1[,c(1:32)]

colnames(onset_states_1) <- states

means_ar_hosp_death_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:32) 

onset_states_2 <- hospitalisation_to_death_epoch_2[,c(1:32)]

colnames(onset_states_2) <- states

means_ar_hosp_death_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:32) 

onset_states_3 <- hospitalisation_to_death_epoch_3[,c(1:32)]

colnames(onset_states_3) <- states

means_ar_hosp_death_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:32) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_ar <- rbind(mean_ar,means_ar_hosp_death_epoch1,means_ar_hosp_death_epoch2,means_ar_hosp_death_epoch3)

#plot

ar_hosp_death_epoch_1 <- ggplot(combined) +
  aes(x=epoch_1, y=reorder(State,-epoch_1)) +
  geom_posterior(
    aes(color="orange"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from hospitalisation-to-death') +
  xlim(5,20)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

ar_hosp_death_epoch_1 <- ar_hosp_death_epoch_1 +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
ar_hosp_death_epoch_1 <- ar_hosp_death_epoch_1 +   geom_posterior(
  aes(x=epoch_3, y=State, color="blue"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

ar_hosp_death_epoch_1 <- ar_hosp_death_epoch_1 + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                                      values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                                      labels = c("Epoch 1", "Epoch 2", "Epoch 3"))


#Onset to hospitalisation  

hosp_state_exclude <- c('La Pampa', 'Salta', 'Santa Cruz', 'Formosa', 'Tucumán', 'San Luis','Chubut','Misiones','Chaco')

argentina_variables_hosp <- filter(argentina_variables, !State %in% hosp_state_exclude)


states = unique(argentina_variables_hosp$State)
states <- sort(states)

onset_to_hospitalisation_epoch_1 = fread('fitting_outputs/onset-to-hospitalisation_epoch_1_ar-samples-lognormal.csv')
onset_to_hospitalisation_epoch_2 = fread('fitting_outputs/onset-to-hospitalisation_epoch_2_ar-samples-lognormal.csv')
onset_to_hospitalisation_epoch_3 = fread('fitting_outputs/onset-to-hospitalisation_epoch_3_ar-samples-lognormal.csv')

onset_to_hospitalisation_epoch_1 <- exp((onset_to_hospitalisation_epoch_1[,1:11]+0.5)*(onset_to_hospitalisation_epoch_1[,12:22]^2)) 
onset_to_hospitalisation_epoch_2 <- exp((onset_to_hospitalisation_epoch_2[,1:11]+0.5)*(onset_to_hospitalisation_epoch_2[,12:22]^2)) 
onset_to_hospitalisation_epoch_3 <- exp((onset_to_hospitalisation_epoch_3[,1:11]+0.5)*(onset_to_hospitalisation_epoch_3[,12:22]^2)) 

#get 
onset_states_1 <- onset_to_hospitalisation_epoch_1[,c(1:11)]

colnames(onset_states_1) <- states

means_ar_onset_hosp_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:11) 

onset_states_2 <- onset_to_hospitalisation_epoch_2[,c(1:11)]

colnames(onset_states_2) <- states

means_ar_onset_hosp_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:11) 

onset_states_3 <- onset_to_hospitalisation_epoch_3[,c(1:11)]

colnames(onset_states_3) <- states

means_ar_onset_hosp_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:11) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_ar_hosp <- rbind(means_ar_onset_hosp_epoch1,means_ar_onset_hosp_epoch2,means_ar_onset_hosp_epoch3)

#plot

ar_hospitalisation <- ggplot(combined) +
  aes(x=epoch_3, y=reorder(State,-epoch_3)) +
  geom_posterior(
    aes(color="blue"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-hospitalisation') +
  xlim(0,25)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

ar_hospitalisation <- ar_hospitalisation +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
ar_hospitalisation <- ar_hospitalisation +   geom_posterior(
  aes(x=epoch_1, y=State, color="orange"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

ar_hospitalisation <- ar_hospitalisation + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                                values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                                labels = c("Epoch 1", "Epoch 2", "Epoch 3"))

ggarrange(
  ar_diagnosis,ar_hospitalisation,ar_death_epoch_1, labels = c("A", "B", 'C'),
  common.legend = TRUE
)

write.csv(mean_ar_diagnosis, file = "results/mean_ar_diagnosis.csv", row.names = TRUE)
write.csv(mean_ar_hosp, file = "results/mean_ar_hosp.csv", row.names = TRUE)
write.csv(mean_ar_death, file = "results/mean_ar_death.csv", row.names = TRUE)

#Colombia

#Load data

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

#Onset to Diagnosis 

colombia_variables <-  colombia %>% 
  summarise('onset-to-confirmed' = delta_symptom_onset_to_confirmed,
            'onset-to-death' = delta_symptom_onset_to_death,
            'onset-to-hospitalisation' = delta_symptom_onset_to_hospitalisation,
            State = ADM1_ES,
            Date =  events.confirmed.date )

#Onset to Diagnosis 

onset_to_diagnosis_epoch_1 = fread('fitting_outputs/onset-to-diagnosis_epoch_1_co-samples-lognormal.csv')
onset_to_diagnosis_epoch_2 = fread('fitting_outputs/onset-to-diagnosis_epoch_2_co-samples-lognormal.csv')
onset_to_diagnosis_epoch_3 = fread('fitting_outputs/onset-to-diagnosis_epoch_3_co-samples-lognormal.csv')

onset_to_diagnosis_epoch_1 <- exp((onset_to_diagnosis_epoch_1[,1:30]+0.5)*(onset_to_diagnosis_epoch_1[,31:60]^2)) 
onset_to_diagnosis_epoch_2 <- exp((onset_to_diagnosis_epoch_2[,1:30]+0.5)*(onset_to_diagnosis_epoch_2[,31:60]^2)) 
onset_to_diagnosis_epoch_3 <- exp((onset_to_diagnosis_epoch_3[,1:30]+0.5)*(onset_to_diagnosis_epoch_3[,31:60]^2)) 

states = unique(columbia_epoch1_onset_to_diagnosis$State)
states <- sort(states)

#get 
onset_states_1 <- onset_to_diagnosis_epoch_1[,c(1:30)]

colnames(onset_states_1) <- states

means_co_onset_diagnosis_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:30) 

onset_states_2 <- onset_to_diagnosis_epoch_2[,c(1:30)]

colnames(onset_states_2) <- states

means_co_onset_diagnosis_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:30) 

onset_states_3 <- onset_to_diagnosis_epoch_3[,c(1:30)]

colnames(onset_states_3) <- states

means_co_onset_diagnosis_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:30) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_co_diagnosis <- rbind(means_co_onset_diagnosis_epoch1,means_co_onset_diagnosis_epoch2,means_co_onset_diagnosis_epoch3)

#plot

co_diagnosis <- ggplot(combined) +
  aes(x=epoch_3, y=reorder(State,-epoch_3)) +
  geom_posterior(
    aes(color="blue"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-diagnosis') +
  xlim(0,75)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

co_diagnosis <- co_diagnosis +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
co_diagnosis <- co_diagnosis +   geom_posterior(
  aes(x=epoch_1, y=State, color="orange"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)

co_diagnosis <- co_diagnosis + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                  values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                  labels = c("Epoch 1", "Epoch 2", "Epoch 3"))

#onset_to_death



states = unique(columbia_epoch1_onset_to_death$State)
states <- sort(states)

#merge all tables together

onset_to_death_epoch_1 = fread('fitting_outputs/onset-to-death_epoch_1_co-samples-gamma.csv')
onset_to_death_epoch_2 = fread('fitting_outputs/onset-to-death_epoch_2_co-samples-gamma.csv')
onset_to_death_epoch_3 = fread('fitting_outputs/onset-to-death_epoch_3_co-samples-gamma.csv')

onset_to_death_epoch_1 <- onset_to_death_epoch_1[,1:28] / onset_to_death_epoch_1[,29:56] 
onset_to_death_epoch_2 <- onset_to_death_epoch_2[,1:28] / onset_to_death_epoch_2[,29:56] 
onset_to_death_epoch_3 <- onset_to_death_epoch_3[,1:28] / onset_to_death_epoch_3[,29:56] 

#get 

onset_states_1 <- onset_to_death_epoch_1[,c(1:28)]

colnames(onset_states_1) <- states

means_co_onset_death_epoch1 <- colMeans(onset_states_1)

onset_states_1 <- onset_states_1 %>% 
  pivot_longer(cols = 1:28) 

onset_states_2 <- onset_to_death_epoch_2[,c(1:28)]

colnames(onset_states_2) <- states

means_co_onset_death_epoch2 <- colMeans(onset_states_2)

onset_states_2 <- onset_states_2 %>% 
  pivot_longer(cols = 1:28) 

onset_states_3 <- onset_to_death_epoch_3[,c(1:28)]

colnames(onset_states_3) <- states

means_co_onset_death_epoch3 <- colMeans(onset_states_3)

onset_states_3 <- onset_states_3 %>% 
  pivot_longer(cols = 1:28) 

combined <- cbind(onset_states_1,onset_states_2,onset_states_3)

colnames(combined) <- c("State","epoch_1","rm1","epoch_2","rm2","epoch_3")

combined <- dplyr:: select (combined, c(1,2,4,6))

#combine means 

mean_co_death <-rbind(means_co_onset_death_epoch1,means_co_onset_death_epoch2,means_co_onset_death_epoch3)

#plot

co_death_epoch_1 <- ggplot(combined) +
  aes(x=epoch_1, y=reorder(State,-epoch_1)) +
  geom_posterior(
    aes(color="orange"),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="State", x='Delay from symtom onset-to-death') +
  xlim(5,75)+
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right")

co_death_epoch_1 <- co_death_epoch_1 +
  geom_posterior(
    aes(x=epoch_2, y=State, color="purple" ),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  )
co_death_epoch_1 <- co_death_epoch_1 +   geom_posterior(
  aes(x=epoch_3, y=State, color="blue"),
  midline=NULL,
  mirror=TRUE,
  fill="#FFFFFF",
  draw_sd=FALSE,
  interval_type="hdi",
  vjust=0,
  position=position_spread(
    reverse=TRUE, # order of spreaded groups within panels
    padding=0.6, # shrink heights of distributions
    height=2 # scale by heights within panels
  ),
  adjust=1.5,
)



co_death_epoch_1 <- co_death_epoch_1 + scale_color_manual("Epoch", breaks=c("orange", "purple", "blue"),
                                                          values=c("orange"="orange","purple"="purple","blue"="blue"),
                                                          labels = c("Epoch 1", "Epoch 2", "Epoch 3"))
ggarrange(
  co_diagnosis,co_death_epoch_1, labels = c("A", "B"),
  common.legend = TRUE
)

write.csv(mean_co_diagnosis, file = "results/mean_co_diagnosis.csv", row.names = TRUE)
write.csv(mean_co_death, file = "results/mean_co_death.csv", row.names = TRUE)


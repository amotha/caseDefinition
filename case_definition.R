## PROGRAM: malaria case definition 
## PURPOSE: Determine parasitemia threshold for malaria case definition
#             To implement the method of Smith et. al 1994; 
#              "Attributable fraction estimates and case definitions for malaria in endemic areas." 

## WRITTEN BY: Amos Thairu
## DATE: June 2014
## UPDATED BY: Amos Thairu
## UPDATED: May 2015

## @knitr junjuParasitemiaThresholds.1
library(ggplot2)
library(glm2)
library(plyr)
library(boot)
library(dplyr)

#import the malaria dataset
#read in the functions in the script 'case_definition_functions'

#specify candidate parasite density thresholds
thresh=c(1,100,500,800,900,1000,1100,1200,1300,1500,1800,2000,2300,2500,2800,3000,3500,4000,5000,6000,7000,8000,9000)

#generate attributable fraction, sensitivity and specificity at the various candidate thresholds
result <- ddply(maldata, .(agecat), specSens)

result$sensitivity <- unlist(result$sensitivity)  
result$specificity <- unlist(result$specificity)
result$threshold <- unlist(result$threshold)
result$AF <- unlist(result$AF)

#plot the result;sensitivity and specificity plots over age categories 
#define breaks to use in the plot
seq.plot <- c(1,100,500,1000,2000,3000,4000,5000,6000,7000,8000,9000) 
result.plot <- result[result$threshold %in% seq.plot,]

ggplot(result.plot, aes(threshold,specificity,color = "Specificity")) + geom_line() +
  geom_line(aes(y=sensitivity,colour="Sensitivity")) +
  xlab(expression(paste("Parasites / ", mu,"L"))) + ylab("Percentage") +
  scale_colour_discrete(name="") +
  facet_wrap(~ agecat, ncol=2) +
  xlim(0, 7500) + ylim(70, 100) +
  theme(panel.background=element_blank(), panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(),
        axis.line=element_line(colour="light grey"), axis.text=element_text(colour='black'))


#bootstrap to calculate the confidence intervals for the MAFs; 1000 bootstrap samples
# boot.result <- boot.AF(maldata)
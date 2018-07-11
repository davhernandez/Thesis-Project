# setup ---------------------------

rm(list = ls())
library(dplyr)
library(ggplot2)

# CCCFRP dates --------------------------------------

#importing CCFRP data
sample_date <- read.csv("~/Desktop/Thesis/Raw Data/CCFRP/Hernandez_BlueRF.csv")
#selecting date for CCFRP
sample_date$Month <- match(sample_date$Month, month.name)
#a new column that is a character vector combining the month, day, and year columns into the full date
sample_date <- mutate(sample_date, Date = paste(sample_date$Month, sample_date$Day, sample_date$Year, sep = "/"))
#since the Date column is character vector, it needs to be converted to an object type Date
sample_date$Date <- as.Date(sample_date$Date, "%m/%d/%Y")
#add a column Frequency to tell the number of observations. The column will be filled with 1 because each row is for an individual fish
sample_date <- mutate(sample_date, Frequency = 1)
#add a column for the source of the data
sample_date <- mutate(sample_date, source = 'CCFRP')

# PISCO dates -----------------------

#importing PISCO data
smys <- read.csv("~/Desktop/Thesis/Raw Data/PISCO/UCSB_FISH.csv")
#filter out all fish that aren't BRF
smys = subset(smys, classcode == 'SMYS')
#a new column that is a character vector combining the month, day, and year columns into the full date
smys <- mutate(smys, Date = paste(smys$month, smys$day, smys$year, sep = "/"))
#since the Date column is character vector, it needs to be converted to an object type Date
smys$Date <- as.Date(smys$Date, "%m/%d/%Y")
#rename 'count' to 'Frequnecy' so that it matchs the column name of 'sample_date'
smys <- rename(smys, Frequency = count)
#add a column for the source of the data
smys <- mutate(smys, source = 'PISCO')

# combining data -----------------------

#select CCFRP dates and add a column for frequency. All frequencies = 1
#select PISCO dates and count column
#join both matricies
joined_dates <- rbind(sample_date[, c("Date", "Frequency", "source")], smys[,c("Date", "Frequency", "source")])
#stacked bar plot of when each fish was sampled
ggplot(joined_dates, aes(x = Date, y = Frequency, fill = source)) +
  geom_bar(stat='identity')
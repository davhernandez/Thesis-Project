# setup ---------------------------

rm(list = ls())
library(dplyr)
library(ggplot2)
library(lubridate)

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

# plotting just month and day ---------------------------
#this combines all of the data across years to look at the trend in what time of year the samples were taken
#this section uses lubridate package to reach its goal. The previous section did it in base R
sample_date <- mutate(sample_date, Month_Day = paste(month(sample_date$Month, label = TRUE), sample_date$Day, sep = "-"))

smys <- mutate(smys, Month_Day = paste(month(smys$month, label = TRUE), smys$day, sep = "-"))

joined_months <- rbind(sample_date[, c("Month_Day", "Frequency", "source")], smys[,c("Month_Day", "Frequency", "source")])

#collapsing all matching data points together
joined_months <- joined_months %>%
  na.omit %>%
  group_by(Month_Day, source) %>%
  summarise(Frequency = sum(Frequency))

#what if I mutate and extract the name of the month `month(data, label = TRUE)` and the day number and then forgo the as.Date?

ggplot(joined_months, aes(x = Month_Day, y = Frequency, fill=source)) +
  geom_bar(stat="identity")

#plotting just the months ---------------------------------------------------------------------

sample_date <- mutate(sample_date, Months = month(sample_date$Month, label = TRUE))

smys <- mutate(smys, Months = month(smys$month, label = TRUE))

joined_months <- rbind(sample_date[, c("Months", "Frequency", "source")], smys[,c("Months", "Frequency", "source")])

#collapsing all matching data points together
joined_months <- joined_months %>%
  na.omit %>%
  group_by(Months, source) %>%
  summarise(Frequency = sum(Frequency))

ggplot(joined_months, aes(x = Months, y = Frequency, fill=source)) +
  geom_bar(stat="identity") +
  ggtitle("Samples based on month collected")

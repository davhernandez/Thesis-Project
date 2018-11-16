#This script does the same exact thing as `Combined_data.R` as of 8/31/18. The only difference is that it only pulls the data for July, August, and September

# setup ---------------------------

#rm(list = ls())
library(dplyr)
library(geosphere)
library(ggplot2)
library(scales)
library(plotly) #for making contour plots
library(webshot) #for exporting plotly output
library(ggridges)

#the first part of this code is lifted from CCFRP_Mixed_CPUE.R
#The second part is lifted from PISCO_abundance.R
# CCCFRP abundance index --------------------------------------
ccfrp_raw <- read.csv("~/Desktop/Thesis/Raw Data/CCFRP/Hernandez_BlueRF.csv")
BothSpecies <- ccfrp_raw
#renaming the columns
BothSpecies <- rename(BothSpecies, fish_ID = Fish.ID, species_code = Species.Code, common_name = Common.Name, size = Length..cm., cell_ID = Cell.ID, site_type = Site..MPA..REF., start_depth = Start.Depth..ft., end_depth = End.Depth..ft., start_time = Start.Time, end_time = End.Time, drift_ID = Drift.ID)
# selecting the columns that will be used
BothSpecies <- select(BothSpecies, IDCell.per.Trip, size, Area, Year, cell_ID:end_time, drift_ID, species_code)
# this line was added to the CCFRP_CPUE file to get rid of NAs in the Lat-Lon columns
BothSpecies <- filter(BothSpecies, ST_LatDD != 'NA', ST_LonDD != 'NA', End_LatDD != 'NA', End_LonDD != 'NA')
#using distCosine from the geosphere package to calculate the distance of each drift based on the start and end lat/long. The lat/longs are concatinated with cbind() instead of just c() so that mutate works on single row at a time. The output is in meters.
BothSpecies <- mutate(BothSpecies, distance_meters = distCosine(p1 = cbind(ST_LonDD, ST_LatDD), p2 = cbind(End_LonDD, End_LatDD)))
#converting the start_time and end_time columns into POSIXct format, which is necessary for the difftime() function
BothSpecies$start_time <- as.POSIXct(BothSpecies$start_time, format = "%H:%M:%S")
BothSpecies$end_time <- as.POSIXct(BothSpecies$end_time, format = "%H:%M:%S")
#creating a new column (total_time) that is the total amount of time (in minutes) spent on the drift
BothSpecies <- mutate(BothSpecies, total_time = difftime(time1 = BothSpecies$end_time, time2 = BothSpecies$start_time, units = "mins"))
# eliminating the catch entries that have 'NA' for size
# the comma at the end of the [(...),] is important because you are defining the rows of BLU_only$length to queries fro NAs
BothSpecies <- BothSpecies[!is.na(BothSpecies$size),]
# below is another way of eliminating the NAs from the row 'size'
#BLU_only <- filter(BLU_only, length != 'NA')

# eliminating all entries that don't have a start or end depth
BothSpecies <- filter(BothSpecies, start_depth != 'NA', end_depth != "NA")
# convert depths from ft to m
BothSpecies$start_depth <- BothSpecies$start_depth * 0.3048
BothSpecies$end_depth <- BothSpecies$end_depth * 0.3048
# rounding down the depth columns
# this is important for calculating the correct dpeths for the Catch matrix
# see 'Catch' section for explanation of the rounding
BothSpecies$start_depth <- floor(BothSpecies$start_depth)
BothSpecies$end_depth <- floor(BothSpecies$end_depth)
#Filter for the months July, August, September
#the only way to filter with multiple character values is to put them in a variable and then use `%in%``
summer_months <- c("July", "August", "September")
BothSpecies <- filter(BothSpecies, Month %in% summer_months)
rm(list = "summer_months")

# Effort -------------------------------------------------------------------------------
if (exists("effort") == TRUE) {
  print("Effort calculation already exists in the environment")
} else {
  # importing the dataset that contains all drifts CCFRP has conducted
  all_drifts_raw <- read.csv("~/Desktop/Thesis/Raw Data/CCFRP/Hernandez_Drifts.csv")
  # loop for calculating the effort
  # this is a function so that the variables remain local instead of being stored globally
  effort_function <- function(){
    effort_data <- all_drifts_raw
    #filtering for just the summer months
    #similar to `summer_months`, I am creating a vector to match Drift.ID column. In this column, the first two numbers represent the month of the sample
    summer_integers <- c("07", "08", "09")
    effort_data <- filter(effort_data, substr(Drift.ID, 4,5) %in% summer_integers)
    
    # renmaing columns
    effort_data <- dplyr::rename(all_drifts_raw, start_time = Start.Time, end_time = End.Time, drift_time_hrs = Drift.time_hrs, total_andjusted_angler_hrs = total.adjusted.angler.hrs, drift_length = Drift.Length..m., start_depth = Start.Depth..ft., end_depth = End.Depth..ft.)
    # eliminating the depths entries that have 'NA'
    effort_data <- filter(effort_data, start_depth != 'NA', end_depth != 'NA')
    # converting start and end depth to meters
    # ththe meters conversion is combined with rounding the resulting number down to the nearest whole meter by using the floor() function
    # I am rounding down because I am binning per meter. So 30.2m counts as falling in the 30-31m bin
    effort_data$start_depth <- floor(effort_data$start_depth * 0.3048)
    effort_data$end_depth <- floor(effort_data$end_depth * 0.3048)
    # The point of this section is to create a vector that contains the range of depth values and add up how many angler hours were spent fishing at each depth.
    
    
    #effort_function(){ where it was originally placed. It was moved up to scope 'effort_data'
    
    # finding the minimum depth that was fished at
    min_depth_effort <- min(c(effort_data$start_depth, effort_data$end_depth))
    max_depth_effort <- max(c(effort_data$start_depth, effort_data$end_depth))
    
    #initialize effort matrix
    #making a new empty vector that is the length of the max depth - min depth
    effort <- rep(0, (max_depth_effort - min_depth_effort +1))
    # naming the vector so that entries can be input by the name/value of depth
    names(effort) <- c(min_depth_effort:max_depth_effort)
    
    for(x in 1:nrow(effort_data)){ #this will run through every row in the effort_data matrix
      # x counts the row of the effort_data matrix that is being currently worked on
      # determine if the start or end depth is deeper and save the deeper value
      # pmax stands for parallel maximum, which compares two columns and finds the maximum value
      deeper <- pmax(effort_data$start_depth[x], effort_data$end_depth[x])
      # determine if the start or end depth is shallower and save the shallower depth
      shallower <-pmin(effort_data$start_depth[x], effort_data$end_depth[x])
      # calls the position in the effort vector named for the depth range of the current row, then divides the total adjusted angler hours for that drift by the range of depths. This evenly splits the fishing time between all depths in the drift. The time is then input for each depth that was fished.
      # In order to call the correct name of the column the code [names(effort) %in% shallower:deeper] searches for the names in the vector effort that perfectly match the names in the series shallower:deeper
      effort[names(effort) %in% shallower:deeper] <- effort[names(effort) %in% shallower:deeper] + (effort_data$total_andjusted_angler_hrs[x]/(deeper- shallower + 1))
      # time keeper
      print(x)
    }
    #coercing it into a two column tidy dataframe
    tidy_effort <- data.frame(min_depth_effort:max_depth_effort, effort)
    colnames(tidy_effort) <- c("depth", "effort")
    #tidy_effort$effort <- effort
    # after the for loop finished, return the final effort vector
    return(tidy_effort)
  }
  # runs the effort_function and stores it in the effort vector
  effort <- effort_function()
}

#bar graph of effort
#ggplot(effort, aes(x = depth, weight = effort)) +
#  geom_bar() +
#  ggtitle("CCFRP sampling effort") +
#  ylab("Effort (minutes)") +
#  xlab("Depth")


# tidy function ---------------------------------

#control for CPUE that the site contributes
#matrix of CPUE per depth per site
# divide the CPUE by the row sum before adding it to total catch matrix

tidy <- function(catch, year=FALSE){
  
  # initialize variables
  # max and min depth
  min_depth <- min(c(catch$start_depth, catch$end_depth))
  max_depth <- max(c(catch$start_depth, catch$end_depth))
  # max and min size
  max_size <- max(catch$size)
  min_size <- min(catch$size)
  
  # initialize marticies
  # total CPUE
  total_CPUE <- matrix(0L, nrow = (max_depth - min_depth + 1), ncol = (max_size - min_size + 1))
  colnames(total_CPUE) <- c(dQuote(min_size:max_size))
  rownames(total_CPUE) <- c(dQuote(min_depth:max_depth))
  # construct total effort for each depth
  # effort <- effort_function()
  
  # tidy dataframe (site, size, depth, CPUE)
  tidy_data <- data.frame(row.names = c("site", "size", "depth", "catch", "year"))
  
  
  #this is the section of code that Duncan Lang vectorized in the Beyond R Basics workshop. It replaces the old loop.
  mn = pmin(catch$start_depth, catch$end_depth)
  mx = pmax(catch$start_depth, catch$end_depth)    
  nrep = mx - mn + 1
  tidy_data = data.frame(site = rep(substr(catch$IDCell.per.Trip, 0, 3), nrep),
                         size = rep(catch$size,  nrep),
                         depth = unlist(mapply(`:`, mn, mx, SIMPLIFY = FALSE)),
                         year = rep(catch$Year, nrep))
  tidy_data$catch = 1
  
  #Duncan added this line, which doesn't work because there is no `tmp` variable. If it did, it would call rbind for all parts of the list `tmp`. However the tidy_data df already is combined without the need for a do.call
  #tidy_data = do.call(rbind, tmp)
  
  # tallies the number of columns that match the 'group_by' call and collapses them down to one row
  if(year==TRUE){
    tidy_data <- tidy_data %>%
      group_by(site, size, depth, year) %>%
      tally()
  }else{
    tidy_data <- tidy_data %>%
      group_by(site, size, depth) %>%
      tally()
  }
  # matches the depth column in 'effort' to the depth column in 'tidy_data' and inputs the associated effort value as a new column in tidy-data
  tidy_data$effort <- effort[match(tidy_data$depth, effort$depth), 2]
  
  # Trying to control for the posibility that one site contributed the majority of the CPUE to the grand total, thereby masking the data from the other sites.
  # The way that it is being calculated here makes the assumption that, since they fished the same number of tansects, they fished for the same amount of time at each location
  # An ANOVA confirms that total drift times are not significantly different across sites
  # After calculating CPUE, total CPUE for each site is added up.
  # Then CPUE is divide by the site total
  
  # add a row that calculates CPUE based on total effort
  tidy_data <- mutate(tidy_data, CPUE = n/effort)
  # new variable that calculates the total CPUE per site
  CPUE_by_site <- tidy_data %>%
    group_by(site) %>%
    summarize(site_total_CPUE = sum(CPUE))
  # merge CPUE_by_site with tidy_data by matching the name 'site'
  tidy_data <- merge(tidy_data, CPUE_by_site, by='site')
  tidy_data <- mutate(tidy_data, normalized_CPUE = CPUE/site_total_CPUE)
  
  
  return(tidy_data)
}

#'tidy_data' will be used to form the conversion estimate
tidy_data <- tidy(catch = BothSpecies, year = FALSE)
#tidy_data_year' is what will be combined after the conversion estimate is calculated
tidy_data_year <- tidy(catch = BothSpecies, year = TRUE)

# CCFRP final matrix ------------------------------------------------------------------------------
# matrix of mixed BLU and DEA with uncontrolled CPUE
# summing CPUE across all depths
tidy_data <- aggregate(CPUE ~ size, data = tidy_data, FUN = sum)
tidy_data <- rename(tidy_data, count = CPUE)
tidy_data_year <- aggregate(CPUE ~ size + year + depth, data = tidy_data_year, FUN = sum)
tidy_data_year <- rename(tidy_data_year, count = CPUE)

# PISCO section -----------------------
# all the PISCO data was lifted from PISCO_abundance

#importing PISCO data
pisco_raw <- read.csv("~/Desktop/Thesis/Raw Data/PISCO/UCSB_FISH.csv")
#I'm pretty sure that `myst_only` is superfluous now
# selecting only the blue rockfish entries
#myst_only <- filter(pisco_raw, classcode == 'SMYS')
# filtering out rows with NA
#myst_only <- filter(myst_only, count != 'NA', fish_tl != 'NA', depth != 'NA')
# rounding down the tail length
#myst_only$fish_tl <- floor(myst_only$fish_tl)
# rounding down the depth
#myst_only$depth <- floor(myst_only$depth)
#group by year, site, side, zone, depth

# subsetting code with the help of Duncan Lang
d = pisco_raw

smys = subset(d, classcode == 'SMYS')
pisco_summer <- c(7,8,9)
smys <- filter(smys, month %in% pisco_summer)
smys$depth <- floor(smys$depth)
smys$length = factor(cut(smys$fish_tl, seq(0, 50, by = 1), labels = seq(0.5, 49.5, by = 1)))
bb = with(smys, by(smys, list(year, site, side, zone, depth),
                   function(d) {
                     tt = table(d$length)
                     ans = data.frame(count = as.integer(tt), length = names(tt))
                     cbind(ans, d[rep(1, nrow(ans)), c("year", "site", "side", "zone", "depth")])
                   }))

bb2 = bb[! sapply(bb, is.null) ]
bb3 = do.call(rbind, bb2)

bb4 = bb3[bb3$count > 0, ]
bb4$length <- floor(as.numeric(levels(bb4$length))[bb4$length])
#rounds it back down from intervals of 15.5, 16.5, 17.5, etc. to integers
bb5 <- aggregate(count ~ depth + length + year, data = bb4, FUN = sum)

# sampling effort -------------------------------------------------

#counts the number of transects sampled at each depth
# output is converted to as.data.frame b/c dplyr automatically coerces it into tbl_df

sampling_effort <- pisco_raw
sampling_effort <- filter(sampling_effort, month %in%pisco_summer)
sampling_effort$fish_tl <- floor(sampling_effort$fish_tl)
sampling_effort$depth <- floor(sampling_effort$depth)
rm(list = "pisco_summer")

sampling_effort <- as.data.frame(
  sampling_effort %>%
    group_by(year, site, side, zone, depth) %>%
    summarize(count = n_distinct(transect)) %>%
    group_by(depth) %>%
    summarize(count = sum(count)) %>%
    dplyr::rename(samples = count)
)
# rename() is in both 'dpylr' and 'reshape'. 'dplyr::' tells r which package to use for that function from

sampling_effort <- na.omit(sampling_effort)

#bar graph of PISCO effort
ggplot(sampling_effort, aes(x = depth, weight = samples)) +
  geom_bar() +
  ggtitle("PISCO sampling effort") +
  xlab("Depth") +
  ylab("Number of samples")

# counts controlled for effort ----------------------------------------------

# each 'count' is needs to be divided by the 'samples' row of sampling_effort
# you have to match up the depth values of bb5 and sampling_effort to makes sure that the correct value for 'samples is used
# I'm guessing that it has something to do with sapply

control_for_effort <- function(counts, samples){
  # counts is the dataframe containing frequency and the associated depth
  # samples is the dataframe containing the number of samples taken at each depth
  
  #the for loop iterates through the rows of 'couts' and divides the frequency by the number of samples at the depth that matches both counts[i,] and samples$depth
  for(i in 1:nrow(counts)){
    print(samples[samples$depth == counts[i,1] , 2])
    counts[i,4] <- counts[i,4] / samples[samples$depth == counts[i,1] , 2]
  }
  # normalize values to 1
  # counts$count is calling the 'count' column in the 'counts' dataframe  
  counts$count <- counts$count / min(counts$count)
  return(counts)
}
bb6 <- control_for_effort(counts = bb5, samples = sampling_effort)
bb6 <- rename(bb6, size = length)
bb7 <- aggregate(count ~ size, data = bb6, FUN = sum)

#KS Test ---------------------------------------------------

ks.test(x = bb7$count[which(bb7$size>=19 & bb7$size<=27)], y = tidy_data$count[which(tidy_data$size>=19 & tidy_data$size<=27)])

#D = 1
#p-value = 4.114e-05
#alternative hypothesis : two-sided
#p<0.05 fails to provide evidence that the distributions are different. However, it does not provide statistical evidence for the sameness of the distribution.

# conversion estimate -------------------------

conversion_estimate <- function(base_data, changing_data, size_range, estimate, incriment = 0.5){
  #`trial` will test the estimate and find the difference between the two histograms with the current estimate
  trial <- function(){
    #initialize two vectors with zeros
    base <- rep(0, length(size_range))
    change <- rep(0, length(size_range))
    #use "size_range" for the names of the vector
    names(base) <- size_range
    names(change) <- size_range
    #convert data with starting conversion estimate
    converted_data <- changing_data
    converted_data$count <- changing_data$count * estimate
    #match data into the named lengths
    mb <- match(names(base), as.character(base_data$size))
    base[which(!is.na(mb))]    <- base_data$count[mb[!is.na(mb)]]
    mc <- match(names(change), as.character(changing_data$size))
    change[which(!is.na(mc))]    <- converted_data$count[mc[!is.na(mc)]]
    #subtract changed data from base data
    difference <- 0
    difference <- base - change
    return(abs(sum(difference)))
  }
  sample1 <- trial()
  estimate <- estimate - incriment
  sample2 <- trial()
  if(sample1>sample2){
    while(sample1>sample2){
      sample1 <- sample2
      estimate <- estimate - incriment
      sample2 <- trial()
    }
    return(estimate)
  } else if(sample1 < sample2){
    estimate <- estimate + 2*incriment
    sample2 <- trial()
    while(sample1 > sample2){
      sample1 <- sample2
      estimate <- estimate + incriment
      sample2 <- trial()
    }
    return(estimate)
  } else {
    return("oops")
  }
}
converter <- conversion_estimate(base_data = bb7, changing_data = tidy_data, size_range = 19:27, estimate = 23, incriment = 0.0001)

#the conversion estimate for size 19-27 is 22.9537000000001!!!!!!

#including 19cm changes it to 57.8662
#however, when plotting the mortality regression, it doesn't change much. When using the conversion 20-38 and plotting log[>=19] the slope is -0.2852. While using the conversion 19-38 plotting log[>=19] the slope is -0.2846. Using the first conversion, and plotting log[>=20] slope = -0.2925. So including 19cm in the regression does change it, but the estimate itself isn't the important part.

# combining data with estimate --------------------------------------------------------------

CCFRP_converted <- tidy_data_year
CCFRP_converted$count <- tidy_data_year$count * converter
combined_dataset <- rbind(bb6, CCFRP_converted) %>%
  group_by(size, year, depth) %>%
  summarise(count = sum(count))
#after using group_by, there is a weird error where it gives multiple warnings: "Unknown column 'depth'". This is because group_by is making it both a tibble and a dataframe. To get rid of the error, just make `combined_dataset` a dataframe
combined_dataset <- data.frame(combined_dataset)

# combined size distribution ----------------------------------------------
ggplot(combined_dataset, aes(x=size, weight=count)) +
  geom_bar() +
  ggtitle("Combined data size distribution (summer)")

# regression mortality slope ----------------------------------------------------------------------------
# plot the regression of log[20-45]
# the slope is mortality of the population
# mortality is combined natural and fishing mortality

#subset data
log_combined <- combined_dataset
log_combined$count <- log(combined_dataset$count)
plot(count ~ size, data = log_combined[which(log_combined$size>=19),], ylab = "log[count]")
abline(lm(count ~ size, data = log_combined[which(log_combined$size>=19),]))
title('Log[19-45cm] size distribution. Slope = -0.2844')

#same regression, except having two different lines. One for 19-33, the other for 34-45.

#these two outputs tell you the equation for the two linear regressions
#lm(count ~ size, data = log_combined[which(log_combined$size>=19 & log_combined<=33),])
#lm(count ~ size, data = log_combined[which(log_combined$size>=34 & log_combined<=45),])

ggplot(data=log_combined[which(log_combined$size>=19 & log_combined$size<=45),] , aes(x=size, y = count)) +
  geom_point() +
  geom_smooth(data=log_combined[which(log_combined$size>=19 & log_combined$size<=33),], method = "lm", se=FALSE) +
  geom_smooth(data=log_combined[which(log_combined$size>=34 & log_combined$size<=45),], method = "lm", se=FALSE) +
  ggtitle("log[19-33cm] & log[34-45cm] size distribution (summer)") +
  geom_text(label="y= -0.06872x + 8.13579", aes(x=30,y=7.5)) +
  geom_text(label="y= -0.5425x + 23.5469", aes(x=40,y=5))


#plotting PISCO and CCRFP 19-25cm on same plot
#trying to figure out if the mortality of the two differs substantially

#sample code for getting ggplot regression
#ggplot(mydata, aes(x=tb, y=ts, col=pop)) + geom_point() +
#  geom_smooth(method="lm", se=FALSE)

#combine pisco and converted ccfrp into one df with a column that denotes which source of data each comes from
log_pisco <- bb7
log_pisco$source <- "PISCO"
log_ccfrp <- tidy_data
log_ccfrp$count <- tidy_data$count * converter
log_ccfrp$source <- "CCFRP"
log_combined <- rbind(log_ccfrp,log_pisco)
log_combined$count <- log(log_combined$count)
#clean-up so that only `log_combined` remains
rm(list = "log_pisco", "log_ccfrp")

#the size 
ggplot(data=log_combined[which(log_combined$size>=19 & log_combined$size<=27),] , aes(x=size, y = count, col=source)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  ggtitle("log[19-27cm] regression")


# regression count vs size ----------------------------
#run a regression for count explained by size
# check the residuals for depth and year to see how much of the variation is explained by those factors

fit1 <- lm(count ~ size, data = combined_dataset)
fit2 <- lm(count ~ size + depth, data = combined_dataset)
fit3 <- lm(count ~ size+ year, data = combined_dataset)
fit4 <- lm(count ~ size + depth + year, data = combined_dataset)

# setup ---------------------------

rm(list = ls())
library(dplyr)
library(geosphere)
library(ggplot2)
library(scales)
library(plotly) #for making contour plots
library(webshot) #for exporting plotly output

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

# Effort -------------------------------------------------------------------------------
if (exists("effort") == TRUE) {
  print("Effort calculation already exists in the environment")
} else {
  # importing the dataset that contains all drifts CCFRP has conducted
  all_drifts_raw <- read.csv("~/Desktop/Thesis/Raw Data/CCFRP/Hernandez_Drifts.csv")
  
  # loop for calculating the effort
  # this is a function so that the variables remain local instead of being stored globally
  effort_function <- function(){
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

tidy <- function(catch, drifts){
  # intialize matricies and variables in this area
  
  #!!!!!!!!
  # These are placeholders to check the code.
  # REMOVE BEFORE RUNNING THE FUNCTION
  #catch <- BothSpecies
  #drifts <- all_drifts_raw
  #!!!!!!!
  
  # initialize variables
  # max and min depth
  min_depth <- min(c(catch$start_depth, catch$end_depth))
  max_depth <- max(c(catch$start_depth, catch$end_depth))
  # max and min size
  max_size <- max(catch$size)
  min_size <- min(catch$size)
  # three letter site names
  sites <- unique(x = substr(drifts$IDCell.per.Trip, 0, 3))
  
  # initialize marticies
  # total CPUE
  total_CPUE <- matrix(0L, nrow = (max_depth - min_depth + 1), ncol = (max_size - min_size + 1))
  colnames(total_CPUE) <- c(dQuote(min_size:max_size))
  rownames(total_CPUE) <- c(dQuote(min_depth:max_depth))
  # site-specific CPUE
  CPUEperSite <- matrix(0L, nrow = (length(sites)), ncol = (max_size - min_size + 2))
  colnames(CPUEperSite) <- c("Site", dQuote(min_size:max_size))
  CPUEperSite[,1] <- sites
  # construct total effort for each depth
  # effort <- effort_function()
  
  # tidy dataframe (site, size, depth, CPUE)
  tidy_data <- data.frame(row.names = c("site", "size", "depth", "catch"))
  
  
  #loops
  # calculate CPUE per site
  for(x in 1:nrow(catch)){ #this will run through every row in the BLU_only_depths matrix
    # x counts the row of the BLU_only_depths matrix that is being currently worked on
    # determine if the start or end depth is deeper and save the deeper value
    # pmax stands for parallel maximum, which compares two columns and finds the maximum value
    deeper <- pmax(catch$start_depth[x], catch$end_depth[x])
    # determine if the start or end depth is shallower and save the shallower depth
    shallower <-pmin(catch$start_depth[x], catch$end_depth[x])
    # new dataframe to collect all the data in this iteration of the loop
    this_fish <- data.frame(site = c(rep(substr(catch$IDCell.per.Trip[x], 0, 3), times = (deeper - shallower +1))), size = c(rep(catch$size[x], times = (deeper - shallower +1))), depth = c(shallower:deeper), catch = c(rep(1, times = (deeper - shallower +1))))
    # Add 1 to the column associated with the fish's length and all the rows in the range of depths fished for that drift
    tidy_data <- rbind(tidy_data, this_fish)
    # time keeper
    print(x)
  }
  
  # tallies the number of columns that match the 'group_by' call and collapses them down to one row
  tidy_data <- tidy_data %>%
    group_by(site, size, depth) %>%
    tally()
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

tidy_data <- tidy(catch = BothSpecies, drifts = all_drifts_raw)

# CCFRP final matrix ------------------------------------------------------------------------------
# matrix of mixed BLU and DEA with uncontrolled CPUE
# summing CPUE across all depths
tidy_data <- aggregate(CPUE ~ size, data = tidy_data, FUN = sum)
tidy_data <- rename(tidy_data, count = CPUE)

# PISCO section -----------------------
# all the PISCO data was lifted from PISCO_abundance

#importing PISCO data
pisco_raw <- read.csv("~/Desktop/Thesis/Raw Data/PISCO/UCSB_FISH.csv")
# selecting only the blue rockfish entries
myst_only <- filter(pisco_raw, classcode == 'SMYS')
# filtering out rows with NA
myst_only <- filter(myst_only, count != 'NA', fish_tl != 'NA', depth != 'NA')
# rounding down the tail length
myst_only$fish_tl <- floor(myst_only$fish_tl)
# rounding down the depth
myst_only$depth <- floor(myst_only$depth)
#group by year, site, side, zone, depth

# subsetting code with the help of Duncan Lang
d = pisco_raw

smys = subset(d, classcode == 'SMYS')
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
bb5 <- aggregate(count ~ depth + length, data = bb4, FUN = sum)

# sampling effort -------------------------------------------------

#counts the number of transects sampled at each depth
# output is converted to as.data.frame b/c dplyr automatically coerces it into tbl_df

sampling_effort <- pisco_raw
sampling_effort$fish_tl <- floor(sampling_effort$fish_tl)
sampling_effort$depth <- floor(sampling_effort$depth)

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
    counts[i,3] <- counts[i,3] / samples[samples$depth == counts[i,1] , 2]
  }
  # normalize values to 1
  # counts$count is calling the 'count' column in the 'counts' dataframe  
  counts$count <- counts$count / min(counts$count)
  return(counts)
}
bb6 <- control_for_effort(counts = bb5, samples = sampling_effort)
bb7 <- aggregate(count ~ length, data = bb6, FUN = sum)
bb7 <- rename(bb7, size = length)

# Ratio of PISCO to CCFRP in overlap ------------------------------

#trying to get a sense of what the conversion ratio should be
# PISCO counts of size 19-25cm divided by CCFRP 19-25cm counts
sum(bb7$count[bb7$size >= 19 & bb7$size <= 25]) / sum(tidy_data$count[tidy_data$size >= 19 & tidy_data$size<= 25])
# ratio is 67.91852

#PISCO 20-25/CCFRP 20-25
sum(bb7$count[bb7$size >= 20 & bb7$size <= 25]) / sum(tidy_data$count[tidy_data$size >= 20 & tidy_data$size<= 25])

#PISCO 21-25/CCFRP 21-25
sum(bb7$count[bb7$size >= 21 & bb7$size <= 25]) / sum(tidy_data$count[tidy_data$size >= 21 & tidy_data$size<= 25])
# ratio is 50.525

#PISCO 20-26/CCFRP 20-26
sum(bb7$count[bb7$size >= 20 & bb7$size <= 26]) / sum(tidy_data$count[tidy_data$size >= 20 & tidy_data$size<= 26])
# ratio is 47.28337

#PISCO 21-26/CCFRP 21-26
sum(bb7$count[bb7$size >= 21 & bb7$size <= 26]) / sum(tidy_data$count[tidy_data$size >= 21 & tidy_data$size<= 26])
# ratio is 49.69914

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
converter <- conversion_estimate(base_data = bb7, changing_data = tidy_data, size_range = 19:25, estimate = 50, incriment = 0.0001)

#the conversion estimate for size 20-38 is 47.1007!!!!!!
#including 19cm changes it to 57.8662
#however, when plotting the mortality regression, it doesn't change much. When using the conversion 20-38 and plotting log[>=19] the slope is -0.2852. While using the conversion 19-38 plotting log[>=19] the slope is -0.2846. Using the first conversion, and plotting log[>=20] slope = -0.2925. So including 19cm in the regression does change it, but the estimate itself isn't the important part.

# combining data with estimate --------------------------------------------------------------

CCFRP_converted <- tidy_data
CCFRP_converted$count <- tidy_data$count * converter
combined_dataset <- rbind(bb7, CCFRP_converted) %>%
  group_by(size) %>%
  summarise(count = sum(count))

# combined size distribution ----------------------------------------------
ggplot(combined_dataset, aes(x=size, weight=count)) +
  geom_bar() +
  ggtitle("Combined data size distribution")

#data for stacked histogram
stack_PISCO <- bb7
stack_PISCO$source <- "PISCO"
stack_CCFRP <- CCFRP_converted
stack_CCFRP$source <- "CCFRP"
stacked_dist <- rbind(stack_CCFRP,stack_PISCO)

#stacked size distribution histogram
ggplot(stacked_dist, aes(x=size, y=count, fill = source)) +
  geom_bar(stat = "identity") +
  ggtitle("BRF size distribution") +
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "black"),
        panel.grid.minor = element_line(colour = "black"))
rm(list = "stack_PISCO","stack_CCFRP","stacked_dist")

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
  ggtitle("log[19-33cm] & log[34-45cm] size distribution") +
  geom_text(label="y= -0.07917x + 8.71043", aes(x=30,y=7.5)) +
  geom_text(label="y= -0.552x + 24.254", aes(x=40,y=5))


#plotting PISCO and CCRFP 19-25cm on same plot
#trying to figure out if the mortality of the two differs substantially

#sample code for getting ggplot regression
ggplot(mydata, aes(x=tb, y=ts, col=pop)) + geom_point() +
  geom_smooth(method="lm", se=FALSE)

#combine pisco and converted ccfrp into one df with a column that denotes which source of data each comes from
log_pisco <- bb7
log_pisco$source <- "PISCO"
log_ccfrp <- tidy_data
log_ccfrp$count <- tidy_data$count * converter
log_ccfrp$source <- "CCFRP"
log_combined <- rbind(log_ccfrp,log_pisco)
log_combined$count <- log(log_combined$count)
#clean-up so that only `log_combined` remains
rm(list = "log_pisco", "log_ccrfp")

#the size 
ggplot(data=log_combined[which(log_combined$size>=20 & log_combined$size<=27),] , aes(x=size, y = count, col=source)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  ggtitle("log[20-27cm] regression")

# heatmap combined data ---------------------------------------------------
#putting the combined data into a matrix that can be used to display a heatmap

tidy_heatmap <- tidy(catch = BothSpecies, drifts = all_drifts_raw)
tidy_heatmap$CPUE <- tidy_heatmap$CPUE * converter
tidy_heatmap <- rename(tidy_heatmap, count = CPUE) %>%
    select(size, depth, count)
bb6_heatmap <- rename(bb6, size = length)
#combined the two datasets
combined_heatmap_data <- rbind(tidy_heatmap, bb6_heatmap)

#initializing matrix to be populated
combined_heatmap_matrix <- matrix(0L,
       ncol = max(combined_heatmap_data$size) - min(combined_heatmap_data$size) + 1,
       nrow = max(combined_heatmap_data$depth) - min(combined_heatmap_data$depth) + 1)
#name the columns based on size
colnames(combined_heatmap_matrix) <- c(dQuote(min(combined_heatmap_data$size):max(combined_heatmap_data$size)))
#name the rows based on depth
rownames(combined_heatmap_matrix) <- c(dQuote(min(combined_heatmap_data$depth):max(combined_heatmap_data$depth)))

for(i in 1:nrow(combined_heatmap_data)){
  #add count to the count in the named column and row
  combined_heatmap_matrix[dQuote(combined_heatmap_data$depth[i]), dQuote(combined_heatmap_data$size[i])] <- combined_heatmap_matrix[dQuote(combined_heatmap_data$depth[i]), dQuote(combined_heatmap_data$size[i])] + combined_heatmap_data$count[i]
  #start with how to call the named column or row
  #then figure out how to add onto the value already there.
}

#clean up the environment
rm(list = "tidy_heatmap", "bb6_heatmap", "combined_heatmap_data")

heatmap(combined_heatmap_matrix, Rowv = NA, Colv = NA, margins = c(3,3), labRow = 1:42, labCol = 2:49, main = "Combined Heatmap", ylab = "Depth (m)", xlab = "Size (cm)", cexRow = 1.2, cexCol = 1.2)

# ggplot heatmap ------------------------------------------

tidy_heatmap <- tidy(catch = BothSpecies, drifts = all_drifts_raw)
tidy_heatmap$CPUE <- tidy_heatmap$CPUE * converter
tidy_heatmap <- rename(tidy_heatmap, count = CPUE)
tidy_heatmap <- select(tidy_heatmap, size, depth, count)

bb6_heatmap <- rename(bb6, size = length)

combined_heatmap_data <- rbind(tidy_heatmap, bb6_heatmap)
combined_heatmap_data <- combined_heatmap_data %>%
  group_by(size, depth) %>%
  summarize(count = sum(count))
rm(list = "tidy_heatmap", "bb6_heatmap")
#makes a range of values from 0 to the top of your scale
main_val = seq(0,600, length=8)
#selecting the color range that you want to use. This must be equal to the length of main_val
mycol = c("blue","cyan","cadetblue1", "aquamarine", "yellow", "orange", "orangered", "red")

#making the heatmap again, but this time with ggplot to control the tile scale better
ggplot(data = combined_heatmap_data[which(combined_heatmap_data$count<=600),], aes(x=size,y=depth,fill=count)) +
  geom_tile()+
  scale_fill_gradientn(colours = mycol, values = scales::rescale(main_val), breaks = c(0,200,400,600,800)) + #begin and end controls the color pallete. Breaks sets the value range
  ggtitle("Heatmap of Combined data (cutoff @ 600)") +
  #sets a bounding box around PISCO and CCFRP data. the `+0.5` or `-0.5` is to set the offset
  geom_rect(xmin=2-0.5, xmax=49+0.5, ymin=1-0.5, ymax=30+0.5, color="black",alpha=0.5, fill=NA, size = 0.3) +
  geom_rect(xmin=7-0.5, xmax=45+0.5, ymin=5-0.5, ymax=42+0.5, color="black",alpha=0.5, fill=NA, size = 0.3)


#ggplot heatmap control across size ----------------------------------------

#normalize each size value so that you can compare two sizes at the same depth
#ratio is the count for the size class that was most abundant divided by the count for a size class.
combined_dataset$ratio <- sapply(1:nrow(combined_dataset), function(x) max(combined_dataset$count)/combined_dataset$count[x])
size_and_ratio <- as.data.frame(cbind(combined_dataset$size, combined_dataset$ratio))
size_and_ratio<- rename(size_and_ratio, size = V1, ratio = V2)
size_and_ratio <- merge(size_and_ratio,combined_heatmap_data, by ="size")
#multiple each count by the ratio
for(i in 1:nrow(size_and_ratio)){
  size_and_ratio$count[i] <- size_and_ratio$count[i] * size_and_ratio$ratio[i]
}

#plotting it with the entire range of values
ggplot(data = size_and_ratio, aes(x=size,y=depth, fill= count)) +
  geom_tile() +
  scale_fill_gradientn(colours = mycol, values = scales::rescale(main_val), breaks = c(0,200,400,600,800)) + #begin and end controls the color pallete. Breaks sets the value range
  ggtitle("Heatmap of normalized sizes (entire range)") +
  #sets a bounding box around PISCO and CCFRP data. the `+0.5` or `-0.5` is to set the offset
  geom_rect(xmin=2-0.5, xmax=49+0.5, ymin=1-0.5, ymax=30+0.5, color="black",alpha=0.5, fill=NA, size = 0.3) +
  geom_rect(xmin=7-0.5, xmax=45+0.5, ymin=5-0.5, ymax=42+0.5, color="black",alpha=0.5, fill=NA, size = 0.3)

#using a cutoff of mean+2sd for better resolution
ggplot(data = size_and_ratio[which(size_and_ratio$count<=828),], aes(x=size,y=depth, fill= count)) +
  geom_tile() +
  scale_fill_gradientn(colours = mycol, values = scales::rescale(main_val), breaks = c(0,200,400,600,800)) + #begin and end controls the color pallete. Breaks sets the value range
  ggtitle("Heatmap of normalized sizes (cutoff at mean+2sd, 828)") +
  #sets a bounding box around PISCO and CCFRP data. the `+0.5` or `-0.5` is to set the offset
  geom_rect(xmin=2-0.5, xmax=49+0.5, ymin=1-0.5, ymax=30+0.5, color="black",alpha=0.5, fill=NA, size = 0.3) +
  geom_rect(xmin=7-0.5, xmax=45+0.5, ymin=5-0.5, ymax=42+0.5, color="black",alpha=0.5, fill=NA, size = 0.3)


#ggplot heatmap controlled across size w/ regression line ----------------------------------------
#the previous section controlled across size by making a ratio of the count at a size to the maximum count (7cm)
#the problem is that it is using noisy data to make a noisy conversion
#this will plot a flat line across 1-19cm, a regression line 20-33cm, and a regression line34-45cm
#the new ratio will be the point on the line at a size divided by maximum count (7cm)
#for the line @ 1-19cm, it will be the value of the regression line @ 20cm

#fit regression line for 34-40cm
#equation is y = -114.3x + 4460.9
lm(combined_dataset$count[33:39] ~ combined_dataset$size[33:39])


#fit regression line for 20-33cm
#equation is y = -105.6x + 4404.5
#lm(combined_dataset$count[19:32] ~ combined_dataset$size[19:32])

#value for the flat line is y = -105.6x + 4404.5
#where x = 20
# value = 2292.5

#plot of the original datasets with the regression lines on top
ggplot(combined_dataset, aes(x =  size, weight = count)) +
  geom_bar() +
  geom_segment(x = 2, y = 2292.5, xend = 19, yend = 2292.5, color = "Red", size = 0.1) +
  geom_segment(x = 20, y = 2292.5, xend = 33, yend = 919.7, color = "Red", size = 0.1) +
  geom_segment(x = 34, y = 574.7, xend = 40, yend = -111.1, color = "Red", size = 0.1) +
  ggtitle("Size distribution with regression lines for ratios")


#normalizing the count data using the regression lines
normalized_by_regression <- combined_dataset
normalized_by_regression <- mutate(normalized_by_regression, ratio = 0)
#this code gets a bit jenk b/c there are no observations for 1cm or 43 cm, so I have to add to i in order to have the right size
for(i in 1:44){
  if(i+1 <= 20){ #i+1 because 1cm is missing. e.g 15cm is in row 14
    normalized_by_regression$ratio[i] <- 5039.7754905/2292.5
  } else if(i+1 > 20 & i+1 < 34){
    normalized_by_regression$ratio[i] <- 5039.7754905/(-105.6*(i+1)+4404.5)
  } else if (i+1 >=34 & i+1 < 43){
    normalized_by_regression$ratio[i] <- 5039.7754905/(-136.8*(i+1)+5263)
  } else if(i==42 || i==43){
    normalized_by_regression$ratio[i] <- 5039.7754905/(-136.8*(i+2)+5263) #i+2 is the size in row 42&43 because 1cm and 43cm are missing
  } else {
    normalized_by_regression$ratio[i] <- 5039.7754905/(-136.8*(49)+5263)
  }
}

rm(list = "i")

#fiddling with ratios that are out of line
#40-49cm spit out negative values
#so I decided to use the absolute value of the ratio
normalized_by_regression$ratio[42:44] <- abs(normalized_by_regression$ratio[42:44])
#ratios for 39-40cm are incosistent with the trend, so we decided to just take the ratio for 38cm for both
normalized_by_regression$ratio[38:39] <- normalized_by_regression$ratio[37]

#eliminate the `count` column so that it can be replaced later
normalized_by_regression$count <- NULL


#associating ratio with each size at each depth
normalized_by_regression <- merge(normalized_by_regression,combined_heatmap_data, by ="size")

#Normalize the ratios so that the ratio for 2-19cm =1 instead of 2.19
#essentially, dividing all ratios by the ratio for 2-19cm
normalized_by_regression$ratio <- normalized_by_regression$ratio/normalized_by_regression$ratio[1]

#multiply count by the ratio
for(i in 1:nrow(normalized_by_regression)){
  normalized_by_regression$count[i] <- normalized_by_regression$count[i] * normalized_by_regression$ratio[i]
}

#plot of the conversion ratios shows that the ratios above 40cm are off the chart
#So the heatmap only includes <=40cm
#using a cutoff of mean+2sd for resolution
#something broke when reproducing the normalized_by_regression plots. The mean+2sd is not 112
#!!!!!recheck mean+2sd!!!!!
ggplot(data = normalized_by_regression[which(normalized_by_regression$count<=112 & normalized_by_regression$size<=40),],
  aes(x=size,y=depth, fill= count)) +
  geom_tile() +
  scale_fill_gradientn(colours = mycol, values = scales::rescale(main_val), breaks = c(0,20,40,60,80, 100, 120)) + #begin and end controls the color pallete. Breaks sets the value range
  ggtitle("Heatmap of regression-normalized sizes below 40cm (cutoff at mean+2sd, 112)") +
  #sets a bounding box around PISCO and CCFRP data. the `+0.5` or `-0.5` is to set the offset
  geom_rect(xmin=2-0.5, xmax=49+0.5, ymin=1-0.5, ymax=30+0.5, color="black",alpha=0.5, fill=NA, size = 0.3) +
  geom_rect(xmin=7-0.5, xmax=45+0.5, ymin=5-0.5, ymax=42+0.5, color="black",alpha=0.5, fill=NA, size = 0.3)


#plotting only depths <=30m
ggplot(data = normalized_by_regression[which(normalized_by_regression$depth>=30 & normalized_by_regression$count<=521),], aes(x=size,y=depth, fill= count)) +
  geom_tile() +
  scale_fill_gradientn(colours = mycol, values = scales::rescale(main_val), breaks = c(0,100,200,300)) + #begin and end controls the color pallete. Breaks sets the value range
  ggtitle("Regression-normalized sizes in deeper water(cutoff at mean+2sd, 521)")



#rescaled the data so that the plot scale is 0-1
normalized_rescaled <- normalized_by_regression [which(normalized_by_regression$size<=40),]
#mean+2sd is 229.1472
#!!!!!!! recheck mean+2sd calculation b/c you plotted a new right-hand regression!!!!!!
normalized_rescaled <- normalized_rescaled [which(normalized_rescaled$count<=230),]
normalized_rescaled$count <- normalized_rescaled$count/max(normalized_rescaled$count)
#change the name of $count to $Relative abundance
normalized_rescaled <- rename(normalized_rescaled, Relative_abundance = count)

#plotting with the normalized, rescaled data
ggplot(data = normalized_rescaled, aes(x=size,y=depth, fill= Relative_abundance)) +
  geom_tile() +
  scale_fill_gradientn(colours = mycol, values = scales::rescale(main_val), breaks = c(0,0.2,0.4,0.6,0.8,1)) + #begin and end controls the color pallete. Breaks sets the value range
  ggtitle("Relative abundance of BRF with regression-normalized sizes (mean+2sd)")


# contour plot ---------------------------------------------------------
#trying to get a contour plot of the heatmap

#using plotly package to get a contour plot
#!!!!!!RECHECK THE MEAN+2SD CALCULATION!!!!!!!!!!!!!!
contour_plot <- plot_ly(data = normalized_by_regression[which(normalized_by_regression$count<=112 & normalized_by_regression$size<=40),],
        type = 'contour',
        x= ~size, y= ~depth, z= ~count, #sets
        colorscale = 'Jet', #changes the color scale
        autocontour = F, #turns of automatice contoursing so that you can set parameters in next line
        contours = list(start = 0, end = 112, size = 7) #sets the range and size of the contour values
        )
#output of the contour plot with a title
contour_heatmap <- layout(contour_plot, title = "Normalized combined data <= 40cm (mean+2sd)")
#install_phantomjs() command line tool necessary for exporting
export(contour_heatmap, file = "Contour heatmap3.png")

# ratios and size distribution --------------------

#the ratios used to to normalized the size distributions for the ggplot heatmaps and contour map
ratios <- normalized_by_regression %>%
  group_by(size, ratio) %>%
  filter(size <=40) %>%
  summarize(count = sum(count))

#bar graph of the ratios
ggplot(ratios, aes(x = size, weight = ratio)) +
  geom_bar() +
  ggtitle("Size distribution conversion ratios") +
  ylab("Conversion ratio") +
  xlab("Size")

#bar graph of the size distribution after normalizing with the ratios
ggplot(ratios, aes(x = size, weight = count/2292.5)) +
  geom_bar() +
  ggtitle("Regression normalized size distribution") +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red")

#ggplot heatmap normalized to 1 ----------------------------------------------
#Loo wants a heatmap with the midpoint = 1.0
#He wants normalized data to show what is more or less abundant
#In theory, this works by taking the mean when all the depths are combined
#The problem is that everything would be below that line on the heatmap.
#So I am going to try something different
#Take the data and only keep mean+2sd
#Then find the mean of the mean+2sd data
#That new mean=1.0 on the scale
#plot the data as a heatmap

#new vector that only has mean+2sd from `normalized_by_regression`
heatmap_scale1 <- normalized_by_regression[which(normalized_by_regression$count <= 338.75 & normalized_by_regression$size<40),]
#normalize count with x-min(x)/max(x)-min(x)
  #this normalization creates a value range 0-1
  #multiply that formula by 2 so that the midpoint is 1 instead of 0.5
heatmap_scale1$count <- 2*(heatmap_scale1$count-min(heatmap_scale1$count))/(max(heatmap_scale1$count)-min(heatmap_scale1$count))

#plot the normalized data
ggplot(data = heatmap_scale1, aes(x=size, y=depth, fill=count)) +
  geom_tile() +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  theme_dark() +
  ggtitle("Something interesing this way comes")

# kelp vs everything and depth -------------------------------------------------------------------
#Loo wants two plots that show the frequency of observations in at each depth in- and outside of kep forest
#collapse counts so that they are summed at each depth

#PISCO
kelp_depths <- aggregate(count ~ depth, data = bb6, FUN = sum)

ggplot(data = kelp_depths, aes(x=depth, weight=count)) +
  geom_bar() +
  ggtitle("Kelp forest observations")

#CCFRP
everywhere_depths <- aggregate(CPUE ~ depth, data = tidy(catch = BothSpecies, drifts = all_drifts_raw), FUN = sum)

ggplot(data = everywhere_depths, aes(x=depth, weight=CPUE*47.1007)) +
  geom_bar() +
  ggtitle("Outside kelp forest observations")

#combined plot of inside and outside
ggplot() +
  geom_bar(data = kelp_depths, aes(x=depth, weight=count, fill="inside")) +
  geom_bar(data = everywhere_depths, aes(x=depth, weight=CPUE*47.1007, fill="outside")) +
  ggtitle("Index of abundance by depth") +
  scale_fill_manual(values = c(inside="blue", outside="red"), name = "location", labels =c("Inside kelp forest", "Outside kelp forest"))

# von Bertalanffy Growth -----------------------------------------------------------

library(FSA)
library(nlstools)

vonbertgrowth <- function(ta,t0,Linf,k){
  Lt = data.frame(size = NA, age = NA)
  for (i in 1:ta){
    Lt[i,1] = Linf * (1-exp(-k*(i-t0)))
  }
  Lt$age <- 1:ta
  return(Lt)
}

blueRockfish_vbg_f <- vonbertgrowth(ta=50,t0=-1.145,Linf=38.150142,k=0.172)
#citations for these parameters are in the email that Lauren sent
#ta?

ggplot(data = blueRockfish_vbg_f, aes(y = size, x = age)) +
  geom_smooth()

smooth_vals = predict(loess(age ~ size, data = blueRockfish_vbg_f), newdata = expand.grid(size=12:38))
diff(smooth_vals)
# gives the point value for the smoothed line at each integer value of x (ie size at each year of age)

#cohort sizes ---------------------------------------------------------
#this is pulled from Lauren's code
#the object of this is to get a graph that shows what size range is for each age

#I put it in a function to scope all the variable to just this instance and not have to worry about extra variables in the global environment
mean_size_at_age <- function(){
##blue rockfish
##derivative of how population changes in time based on von bertalanffy growth.  i.e., this reflects the size distribution at equilibrium under constant recruitment. 
eqn=function(l,y,parms) {
  N1=y[1]
  k=parms[1] 
  M=parms[2] 
  Linf=parms[3] 
  
  dN=numeric(1); # Vector to hold the integrands
  dN[1]=-N1*((M-k)/(k*(Linf-l)))
  return(list(dN))
}

##If need initial abundance N1 to match the mockdata
# F0.05 <- read.csv("/R_stuff/R_Research/MPAs/gopher_rockfish_sizes_F0.05.csv")
# N1 <- sum(F0.05$Year1)
# y=c(N1=20); 
##

##blue rockfish parms
y=c(N1=10); 
parms=c(k=0.172, M=0.14, Linf=38.15);
l=seq(1,parms[3]-0.01,length=10000); #sequence of lengths, previously 10000; another option: length=round(parms[3]

#Integrate the derivative to get abundance N as a function of length l
out=ode(y,l,eqn,parms) #out[,2] are the abundances at sizes 1-34
t0=-0.5
max_age=35 #was 35; try to lower this to 10, 15, 20 to reduce # of cohorts
bin_w <- 1

#Plot the integral of abundance N as a function of length l
plot(out, col=1, xlab = "Length", ylab = "Abundance (N)", main = "Gopher rockfish", type = "p")

colnames(out) <- c("length", "N")
abundance_at_size <- data.frame(out)

#Convert ages to mean size-at-age based on von bertalanffy eq.
#t0 = -0.5, so x0 = 3.6375
Linf <- parms[3]
k <- parms[1]

size_at_age <- vector()
for (j in 0:max_age){
  size_at_age[j+1] <- Linf*(1-exp(-k*(j-t0)))
}

mean_size <- vector()
lower_size <- vector()
upper_size <- vector()
sd <- vector()

for (s in 1:length(size_at_age)){
  mean_size[s] <- size_at_age[s]
  sd[s] <- (mean_size[s])*.1
  lower_size[s] <- mean_size[s]-(3*sd[s])
  upper_size[s] <- mean_size[s]+(3*sd[s])
  #upper_size[s] <- ifelse(upper_size[s] > Linf, 34, upper_size[s])
}

blue_sizes <- data.frame(mean_size = size_at_age, sd = sd, lower_size = lower_size, upper_size = upper_size)

blue_sizes <- blue_sizes %>% mutate_all(.,funs(round(.,digits = 2)))

#return(blue_sizes)
}
blue_sizes <- mean_size_at_age()

#graphing the size range
#ggplot with geom_pointrange(x = age, y = mean_size, ymin = blue_sizes$lower_size[i], yman = blue_sizes$upper_size[i])
ggplot(data=blue_sizes, aes(x=1:36, y=mean_size, ymin=lower_size, ymax = upper_size)) +
  geom_pointrange() +
  xlab("Age") +
  ylab("Mean size") +
  ggtitle("Mean size at age w/ 3 sd")

#Solving for F -----------------------------------

#regression slopes are 19-33cm = -0.07917 and 34-45cm = -0.552
#solve for F at each size using slope = -(k + M +F)/((Linf - l)*k)
#where parameters are based on Lauren's model
#k = 0.172
#M = 0.14
#Linf = 38.15
#l = a size of fish

#initialize dataframe
mortality <- data.frame(matrix(NA, ncol=2, nrow=length((19:45)+1)))
colnames(mortality) <- c("size", "F")
mortality[1:27,1] <- c(19:45)
#calculated F for each size class using the equation where you solve for F
for(i in 19:45){
  if(i < 34){
  mortality$F[which(mortality$size==i)] <- 0.01361724*(38.15-i)-0.312
  } else {
  mortality$F[which(mortality$size==i)] <- 0.094944*(38.15-i)-0.312
  }
}


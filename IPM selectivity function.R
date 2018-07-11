#determining if the IPM model was correct

rm(list = ls())
library(dplyr)
library(ggplot2)
library(geosphere)

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

# PISCO/(PISCO+CCFRP) by size -----------------------------------------------------------

#the combined data has some sizes that aren't in PISCO, so start with combined_dataset to get all sizes
IPM <- combined_dataset
#join the combined data with PISCO data by matching 'size'row. This allows join even though they are not the same length
IPM <- left_join(IPM,bb7, by='size')
#remove the 'count' column from the combined data
IPM$count.x <- NULL
#replace NA with 0
IPM[is.na(IPM)] <- 0
#left_join renamed the columns b/c there were two 'count' columns. So this is renaming it back to 'count'
IPM <- rename(IPM, count = count.y)
#PISCO/(PISCO+CCFRP) by size
IPM$count <-IPM$count/combined_dataset$count

ggplot(IPM, aes(x=size, y =count)) +
  geom_bar(stat= 'identity')

#count density 1-45cm. Had to manually input 0s for the missing size classses
IPM_density <- c(0, IPM$count[1:41], 0,1,0,0,0,0,1)
#mean(IPM_density) = 0.5454389
#sd(IPM_density) = 0.3727877
#an alternative when you eliminate the outliers are the right tail
IPM_density <- c(0, IPM$count[1:37], 0, IPM$count[39:41],0,0,0,0,0,0,0)
#mean(IPM_density) = 0.4867295
#sd(IPM_density) = 0.3783512
IPM_curve <- qnorm(IPM_density, 0.4867295, 0.3783512)
plot(1:49, IPM_curve, xlab = "Size (cm)", ylab = "Probability of being observed", main = "CDF for IPM estimation")

qplot(1:49, IPM_curve, geom = "smooth"
      )

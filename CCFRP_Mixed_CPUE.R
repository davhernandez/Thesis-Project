# Setup ----------------------------

rm(list = ls())
library(dplyr)
library(geosphere)
library(ggplot2)
library(gridExtra)

# raw data -----------------------------------------------------------

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

# Heatmap Matrix ---------------------------------------------------------
  
  # putting data info into a matrix for a heatmap
heatmap_matrix <- function(data_input, CPUE_column){
  CPUE_matrix <- matrix(0L, nrow = (max(data_input$depth) - min(data_input$depth) + 1), ncol = (max(data_input$size) - min(data_input$size) + 1))
  # naming columns and rows
  rownames(CPUE_matrix) <- c(dQuote(min(data_input$depth):max(data_input$depth)))
  colnames(CPUE_matrix) <- c(dQuote(min(data_input$size):max(data_input$size)))
  
  for(i in 1:nrow(data_input)){
    # add the row's value of CPUE_matrix to the current value in the matrix
    CPUE_matrix[dQuote(data_input$depth[i]), dQuote(data_input$size[i])] <- data_input[i,eval(CPUE_column)] +       CPUE_matrix[dQuote(data_input$depth[i]), dQuote(data_input$size[i])]
    }
  return(CPUE_matrix)
  }

# heatmap Mixed Raw CPUE ------------------------------------------------------------------------------
# heatmap of mixed BLU and DEA with uncontrolled CPUE
uncontrolled_CPUE <- heatmap_matrix(data_input = tidy_data, CPUE_column = "CPUE")
heatmap(uncontrolled_CPUE, Rowv = NA, Colv = NA, margins = c(3,3), labRow = 5:54, labCol = 7:45, main = "CCFRP CPUE (BLU & DEA)", ylab = "Depth (m)", xlab = "Size (cm)", cexRow = 1.2, cexCol = 1.2)

# heatmap controlled CPUE -----------------------------------------------------------------------------

# heatmap for the CPUE that is controlled for amount of effort at each site
controlled_CPUE <- heatmap_matrix(data_input = tidy_data, CPUE_column = "normalized_CPUE")
heatmap(controlled_CPUE, Rowv = NA, Colv = NA, margins = c(3,3), labRow = 5:54, labCol = 7:45, main = "CCFRP CPUE Controlled for Site CPUE", ylab = "Depth (m)", xlab = "Size (cm)", cexRow = 1.2, cexCol = 1.2)


# Plot of normailized CPUE summed across depths and another summed across size -------------------------------
  # CPUE summed for each size
  ggplot(tidy_data, aes(x = size, y = normalized_CPUE)) +
    geom_bar(stat = "identity") +
    ylab("CPUE")
  #CPUE summed for each depth
  ggplot(tidy_data, aes(x = depth, y = normalized_CPUE)) +
    geom_bar(stat = "identity") +
    ylab("CPUE")
  
# Plot N for each size ------------------------------------------------------
  # simple distribution of catch
  # counting the total number of fish caught for each size
  total_n <- tidy_data %>%
    group_by(size) %>%
    summarize(n = sum(n))
  # plot the number of fish caught at each size
  ggplot(total_n, aes(x = size, y = n)) +
    geom_bar(stat = "identity")
  # now check these distributions for each site
  total_n <- tidy_data %>%
  group_by(site, size) %>%
    summarize(count = sum(n))
  # this plot doesnot work yet. It may be because tidy_data is making sites with 'NA'. It may be something else
  ggplot(total_n, aes(x = size, y = count)) +
    geom_bar(stat = "identity") +
    facet_grid(~site)

# BRF only heatmap CPUE ------------------------------------------------------
  # They started distinguishing between BRF and Deacon in 2011
    # select only BRF caught in 2011 or later
  BRF_only <- filter(BothSpecies, species_code == "BLU", Year >= 2011)
  # run through tidy function
  tidy_BRF <- tidy(catch = BRF_only, drifts = all_drifts_raw)
  #decide whether or not to control for effort
    #decided against it
  # adding data to a matrix for the heatmap
  BRF_heatmap <- heatmap_matrix(data_input = tidy_BRF, CPUE_column = "CPUE")
  # making the heatmap
  heatmap(BRF_heatmap, Rowv = NA, Colv = NA, margins = c(3,3), labRow = 5:42, labCol = 7:41, main = "CCFRP CPUE for Blue RF 2011-2016", ylab = "Depth (m)", xlab = "Size (cm)", cexRow = 1.2, cexCol = 1.2)
  #BRF_only size distribution plot
  
# Deacon only heatmap CPUE ---------------------------------------------------
  # They started distinguishing between BRF and Deacon in 2011
  # select only BRF caught in 2011 or later
  DEA_only <- filter(BothSpecies, species_code == "DEA", Year >= 2011)
  # run through tidy function
  tidy_DEA <- tidy(catch = DEA_only, drifts = all_drifts_raw)
  #decide whether or not to control for effort
  #decided against it
  # adding data to a matrix for the heatmap
  DEA_heatmap <- heatmap_matrix(data_input = tidy_DEA, CPUE_column = "CPUE")
  # making the heatmap
  heatmap(DEA_heatmap, Rowv = NA, Colv = NA, margins = c(3,3), labRow = 5:42, labCol = 11:38, main = "CCFRP CPUE for Deacon RF 2011-2016", ylab = "Depth (m)", xlab = "Size (cm)", cexRow = 1.2, cexCol = 1.2)
  
  
# normalized CPUE ------------------------------------------------------------------------------
  
  # normalize each cloumn (size) so that the heatmap shows a better depiction of depth distribution within each size
size_normalization <- function(input_CPUE){
  #initialize the output
  output_matrix <- matrix(data = 0L, nrow = nrow(input_CPUE), ncol = ncol(input_CPUE))
  
  for(i in 1:ncol(input_CPUE)){
    #calculate mean and sd for the current column
    k <- mean(input_CPUE[,i])
    l <- sd(input_CPUE[,i])
    
    for(j in 1:nrow(input_CPUE)){
      #normalize the current cell
      output_matrix[j,i] <- (input_CPUE[j,i] - k)/l
      }
    }
  return(output_matrix)
  }
  
# normalizing raw CPUE --------------------------------------------------------
  
  nrm_uncontrolled_CPUE <- size_normalization(input_CPUE = uncontrolled_CPUE)
  heatmap(nrm_uncontrolled_CPUE, Rowv = NA, Colv = NA, margins = c(3,3), labRow = 5:54, labCol = 7:45, main = "CCFRP CPUE Normalized by size", ylab = "Depth (m)", xlab = "Size (cm)", cexRow = 1.2, cexCol = 1.2)
  
  # look at BothSpecies for the catch
  # figure out if you divide catch by effort per site or total effort for the depths fished
    # so are you looking at catch per effort at a location or catch per effort in its entirety
  # find the row sums from the site-specific CPUE and save it
  
  # one is catch per unit effort at a location
  #the other is catch per total effort, controlling for the amount of catch that each site contributes to the total CPUE
  # you want the second one, which is catch/total effort for the depths, then divided by site CPUE

# CPUE distribution for size --------------------------------------------
p1<-ggplot(aggregate(CPUE ~ size, data = tidy_data, FUN = sum),
    aes(x = size, y = CPUE)) +
    geom_bar(stat = "identity") +
    ggtitle("Mixed CPUE by size")
p2<-ggplot(aggregate(CPUE ~ size, data = tidy_data, FUN = sum),
    aes(x = size, y = log(CPUE))) +
    geom_bar(stat = "identity") +
    ggtitle("log Mixed CPUE by size")
p3<-ggplot(aggregate(CPUE ~ size, data = tidy_BRF, FUN = sum),
    aes(x = size, y = CPUE)) +
    geom_bar(stat = "identity") +
    ggtitle("BRF CPUE by size")
p4<-ggplot(aggregate(CPUE ~ size, data = tidy_BRF, FUN = sum),
    aes(x = size, y = log(CPUE))) +
    geom_bar(stat = "identity") +
    ggtitle("log BRF CPUE by size")
p5<-ggplot(aggregate(CPUE ~ size, data = tidy_DEA, FUN = sum),
    aes(x = size, y = CPUE)) +
    geom_bar(stat = "identity") +
    ggtitle("DEA CPUE by size")
p6<-ggplot(aggregate(CPUE ~ size, data = tidy_DEA, FUN = sum),
    aes(x = size, y = log(CPUE))) +
    geom_bar(stat = "identity") +
    ggtitle("log DEA CPUE by size")
grid.arrange(grobs = list(p1, p2, p3, p4, p5, p6))
rm(list = "p1", "p2", "p3", "p4", "p5", "p6")

# size distributions for BRF and DEA separately --------------------------

tidy_BRF <- aggregate(CPUE ~ size, data = tidy_BRF, FUN = sum)
ggplot(data=tidy_BRF, aes(x = size, weight = CPUE)) +
  geom_bar() +
  ylab("CPUE") +
  ggtitle("CCFRP Blues only (2011-2016)")

tidy_DEA <- aggregate(CPUE ~ size, data = tidy_DEA, FUN = sum)
ggplot(data=tidy_DEA, aes(x = size, weight = CPUE)) +
  geom_bar() +
  ylab("CPUE") +
  ggtitle("CCFRP Deacon only (2011-2016)")

#both on the same plot
ggplot() +
  geom_bar(data=tidy_BRF, aes(x = size, weight = CPUE, fill = "BRF")) +
  geom_bar(data=tidy_DEA, aes(x = size, weight = CPUE, fill = "DEA")) +
  ylab("CPUE") +
  ggtitle("CCFRP Blue vs Deacon size distribution (2011-2016)") +
  scale_fill_manual(values = c(BRF="blue", DEA="red"), name = "species", labels =c("Blue Rockfish", "Deacon Rockfish"))

# sandbox ---------------------------------------------------------------

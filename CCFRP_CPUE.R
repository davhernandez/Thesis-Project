# setup ------------------------------------------------------------------
library(dplyr)
library(geosphere)
library(ggplot2)
library(ggjoy)
library(matrixStats)
library(cowplot)
library(gridExtra)
#importing the CCFRP dataset
#the next two lines are how to read the excel file
#library(readxl)
#ccfrp_raw <- read_excel("~/Desktop/Thesis/Raw Data/CCFRP/Hernandez_BlueRF.xlsx")
ccfrp_raw <- read.csv("~/Desktop/Thesis/Raw Data/CCFRP/Hernandez_BlueRF.csv")


# Catch data manipulation ------------------------------------------------------

# removing entries for Deacon Rockfish (the subspecies that looks like blue rockfish)
BLU_only <- filter(ccfrp_raw, Species.Code == "BLU")
#renaming the columns
BLU_only <- rename(BLU_only, fish_ID = Fish.ID, species_code = Species.Code, common_name = Common.Name, length = Length..cm., cell_ID = Cell.ID, site_type = Site..MPA..REF., start_depth = Start.Depth..ft., end_depth = End.Depth..ft., start_time = Start.Time, end_time = End.Time, drift_ID = Drift.ID)
# selecting the columns that will be used
BLU_only <- select(BLU_only, IDCell.per.Trip, length, Area, Year, cell_ID:end_time, drift_ID)
#using distCosine from the geosphere package to calculate the distance of each drift based on the start and end lat/long. The lat/longs are concatinated with cbind() instead of just c() so that mutate works on single row at a time. The output is in meters.
BLU_only <- mutate(BLU_only, distance_meters = distCosine(p1 = cbind(ST_LonDD, ST_LatDD), p2 = cbind(End_LonDD, End_LatDD)))
#converting the start_time and end_time columns into POSIXct format, which is necessary for the difftime() function
BLU_only$start_time <- as.POSIXct(BLU_only$start_time, format = "%H:%M:%S")
BLU_only$end_time <- as.POSIXct(BLU_only$end_time, format = "%H:%M:%S")
#creating a new column (total_time) that is the total amount of time (in minutes) spent on the drift
BLU_only <- mutate(BLU_only, total_time = difftime(time1 = BLU_only$end_time, time2 = BLU_only$start_time, units = "mins"))
# eliminating the catch entries that have 'NA' for length
# the comma at the end of the [(...),] is important because you are defining the rows of BLU_only$length to queries fro NAs
BLU_only <- BLU_only[!is.na(BLU_only$length),]
# below is another way of eliminating the NAs from the row 'length'
#BLU_only <- filter(BLU_only, length != 'NA')

# eliminating all entries that don't have a start or end depth
BLU_only_depths <- filter(BLU_only, start_depth != 'NA', end_depth != "NA")
# convert depths from ft to m
BLU_only_depths$start_depth <- BLU_only_depths$start_depth * 0.3048
BLU_only_depths$end_depth <- BLU_only_depths$end_depth * 0.3048
# rounding down the depth columns
# this is important for calculating the correct dpeths for the Catch matrix
# see 'Catch' section for explanation of the rounding
BLU_only_depths$start_depth <- floor(BLU_only_depths$start_depth)
BLU_only_depths$end_depth <- floor(BLU_only_depths$end_depth)


# simple frequency distribution -----------------------------------------------

#plots a simple histogram of the frequency of observations for each length
hist(BLU_only$length)
#the same histogram achieved with ggplot and bining to 1cm
ggplot(BLU_only, aes(length)) + geom_histogram(binwidth = 1)


#checking the depth distrbution -----------------------------------------------

# making a joy plot to see if there is a large spread in the depths that they fished at for both start and end depths
#new blank matrix with a column for the name of the factor (start/end) and a column for the depth in meters
small_only <- data.frame(data=NA, nrow = 2*10582, ncol = 2)
#applying column nmaes to the blank matrix
colnames(small_only) <- c("start_end", "depth")
#making the first half of the entries assigned to the catergorical variable 'start'
small_only[1:10570,1] <- "start"
# inserting the values for the start depths from the BLU_only_depth matrix
small_only[1:10570,2] <- BLU_only_depths$start_depth
# making the second half of the entries assigned to the catergorical variable 'end'
small_only[10571:21140,1] <- "end"
# inserting the values for the end depths from the BLU_only_depths matrix
small_only[10571:21140, 2] <- BLU_only_depths$end_depth

#scales for ggjoy axes
# find the minimum depth and round it to the nearest whole number
mins <- round(min(small_only[,2]), digits = 0)
# find the maximum dpeth and round it to nearest whole number
maxs <- round(max(small_only[,2]), digits = 0)

#converting the matrix to a dataframe since ggplot cannot use matrices
small_only <- as.data.frame(small_only)

#the plot
ggplot(small_only, aes(x = depth, y = start_end, group = start_end)) + # "group=" is necessary for ggplot to distinguish how to group the data since start_end is a categorical variable
  geom_joy() + #make it a joy plot
  scale_x_continuous(breaks = c(mins, 10, 15, 20, 25, 30, 35, maxs)) + #labeling the x-axis values
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) # removes the y-axis title and the black tick mark next to "start" and "end on the y-axis


# Catch matrix is below the effort calculation -----------------------------------------------------------------------

# the catch calculation was moved to be below the effort calculation because effort fished over a wider range of depths. In order for the column by vector division to work, the catch matrix needed to have the same range of depths as the effort vector

# Effort data manipulation ----------------------------------------------------------------------------------

# importing the dataset that contains all drifts CCFRP has conducted
all_drifts_raw <- read.csv("~/Desktop/Thesis/Raw Data/CCFRP/Hernandez_Drifts.csv")
# renmaing columns
effort_data <- dplyr::rename(all_drifts_raw, start_time = Start.Time, end_time = End.Time, drift_time_hrs = Drift.time_hrs, total_andjusted_angler_hrs = total.adjusted.angler.hrs, drift_length = Drift.Length..m., start_depth = Start.Depth..ft., end_depth = End.Depth..ft.)
# eliminating the depths entries that have 'NA'
effort_data <- filter(effort_data, start_depth != 'NA', end_depth != 'NA')
# converting start and end depth to meters
# ththe meters conversion is combined with rounding the resulting number down to the nearest whole meter by using the floor() function
# I am rounding down because I am binning per meter. So 30.2m counts as falling in the 30-31m bin
effort_data$start_depth <- floor(effort_data$start_depth * 0.3048)
effort_data$end_depth <- floor(effort_data$end_depth * 0.3048)


# Effort Calculation ------------------------------------------------------------

# The point of this section is to create a vector that contains the range of depth values and add up how many angler hours were spent fishing at each depth.
# finding the minimum depth that was fished at
min_depth_effort <- min(c(effort_data$start_depth, effort_data$end_depth))
max_depth_effort <- max(c(effort_data$start_depth, effort_data$end_depth))
#making a new empty vector that is the length of the max depth - min depth
effort <- rep(0, (max_depth_effort - min_depth_effort +1))
# naming the vector so that entries can be input by the name/value of depth
names(effort) <- c(min_depth_effort:max_depth_effort)

# check to see if the names were assigned correctly
# str(effort)

# loop for calculating the effort
# this is a function so that the variables remain local instead of being stored globally
effort_function <- function(){
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
  # after the for loop finished, return the final effort vector
  return(effort)
}

# runs the effort_function and stores it in the effort vector
effort <- effort_function()

# Catch Calculation -----------------------------------------------------------------------------

# creating a matrix for the catch part of CPUE
# this matrix will be divided, column by column, by the the Effort matrix
# Use BLU_only_depths matrix to ensure that depths are in m
# start and end depths will be rounded down because 30.2m counts for the full depth range of 30-31m


# finding the range of fish lengths to make an appropriately sized matrix
# this is the smallest size class of fish that was caught
min_length <- min(BLU_only_depths$length, na.rm = TRUE)
# this is the largest size class of fish that was caught
max_length <- max(BLU_only_depths$length, na.rm = TRUE)
# find the minimum depth that a fish was caught at
min_depth_catch <- min(c(BLU_only_depths$start_depth, BLU_only_depths$end_depth))
# finding the maximum depth that a fish was caught at
max_depth_catch <- max(c(BLU_only_depths$start_depth, BLU_only_depths$end_depth))

# a blank matrix with enough columns to be named min_length:max_length and enough rows to be named min_depth:max_depth
# this matrix will contain all the catch data to ranges over all size classes caught and all depths fished
# starting out with all zeros instead of NA so that later calculations won't run into problems with having non-integer values
# 0L means it is an interger, not a numeric. It uses half the memory (4 bytes instead of 8) which will make loops and calculations faster
catch_matrix <- matrix(0L, nrow = (max_depth_effort - min_depth_effort + 1), ncol = (max_length - min_length + 1))
# name the columns with the range of size classes from the smallest fish to the largest fish caught
# dQuote() adds double quotations around the value called
colnames(catch_matrix) <- c(dQuote(min_length:max_length))
# name the rows with the range of values from the shallowest depths to the deepest depth fished for the entire project. Since catch data for BRF only ranges 5-42m, I used the effort data range, which is 5-54m
rownames(catch_matrix) <- c(dQuote(min_depth_effort:max_depth_effort))

# going row by row, add 1 to the column of the same length for the entire range of depths that were fished
# dQuote is being used to add quotations around the matrix value so that the command is calling the name of the row or column instead of the row or column number
# so that "25" is the column named "25" instead of column number 25 (ie [,25])
# this is a proof of concept piece of code that adds 1 to the column which is associated with the first fish length and the depths of the drift that it was caught at
# If it works, there should be 1s in column 25 (the fish is length: 25cm) and rows 12-13 (start depth: 12, end depth: 13)
# catch_matrix[dQuote(BLU_only_depths$start_depth[1]:BLU_only_depths$end_depth[1]), dQuote(BLU_only_depths$length[1])] <- catch_matrix[dQuote(BLU_only_depths$start_depth[1]:BLU_only_depths$end_depth[1]), dQuote(BLU_only_depths$length[1])] + 1

# loop to fill in the catch matrix in the same manner as above
# I am calling this as a function so that the variables x, deeper, and shallower remain contained in the local environment of the function and do not become variables in the global environment
catch_function <- function(){
  for(x in 1:nrow(BLU_only_depths)){ #this will run through every row in the BLU_only_depths matrix
    # x counts the row of the BLU_only_depths matrix that is being currently worked on
    # determine if the start or end depth is deeper and save the deeper value
    # pmax stands for parallel maximum, which compares two columns and finds the maximum value
    deeper <- pmax(BLU_only_depths$start_depth[x], BLU_only_depths$end_depth[x])
    # determine if the start or end depth is shallower and save the shallower depth
    shallower <-pmin(BLU_only_depths$start_depth[x], BLU_only_depths$end_depth[x])
    # Add 1 to the column associated with the fish's length and all the rows in the range of depths fished for that drift
    catch_matrix[dQuote(shallower:deeper), dQuote(BLU_only_depths$length[x])] <- catch_matrix[dQuote(shallower:deeper), dQuote(BLU_only_depths$length[x])] + 1
    # time keeper
    print(x)
  }
  # returns the completed catch matrix
  return(catch_matrix)
}

# running the catch_function and saving the results to the catch_matrix
catch_matrix <- catch_function()


# CPUE Calculation ----------------------------------------------------------------------

# divide each column of the catch_matrix by the effort vector to calculate CPUE
CPUE <- apply(catch_matrix, 2, function(x) (x/effort))

# simple heatmap to look at the distribution
# cexRow and cexCol scale the size of the axis label text
heatmap(CPUE, Rowv = NA, Colv = NA, margins = c(3,3), labRow = 5:54, labCol = 7:45, main = "CCFRP CPUE", ylab = "Depth (m)", xlab = "Size (cm)", cexRow = 1.2, cexCol = 1.2)


# binned plots/finding a length inflection point ----------------------------------------------------
# by binning fish to small or large sizes and varying what size cosititutes small or large, I hope to find that there is a size at which depth preference changes drastically

#decide on the size for the bins. In this case, 'small' is 7-27cm and 'large' is 28-45cm
# average all CPUE across the size range for 'small'
# dQuote adds quotes around the numbers so that you are calling the column name (fish size) istead of a numbered column
small_bin <- rowMeans(x = CPUE[,dQuote(7:27)])
# average all CPUE across the size range for 'large'
large_bin <- rowMeans(x = CPUE[,dQuote(28:45)])
# combine bins into a dataframe
binned_CPUE <- data.frame(Small = small_bin, Large = large_bin, Depth = 5:54)
# plot the two bins
# the separate geom_line() allow you to plot multiple columns on the same plot
ggplot(binned_CPUE,
  aes(x= Depth)) +
    geom_line(aes(y = Small, colour = "Small")) +
    geom_line(aes(y = Large, colour = "Large")) +
      ylab("CPUE")

# looping it!
binned_plots <- function(){
  #i iterates from 1 more than the shortest length(7cm) to 1 less than the longest length (45cm) because the column call needs two dimensions (e.g. 7:8 is two dimensions, but 7:7 isn't). Using "i in 8:44" gives both bounds of 7:8 (for small bin) and 44:45 (for large bin)
  j <- 0 #initializes 'j', which is the position that a plot is put in for plotlist
  plotlist <- list() #initializing blank list
  for(i in 8:43){
    print(i) #tracks the progress of the loop
    j <- j+1 #increases j with each loop
    small_bin <- rowMeans( x = CPUE[, dQuote(7:i)]) #averages CPUE for the small sizes
    large_bin <- rowMeans(x= CPUE[, dQuote(i:45)])  #averages CPUE for the large sizes
    binned_CPUE <- data.frame(Small = small_bin, Large = large_bin, Depth = 5:54) #makes a new dataframe with three columns: Small, Large, and Depth)
   #print("x")
  
    #prints the plot for each run og the loop  
  plotlist[[j]]<- (ggplot(binned_CPUE,
                aes(x= Depth)) +
           geom_line(aes(y = Small)) + #to add a legend, make aes(y = Small, colour = "Small")
           geom_line(aes(y = Large)) + #to add a legend, make aes(y = Large, colour = "Large")
           ylab("CPUE") +
           labs(title = paste("split at", i, "cm"))
  )
  }
  # calls the grid.arrange function from gridExtra package and applies it to plotlist. do.call expects to pass the function over a list of values
  do.call(grid.arrange, plotlist)
}

binned_plots()

# binned for maturity --------------------------------------------------------------

# making a plot of binned CPUE for juvenile, maturing, and mature at each depth.
# I tried dividing the bins by the number of size classes used to control for bin width, but this doesn't seem to make a difference in the overall shape of the lines. It just changes the height.
# Since this plot is looking for trends at different depths, magnitude of CPUE is not important
#Therefore, controlling for bin width is not important

# average CPUE across all juvenile fish
juvenile_bin <- rowMeans(x = CPUE[,dQuote(7:25)])
# average CPUE across all maturing fish
maturing_bin <- rowMeans(x = CPUE[,dQuote(26:31)])
# average CPUE across all mature fish
mature_bin <- rowMeans (x = CPUE [, dQuote((32:45))])
# combine bins into a dataframe
binned_CPUE <- data.frame(Juvenile = juvenile_bin, Maturing = maturing_bin, Mature = mature_bin, Depth = 5:54)
# plot the two bins
# the separate geom_line() allow you to plot multiple columns on the same plot
ggplot(binned_CPUE,
       aes(x= Depth)) +
  geom_line(aes(y = Juvenile,
                colour = "Juvenile (7-25cm)")) +
  geom_line(aes(y = Maturing,
                colour = "Maturing (26-31cm)")) +
  geom_line(aes(y = Mature,
                colour = "Mature (32-45cm)")) +
  labs(x = "Depth", y = "CPUE", title = "Binned CPUE") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 17)) +
  # centers the title
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))
  # changes the size of the axis labels and axis titles

# plotting size vs CPUE at various depths --------------------------------------------------

size_CPUE_by_depth <- function(){
  #i iterates from 1 more than the shortest length(7cm) to 1 less than the longest length (45cm) because the column call needs two dimensions (e.g. 7:8 is two dimensions, but 7:7 isn't). Using "i in 8:44" gives both bounds of 7:8 (for small bin) and 44:45 (for large bin)
  j <- 0 #initializes 'j', which is the position that a plot is put in for plotlist
  plotlist <- list() #initializing blank list
  for(i in 5:10){
    print(i) #tracks the progress of the loop
    j <- j+1 #increases j with each loop
    #prints the plot for each run og the loop  
    plotlist[[j]]<- (ggplot(as.data.frame(CPUE), aes(x= i)) +
                       geom_histogram(binwidth = 1) +
                       ylab("CPUE") +
                       xlab("Size") +
                       labs(title = paste("Depth:", i, "m"))
    )
  }
  # calls the grid.arrange function from gridExtra package and applies it to plotlist. do.call expects to pass the function over a list of values
  do.call(grid.arrange, plotlist)
}

size_CPUE_by_depth()

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

ggplot(data = blueRockfish_vbg_f, aes(y = age, x = size)) +
  geom_smooth()

smooth_vals = predict(loess(age ~ size, data = blueRockfish_vbg_f), newdata = expand.grid(size=12:38))
diff(smooth_vals)
# gives the point value for the smoothed line at each integer value of x (ie size at each year of age)


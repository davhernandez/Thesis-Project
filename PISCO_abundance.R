# setup ------------------------------------------------------------------

library(dplyr)
library(geosphere)
library(ggplot2)
library(rethinking)
library(pscl)
library(vegan)
library(mgcv)
library(reshape)
library(gridExtra)
#importing PISCO data
pisco_raw <- read.csv("~/Desktop/Thesis/Raw Data/PISCO/UCSB_FISH.csv")

# data manipulation ------------------------------------------------------

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
#tallies the number of rows for each unique depth-length combination
ggplot(bb5, aes(x = length, y = count, group = depth, color = depth)) +
    geom_line()

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
# View(sampling_effort)

ggplot(sampling_effort, aes(depth, samples)) +
  geom_bar(stat = "identity") +
  labs(x = "Depth",
       y = "Samples",
       title = "Sampling Effort") +
  theme_minimal() +
  # a default theme that will remove that changes girdlines and background
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 17)) +
  #centers the title
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))
  # changes the size of the axis labels and axis titles
  
# plotting bb5 ---------------------------------------------
# make a plot of depth x length where count freuqncy is also displayed. It will be similar to a heatmap?
#change bb5 so that it is a df with rows of 'depth', columns of 'length', and values of 'count' using the 'cast' function in 'reshape' package
heatmap_transform <- as.matrix(cast(bb5, depth~length, value="count", fill = 0))
# 'fill = 0' will replace 'NA' with 0
#View(heatmap_transform)
heatmap(heatmap_transform, Rowv = NA, Colv = NA, main = "PISCO Abundance", ylab = "Depth (m)", xlab = "Size (cm)", margins = c(3,3), cexRow = 1.5, cexCol = 1.5)

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

heatmap_effort <- as.matrix(cast(bb6, depth~length, value="count", fill = 0))
# 'fill = 0' will replace 'NA' with 0
heatmap(x = heatmap_effort, Rowv = NA, Colv = NA, main = "PISCO Abundance", xlab = "Size (cm)", margins = c(3,3), cexRow = 1.5, cexCol = 1.5)


#WTF why is the heatmap for the data controlled for effort the same as the heatmap that isn't controlled?!?!?!?!


# selecting only coastal sites ------------------------
  # ie removing all sites around the channel islands

# Pseudo code
  # use grep(), select dataframe[grep(pattern = "site1 | site 2 | site3", x = dataframe$site),]
  # placing grep() in the 'x' position of [x,y] ensures that it checks row values
  # 'x = dataframe$site' looks for these row values in the 'site' column
  # invert = TRUE simply inverts the selection so that it excludes everything that matches the pattern. Just like doing !=

# selects all row values in the 'site' column that don't partially match ANACAPA, CAT,SBI, SCI, SMI, or SRI
coastal_myst <- bb4[grep(pattern = "ANACAPA|CAT|SBI|SCI|SMI|SRI|SNI", x = bb4$site, invert = TRUE),]


# poisson model -------------------------------------------------------
#first check for overdispersion
  # is the mean of the length/count data much smaller than the variance
  # conduct this test by doing mean/var
  # create a vector that places each length as a multiple of the count. Then do mean and variance from that
overdispersion_test  <- function(input){
    string<-c()
    for(i in 1:nrow(input)){
      y <- rep(input$fish_tl[i], times = input$count[i])
      string<-c(string, y)
    }
    return(mean(string)/var(string))
}
# running the overdispersion test
overdispersion_test(myst_only)
# overdispersion is not indicated. Mean/var = 0.2515809
# build poisson model based on rethinking package
model1 <- glm(depth ~ length, family = 'poisson', data =coastal_myst)
summary(model1)
model2 <- glm(depth ~ length * count, family = 'poisson', data = coastal_myst)
summary(model2)

# a quick histogram of the size distributions --------------------------------
ggplot(myst_only, aes(x=fish_tl, weight=count)) + geom_histogram(binwidth = 1)

# presence absence binomial model ----------------------------------------------
# GOAL: create a binomial model of the data based on the presence or absence of s.mystinus at a site in a given year.
#creating a new column for presence/absence
#pisco_raw["presence"] <- NA
#if SMYS is i the classcode column, it gets a '1'. If not, it get s a 0
#pisco_raw$presence <- ifelse(pisco_raw$classcode == 'SMYS', 1, 0)
# sum the prsence absence columns grouped into the heirarchical tree of year, site, side, and zone
#presence_absence <- aggregate(cbind(presence) ~ year+site+side+zone+depth, data=pisco_raw, FUN = sum)
# change any value >0 in the presnece column to a '1'
#presence_absence$presence[presence_absence$presence>0] = 1
#range(presence_absence$presence)
#mean(presence_absence$presence)
# glm for presence ~ depth
#presence_myst_glm <- glm(presence~depth, family = "binomial", data = presence_absence)
#summary(presence_myst_glm)
#plot(presence~depth, data = presence_absence)
#presence_myst_glm_pred <- predict(presence_myst_glm, newdata = presence_absence, type = 'response')
#lines(presence_absence$depth, presence_myst_glm_pred, lwd=2, col=2)
#lot(P.magellanicus.Pabs~Didemnum,data=all.new)
#P.magellanicus.Pabs.pred=predict(P.magellanicus.Pabs.glm,newdata=new.data, type='response')
#lines(new.data$Didemnum,P.magellanicus.Pabs.pred, lwd=2,col=2)

# catch by year --------------------------------------------------

# the purpose of this section is to graph catch numbers for each year to see if there is a dominate cohort
# loop that plots catch number and frequency in each year
plot_list <- list()
for(i in 1999:2009){
  dummy <- subset(myst_only, year == i)
  plot_list[[i]] <- ggplot(dummy, aes(x = fish_tl, y = count)) +
    geom_bar(stat= "identity") +
    xlab("size (cm)") +
    ylab("Number of sightings") +
    ggtitle(paste(i))
    
}
plot_list <- plot_list[1999:2009]
do.call(grid.arrange, plot_list)

rm(list ="i", "plot_list")

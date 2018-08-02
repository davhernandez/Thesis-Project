library(ggplot2)
library(ggmap)

#import raw data with lat lon
PISCO_sites <- read.csv("~/Desktop/Thesis/Raw Data/PISCO/PISCO Site Codes.csv")
CCFRP_sites <- read.csv("~/Desktop/Thesis/Raw Data/CCFRP/Hernandez_Drifts.csv")

#make dataframes with just lat lon
#the PISCO data table misnamed their lat lon columns
PISCO_gps <- as.data.frame(cbind(as.numeric(PISCO_sites$lon_wgs84), as.numeric(PISCO_sites$lat_wgs84), "PISCO"))
CCFRP_gps <- as.data.frame(cbind(CCFRP_sites$ST_LatDD, CCFRP_sites$ST_LonDD, "CCFRP"))
gps_points <- rbind(PISCO_gps, CCFRP_gps)
colnames(gps_points) <- c("lat", "lon", "survey")
#coercing the lat-lon columns from being type factor to type numeric
gps_points$lat <- as.numeric(as.character(gps_points$lat))
gps_points$lon <- as.numeric(as.character(gps_points$lon))
gps_points$survey <- as.character(gps_points$survey)
#removing duplicated gps points
gps_points <- unique(gps_points)
#removing NAs
gps_points <- na.omit(gps_points)
#remove the data points that look to be incorreclt entered
gps_points <- gps_points[which(gps_points$lon>-122.9),]
gps_points <- gps_points[which(gps_points$lat != 35.8871 & gps_points$lon != -120.89602),]
#cleaning up unused items
rm(list = "PISCO_gps", "CCFRP_gps")

#download the map
sitemap <- get_map(location = c(lon = mean(c(max(gps_points$lon), min(gps_points$lon))), lat = mean(c(max(gps_points$lat), min(gps_points$lat)))), zoom = 7, maptype = "satellite", scale = 2)

#plot gps points on the map
ggmap(sitemap) +
  geom_point(data = gps_points, aes(x = lon, y = lat, colour = survey)) +
  scale_color_discrete(name = "Dataset") +
  ggtitle("Location of sampling")

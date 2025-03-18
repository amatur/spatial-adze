# load required packages
library(ggplot2)
library(mapdata)

height <- 1
width <- 1

getX <- function(lon) {
  return(-(width/180)*(90-lon))
}
getY <- function(lat) {
  return((height/360.0)*(180+lat))
}

getLon <- function(x) {
  lon <- (180 * x + 180 * width / 2) / width
  return(lon)
}

getLat <- function(y) {
  lat <- (360.0 * y - height * 180) / height
  return(lat)
}


# Read the tab-separated file
data <- read.table("~/code/spatial-data/ara.tsv", header = TRUE, sep = "\t")
# Filter the "lat" column to have values less than 50
# data <- data[data$longitude > -10 ,  ]
# data <- data[data$longitude < 15  , ]
# data <- data[data$latitude < 60  , ]
# data <- data[data$latitude > 35  , ]
data <- data[data$longitude > -10 ,  ]
data <- data[data$longitude < 100  , ]
data <- data[data$latitude < 60  , ]
data <- data[data$latitude > -30  , ]

# Extract the column labeled "lat"
lat_column <- data$latitude
lat <- data$latitude
lon <- data$longitude
value <- data$id

# 
# x <- (getLon(getX(data$longitude)))
# y <- (getLat(getY(data$latitude)))
# plot(x, y, pch = 19, col = "blue", xlab = "x", ylab = "y", main = "Scatterplot of x vs y")
# # Add grid lines
# abline(h = seq(35, 60, by = 5), col = "lightgray", lty = "dotted")
# abline(v = seq(-10, 15, by = 5), col = "lightgray", lty = "dotted")
# 






# create example data with coordinates and values
# lon <- c(9.481544, 2.352222, -74.005973, 139.650312)
# lat <- c(51.312801, 48.856613, 40.712776, 35.676191)
# value <- c(12, 15, 19, 30)
foo1 <- data.frame(lon, lat, value)
foo1

# get world map data
world_map <- map_data("world")
head(world_map,10)

# create heat map using ggplot2
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "white", color = "grey") +
  geom_point(data = foo1, aes(x = lon, y = lat, fill = value), size = 1, shape = 21) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Longitude", y = "Latitude", title = "Example Heat Map")
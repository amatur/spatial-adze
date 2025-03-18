library(ggplot2)

MINX=-10
MAXX=15
MINY=35
MAXY=60
GRIDX=GRIDY=5

MINX=-10
MAXX=100
MINY=-30
MAXY=60
GRIDX=9
GRIDY=9

# Read the tab-separated file
data <- read.table("/Users/amatur/code/spatial-adze/src/output9.tsv", header = TRUE, sep = "\t")
data <- data[data$mean_alpha > 0,  ]

# Sample data
df <- data.frame(x = data$x,
                 y = data$y,
                 mean_alpha = data$mean_alpha)

# Define color scale
color_scale <- scale_color_gradient(low = "blue", high = "red")

# Plot
ggplot(df, aes(x = x, y = y, color = mean_alpha)) +
  geom_point() +
  scale_x_continuous(breaks = seq(MINX, MAXX, by =(MAXX-MINX)/GRIDX) ) +
  scale_y_continuous(breaks = seq(MINY, MAXY, by = (MAXY-MINY)/GRIDY) ) +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = "white")) +
  color_scale


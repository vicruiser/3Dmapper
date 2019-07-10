#Load libraries
library(mathart)
library(ggart)
library(ggforce)
library(Rcpp)
library(tidyverse)

# Example

points <- mathart::points

result <- kdtree(points, minmax = TRUE)

p <- ggplot() +
  geom_segment(aes(x, y, xend = xend, yend = yend), result) +
  coord_equal() +
  xlim(0, 10000) + ylim(0, 10000) +
  theme_blankcanvas(margin_cm = 0)


####Working directory----
setwd("SetYourPATH/DEanalysis/")
getwd()
#----

####Libraries----
library(ggplot2)
library(ggrepel)
library(tidyverse)
#----

####Data----
data <- read.csv2("./data/chick_cor.csv")
data$pattern <- as.factor(as.character(data$pattern))
data$GRAN <- as.factor(as.character(data$GRAN)) 
#----

####Axis---- 
# constants 
axis_begin_x  <- round(min(data$logFC.x))
axis_end_x  <- round(max(data$logFC.x))
total_ticks_x <- abs((axis_begin_x)) + abs((axis_end_x)) + 1

axis_begin_y  <- -round(min(data$logFC.y))
axis_end_y  <- round(max(data$logFC.x))
total_ticks_y <- abs((axis_begin_y)) + abs((axis_end_y)) + 1
# chart junk data
tick_frame_x <- 
  data.frame(ticks = seq(axis_begin_x, axis_end_x, length.out = total_ticks_x), 
             zero=0) %>%
  subset(ticks != 0)
  tick_frame_y <- 
  data.frame(ticks = seq(axis_begin_y, axis_end_y, length.out = total_ticks_y), 
             zero=0) %>%
  subset(ticks != 0)

lab_frame_x <- data.frame(lab = seq(axis_begin_x, axis_end_x),
                          zero = 0
) %>%
  subset(lab != 0)

lab_frame_y <- data.frame(lab = seq(axis_begin_y, axis_end_y),
                          zero = 0) %>%
  subset(lab != 0)

tick_sz_x <- (tail(lab_frame_y$lab, 1) -  lab_frame_y$lab[1]) / -90
tick_sz_y <- (tail(lab_frame_y$lab, 1) -  lab_frame_y$lab[1]) / -90
#----

####Plotting----
ggplot(data, aes(x=logFC.x, y=logFC.y))+
  geom_point(aes(size=GRAN,col=pattern))+
  scale_size_manual(values = c(6,10))+
  scale_color_manual(values = c("#898989","#ff0101","#0000e7"))+
  geom_point(data = data[data$GRAN == "SI",],
             pch = 21, fill=NA, size = 10, colour = "black", stroke=2.5)+
  
  #----

####Add ons----
# y axis line
geom_segment(x = 0, xend = 0, 
             y = lab_frame_y$lab[1], yend = tail(lab_frame_y$lab, 1),
             size = 1) +
  # x axis line
  geom_segment(y = 0, yend = 0, 
               x = lab_frame_x$lab[1], xend = tail(lab_frame_x$lab, 1),
               size = 1) +
  # x ticks
  geom_segment(data = tick_frame_x, 
               aes(x = ticks, xend = ticks, 
                   y = zero, yend = zero + tick_sz_x),
               size = 1) +
  # y ticks
  geom_segment(data = tick_frame_y, 
               aes(x = zero, xend = zero + tick_sz_y, 
                   y = ticks, yend = ticks),
               size = 1) + 
  # labels
  geom_text(data=lab_frame_x, aes(x=lab, y=zero, label=lab),vjust=1.8, size = 14) +
  geom_text(data=lab_frame_y, aes(x=zero, y=lab, label=lab),hjust=2, size = 14) +
  #----  

####Finishing the plot----
  theme_void()+
  theme(legend.position = "none")+

#Title
  geom_text(mapping = aes(label="Gallus gallus",x=-3,y=5,fontface=3),
            size = 18)
#----  

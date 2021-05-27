setwd("/home/frankaylward/recruitment_plots")
x <- read.table(file="testlastout.finalcoords.txt", header=T)
bounds <- read.table("testlastout.contig_boundaries.txt")

library(ggplot2)
ggplot(x) + geom_point(aes(x=coord, y=percid), colour="dodgerblue", alpha=0.3) + theme_bw() + geom_vline(xintercept = bounds$V2, linetype="solid", color = "grey50", size=0.5)

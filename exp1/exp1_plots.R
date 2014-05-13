# Chart Maker for SVM Experiment 1
# N. DiQuattro - Spring 2014

# Set-up ------------------------------------------------------------------

# Clear all
rm(list = ls())

# Set working directory
if (Sys.info()['sysname'] == "Darwin") {
  setwd("~/Dropbox/analysis/sac_svm/Rscripts/exp1")
} else {
  setwd('C:\\Dropbox\\analysis\\sac_svm\\Rscripts\\exp1')
}

# Load Libraries
library(ggplot2)

# Mean DV Plots --------------------------------------------------------------

# Load data
sdf <- read.csv("svm1_1sac_dat.txt", header = T)
  # Make Text Vectors
  sdf$dtype <- NA
  sdf$sdest <- NA
  sdf <- within(sdf,{
    dtype[ttype==1] = "Dissimilar"
    dtype[ttype==2] = "Similiar"
    sdest[endia==3] = "Distractor"
    sdest[endia==2] = "Target"
  }) 

# Find Subject Level means of DVs
mdat <- aggregate(cbind(slat, samp, fixdur, rt) ~ dtype + sdest, mean, data=sdf)

# Set-up shared plot theme
nogrid <- theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.title=element_blank(),
                )

poster.theme <- theme(legend.title=element_blank(),
                    legend.position = c(0, 1),
                    legend.justification = c(0, 1),
                    legend.text = element_text(size = 18),
                    axis.text = element_text(size = 24),
                    axis.title = element_text(size = 32),
                    axis.title.y = element_text(vjust = 0.3),
                    axis.title.x = element_text(vjust = 0.3)
                    )

# Set up shared Plot Values for DV means
pvals <- list(xlab("Saccade Destination"),
              stat_summary(fun.y=mean, geom="bar",
                           position=position_dodge(width=.9)),
              stat_summary(fun.data = mean_cl_normal, geom = "linerange",
                           position=position_dodge(width=0.9), show_guide=F),
              poster.theme
              )

# Make Plots
latplot <- ggplot(sdf,aes(sdest,slat,fill=dtype))+
                  ylab("Saccade Latency (ms)")+
                  coord_cartesian(ylim=c(210,280))+
                  pvals

ampplot <- ggplot(sdf,aes(sdest,samp,fill=dtype))+
                  ylab("Saccade Amplitude (Degs)")+
                  coord_cartesian(ylim=c(3.8,5.1))+
                  pvals

fixplot <- ggplot(sdf,aes(sdest,fixdur,fill=dtype))+
                  ylab("Fixation Duration (ms)")+
                  coord_cartesian(ylim=c(70,227))+
                  pvals

rtplot  <- ggplot(sdf,aes(sdest,rt,fill=dtype))+
                  ylab("Reaction Time (ms)")+
                  coord_cartesian(ylim=c(815,1082))+
                  pvals

# Display plots
latplot
ampplot
fixplot
rtplot

# Save Plots
psize <- 20
# ggsave("figs/lat.tiff", latplot, height = psize, width = psize, units = "cm",
#        dpi = 600)
# ggsave("figs/amp.tiff", ampplot, height = psize, width = psize, units = "cm",
#        dpi = 600)
# ggsave("figs/fix.tiff", fixplot, height = psize, width = psize, units = "cm",
#        dpi = 600)
#ggsave("figs/rt.png",rtplot)

# Classifier Plot ---------------------------------------------------------

# Load data
cdat <- read.csv("svm1_svmc_dat.txt", header=T)
meana <- data.frame(macc=mean(cdat$acc))  # Find mean

# Make plot
cplot <- ggplot(cdat, aes(snum,acc))+
                labs(x = "Subject", y = "Classifier Accuracy")+
                geom_bar(stat = "identity", alpha = .6)+
                geom_hline(aes(yintercept = .5, size = 1), linetype = "dotted",
                           size = 1)+
                geom_hline(aes(yintercept = mean(cdat$acc), size = .5))+
                coord_cartesian(ylim = c(.4, 1))+
                annotate("text", x = 22, y = mean(cdat$acc)-.01,
                         label = "Group Mean", color = "black", fontface = 2)+
                poster.theme

# Display and Save
cplot

# ggsave("figs/classacc.tiff",cplot, height = psize, width = psize, units = "cm",
#        dpi = 600)

# Saccade Proportion Plot -------------------------------------------------
  
# Load Data
spd <- read.csv("svm1_prop_dat.txt", header=T)
  # Make Text Vectors
  spd$dtype <- NA
  spd$sdest <- NA
  spd <- within(spd,{
    dtype[ttype==1] = "Dissimilar"
    dtype[ttype==2] = "Similiar"
    sdest[endia==3] = "Distractor"
    sdest[endia==2] = "Target"
  })

# Make Plot
prplot <- ggplot(spd, aes(sdest, props, fill=dtype))+
                 ylab("Proportion of First Saccade")+
                 pvals

#Display and Save
prplot

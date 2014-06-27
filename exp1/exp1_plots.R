# Chart Maker for SVM Experiment 1
# N. DiQuattro - Spring 2014 - Hi world!

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
library(dplyr)

# Mean DV Plots --------------------------------------------------------------

# Load data
sdf <- read.csv("svm1_1sac_dat.txt", header = T)
  # Make Text Vectors
  sdf$dtype <- NA
  sdf$sdest <- NA
  sdf <- within(sdf,{
    dtype[ttype==1] = "Dissimilar"
    dtype[ttype==2] = "Similar"
    dtype[ttype==3] = "Salient"
    sdest[endia==3] = "Distractor"
    sdest[endia==2] = "Target"
  }) 

# acc.dat <- read.csv("svm1_acc_dat.txt", header = TRUE)
# 
# sdf <- left_join(sdf, acc.dat, by=c("sub", "ttype", "endia"))

# Paper stats -------------------------------------------------------------

# Dissimilar t-tests
# t.test(rt ~ dtype, sdf, dtype!="Similar"&endia==3&sub!="10_tn", paired=TRUE)
# t.test(cor.y ~ dtype, sdf, dtype!="Similar"&endia==3&sub!="10_tn", paired=TRUE)

# # Behavoir tests
# t.test(rt ~ dtype, sdf, endia==3, paired=TRUE)
# t.test(cor.y ~ dtype, sdf, endia==3, paired=TRUE)

# Plot --------------------------------------------------------------------

# Set-up shared plot theme
# nogrid <- theme(panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.background = element_blank(),
#                 axis.line = element_line(colour = "black"),
#                 legend.title=element_blank()
#                 )
# 
# poster.theme <- theme(legend.title=element_blank(),
#                     legend.position = c(0, 1),
#                     legend.justification = c(0, 1),
#                     legend.text = element_text(size = 18),
#                     axis.text = element_text(size = 24),
#                     axis.title = element_text(size = 32),
#                     axis.title.y = element_text(vjust = 0.3),
#                     axis.title.x = element_text(vjust = 0.3)
#                     )

apa.theme <- theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.title=element_blank(),
                   axis.title = element_text(size = 12),
                   axis.title.y = element_text(vjust = 1),
                   axis.text = element_text(size = 10),
                   axis.ticks.x = element_blank() )

# Set up shared Plot Values for DV means
dcolors <- c("#984ea3", "#4daf4a")
pvals <- list(xlab("Distractor Type"),
              stat_summary(fun.y=mean, geom="bar"),
              stat_summary(fun.data = mean_cl_normal, geom = "linerange",
                           show_guide=FALSE),
              scale_fill_manual(values=dcolors, guide=FALSE),
              apa.theme)

# Subset data
pdat <- subset(sdf, sdest=="Distractor")

# Make Plots
latplot <- ggplot(pdat, aes(dtype, slat, fill=dtype))+
                  ylab("Saccade Latency (ms)")+
                  coord_cartesian(ylim=c(200, 275))+
                  scale_y_continuous(breaks=seq(200,275,10))+
                  pvals

ampplot <- ggplot(pdat,aes(dtype, samp, fill=dtype))+
                  ylab("Saccade Amplitude (Degs)")+
                  coord_cartesian(ylim=c(3.75, 5))+
                  scale_y_continuous(breaks=seq(3.75, 5, .2))+
                  pvals

fixplot <- ggplot(pdat, aes(dtype, fixdur, fill=dtype))+
                  ylab("Fixation Duration (ms)")+
                  coord_cartesian(ylim=c(50, 225))+
                  scale_y_continuous(breaks=seq(50, 225, 20))+
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
psize <- 3
ggsave("figs/lat_exp1.tiff", latplot, height=psize, width=psize, units="in",
       dpi=600)
ggsave("figs/amp_exp1.tiff", ampplot, height=psize, width=psize, units="in",
       dpi=600)
ggsave("figs/fix_exp1.tiff", fixplot,height=psize, width=psize, units="in",
       dpi=600)
# ggsave("figs/rt.png",rtplot)

# Classifier Plot ---------------------------------------------------------

# Load data
cdat <- read.csv("svm1_svmc_dat.txt", header=T)
con.dat <- read.csv("svm1_svmcontrol_dat.txt", header=T)
  con.dat <- con.dat[-1,]  # Remove subject without Sal captures


# Find means
cl.mn <- summarise(cdat,
            macc = mean(acc),
            sacc = sd(acc)
            )

conl.mn <- summarise(con.dat,
                   macc = mean(acc),
                   sacc = sd(acc) )

# Compare subject means to chance level
wt <- wilcox.test(cdat$acc, mu=.5, exact=FALSE)
  wt.z <- qnorm(wt$p.value/2)

wtc <- wilcox.test(con.dat$acc, mu=.5, exact=FALSE)
wtc.z <- qnorm(wtc$p.value/2)

# Weight Tests
wilcox.test(abs(cdat$slatw), abs(cdat$sampw), paired=TRUE, exact=FALSE)


# Make plot
cplot <- ggplot(cdat, aes(snum, acc))+
                labs(x="Subject", y="Classifier Accuracy")+
                geom_bar(stat="identity", alpha=.6)+
                geom_hline(aes(yintercept=.5, size=1), linetype="dotted", size=1)+
                geom_hline(aes(yintercept=mean(cdat$acc), size=.1))+
                coord_cartesian(ylim=c(.4, 1))+
                scale_x_discrete(breaks=seq(1, 24, 4))+
                apa.theme

# Display and Save
cplot

ggsave("figs/classacc_paper.tiff", cplot, height = 3, width = 3, units = "in",
       dpi = 600)

# Saccade Proportion Plot -------------------------------------------------
  
# # Load Data
# spd <- read.csv("svm1_prop_dat.txt", header=T)
#   # Make Text Vectors
#   spd$dtype <- NA
#   spd$sdest <- NA
#   spd <- within(spd,{
#     dtype[ttype==1] = "Dissimilar"
#     dtype[ttype==2] = "Similar"
#     sdest[endia==3] = "Distractor"
#     sdest[endia==2] = "Target"
#   })
# 
# # Find group means
# sp.mn <- spd %>%
#           group_by(ttype, endia) %>%
#           summarise(
#             props = mean(props))
# 
# # Stats
# spt  <- t.test(props ~ dtype, spd, endia==3, paired = TRUE)
# dist <- t.test(spd$props[spd$endia==3&spd$dtype=="Dissimilar"])
# simt <- t.test(spd$props[spd$endia==3&spd$dtype=="Similar"])
# 
# # Make Plot
# prplot <- ggplot(spd, aes(sdest, props, fill=dtype))+
#                  ylab("Proportion of First Saccade")+
#                  pvals
# 
# #Display and Save
# prplot

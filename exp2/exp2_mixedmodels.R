# SVM 3 - Mixed Models
# N. DiQuattro - Spring 2014

# Set-up ------------------------------------------------------------------

# Clear all
rm(list = ls())

# Set working directory
if (Sys.info()['sysname'] == "Darwin") {
  setwd("~/Dropbox/analysis/sac_svm/Rscripts/exp3")
} else {
  setwd('C:\\Dropbox\\analysis\\sac_svm\\Rscripts\\exp3')
}

# Load Libraries
#library(lme4)
#library(lmerTest)
library(ggplot2)
library(dplyr)

# Load Data ---------------------------------------------------------------
soast = c("-100ms", "-68ms", "-33ms", "0ms")
dtypest = c("Similar", "Dissimilar")
nendiast = c("Fixation", "Target", "Distractor")

sdf <- read.csv("svm3_1sac_dat2.txt", header = T) %>%
        mutate(
          dtype = dtypest[ttype],
          soast = soast[soa],
          soa = soa - 1,
          nendia = nendiast[nendia],
          nendia = ifelse(is.na(nendia), "NoSac2", nendia) ) %>%
        filter(samp < 7) # Remove trials with weird amplitdues

# Behavioral Stats --------------------------------------------------------

# # Calculate means
bdat <- sdf %>%
          group_by(dtype, soast, sub) %>%
          summarise(
            acc = mean(cor, na.rm = TRUE),
            rt = mean(rt[cor==1], na.rm = TRUE) )
          
# 
# obdat <- summarise(bdat,
#                    acc = mean(acc),
#                    rt = mean(rt) ) %>%
#          summarise(
#            acc = mean(acc),
#            rt = mean(rt) )
# 
# bdat <- as.data.frame(bdat)
# 
# plot it to check it out
ggplot(bdat, aes(soast, acc, color=dtype, group=dtype)) + 
        stat_summary(fun.y=mean, geom="line", size = 2) +
        stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size=2)
# 
# # Accuracy
#   # mixed model
#   acc.mod <- lmer(acc ~ dtype*soa + (dtype*soa|sub), bdat)
#     acc.sum <- summary(acc.mod)
#     acc.aov <- anova(acc.mod)
# 
#   # Normal RM ANOVA
#   bdat$soaF <- as.factor(bdat$soa)  # make factor vector
# 
#   acc.aovrm <- aov(acc ~ dtype*soaF + Error(sub/(dtype*soaF)), bdat)
#   acc.aovrm.sum <- summary(acc.aovrm)
# 
# # Reaction Time
#   # RM ANOVA
#   rt.aovrm <- aov(rt ~ dtype*soaF + Error(sub/(dtype*soaF)), bdat)
#     rt.aov.rm.sum <- summary(rt.aovrm)


# Summary Stats for Eye Data ----------------------------------------------

sdat <- sdf %>%
          filter(cor == 1) %>%
          group_by(dtype, soast, sub) %>%
          summarise(
            slat   = mean(slat),
            samp   = mean(samp),
            fixdur = mean(fixdur, na.rm=TRUE) )

# Shared Plot Values
  # Theme Creation
#   nogrid <- theme(panel.grid.major = element_blank(),
#                   panel.grid.minor = element_blank(),
#                   panel.background = element_blank(),
#                   axis.line = element_line(colour = "black"),
#                   legend.title=element_blank())
# 
#   poster.theme <- theme(legend.title=element_blank(),
#                         legend.position = c(0, 1),
#                         legend.justification = c(0, 1),
#                         legend.text = element_text(size = 18),
#                         axis.text = element_text(size = 24),
#                         axis.title = element_text(size = 32),
#                         axis.title.y = element_text(vjust = 0.3),
#                         axis.title.x = element_text(vjust = 0.3) )

  apa.theme <- theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     legend.title=element_blank(),
                     legend.position="none",
                     axis.title = element_text(size = 12),
                     axis.title.y = element_text(vjust = 1),
                     axis.text = element_text(size = 10),
                     axis.ticks.x = element_blank() )


  # Plot values
  dcolors <- c("#984ea3", "#4daf4a")
  pvals <- list(scale_x_continuous(labels=c("-100", "-68", "-33", "0"),
                                   name="Distractor/Target SOA (ms)"),
                stat_summary(fun.y=mean, geom="line", size=.7),
                stat_summary(fun.data=mean_cl_normal, geom="pointrange",
                              size=.7),
                scale_color_manual(values=dcolors) )

# Make Plots
lat.plot <- ggplot(sdat, aes(soa, slat, color=dtype))+
                  ylab("Saccade Latency (ms)")+
                  coord_cartesian(ylim=c(200, 275))+
                  scale_y_continuous(breaks=seq(200,275,10))+
                  pvals+
                  apa.theme

amp.plot <- ggplot(sdat, aes(soa, samp, color=dtype))+
                  ylab("Saccade Amplitude (deg)")+
                  coord_cartesian(ylim=c(3.75, 5))+
                  scale_y_continuous(breaks=seq(3.75, 5, .2))+
                  pvals+
                  apa.theme

fix.plot <- ggplot(sdat, aes(soa, fixdur, color=dtype))+
                  ylab("Fixation Duration (ms)")+
                  coord_cartesian(ylim=c(50, 225))+
                  scale_y_continuous(breaks=seq(50, 225, 20))+
                  pvals+
                  apa.theme

# Display and Save
lat.plot
amp.plot
fix.plot

# Try out a desnity plot
# lat.den <- ggplot(sdf, aes(x=rt, fill=dtype, alpha=.6))+
#                   geom_density()+
#                   facet_wrap(~soa)

# kPsize <- 3
# ggsave("figs/lat.tiff", lat.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)
# ggsave("figs/amp.tiff", amp.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)
# ggsave("figs/fix.tiff", fix.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)

# kPsize <- 6
# ggsave("figs/lat_small.tiff", lat.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)
# ggsave("figs/amp_small.tiff", amp.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)
# ggsave("figs/fix_small.tiff", fix.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)


# Fit mixed models --------------------------------------------------------
# # Remove incorrect trials
# sdf <- subset(sdf, cor==1)
# 
# # Increase interations
# iters <- lmerControl(optCtrl=list(maxfun=1000000))
# #, control=iters
# 
# # Latency Models
# sdf$dtype <- relevel(sdf$dtype, "Similar")
# slat.mod <- lmer(slat ~ dtype*soa + (dtype*soa|sub), sdf)
#   slat.sum <- summary(slat.mod)
#   slat.aov <- anova(slat.mod)
# 
#   # Diagnostics
#   plot(slat.mod)
#   qqnorm(residuals(slat.mod))
#   qqline(residuals(slat.mod))
#   hist(residuals(slat.mod))
#   ggplot(slat.mod, aes(soa, .resid, color=dtype))+
#     geom_point(position="jitter")+
#     geom_hline(yintercept=0)
# 
# # Amplitude Models
# amp.mod <- lmer(samp ~ dtype*soa + (dtype*soa|sub), sdf)
#   amp.sum <- summary(amp.mod)
#   amp.aov <- anova(amp.mod)
# 
#   # Diagnostics
#   plot(amp.mod)
#   qqnorm(residuals(amp.mod))
#   qqline(residuals(amp.mod))
#   hist(residuals(amp.mod))
#   ggplot(amp.mod, aes(soa, .resid, color=dtype))+
#     geom_point(position="jitter")+
#     geom_hline(yintercept=0)
# 
# # Fixation Duration model
# fix.mod <- lmer(fixdur ~ dtype*soa + (dtype*soa|sub), sdf)
#   fix.sum <- summary(fix.mod)
#   fix.aov <- anova(fix.mod)
#   
#   # Diagnostics
#   plot(fix.mod)
#   qqnorm(residuals(fix.mod))
#   qqline(residuals(fix.mod))
#   hist(residuals(fix.mod))
#   ggplot(fix.mod, aes(soa, .resid, color=dtype))+
#     geom_point(position="jitter")+
#     geom_hline(yintercept=0)


# Second saccade stuff ----------------------------------------------------

dat2 <- sdf %>%
          filter(cor == 1, endia == 3) %>%
          group_by(dtype, soast, nendia, sub) %>%
          summarise(
            slat = mean(slat),
            samp = mean(samp),
            fix  = mean(fixdur, na.rm=TRUE),
            tnum = n()) %>%
          filter(nendia != "Distractor") %>%
          filter(nendia != "Fixation")

ggplot(dat2, aes(soast, slat, color=dtype))+
  stat_summary(fun.data=mean_cl_normal, geom="pointrange", size=1)+
  facet_wrap(~nendia)

fmod <- lmer(fixdur ~ ttype*soa + (ttype*soa|sub), dat2)
summary(fmod)

# Saccade Proportion Data -------------------------------------------------

# Load Data
soast = c("-100ms", "-68ms", "-33ms", "0ms")
dtypest = c("Similar", "Dissimilar")
sprop <- read.csv("svm3_sacprop_dat.txt", header=TRUE,
                  stringsAsFactors=FALSE) %>%
          mutate(
            dtype = dtypest[ttype],
            soast = soast[soa],
            soa = soa - 1)

# Find means
pmeans <- sprop %>%
            group_by(dtype, soa, endia) %>%
            summarise(
              mpro = mean(props))

# plot
ggplot(pmeans, aes(soa, mpro, group=dtype, color=dtype)) + geom_line()

ggplot(sprop[sprop$endia==2,], aes(soa, props, color=dtype))+
  ylab("Proportion of Capture")+
 # scale_y_continuous(breaks=seq(0,1,.1))+
  pvals+
  facet_wrap(~endia)
 # apa.theme
  #theme(legend.position="bottom")

# ggsave("figs/props.tiff", pro.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)

# # Try stats
# prop.mod <- lmer(props ~ dtype*soa + (dtype*soa|sub), sprop)
# 
# # RM ANOVA
# prop.aov <- aov(props ~ dtype*soa + Error(sub / (dtype*soa)), sprop)

# Correlation with eye mets
cordat <- left_join(sprop, sdat, by=c("sub", "dtype", "soast")) %>%
          select(sub, dtype, soast, props, slat, samp, fixdur) %>%
          mutate(
            soast = factor(soast, c("-100ms", "-68ms", "-33ms", "0ms")))

# Scatter plots with fits
ggplot(cordat, aes(samp, fixdur)) +
  geom_point()+
  geom_smooth(method=lm)+
  facet_grid(dtype ~ soast)

cor.test(~fixdur + samp, cordat)

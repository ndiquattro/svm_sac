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
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)

# Load Data ---------------------------------------------------------------

sdf <- read.csv("svm3_1sac_dat.txt", header = T)
  # Make Text Vectors
  sdf$dtype[sdf$ttype==1] <- "Similar"
  sdf$dtype[sdf$ttype==2] <- "Dissimilar"

# Define factors as factors and center SOA to -103ms
sdf <- within(sdf,{
              dtype = factor(dtype, levels=c("Dissimilar","Similar"))
              sub = factor(sub)
              soa = soa - 1
})

# Outlier Removal
sdf <- subset(sdf, samp<7)  # Remove trials with weird amplitdues


# Behavioral Stats --------------------------------------------------------

# Calculate means
bdat <- sdf %>%
          group_by(dtype, soa, sub) %>%
          summarise(
            acc = mean(cor, na.rm = TRUE),
            rt = mean(rt[cor==1], na.rm = TRUE) )

obdat <- summarise(bdat,
                   acc = mean(acc),
                   rt = mean(rt) ) %>%
         summarise(
           acc = mean(acc),
           rt = mean(rt) )

bdat <- as.data.frame(bdat)

# plot it to check it out
ggplot(bdat, aes(soa, acc, color=dtype)) + 
        stat_summary(fun.y=mean, geom="line", size = 2) +
        stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size=2)

# Accuracy
  # mixed model
  acc.mod <- lmer(acc ~ dtype*soa + (dtype*soa|sub), bdat)
    acc.sum <- summary(acc.mod)
    acc.aov <- anova(acc.mod)

  # Normal RM ANOVA
  bdat$soaF <- as.factor(bdat$soa)  # make factor vector

  acc.aovrm <- aov(acc ~ dtype*soaF + Error(sub/(dtype*soaF)), bdat)
  acc.aovrm.sum <- summary(acc.aovrm)

# Reaction Time
  # RM ANOVA
  rt.aovrm <- aov(rt ~ dtype*soaF + Error(sub/(dtype*soaF)), bdat)
    rt.aov.rm.sum <- summary(rt.aovrm)


# Summary Stats for Eye Data ----------------------------------------------

# Remove incorrect trials
sdf <- subset(sdf, cor==1)

# Find means
sdat <- aggregate(cbind(slat,samp,fixdur) ~ dtype+soa+sub, FUN=mean, sdf)
gdat <- aggregate(cbind(slat,samp,fixdur) ~ dtype+soa, FUN=mean, sdat)

# Shared Plot Values
  # Theme Creation
  nogrid <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.title=element_blank())

  poster.theme <- theme(legend.title=element_blank(),
                        legend.position = c(0, 1),
                        legend.justification = c(0, 1),
                        legend.text = element_text(size = 18),
                        axis.text = element_text(size = 24),
                        axis.title = element_text(size = 32),
                        axis.title.y = element_text(vjust = 0.3),
                        axis.title.x = element_text(vjust = 0.3))


  # Plot values
  pvals <- list(scale_x_continuous(labels=c("-100","-68","-33","0"),
                                   name="Distractor -> Target SOA"),
                stat_summary(fun.y=mean, geom="line", size = 2),
                stat_summary(fun.data = mean_cl_normal, geom = "pointrange",
                             size=2),
                poster.theme)

# Make Plots
lat.plot <- ggplot(sdat, aes(soa, slat, color=dtype))+
                  ylab("Saccade Latency (ms)") + pvals

amp.plot <- ggplot(sdat, aes(soa, samp, color=dtype))+
                  ylab("Saccade Amplitude (deg)") + pvals+
                  theme(legend.position = c(-1,-1))

fix.plot <- ggplot(sdat, aes(soa, fixdur, color=dtype))+
                  ylab("Fixation Duration (ms)") + pvals+
                  theme(legend.position = c(-1,-1))

# Display and Save
lat.plot
amp.plot
fix.plot

# Try out a desnity plot
lat.den <- ggplot(sdf, aes(x=rt, fill=dtype, alpha=.6))+
                  geom_density()+
                  facet_wrap(~soa)

# kPsize <- 20
# ggsave("figs/lat.tiff", lat.plot, height = kPsize, width = kPsize, units = "cm",
#        dpi = 600)
# ggsave("figs/amp.tiff", amp.plot, height = kPsize, width = kPsize, units = "cm",
#        dpi = 600)
# ggsave("figs/fix.tiff", fix.plot, height = kPsize, width = kPsize, units = "cm",
#        dpi = 600)

# kPsize <- 6
# ggsave("figs/lat_small.tiff", lat.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)
# ggsave("figs/amp_small.tiff", amp.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)
# ggsave("figs/fix_small.tiff", fix.plot, height = kPsize, width = kPsize, units = "in",
#        dpi = 600)


# Fit mixed models --------------------------------------------------------
# Remove incorrect trials
sdf <- subset(sdf, cor==1)

# Increase interations
iters <- lmerControl(optCtrl=list(maxfun=1000000))
#, control=iters

# Latency Models
sdf$dtype <- relevel(sdf$dtype, "Similar")
slat.mod <- lmer(slat ~ dtype*soa + (dtype*soa|sub), sdf)
  slat.sum <- summary(slat.mod)
  slat.aov <- anova(slat.mod)

  # Diagnostics
  plot(slat.mod)
  qqnorm(residuals(slat.mod))
  qqline(residuals(slat.mod))
  hist(residuals(slat.mod))
  ggplot(slat.mod, aes(soa, .resid, color=dtype))+
    geom_point(position="jitter")+
    geom_hline(yintercept=0)

# Amplitude Models
amp.mod <- lmer(samp ~ dtype*soa + (dtype*soa|sub), sdf)
  amp.sum <- summary(amp.mod)
  amp.aov <- anova(amp.mod)

  # Diagnostics
  plot(amp.mod)
  qqnorm(residuals(amp.mod))
  qqline(residuals(amp.mod))
  hist(residuals(amp.mod))
  ggplot(amp.mod, aes(soa, .resid, color=dtype))+
    geom_point(position="jitter")+
    geom_hline(yintercept=0)

# Fixation Duration model
fix.mod <- lmer(fixdur ~ dtype*soa + (dtype*soa|sub), sdf)
  fix.sum <- summary(fix.mod)
  fix.aov <- anova(fix.mod)
  
  # Diagnostics
  plot(fix.mod)
  qqnorm(residuals(fix.mod))
  qqline(residuals(fix.mod))
  hist(residuals(fix.mod))
  ggplot(fix.mod, aes(soa, .resid, color=dtype))+
    geom_point(position="jitter")+
    geom_hline(yintercept=0)


# Saccade Proportion Data -------------------------------------------------

# Load Data
sprop <- read.csv("svm3_sacprop_dat.txt", header=TRUE, stringsAsFactors=FALSE)
  # Make Text Vectors
  sprop$dtype[sprop$ttype==1] <- "Similar"
  sprop$dtype[sprop$ttype==2] <- "Dissimilar"

# Set data up
sprop <- within(sprop, {
                  dtype = factor(dtype, levels=c("Dissimilar", "Similar"))
                  soa = factor(soa)
                  sub = factor(sub)
                  }
                )

# Find means
pmeans <- sprop %>%
            group_by(dtype, soa) %>%
            summarise(
              mpro = mean(props))

# plot
ggplot(pmeans, aes(soa, mpro, group=dtype, color=dtype)) + geom_line()

# Try stats
prop.mod <- lmer(props ~ dtype*soa + (dtype*soa|sub), sprop)

# RM ANOVA
prop.aov <- aov(props ~ dtype*soa + Error(sub / (dtype*soa)), sprop)

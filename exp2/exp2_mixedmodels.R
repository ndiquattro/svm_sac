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

# Functions ---------------------------------------------------------------

ForDiffMaker <- function(k) {
  # Makes a contrast matrix for forward difference coding used in regression
  # Args:
  #   k: number of levels in factor
  # Returns:
  #   cmat: a matrix
  cmat = matrix(nrow = k, ncol = k -1)  # Create an empty matrix
    for (i in 1:ncol(cmat)) {
      cmat[1:i,i] = k-i
      cmat[(i+1):nrow(cmat),i] = i * -1
    }
  cmat = cmat / k
  return(cmat)
 }   
OutlierMarker = function(df, form, cut) {
  # Finds outliers in a dataframe
  # Args:
  #   df: data frame
  #   form: formula for subsetting the data
  #   cut: Number of Standard Deviations to set threshold at
  # Returns:
  #   lv: logical vector of outlier presence
  
  # Decompose formula
  fvars <- all.vars(form)
  
  # Find means and SDs based on formula
  mn.dat <- aggregate(form, df, mean)
  sd.dat <- aggregate(form, df, sd)
    # Merge stat data
    outs <- merge(mn.dat, sd.dat, by=fvars[2:length(fvars)], suffixes=c(".mn",
                                                                        ".sd"))

  # Make upper and lower thresholds based on cut
  odim <- dim(outs)
  outs$lpc = outs[,odim[2]-1] - (outs[,odim[2]] * cut)
  outs$upc = outs[,odim[2]-1] + (outs[,odim[2]] * cut)

  # Make index for orignal dataframe where values are above or below thresholds
  cnames <- colnames(outs)
  lv <- logical(length=dim(df)[1])
  
  for (a in 1:odim[1]) {
    cur.row <- outs[a,]
    cur.log <- with(df, {
                    cnames[1] == cur.row[cnames[1]] &
                    cnames[2] == cur.row[cnames[2]] &
                    cnames[3] == cur.row[cnames[3]] &(
                    fvars[1] < cur.row$lpc |
                    fvars[1] > cur.row$upc)
    })
    #cur.log <- df[,cnames[1]] == currow[,cnames[1]]&
    #df[,cnames[2]] == currow[,cnames[2]] &df[,cnames[3]] == currow[,cnames[3]]&
    #(df[,fvars[1]] < currow$lpc | df[,fvars[1]] > currow$upc)
    lv <- lv + cur.log
  }
  
  return(as.logical(lv))
 }

# Load Data ---------------------------------------------------------------

sdf <- read.csv("svm3_1sac_dat.txt", header = T)
  # Make Text Vectors
  sdf$dtype[sdf$ttype==1] <- "Similar"
  sdf$dtype[sdf$ttype==2] <- "Dissimilar"

# Define factors as factors and center SOA to -103ms
sdf <- within(sdf,{
              dtype = factor(dtype,levels=c("Dissimilar","Similar"))
              sub = factor(sub)
              soa = soa - 1
})

# Outlier Removal
sdf <- subset(sdf,samp<7)  # Remove trials with weird amplitdues
fout <- OutlierMarker(sdf, fixdur~dtype+soa+sub, 2.5) # Fix duraiton outliers
sout <- OutlierMarker(sdf, slat~dtype+soa+sub, 2.5)  # Saccade latency outliers
#sdf <- subset(sdf, !sout | !fout)  # Remove outliers

# Plot Means --------------------------------------------------------------

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
                                   name="Distractor -> Target SOA (ms)"),
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

kPsize <- 20
ggsave("figs/lat.tiff", lat.plot, height = kPsize, width = kPsize, units = "cm",
       dpi = 600)
ggsave("figs/amp.tiff", amp.plot, height = kPsize, width = kPsize, units = "cm",
       dpi = 600)
ggsave("figs/fix.tiff", fix.plot, height = kPsize, width = kPsize, units = "cm",
       dpi = 600)

# Fit mixed models --------------------------------------------------------

# Increase interations
iters <- lmerControl(optCtrl=list(maxfun=1000000))

# Latency Models
slat.mod <- lmer(slat ~ dtype*soa + (dtype*soa|sub), sdf, control=iters)
  slat.sum <- summary(slat.mod)

# Amplitude Models
amp.mod <- lmer(samp ~ dtype*soa + (dtype*soa|sub), sdf, control=iters)
  amp.sum <- summary(amp.mod)
#ampnull <- lmer(samp ~ 1 + (dtype*soa|sub), sdf, control=iters)  # Null model

# Fixation Duration model
fix.mod <- lmer(fixdur ~ dtype*soa + (dtype*soa|sub), sdf, control=iters)
  fix.sum <- summary(fix.mod)

# Parameter Confidence Intervals ------------------------------------------
  # Talked to Jon about this method, technically legit, but not ideal for
  # comparing when you have different number of samples per condition.
  # best thing is to flip reference level and re-run the model.

#   # Find CIs for each model
#   cimeth <- "Wald"  # Use "boot" for best
#   slat.ci <- confint.merMod(slat.mod, method=cimeth)
#   amp.ci <- confint.merMod(amp.mod, method=cimeth)
#   fix.ci <- confint.merMod(fix.mod, method=cimeth)
# 
#   # CI plot function
#   PlotCI <- function(mod, cidf, dv){
#     # Plots model parameters with Confidence Intervals
#     # NOTE: There is specific code to this experiment in here.
#     # Args:
#     #   mod: Fitted lmer model
#     #   cidf: results from CI function
#     #   dv: Name of Depdent Variable
#     # Returns:
#     #   A ggplot
# 
#     # Get fixed effects and merge with CIs
#     fx <- fixef(mod)
#     asig <- merge(fx, cidf, by=0)  # 0 for by means by row names
#       colnames(asig) <- c("Param","Mean","Low","Up")
#       asig$comp <- c("Distractor Type","Distractor Type","SOA","SOA")
#       asig$Param <- as.character(asig$Param)
# 
#     # Calculate values by adding Means to CIs
#     asig[2,2:4] <- asig[1,2] + asig[2,2:4]
#     asig[3,2:4] <- asig[4,2] + asig[3,2:4]
# 
#     # Create plot
#     ldat <- data.frame(z=0,comp="SOA")
#     ci.plot <- ggplot(asig, aes(x=0, y=Mean, color=Param, ymin=Low, ymax=Up))+
#                       scale_color_discrete(labels=c("Dis Intercept","Sim Intercept","Sim Slope","Dis Slope"))+
#                       geom_pointrange(size=2.5, alpha=.5)+
#                       geom_hline(aes(yintercept=z), ldat)+
#                       facet_wrap(~comp, scales="free")+ 
#                       labs(y=dv,x=NULL)+
#                       theme(legend.title=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
#     return(ci.plot)
#   }
# 
#   # Plot CIs for each model
#   PlotCI(slat.mod, slat.ci, "Saccade Latency (ms)")
#   PlotCI(amp.mod, amp.ci, "Saccade Amplitude (degs)")
#   PlotCI(fix.mod, fix.ci, "Fixation Duration (ms)")

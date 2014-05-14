# sac_svm 5
# N. DiQuattro - Spring 2014

# Oms SOA Key
#   1_0: Target and Sim
#   2_0: Target and Dsim
#   3_0: Sim and Dsim
#   4_0: Dsim and Dsim

# Set-up ------------------------------------------------------------------

# Clear all
rm(list = ls())

# Set working directory
if (Sys.info()['sysname'] == "Darwin"){
  setwd("~/Dropbox/analysis/sac_svm/Rscripts/exp5")
}else{
  setwd('C:\\Dropbox\\analysis\\sac_svm\\Rscripts\\exp5')
}

# Load Libraries
library(ggplot2)

# Load Data ---------------------------------------------------------------

dat <- read.csv("svm_5_dat.txt", header = TRUE, na.strings = c(".", "NA"),
                stringsAsFactors = FALSE, strip.white = T, 
                col.names = c("sub", "tnum", "pupil", "on1", "on2", "soa",
                              "tloc", "dcol", "tcol", "cor", "rt", "sidx",
                              "endia", "stia", "samp", "slat", "fixdur",
                              "pendia", "pstia"))

# Clean up data
dat <- within(dat, {
#               endia = substr(endia, 1, nchar(endia)-1)  # Remove extra space 
#               stia  = substr(stia, 1, nchar(stia)-1)
#               pendia = substr(pendia, 1, nchar(pendia)-1)
#               pstia  = substr(pstia, 1, nchar(pstia)-1)
              # Set no target trials to correct if no response
              cor[on1 > 1 & on2 > 1 & rt == 0] = 1
              cor[on1 > 1 & on2 == 0 & rt == 0] = 1
              soa[soa==1] = "-68ms"  # Make the soa vector informative
              soa[soa==0] = "0ms"
            })

# Collapse across some conditions
dat$endia[dat$endia=="OpDcol"] = "Dcol"

# Remove subjects
dat <- subset(dat, !sub %in% c("24_lp",  # No D sacs
                               "25_rs"   # No D sacs
                               ))

# Subject and Group Means -------------------------------------------------

# Make datasets for analysis
bhav <- subset(dat, sidx==1 & pupil==0 & samp < 7 & endia != "fix")
fsac <- subset(dat, sidx==1 & pupil==0 & samp < 7 & endia != "fix" & cor==1)
sac2 <- subset(dat, sidx==2 & pupil==0 & cor==1 & pstia=="fix" &
                    (pendia=="Distractor" | pendia=="On1"))

# Behavioral means
oacc.mn <- aggregate(cor ~ sub, bhav, mean)  # Overall Accuracy
sub.acc <- aggregate(cor ~ on1 + on2 + endia + soa + sub, bhav, mean)
  grp.acc <- aggregate(cor ~ on1 + on2 + endia + soa, sub.acc, mean)

# Subject level means of DVs
sub.mn <- aggregate(cbind(slat, samp, fixdur) ~ on1 + on2 + endia + soa + sub,
                    fsac, mean)
  # Make vectors for plotting
  sub.mn$cond   = paste(sub.mn$on1, sub.mn$on2, sub.mn$endia, sep="_")
  sub.mn$onsets = paste(sub.mn$on1, sub.mn$on2, sep="_")

# Group means
grp.mn <- aggregate(cbind(slat, samp, fixdur) ~ cond, sub.mn, mean)

# Trial Numbers
sub.ln <- aggregate(slat ~ on1 + on2 + endia + soa + sub, fsac, length)
  sub.ln <- within(sub.ln, {
                   cond   = paste(on1, on2, endia, sep="_")
                   onsets = paste(on1, on2, sep="_")
                 })

# Saccade Proportions
num <- aggregate(slat ~ on1 + on2 + endia + sub + soa, fsac, length)
  num$cond <- paste(num$on1, num$on2, num$endia, sep="_")
  num$onsets <- paste(num$on1, num$on2, sep="_")

den <- aggregate(slat ~ on1 + on2 + sub, fsac, length)
  # Give up and write a loop
  num$props = NA
  for (a in 1:dim(num)[1]) {
    num$props[a] = num[a, "slat"] / den[den$on1==num[a, "on1"] &
                                    den$on2==num[a, "on2"]& 
                                    den$sub==num[a, "sub"], "slat"]
  }

# Find 2nd saccade proportions
num2 <- aggregate(slat ~ on1 + on2 + endia + sub + soa, sac2, length)
  num2$cond <- paste(num2$on1, num2$on2, num2$endia, sep="_")
  num2$onsets <- paste(num2$on1, num2$on2, sep="_")

den2 <- aggregate(slat ~ on1 + on2 + sub, sac2, length)
  # Give up and write a loop
  num2$props = NA
  for (a in 1:dim(num2)[1]) {
    num2$props[a] = num2[a, "slat"] / den2[den2$on1==num2[a, "on1"] &
                                        den2$on2==num2[a, "on2"]& 
                                        den2$sub==num2[a, "sub"], "slat"]
}      
# Stats -------------------------------------------------------------------
  
# Subset for Stats
sdat <- subset(sub.mn, cond == "3_1_Distractor" | cond == "3_3_On1" |
                       cond == "2_1_Distractor" | cond == "2_2_On1")

# Anovas
  slat.mod <- aov(slat ~ cond + Error(sub/cond), sdat)
    summary(slat.mod)
  samp.mod <- aov(samp ~ cond + Error(sub/cond), sdat)
    summary(samp.mod)
  fix.mod <- aov(fixdur ~ cond + Error(sub/cond), sdat)
    summary(fix.mod)

  # Post-hoc t-tests
  slat.t <- with(sdat, pairwise.t.test(slat, cond, "bon", paired=TRUE))
  amp.t <- with(sdat, pairwise.t.test(samp, cond, "bon", paired=TRUE))
  fix.t <- with(sdat, pairwise.t.test(fixdur, cond, "bon", paired=TRUE))

# Plot it! ----------------------------------------------------------------

# Which data for plotting?
plot.dat <- sdat
  # Make vector for plotting based on Sac destination
  plot.dat$fillon <- paste(plot.dat$on1)
    plot.dat$fillon[plot.dat$fillon=="4"] = "3"
    plot.dat$fillon[plot.dat$fillon=="3" & plot.dat$endia=="Sim"] = "2"
    plot.dat$fillon[plot.dat$fillon=="2" & 
                    plot.dat$on2 == 0 &
                    plot.dat$endia=="Distractor"] = "3"
    plot.dat$fillon[plot.dat$fillon=="1"] = "2"
  # Make factor to so we can define level order
  plot.dat$fillon = factor(plot.dat$fillon, levels=c("3", "2"))

  # Try something tricky to get x labels to be right
  plot.dat <- within(plot.dat, {
    cond[cond=="2_1_Distractor"] = "Target"
    cond[cond=="2_2_On1"] = "Similar"
    cond[cond=="3_1_Distractor"] = "Target "
    cond[cond=="3_3_On1"] = "Dissimilar"
    cond[cond=="1_0_Distractor"] = "Target  "
    cond[cond=="4_0_Dcol"] = "Dissimilar "
    cond = factor(cond, levels=c("Target", "Similar", "Target ", "Dissimilar",
                                 "Target  ", "Dissimilar "))
  })

# Shared Pvals
ptheme <- theme(axis.text.x = element_text(angle=-45, hjust=0, size=12,
                                           color="black"),
                axis.title.x = element_blank())

poster.theme <- theme(legend.position = c(0, 1),
                      legend.justification = c(0, 1),
                      legend.text = element_text(size = 18),
                      axis.text = element_text(size = 24),
                      axis.title = element_text(size = 32),
                      axis.title.y = element_text(vjust = 0.3),
                      axis.title.x = element_text(vjust = 0.3),
                      #axis.text.x = element_text(angle=-45, hjust=0),
                      axis.title.x = element_blank())

pvals <- list(stat_summary(fun.y=mean, geom="bar", position="dodge"),
              stat_summary(fun.data=mean_cl_normal, geom="linerange"),
              scale_x_discrete("Second Onset"),
              scale_fill_discrete("Saccade Destination", labels=c("Dissimilar",
                                                                  "Similar")),
              poster.theme)

# Diagnostic Plot
tnum.plot <- ggplot(sub.ln, aes(cond, slat, fill=onsets))+
                    scale_y_continuous("Mean Number of Trials", breaks=1:40)+
                    geom_boxplot()+ ptheme
                    #pvals

# Saccade Proportion plot
spro.plot <- ggplot(num, aes(cond, props, fill=onsets))+
              scale_y_continuous("Proportion of Saccades", breaks=seq(0,1,.1))+
              pvals + ptheme

# 2nd saccades proportion plot
pdat2 <- subset(num2, cond == "3_1_Target" | cond == "3_3_On2" |
                      cond=="2_2_On2")
ggplot(pdat2, aes(cond, props, fill=onsets)) + labs(title="2nd Saccades")+
  scale_y_continuous("Proportion of Saccades", breaks=seq(0,1,.1))+
  pvals

# Eye DVs
slat.plot <- ggplot(plot.dat, aes(cond, slat, fill=fillon))+
                    ylab("Saccade Latency (ms)")+
                    coord_cartesian(ylim=c(170,290))+
                    annotate("text", 2, 270, label="*", size=24)+
                    annotate("segment", x=1, xend=3, y=269, yend=269, size=1.5)+      
                    pvals

amp.plot <- ggplot(plot.dat, aes(cond, samp, fill=fillon))+
                  ylab("Saccade Amplitude (deg)")+
                  coord_cartesian(ylim=c(3, 5.1))+
                  pvals+
                  theme(legend.position = c(-1,-1))+
                  annotate("segment", x=c(.9, 2.5, 2.5), xend=c(2.5, 2.5, 4.1),
                           y=c(4.9, 4.9, 4.6), yend=c(4.9, 4.6, 4.6), size=1.5)+      
                  annotate("text", 2.5, 4.95, label="*", size=24)

fix.plot  <- ggplot(plot.dat, aes(cond, fixdur, fill=fillon))+
                    ylab("Fixation Duration (ms)")+
                    pvals+
                    theme(legend.position = c(-1,-1))+
                    annotate("text", 3.4, 175, label="*", size=24)+
                    annotate("text", 2.7, 200, label="*", size=24, alpha=.4)+
                    annotate("text", 2.7, 208, label="X", size=12, color="red",
                             alpha=.6)+
                    annotate("segment", x=1.9, xend=3.9, y=200, yend=200,
                             size=1.5, alpha=.4)        
  

# Display and save
slat.plot
amp.plot
fix.plot

kPsize <- 20
ggsave("figs/lat.tiff", slat.plot, height=kPsize, width=kPsize, units="cm",
       dpi = 600)
ggsave("figs/amp.tiff", amp.plot, height = kPsize, width = kPsize,
       units = "cm",
       dpi = 600)
ggsave("figs/fix.tiff", fix.plot, height = kPsize, width = kPsize,
       units = "cm",
       dpi = 600)

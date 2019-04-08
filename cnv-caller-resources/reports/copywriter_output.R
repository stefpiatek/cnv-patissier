library(ggplot2)
library(dplyr)
library(plotROC)
library(here)
library(ggsci)

load(here::here("../../output/actual_copywriter_threshold.RData"))

ggplot(copywriter_data, aes(d = roc_outcome, m = seg_mean)) +
  geom_abline(slope=1, intercept=0, colour="grey") +
  geom_roc(linealpha = 0.3, colour = "#005aac", cutoffs.at = c(0:10 * 0.1, 2)) +
  theme_bw() +
  ylab("Sensitivity") +
  xlab("1 - Specificity") + 
  ggsave(here::here("../../output/plots_tables/copywriter_threshold.pdf"), height = 5, width = 5) +
  ggsave(here::here("../../output/plots_tables/copywriter_threshold.png"), height = 5, width = 5)
  

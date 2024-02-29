library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
preds <- args[1]
outputDir <- args[2]
classifyName <- args[3]

data<-read_tsv(preds)

mb_adj<-1000000
data %>%
  dplyr::mutate(window_middle=(window_start+window_end)/2) %>%
  ggplot() + 
  aes(x=window_middle/mb_adj,y=`prob(sweep)`) +
  geom_point(size=1) +
  xlab("chromosome position (Mb)") +
  theme_bw() +
  ylab("sweep probability")
ggsave(paste0(outputDir,"/plots/",classifyName,"classification_plot.png"))

